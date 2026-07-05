#!/usr/bin/env python3
"""
run_gen_fixings_chunked.py -- divide-and-conquer runner for the gen_fixings scan.
=================================================================================

WHY: a single gen_fixings run over the full bronze CSV at calibration_n=10000
produced a 5-8 GB output CSV and a ~5.9 GB heap (12 OMP threads each holding a
10000-point cloud); RAM hit 94%, swap thrashed, and the host froze (see
`docs`/the jun-2026 incident). This driver runs the bronze rows in SEQUENTIAL
chunks, each a SEPARATE orchestrator subprocess, so:

  - heap is fully released between parts (fresh process per chunk), and
  - each part writes a bounded (~1.4 GB) output instead of one 8 GB file,

keeping peak memory well under RAM so calibration_n=10000 actually COMPLETES.

It does NOT change any C++/orchestrator code -- it is a thin loop around the
existing `python -m dihiggs.app.orchestrator` CLI, plus a concat of the per-chunk
outputs. Run it under ./safe_run.sh for the cgroup backstop.

KEY CORRECTNESS POINTS
----------------------
* RNG alignment: the engine seeds each row as ``rng_seed + local_row_index``
  (GenScanWithFixings.cpp). Each chunk is therefore launched with
  ``--rng-seed = base_seed + chunk_global_start_row`` so per-row seeds equal
  ``base_seed + global_index`` -- the chunked output is content-identical to the
  unsplit run (row order differs, which nothing downstream depends on).
* Resume gate: the orchestrator's own skip logic keys on "output CSV exists",
  NOT on completion, so a chunk killed mid-write leaves a truncated-but-nonempty
  CSV. We instead gate on the orchestrator's completion signal
  ``scan_meta.json`` ``event == "done"`` (written only on success):
    - done            -> skip (no subprocess)
    - incomplete      -> re-run WITH --force (so the orchestrator overwrites the
                         truncated file instead of skipping it)
    - absent          -> run normally
  If any chunk is not "done" after its run, we ABORT before concat so
  silver_all.csv is never built from a partial set. Re-running resumes.

USAGE
-----
  python3 scripts/run_gen_fixings_chunked.py \
      --bronze-csv data/lam1eq1_allchris_var1000_2026jun/concat/bronze_all.csv \
      --calibration-n 10000 --rows-per-chunk 100 \
      --out-root data/lam1eq1_allchris_var10000_2026jun/chunked \
      --campaign lam1eq1_allchris_var10000_2026jun

  # under the resource backstop, plotting on completion:
  ./safe_run.sh --name var10000_chunked \
      --then "python3 scripts/compare_chris_vs_2hdmc_lam1eq1.py \
                --silver-csv data/lam1eq1_allchris_var10000_2026jun/chunked/silver_all.csv \
                --outdir data/lam1eq1_allchris_var10000_2026jun/chunked/comparison" \
      -- python3 scripts/run_gen_fixings_chunked.py ... --no-plot
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_BRONZE = "data/lam1eq1_allchris_var1000_2026jun/concat/bronze_all.csv"
DEFAULT_EXEC = "dihiggs/app/GenScanWithFixings"
DEFAULT_PLOTTER = "scripts/compare_chris_vs_2hdmc_lam1eq1.py"
LAKE_NAME = "dihiggs_lake"
RUN_NAME = "run01"


# ---------------------------------------------------------------------------
@dataclass
class Chunk:
    index: int
    start_row: int       # global 0-based index of this chunk's first data row
    n_rows: int
    csv_path: Path

    @property
    def tag(self) -> str:
        return f"chunk_{self.index:02d}"


def log(msg: str) -> None:
    print(f"[chunked] {msg}", flush=True)


# ---------------------------------------------------------------------------
def split_bronze(bronze_path: Path, work_dir: Path, rows_per_chunk: int) -> List[Chunk]:
    """Split a bronze CSV (header + data rows) into contiguous chunk files."""
    lines = bronze_path.read_text(encoding="utf-8").splitlines(keepends=True)
    if len(lines) < 2:
        sys.exit(f"[chunked] ERROR: bronze CSV has no data rows: {bronze_path}")
    header, data = lines[0], lines[1:]
    n = len(data)
    work_dir.mkdir(parents=True, exist_ok=True)

    chunks: List[Chunk] = []
    for i, start in enumerate(range(0, n, rows_per_chunk)):
        slice_rows = data[start:start + rows_per_chunk]
        if not slice_rows:
            continue  # never emit an empty chunk
        chunk = Chunk(index=i, start_row=start, n_rows=len(slice_rows),
                      csv_path=work_dir / f"chunk_{i:02d}.csv")
        with chunk.csv_path.open("w", encoding="utf-8") as fh:
            fh.write(header if header.endswith("\n") else header + "\n")
            fh.writelines(slice_rows)
        chunks.append(chunk)
    log(f"split {n} data rows into {len(chunks)} chunk(s) of <= {rows_per_chunk} rows "
        f"-> {work_dir}")
    return chunks


def chunk_outdir(out_root: Path, chunk: Chunk) -> Path:
    return out_root / chunk.tag / "lake"


def _find_one(outdir: Path, name_glob: str) -> Optional[Path]:
    """Return the single match for a recursive glob under outdir, else None.
    Fails loudly on >1 match (a stray/duplicated run would corrupt the concat)."""
    matches = sorted((outdir / LAKE_NAME).glob(f"**/{name_glob}"))
    if len(matches) > 1:
        sys.exit(f"[chunked] ERROR: expected 1 '{name_glob}' under {outdir}, "
                 f"found {len(matches)}: {matches}")
    return matches[0] if matches else None


def chunk_status(outdir: Path) -> tuple[str, Optional[Path]]:
    """Classify a chunk by the orchestrator's completion signal.

    Returns (status, silver_path) where status is 'done' | 'incomplete' | 'absent'.
    'done' requires scan_meta.json event=="done" AND a non-empty silver CSV.
    """
    silver = _find_one(outdir, "scan_tb_*.csv")
    meta = _find_one(outdir, "scan_meta.json")
    if meta is not None:
        try:
            event = json.loads(meta.read_text(encoding="utf-8")).get("event")
        except (OSError, json.JSONDecodeError):
            event = None
        if event == "done" and silver is not None and silver.stat().st_size > 0:
            return "done", silver
    if silver is not None:
        return "incomplete", silver
    return "absent", None


def run_chunk(chunk: Chunk, args: argparse.Namespace, force: bool) -> int:
    """Launch one orchestrator subprocess for this chunk (fresh process)."""
    outdir = chunk_outdir(Path(args.out_root), chunk)
    seed = args.base_rng_seed + chunk.start_row
    cmd = [
        sys.executable, "-m", "dihiggs.app.orchestrator",
        "--engine", "gen_fixings",
        "--exec", args.exec_path,
        "--campaign", args.campaign,
        "--run-name", RUN_NAME,
        "--lake-name", LAKE_NAME,
        "--bronze-csv", str(chunk.csv_path),
        "--calibration-n", str(args.calibration_n),
        "--calibration-frac", repr(args.calibration_frac),
        "--rng-seed", str(seed),
        "--tanbeta", "10",
        "--outdir", str(outdir),
    ]
    if force:
        cmd.append("--force")
    if args.dry_run:
        cmd.append("--dry-run")
    log(f"{chunk.tag}: rows[{chunk.start_row}:{chunk.start_row + chunk.n_rows}] "
        f"seed={seed} force={force} -> {outdir}")
    return subprocess.run(cmd, cwd=REPO_ROOT).returncode


def concat_silvers(silver_paths: List[Path], dest: Path) -> int:
    """Concatenate per-chunk silver CSVs into one (header once, data appended)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    header: Optional[str] = None
    total_rows = 0
    with dest.open("w", encoding="utf-8") as out:
        for sp in silver_paths:
            with sp.open("r", encoding="utf-8") as fh:
                this_header = fh.readline()
                if not this_header:
                    continue
                if header is None:
                    header = this_header
                    out.write(header)
                elif this_header != header:
                    sys.exit(f"[chunked] ERROR: header mismatch in {sp}\n"
                             f"  expected: {header!r}\n  got:      {this_header!r}")
                for line in fh:
                    out.write(line)
                    total_rows += 1
    log(f"concatenated {len(silver_paths)} chunk(s) -> {dest} ({total_rows} data rows)")
    return total_rows


def maybe_plot(silver_all: Path, out_root: Path) -> None:
    comp_dir = out_root / "comparison"
    cmd = [sys.executable, DEFAULT_PLOTTER,
           "--silver-csv", str(silver_all), "--outdir", str(comp_dir)]
    log(f"plotting -> {comp_dir}")
    rc = subprocess.run(cmd, cwd=REPO_ROOT).returncode
    if rc != 0:
        sys.exit(f"[chunked] ERROR: plotter exited {rc}")


# ---------------------------------------------------------------------------
def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bronze-csv", default=DEFAULT_BRONZE,
                   help="Bronze CSV to split and scan (header + data rows).")
    p.add_argument("--campaign", required=True, help="Campaign label for the runs.")
    p.add_argument("--out-root", required=True,
                   help="Root dir for per-chunk lakes, silver_all.csv, comparison/.")
    p.add_argument("--calibration-n", type=int, required=True,
                   help="gen_fixings calibration candidates per bronze row.")
    p.add_argument("--calibration-frac", type=float, default=0.1,
                   help="Calibration variation fraction (default 0.10).")
    p.add_argument("--base-rng-seed", type=int, default=0,
                   help="Base RNG seed; chunk i uses base + global_start_row.")
    grp = p.add_mutually_exclusive_group()
    grp.add_argument("--rows-per-chunk", type=int, default=100,
                     help="Bronze rows per chunk (default 100 -> K=6 on 600 rows).")
    grp.add_argument("--chunks", type=int, default=None,
                     help="Number of chunks K (alternative to --rows-per-chunk).")
    p.add_argument("--exec", dest="exec_path", default=DEFAULT_EXEC,
                   help="Path to the GenScanWithFixings binary.")
    p.add_argument("--work-dir", default=None,
                   help="Where chunk CSVs are written (default <out-root>/_chunks).")
    plot = p.add_mutually_exclusive_group()
    plot.add_argument("--plot", dest="plot", action="store_true", default=True,
                      help="Run the comparison plotter on silver_all.csv (default).")
    plot.add_argument("--no-plot", dest="plot", action="store_false",
                      help="Skip plotting (e.g. when safe_run --then does it).")
    p.add_argument("--dry-run", action="store_true",
                   help="Pass --dry-run to each chunk (no C++), skip concat/plot.")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    bronze = Path(args.bronze_csv).expanduser().resolve()
    if not bronze.exists():
        sys.exit(f"[chunked] ERROR: bronze CSV not found: {bronze}")
    out_root = Path(args.out_root).expanduser().resolve()
    work_dir = Path(args.work_dir).expanduser().resolve() if args.work_dir \
        else out_root / "_chunks"

    n_data = max(0, len(bronze.read_text(encoding="utf-8").splitlines()) - 1)
    rows_per_chunk = (math.ceil(n_data / args.chunks)
                      if args.chunks else args.rows_per_chunk)
    if rows_per_chunk < 1:
        sys.exit("[chunked] ERROR: rows-per-chunk resolved to < 1")

    chunks = split_bronze(bronze, work_dir, rows_per_chunk)

    if args.dry_run:
        for chunk in chunks:
            run_chunk(chunk, args, force=False)
        log("dry-run complete (no concat/plot).")
        return 0

    # --- run each chunk sequentially, gating on the done-signal --------------
    silver_paths: List[Path] = []
    for chunk in chunks:
        outdir = chunk_outdir(out_root, chunk)
        status, silver = chunk_status(outdir)
        if status == "done":
            log(f"{chunk.tag}: skip (already done) -> {silver}")
            silver_paths.append(silver)  # type: ignore[arg-type]
            continue
        rc = run_chunk(chunk, args, force=(status == "incomplete"))
        status, silver = chunk_status(outdir)
        if rc != 0 or status != "done" or silver is None:
            log(f"{chunk.tag}: FAILED (rc={rc}, status={status}). Aborting before "
                f"concat so silver_all.csv is never partial. Re-run to resume.")
            log("  If a single chunk keeps OOM-ing, more (smaller) chunks will "
                "NOT help: peak heap is set by thread count (per-thread 2HDMC "
                "state), not row count. Lower --threads (e.g. safe_run "
                "--reserve); --calibration-n mainly controls output size.")
            return 1
        silver_paths.append(silver)

    # --- concat + optional plot ---------------------------------------------
    silver_all = out_root / "silver_all.csv"
    concat_silvers(silver_paths, silver_all)
    if args.plot:
        maybe_plot(silver_all, out_root)
    log(f"DONE. {len(silver_paths)} chunk(s) -> {silver_all}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
