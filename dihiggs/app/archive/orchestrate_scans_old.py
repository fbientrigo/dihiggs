#!/usr/bin/env python3
"""
Orchestrate 2HDM parameter scans and persist results in a CERNBox-synced "mini data lake".

This script is designed to be:
- **Resumable**: skips existing outputs unless `--force` is used.
- **Multi-machine friendly**: each run writes into a unique run directory (`run_id`),
  so different machines won't overwrite each other.
- **Analysis-ready**: writes a `run_manifest.json` capturing *all* fixed hyperparameters
  and scan-grid settings (including ranges and point counts).

Folder layout (default):
    /mnt/c/Users/Asus/cern_db/dihiggs_lake/
        campaign=<campaign>/
            fixed_sinba=<...>_l6=<...>_l7=<...>_mA=<...>/
                run_<run_id>_host=<host>_git=<sha>/
                    run_manifest.json
                    orchestrator.log
                    task_summary.jsonl
                    tb_<...>/
                        scan_tb_<...>.csv

Notes
-----
- CVMFS (/cvmfs) is read-only and should not be used for outputs.
- For safety/portability, avoid writing outputs inside your repo; instead point `--outdir`
  to a CERNBox-synced folder (default already does this).
- This script does NOT use sudo. Install prefixes should be writable by the user.

Examples
--------
Basic (defaults to CERNBox outdir):
    python3 orchestrate_scans.py --exec ./PhysScanWithFixings --campaign scan_999

Choose fixed hyperparameters and tan(beta) list:
    python3 orchestrate_scans.py \
        --campaign scan_999 \
        --sin-ba 0.995 --lambda6 0.10 --lambda7 0.00 --mA 300 \
        --tanbeta 10000,15000,20000

Customize scan grid (mphi and lambda1 ranges, N points):
    python3 orchestrate_scans.py \
        --mphi-min 130 --mphi-max 290 --n-mphi 15 \
        --lam1-min 0 --lam1-max 12 --n-lam1 666

Force overwrite of existing CSVs:
    python3 orchestrate_scans.py --force

Dry-run (print commands, create folders + manifest, but don't execute):
    python3 orchestrate_scans.py --dry-run
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import socket
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# -----------------------------------------------------------------------------
# Defaults (safe + convenient for WSL + CERNBox Desktop on Windows)
# -----------------------------------------------------------------------------
DEFAULT_CERNBOX_OUTDIR = "/mnt/c/Users/Asus/cern_db"
DEFAULT_LAKE_DIRNAME = "dihiggs_lake"

# Original defaults (kept as defaults, but now configurable via CLI)
DEFAULT_MPHI_MIN = 130.0
DEFAULT_MPHI_MAX = 290.0
DEFAULT_N_MPHI = 15

DEFAULT_LAM1_MIN = 0.0
DEFAULT_LAM1_MAX = 12.0
DEFAULT_N_LAM1 = 666

DEFAULT_MA = 300.0
DEFAULT_SIN_BA = 0.995
DEFAULT_L6 = 0.1
DEFAULT_L7 = 0.0

DEFAULT_TANBETA_LIST = [10000.0, 15000.0, 20000.0]


# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------
def utc_now_iso() -> str:
    """Return current UTC timestamp in ISO-8601 format."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def utc_now_compact() -> str:
    """Return compact UTC timestamp suitable for folder naming."""
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def log_message(msg: str, log_file: Path) -> None:
    """
    Write a timestamped message to stdout and append to a log file.

    Parameters
    ----------
    msg:
        Message to log.
    log_file:
        Path to the orchestrator log file.
    """
    ts = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    full = f"{ts} {msg}"
    print(full, flush=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("a", encoding="utf-8") as f:
        f.write(full + "\n")


def safe_write_json(path: Path, payload: Dict[str, Any]) -> None:
    """
    Atomically write a JSON file.

    Parameters
    ----------
    path:
        Destination JSON path.
    payload:
        JSON-serializable dictionary.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=False)
        f.write("\n")
    tmp.replace(path)


def append_jsonl(path: Path, record: Dict[str, Any]) -> None:
    """
    Append a single JSON record to a JSON Lines file.

    Parameters
    ----------
    path:
        JSONL file path.
    record:
        JSON-serializable dictionary record.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record) + "\n")


def parse_csv_floats(s: str) -> List[float]:
    """
    Parse a comma-separated list of floats.

    Examples
    --------
    "100,200,500" -> [100.0, 200.0, 500.0]
    """
    if not s.strip():
        return []
    parts = [p.strip() for p in s.split(",")]
    out: List[float] = []
    for p in parts:
        if p == "":
            continue
        out.append(float(p))
    return out


def sanitize_for_path(s: str) -> str:
    """
    Sanitize a string for safe path components.

    Replaces:
    - '.' -> 'p'
    - '-' -> 'm'
    - any other non-alnum/_/= -> '_'
    """
    s = s.replace(".", "p").replace("-", "m")
    s = re.sub(r"[^A-Za-z0-9_=]+", "_", s)
    return s.strip("_")


def format_float_tag(x: float, ndp: int = 4) -> str:
    """
    Format a float for stable tagging in folder names.

    Parameters
    ----------
    x:
        Value to format.
    ndp:
        Decimal places.

    Returns
    -------
    str
        A compact, deterministic string representation.
    """
    # Use fixed decimals for stability in folder naming
    return f"{x:.{ndp}f}"


def format_tanbeta_folder(tb: float) -> Tuple[str, str]:
    """
    Create a folder name and file tag for tan(beta).

    Handles integers with zero padding (tb_00010), and non-integers safely.

    Parameters
    ----------
    tb:
        tan(beta) value.

    Returns
    -------
    folder_name:
        Directory name for this tan(beta).
    tag:
        Filename tag for this tan(beta).
    """
    if abs(tb - round(tb)) < 1e-12 and tb >= 0:
        tb_int = int(round(tb))
        # Pad at least 5 digits; if larger, it will expand naturally.
        folder = f"tb_{tb_int:05d}"
        tag = str(tb_int)
        return folder, tag

    # Non-integer or negative: safe representation
    raw = f"{tb:.6g}"
    safe = sanitize_for_path(raw)
    folder = f"tb_{safe}"
    tag = safe
    return folder, tag


def detect_git_info(start_dir: Path) -> Dict[str, Optional[str]]:
    """
    Best-effort git metadata (commit + dirty status).

    Parameters
    ----------
    start_dir:
        A directory inside (or near) the git repo.

    Returns
    -------
    dict
        Keys: "commit", "commit_short", "is_dirty", "repo_root"
    """
    def run_git(args: List[str]) -> Optional[str]:
        try:
            r = subprocess.run(
                ["git"] + args,
                cwd=str(start_dir),
                capture_output=True,
                text=True,
                check=False,
            )
            if r.returncode != 0:
                return None
            return r.stdout.strip()
        except Exception:
            return None

    repo_root = run_git(["rev-parse", "--show-toplevel"])
    commit = run_git(["rev-parse", "HEAD"])
    dirty = run_git(["status", "--porcelain"])

    commit_short = commit[:7] if commit else None
    is_dirty = None if dirty is None else ("yes" if dirty != "" else "no")

    return {
        "repo_root": repo_root,
        "commit": commit,
        "commit_short": commit_short,
        "is_dirty": is_dirty,
    }


@dataclass(frozen=True)
class ScanGrid:
    """Container for scan-grid parameters."""
    mphi_min: float
    mphi_max: float
    n_mphi: int
    lam1_min: float
    lam1_max: float
    n_lam1: int


@dataclass(frozen=True)
class FixedParams:
    """Container for fixed (experiment-level) hyperparameters."""
    mA: float
    sin_ba: float
    lambda6: float
    lambda7: float


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    p = argparse.ArgumentParser(
        description="Orchestrator for 2HDM scans (C++ driver) with CERNBox data-lake output."
    )

    # Executable + runtime
    p.add_argument(
        "--exec",
        type=str,
        default="./PhysScanWithFixings",
        help="Path to the compiled C++ executable.",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=None,
        help="OMP_NUM_THREADS (optional). If omitted, use default (all cores).",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing CSV outputs (otherwise skip).",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Create folders + manifest, print commands, but do not execute the C++ binary.",
    )

    # Output organization
    p.add_argument(
        "--outdir",
        type=str,
        default=DEFAULT_CERNBOX_OUTDIR,
        help=(
            "Root output directory (default points to CERNBox-synced folder in WSL). "
            f"Default: {DEFAULT_CERNBOX_OUTDIR}"
        ),
    )
    p.add_argument(
        "--lake-name",
        type=str,
        default=DEFAULT_LAKE_DIRNAME,
        help=f"Subfolder under outdir used as data-lake root. Default: {DEFAULT_LAKE_DIRNAME}",
    )
    p.add_argument(
        "--campaign",
        type=str,
        default="scan",
        help="Human label for this scan campaign (e.g. scan_999).",
    )
    p.add_argument(
        "--run-name",
        type=str,
        default=None,
        help=(
            "Optional run name. If omitted, uses UTC timestamp + hostname (+ git short sha if available)."
        ),
    )

    # Fixed hyperparameters (experiment-level)
    p.add_argument("--mA", type=float, default=DEFAULT_MA, help="Fixed mA (also mHp).")
    p.add_argument("--sin-ba", type=float, default=DEFAULT_SIN_BA, help="Fixed sin(b-a).")
    p.add_argument("--lambda6", type=float, default=DEFAULT_L6, help="Fixed lambda_6.")
    p.add_argument("--lambda7", type=float, default=DEFAULT_L7, help="Fixed lambda_7.")

    # Scan grid
    p.add_argument("--mphi-min", type=float, default=DEFAULT_MPHI_MIN, help="m_phi min.")
    p.add_argument("--mphi-max", type=float, default=DEFAULT_MPHI_MAX, help="m_phi max.")
    p.add_argument("--n-mphi", type=int, default=DEFAULT_N_MPHI, help="Number of m_phi points.")

    p.add_argument("--lam1-min", type=float, default=DEFAULT_LAM1_MIN, help="lambda_1 min.")
    p.add_argument("--lam1-max", type=float, default=DEFAULT_LAM1_MAX, help="lambda_1 max.")
    # Keep backward compatibility with your old flag spelling:
    p.add_argument(
        "--n-lam1",
        "--N-lam1",
        dest="n_lam1",
        type=int,
        default=DEFAULT_N_LAM1,
        help="Number of lambda_1 points. Default: 666. Use 10 for tests.",
    )

    # TanBeta list
    p.add_argument(
        "--tanbeta",
        type=str,
        default=",".join(str(int(x)) for x in DEFAULT_TANBETA_LIST),
        help="Comma-separated list of tan(beta) values, e.g. '100,500,1000'.",
    )

    return p


# -----------------------------------------------------------------------------
# Main orchestrator
# -----------------------------------------------------------------------------
def build_run_dir(
    outdir: Path,
    lake_name: str,
    campaign: str,
    fixed: FixedParams,
    run_name: str,
) -> Path:
    """
    Build the run directory path.

    The directory is split into:
      - campaign directory
      - fixed-parameter directory
      - run directory

    Returns
    -------
    Path
        Full path to run directory.
    """
    lake_root = outdir / lake_name
    campaign_dir = lake_root / f"campaign={sanitize_for_path(campaign)}"

    fixed_dir_name = (
        f"fixed_sinba={sanitize_for_path(format_float_tag(fixed.sin_ba, 4))}"
        f"_l6={sanitize_for_path(format_float_tag(fixed.lambda6, 4))}"
        f"_l7={sanitize_for_path(format_float_tag(fixed.lambda7, 4))}"
        f"_mA={sanitize_for_path(format_float_tag(fixed.mA, 1))}"
    )
    fixed_dir = campaign_dir / fixed_dir_name

    run_dir = fixed_dir / sanitize_for_path(run_name)
    return run_dir


def parse_cpp_stats(stdout: str) -> Tuple[Optional[int], Optional[str]]:
    """
    Parse key stats from the C++ stdout.

    Looks for:
      - "Total Attempts: <int>"
      - "TRIPLE_OK_POINTS <something>"

    Returns
    -------
    attempts:
        Parsed integer attempts if found, else None.
    triple_ok:
        Parsed string for TRIPLE_OK_POINTS if found, else None.
    """
    attempts: Optional[int] = None
    triple_ok: Optional[str] = None

    for line in stdout.splitlines():
        if "TRIPLE_OK_POINTS" in line:
            parts = line.split()
            if parts:
                triple_ok = parts[-1]
        if "Total Attempts:" in line:
            try:
                attempts = int(line.split(":")[-1].strip())
            except Exception:
                pass

    return attempts, triple_ok


def main() -> int:
    args = build_parser().parse_args()

    # Resolve executable path
    executable = Path(args.exec).expanduser().resolve()
    if not executable.exists():
        print(f"[ERROR] Executable not found: {executable}", file=sys.stderr)
        return 2

    # Prepare scan settings
    grid = ScanGrid(
        mphi_min=args.mphi_min,
        mphi_max=args.mphi_max,
        n_mphi=args.n_mphi,
        lam1_min=args.lam1_min,
        lam1_max=args.lam1_max,
        n_lam1=args.n_lam1,
    )
    fixed = FixedParams(
        mA=args.mA,
        sin_ba=args.sin_ba,
        lambda6=args.lambda6,
        lambda7=args.lambda7,
    )
    tanbeta_list = parse_csv_floats(args.tanbeta)
    if not tanbeta_list:
        print("[ERROR] Empty --tanbeta list.", file=sys.stderr)
        return 2

    # Environment (OpenMP)
    env = os.environ.copy()
    if args.threads is not None:
        env["OMP_NUM_THREADS"] = str(args.threads)

    # Metadata (host + git)
    host = socket.gethostname()
    git = detect_git_info(Path.cwd())

    # run name
    if args.run_name:
        run_name = args.run_name
    else:
        # Example: run_20260126T182031Z_host=zephyrus_git=abc1234
        sha = git.get("commit_short") or "nogit"
        run_name = f"run_{utc_now_compact()}_host={host}_git={sha}"

    # Output directories
    outdir = Path(args.outdir).expanduser().resolve()
    run_dir = build_run_dir(
        outdir=outdir,
        lake_name=args.lake_name,
        campaign=args.campaign,
        fixed=fixed,
        run_name=run_name,
    )
    run_dir.mkdir(parents=True, exist_ok=True)

    log_path = run_dir / "orchestrator.log"
    jsonl_path = run_dir / "task_summary.jsonl"
    manifest_path = run_dir / "run_manifest.json"

    # Start logging
    log_message("[INIT] Starting scan orchestration.", log_path)
    log_message(f"[PATH] run_dir = {run_dir}", log_path)
    log_message(f"[EXEC] executable = {executable}", log_path)

    if args.threads is not None:
        log_message(f"[CONF] OMP_NUM_THREADS = {args.threads}", log_path)
    log_message(f"[CONF] dry_run = {args.dry_run}", log_path)
    log_message(f"[CONF] force   = {args.force}", log_path)

    log_message(
        f"[CONF] fixed: mA={fixed.mA}, sin(b-a)={fixed.sin_ba}, l6={fixed.lambda6}, l7={fixed.lambda7}",
        log_path,
    )
    log_message(
        f"[CONF] grid: m_phi=[{grid.mphi_min}, {grid.mphi_max}] (N={grid.n_mphi}), "
        f"lam1=[{grid.lam1_min}, {grid.lam1_max}] (N={grid.n_lam1})",
        log_path,
    )
    log_message(f"[CONF] tanbeta_list = {tanbeta_list}", log_path)

    # Write initial manifest
    manifest: Dict[str, Any] = {
        "created_utc": utc_now_iso(),
        "campaign": args.campaign,
        "run_name": run_name,
        "paths": {
            "outdir": str(outdir),
            "lake_name": args.lake_name,
            "run_dir": str(run_dir),
            "log": str(log_path),
            "task_summary_jsonl": str(jsonl_path),
        },
        "host": {
            "hostname": host,
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "user": os.environ.get("USER") or os.environ.get("USERNAME"),
        },
        "git": git,
        "runtime": {
            "executable": str(executable),
            "dry_run": args.dry_run,
            "force": args.force,
            "omp_num_threads": args.threads,
        },
        "fixed_params": {
            "mA": fixed.mA,
            "sin_ba": fixed.sin_ba,
            "lambda6": fixed.lambda6,
            "lambda7": fixed.lambda7,
        },
        "scan_grid": {
            "mphi_min": grid.mphi_min,
            "mphi_max": grid.mphi_max,
            "n_mphi": grid.n_mphi,
            "lam1_min": grid.lam1_min,
            "lam1_max": grid.lam1_max,
            "n_lam1": grid.n_lam1,
            "tanbeta_list": tanbeta_list,
        },
        "summary": {},  # filled at end
    }
    safe_write_json(manifest_path, manifest)
    log_message(f"[MANI] Wrote manifest: {manifest_path}", log_path)

    # Run loop
    total_tasks = len(tanbeta_list)
    grand_t0 = time.time()

    total_attempts = 0
    tasks_ok = 0
    tasks_skipped = 0
    tasks_failed = 0

    for idx, tb in enumerate(tanbeta_list, start=1):
        folder_name, tb_tag = format_tanbeta_folder(tb)
        current_dir = run_dir / folder_name
        current_dir.mkdir(parents=True, exist_ok=True)

        output_csv = current_dir / f"scan_tb_{tb_tag}.csv"

        # resumability
        if output_csv.exists() and not args.force:
            tasks_skipped += 1
            log_message(f"[SKIP] ({idx}/{total_tasks}) tanbeta={tb} exists: {output_csv}", log_path)
            append_jsonl(
                jsonl_path,
                {
                    "event": "skip",
                    "utc": utc_now_iso(),
                    "index": idx,
                    "total": total_tasks,
                    "tanbeta": tb,
                    "output_csv": str(output_csv),
                    "reason": "exists_and_force_false",
                },
            )
            continue

        log_message(f"[RUN ] ({idx}/{total_tasks}) tanbeta={tb} -> {output_csv.name}", log_path)

        cmd = [
            str(executable),
            f"{grid.mphi_min:.4f}",
            f"{grid.mphi_max:.4f}",
            str(grid.n_mphi),
            f"{grid.lam1_min:.4f}",
            f"{grid.lam1_max:.4f}",
            str(grid.n_lam1),
            f"{fixed.mA:.4f}",
            f"{fixed.sin_ba:.4f}",
            f"{tb:.4f}",
            f"{fixed.lambda6:.4f}",
            f"{fixed.lambda7:.4f}",
            str(output_csv),
        ]

        if args.dry_run:
            log_message(f"[DRY ] cmd = {' '.join(cmd)}", log_path)
            append_jsonl(
                jsonl_path,
                {
                    "event": "dry_run",
                    "utc": utc_now_iso(),
                    "index": idx,
                    "total": total_tasks,
                    "tanbeta": tb,
                    "cmd": cmd,
                    "output_csv": str(output_csv),
                },
            )
            continue

        t0 = time.time()
        try:
            result = subprocess.run(cmd, env=env, capture_output=True, text=True)
            elapsed = time.time() - t0

            attempts, triple_ok = parse_cpp_stats(result.stdout)

            record: Dict[str, Any] = {
                "utc": utc_now_iso(),
                "index": idx,
                "total": total_tasks,
                "tanbeta": tb,
                "output_csv": str(output_csv),
                "elapsed_sec": elapsed,
                "returncode": result.returncode,
                "attempts": attempts,
                "triple_ok_points": triple_ok,
            }

            if result.returncode == 0:
                tasks_ok += 1
                if attempts is not None:
                    total_attempts += attempts
                log_message(
                    f"[DONE] tanbeta={tb} in {elapsed:.2f}s | attempts={attempts} | tripleOK={triple_ok}",
                    log_path,
                )
                record["event"] = "done"
            else:
                tasks_failed += 1
                log_message(f"[FAIL] tanbeta={tb} code={result.returncode}", log_path)
                # Keep stderr in a file for debugging (avoid polluting orchestrator.log)
                stderr_path = current_dir / "stderr.log"
                stdout_path = current_dir / "stdout.log"
                stderr_path.write_text(result.stderr or "", encoding="utf-8")
                stdout_path.write_text(result.stdout or "", encoding="utf-8")
                log_message(f"[FAIL] wrote stdout/stderr logs to: {stdout_path.name}, {stderr_path.name}", log_path)
                record["event"] = "fail"
                record["stderr_path"] = str(stderr_path)
                record["stdout_path"] = str(stdout_path)

            append_jsonl(jsonl_path, record)

        except Exception as e:
            tasks_failed += 1
            elapsed = time.time() - t0
            log_message(f"[CRASH] tanbeta={tb} exception={e}", log_path)
            append_jsonl(
                jsonl_path,
                {
                    "event": "crash",
                    "utc": utc_now_iso(),
                    "index": idx,
                    "total": total_tasks,
                    "tanbeta": tb,
                    "output_csv": str(output_csv),
                    "elapsed_sec": elapsed,
                    "exception": repr(e),
                    "cmd": cmd,
                },
            )

    # Final summary
    grand_elapsed = time.time() - grand_t0
    log_message("-" * 70, log_path)
    log_message("[FIN ] Orchestration complete.", log_path)
    log_message(f"       Total wall time: {grand_elapsed:.2f} s ({grand_elapsed/60:.2f} min)", log_path)
    log_message(f"       Tasks: ok={tasks_ok}, skipped={tasks_skipped}, failed={tasks_failed}", log_path)
    log_message(f"       Total attempts (parsed): {total_attempts}", log_path)

    if total_attempts > 0 and grand_elapsed > 0:
        sec_per_point = grand_elapsed / total_attempts
        points_per_hour = total_attempts / (grand_elapsed / 3600.0)
        log_message(f"       Avg speed: {sec_per_point:.6f} s/attempt", log_path)
        log_message(f"       Throughput: {int(points_per_hour)} attempts/hour", log_path)
    else:
        log_message("       No new attempts parsed (all skipped/dry-run/failed?).", log_path)

    # Update manifest with summary
    manifest["summary"] = {
        "finished_utc": utc_now_iso(),
        "wall_time_sec": grand_elapsed,
        "tasks_total": total_tasks,
        "tasks_ok": tasks_ok,
        "tasks_skipped": tasks_skipped,
        "tasks_failed": tasks_failed,
        "total_attempts_parsed": total_attempts,
    }
    safe_write_json(manifest_path, manifest)
    log_message(f"[MANI] Updated manifest summary: {manifest_path}", log_path)
    log_message("-" * 70, log_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
