#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


REQUIRED_HEADER_COLUMNS = {
    "m_phi",
    "mA",
    "lambda6",
    "lambda7",
    "m12",
    "sin_ba",
    "tan_beta",
    "positivity_ok",
    "unitarity_ok",
    "perturbativity_ok",
}


@dataclass
class CsvCandidate:
    path: str
    campaign: str
    size_bytes: int
    mtime: float
    header_match: bool
    scan_tb_like: bool


@dataclass
class RunResult:
    run_index: int
    campaign: str
    input_csv: str
    mtime: float
    size_bytes: int
    mode: str
    requested_samples: int
    actual_samples_written: str
    input_rows_parsed: str
    parseable_rows: str
    eligible_rows: str
    baseline_failures: str
    probe_failures: str
    max_abs_error: str
    max_rel_error: str
    max_scaled_error: str
    warning_rate_abs: str
    warning_rate_rel: str
    warning_rate_scaled: str
    elapsed_sec: float
    returncode: int
    output_csv: str
    stdout_log: str
    stderr_log: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark CalcLam1ScanFromLake across official lake campaigns"
    )
    parser.add_argument("--lake-root", required=True, type=Path)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)

    parser.add_argument("--mode", default="triple_ok")
    parser.add_argument("--n-campaigns", type=int, default=10)
    parser.add_argument("--n-files-per-campaign", type=int, default=1)
    parser.add_argument("--samples-per-file", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=123456)
    parser.add_argument("--campaign-regex", default="")
    parser.add_argument("--max-files-scan", type=int, default=5000)
    parser.add_argument(
        "--require-header-match",
        dest="require_header_match",
        action="store_true",
        default=True,
        help="Require PhysScanWithFixings header columns",
    )
    parser.add_argument(
        "--no-require-header-match",
        dest="require_header_match",
        action="store_false",
        help="Allow non-matching header files",
    )
    return parser.parse_args()


def sanitize_name(value: str, max_len: int = 80) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    if not cleaned:
        cleaned = "item"
    return cleaned[:max_len]


def campaign_from_path(path: Path) -> str:
    for part in path.parts:
        if part.startswith("campaign="):
            return part
    return "campaign=unknown"


def read_csv_header(path: Path) -> list[str]:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            first = handle.readline().strip()
    except OSError:
        return []
    if not first:
        return []
    return [col.strip() for col in first.split(",")]


def header_matches_schema(columns: Iterable[str]) -> bool:
    colset = set(columns)
    return REQUIRED_HEADER_COLUMNS.issubset(colset)


def discover_csv_candidates(lake_root: Path, max_files_scan: int) -> tuple[list[CsvCandidate], int]:
    candidates: list[CsvCandidate] = []
    scanned_csv = 0
    for root, dirs, files in os.walk(lake_root):
        dirs.sort()
        files.sort()
        for filename in files:
            if not filename.lower().endswith(".csv"):
                continue
            scanned_csv += 1
            full = Path(root) / filename
            try:
                stat = full.stat()
            except OSError:
                continue
            columns = read_csv_header(full)
            candidate = CsvCandidate(
                path=str(full),
                campaign=campaign_from_path(full),
                size_bytes=stat.st_size,
                mtime=stat.st_mtime,
                header_match=header_matches_schema(columns),
                scan_tb_like=bool(re.fullmatch(r"scan_tb_.*\.csv", filename)),
            )
            candidates.append(candidate)
            if scanned_csv >= max_files_scan:
                return candidates, scanned_csv
    return candidates, scanned_csv


def select_candidates(
    candidates: list[CsvCandidate],
    n_campaigns: int,
    n_files_per_campaign: int,
    campaign_regex: str,
    require_header_match: bool,
) -> list[CsvCandidate]:
    filtered = candidates
    if require_header_match:
        filtered = [c for c in filtered if c.header_match]
    if campaign_regex:
        pattern = re.compile(campaign_regex)
        filtered = [c for c in filtered if pattern.search(c.campaign)]

    by_campaign: dict[str, list[CsvCandidate]] = {}
    for candidate in filtered:
        by_campaign.setdefault(candidate.campaign, []).append(candidate)

    campaigns = sorted(
        by_campaign.keys(),
        key=lambda campaign: (
            -max(item.mtime for item in by_campaign[campaign]),
            campaign,
        ),
    )
    selected_campaigns = campaigns[:n_campaigns]

    selected: list[CsvCandidate] = []
    for campaign in selected_campaigns:
        ranked = sorted(
            by_campaign[campaign],
            key=lambda c: (
                0 if c.header_match else 1,
                0 if c.scan_tb_like else 1,
                -c.size_bytes,
                -c.mtime,
                c.path,
            ),
        )
        selected.extend(ranked[:n_files_per_campaign])
    return selected


def derive_seed(global_seed: int, run_index: int, input_csv: str) -> int:
    digest = hashlib.sha256(f"{global_seed}:{run_index}:{input_csv}".encode("utf-8")).digest()
    value = int.from_bytes(digest[:8], byteorder="big", signed=False)
    return value % (2**31 - 1)


def parse_metric(text: str, patterns: list[str]) -> str:
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE | re.MULTILINE)
        if match:
            value = match.group(1).strip()
            value = value.rstrip("%")
            return value
    return ""


def parse_validator_stdout(stdout_text: str) -> dict[str, str]:
    return {
        "input_rows_parsed": parse_metric(
            stdout_text,
            [
                r"^\s*(?:input\s*rows\s*parsed|input\s*rows|rows\s*parsed)\s*[:=]\s*([0-9]+)",
            ],
        ),
        "parseable_rows": parse_metric(
            stdout_text,
            [
                r"^\s*(?:parseable\s*rows|rows\s*parseable)\s*[:=]\s*([0-9]+)",
            ],
        ),
        "eligible_rows": parse_metric(
            stdout_text,
            [
                r"^\s*(?:eligible\s*rows(?:\s*after\s*filtering)?)\s*[:=]\s*([0-9]+)",
            ],
        ),
        "requested_samples": parse_metric(
            stdout_text,
            [
                r"^\s*(?:requested\s*samples)\s*[:=]\s*([0-9]+)",
            ],
        ),
        "actual_samples_written": parse_metric(
            stdout_text,
            [
                r"^\s*(?:actual\s*samples\s*written|samples\s*written)\s*[:=]\s*([0-9]+)",
            ],
        ),
        "baseline_failures": parse_metric(
            stdout_text,
            [r"^\s*(?:baseline\s*failures?)\s*[:=]\s*([0-9]+)"],
        ),
        "probe_failures": parse_metric(
            stdout_text,
            [r"^\s*(?:probe\s*failures?)\s*[:=]\s*([0-9]+)"],
        ),
        "max_abs_error": parse_metric(
            stdout_text,
            [r"^\s*(?:max\s*(?:absolute|abs)\s*error)\s*[:=]\s*([-+0-9.eE]+)"],
        ),
        "max_rel_error": parse_metric(
            stdout_text,
            [r"^\s*(?:max\s*(?:relative|rel)\s*error)\s*[:=]\s*([-+0-9.eE]+)"],
        ),
        "max_scaled_error": parse_metric(
            stdout_text,
            [r"^\s*(?:max\s*scaled\s*error)\s*[:=]\s*([-+0-9.eE]+)"],
        ),
        "warning_rate_abs": parse_metric(
            stdout_text,
            [r"^\s*(?:warning\s*rate\s*(?:abs|absolute))\s*[:=]\s*([-+0-9.eE%]+)"],
        ),
        "warning_rate_rel": parse_metric(
            stdout_text,
            [r"^\s*(?:warning\s*rate\s*rel(?:ative)?)\s*[:=]\s*([-+0-9.eE%]+)"],
        ),
        "warning_rate_scaled": parse_metric(
            stdout_text,
            [r"^\s*(?:warning\s*rate\s*scaled)\s*[:=]\s*([-+0-9.eE%]+)"],
        ),
    }


def run_validator(
    binary: Path,
    mode: str,
    candidate: CsvCandidate,
    samples_per_file: int,
    seed: int,
    run_dir: Path,
    run_index: int,
) -> RunResult:
    run_dir.mkdir(parents=True, exist_ok=True)
    output_csv = run_dir / "validator_output.csv"
    stdout_log = run_dir / "stdout.txt"
    stderr_log = run_dir / "stderr.txt"

    command = [
        str(binary),
        candidate.path,
        str(output_csv),
        str(samples_per_file),
        str(seed),
        mode,
    ]

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"

    t0 = time.perf_counter()
    proc = subprocess.run(command, capture_output=True, text=True, env=env, check=False)
    elapsed = time.perf_counter() - t0

    stdout_log.write_text(proc.stdout, encoding="utf-8")
    stderr_log.write_text(proc.stderr, encoding="utf-8")

    parsed = parse_validator_stdout(proc.stdout)
    actual_samples = parsed.get("actual_samples_written", "")

    metadata = {
        "command": command,
        "returncode": proc.returncode,
        "elapsed_sec": elapsed,
        "seed": seed,
        "mode": mode,
    }
    (run_dir / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    return RunResult(
        run_index=run_index,
        campaign=candidate.campaign,
        input_csv=candidate.path,
        mtime=candidate.mtime,
        size_bytes=candidate.size_bytes,
        mode=mode,
        requested_samples=samples_per_file,
        actual_samples_written=actual_samples,
        input_rows_parsed=parsed.get("input_rows_parsed", ""),
        parseable_rows=parsed.get("parseable_rows", ""),
        eligible_rows=parsed.get("eligible_rows", ""),
        baseline_failures=parsed.get("baseline_failures", ""),
        probe_failures=parsed.get("probe_failures", ""),
        max_abs_error=parsed.get("max_abs_error", ""),
        max_rel_error=parsed.get("max_rel_error", ""),
        max_scaled_error=parsed.get("max_scaled_error", ""),
        warning_rate_abs=parsed.get("warning_rate_abs", ""),
        warning_rate_rel=parsed.get("warning_rate_rel", ""),
        warning_rate_scaled=parsed.get("warning_rate_scaled", ""),
        elapsed_sec=elapsed,
        returncode=proc.returncode,
        output_csv=str(output_csv),
        stdout_log=str(stdout_log),
        stderr_log=str(stderr_log),
    )


def write_summary_csv(path: Path, rows: list[RunResult]) -> None:
    fieldnames = [
        "run_index",
        "campaign",
        "input_csv",
        "mtime",
        "size_bytes",
        "mode",
        "requested_samples",
        "actual_samples_written",
        "input_rows_parsed",
        "parseable_rows",
        "eligible_rows",
        "baseline_failures",
        "probe_failures",
        "max_abs_error",
        "max_rel_error",
        "max_scaled_error",
        "warning_rate_abs",
        "warning_rate_rel",
        "warning_rate_scaled",
        "elapsed_sec",
        "returncode",
        "output_csv",
        "stdout_log",
        "stderr_log",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def to_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def write_summary_md(path: Path, args: argparse.Namespace, selected: list[CsvCandidate], rows: list[RunResult]) -> None:
    campaigns_tested = len({row.campaign for row in rows})
    files_tested = len(rows)
    total_eligible = sum(to_int(row.eligible_rows) for row in rows)
    total_samples = sum(to_int(row.actual_samples_written) for row in rows)
    total_baseline_fail = sum(to_int(row.baseline_failures) for row in rows)
    total_probe_fail = sum(to_int(row.probe_failures) for row in rows)
    max_abs = max((to_float(row.max_abs_error) for row in rows), default=0.0)
    max_rel = max((to_float(row.max_rel_error) for row in rows), default=0.0)
    max_scaled = max((to_float(row.max_scaled_error) for row in rows), default=0.0)
    nonzero_returncodes = sum(1 for row in rows if row.returncode != 0)

    lines: list[str] = []
    lines.append("# CalcLam1ScanFromLake official-lake benchmark")
    lines.append("")
    lines.append("## Configuration")
    lines.append("")
    lines.append(f"- lake_root: `{args.lake_root}`")
    lines.append(f"- binary: `{args.binary}`")
    lines.append(f"- mode: `{args.mode}`")
    lines.append(f"- n_campaigns: `{args.n_campaigns}`")
    lines.append(f"- n_files_per_campaign: `{args.n_files_per_campaign}`")
    lines.append(f"- samples_per_file: `{args.samples_per_file}`")
    lines.append(f"- seed: `{args.seed}`")
    lines.append(f"- campaign_regex: `{args.campaign_regex or '<none>'}`")
    lines.append(f"- require_header_match: `{args.require_header_match}`")
    lines.append("")

    lines.append("## Selected official files")
    lines.append("")
    for candidate in selected:
        lines.append(f"- `{candidate.campaign}` — `{candidate.path}`")
    lines.append("")

    lines.append("## Results")
    lines.append("")
    lines.append("| campaign | returncode | eligible_rows | actual_samples_written | baseline_failures | probe_failures | max_abs_error | max_rel_error | max_scaled_error | elapsed_sec |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        lines.append(
            "| "
            f"{row.campaign} | {row.returncode} | {row.eligible_rows or ''} | {row.actual_samples_written or ''} | "
            f"{row.baseline_failures or ''} | {row.probe_failures or ''} | {row.max_abs_error or ''} | "
            f"{row.max_rel_error or ''} | {row.max_scaled_error or ''} | {row.elapsed_sec:.3f} |"
        )
    lines.append("")

    lines.append("## Aggregate summary")
    lines.append("")
    lines.append(f"- campaigns tested: **{campaigns_tested}**")
    lines.append(f"- files tested: **{files_tested}**")
    lines.append(f"- total eligible rows: **{total_eligible}**")
    lines.append(f"- total actual samples written: **{total_samples}**")
    lines.append(f"- total baseline failures: **{total_baseline_fail}**")
    lines.append(f"- total probe failures: **{total_probe_fail}**")
    lines.append(f"- runs with non-zero return code: **{nonzero_returncodes}**")
    lines.append(f"- max(max_abs_error): **{max_abs:.12g}**")
    lines.append(f"- max(max_rel_error): **{max_rel:.12g}**")
    lines.append(f"- max(max_scaled_error): **{max_scaled:.12g}**")
    lines.append("")

    if nonzero_returncodes > 0:
        interpretation = (
            "One or more validator runs returned non-zero exit codes, so this benchmark should be treated as "
            "an execution audit rather than a physics-validation result. Verify the binary path and CLI contract "
            "before drawing robustness conclusions."
        )
    elif total_baseline_fail == 0 and total_probe_fail == 0:
        interpretation = (
            "In this tested sample from multiple official campaigns, no baseline/probe failures were "
            "reported. These results support consistency for the exercised files and sample counts, "
            "without claiming full parameter-space coverage."
        )
    else:
        interpretation = (
            "The benchmark covered multiple official campaigns and observed non-zero failures in the tested "
            "sample. This indicates specific file/sample combinations requiring follow-up inspection before "
            "broader robustness claims."
        )
    lines.append("## Interpretation")
    lines.append("")
    lines.append(interpretation)
    lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    lake_root = args.lake_root.resolve()
    binary = args.binary.resolve()
    outdir = args.outdir.resolve()

    if not lake_root.exists() or not lake_root.is_dir():
        print(f"ERROR: lake root does not exist or is not a directory: {lake_root}", file=sys.stderr)
        return 2
    if not binary.exists() or not binary.is_file():
        print(f"ERROR: binary does not exist or is not a file: {binary}", file=sys.stderr)
        return 2

    outdir.mkdir(parents=True, exist_ok=True)
    runs_dir = outdir / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)

    candidates, scanned_csv = discover_csv_candidates(lake_root, args.max_files_scan)
    selected = select_candidates(
        candidates,
        n_campaigns=args.n_campaigns,
        n_files_per_campaign=args.n_files_per_campaign,
        campaign_regex=args.campaign_regex,
        require_header_match=args.require_header_match,
    )

    print(f"Discovered {len(candidates)} CSV candidates (scanned_csv={scanned_csv}, max_files_scan={args.max_files_scan}).")
    print(f"Selected {len(selected)} files across up to {args.n_campaigns} campaigns.")
    for idx, candidate in enumerate(selected):
        print(f"[{idx:03d}] {candidate.campaign} :: {candidate.path}")

    selected_payload = {
        "config": {
            "lake_root": str(lake_root),
            "binary": str(binary),
            "mode": args.mode,
            "n_campaigns": args.n_campaigns,
            "n_files_per_campaign": args.n_files_per_campaign,
            "samples_per_file": args.samples_per_file,
            "seed": args.seed,
            "campaign_regex": args.campaign_regex,
            "require_header_match": args.require_header_match,
            "max_files_scan": args.max_files_scan,
        },
        "selected": [asdict(candidate) for candidate in selected],
    }
    (outdir / "selected_inputs.json").write_text(json.dumps(selected_payload, indent=2), encoding="utf-8")

    if not selected:
        print("ERROR: no benchmark files selected.", file=sys.stderr)
        return 3

    rows: list[RunResult] = []
    for run_index, candidate in enumerate(selected):
        run_seed = derive_seed(args.seed, run_index, candidate.path)
        base_name = sanitize_name(f"{run_index:03d}_{candidate.campaign}_{Path(candidate.path).stem}")
        run_dir = runs_dir / base_name
        print(f"Running [{run_index+1}/{len(selected)}]: seed={run_seed} input={candidate.path}")
        result = run_validator(
            binary=binary,
            mode=args.mode,
            candidate=candidate,
            samples_per_file=args.samples_per_file,
            seed=run_seed,
            run_dir=run_dir,
            run_index=run_index,
        )
        rows.append(result)

    write_summary_csv(outdir / "benchmark_summary.csv", rows)
    write_summary_md(outdir / "benchmark_summary.md", args, selected, rows)
    print(f"Wrote summary CSV: {outdir / 'benchmark_summary.csv'}")
    print(f"Wrote summary MD : {outdir / 'benchmark_summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
