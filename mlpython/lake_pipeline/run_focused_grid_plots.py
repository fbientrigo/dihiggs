#!/usr/bin/env python3
"""Run existing plotting scripts on focused-grid slices only.

Reads ``manifest.csv`` and, for each ``slice.parquet``, builds plotting
commands for the existing scripts.  It NEVER targets the full
``silver_all.parquet`` and refuses any input that is not a focused
``slice.parquet`` (hard safety guard against the unbounded global family mode).

Plotting logic is not duplicated here; this only orchestrates subprocess calls.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PRIMARY_SCRIPT = "paper_like_three_br_ctau_lambda1.py"
MPHI_SCRIPT = "paper_like_mphi_vs_br_patched.py"
CTAU_SCRIPT = "ctau_vs_br_multislice_patched.py"

PLOT_MANIFEST_FIELDS = [
    "slice_path", "plot_output_dir", "script", "command", "returncode",
    "status", "stdout_path", "stderr_path", "created_utc",
]


def fmt(value: float) -> str:
    return f"{value:g}"


def assert_focused_slice(slice_path: Path) -> None:
    """Refuse anything that is not a focused slice.parquet."""
    name = slice_path.name
    if name != "slice.parquet" or "silver_all" in str(slice_path):
        raise SystemExit(
            f"[plots] REFUSING non-focused input: {slice_path}. "
            "Only focused 'slice.parquet' files are allowed."
        )


def plot_subdir(output_root: Path, row: dict, script_stem: str) -> Path:
    return (
        output_root
        / f"mA={fmt(float(row['mA_target']))}"
        / f"lambda6={fmt(float(row['lambda6']))}"
        / f"tan_beta={fmt(float(row['tan_beta']))}"
        / script_stem
    )


def build_commands(row: dict, output_root: Path, args) -> list[tuple[str, Path, list[str]]]:
    slice_path = Path(row["slice_path"])
    assert_focused_slice(slice_path)
    cmds: list[tuple[str, Path, list[str]]] = []

    primary_dir = plot_subdir(output_root, row, Path(PRIMARY_SCRIPT).stem)
    primary_cmd = [
        sys.executable, str(SCRIPT_DIR / PRIMARY_SCRIPT),
        "--input", str(slice_path),
        "--output-dir", str(primary_dir),
        "--mA", fmt(float(row["mA_target"])),
        "--tan-beta", fmt(float(row["tan_beta"])),
        "--lambda6", fmt(float(row["lambda6"])),
        "--lambda7", fmt(float(row["lambda7"])),
    ]
    if args.sin_ba is not None:
        primary_cmd += ["--sin-ba", fmt(float(args.sin_ba))]
    cmds.append((PRIMARY_SCRIPT, primary_dir, primary_cmd))

    if args.with_mphi:
        mphi_dir = plot_subdir(output_root, row, Path(MPHI_SCRIPT).stem)
        cmds.append((MPHI_SCRIPT, mphi_dir, [
            sys.executable, str(SCRIPT_DIR / MPHI_SCRIPT),
            "--input", str(slice_path),
            "--output-dir", str(mphi_dir),
            "--mA", fmt(float(row["mA_target"])),
            "--tan-beta", fmt(float(row["tan_beta"])),
            "--lambda6", fmt(float(row["lambda6"])),
            "--lambda7", fmt(float(row["lambda7"])),
        ]))

    if args.with_ctau:
        ctau_dir = plot_subdir(output_root, row, Path(CTAU_SCRIPT).stem)
        cmds.append((CTAU_SCRIPT, ctau_dir, [
            sys.executable, str(SCRIPT_DIR / CTAU_SCRIPT),
            "--input", str(slice_path),
            "--output-dir", str(ctau_dir),
        ]))
    return cmds


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", required=True, help="Path to slices manifest.csv")
    p.add_argument("--output-root", required=True,
                   help="Root directory for plot outputs")
    p.add_argument("--dry-run", action="store_true",
                   help="Print commands without executing")
    p.add_argument("--limit", type=int, default=None,
                   help="Process at most N slices (smoke testing)")
    p.add_argument("--sin-ba", type=float, default=None,
                   help="Optional sin_ba cut for the primary script "
                        "(default: let the script use its own default of 1.0)")
    p.add_argument("--with-mphi", action="store_true",
                   help="Also run paper_like_mphi_vs_br_patched.py")
    p.add_argument("--with-ctau", action="store_true",
                   help="Also run ctau_vs_br_multislice_patched.py")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    manifest_path = Path(args.manifest)
    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    with manifest_path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if args.limit is not None:
        rows = rows[: args.limit]

    created = datetime.now(timezone.utc).isoformat()
    records: list[dict] = []
    for row in rows:
        for script_name, plot_dir, cmd in build_commands(row, output_root, args):
            command_str = " ".join(cmd)
            record = {
                "slice_path": row["slice_path"],
                "plot_output_dir": str(plot_dir),
                "script": script_name,
                "command": command_str,
                "returncode": None,
                "status": "dry_run",
                "stdout_path": None,
                "stderr_path": None,
                "created_utc": created,
            }
            if args.dry_run:
                print(f"[dry-run] {command_str}")
            else:
                plot_dir.mkdir(parents=True, exist_ok=True)
                stdout_path = plot_dir / "plot_stdout.txt"
                stderr_path = plot_dir / "plot_stderr.txt"
                proc = subprocess.run(cmd, capture_output=True, text=True)
                stdout_path.write_text(proc.stdout)
                stderr_path.write_text(proc.stderr)
                record["returncode"] = proc.returncode
                record["status"] = "ok" if proc.returncode == 0 else "failed"
                record["stdout_path"] = str(stdout_path)
                record["stderr_path"] = str(stderr_path)
                print(f"[plots] {record['status']} ({proc.returncode}): {command_str}")
            records.append(record)

    records.sort(key=lambda r: (r["slice_path"], r["script"]))
    plot_csv = output_root / "plot_manifest.csv"
    plot_json = output_root / "plot_manifest.json"
    with plot_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=PLOT_MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerows(records)
    with plot_json.open("w") as f:
        json.dump(records, f, indent=2, sort_keys=True, default=str)
        f.write("\n")

    print(f"[plots] wrote {plot_csv}")
    print(f"[plots] wrote {plot_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
