#!/usr/bin/env python3
"""
run_chris_existing_plots.py
============================
Minimal wrapper to make Christopher's production run reproducibly plottable.

Converts silver_all.csv to parquet and invokes existing mlpython/lake_pipeline plotters
without duplicating their logic.

Usage:
    python scripts/run_chris_existing_plots.py \
      --input data/lam1eq1_allchris_var10000_2026jun/chunked/silver_all.csv \
      --output-root data/lam1eq1_allchris_var10000_2026jun/chunked/plots_chris

Or with defaults:
    python scripts/run_chris_existing_plots.py
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

import polars as pl

logging.basicConfig(level=logging.INFO, format="[%(name)s] %(message)s")
logger = logging.getLogger("run_chris_plots")

REPO_ROOT = Path(__file__).resolve().parent.parent
MLPYTHON_DIR = REPO_ROOT / "mlpython" / "lake_pipeline"

DEFAULT_INPUT_CSV = (
    REPO_ROOT / "data" / "lam1eq1_allchris_var10000_2026jun" / "chunked" / "silver_all.csv"
)
DEFAULT_OUTPUT_ROOT = (
    REPO_ROOT / "data" / "lam1eq1_allchris_var10000_2026jun" / "chunked" / "plots_chris"
)


def csv_to_parquet(csv_path: Path, parquet_path: Path) -> None:
    """Convert CSV to parquet with zstd compression."""
    logger.info("Reading CSV: %s", csv_path)
    df = pl.read_csv(str(csv_path))
    logger.info("Writing parquet (%d rows) with zstd compression: %s", len(df), parquet_path)
    df.write_parquet(str(parquet_path), compression="zstd")


def run_plotter(script_name: str, args: list[str]) -> int:
    """Invoke a plotter script and return exit code."""
    script_path = MLPYTHON_DIR / script_name
    # Use sys.executable so subprocesses share the same Python env (polars, matplotlib, etc.)
    cmd = [sys.executable, str(script_path)] + args
    logger.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, cwd=str(REPO_ROOT))
    return result.returncode


def _write_readme(
    input_csv: Path,
    parquet_file: Path,
    output_root: Path,
    mphi_br_dir: Path,
    mphi_br_lambda1_dir: Path,
    ctau_br_dir: Path,
) -> None:
    py = sys.executable
    lines = [
        "# Christopher production run — plot report",
        "",
        "## Input",
        f"- CSV: `{input_csv}`",
        f"- Parquet (zstd): `{parquet_file}`",
        "",
        "## Scripts invoked",
        f"1. `mlpython/lake_pipeline/paper_like_mphi_vs_br_patched.py` (color by tan_beta)",
        f"2. `mlpython/lake_pipeline/paper_like_mphi_vs_br_patched.py` (color by lambda1)",
        f"3. `mlpython/lake_pipeline/ctau_vs_br_multislice_patched.py`",
        "",
        "## Output directories",
        f"- `{mphi_br_dir}` — m_phi vs BR, colored by tan_beta",
        f"- `{mphi_br_lambda1_dir}` — m_phi vs BR, colored by lambda1",
        f"- `{ctau_br_dir}` — ctau vs BR (multi-slice)",
        "",
        "## Exact commands",
        "```bash",
        f"{py} mlpython/lake_pipeline/paper_like_mphi_vs_br_patched.py \\",
        f"    --input {parquet_file} \\",
        f"    --output-dir {mphi_br_dir} \\",
        "    --color-by tan_beta",
        "",
        f"{py} mlpython/lake_pipeline/paper_like_mphi_vs_br_patched.py \\",
        f"    --input {parquet_file} \\",
        f"    --output-dir {mphi_br_lambda1_dir} \\",
        "    --color-by lambda1",
        "",
        f"{py} mlpython/lake_pipeline/ctau_vs_br_multislice_patched.py \\",
        f"    --input {parquet_file} \\",
        f"    --output-dir {ctau_br_dir} \\",
        "    --color-by both",
        "```",
        "",
        "## Smoke command (re-run from repo root)",
        "```bash",
        "bash mlpython/lake_pipeline/run_chris_existing_plots.sh \\",
        f"  {input_csv} \\",
        f"  {output_root}",
        "```",
    ]
    readme_path = output_root / "README.md"
    readme_path.write_text("\n".join(lines) + "\n")
    logger.info("README written to: %s", readme_path)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Wrapper to plot Christopher's production run using existing mlpython/lake_pipeline tools."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT_CSV,
        help=f"Input CSV (default: {DEFAULT_INPUT_CSV})",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help=f"Output root directory (default: {DEFAULT_OUTPUT_ROOT})",
    )
    args = parser.parse_args()

    input_csv = args.input
    output_root = args.output_root

    if not input_csv.exists():
        logger.error("Input CSV not found: %s", input_csv)
        return 1

    output_root.mkdir(parents=True, exist_ok=True)
    parquet_file = output_root / "silver_all.parquet"

    # Convert CSV to parquet
    csv_to_parquet(input_csv, parquet_file)

    # Set up plot subdirectories
    mphi_br_dir = output_root / "mphi_vs_br"
    mphi_br_lambda1_dir = output_root / "mphi_vs_br_lambda1"
    ctau_br_dir = output_root / "ctau_vs_br"

    mphi_br_dir.mkdir(parents=True, exist_ok=True)
    mphi_br_lambda1_dir.mkdir(parents=True, exist_ok=True)
    ctau_br_dir.mkdir(parents=True, exist_ok=True)

    # Plotter 1: m_phi vs br_gaga (colored by tan_beta, family mode)
    # silver_all is already positivity/unitarity/perturbativity-filtered, so --apply-phys-filter
    # is intentionally omitted (it triggers a column-not-found bug in family mode).
    logger.info("=== Plotting m_phi vs br_gaga (family mode, color by tan_beta) ===")
    rc = run_plotter(
        "paper_like_mphi_vs_br_patched.py",
        [
            "--input", str(parquet_file),
            "--output-dir", str(mphi_br_dir),
            "--color-by", "tan_beta",
        ],
    )
    if rc != 0:
        logger.error("paper_like_mphi_vs_br_patched.py (tan_beta mode) failed with rc=%d", rc)
        return rc

    # Plotter 2: m_phi vs br_gaga (colored by lambda1)
    logger.info("=== Plotting m_phi vs br_gaga (family mode, color by lambda1) ===")
    rc = run_plotter(
        "paper_like_mphi_vs_br_patched.py",
        [
            "--input", str(parquet_file),
            "--output-dir", str(mphi_br_lambda1_dir),
            "--color-by", "lambda1",
        ],
    )
    if rc != 0:
        logger.error("paper_like_mphi_vs_br_patched.py (lambda1 mode) failed with rc=%d", rc)
        return rc

    # Plotter 3: ctau vs br_gaga (multi-slice)
    logger.info("=== Plotting ctau vs br_gaga (multi-slice) ===")
    rc = run_plotter(
        "ctau_vs_br_multislice_patched.py",
        [
            "--input", str(parquet_file),
            "--output-dir", str(ctau_br_dir),
            "--color-by", "both",
        ],
    )
    if rc != 0:
        logger.error("ctau_vs_br_multislice_patched.py failed with rc=%d", rc)
        return rc

    # Write human-readable report
    _write_readme(input_csv, parquet_file, output_root, mphi_br_dir, mphi_br_lambda1_dir, ctau_br_dir)

    # Write summary report
    summary = {
        "input_csv": str(input_csv),
        "converted_parquet": str(parquet_file),
        "output_root": str(output_root),
        "plotters": [
            {
                "name": "paper_like_mphi_vs_br_patched.py (tan_beta)",
                "script": "mlpython/lake_pipeline/paper_like_mphi_vs_br_patched.py",
                "output_dir": str(mphi_br_dir),
                "command": f"python3 mlpython/lake_pipeline/paper_like_mphi_vs_br_patched.py --input {parquet_file} --output-dir {mphi_br_dir} --color-by tan_beta",
            },
            {
                "name": "paper_like_mphi_vs_br_patched.py (lambda1)",
                "script": "mlpython/lake_pipeline/paper_like_mphi_vs_br_patched.py",
                "output_dir": str(mphi_br_lambda1_dir),
                "command": f"python3 mlpython/lake_pipeline/paper_like_mphi_vs_br_patched.py --input {parquet_file} --output-dir {mphi_br_lambda1_dir} --color-by lambda1",
            },
            {
                "name": "ctau_vs_br_multislice_patched.py",
                "script": "mlpython/lake_pipeline/ctau_vs_br_multislice_patched.py",
                "output_dir": str(ctau_br_dir),
                "command": f"python3 mlpython/lake_pipeline/ctau_vs_br_multislice_patched.py --input {parquet_file} --output-dir {ctau_br_dir} --color-by both",
            },
        ],
    }

    summary_path = output_root / "plotting_summary.json"
    logger.info("Writing summary to: %s", summary_path)
    summary_path.write_text(json.dumps(summary, indent=2))

    logger.info("=== All plotting complete ===")
    logger.info("Output root: %s", output_root)
    logger.info("  m_phi vs br_gaga (tan_beta): %s", mphi_br_dir)
    logger.info("  m_phi vs br_gaga (lambda1): %s", mphi_br_lambda1_dir)
    logger.info("  ctau vs br_gaga: %s", ctau_br_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
