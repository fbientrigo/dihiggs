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

    plotters = [
        ("paper_like_mphi_vs_br_patched.py", mphi_br_dir, "tan_beta"),
        ("paper_like_mphi_vs_br_patched.py", mphi_br_lambda1_dir, "lambda1"),
        ("ctau_vs_br_multislice_patched.py", ctau_br_dir, "both"),
    ]
    for script, out_dir, color_by in plotters:
        logger.info("=== Plotting %s (color by %s) ===", script, color_by)
        rc = run_plotter(
            script,
            ["--input", str(parquet_file), "--output-dir", str(out_dir), "--color-by", color_by],
        )
        if rc != 0:
            logger.error("%s failed with rc=%d", script, rc)
            return rc

    logger.info("=== All plotting complete ===")
    logger.info("Output root: %s", output_root)
    logger.info("  m_phi vs br_gaga (tan_beta): %s", mphi_br_dir)
    logger.info("  m_phi vs br_gaga (lambda1): %s", mphi_br_lambda1_dir)
    logger.info("  ctau vs br_gaga: %s", ctau_br_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
