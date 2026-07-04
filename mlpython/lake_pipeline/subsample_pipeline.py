#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
subsample_pipeline.py
=====================
Create a physically filtered subsample from a parquet based on BR/ctau thresholds,
summarize relevant hyperparameters, and optionally run the full plotting pipeline
(ctau, mphi, parallel) inside a dedicated output folder.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

import polars as pl

HBAR = 6.582119569e-25  # GeV s
C_SPEED = 2.99792458e11  # mm/s

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def _stream_collect(lf: pl.LazyFrame) -> pl.DataFrame:
    try:
        return lf.collect(engine="streaming")
    except Exception:
        return lf.collect()


def resolve_column(requested: str, schema_names: list[str], fallback_aliases: list[str] | None = None) -> str:
    if requested in schema_names:
        return requested

    lower_map = {c.lower(): c for c in schema_names}
    if requested.lower() in lower_map:
        return lower_map[requested.lower()]

    if fallback_aliases:
        for alias in fallback_aliases:
            if alias.lower() in lower_map:
                resolved = lower_map[alias.lower()]
                logger.info("Resolved requested '%s' to '%s'", requested, resolved)
                return resolved

    raise ValueError(f"Could not resolve column '{requested}' in schema.")


def _flag_true_expr(col_name: str) -> pl.Expr:
    as_num = pl.col(col_name).cast(pl.Float64, strict=False)
    as_txt = (
        pl.col(col_name)
        .cast(pl.Utf8, strict=False)
        .str.strip_chars()
        .str.to_lowercase()
    )
    return (as_num >= 0.5) | as_txt.is_in(["1", "1.0", "true", "t"])


def _try_run(command: list[str], cwd: Path) -> None:
    logger.info("Running: %s", " ".join(command))
    subprocess.run(command, cwd=str(cwd), check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build and analyze a BR/ctau threshold subsample.")
    parser.add_argument("--input", type=str, required=True, help="Input parquet (preferably phys-only)")
    parser.add_argument("--output-root", type=str, default="./subsample", help="Root output directory")
    parser.add_argument("--br-col", type=str, default="br_gaga", help="BR column used for threshold")
    parser.add_argument("--br-min", type=float, required=True, help="Keep rows with BR >= br-min")
    parser.add_argument("--ctau-max", type=float, required=True, help="Keep rows with c_tau_mm <= ctau-max")
    parser.add_argument("--apply-phys-filter", action="store_true", help="Apply positivity/unitarity/perturbativity if available")
    parser.add_argument("--run-pipeline", action="store_true", help="Run ctau/mphi/parallel scripts on the generated subsample")
    parser.add_argument("--point-size", type=float, default=2.0, help="Point size forwarded to plotting scripts")
    parser.add_argument("--point-alpha", type=float, default=0.5, help="Point alpha forwarded to plotting scripts")

    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input parquet not found: {input_path}")

    work_dir = Path(__file__).resolve().parent
    output_root = Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    schema_names = pl.scan_parquet(input_path).collect_schema().names()

    col_mphi = resolve_column("m_phi", schema_names, ["mphi"])
    col_ma = resolve_column("mA", schema_names)
    col_l6 = resolve_column("lambda6", schema_names, ["lambda_6", "lam6"])
    col_tb = resolve_column("tan_beta", schema_names, ["tanbeta"])
    col_w = resolve_column("total_width", schema_names, ["total_decay_width", "w_total"])
    col_br = resolve_column(args.br_col, schema_names, ["branching_ratio", "br_gg"])

    optional_cols = []
    for candidate in ["sin_ba", "lambda7", "lam1", "br_bb", "br_tautau", "br_zz", "br_ww"]:
        if candidate in schema_names:
            optional_cols.append(candidate)

    needed_cols = [col_mphi, col_ma, col_l6, col_tb, col_w, col_br] + optional_cols

    flag_cols = [c for c in ["positivity_ok", "unitarity_ok", "perturbativity_ok"] if c in schema_names]
    if args.apply_phys_filter:
        needed_cols.extend([c for c in flag_cols if c not in needed_cols])

    lf = pl.scan_parquet(input_path).select(needed_cols)

    if args.apply_phys_filter and flag_cols:
        expr = None
        for c in flag_cols:
            e = _flag_true_expr(c)
            expr = e if expr is None else (expr & e)
        if expr is not None:
            lf = lf.filter(expr)

    lf = lf.with_columns(
        pl.when(pl.col(col_w).is_not_null() & (~pl.col(col_w).is_nan()) & (pl.col(col_w) > 0))
        .then((HBAR * C_SPEED) / pl.col(col_w).cast(pl.Float64))
        .otherwise(None)
        .alias("c_tau_mm")
    )

    lf_sub = lf.filter(
        pl.col(col_br).is_not_null()
        & pl.col("c_tau_mm").is_not_null()
        & (pl.col(col_br) >= float(args.br_min))
        & (pl.col("c_tau_mm") <= float(args.ctau_max))
    )

    total_rows = int(_stream_collect(lf.select(pl.len().alias("n")))["n"].item())
    sub_rows = int(_stream_collect(lf_sub.select(pl.len().alias("n")))["n"].item())

    logger.info("Rows before thresholds: %d", total_rows)
    logger.info("Rows after thresholds:  %d", sub_rows)

    if sub_rows == 0:
        raise RuntimeError("No rows survive thresholds; adjust --br-min / --ctau-max.")

    sub_df = _stream_collect(lf_sub)
    subsample_parquet = output_root / "subsample.parquet"
    sub_df.write_parquet(subsample_parquet)

    interesting_cols = [
        c
        for c in [col_ma, col_l6, col_tb, "sin_ba", "lambda7", col_mphi, col_br, "c_tau_mm", col_w, "br_bb", "br_tautau"]
        if c in sub_df.columns
    ]

    sub_df.select(interesting_cols).write_csv(output_root / "subsample_points.csv")

    unique_hyperparams = {
        "mA_values": sub_df.select(col_ma).drop_nulls().unique().sort(col_ma)[col_ma].to_list(),
        "lambda6_values": sub_df.select(col_l6).drop_nulls().unique().sort(col_l6)[col_l6].to_list(),
        "tan_beta_values": sub_df.select(col_tb).drop_nulls().unique().sort(col_tb)[col_tb].to_list(),
    }
    if "sin_ba" in sub_df.columns:
        unique_hyperparams["sin_ba_values"] = sub_df.select("sin_ba").drop_nulls().unique().sort("sin_ba")["sin_ba"].to_list()
    if "lambda7" in sub_df.columns:
        unique_hyperparams["lambda7_values"] = sub_df.select("lambda7").drop_nulls().unique().sort("lambda7")["lambda7"].to_list()

    summary = {
        "input_file": str(input_path),
        "subsample_file": str(subsample_parquet),
        "thresholds": {"br_col": col_br, "br_min": float(args.br_min), "ctau_max_mm": float(args.ctau_max)},
        "rows_before": total_rows,
        "rows_after": sub_rows,
        "resolved_columns": {
            "m_phi": col_mphi,
            "mA": col_ma,
            "lambda6": col_l6,
            "tan_beta": col_tb,
            "total_width": col_w,
            "br_col": col_br,
        },
        "unique_hyperparams": unique_hyperparams,
        "apply_phys_filter": bool(args.apply_phys_filter),
    }

    if args.run_pipeline:
        ctau_out = output_root / "ctau_plots"
        mphi_out = output_root / "paper_plots_mphi_br"
        parallel_out = output_root / "parallel_plots"

        _try_run(
            [
                sys.executable,
                "ctau_vs_br_multislice_patched.py",
                "--input",
                str(subsample_parquet),
                "--output-dir",
                str(ctau_out),
                "--br-col",
                col_br,
                "--color-by",
                "both",
                "--point-size",
                str(args.point_size),
                "--point-alpha",
                str(args.point_alpha),
            ],
            cwd=work_dir,
        )

        _try_run(
            [
                sys.executable,
                "paper_like_mphi_vs_br_patched.py",
                "--input",
                str(subsample_parquet),
                "--output-dir",
                str(mphi_out),
                "--mphi-col",
                col_mphi,
                "--br-col",
                col_br,
                "--point-size",
                str(args.point_size),
                "--point-alpha",
                str(args.point_alpha),
            ],
            cwd=work_dir,
        )

        _try_run(
            [
                sys.executable,
                "parallel_coordinates_roi.py",
                "--input",
                str(subsample_parquet),
                "--output-dir",
                str(parallel_out),
                "--br-threshold",
                str(args.br_min),
                "--ctau-threshold",
                str(args.ctau_max),
                "--background-alpha",
                "0.03",
                "--roi-alpha",
                "0.8",
            ],
            cwd=work_dir,
        )

        summary["pipeline_outputs"] = {
            "ctau_plots": str(ctau_out),
            "paper_plots_mphi_br": str(mphi_out),
            "parallel_plots": str(parallel_out),
        }

    summary_path = output_root / "subsample_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    logger.info("Subsample saved to: %s", subsample_parquet)
    logger.info("Summary written to: %s", summary_path)


if __name__ == "__main__":
    main()
