#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ctau_vs_br_multislice_final.py
==============================
Memory-safer replacement for the original `ctau_vs_br_multislice.py`.

What changed
------------
- Keeps the subset behavior and figure style that already worked.
- Avoids reading the full parquet into pandas at once.
- Reads only the columns needed for the current task.
- Computes c_tau lazily with Polars instead of copying a huge pandas DataFrame.
- Skips redundant internal physical filtering by default when the input already
  looks like a phys-only parquet.
- Materializes only one plotting slice at a time into pandas.
- Fixes the wrapper/call mismatch when `--max-combinations` is active.

Recommended usage
-----------------
Subset (already filtered by your workflow):
    python ctau_vs_br_multislice_final.py --input temp_subspace.parquet --color-by both

Full phys-only parquet:
    python ctau_vs_br_multislice_final.py \
        --input /home/fabi/cern_db/dihiggs_consolidated/dihiggs_lake_phys_only.parquet \
        --color-by both

Raw parquet (apply the 3 *_ok filters inside the script):
    python ctau_vs_br_multislice_final.py \
        --input /home/fabi/cern_db/dihiggs_consolidated/dihiggs_lake.parquet \
        --color-by both \
        --apply-phys-filter
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl

# ------------------------------------------------------------------------------
# Set up plotting style (paper-friendly)
# ------------------------------------------------------------------------------
plt.style.use("default")
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 18,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 12,
        "legend.frameon": True,
        "legend.edgecolor": "black",
        "legend.fancybox": False,
        "figure.titlesize": 20,
        "figure.figsize": (8, 6),
        "lines.linewidth": 2.0,
        "axes.linewidth": 1.5,
        "xtick.major.width": 1.5,
        "ytick.major.width": 1.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    }
)

# ------------------------------------------------------------------------------
# Config & Setup
# ------------------------------------------------------------------------------
DEFAULT_PARQUET = Path(__file__).parent / "temp_subspace.parquet"
OUTPUT_DIR = Path(__file__).parent / "ctau_plots"

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

HBAR = 6.582119569e-25  # GeV s
C_SPEED = 2.99792458e11  # mm/s


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

    raise ValueError(f"Could not resolve column '{requested}' in the schema. Available: {schema_names}")


def _stream_collect(lf: pl.LazyFrame) -> pl.DataFrame:
    try:
        return lf.collect(engine="streaming")
    except TypeError:
        return lf.collect(streaming=True)
    except Exception:
        return lf.collect()


def _flag_true_expr(col_name: str) -> pl.Expr:
    as_num = pl.col(col_name).cast(pl.Float64, strict=False)
    as_txt = (
        pl.col(col_name)
        .cast(pl.Utf8, strict=False)
        .str.strip_chars()
        .str.to_lowercase()
    )
    return (as_num >= 0.5) | as_txt.is_in(["1", "1.0", "true", "t"])


def generate_family_plot(
    df: pd.DataFrame,
    fixed_params: dict,
    varying_col: str,
    varying_vals: list,
    x_col: str,
    y_col: str,
    output_stem: str,
    point_size: float,
    point_alpha: float,
) -> bool:
    fig, ax = plt.subplots()
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(varying_vals)))
    groups_plotted = 0

    for idx, v_val in enumerate(sorted(varying_vals)):
        df_curve = df[df[varying_col] == v_val]
        if len(df_curve) == 0:
            logger.warning("  -> Skipping %s=%s: empty group.", varying_col, v_val)
            continue

        if varying_col == "tan_beta":
            val_label = f"$\\tan\\beta = {v_val:g}$"
        elif varying_col == "lambda6":
            val_label = f"$\\lambda_6 = {v_val:g}$"
        else:
            val_label = f"${varying_col} = {v_val:g}$"

        raw_x = df_curve[x_col].to_numpy()
        raw_y = df_curve[y_col].to_numpy()
        valid_scatter = (~np.isnan(raw_x)) & (~np.isnan(raw_y))
        if not valid_scatter.any():
            logger.warning("  -> Skipping %s=%s: no valid points after NaN filtering.", varying_col, v_val)
            continue

        ax.scatter(
            raw_x[valid_scatter],
            raw_y[valid_scatter],
            color=colors[idx],
            alpha=point_alpha,
            s=point_size,
            edgecolors="none",
            label=val_label,
            zorder=2,
        )
        groups_plotted += 1

    if groups_plotted == 0:
        logger.warning("No valid point groups plotted for this slice combination.")
        plt.close(fig)
        return False

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$c\tau$ [mm]")
    br_display = r"BR($\phi \to \gamma\gamma$)" if "gaga" in y_col.lower() else f"BR ({y_col})"
    ax.set_ylabel(br_display)

    title_parts = []
    if "mA" in fixed_params:
        title_parts.append(f"$m_A = {fixed_params['mA']:g}$")
    if "lambda6" in fixed_params:
        title_parts.append(f"$\\lambda_6 = {fixed_params['lambda6']:g}$")
    if "tan_beta" in fixed_params:
        title_parts.append(f"$\\tan\\beta = {fixed_params['tan_beta']:g}$")

    ax.set_title(r"$\quad$".join(title_parts))
    ax.grid(True, linestyle="--", alpha=0.4, which="both")
    ax.legend(loc="best")
    fig.tight_layout()

    png_path = OUTPUT_DIR / f"{output_stem}.png"
    pdf_path = OUTPUT_DIR / f"{output_stem}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return True


def run_analysis_lazy(
    lf: pl.LazyFrame,
    br_col: str,
    color_mode: str,
    summary: dict,
    point_size: float,
    point_alpha: float,
) -> None:
    logger.info("\n=========================================")
    logger.info("Running Variant: Color by %s", color_mode)
    logger.info("=========================================")

    if color_mode == "lambda6":
        fixed_cols = ["mA", "tan_beta"]
    elif color_mode == "tan_beta":
        fixed_cols = ["mA", "lambda6"]
    else:
        raise ValueError("Invalid color_mode")

    unique_colors = _stream_collect(lf.select(color_mode).drop_nulls().unique().sort(color_mode))[color_mode].to_list()
    logger.info("Unique values for %s: %d", color_mode, len(unique_colors))
    summary[f"unique_{color_mode}_count"] = int(len(unique_colors))

    if len(unique_colors) < 2:
        logger.warning("[!] The dataset only has %d unique value for '%s' (%s).", len(unique_colors), color_mode, unique_colors)
        logger.warning("    Therefore, plots differentiating by '%s' are not highly informative.", color_mode)
        logger.warning("    They will be plotted as a single curve. To get a family of curves, consider using the full dataset.")
    elif len(unique_colors) > 15:
        logger.warning("    Too many %s curves (%d). Sampling 10 evenly spaced values for plotting readability.", color_mode, len(unique_colors))
        sorted_vals = sorted(unique_colors)
        indices = np.linspace(0, len(sorted_vals) - 1, 10, dtype=int)
        unique_colors = [sorted_vals[i] for i in indices]

    clean_lf = lf.filter(
        pl.col("c_tau_mm").is_not_null()
        & pl.col(br_col).is_not_null()
        & pl.col(color_mode).is_not_null()
        & pl.all_horizontal([pl.col(c).is_not_null() for c in fixed_cols])
    )

    unique_params = _stream_collect(clean_lf.select(fixed_cols).unique().sort(fixed_cols)).to_pandas()
    logger.info("Found %d unique combinations of fixed parameters %s.", len(unique_params), fixed_cols)

    var_summary = {
        "mode": color_mode,
        "total_combinations_investigated": int(len(unique_params)),
        "plots_generated": [],
        "discarded_combinations": [],
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for _, row in unique_params.iterrows():
        fixed_dict = {col: float(row[col]) for col in fixed_cols}
        logger.info("  -> Processing slice: %s", fixed_dict)

        expr = None
        for col, val in fixed_dict.items():
            e = pl.col(col) == val
            expr = e if expr is None else (expr & e)

        df_slice = _stream_collect(clean_lf.filter(expr).select(["c_tau_mm", br_col, color_mode])).to_pandas()

        if len(df_slice) == 0:
            var_summary["discarded_combinations"].append({**fixed_dict, "reason": "empty slice"})
            continue

        actual_varying = [c for c in unique_colors if c in df_slice[color_mode].values]
        stem_parts = [f"{k}{v:g}".replace(".", "p") for k, v in fixed_dict.items()]
        output_stem = f"ctau_vs_{br_col}_color_{color_mode}_" + "_".join(stem_parts)

        success = generate_family_plot(
            df_slice,
            fixed_dict,
            color_mode,
            actual_varying,
            "c_tau_mm",
            br_col,
            output_stem,
            point_size,
            point_alpha,
        )
        if success:
            var_summary["plots_generated"].append(
                {
                    "fixed_params": fixed_dict,
                    "curves_plotted": actual_varying,
                    "points_in_slice": int(len(df_slice)),
                    "file_stem": output_stem,
                }
            )
        else:
            var_summary["discarded_combinations"].append({**fixed_dict, "reason": "no valid point groups"})

    summary["variants"].append(var_summary)


def main() -> None:
    global OUTPUT_DIR
    parser = argparse.ArgumentParser(description="Generate ctau vs BR plots sliced by lambda6 or tan_beta.")
    parser.add_argument("--input", type=str, default=str(DEFAULT_PARQUET), help="Path to the parquet file (subset or full lake).")
    parser.add_argument("--output-dir", type=str, default=str(OUTPUT_DIR), help="Directory where plots and summary are written")
    parser.add_argument("--br-col", type=str, default="br_gaga", help="Name of the Branching Ratio column")
    parser.add_argument("--color-by", type=str, choices=["lambda6", "tan_beta", "both"], default="both", help="Which parameter to vary inside a single plot.")
    parser.add_argument(
        "--apply-phys-filter",
        action="store_true",
        help="Apply positivity/unitarity/perturbativity filtering inside this script. Leave OFF for phys-only parquet or subset parquet.",
    )
    parser.add_argument(
        "--max-combinations",
        type=int,
        default=None,
        help="Optional safety limit for number of fixed-parameter combinations processed in each variant.",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=2.0,
        help="Scatter marker size in points^2.",
    )
    parser.add_argument(
        "--point-alpha",
        type=float,
        default=0.5,
        help="Scatter alpha in [0, 1].",
    )
    args = parser.parse_args()
    OUTPUT_DIR = Path(args.output_dir)

    input_path = Path(args.input)
    if not input_path.exists():
        logger.error("Input file not found: %s", input_path)
        return

    logger.info("Loading data lazily from %s", input_path)
    schema = pl.scan_parquet(input_path).collect_schema()
    schema_names = schema.names()

    br_col = resolve_column(args.br_col, schema_names, ["branching_ratio", "br_gg"])
    width_col = resolve_column("total_width", schema_names, ["total_decay_width", "w_total"])
    logger.info("Resolved target columns: BR='%s', Width='%s'", br_col, width_col)

    needed_cols = [br_col, width_col, "mA", "lambda6", "tan_beta"]
    flag_cols = ["positivity_ok", "unitarity_ok", "perturbativity_ok"]
    present_flags = [c for c in flag_cols if c in schema_names]
    if args.apply_phys_filter:
        needed_cols.extend([c for c in present_flags if c not in needed_cols])

    lf = pl.scan_parquet(input_path).select(needed_cols)

    if args.apply_phys_filter:
        if present_flags:
            logger.info("Applying mandatory theoretical filters: %s", present_flags)
            expr = None
            for c in present_flags:
                e = _flag_true_expr(c)
                expr = e if expr is None else (expr & e)
            lf = lf.filter(expr)
            rows_after = _stream_collect(lf.select(pl.len().alias("n")))["n"].item()
            logger.info("Rows remaining after filtering: %d", rows_after)
        else:
            logger.info("No theoretical filters found in schema. Proceeding without filtering.")
    else:
        if "phys_only" in input_path.name:
            logger.info("Input looks like a phys-only parquet; skipping redundant internal physical filtering.")
        else:
            logger.info("Skipping internal physical filter. Use --apply-phys-filter only for raw parquet inputs.")

    logger.info("Calculating ctau [mm] using convention: c_tau = hbar * c / %s", width_col)
    lf = lf.with_columns(
        pl.when(pl.col(width_col).is_not_null() & (~pl.col(width_col).is_nan()) & (pl.col(width_col) > 0))
        .then((HBAR * C_SPEED) / pl.col(width_col).cast(pl.Float64))
        .otherwise(None)
        .alias("c_tau_mm")
    )

    counts_df = _stream_collect(
        lf.select(
            pl.len().alias("rows"),
            pl.col("c_tau_mm").is_not_null().sum().alias("valid_ctau"),
        )
    )
    total_rows = int(counts_df["rows"].item())
    valid_ctau = int(counts_df["valid_ctau"].item())
    logger.info("Successfully computed ctau for %d rows. (%d rows were invalid/NaN)", valid_ctau, total_rows - valid_ctau)

    if valid_ctau == 0:
        logger.error("No valid ctau values computed. Exiting.")
        return

    summary: dict[str, object] = {
        "input_file": str(input_path),
        "apply_phys_filter": bool(args.apply_phys_filter),
        "total_rows_read": total_rows,
        "valid_ctau_rows": valid_ctau,
        "variants": [],
    }

    if args.max_combinations is not None:
        logger.info("Safety cap active: max combinations per variant = %d", args.max_combinations)

        def capped_run_analysis_lazy(
            lf: pl.LazyFrame,
            br_col: str,
            color_mode: str,
            summary: dict,
            point_size: float,
            point_alpha: float,
        ) -> None:
            logger.info("\n=========================================")
            logger.info("Running Variant: Color by %s", color_mode)
            logger.info("=========================================")

            if color_mode == "lambda6":
                fixed_cols = ["mA", "tan_beta"]
            else:
                fixed_cols = ["mA", "lambda6"]

            unique_colors = _stream_collect(lf.select(color_mode).drop_nulls().unique().sort(color_mode))[color_mode].to_list()
            logger.info("Unique values for %s: %d", color_mode, len(unique_colors))
            summary[f"unique_{color_mode}_count"] = int(len(unique_colors))
            if len(unique_colors) > 15:
                sorted_vals = sorted(unique_colors)
                indices = np.linspace(0, len(sorted_vals) - 1, 10, dtype=int)
                unique_colors = [sorted_vals[i] for i in indices]

            clean_lf = lf.filter(
                pl.col("c_tau_mm").is_not_null()
                & pl.col(br_col).is_not_null()
                & pl.col(color_mode).is_not_null()
                & pl.all_horizontal([pl.col(c).is_not_null() for c in fixed_cols])
            )
            unique_params = _stream_collect(clean_lf.select(fixed_cols).unique().sort(fixed_cols)).to_pandas().head(args.max_combinations)
            logger.info("Found %d unique combinations of fixed parameters %s (capped).", len(unique_params), fixed_cols)

            var_summary = {
                "mode": color_mode,
                "total_combinations_investigated": int(len(unique_params)),
                "plots_generated": [],
                "discarded_combinations": [],
            }
            OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

            for _, row in unique_params.iterrows():
                fixed_dict = {col: float(row[col]) for col in fixed_cols}
                logger.info("  -> Processing slice: %s", fixed_dict)
                expr = None
                for col, val in fixed_dict.items():
                    e = pl.col(col) == val
                    expr = e if expr is None else (expr & e)
                df_slice = _stream_collect(clean_lf.filter(expr).select(["c_tau_mm", br_col, color_mode])).to_pandas()
                if len(df_slice) == 0:
                    var_summary["discarded_combinations"].append({**fixed_dict, "reason": "empty slice"})
                    continue
                actual_varying = [c for c in unique_colors if c in df_slice[color_mode].values]
                stem_parts = [f"{k}{v:g}".replace(".", "p") for k, v in fixed_dict.items()]
                output_stem = f"ctau_vs_{br_col}_color_{color_mode}_" + "_".join(stem_parts)
                success = generate_family_plot(
                    df_slice,
                    fixed_dict,
                    color_mode,
                    actual_varying,
                    "c_tau_mm",
                    br_col,
                    output_stem,
                    point_size,
                    point_alpha,
                )
                if success:
                    var_summary["plots_generated"].append({
                        "fixed_params": fixed_dict,
                        "curves_plotted": actual_varying,
                        "points_in_slice": int(len(df_slice)),
                        "file_stem": output_stem,
                    })
                else:
                        var_summary["discarded_combinations"].append({**fixed_dict, "reason": "no valid point groups"})
            summary["variants"].append(var_summary)

        runner = capped_run_analysis_lazy
    else:
        runner = run_analysis_lazy

    if args.color_by in ["lambda6", "both"]:
        runner(lf, br_col, "lambda6", summary, args.point_size, args.point_alpha)
    if args.color_by in ["tan_beta", "both"]:
        runner(lf, br_col, "tan_beta", summary, args.point_size, args.point_alpha)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_path = OUTPUT_DIR / "ctau_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    logger.info("\nAnalysis complete! Summary saved to %s", summary_path)
    logger.info("Figures saved in %s/", OUTPUT_DIR)


if __name__ == "__main__":
    main()
