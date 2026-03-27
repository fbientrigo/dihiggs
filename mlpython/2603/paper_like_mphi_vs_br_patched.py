#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
paper_like_mphi_vs_br_patched.py
================================
Memory-safer replacement for the original `paper_like_mphi_vs_br.py`.

What changed
------------
- Keeps the subset behavior and plotting style that already worked.
- Avoids loading the full parquet into pandas in one shot.
- Uses Polars lazy scanning to read only the columns that are needed.
- Builds the list of (mA, lambda6) slices lazily.
- Materializes only one slice at a time into pandas for plotting.
- Optional physical filtering can be enabled explicitly if the input parquet is raw.

Recommended usage
-----------------
Subset (already filtered by your workflow):
    python paper_like_mphi_vs_br_patched.py --input temp_subspace.parquet

Full phys-only parquet:
    python paper_like_mphi_vs_br_patched.py \
        --input /home/fabi/cern_db/dihiggs_consolidated/dihiggs_lake_phys_only.parquet

Raw parquet (only if you explicitly want this script to apply the 3 *_ok filters):
    python paper_like_mphi_vs_br_patched.py \
        --input /home/fabi/cern_db/dihiggs_consolidated/dihiggs_lake.parquet \
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
OUTPUT_DIR = Path(__file__).parent / "paper_plots_mphi_br"

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def resolve_column(requested: str, schema_names: list[str], fallback_aliases: list[str] | None = None) -> str:
    """Find the exact column name in the schema, with optional fallback heuristics."""
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
    """Collect with streaming engine when available; fall back safely otherwise."""
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


def process_slice(df_slice: pd.DataFrame, mphi_col: str, br_col: str) -> pd.DataFrame:
    """Aggregate duplicate m_phi points to mean/std/count for cleaner curves."""
    return (
        df_slice.groupby(mphi_col)
        .agg(
            mean_br=(br_col, "mean"),
            std_br=(br_col, "std"),
            count_br=(br_col, "count"),
        )
        .reset_index()
        .sort_values(by=mphi_col)
    )


def generate_family_plot(
    df: pd.DataFrame,
    ma_val: float,
    l6_val: float,
    tb_vals: list[float],
    mphi_col: str,
    br_col: str,
    output_stem: str,
) -> bool:
    """Plot m_phi vs BR for a fixed (mA, lambda6) and multiple tan_beta values."""
    fig, ax = plt.subplots()
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(tb_vals)))
    curves_plotted = 0

    for idx, tb in enumerate(sorted(tb_vals)):
        df_curve = df[df["tan_beta"] == tb]
        if len(df_curve) < 2:
            logger.warning("  -> Skipping tan_beta=%s: only %d points.", tb, len(df_curve))
            continue

        df_agg = process_slice(df_curve, mphi_col, br_col)
        if len(df_agg) < 2:
            logger.warning("  -> Skipping tan_beta=%s: not enough aggregated points.", tb)
            continue

        m_phi = df_agg[mphi_col].to_numpy()
        br_mean = df_agg["mean_br"].to_numpy()
        br_std = df_agg["std_br"].to_numpy()
        tb_label = f"$\\tan\\beta = {tb:g}$"

        raw_mphi = df_curve[mphi_col].to_numpy()
        raw_br = df_curve[br_col].to_numpy()
        valid = (~np.isnan(raw_mphi)) & (~np.isnan(raw_br))
        ax.scatter(raw_mphi[valid], raw_br[valid], color=colors[idx], alpha=0.15, s=15, edgecolors="none", zorder=1)
        ax.plot(m_phi, br_mean, color=colors[idx], label=tb_label, zorder=2)

        valid_std = ~np.isnan(br_std)
        if valid_std.any():
            ax.fill_between(
                m_phi[valid_std],
                br_mean[valid_std] - br_std[valid_std],
                br_mean[valid_std] + br_std[valid_std],
                color=colors[idx],
                alpha=0.2,
                zorder=1,
            )
        curves_plotted += 1

    if curves_plotted == 0:
        logger.warning("No valid curves plotted for this slice combination.")
        plt.close(fig)
        return False

    ax.set_yscale("log")
    ax.set_xlabel(r"$m_\phi$ [GeV]")
    br_display = r"BR($\phi \to \gamma\gamma$)" if "gaga" in br_col.lower() else f"BR ({br_col})"
    ax.set_ylabel(br_display)
    ax.set_title(f"$m_A = {ma_val:g}$ GeV,  $\\lambda_6 = {l6_val:g}$")
    ax.grid(True, linestyle="--", alpha=0.4, which="both")
    ax.legend()
    fig.tight_layout()

    png_path = OUTPUT_DIR / f"{output_stem}.png"
    pdf_path = OUTPUT_DIR / f"{output_stem}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate paper-quality mphi vs BR plots.")
    parser.add_argument("--input", type=str, default=str(DEFAULT_PARQUET), help="Path to the parquet file (subset or full lake).")
    parser.add_argument("--mphi-col", type=str, default="m_phi", help="Name of the m_phi column")
    parser.add_argument("--br-col", type=str, default="br_gaga", help="Name of the Branching Ratio column")
    parser.add_argument(
        "--apply-phys-filter",
        action="store_true",
        help="Apply positivity/unitarity/perturbativity filters inside this script. Leave OFF for phys-only parquet or subset parquet.",
    )
    parser.add_argument(
        "--max-slices",
        type=int,
        default=None,
        help="Optional safety limit for number of (mA, lambda6) slices to process.",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        logger.error("Input file not found: %s", input_path)
        logger.info("If you meant to use the full lake, point --input to it.")
        return

    logger.info("Loading data lazily from %s", input_path)
    schema = pl.scan_parquet(input_path).collect_schema()
    schema_names = schema.names()
    logger.info("Available columns: %d", len(schema_names))

    mphi_col = resolve_column(args.mphi_col, schema_names, ["mphi"])
    br_col = resolve_column(args.br_col, schema_names, ["branching_ratio", "br_gg"])

    for req in ["mA", "lambda6", "tan_beta"]:
        if req not in schema_names:
            logger.error("Required parameter column '%s' not found in schema!", req)
            return

    logger.info("Resolved target columns: X='%s', Y='%s'", mphi_col, br_col)

    needed_cols = [mphi_col, br_col, "mA", "lambda6", "tan_beta"]
    present_flags = [c for c in ["positivity_ok", "unitarity_ok", "perturbativity_ok"] if c in schema_names]
    if args.apply_phys_filter:
        needed_cols.extend([c for c in present_flags if c not in needed_cols])

    lf = pl.scan_parquet(input_path).select(needed_cols)

    if args.apply_phys_filter:
        if present_flags:
            logger.info("Applying internal physical filter: %s", present_flags)
            expr = None
            for c in present_flags:
                e = _flag_true_expr(c)
                expr = e if expr is None else (expr & e)
            lf = lf.filter(expr)
        else:
            logger.warning("--apply-phys-filter requested, but no *_ok columns were found.")
    else:
        logger.info("Skipping internal physical filter. Assumes input is already the desired dataset.")

    lf = lf.filter(pl.col(mphi_col).is_not_null() & pl.col(br_col).is_not_null())

    unique_params_lf = lf.select(["mA", "lambda6"]).unique().sort(["mA", "lambda6"])
    unique_params = _stream_collect(unique_params_lf).to_pandas()

    if args.max_slices is not None:
        unique_params = unique_params.head(args.max_slices)

    logger.info("Found %d unique (mA, lambda6) combinations.", len(unique_params))

    summary: dict[str, object] = {
        "input_file": str(input_path),
        "apply_phys_filter": bool(args.apply_phys_filter),
        "total_combinations_investigated": int(len(unique_params)),
        "plots_generated": [],
        "discarded_combinations": [],
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for _, row in unique_params.iterrows():
        ma_val = float(row["mA"])
        l6_val = float(row["lambda6"])
        logger.info("\nProcessing slice: mA=%s, lambda6=%s", ma_val, l6_val)

        df_slice = _stream_collect(
            lf.filter((pl.col("mA") == ma_val) & (pl.col("lambda6") == l6_val)).select([mphi_col, br_col, "tan_beta"])
        ).to_pandas()

        if len(df_slice) == 0:
            logger.warning("  -> Slice is empty after collection. Discarding.")
            summary["discarded_combinations"].append({"mA": ma_val, "lambda6": l6_val, "reason": "empty"})
            continue

        unique_tb = df_slice["tan_beta"].dropna().unique().tolist()
        logger.info("  -> Found %d unique tan_beta values: %s", len(unique_tb), unique_tb)
        logger.info("  -> Total points in this slice: %d", len(df_slice))

        if len(unique_tb) > 15:
            logger.warning("  -> Too many tan_beta curves (>15). Taking 10 linearly spaced values to avoid clutter.")
            sorted_tb = sorted(unique_tb)
            indices = np.linspace(0, len(sorted_tb) - 1, 10, dtype=int)
            unique_tb = [sorted_tb[i] for i in indices]

        output_stem = f"mphi_vs_{args.br_col}_mA{ma_val:g}_l6_{l6_val:g}".replace(".", "p")
        success = generate_family_plot(df_slice, ma_val, l6_val, unique_tb, mphi_col, br_col, output_stem)

        if success:
            summary["plots_generated"].append(
                {
                    "mA": ma_val,
                    "lambda6": l6_val,
                    "tan_beta_curves": unique_tb,
                    "points_in_slice": int(len(df_slice)),
                    "file_stem": output_stem,
                }
            )
        else:
            summary["discarded_combinations"].append(
                {"mA": ma_val, "lambda6": l6_val, "reason": "no valid curves generated (<2 pts per tb)"}
            )

    summary_path = OUTPUT_DIR / "summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    logger.info("\nAnalysis complete! Summary saved to %s", summary_path)
    logger.info("Figures saved in %s/", OUTPUT_DIR)


if __name__ == "__main__":
    main()
