#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
paper_like_mphi_vs_br_patched.py
================================
Dual-mode plotting utility:
1) Family mode (default): m_phi vs BR for each (mA, lambda6), colored by tan_beta.
2) Fixed-cut mode: single plot after applying fixed cuts on
   (sin_ba, tan_beta, mA, lambda6, lambda7).

This keeps compatibility with both:
- run_pipeline.sh (family mode)
- previous fixed-cut invocations.
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


plt.style.use("default")
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 18,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 11,
        "legend.frameon": True,
        "legend.edgecolor": "black",
        "legend.fancybox": False,
        "figure.titlesize": 20,
        "figure.figsize": (8, 6),
        "axes.linewidth": 1.5,
        "xtick.major.width": 1.5,
        "ytick.major.width": 1.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    }
)

DEFAULT_PARQUET = Path(__file__).parent / "temp_subspace.parquet"
OUTPUT_DIR = Path(__file__).parent / "paper_plots_mphi_br"

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def _stream_collect(lf: pl.LazyFrame) -> pl.DataFrame:
    try:
        return lf.collect(engine="streaming")
    except TypeError:
        return lf.collect(streaming=True)
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

    raise ValueError(f"Could not resolve column '{requested}' in the schema.")


def _flag_true_expr(col_name: str) -> pl.Expr:
    as_num = pl.col(col_name).cast(pl.Float64, strict=False)
    as_txt = (
        pl.col(col_name)
        .cast(pl.Utf8, strict=False)
        .str.strip_chars()
        .str.to_lowercase()
    )
    return (as_num >= 0.5) | as_txt.is_in(["1", "1.0", "true", "t"])


def _plot_scatter_family(
    df: pd.DataFrame,
    ma_val: float,
    l6_val: float,
    tb_vals: list[float],
    mphi_col: str,
    br_col: str,
    point_size: float,
    point_alpha: float,
    output_stem: str,
) -> bool:
    fig, ax = plt.subplots()
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(tb_vals)))
    groups_plotted = 0

    for idx, tb in enumerate(sorted(tb_vals)):
        dfi = df[df["tan_beta"] == tb]
        if len(dfi) == 0:
            continue
        x = dfi[mphi_col].to_numpy()
        y = dfi[br_col].to_numpy()
        valid = np.isfinite(x) & np.isfinite(y)
        if not valid.any():
            continue
        ax.scatter(
            x[valid],
            y[valid],
            color=colors[idx],
            alpha=point_alpha,
            s=point_size,
            edgecolors="none",
            label=f"$\\tan\\beta={tb:g}$",
            zorder=2,
        )
        groups_plotted += 1

    if groups_plotted == 0:
        plt.close(fig)
        return False

    ax.set_yscale("log")
    ax.set_xlabel(r"$m_\phi$ [GeV]")
    br_display = r"BR($\phi \to \gamma\gamma$)" if "gaga" in br_col.lower() else f"BR ({br_col})"
    ax.set_ylabel(br_display)
    ax.set_title(f"$m_A={ma_val:g}$ GeV, $\\lambda_6={l6_val:g}$")
    ax.grid(True, linestyle="--", alpha=0.4, which="both")
    ax.legend(loc="best")
    fig.tight_layout()

    png_path = OUTPUT_DIR / f"{output_stem}.png"
    pdf_path = OUTPUT_DIR / f"{output_stem}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return True


def _plot_fixed_mode(
    df: pd.DataFrame,
    args: argparse.Namespace,
    mphi_col: str,
    br_col: str,
) -> tuple[Path, Path]:
    fig, ax = plt.subplots()
    x = df[mphi_col].to_numpy()
    y = df[br_col].to_numpy()
    valid = np.isfinite(x) & np.isfinite(y)
    if not valid.any():
        raise RuntimeError("No valid points in fixed-cut mode after NaN filtering.")

    ax.scatter(x[valid], y[valid], s=args.point_size, alpha=args.point_alpha, color="#d62728", edgecolors="none")
    ax.set_yscale("log")
    ax.set_xlabel(r"$m_\phi$ [GeV]")
    ax.set_ylabel(r"BR($\phi \to \gamma\gamma$)" if "gaga" in br_col.lower() else f"BR ({br_col})")
    ax.grid(True, linestyle="--", alpha=0.4, which="both")

    title = (
        rf"$\sin(\beta-\alpha)={args.sin_ba:g}$, "
        rf"$\tan\beta={args.tan_beta:g}$, "
        rf"$m_A={args.mA:g}$, "
        rf"$\lambda_6={args.lambda6:g}$, "
        rf"$\lambda_7={args.lambda7:g}$"
    )
    ax.set_title(title)
    fig.tight_layout()

    output_stem = (
        f"mphi_vs_{br_col}_sinba{args.sin_ba:g}_tb{args.tan_beta:g}_"
        f"mA{args.mA:g}_l6{args.lambda6:g}_l7{args.lambda7:g}"
    ).replace(".", "p")

    png_path = OUTPUT_DIR / f"{output_stem}.png"
    pdf_path = OUTPUT_DIR / f"{output_stem}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, pdf_path


def _run_family_mode(lf: pl.LazyFrame, args: argparse.Namespace, mphi_col: str, br_col: str) -> dict[str, object]:
    logger.info("Running family mode (pipeline-compatible).")

    needed = [mphi_col, br_col, "mA", "lambda6", "tan_beta"]
    lf = lf.select(needed)

    if args.apply_phys_filter and args.present_flags:
        expr = None
        for c in args.present_flags:
            e = _flag_true_expr(c)
            expr = e if expr is None else (expr & e)
        lf = lf.filter(expr)

    lf = lf.filter(pl.col(mphi_col).is_not_null() & pl.col(br_col).is_not_null())

    unique_params = _stream_collect(lf.select(["mA", "lambda6"]).unique().sort(["mA", "lambda6"]))
    unique_df = unique_params.to_pandas()
    if args.max_slices is not None:
        unique_df = unique_df.head(args.max_slices)

    summary: dict[str, object] = {
        "mode": "family",
        "total_combinations_investigated": int(len(unique_df)),
        "plots_generated": [],
        "discarded_combinations": [],
    }

    for _, row in unique_df.iterrows():
        ma_val = float(row["mA"])
        l6_val = float(row["lambda6"])

        df_slice = _stream_collect(
            lf.filter((pl.col("mA") == ma_val) & (pl.col("lambda6") == l6_val)).select([mphi_col, br_col, "tan_beta"])
        ).to_pandas()

        if len(df_slice) == 0:
            summary["discarded_combinations"].append({"mA": ma_val, "lambda6": l6_val, "reason": "empty"})
            continue

        tb_vals = df_slice["tan_beta"].dropna().unique().tolist()
        if len(tb_vals) > 20:
            sorted_tb = sorted(tb_vals)
            idxs = np.linspace(0, len(sorted_tb) - 1, 12, dtype=int)
            tb_vals = [sorted_tb[i] for i in idxs]

        output_stem = f"mphi_vs_{br_col}_mA{ma_val:g}_l6_{l6_val:g}".replace(".", "p")
        ok = _plot_scatter_family(
            df_slice,
            ma_val,
            l6_val,
            tb_vals,
            mphi_col,
            br_col,
            args.point_size,
            args.point_alpha,
            output_stem,
        )

        if ok:
            summary["plots_generated"].append(
                {
                    "mA": ma_val,
                    "lambda6": l6_val,
                    "tan_beta_groups": tb_vals,
                    "points_in_slice": int(len(df_slice)),
                    "file_stem": output_stem,
                }
            )
        else:
            summary["discarded_combinations"].append(
                {"mA": ma_val, "lambda6": l6_val, "reason": "no valid point groups"}
            )

    return summary


def _run_fixed_mode(lf: pl.LazyFrame, args: argparse.Namespace, mphi_col: str, br_col: str) -> dict[str, object]:
    logger.info("Running fixed-cut mode (legacy-compatible).")

    required = ["sin_ba", "tan_beta", "mA", "lambda6", "lambda7", mphi_col, br_col]
    df = _stream_collect(lf.select(required).filter(
        (pl.col("sin_ba") == float(args.sin_ba))
        & (pl.col("tan_beta") == float(args.tan_beta))
        & (pl.col("mA") == float(args.mA))
        & (pl.col("lambda6") == float(args.lambda6))
        & (pl.col("lambda7") == float(args.lambda7))
    )).to_pandas()

    if df.empty:
        raise RuntimeError("No rows found for fixed cuts. Check --tan-beta/--mA/--lambda6/--lambda7/--sin-ba values.")

    png_path, pdf_path = _plot_fixed_mode(df, args, mphi_col, br_col)
    return {
        "mode": "fixed",
        "n_rows_after_cuts": int(len(df)),
        "files": {"png": str(png_path), "pdf": str(pdf_path)},
    }


def main() -> None:
    global OUTPUT_DIR
    parser = argparse.ArgumentParser(description="Generate m_phi vs BR plots (family mode + fixed-cut compatibility).")
    parser.add_argument("--input", type=str, default=str(DEFAULT_PARQUET), help="Path to parquet")
    parser.add_argument("--output-dir", type=str, default=str(OUTPUT_DIR), help="Directory where plots and summary are written")
    parser.add_argument("--mphi-col", type=str, default="m_phi", help="Name of m_phi column")
    parser.add_argument("--br-col", type=str, default="br_gaga", help="Name of BR column")
    parser.add_argument("--apply-phys-filter", action="store_true", help="Apply *_ok filters when present")
    parser.add_argument("--max-slices", type=int, default=None, help="Optional cap for (mA, lambda6) slices in family mode")
    parser.add_argument("--point-size", type=float, default=2.0, help="Scatter marker size")
    parser.add_argument("--point-alpha", type=float, default=0.5, help="Scatter marker alpha")

    # Legacy/fixed-cut compatibility args.
    parser.add_argument("--sin-ba", type=float, default=1.0, help="Fixed sin_ba cut (fixed mode)")
    parser.add_argument("--tan-beta", type=float, default=None, help="Fixed tan_beta cut (fixed mode)")
    parser.add_argument("--mA", type=float, default=None, help="Fixed mA cut (fixed mode)")
    parser.add_argument("--lambda6", type=float, default=None, help="Fixed lambda6 cut (fixed mode)")
    parser.add_argument("--lambda7", type=float, default=None, help="Fixed lambda7 cut (fixed mode)")

    args = parser.parse_args()
    OUTPUT_DIR = Path(args.output_dir)

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    schema = pl.scan_parquet(input_path).collect_schema()
    schema_names = schema.names()

    mphi_col = resolve_column(args.mphi_col, schema_names, ["mphi"])
    br_col = resolve_column(args.br_col, schema_names, ["branching_ratio", "br_gg"])

    present_flags = [c for c in ["positivity_ok", "unitarity_ok", "perturbativity_ok"] if c in schema_names]
    args.present_flags = present_flags

    needed_base = [mphi_col, br_col, "mA", "lambda6", "tan_beta"]
    for c in ["sin_ba", "lambda7"]:
        if c in schema_names and c not in needed_base:
            needed_base.append(c)
    for c in present_flags:
        if c not in needed_base:
            needed_base.append(c)

    lf = pl.scan_parquet(input_path).select(needed_base)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Fixed mode is enabled only when all fixed parameters are explicitly provided.
    use_fixed_mode = all(v is not None for v in [args.tan_beta, args.mA, args.lambda6, args.lambda7])

    if use_fixed_mode:
        result = _run_fixed_mode(lf, args, mphi_col, br_col)
    else:
        result = _run_family_mode(lf, args, mphi_col, br_col)

    summary = {
        "input_file": str(input_path),
        "output_dir": str(OUTPUT_DIR),
        "resolved_columns": {"mphi_col": mphi_col, "br_col": br_col},
        "apply_phys_filter": bool(args.apply_phys_filter),
        "result": result,
    }

    summary_path = OUTPUT_DIR / "summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    logger.info("Done. Summary written to %s", summary_path)


if __name__ == "__main__":
    main()
