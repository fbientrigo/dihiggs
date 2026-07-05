#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
paper_like_mphi_vs_br_patched.py
================================
Dual-mode plotting utility:
1) Family mode (default): m_phi vs BR for each (mA, lambda6), colored by tan_beta.
2) Fixed-cut mode: single plot after applying fixed cuts on
   (sin_ba, tan_beta, mA, lambda6, lambda7).

Extended with optional lambda1-aware filtering/coloring while preserving
legacy defaults.
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
from matplotlib.lines import Line2D


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


def _parse_float_bins(raw_bins: str | None) -> list[float] | None:
    if not raw_bins:
        return None
    vals = [float(x.strip()) for x in raw_bins.split(",") if x.strip()]
    if len(vals) < 2:
        raise ValueError("--lambda1-bins requires at least two comma-separated edges.")
    if any(not np.isfinite(v) for v in vals):
        raise ValueError("--lambda1-bins contains non-finite values.")
    if any(vals[i] >= vals[i + 1] for i in range(len(vals) - 1)):
        raise ValueError("--lambda1-bins must be strictly increasing.")
    return vals


def _format_num_for_stem(v: float) -> str:
    return f"{v:g}".replace("-", "m").replace(".", "p")


def _build_suffix_parts(args: argparse.Namespace, include_lambda1_filters: bool) -> list[str]:
    parts: list[str] = []
    if include_lambda1_filters:
        lo = "min" if args.lambda1_min is None else _format_num_for_stem(float(args.lambda1_min))
        hi = "max" if args.lambda1_max is None else _format_num_for_stem(float(args.lambda1_max))
        parts.append(f"lam1_{lo}_{hi}")
    if args.color_by == "lambda1":
        parts.append("color_lambda1")
    if args.lambda1_bins:
        parts.append("lambda1_bins")
    elif args.lambda1_quantile_bins:
        parts.append(f"lambda1_qbins{int(args.lambda1_quantile_bins)}")
    return parts


def _add_lambda1_stats(record: dict[str, object], df: pd.DataFrame, lambda1_col: str | None) -> None:
    if lambda1_col is None or lambda1_col not in df.columns or df.empty:
        return
    s = pd.to_numeric(df[lambda1_col], errors="coerce").dropna()
    if len(s) == 0:
        return
    record.update(
        {
            "lambda1_min": float(s.min()),
            "lambda1_max": float(s.max()),
            "lambda1_q05": float(np.quantile(s, 0.05)),
            "lambda1_q50": float(np.quantile(s, 0.50)),
            "lambda1_q95": float(np.quantile(s, 0.95)),
        }
    )


def _downsample_df(df: pd.DataFrame, max_points: int | None, seed: int) -> tuple[pd.DataFrame, bool]:
    if max_points is None or len(df) <= max_points:
        return df, False
    sampled = df.sample(n=max_points, random_state=seed)
    return sampled, True


def _plot_family_tanbeta(
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


def _plot_family_lambda1(
    df: pd.DataFrame,
    ma_val: float,
    l6_val: float,
    mphi_col: str,
    br_col: str,
    lambda1_col: str,
    point_size: float,
    point_alpha: float,
    output_stem: str,
    lambda1_bins: list[float] | None,
    lambda1_cmap: str,
) -> bool:
    fig, ax = plt.subplots()
    x = pd.to_numeric(df[mphi_col], errors="coerce").to_numpy()
    y = pd.to_numeric(df[br_col], errors="coerce").to_numpy()
    l1 = pd.to_numeric(df[lambda1_col], errors="coerce").to_numpy()

    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(l1)
    if not valid.any():
        plt.close(fig)
        return False

    x = x[valid]
    y = y[valid]
    l1 = l1[valid]

    if lambda1_bins is None:
        sc = ax.scatter(
            x,
            y,
            c=l1,
            cmap=lambda1_cmap,
            alpha=point_alpha,
            s=point_size,
            edgecolors="none",
            zorder=2,
        )
        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label(r"$\lambda_1$")
    else:
        bin_idx = np.digitize(l1, lambda1_bins, right=False) - 1
        nbins = len(lambda1_bins) - 1
        mask = (bin_idx >= 0) & (bin_idx < nbins)
        if not mask.any():
            plt.close(fig)
            return False
        x, y, l1, bin_idx = x[mask], y[mask], l1[mask], bin_idx[mask]
        cmap = plt.get_cmap(lambda1_cmap, nbins)
        colors = cmap(np.arange(nbins))
        labels: list[str] = []
        handles: list[Line2D] = []
        for i in range(nbins):
            m = bin_idx == i
            if not m.any():
                continue
            ax.scatter(
                x[m],
                y[m],
                color=colors[i],
                alpha=point_alpha,
                s=point_size,
                edgecolors="none",
                zorder=2,
            )
            lo, hi = lambda1_bins[i], lambda1_bins[i + 1]
            labels.append(f"[{lo:g}, {hi:g})")
            handles.append(Line2D([0], [0], marker="o", color="none", markerfacecolor=colors[i], markersize=6, label=labels[-1]))
        if handles:
            ax.legend(handles=handles, title=r"$\lambda_1$ bins", loc="best")

    ax.set_yscale("log")
    ax.set_xlabel(r"$m_\phi$ [GeV]")
    br_display = r"BR($\phi \to \gamma\gamma$)" if "gaga" in br_col.lower() else f"BR ({br_col})"
    ax.set_ylabel(br_display)

    n_tb = int(pd.Series(df["tan_beta"]).dropna().nunique()) if "tan_beta" in df.columns else 0
    subtitle = f"; tan_beta groups={n_tb}" if n_tb > 1 else ""
    ax.set_title(f"$m_A={ma_val:g}$ GeV, $\\lambda_6={l6_val:g}${subtitle}")
    ax.grid(True, linestyle="--", alpha=0.4, which="both")
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
    lambda1_col: str | None,
    output_stem_suffix: str,
) -> tuple[Path, Path]:
    fig, ax = plt.subplots()
    x = pd.to_numeric(df[mphi_col], errors="coerce").to_numpy()
    y = pd.to_numeric(df[br_col], errors="coerce").to_numpy()

    if args.color_by == "lambda1":
        if lambda1_col is None:
            raise RuntimeError("--color-by lambda1 requested, but no lambda1/lam1/lambda_1 column is available.")
        l1 = pd.to_numeric(df[lambda1_col], errors="coerce").to_numpy()
        valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(l1)
    else:
        valid = np.isfinite(x) & np.isfinite(y)

    if not valid.any():
        raise RuntimeError("No valid points in fixed-cut mode after NaN filtering.")

    x, y = x[valid], y[valid]
    if args.color_by == "lambda1":
        l1 = l1[valid]
        if args.lambda1_bins_values is None:
            sc = ax.scatter(x, y, c=l1, cmap=args.lambda1_cmap, s=args.point_size, alpha=args.point_alpha, edgecolors="none")
            cbar = fig.colorbar(sc, ax=ax)
            cbar.set_label(r"$\lambda_1$")
        else:
            bins = args.lambda1_bins_values
            idx = np.digitize(l1, bins, right=False) - 1
            nbins = len(bins) - 1
            mask = (idx >= 0) & (idx < nbins)
            if not mask.any():
                raise RuntimeError("No points fall inside provided --lambda1-bins in fixed mode.")
            x, y, idx = x[mask], y[mask], idx[mask]
            cmap = plt.get_cmap(args.lambda1_cmap, nbins)
            handles: list[Line2D] = []
            for i in range(nbins):
                m = idx == i
                if not m.any():
                    continue
                color = cmap(i)
                ax.scatter(x[m], y[m], color=color, s=args.point_size, alpha=args.point_alpha, edgecolors="none")
                lo, hi = bins[i], bins[i + 1]
                handles.append(Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markersize=6, label=f"[{lo:g}, {hi:g})"))
            if handles:
                ax.legend(handles=handles, title=r"$\lambda_1$ bins", loc="best")
    else:
        ax.scatter(x, y, s=args.point_size, alpha=args.point_alpha, color="#d62728", edgecolors="none")

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
        f"mA{args.mA:g}_l6{args.lambda6:g}_l7{args.lambda7:g}{output_stem_suffix}"
    ).replace(".", "p")

    png_path = OUTPUT_DIR / f"{output_stem}.png"
    pdf_path = OUTPUT_DIR / f"{output_stem}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, pdf_path


def _run_family_mode(lf: pl.LazyFrame, args: argparse.Namespace, mphi_col: str, br_col: str, lambda1_col: str | None) -> dict[str, object]:
    logger.info("Running family mode (pipeline-compatible).")

    needed = [mphi_col, br_col, "mA", "lambda6", "tan_beta"]
    if lambda1_col and lambda1_col not in needed:
        needed.append(lambda1_col)
    lf = lf.select(needed)

    if args.apply_phys_filter and args.present_flags:
        expr = None
        for c in args.present_flags:
            e = _flag_true_expr(c)
            expr = e if expr is None else (expr & e)
        lf = lf.filter(expr)

    lf = lf.filter(pl.col(mphi_col).is_not_null() & pl.col(br_col).is_not_null())

    rows_before_lambda1 = int(_stream_collect(lf.select(pl.len().alias("n"))).item(0, 0))
    if args.use_lambda1_filters:
        if lambda1_col is None:
            raise ValueError("Lambda1 filters requested, but no lambda1-compatible column (lambda1/lam1/lambda_1) exists.")
        if args.lambda1_min is not None:
            lf = lf.filter(pl.col(lambda1_col) >= float(args.lambda1_min))
        if args.lambda1_max is not None:
            lf = lf.filter(pl.col(lambda1_col) <= float(args.lambda1_max))
    rows_after_lambda1 = int(_stream_collect(lf.select(pl.len().alias("n"))).item(0, 0))

    if args.color_by == "lambda1" and lambda1_col is None:
        raise ValueError("--color-by lambda1 requested, but no lambda1-compatible column exists.")

    unique_params = _stream_collect(lf.select(["mA", "lambda6"]).unique().sort(["mA", "lambda6"]))
    unique_df = unique_params.to_pandas()
    if args.max_slices is not None:
        unique_df = unique_df.head(args.max_slices)

    summary: dict[str, object] = {
        "mode": "family",
        "color_by": args.color_by,
        "rows_before_lambda1_filter": rows_before_lambda1,
        "rows_after_lambda1_filter": rows_after_lambda1,
        "total_combinations_investigated": int(len(unique_df)),
        "plots_generated": [],
        "discarded_combinations": [],
    }

    suffix = ""
    suffix_parts = _build_suffix_parts(args, args.use_lambda1_filters)
    if suffix_parts:
        suffix = "_" + "_".join(suffix_parts)

    for _, row in unique_df.iterrows():
        ma_val = float(row["mA"])
        l6_val = float(row["lambda6"])

        select_cols = [mphi_col, br_col, "tan_beta"]
        if lambda1_col and lambda1_col not in select_cols:
            select_cols.append(lambda1_col)

        df_slice = _stream_collect(
            lf.filter((pl.col("mA") == ma_val) & (pl.col("lambda6") == l6_val)).select(select_cols)
        ).to_pandas()

        if len(df_slice) == 0:
            summary["discarded_combinations"].append({"mA": ma_val, "lambda6": l6_val, "reason": "empty"})
            continue

        rows_in_slice_before_downsampling = int(len(df_slice))
        df_plot, downsample_applied = _downsample_df(df_slice, args.max_points_per_slice, args.random_seed)

        output_stem = f"mphi_vs_{br_col}_mA{ma_val:g}_l6_{l6_val:g}{suffix}".replace(".", "p")

        if args.color_by == "tan_beta":
            tb_vals = df_plot["tan_beta"].dropna().unique().tolist()
            if len(tb_vals) > 20:
                sorted_tb = sorted(tb_vals)
                idxs = np.linspace(0, len(sorted_tb) - 1, 12, dtype=int)
                tb_vals = [sorted_tb[i] for i in idxs]
            ok = _plot_family_tanbeta(
                df_plot,
                ma_val,
                l6_val,
                tb_vals,
                mphi_col,
                br_col,
                args.point_size,
                args.point_alpha,
                output_stem,
            )
        else:
            ok = _plot_family_lambda1(
                df_plot,
                ma_val,
                l6_val,
                mphi_col,
                br_col,
                lambda1_col=lambda1_col,
                point_size=args.point_size,
                point_alpha=args.point_alpha,
                output_stem=output_stem,
                lambda1_bins=args.lambda1_bins_values,
                lambda1_cmap=args.lambda1_cmap,
            )

        if ok:
            rec: dict[str, object] = {
                "mA": ma_val,
                "lambda6": l6_val,
                "points_in_slice": rows_in_slice_before_downsampling,
                "points_plotted": int(len(df_plot)),
                "downsampling_applied": bool(downsample_applied),
                "file_stem": output_stem,
                "tan_beta_unique_count": int(df_slice["tan_beta"].dropna().nunique()),
            }
            if args.color_by == "tan_beta":
                rec["tan_beta_groups"] = sorted(df_plot["tan_beta"].dropna().unique().tolist())
            _add_lambda1_stats(rec, df_plot, lambda1_col)
            summary["plots_generated"].append(rec)
        else:
            summary["discarded_combinations"].append(
                {"mA": ma_val, "lambda6": l6_val, "reason": "no valid points"}
            )

    return summary


def _run_fixed_mode(lf: pl.LazyFrame, args: argparse.Namespace, mphi_col: str, br_col: str, lambda1_col: str | None) -> dict[str, object]:
    logger.info("Running fixed-cut mode (legacy-compatible).")

    required = ["sin_ba", "tan_beta", "mA", "lambda6", "lambda7", mphi_col, br_col]
    if lambda1_col and lambda1_col not in required:
        required.append(lambda1_col)

    lf_fixed = lf.select(required).filter(
        (pl.col("sin_ba") == float(args.sin_ba))
        & (pl.col("tan_beta") == float(args.tan_beta))
        & (pl.col("mA") == float(args.mA))
        & (pl.col("lambda6") == float(args.lambda6))
        & (pl.col("lambda7") == float(args.lambda7))
    )

    rows_before_lambda1 = int(_stream_collect(lf_fixed.select(pl.len().alias("n"))).item(0, 0))

    if args.use_lambda1_filters:
        if lambda1_col is None:
            raise ValueError("Lambda1 filters requested, but no lambda1-compatible column (lambda1/lam1/lambda_1) exists.")
        if args.lambda1_min is not None:
            lf_fixed = lf_fixed.filter(pl.col(lambda1_col) >= float(args.lambda1_min))
        if args.lambda1_max is not None:
            lf_fixed = lf_fixed.filter(pl.col(lambda1_col) <= float(args.lambda1_max))

    rows_after_lambda1 = int(_stream_collect(lf_fixed.select(pl.len().alias("n"))).item(0, 0))

    df = _stream_collect(lf_fixed).to_pandas()

    if df.empty:
        raise RuntimeError("No rows found for fixed cuts. Check fixed cuts and lambda1 filters.")

    df_plot, downsample_applied = _downsample_df(df, args.max_points_per_slice, args.random_seed)

    suffix_parts = _build_suffix_parts(args, args.use_lambda1_filters)
    output_stem_suffix = ""
    if suffix_parts:
        output_stem_suffix = "_" + "_".join(suffix_parts)

    png_path, pdf_path = _plot_fixed_mode(df_plot, args, mphi_col, br_col, lambda1_col, output_stem_suffix)

    rec: dict[str, object] = {
        "mode": "fixed",
        "color_by": args.color_by,
        "rows_before_lambda1_filter": rows_before_lambda1,
        "rows_after_lambda1_filter": rows_after_lambda1,
        "n_rows_after_cuts": int(len(df)),
        "points_plotted": int(len(df_plot)),
        "downsampling_applied": bool(downsample_applied),
        "files": {"png": str(png_path), "pdf": str(pdf_path)},
    }
    _add_lambda1_stats(rec, df_plot, lambda1_col)
    return rec


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

    # Lambda1 options (backward compatible defaults).
    parser.add_argument("--lambda1-col", type=str, default="lambda1", help="Preferred lambda1 column name (aliases: lam1, lambda_1)")
    parser.add_argument("--lambda1-min", type=float, default=None, help="Optional lower filter on lambda1")
    parser.add_argument("--lambda1-max", type=float, default=None, help="Optional upper filter on lambda1")
    parser.add_argument("--lambda1-bins", type=str, default=None, help="Optional comma-separated lambda1 bin edges, e.g. 0.2,1,2,3")
    parser.add_argument("--color-by", type=str, choices=["tan_beta", "lambda1"], default="tan_beta", help="Color variable for scatter")
    parser.add_argument("--lambda1-cmap", type=str, default="plasma", help="Colormap for lambda1 color mode")
    parser.add_argument("--lambda1-quantile-bins", type=int, default=None, help="If set and --lambda1-bins omitted, use global quantile bins")
    parser.add_argument("--max-points-per-slice", type=int, default=None, help="Optional deterministic per-slice downsampling cap")
    parser.add_argument("--random-seed", type=int, default=42, help="Deterministic random seed for downsampling")

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

    if args.lambda1_quantile_bins is not None and args.lambda1_quantile_bins < 2:
        raise ValueError("--lambda1-quantile-bins must be >= 2.")

    schema = pl.scan_parquet(input_path).collect_schema()
    schema_names = schema.names()

    mphi_col = resolve_column(args.mphi_col, schema_names, ["mphi"])
    br_col = resolve_column(args.br_col, schema_names, ["branching_ratio", "br_gg"])

    lambda1_related_requested = any(
        [
            args.lambda1_min is not None,
            args.lambda1_max is not None,
            args.lambda1_bins is not None,
            args.lambda1_quantile_bins is not None,
            args.color_by == "lambda1",
            args.lambda1_col.lower() != "lambda1",
        ]
    )
    args.use_lambda1_filters = (args.lambda1_min is not None) or (args.lambda1_max is not None)

    lambda1_col: str | None = None
    if lambda1_related_requested:
        lambda1_col = resolve_column(args.lambda1_col, schema_names, ["lam1", "lambda_1", "lambda1"])

    present_flags = [c for c in ["positivity_ok", "unitarity_ok", "perturbativity_ok"] if c in schema_names]
    args.present_flags = present_flags

    needed_base = [mphi_col, br_col, "mA", "lambda6", "tan_beta"]
    for c in ["sin_ba", "lambda7"]:
        if c in schema_names and c not in needed_base:
            needed_base.append(c)
    for c in present_flags:
        if c not in needed_base:
            needed_base.append(c)
    if lambda1_col and lambda1_col not in needed_base:
        needed_base.append(lambda1_col)

    lf = pl.scan_parquet(input_path).select(needed_base)

    args.lambda1_bins_values = _parse_float_bins(args.lambda1_bins)

    # If explicit bins are absent and quantile bins requested, compute global edges once.
    # Global is chosen as the safer/simple behavior to keep plot-to-plot color meaning consistent.
    if args.lambda1_bins_values is None and args.lambda1_quantile_bins is not None:
        if lambda1_col is None:
            raise ValueError("--lambda1-quantile-bins requested but no lambda1-compatible column exists.")
        q = np.linspace(0.0, 1.0, int(args.lambda1_quantile_bins) + 1)
        qvals = _stream_collect(
            lf.select([pl.col(lambda1_col).quantile(float(v)).alias(f"q{int(v*1000):03d}") for v in q])
        ).to_pandas().iloc[0].to_numpy(dtype=float)
        qvals = np.unique(qvals[np.isfinite(qvals)])
        if len(qvals) >= 2:
            args.lambda1_bins_values = qvals.tolist()
        else:
            logger.warning("Could not derive valid global quantile bins for lambda1; falling back to continuous coloring.")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    use_fixed_mode = all(v is not None for v in [args.tan_beta, args.mA, args.lambda6, args.lambda7])

    if use_fixed_mode:
        result = _run_fixed_mode(lf, args, mphi_col, br_col, lambda1_col)
    else:
        result = _run_family_mode(lf, args, mphi_col, br_col, lambda1_col)

    summary = {
        "input_file": str(input_path),
        "output_dir": str(OUTPUT_DIR),
        "resolved_columns": {
            "mphi_col": mphi_col,
            "br_col": br_col,
            "lambda1_col": lambda1_col,
        },
        "apply_phys_filter": bool(args.apply_phys_filter),
        "color_by": args.color_by,
        "lambda1_filter_min": args.lambda1_min,
        "lambda1_filter_max": args.lambda1_max,
        "lambda1_bins": args.lambda1_bins_values,
        "lambda1_quantile_bins": args.lambda1_quantile_bins,
        "max_points_per_slice": args.max_points_per_slice,
        "random_seed": args.random_seed,
        "result": result,
    }

    summary_path = OUTPUT_DIR / "summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    logger.info("Done. Summary written to %s", summary_path)


if __name__ == "__main__":
    main()
