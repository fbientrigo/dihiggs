#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import polars as pl

plt.style.use("default")
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 13,
        "axes.labelsize": 15,
        "axes.titlesize": 16,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 11,
        "legend.frameon": True,
        "legend.edgecolor": "black",
        "legend.fancybox": False,
        "figure.figsize": (9, 10),
        "axes.linewidth": 1.4,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    }
)

HBAR_C_GEV_MM = 1.973269804e-13

STYLE_PRESETS = {
    "paper": {"point_size": 8.0, "point_alpha": 0.55, "guide_linewidth": 1.2, "label_fontsize": 10.0, "fig_width": 9.0, "fig_height": 10.0, "grid_alpha": 0.25},
    "talk": {"point_size": 16.0, "point_alpha": 0.75, "guide_linewidth": 1.8, "label_fontsize": 13.0, "fig_width": 11.0, "fig_height": 6.5, "grid_alpha": 0.30},
    "draft": {"point_size": 10.0, "point_alpha": 0.65, "guide_linewidth": 1.3, "label_fontsize": 10.5, "fig_width": 9.5, "fig_height": 10.0, "grid_alpha": 0.35},
    "compact": {"point_size": 6.5, "point_alpha": 0.50, "guide_linewidth": 1.0, "label_fontsize": 9.0, "fig_width": 7.5, "fig_height": 8.0, "grid_alpha": 0.22},
}

CHANNEL_COLORS = {
    "default": {"br_gaga": "tab:blue", "br_Zga": "tab:orange", "br_bb": "tab:green"},
    "grayscale": {"br_gaga": "0.15", "br_Zga": "0.45", "br_bb": "0.70"},
    "classic": {"br_gaga": "#1f77b4", "br_Zga": "#d62728", "br_bb": "#2ca02c"},
}

ATOL = {"mA": 1e-8, "tan_beta": 1e-3, "lambda6": 1e-7, "lambda7": 1e-9, "sin_ba": 1e-9}
logger = logging.getLogger(__name__)


def _stream_collect(lf: pl.LazyFrame) -> pl.DataFrame:
    try:
        return lf.collect(engine="streaming")
    except TypeError:
        return lf.collect(streaming=True)
    except Exception:
        return lf.collect()


def _resolve_column(requested: str, schema_names: list[str], aliases: list[str] | None = None) -> str:
    if requested in schema_names:
        return requested
    lower_map = {c.lower(): c for c in schema_names}
    if requested.lower() in lower_map:
        return lower_map[requested.lower()]
    for a in aliases or []:
        if a.lower() in lower_map:
            resolved = lower_map[a.lower()]
            logger.info("Resolved '%s' -> '%s'", requested, resolved)
            return resolved
    raise ValueError(f"Could not resolve column '{requested}'.")


def _parse_pair(raw: str | None, arg_name: str) -> tuple[float, float] | None:
    if not raw:
        return None
    vals = [float(x.strip()) for x in raw.split(",") if x.strip()]
    if len(vals) != 2:
        raise ValueError(f"{arg_name} requires exactly two comma-separated values.")
    lo, hi = vals
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        raise ValueError(f"{arg_name} requires finite values with high > low.")
    return float(lo), float(hi)


def _parse_float_list(raw: str | None) -> list[float]:
    if not raw:
        return []
    vals = [float(x.strip()) for x in raw.split(",") if x.strip()]
    return [float(v) for v in vals if np.isfinite(v)]


def _float_tol_cut(col: str, target: float, atol: float) -> pl.Expr:
    return (pl.col(col).cast(pl.Float64) - pl.lit(float(target))).abs() <= float(atol)


def _series_stats(x: np.ndarray) -> dict[str, float | None]:
    vals = x[np.isfinite(x)]
    if vals.size == 0:
        return {"min": None, "q05": None, "q50": None, "q95": None, "max": None}
    return {
        "min": float(np.min(vals)),
        "q05": float(np.quantile(vals, 0.05)),
        "q50": float(np.quantile(vals, 0.50)),
        "q95": float(np.quantile(vals, 0.95)),
        "max": float(np.max(vals)),
    }


def _best_point(df: pd.DataFrame, metric_col: str, mode: str) -> dict[str, float] | None:
    if df.empty or metric_col not in df.columns:
        return None
    s = pd.to_numeric(df[metric_col], errors="coerce")
    d = df.loc[np.isfinite(s.to_numpy())].copy()
    if d.empty:
        return None
    row = d.loc[s.loc[d.index].idxmax()] if mode == "max" else d.loc[s.loc[d.index].idxmin()]
    return {
        "m_phi": float(row["m_phi"]),
        "lambda1": float(row["lambda1"]),
        "br_gaga": float(row["br_gaga"]),
        "br_Zga": float(row["br_Zga"]),
        "br_bb": float(row["br_bb"]),
        "ctau_mm": float(row["ctau_mm"]),
        "m12": float(row["m12"]) if "m12" in row and np.isfinite(row["m12"]) else None,
        metric_col: float(row[metric_col]),
    }


def _channel_guide_curve(x: np.ndarray, y: np.ndarray, bins: int, min_count: int) -> tuple[np.ndarray, np.ndarray]:
    valid = np.isfinite(x) & np.isfinite(y) & (y > 0)
    if not valid.any():
        return np.array([], dtype=float), np.array([], dtype=float)
    xv = x[valid]
    yv = y[valid]
    if xv.size < max(2, int(min_count)):
        return np.array([], dtype=float), np.array([], dtype=float)
    edges = np.linspace(float(np.nanmin(xv)), float(np.nanmax(xv)), int(max(2, bins)) + 1)
    idx = np.digitize(xv, edges) - 1
    idx[idx == len(edges) - 1] = len(edges) - 2
    gx, gy = [], []
    for i in range(len(edges) - 1):
        m = idx == i
        if int(np.count_nonzero(m)) < int(min_count):
            continue
        xb = xv[m]
        yb = yv[m]
        yb = yb[np.isfinite(yb) & (yb > 0)]
        if yb.size < int(min_count):
            continue
        gx.append(float(np.median(xb)))
        gy.append(float(np.exp(np.median(np.log(yb)))))
    if len(gx) < 2:
        return np.array([], dtype=float), np.array([], dtype=float)
    order = np.argsort(np.array(gx))
    return np.array(gx)[order], np.array(gy)[order]


def _pick_targets(df: pd.DataFrame, strategy: str, manual_targets: list[float], default_targets: list[float]) -> list[float]:
    lam = pd.to_numeric(df["lambda1"], errors="coerce").to_numpy()
    lam = lam[np.isfinite(lam)]
    if strategy == "manual":
        return manual_targets if manual_targets else default_targets
    if lam.size == 0:
        return []
    if strategy == "quantiles":
        qs = np.quantile(lam, [0.2, 0.5, 0.8])
        return [float(x) for x in qs]

    d = df.copy()
    d["score"] = d["br_gaga"] / np.sqrt(d["br_bb"] + 1e-12)
    g = _best_point(d, "br_gaga", "max")
    p = _best_point(d, "score", "max")
    b = _best_point(d, "br_bb", "min")
    out = []
    for t in [g, p, b]:
        if t is not None and t.get("lambda1") is not None:
            out.append(float(t["lambda1"]))
    dedup = []
    for v in out:
        if not any(abs(v - u) <= 1e-9 for u in dedup):
            dedup.append(v)
    if len(dedup) < 3:
        qs = [float(x) for x in np.quantile(lam, [0.2, 0.5, 0.8])]
        for q in qs:
            if not any(abs(q - u) <= 1e-9 for u in dedup):
                dedup.append(q)
    return dedup[:3]


def _style_list(generate_style_variants: bool, style_preset: str) -> list[str]:
    return ["paper", "talk", "compact"] if generate_style_variants else [style_preset]


def _subtitle_m12_value(df: pd.DataFrame, mode: str, m12_stats: dict[str, float | None], best_gaga: dict[str, float] | None, best_pareto: dict[str, float] | None) -> float | None:
    if mode == "best_gaga" and best_gaga is not None:
        return best_gaga.get("m12")
    if mode == "best_pareto" and best_pareto is not None:
        return best_pareto.get("m12")
    return m12_stats.get("q50")


def _plot_and_save(df: pd.DataFrame, args: argparse.Namespace, out_dir: Path, style_name: str, fixed_mode: bool, lambda1_target: float | None, rows_meta: dict, m12_label_tex: str) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    preset = STYLE_PRESETS[style_name]
    fig_width = float(args.fig_width) if args.fig_width is not None else float(preset["fig_width"])
    fig_height = float(args.fig_height) if args.fig_height is not None else float(preset["fig_height"])
    point_size = float(args.point_size) if args.point_size is not None else float(preset["point_size"])
    point_alpha = float(args.point_alpha) if args.point_alpha is not None else float(preset["point_alpha"])
    channel_label_fontsize = float(args.channel_label_fontsize) if args.channel_label_fontsize is not None else float(preset["label_fontsize"])
    guide_linewidth = float(args.channel_guide_linewidth) if args.channel_guide_linewidth is not None else float(preset["guide_linewidth"])
    grid_alpha = float(args.grid_alpha) if args.grid_alpha is not None else float(preset["grid_alpha"])
    ctau_band = _parse_pair(args.ctau_band, "--ctau-band")

    x_all = pd.to_numeric(df["m_phi"], errors="coerce").to_numpy()
    finite_x = x_all[np.isfinite(x_all)]
    x_min_data = float(np.nanmin(finite_x))
    x_max_data = float(np.nanmax(finite_x))
    x_min = float(args.mphi_min_plot) if args.mphi_min_plot is not None else x_min_data
    x_max = float(args.mphi_max_plot) if args.mphi_max_plot is not None else x_max_data
    x_span = max(x_max - x_min, 1e-9)
    xpad = float(args.channel_label_xpad_frac) * x_span

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(fig_width, fig_height), gridspec_kw={"hspace": 0.06})

    channels = [
        ("br_gaga", "o", r"$\mathrm{BR}(H_2\to\gamma\gamma)$", "-"),
        ("br_Zga", "^", r"$\mathrm{BR}(H_2\to Z\gamma)$", "--"),
        ("br_bb", "s", r"$\mathrm{BR}(H_2\to b\bar b)$", ":"),
    ]

    channel_color_map = CHANNEL_COLORS[args.channel_color_mode]
    lam1_values = pd.to_numeric(df["lambda1"], errors="coerce").to_numpy()
    lam1_min = float(np.nanmin(lam1_values))
    lam1_max = float(np.nanmax(lam1_values))
    norm = Normalize(vmin=lam1_min, vmax=lam1_max)
    cmap = plt.get_cmap(args.lambda1_cmap)

    if args.channel_end_labels:
        if args.channel_label_side == "right":
            ax1.set_xlim(x_min, x_max + xpad * 5.0)
        else:
            ax1.set_xlim(x_min - xpad * 5.0, x_max)
    else:
        ax1.set_xlim(x_min, x_max)

    channel_label_positions = {}
    for col, marker, direct_label, line_style in channels:
        x = pd.to_numeric(df["m_phi"], errors="coerce").to_numpy()
        y = pd.to_numeric(df[col], errors="coerce").to_numpy()
        if fixed_mode:
            valid = np.isfinite(x) & np.isfinite(y) & (y > 0)
            ax1.scatter(x[valid], y[valid], c=channel_color_map[col], s=point_size, alpha=point_alpha, marker=marker, edgecolors="none", linewidths=0, zorder=2)
            guide_color = channel_color_map[col]
        else:
            c = pd.to_numeric(df["lambda1"], errors="coerce").to_numpy()
            valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(c) & (y > 0)
            ax1.scatter(x[valid], y[valid], c=c[valid], cmap=cmap, norm=norm, s=point_size, alpha=point_alpha, marker=marker, edgecolors="none", linewidths=0, zorder=2)
            guide_color = "black"

        gx, gy = _channel_guide_curve(x=x, y=y, bins=int(args.channel_guide_bins), min_count=int(args.channel_guide_min_count))
        if gx.size >= 2:
            if args.channel_guides_enabled and float(guide_linewidth) > 0:
                ax1.plot(gx, gy, color=guide_color, linestyle=line_style, linewidth=float(guide_linewidth), alpha=0.85, zorder=float(args.channel_guide_zorder))
            if args.channel_end_labels:
                if args.channel_label_side == "right":
                    lx = float(gx[-1] + xpad)
                    ha = "left"
                    ly_base = float(gy[-1])
                else:
                    lx = float(gx[0] - xpad)
                    ha = "right"
                    ly_base = float(gy[0])
                ly = float(ly_base * float(args.channel_label_y_scale))
                ax1.text(lx, ly, direct_label, fontsize=float(channel_label_fontsize), va="center", ha=ha, clip_on=False, bbox={"facecolor": "white", "alpha": 0.78, "edgecolor": "none", "pad": 0.2}, zorder=5)
                channel_label_positions[col] = {"x": lx, "y": ly, "n_guide_points": int(gx.size)}
            else:
                channel_label_positions[col] = {"x": float(gx[-1]), "y": float(gy[-1]), "n_guide_points": int(gx.size)}
        else:
            channel_label_positions[col] = None

    x2 = pd.to_numeric(df["m_phi"], errors="coerce").to_numpy()
    y2 = pd.to_numeric(df["ctau_mm"], errors="coerce").to_numpy()
    valid2 = np.isfinite(x2) & np.isfinite(y2) & (y2 > 0)
    if fixed_mode:
        ax2.scatter(x2[valid2], y2[valid2], c="0.35", s=point_size * 0.9, alpha=min(0.85, point_alpha + 0.1), marker="o", edgecolors="none", linewidths=0, zorder=2)
        sc2 = None
    else:
        c2 = pd.to_numeric(df["lambda1"], errors="coerce").to_numpy()
        valid2 = valid2 & np.isfinite(c2)
        sc2 = ax2.scatter(x2[valid2], y2[valid2], c=c2[valid2], cmap=cmap, norm=norm, s=point_size, alpha=point_alpha, marker="o", edgecolors="none", linewidths=0, zorder=2)

    if ctau_band is not None:
        band_lo, band_hi = ctau_band
        ax2.axhspan(band_lo, band_hi, color="gray", alpha=0.12, zorder=0)
        x_for_band = x_min + 0.02 * (x_max - x_min)
        y_for_band = float(np.sqrt(band_lo * band_hi))
        ax2.text(x_for_band, y_for_band, args.ctau_band_label, fontsize=max(8.0, channel_label_fontsize - 1.0), color="dimgray", va="center", ha="left", bbox={"facecolor": "white", "alpha": 0.65, "edgecolor": "none", "pad": 0.2})

    df_anno = df.copy()
    df_anno["score"] = df_anno["br_gaga"] / np.sqrt(df_anno["br_bb"] + 1e-12)
    best_gaga = _best_point(df_anno, "br_gaga", "max")
    best_pareto = _best_point(df_anno, "score", "max")
    best_min_bb = _best_point(df_anno, "br_bb", "min")

    ax1.set_yscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"$m_\phi$ [GeV]")
    ax1.set_ylabel("Branching ratio")
    ax2.set_ylabel(r"$c\tau$ [mm]")
    if not args.no_grid:
        ax1.grid(True, linestyle="--", alpha=grid_alpha, which="both")
        ax2.grid(True, linestyle="--", alpha=grid_alpha, which="both")

    if args.br_ymin is not None or args.br_ymax is not None:
        ax1.set_ylim(bottom=float(args.br_ymin) if args.br_ymin is not None else None, top=float(args.br_ymax) if args.br_ymax is not None else None)
    if args.ctau_ymin is not None or args.ctau_ymax is not None:
        ax2.set_ylim(bottom=float(args.ctau_ymin) if args.ctau_ymin is not None else None, top=float(args.ctau_ymax) if args.ctau_ymax is not None else None)

    if fixed_mode and not args.channel_end_labels:
        legend_handles = [Line2D([0], [0], marker="o", color=CHANNEL_COLORS[args.channel_color_mode]["br_gaga"], linestyle="None", markersize=6, label=r"$H_2\to\gamma\gamma$"), Line2D([0], [0], marker="^", color=CHANNEL_COLORS[args.channel_color_mode]["br_Zga"], linestyle="None", markersize=6, label=r"$H_2\to Z\gamma$"), Line2D([0], [0], marker="s", color=CHANNEL_COLORS[args.channel_color_mode]["br_bb"], linestyle="None", markersize=6, label=r"$H_2\to b\bar{b}$")]
        ax1.legend(handles=legend_handles, title="Channel", loc="best")

    if not fixed_mode and sc2 is not None:
        orientation = "vertical" if args.colorbar_location in ("right", "left") else "horizontal"
        cbar = fig.colorbar(sc2, ax=[ax1, ax2], pad=float(args.colorbar_pad), location=args.colorbar_location, fraction=float(args.colorbar_fraction), shrink=float(args.colorbar_shrink), orientation=orientation)
        cbar.set_label(r"$\lambda_1$")

    m12_stats = _series_stats(pd.to_numeric(df["m12"], errors="coerce").to_numpy())
    subtitle_m12 = _subtitle_m12_value(df_anno, args.subtitle_m12_mode, m12_stats, best_gaga, best_pareto)
    m12_nearly_constant = False
    if m12_stats["q05"] is not None and m12_stats["q95"] is not None and m12_stats["q50"] is not None:
        m12_nearly_constant = abs(float(m12_stats["q95"]) - float(m12_stats["q05"])) <= max(1e-6, 1e-3 * max(1.0, abs(float(m12_stats["q50"]))))

    title = rf"$m_A={args.mA:g}$, $\tan\beta={args.tan_beta:g}$, $\lambda_6={args.lambda6:g}$, $\lambda_7={args.lambda7:g}$"
    if fixed_mode and lambda1_target is not None:
        if m12_nearly_constant and subtitle_m12 is not None:
            subtitle = rf"$\lambda_1={lambda1_target:.6g}$, ${m12_label_tex}={subtitle_m12:.6g}\ \mathrm{{GeV}}^2$"
        elif subtitle_m12 is not None:
            subtitle = rf"$\lambda_1={lambda1_target:.6g}$, $m_{{12,\mathrm{{med}}}}^2={subtitle_m12:.6g}\ \mathrm{{GeV}}^2$"
        else:
            subtitle = rf"$\lambda_1={lambda1_target:.6g}$"
        title = title + "\n" + subtitle
    else:
        title = title + rf", $\lambda_1\in[{lam1_min:.4g},{lam1_max:.4g}]$"

    if args.title_suffix:
        title += f"  {args.title_suffix}"

    fig.suptitle(title, y=0.98)
    tight_right = 0.94 if (not fixed_mode) else 0.98
    fig.tight_layout(rect=[0.0, 0.0, tight_right, 0.96])

    stem_base = f"three_br_ctau_mA{args.mA:g}_tb{args.tan_beta:g}_l6{args.lambda6:g}_l7{args.lambda7:g}".replace(".", "p")
    stem = f"{stem_base}_lambda1_{lambda1_target:.6g}".replace(".", "p") if lambda1_target is not None else stem_base

    output_files = []
    png_path = out_dir / f"{stem}.png"
    pdf_path = out_dir / f"{stem}.pdf"
    fig.savefig(png_path, dpi=int(args.dpi), bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    output_files.extend([str(png_path), str(pdf_path)])
    if args.save_svg:
        svg_path = out_dir / f"{stem}.svg"
        fig.savefig(svg_path, bbox_inches="tight")
        output_files.append(str(svg_path))
    plt.close(fig)

    channel_stats = {}
    for k in ["br_gaga", "br_Zga", "br_bb"]:
        arr = pd.to_numeric(df[k], errors="coerce").to_numpy()
        st = _series_stats(arr)
        st["n"] = int(np.isfinite(arr).sum())
        channel_stats[k] = st
    ctau_stats = _series_stats(pd.to_numeric(df["ctau_mm"], errors="coerce").to_numpy())

    lam_stats = _series_stats(pd.to_numeric(df["lambda1"], errors="coerce").to_numpy())

    summary = {
        "input_path": str(args.input),
        "style_preset": style_name,
        "fixed_lambda1_mode": bool(fixed_mode),
        "channel_color_mode": args.channel_color_mode,
        "fixed_cuts": {"mA": args.mA, "tan_beta": args.tan_beta, "lambda6": args.lambda6, "lambda7": args.lambda7, "sin_ba": args.sin_ba, "atol": ATOL},
        "lambda1_target": lambda1_target,
        "lambda1_tolerance": args.lambda1_tolerance,
        "lambda1_stats_after_filter": lam_stats,
        "rows_before_filters": rows_meta["rows_before_filters"],
        "rows_after_fixed_cuts": rows_meta["rows_after_fixed_cuts"],
        "rows_after_lambda1_filter": rows_meta["rows_after_lambda1_filter"],
        "rows_plotted": int(len(df)),
        "m12_semantics": "m12 interpreted as m12^2 based on PARQUET_SCHEMA.md column description",
        "m12_stats": m12_stats,
        "subtitle_m12_mode": args.subtitle_m12_mode,
        "subtitle_m12_value": subtitle_m12,
        "m12_nearly_constant": bool(m12_nearly_constant),
        "m12_at_best_points": {
            "best_gaga": None if best_gaga is None else best_gaga.get("m12"),
            "best_pareto": None if best_pareto is None else best_pareto.get("m12"),
            "min_bb": None if best_min_bb is None else best_min_bb.get("m12"),
        },
        "best_points": {"max_br_gaga": best_gaga, "max_pareto_like": best_pareto, "min_br_bb": best_min_bb},
        "per_channel_statistics": channel_stats,
        "ctau_statistics": ctau_stats,
        "output_files": output_files,
    }

    summary_path = out_dir / "summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    return {
        "summary_path": str(summary_path),
        "rows_plotted": int(len(df)),
        "lambda1_target": lambda1_target,
        "style": style_name,
        "subtitle_m12_value": subtitle_m12,
        "m12_nearly_constant": bool(m12_nearly_constant),
        "output_files": output_files,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Two-panel publication-like plot: BR channels + ctau.")
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--output-dir", type=str, required=True)
    parser.add_argument("--mA", type=float, required=True)
    parser.add_argument("--tan-beta", type=float, required=True)
    parser.add_argument("--lambda6", type=float, required=True)
    parser.add_argument("--lambda7", type=float, required=True)
    parser.add_argument("--sin-ba", type=float, default=1.0)
    parser.add_argument("--lambda1-col", type=str, default="lambda1")
    parser.add_argument("--lambda1-min", type=float, default=None)
    parser.add_argument("--lambda1-max", type=float, default=None)
    parser.add_argument("--lambda1-cmap", type=str, default="plasma")
    parser.add_argument("--max-points-per-slice", type=int, default=50000)
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--style-preset", type=str, choices=["paper", "talk", "draft", "compact"], default="paper")
    parser.add_argument("--fig-width", type=float, default=None)
    parser.add_argument("--fig-height", type=float, default=None)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--br-ymin", type=float, default=None)
    parser.add_argument("--br-ymax", type=float, default=None)
    parser.add_argument("--ctau-ymin", type=float, default=None)
    parser.add_argument("--ctau-ymax", type=float, default=None)
    parser.add_argument("--mphi-min-plot", type=float, default=None)
    parser.add_argument("--mphi-max-plot", type=float, default=None)
    parser.add_argument("--point-alpha", type=float, default=None)
    parser.add_argument("--point-size", type=float, default=None)
    parser.add_argument("--title-suffix", type=str, default="")
    parser.add_argument("--channel-guide-bins", type=int, default=40)
    parser.add_argument("--channel-guide-min-count", type=int, default=5)
    parser.add_argument("--channel-label-fontsize", type=float, default=None)
    parser.add_argument("--channel-label-xpad-frac", type=float, default=0.025)
    parser.add_argument("--channel-guide-linewidth", type=float, default=None)
    parser.add_argument("--channel-guides-enabled", dest="channel_guides_enabled", action="store_true")
    parser.add_argument("--no-channel-guides", dest="channel_guides_enabled", action="store_false")
    parser.add_argument("--channel-guide-zorder", type=float, default=1.0)
    parser.add_argument("--channel-end-labels", dest="channel_end_labels", action="store_true")
    parser.add_argument("--no-channel-end-labels", dest="channel_end_labels", action="store_false")
    parser.add_argument("--channel-label-side", type=str, choices=["right", "left"], default="right")
    parser.add_argument("--channel-label-y-scale", type=float, default=1.0)
    parser.add_argument("--colorbar-location", type=str, choices=["right", "left", "top", "bottom"], default="right")
    parser.add_argument("--colorbar-pad", type=float, default=0.02)
    parser.add_argument("--colorbar-fraction", type=float, default=0.045)
    parser.add_argument("--colorbar-shrink", type=float, default=0.98)
    parser.add_argument("--grid-alpha", type=float, default=None)
    parser.add_argument("--no-grid", action="store_true")
    parser.add_argument("--ctau-band", type=str, default=None)
    parser.add_argument("--ctau-band-label", type=str, default="mm-scale displaced")
    parser.add_argument("--save-svg", action="store_true")

    parser.add_argument("--fixed-lambda1-mode", action="store_true")
    parser.add_argument("--lambda1-targets", type=str, default="4.5,5.5,6.5")
    parser.add_argument("--lambda1-tolerance", type=float, default=1e-3)
    parser.add_argument("--channel-color-mode", choices=["default", "grayscale", "classic"], default="default")
    parser.add_argument("--generate-style-variants", action="store_true")
    parser.add_argument("--lambda1-target-strategy", choices=["manual", "quantiles", "best"], default="manual")
    parser.add_argument("--subtitle-m12-mode", choices=["median", "best_gaga", "best_pareto"], default="median")

    parser.set_defaults(channel_end_labels=True, channel_guides_enabled=True)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    schema = pl.scan_parquet(input_path).collect_schema().names()
    mphi_col = _resolve_column("m_phi", schema, aliases=["mphi", "mH2", "m_h2", "mH_2", "m_H2", "m_H_2"])
    lam1_col = _resolve_column(args.lambda1_col, schema, aliases=["lambda1", "lam1", "lambda_1"])
    total_width_col = _resolve_column("total_width", schema, aliases=["total_decay_width", "w_total"])
    width_gaga_col = _resolve_column("width_gaga", schema, aliases=["width_gg"])
    width_zga_col = _resolve_column("width_Zga", schema, aliases=["width_zga"])
    width_bb_col = _resolve_column("width_bb", schema, aliases=["width_bbar", "width_bbbar"])
    m12_col = _resolve_column("m12", schema, aliases=["m_12", "m12sq", "m12_2"])

    br_gaga_col = next((cand for cand in ["br_gaga", "BR_gaga", "br_gg", "BR_gg"] if cand in schema), None)

    required = [mphi_col, "mA", "tan_beta", "lambda6", "lambda7", "sin_ba", lam1_col, total_width_col, width_gaga_col, width_zga_col, width_bb_col, m12_col]
    if br_gaga_col:
        required.append(br_gaga_col)

    lf_all = pl.scan_parquet(input_path).select(required)
    rows_before = int(_stream_collect(lf_all.select(pl.len().alias("n"))).item(0, 0))

    fixed_expr = (
        _float_tol_cut("mA", args.mA, ATOL["mA"])
        & _float_tol_cut("tan_beta", args.tan_beta, ATOL["tan_beta"])
        & _float_tol_cut("lambda6", args.lambda6, ATOL["lambda6"])
        & _float_tol_cut("lambda7", args.lambda7, ATOL["lambda7"])
        & _float_tol_cut("sin_ba", args.sin_ba, ATOL["sin_ba"])
    )

    lf_fixed = lf_all.filter(fixed_expr)
    rows_after_fixed = int(_stream_collect(lf_fixed.select(pl.len().alias("n"))).item(0, 0))

    if args.lambda1_min is not None:
        lf_fixed = lf_fixed.filter(pl.col(lam1_col).cast(pl.Float64) >= float(args.lambda1_min))
    if args.lambda1_max is not None:
        lf_fixed = lf_fixed.filter(pl.col(lam1_col).cast(pl.Float64) <= float(args.lambda1_max))

    lf_base = (
        lf_fixed
        .with_columns(
            pl.col(mphi_col).cast(pl.Float64).alias("m_phi"),
            pl.col(lam1_col).cast(pl.Float64).alias("lambda1"),
            pl.col(m12_col).cast(pl.Float64).alias("m12"),
            pl.col(total_width_col).cast(pl.Float64).alias("total_width"),
            pl.col(width_zga_col).cast(pl.Float64).alias("width_Zga"),
            pl.col(width_bb_col).cast(pl.Float64).alias("width_bb"),
            pl.col(width_gaga_col).cast(pl.Float64).alias("width_gaga"),
            (pl.col(br_gaga_col).cast(pl.Float64) if br_gaga_col else (pl.col(width_gaga_col).cast(pl.Float64) / pl.col(total_width_col).cast(pl.Float64))).alias("br_gaga"),
        )
        .with_columns(
            (pl.col("width_Zga") / pl.col("total_width")).alias("br_Zga"),
            (pl.col("width_bb") / pl.col("total_width")).alias("br_bb"),
            (pl.lit(HBAR_C_GEV_MM) / pl.col("total_width")).alias("ctau_mm"),
        )
        .filter(pl.col("total_width") > 0)
        .filter(pl.all_horizontal([pl.col("m_phi").is_finite(), pl.col("lambda1").is_finite(), pl.col("m12").is_finite(), pl.col("br_gaga").is_finite(), pl.col("br_Zga").is_finite(), pl.col("br_bb").is_finite(), pl.col("ctau_mm").is_finite(), pl.col("br_gaga") > 0, pl.col("br_Zga") > 0, pl.col("br_bb") > 0, pl.col("ctau_mm") > 0]))
        .select(["m_phi", "lambda1", "m12", "br_gaga", "br_Zga", "br_bb", "ctau_mm"])
    )

    df_base = _stream_collect(lf_base).to_pandas()
    if df_base.empty:
        raise RuntimeError("No rows remain after fixed cuts and finite/positive filtering.")

    if args.max_points_per_slice is not None and len(df_base) > int(args.max_points_per_slice):
        df_base = df_base.sample(n=int(args.max_points_per_slice), random_state=args.random_seed)

    style_names = _style_list(args.generate_style_variants, args.style_preset)

    if not args.fixed_lambda1_mode:
        rows_meta = {"rows_before_filters": rows_before, "rows_after_fixed_cuts": rows_after_fixed, "rows_after_lambda1_filter": len(df_base)}
        m12_label_tex = "m_{12}^{2}"
        result = _plot_and_save(df_base, args, output_dir, style_names[0], False, None, rows_meta, m12_label_tex)
        logger.info("Done non-fixed mode. Summary: %s", result["summary_path"])
        return

    manual_targets = _parse_float_list(args.lambda1_targets)
    default_targets = [4.5, 5.5, 6.5]
    targets = _pick_targets(df_base, args.lambda1_target_strategy, manual_targets, default_targets)
    if len(targets) == 0:
        targets = default_targets

    skipped_targets = []
    aggregate = {"targets_used": [], "skipped_targets": [], "runs": []}
    m12_label_tex = "m_{12}^{2}"

    for target in targets:
        target_dir = output_dir / f"lambda1_{target:.6g}".replace(".", "p")
        target_df = df_base[np.abs(pd.to_numeric(df_base["lambda1"], errors="coerce").to_numpy() - float(target)) <= float(args.lambda1_tolerance)].copy()
        rows_after_lambda1 = int(len(target_df))
        if rows_after_lambda1 < max(5, args.channel_guide_min_count):
            logger.warning("Skipping lambda1 target %.6g: too few rows (%d)", target, rows_after_lambda1)
            skipped_targets.append({"target": float(target), "rows_after_lambda1_filter": rows_after_lambda1})
            continue

        aggregate["targets_used"].append(float(target))
        for style_name in style_names:
            style_dir = target_dir / style_name
            rows_meta = {"rows_before_filters": rows_before, "rows_after_fixed_cuts": rows_after_fixed, "rows_after_lambda1_filter": rows_after_lambda1}
            result = _plot_and_save(target_df, args, style_dir, style_name, True, float(target), rows_meta, m12_label_tex)
            aggregate["runs"].append(result)
            logger.info("Target %.6g style %s rows=%d", target, style_name, result["rows_plotted"])

    aggregate["skipped_targets"] = skipped_targets
    aggregate_path = output_dir / "summary_all_targets.json"
    with open(aggregate_path, "w", encoding="utf-8") as f:
        json.dump(aggregate, f, indent=2)

    logger.info("Done fixed-lambda1 mode. Aggregate summary: %s", aggregate_path)


if __name__ == "__main__":
    main()
