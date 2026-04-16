#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
parallel_coordinates_roi.py
===========================
Generates parallel-coordinates visualizations from a parquet dataset with:
1) Full-sample plot colored by a chosen variable (default: ctau).
2) ROI-highlighted plot for low-ctau and high-BR trajectories.

The implementation is memory-aware:
- scans only required columns with Polars,
- converts to NumPy once,
- draws lines in chunks via LineCollection.

Important: the drawn coordinates are shown with per-axis numeric limits
(raw or log10 mode), not as global min-max normalized values.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib import colormaps
from matplotlib.collections import LineCollection


plt.style.use("default")
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 12,
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 10,
        "figure.figsize": (12, 7),
        "axes.linewidth": 1.2,
    }
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

DEFAULT_PARQUET = Path(__file__).parent / "temp_subspace_l6_tb_high.parquet"
OUTPUT_DIR = Path(__file__).parent / "parallel_plots"

# Canonical available columns in preferred order.
ALL_COLUMNS = ["tan_beta", "lambda6", "m_phi", "lam1", "total_width", "ctau", "br_gaga"]

# Default axis list for plotting. ctau is omitted by default because it is
# encoded in line color in the full-sample plot.
DEFAULT_AXES = ["tan_beta", "lambda6", "m_phi", "lam1", "total_width", "br_gaga"]


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

    raise ValueError(f"Could not resolve column '{requested}' in schema.")


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [item.strip() for item in raw.split(",") if item.strip()]


def transform_for_plot(values: np.ndarray, use_log: bool) -> np.ndarray:
    out = values.astype(np.float64, copy=True)
    if use_log:
        # For strictly positive physics quantities, log10 helps compress decades.
        out[out <= 0] = np.nan
        out = np.log10(out)
    return out


def format_limit(v: float) -> str:
    av = abs(v)
    if av == 0:
        return "0"
    if av >= 1e4 or av < 1e-3:
        return f"{v:.2e}"
    if av >= 100:
        return f"{v:.2f}"
    return f"{v:.4g}"


def build_axis_projection(values_2d: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mins = np.nanmin(values_2d, axis=0)
    maxs = np.nanmax(values_2d, axis=0)
    spans = maxs - mins

    projected = np.empty_like(values_2d, dtype=np.float64)
    for idx in range(values_2d.shape[1]):
        if np.isclose(spans[idx], 0.0):
            projected[:, idx] = 0.5
        else:
            projected[:, idx] = (values_2d[:, idx] - mins[idx]) / spans[idx]

    return projected, mins, maxs


def add_lines_chunked(
    ax: plt.Axes,
    y_values: np.ndarray,
    x_positions: np.ndarray,
    colors: np.ndarray | str,
    alpha: float,
    linewidth: float,
    chunk_size: int,
) -> int:
    if y_values.shape[0] == 0:
        return 0

    added = 0
    for start in range(0, y_values.shape[0], chunk_size):
        stop = min(start + chunk_size, y_values.shape[0])
        chunk = y_values[start:stop]
        segments = np.stack(
            [
                np.broadcast_to(x_positions, chunk.shape),
                chunk,
            ],
            axis=2,
        )
        if isinstance(colors, str):
            chunk_colors: np.ndarray | str = colors
        else:
            chunk_colors = colors[start:stop]
        lc = LineCollection(segments, colors=chunk_colors, linewidths=linewidth, alpha=alpha)
        ax.add_collection(lc)
        added += chunk.shape[0]

    return added


def draw_axis_limits(
    ax: plt.Axes,
    x_positions: np.ndarray,
    axis_labels: list[str],
    mins: np.ndarray,
    maxs: np.ndarray,
    mode_label: str,
) -> None:
    for idx, x in enumerate(x_positions):
        ax.vlines(x, 0.0, 1.0, color="black", linewidth=0.8, alpha=0.35)
        ax.text(x, -0.07, format_limit(mins[idx]), ha="center", va="top", fontsize=8)
        ax.text(x, 1.03, format_limit(maxs[idx]), ha="center", va="bottom", fontsize=8)
        if mode_label == "log":
            ax.text(x, 1.11, "log10", ha="center", va="bottom", fontsize=7, color="#555555")

    ax.set_xticks(x_positions)
    ax.set_xticklabels(axis_labels, rotation=20, ha="right")


def compute_rgba_by_value(values: np.ndarray, cmap_name: str, alpha: float) -> np.ndarray:
    vmin = float(np.nanmin(values))
    vmax = float(np.nanmax(values))
    if np.isclose(vmin, vmax):
        normed = np.full(values.shape, 0.5, dtype=np.float64)
    else:
        normed = (values - vmin) / (vmax - vmin)
    cmap = colormaps.get_cmap(cmap_name)
    rgba = cmap(normed)
    rgba[:, 3] = alpha
    return rgba


def resolve_axes(axes_csv: str, exclude_csv: str) -> list[str]:
    requested = parse_csv_list(axes_csv)
    exclude = set(parse_csv_list(exclude_csv))
    if not requested:
        requested = list(DEFAULT_AXES)

    final_axes = [c for c in requested if c in ALL_COLUMNS and c not in exclude]
    if len(final_axes) < 2:
        raise ValueError("At least two plotting axes are required after applying --axes/--exclude-axes.")
    return final_axes


def main() -> None:
    parser = argparse.ArgumentParser(description="Parallel coordinates for diHiggs parameter space with ROI highlighting.")
    parser.add_argument("--input", type=str, default=str(DEFAULT_PARQUET), help="Path to parquet input")
    parser.add_argument("--output-dir", type=str, default=str(OUTPUT_DIR), help="Directory for figures and summary")
    parser.add_argument(
        "--axes",
        type=str,
        default=",".join(DEFAULT_AXES),
        help="Comma-separated axes to plot in order. Example: tan_beta,lambda6,m_phi,lam1,total_width,br_gaga",
    )
    parser.add_argument(
        "--exclude-axes",
        type=str,
        default="",
        help="Comma-separated axes to remove after --axes selection.",
    )
    parser.add_argument(
        "--color-by",
        type=str,
        default="ctau",
        help="Column used to color full-sample trajectories.",
    )
    parser.add_argument(
        "--plot-modes",
        type=str,
        choices=["raw", "log", "both"],
        default="both",
        help="Generate unscaled plots in raw values, log10 values, or both.",
    )
    parser.add_argument(
        "--log-cols",
        type=str,
        default="total_width,ctau,br_gaga",
        help="Comma-separated columns transformed with log10 when --plot-modes includes log.",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Optional cap on rows to render (head after filtering).",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=50000,
        help="Chunk size for drawing line collections.",
    )
    parser.add_argument(
        "--background-alpha",
        type=float,
        default=0.03,
        help="Alpha for non-ROI trajectories.",
    )
    parser.add_argument(
        "--roi-alpha",
        type=float,
        default=0.8,
        help="Alpha for ROI trajectories.",
    )
    parser.add_argument(
        "--line-width",
        type=float,
        default=0.25,
        help="Line width for trajectories.",
    )
    parser.add_argument(
        "--full-alpha",
        type=float,
        default=0.08,
        help="Line alpha for full-sample colored plot.",
    )
    parser.add_argument(
        "--full-cmap",
        type=str,
        default="viridis",
        help="Colormap for full-sample color encoding.",
    )
    parser.add_argument(
        "--ctau-threshold",
        type=float,
        default=None,
        help="ROI condition: ctau < threshold. If omitted, uses ctau quantile.",
    )
    parser.add_argument(
        "--br-threshold",
        type=float,
        default=None,
        help="ROI condition: br_gaga > threshold. If omitted, uses BR quantile.",
    )
    parser.add_argument(
        "--ctau-quantile",
        type=float,
        default=0.2,
        help="Quantile used for ctau threshold when --ctau-threshold is omitted.",
    )
    parser.add_argument(
        "--br-quantile",
        type=float,
        default=0.8,
        help="Quantile used for BR threshold when --br-threshold is omitted.",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input parquet not found: {input_path}")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    schema = pl.scan_parquet(input_path).collect_schema()
    schema_names = schema.names()

    resolved = {
        "tan_beta": resolve_column("tan_beta", schema_names),
        "lambda6": resolve_column("lambda6", schema_names),
        "m_phi": resolve_column("m_phi", schema_names, ["mphi"]),
        "lam1": resolve_column("lam1", schema_names, ["lambda1"]),
        "total_width": resolve_column("total_width", schema_names, ["total_decay_width", "w_total"]),
        "br_gaga": resolve_column("br_gaga", schema_names, ["branching_ratio", "br_gg"]),
    }

    if args.color_by not in ALL_COLUMNS:
        raise ValueError(f"--color-by must be one of: {ALL_COLUMNS}")

    axes_cols = resolve_axes(args.axes, args.exclude_axes)

    selected_cols = [
        resolved["tan_beta"],
        resolved["lambda6"],
        resolved["m_phi"],
        resolved["lam1"],
        resolved["total_width"],
        resolved["br_gaga"],
    ]

    logger.info("Reading required columns lazily from %s", input_path)
    lf = pl.scan_parquet(input_path).select(selected_cols)

    lf = lf.with_columns(
        pl.when(pl.col(resolved["total_width"]).is_not_null() & (~pl.col(resolved["total_width"]).is_nan()) & (pl.col(resolved["total_width"]) > 0))
        .then((6.582119569e-25 * 2.99792458e11) / pl.col(resolved["total_width"]).cast(pl.Float64))
        .otherwise(None)
        .alias("ctau")
    )

    lf = lf.filter(
        pl.all_horizontal(
            [
                pl.col(resolved["tan_beta"]).is_not_null(),
                pl.col(resolved["lambda6"]).is_not_null(),
                pl.col(resolved["m_phi"]).is_not_null(),
                pl.col(resolved["lam1"]).is_not_null(),
                pl.col(resolved["total_width"]).is_not_null(),
                pl.col("ctau").is_not_null(),
                pl.col(resolved["br_gaga"]).is_not_null(),
            ]
        )
    )

    if args.max_rows is not None and args.max_rows > 0:
        lf = lf.limit(args.max_rows)

    df = _stream_collect(
        lf.select(
            [
                pl.col(resolved["tan_beta"]).cast(pl.Float64).alias("tan_beta"),
                pl.col(resolved["lambda6"]).cast(pl.Float64).alias("lambda6"),
                pl.col(resolved["m_phi"]).cast(pl.Float64).alias("m_phi"),
                pl.col(resolved["lam1"]).cast(pl.Float64).alias("lam1"),
                pl.col(resolved["total_width"]).cast(pl.Float64).alias("total_width"),
                pl.col("ctau").cast(pl.Float64).alias("ctau"),
                pl.col(resolved["br_gaga"]).cast(pl.Float64).alias("br_gaga"),
            ]
        )
    )

    data = df.to_numpy()
    col_index = {name: idx for idx, name in enumerate(ALL_COLUMNS)}

    finite_rows = np.isfinite(data).all(axis=1)
    data = data[finite_rows]

    if data.shape[0] == 0:
        raise RuntimeError("No valid rows remain after filtering finite values.")

    ctau_vals = data[:, col_index["ctau"]]
    br_vals = data[:, col_index["br_gaga"]]

    ctau_threshold = args.ctau_threshold
    if ctau_threshold is None:
        ctau_threshold = float(np.quantile(ctau_vals, args.ctau_quantile))

    br_threshold = args.br_threshold
    if br_threshold is None:
        br_threshold = float(np.quantile(br_vals, args.br_quantile))

    roi_mask = (ctau_vals < ctau_threshold) & (br_vals > br_threshold)

    log_cols = set(parse_csv_list(args.log_cols))
    plot_modes = ["raw", "log"] if args.plot_modes == "both" else [args.plot_modes]
    saved_files: dict[str, dict[str, str]] = {}

    for mode in plot_modes:
        use_log_mode = mode == "log"

        transformed = np.zeros_like(data, dtype=np.float64)
        for idx, col in enumerate(ALL_COLUMNS):
            use_log_col = use_log_mode and (col in log_cols)
            transformed[:, idx] = transform_for_plot(data[:, idx], use_log=use_log_col)

        axis_idx = [col_index[c] for c in axes_cols]
        axis_values = transformed[:, axis_idx]

        color_idx = col_index[args.color_by]
        color_values = transformed[:, color_idx]

        valid = np.isfinite(axis_values).all(axis=1) & np.isfinite(color_values)
        axis_values = axis_values[valid]
        color_values = color_values[valid]
        roi_mode_mask = roi_mask[valid]

        if axis_values.shape[0] == 0:
            logger.warning("Skipping mode '%s': no finite rows after transform/filter.", mode)
            continue

        projected, mins, maxs = build_axis_projection(axis_values)
        x_positions = np.arange(len(axes_cols), dtype=np.float64)

        rgba = compute_rgba_by_value(color_values, cmap_name=args.full_cmap, alpha=args.full_alpha)

        fig_all, ax_all = plt.subplots(figsize=(13, 7))
        n_all = add_lines_chunked(
            ax=ax_all,
            y_values=projected,
            x_positions=x_positions,
            colors=rgba,
            alpha=1.0,
            linewidth=args.line_width,
            chunk_size=args.chunk_size,
        )

        ax_all.set_xlim(x_positions[0], x_positions[-1])
        ax_all.set_ylim(0.0, 1.0)
        ax_all.set_ylabel("Per-axis value range")
        ax_all.set_title(f"Parallel Coordinates (full): color by {args.color_by} [{mode}]")
        ax_all.grid(axis="x", linestyle="--", alpha=0.25)
        draw_axis_limits(ax_all, x_positions, axes_cols, mins, maxs, mode)

        sm = plt.cm.ScalarMappable(cmap=args.full_cmap)
        sm.set_clim(vmin=float(np.nanmin(color_values)), vmax=float(np.nanmax(color_values)))
        cbar = fig_all.colorbar(sm, ax=ax_all, pad=0.02)
        cbar.set_label(f"{args.color_by} ({mode})")

        fig_all.tight_layout()

        full_png = out_dir / f"parallel_full_colorby_{args.color_by}_{mode}.png"
        full_pdf = out_dir / f"parallel_full_colorby_{args.color_by}_{mode}.pdf"
        fig_all.savefig(full_png, dpi=300, bbox_inches="tight")
        fig_all.savefig(full_pdf, bbox_inches="tight")
        plt.close(fig_all)

        bg = projected[~roi_mode_mask]
        fg = projected[roi_mode_mask]

        fig_roi, ax_roi = plt.subplots(figsize=(13, 7))
        n_bg = add_lines_chunked(
            ax=ax_roi,
            y_values=bg,
            x_positions=x_positions,
            colors="#8a8a8a",
            alpha=args.background_alpha,
            linewidth=args.line_width,
            chunk_size=args.chunk_size,
        )
        n_fg = add_lines_chunked(
            ax=ax_roi,
            y_values=fg,
            x_positions=x_positions,
            colors="#d62728",
            alpha=args.roi_alpha,
            linewidth=max(args.line_width, 0.4),
            chunk_size=args.chunk_size,
        )

        ax_roi.set_xlim(x_positions[0], x_positions[-1])
        ax_roi.set_ylim(0.0, 1.0)
        ax_roi.set_ylabel("Per-axis value range")
        ax_roi.set_title(
            f"Parallel Coordinates ROI [{mode}]: ctau < {ctau_threshold:.3e} and br_gaga > {br_threshold:.3e}"
        )
        ax_roi.grid(axis="x", linestyle="--", alpha=0.25)
        draw_axis_limits(ax_roi, x_positions, axes_cols, mins, maxs, mode)
        fig_roi.tight_layout()

        roi_png = out_dir / f"parallel_roi_{mode}.png"
        roi_pdf = out_dir / f"parallel_roi_{mode}.pdf"
        fig_roi.savefig(roi_png, dpi=300, bbox_inches="tight")
        fig_roi.savefig(roi_pdf, bbox_inches="tight")
        plt.close(fig_roi)

        saved_files[mode] = {
            "full_png": str(full_png),
            "full_pdf": str(full_pdf),
            "roi_png": str(roi_png),
            "roi_pdf": str(roi_pdf),
            "rows_full_plot": int(n_all),
            "rows_background": int(n_bg),
            "rows_roi": int(n_fg),
        }

    summary = {
        "input_file": str(input_path),
        "output_dir": str(out_dir),
        "resolved_columns": resolved,
        "available_columns": ALL_COLUMNS,
        "axes_used": axes_cols,
        "color_by": args.color_by,
        "plot_modes": plot_modes,
        "log_transformed_columns": sorted(log_cols),
        "rows_total_after_numeric_filter": int(data.shape[0]),
        "outputs": saved_files,
        "ctau_threshold": float(ctau_threshold),
        "br_threshold": float(br_threshold),
        "ctau_quantile": float(args.ctau_quantile),
        "br_quantile": float(args.br_quantile),
        "background_alpha": float(args.background_alpha),
        "roi_alpha": float(args.roi_alpha),
        "line_width": float(args.line_width),
    }

    summary_path = out_dir / "parallel_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    logger.info("Saved plot modes: %s", sorted(saved_files.keys()))
    logger.info("Summary: %s", summary_path)


if __name__ == "__main__":
    main()
