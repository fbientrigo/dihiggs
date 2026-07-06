#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""parallel_coordinates_population.py
=====================================
Parallel-coordinates plots over the small/big ctau population sample
produced by ``extract_ctau_population.py`` (``population_for_parcoords.parquet``).

Two independent color modes:
* categorical -- small (blue) vs big (red) population, 2-color overlay.
* gradient -- all rows, continuous colormap by ctau_mm value.

Standalone script; duplicates small helpers from ``parallel_coordinates_roi.py``
rather than importing it (no shared-utils module in this repo, and that
script's ROI mask/numeric-only pipeline don't fit a categorical class column).
"""
from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib import colormaps
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

plt.style.use("default")
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 12,
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 10,
        "figure.figsize": (13, 7),
        "axes.linewidth": 1.2,
    }
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")

DEFAULT_AXES = ["mA_target", "tan_beta", "lambda6", "m_phi", "lam1", "total_width", "br_gaga"]


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
    raise ValueError(f"Could not resolve column '{requested}' in schema. Available: {schema_names}")


def parse_csv_list(raw: str | None) -> list[str]:
    if not raw:
        return []
    return [item.strip() for item in raw.split(",") if item.strip()]


def transform_for_plot(values: np.ndarray, use_log: bool) -> np.ndarray:
    out = values.astype(np.float64, copy=True)
    if use_log:
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
            [np.broadcast_to(x_positions, chunk.shape), chunk],
            axis=2,
        )
        chunk_colors = colors if isinstance(colors, str) else colors[start:stop]
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


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Path to population_for_parcoords.parquet")
    p.add_argument("--output-dir", required=True, help="Directory for figures/summary")
    p.add_argument("--axes", type=str, default=",".join(DEFAULT_AXES))
    p.add_argument("--color-mode", choices=["categorical", "gradient", "both"], default="both")
    p.add_argument("--color-by", type=str, default="ctau_mm")
    p.add_argument("--plot-modes", choices=["raw", "log", "both"], default="both")
    p.add_argument("--log-cols", type=str, default="total_width,ctau_mm,br_gaga")
    p.add_argument("--small-color", type=str, default="#4477aa")
    p.add_argument("--big-color", type=str, default="#cc3311")
    p.add_argument("--small-alpha", type=float, default=0.05)
    p.add_argument("--big-alpha", type=float, default=0.5)
    p.add_argument("--full-cmap", type=str, default="viridis")
    p.add_argument("--chunk-size", type=int, default=50000)
    p.add_argument("--line-width", type=float, default=0.3)
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    input_path = Path(args.input)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    schema_names = pl.scan_parquet(input_path).collect_schema().names()
    axes_cols = parse_csv_list(args.axes)
    if len(axes_cols) < 2:
        raise ValueError("At least two --axes are required.")
    for c in axes_cols:
        resolve_column(c, schema_names)
    color_col = resolve_column(args.color_by, schema_names)

    needed_cols = sorted(set(axes_cols + [color_col, "population_class"]))
    df = _stream_collect(pl.scan_parquet(input_path).select(needed_cols))

    class_labels = df["population_class"].to_numpy()
    numeric_cols = [c for c in needed_cols if c != "population_class"]
    data = df.select(numeric_cols).to_numpy()
    col_index = {name: idx for idx, name in enumerate(numeric_cols)}

    finite_rows = np.isfinite(data).all(axis=1)
    data = data[finite_rows]
    class_labels = class_labels[finite_rows]
    if data.shape[0] == 0:
        raise RuntimeError("No finite rows remain after filtering.")

    log_cols = set(parse_csv_list(args.log_cols))
    plot_modes = ["raw", "log"] if args.plot_modes == "both" else [args.plot_modes]
    color_modes = ["categorical", "gradient"] if args.color_mode == "both" else [args.color_mode]

    saved_files: dict[str, dict] = {}

    for mode in plot_modes:
        use_log_mode = mode == "log"
        transformed = np.zeros_like(data, dtype=np.float64)
        for idx, col in enumerate(numeric_cols):
            use_log_col = use_log_mode and (col in log_cols)
            transformed[:, idx] = transform_for_plot(data[:, idx], use_log=use_log_col)

        axis_idx = [col_index[c] for c in axes_cols]
        axis_values = transformed[:, axis_idx]
        color_idx = col_index[color_col]
        color_values = transformed[:, color_idx]

        valid = np.isfinite(axis_values).all(axis=1) & np.isfinite(color_values)
        axis_values_v = axis_values[valid]
        color_values_v = color_values[valid]
        class_labels_v = class_labels[valid]

        if axis_values_v.shape[0] == 0:
            logger.warning("Skipping mode '%s': no finite rows after transform/filter.", mode)
            continue

        projected, mins, maxs = build_axis_projection(axis_values_v)
        x_positions = np.arange(len(axes_cols), dtype=np.float64)

        if "categorical" in color_modes:
            mask_big = class_labels_v == "big"
            fig, ax = plt.subplots(figsize=(13, 7))
            n_small_drawn = add_lines_chunked(
                ax, projected[~mask_big], x_positions, args.small_color,
                args.small_alpha, args.line_width, args.chunk_size,
            )
            n_big_drawn = add_lines_chunked(
                ax, projected[mask_big], x_positions, args.big_color,
                args.big_alpha, max(args.line_width, 0.4), args.chunk_size,
            )
            ax.set_xlim(x_positions[0], x_positions[-1])
            ax.set_ylim(0.0, 1.0)
            ax.set_ylabel("Per-axis value range")
            ax.set_title(f"Parallel Coordinates [{mode}]: small vs big ctau population")
            ax.grid(axis="x", linestyle="--", alpha=0.25)
            draw_axis_limits(ax, x_positions, axes_cols, mins, maxs, mode)
            legend_handles = [
                Line2D([0], [0], color=args.small_color, lw=2, label=f"small (n={int((~mask_big).sum())})"),
                Line2D([0], [0], color=args.big_color, lw=2, label=f"big (n={int(mask_big.sum())})"),
            ]
            ax.legend(handles=legend_handles, loc="upper right")
            fig.tight_layout()
            cat_png = out_dir / f"parallel_categorical_{mode}.png"
            cat_pdf = out_dir / f"parallel_categorical_{mode}.pdf"
            fig.savefig(cat_png, dpi=300, bbox_inches="tight")
            fig.savefig(cat_pdf, bbox_inches="tight")
            plt.close(fig)
            saved_files.setdefault(mode, {})["categorical"] = {
                "png": str(cat_png), "pdf": str(cat_pdf),
                "n_small": int(n_small_drawn), "n_big": int(n_big_drawn),
            }

        if "gradient" in color_modes:
            rgba = compute_rgba_by_value(color_values_v, cmap_name=args.full_cmap, alpha=0.15)
            fig2, ax2 = plt.subplots(figsize=(13, 7))
            n_drawn = add_lines_chunked(
                ax2, projected, x_positions, rgba, 1.0, args.line_width, args.chunk_size,
            )
            ax2.set_xlim(x_positions[0], x_positions[-1])
            ax2.set_ylim(0.0, 1.0)
            ax2.set_ylabel("Per-axis value range")
            ax2.set_title(f"Parallel Coordinates [{mode}]: color by {color_col}")
            ax2.grid(axis="x", linestyle="--", alpha=0.25)
            draw_axis_limits(ax2, x_positions, axes_cols, mins, maxs, mode)
            sm = plt.cm.ScalarMappable(cmap=args.full_cmap)
            sm.set_clim(vmin=float(np.nanmin(color_values_v)), vmax=float(np.nanmax(color_values_v)))
            cbar = fig2.colorbar(sm, ax=ax2, pad=0.02)
            cbar.set_label(f"{color_col} ({mode})")
            fig2.tight_layout()
            grad_png = out_dir / f"parallel_gradient_colorby_{color_col}_{mode}.png"
            grad_pdf = out_dir / f"parallel_gradient_colorby_{color_col}_{mode}.pdf"
            fig2.savefig(grad_png, dpi=300, bbox_inches="tight")
            fig2.savefig(grad_pdf, bbox_inches="tight")
            plt.close(fig2)
            saved_files.setdefault(mode, {})["gradient"] = {
                "png": str(grad_png), "pdf": str(grad_pdf), "n_rows": int(n_drawn),
            }

    summary = {
        "input_file": str(input_path),
        "output_dir": str(out_dir),
        "axes_used": axes_cols,
        "color_by": color_col,
        "color_modes": color_modes,
        "plot_modes": plot_modes,
        "log_transformed_columns": sorted(log_cols),
        "rows_total_after_finite_filter": int(data.shape[0]),
        "outputs": saved_files,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    summary_path = out_dir / "parallel_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write("\n")

    logger.info("Saved color modes: %s, plot modes: %s", color_modes, plot_modes)
    logger.info("Summary: %s", summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
