#!/usr/bin/env python3
"""Histograms of the ctau distribution, small vs big population.

Computes log-spaced bin counts directly via Polars lazy aggregation on
``silver_all.parquet`` -- no per-row materialization, so this is safe and
cheap regardless of the 6,000,000-row input size. Produces two figures:

* categorical: small (ctau < boundary) vs big (ctau >= boundary), 2 colors.
* gradient: single histogram, bars colored by bin-center ctau value.
"""
from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib import colormaps

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
        "figure.figsize": (10, 6),
        "axes.linewidth": 1.2,
    }
)

# shared constant (was 1.973269804e-13); see physics_conventions.py
from physics_conventions import HBAR_C_GEV_MM as HBARC_GEV_MM


def _stream_collect(lf: pl.LazyFrame) -> pl.DataFrame:
    try:
        return lf.collect(engine="streaming")
    except TypeError:
        return lf.collect(streaming=True)
    except Exception:
        return lf.collect()


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Path to silver_all.parquet")
    p.add_argument("--output-dir", required=True, help="Directory for figures/summary")
    p.add_argument("--ctau-min", type=float, default=1e-11)
    p.add_argument("--ctau-max", type=float, default=300.0)
    p.add_argument("--n-bins", type=int, default=120)
    p.add_argument("--boundary", type=float, default=1.0,
                   help="ctau_mm threshold defining the 'big' population")
    p.add_argument("--cmap", type=str, default="viridis")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log_min = math.log10(args.ctau_min)
    log_max = math.log10(args.ctau_max)
    bin_width = (log_max - log_min) / args.n_bins
    edges = np.logspace(log_min, log_max, args.n_bins + 1)
    centers = np.sqrt(edges[:-1] * edges[1:])

    lf = pl.scan_parquet(input_path).select(["total_width"]).with_columns(
        (HBARC_GEV_MM / pl.col("total_width").cast(pl.Float64)).alias("ctau_mm")
    ).with_columns(
        pl.col("ctau_mm").log10().alias("log_ctau")
    ).with_columns(
        ((pl.col("log_ctau") - log_min) / bin_width)
        .floor()
        .cast(pl.Int64)
        .clip(0, args.n_bins - 1)
        .alias("bin_idx")
    )

    binned = _stream_collect(
        lf.group_by("bin_idx").agg(pl.len().alias("n"))
    ).sort("bin_idx")
    counts = np.zeros(args.n_bins, dtype=np.int64)
    for row in binned.iter_rows(named=True):
        counts[int(row["bin_idx"])] = row["n"]

    total_rows = int(counts.sum())

    # exact cross-check totals (not binned)
    exact = _stream_collect(
        lf.select(
            (pl.col("ctau_mm") >= args.boundary).sum().alias("n_big"),
            (pl.col("ctau_mm") < args.boundary).sum().alias("n_small"),
        )
    )
    n_big_exact = int(exact["n_big"].item())
    n_small_exact = int(exact["n_small"].item())

    is_big_bin = centers >= args.boundary

    # --- categorical histogram ---
    fig, ax = plt.subplots()
    ax.bar(centers[~is_big_bin], counts[~is_big_bin],
           width=np.diff(edges)[~is_big_bin], color="#4477aa", alpha=0.85,
           label=f"small (ctau < {args.boundary:g} mm)", align="center")
    ax.bar(centers[is_big_bin], counts[is_big_bin],
           width=np.diff(edges)[is_big_bin], color="#cc3311", alpha=0.85,
           label=f"big (ctau >= {args.boundary:g} mm)", align="center")
    ax.axvline(args.boundary, color="black", linestyle="--", linewidth=1.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$c\tau$ [mm]")
    ax.set_ylabel("Counts")
    ax.set_title("ctau distribution: small vs big population")
    ax.grid(True, which="both", linestyle="--", alpha=0.3)
    ax.legend(loc="upper right")
    fig.tight_layout()
    cat_png = output_dir / "ctau_hist_categorical.png"
    cat_pdf = output_dir / "ctau_hist_categorical.pdf"
    fig.savefig(cat_png, dpi=300, bbox_inches="tight")
    fig.savefig(cat_pdf, bbox_inches="tight")
    plt.close(fig)

    # --- gradient histogram ---
    cmap = colormaps.get_cmap(args.cmap)
    log_centers = np.log10(centers)
    vmin, vmax = log_centers.min(), log_centers.max()
    normed = (log_centers - vmin) / (vmax - vmin)
    bar_colors = cmap(normed)

    fig2, ax2 = plt.subplots()
    ax2.bar(centers, counts, width=np.diff(edges), color=bar_colors, align="center")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel(r"$c\tau$ [mm]")
    ax2.set_ylabel("Counts")
    ax2.set_title("ctau distribution: colored by ctau value")
    ax2.grid(True, which="both", linestyle="--", alpha=0.3)
    sm = plt.cm.ScalarMappable(cmap=args.cmap)
    sm.set_clim(vmin=float(centers.min()), vmax=float(centers.max()))
    cbar = fig2.colorbar(sm, ax=ax2, pad=0.02)
    cbar.set_label(r"$c\tau$ [mm]")
    fig2.tight_layout()
    grad_png = output_dir / "ctau_hist_gradient.png"
    grad_pdf = output_dir / "ctau_hist_gradient.pdf"
    fig2.savefig(grad_png, dpi=300, bbox_inches="tight")
    fig2.savefig(grad_pdf, bbox_inches="tight")
    plt.close(fig2)

    summary = {
        "input_parquet": str(input_path),
        "output_dir": str(output_dir),
        "ctau_min": args.ctau_min,
        "ctau_max": args.ctau_max,
        "n_bins": args.n_bins,
        "boundary": args.boundary,
        "bin_edges": edges.tolist(),
        "bin_centers": centers.tolist(),
        "bin_counts": counts.tolist(),
        "total_rows_binned": total_rows,
        "n_big_exact": n_big_exact,
        "n_small_exact": n_small_exact,
        "categorical_png": str(cat_png),
        "gradient_png": str(grad_png),
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    summary_path = output_dir / "ctau_hist_summary.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write("\n")

    print(f"[hist] total_rows_binned={total_rows} n_big_exact={n_big_exact} n_small_exact={n_small_exact}")
    print(f"[hist] wrote {cat_png}")
    print(f"[hist] wrote {grad_png}")
    print(f"[hist] wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
