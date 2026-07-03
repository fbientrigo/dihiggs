#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
plot_three_brs_fixed.py
=======================
Plot 3 branching ratios (BR_Zga, BR_gaga, BR_bb) in the same figure with different colors.
Fixed-cut mode for specific parameter values.
"""

from __future__ import annotations

import argparse
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
        "figure.figsize": (10, 7),
        "axes.linewidth": 1.5,
        "xtick.major.width": 1.5,
        "ytick.major.width": 1.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    }
)

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


def plot_three_brs(
    df: pd.DataFrame,
    sin_ba: float,
    tan_beta: float,
    mA: float,
    lambda6: float,
    lambda7: float,
    point_size: float = 3.0,
    point_alpha: float = 0.6,
) -> tuple[Path, Path]:
    """
    Plot m_phi vs three branching ratios in the same figure.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with columns ['m_phi', 'width_bb', 'width_Zga', 'width_gaga', 'total_width']
    sin_ba, tan_beta, mA, lambda6, lambda7 : float
        Fixed parameter values for the title
    point_size : float
        Scatter marker size
    point_alpha : float
        Scatter marker transparency
    
    Returns:
    --------
    tuple[Path, Path]
        Paths to saved PNG and PDF files
    """
    
    # Calculate branching ratios
    df = df.copy()
    df['BR_bb'] = df['width_bb'] / df['total_width']
    df['BR_Zga'] = df['width_Zga'] / df['total_width']
    df['BR_gaga'] = df['width_gaga'] / df['total_width']
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Define colors for each BR
    colors = {
        'BR_Zga': '#1f77b4',  # Blue
        'BR_gaga': '#ff7f0e',  # Orange
        'BR_bb': '#2ca02c',    # Green
    }
    
    labels = {
        'BR_Zga': r'$BR(\phi \to Z\gamma)$',
        'BR_gaga': r'$BR(\phi \to \gamma\gamma)$',
        'BR_bb': r'$BR(\phi \to b\bar{b})$',
    }
    
    # Plot each BR
    for br_name, color in colors.items():
        x = df['m_phi'].to_numpy()
        y = df[br_name].to_numpy()
        valid = np.isfinite(x) & np.isfinite(y) & (y > 0)
        
        if not valid.any():
            logger.warning(f"No valid points for {br_name}")
            continue
        
        ax.scatter(
            x[valid],
            y[valid],
            color=color,
            alpha=point_alpha,
            s=point_size,
            edgecolors='none',
            label=labels[br_name],
            zorder=2,
        )
    
    ax.set_yscale('log')
    ax.set_xlabel(r'$m_\phi$ [GeV]')
    ax.set_ylabel('Branching Ratio')
    
    title = (
        rf"$\sin(\beta-\alpha)={sin_ba:g}$, "
        rf"$\tan\beta={tan_beta:g}$, "
        rf"$m_A={mA:g}$ GeV, "
        rf"$\lambda_6={lambda6:g}$, "
        rf"$\lambda_7={lambda7:g}$"
    )
    ax.set_title(title)
    ax.grid(True, linestyle='--', alpha=0.4, which='both')
    ax.legend(loc='best', framealpha=0.95)
    fig.tight_layout()
    
    # Generate output filename
    output_stem = (
        f"three_brs_sinba{sin_ba:g}_tb{tan_beta:g}_"
        f"mA{mA:g}_l6{lambda6:g}_l7{lambda7:g}"
    ).replace(".", "p")
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    png_path = OUTPUT_DIR / f"{output_stem}.png"
    pdf_path = OUTPUT_DIR / f"{output_stem}.pdf"
    
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    logger.info(f"Saved PNG: {png_path}")
    
    fig.savefig(pdf_path, bbox_inches='tight')
    logger.info(f"Saved PDF: {pdf_path}")
    
    plt.close(fig)
    
    return png_path, pdf_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot 3 branching ratios in the same figure for fixed parameters."
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to parquet file"
    )
    parser.add_argument("--sin-ba", type=float, default=1.0, help="Fixed sin_ba cut")
    parser.add_argument("--tan-beta", type=float, required=True, help="Fixed tan_beta cut")
    parser.add_argument("--mA", type=float, required=True, help="Fixed mA cut")
    parser.add_argument("--lambda6", type=float, required=True, help="Fixed lambda6 cut")
    parser.add_argument("--lambda7", type=float, required=True, help="Fixed lambda7 cut")
    parser.add_argument("--point-size", type=float, default=3.0, help="Scatter marker size")
    parser.add_argument("--point-alpha", type=float, default=0.6, help="Scatter marker alpha")
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    logger.info(f"Loading data from {input_path}")
    
    # Define required columns
    required = [
        "m_phi",
        "width_bb",
        "width_Zga",
        "width_gaga",
        "total_width",
        "sin_ba",
        "tan_beta",
        "mA",
        "lambda6",
        "lambda7",
    ]
    
    # Load and filter data
    lf = pl.scan_parquet(input_path).select(required)
    df = _stream_collect(
        lf.filter(
            (pl.col("sin_ba") == float(args.sin_ba))
            & (pl.col("tan_beta") == float(args.tan_beta))
            & (pl.col("mA") == float(args.mA))
            & (pl.col("lambda6") == float(args.lambda6))
            & (pl.col("lambda7") == float(args.lambda7))
        )
    ).to_pandas()
    
    if df.empty:
        raise RuntimeError(
            f"No rows found for:\n"
            f"  sin_ba = {args.sin_ba}\n"
            f"  tan_beta = {args.tan_beta}\n"
            f"  mA = {args.mA}\n"
            f"  lambda6 = {args.lambda6}\n"
            f"  lambda7 = {args.lambda7}"
        )
    
    logger.info(f"Found {len(df)} rows with the specified parameters")
    
    # Plot
    png_path, pdf_path = plot_three_brs(
        df,
        args.sin_ba,
        args.tan_beta,
        args.mA,
        args.lambda6,
        args.lambda7,
        point_size=args.point_size,
        point_alpha=args.point_alpha,
    )
    
    logger.info(f"Done! Output files:")
    logger.info(f"  PNG: {png_path}")
    logger.info(f"  PDF: {pdf_path}")


if __name__ == "__main__":
    main()
