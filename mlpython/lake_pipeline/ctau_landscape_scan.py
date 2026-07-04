#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl

HBAR_C_GEV_MM = 1.973269804e-13
ATOL = {"mA": 1e-8, "tan_beta": 1e-3, "lambda6": 1e-7, "lambda7": 1e-9, "sin_ba": 1e-9}

STYLE_PRESETS = {
    "paper": {"figsize": (8.6, 5.6), "point_size": 12, "alpha": 0.7, "grid_alpha": 0.25},
    "talk": {"figsize": (10.5, 6.8), "point_size": 22, "alpha": 0.8, "grid_alpha": 0.3},
    "draft": {"figsize": (9.2, 5.8), "point_size": 14, "alpha": 0.72, "grid_alpha": 0.3},
    "compact": {"figsize": (7.2, 4.6), "point_size": 10, "alpha": 0.65, "grid_alpha": 0.22},
}

logger = logging.getLogger("ctau_landscape")


def _stream_collect(lf: pl.LazyFrame) -> pl.DataFrame:
    try:
        return lf.collect(engine="streaming")
    except TypeError:
        return lf.collect(streaming=True)
    except Exception:
        return lf.collect()


def _parse_float_list(raw: str | None) -> list[float]:
    if not raw:
        return []
    out = []
    for x in raw.split(","):
        x = x.strip()
        if not x:
            continue
        v = float(x)
        if np.isfinite(v):
            out.append(float(v))
    return out


def _resolve_column(schema_names: list[str], requested: str, aliases: list[str] | None = None) -> str:
    if requested in schema_names:
        return requested
    lower_map = {c.lower(): c for c in schema_names}
    if requested.lower() in lower_map:
        return lower_map[requested.lower()]
    for a in aliases or []:
        if a.lower() in lower_map:
            return lower_map[a.lower()]
    raise ValueError(f"Could not resolve column '{requested}' from schema.")


def _float_tol_cut(col: str, target: float, atol: float) -> pl.Expr:
    return (pl.col(col).cast(pl.Float64) - pl.lit(float(target))).abs() <= float(atol)


def _apply_value_filter(lf: pl.LazyFrame, col: str, values: list[float], atol: float) -> pl.LazyFrame:
    if not values:
        return lf
    expr = None
    for v in values:
        e = _float_tol_cut(col, v, atol)
        expr = e if expr is None else (expr | e)
    return lf.filter(expr)


def _quantile_or_none(arr: np.ndarray, q: float) -> float | None:
    vals = arr[np.isfinite(arr)]
    if vals.size == 0:
        return None
    return float(np.quantile(vals, q))


def _safe_min(arr: np.ndarray) -> float | None:
    vals = arr[np.isfinite(arr)]
    return float(np.min(vals)) if vals.size else None


def _safe_max(arr: np.ndarray) -> float | None:
    vals = arr[np.isfinite(arr)]
    return float(np.max(vals)) if vals.size else None


def _first_existing(cols: list[str], candidates: list[str]) -> str | None:
    lower = {c.lower(): c for c in cols}
    for c in candidates:
        if c in cols:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def _subsample(df: pd.DataFrame, nmax: int, seed: int) -> pd.DataFrame:
    if len(df) <= nmax:
        return df
    return df.sample(n=nmax, random_state=seed)


def _fmt_val(v: float) -> str:
    return f"{float(v):.8g}".replace("-", "m").replace(".", "p")


def _savefig(fig: plt.Figure, stem: Path, save_svg: bool) -> list[str]:
    out = []
    fig.tight_layout()
    p_png = Path(str(stem) + ".png")
    p_pdf = Path(str(stem) + ".pdf")
    fig.savefig(p_png, dpi=220)
    fig.savefig(p_pdf)
    out += [str(p_png), str(p_pdf)]
    if save_svg:
        p_svg = Path(str(stem) + ".svg")
        fig.savefig(p_svg)
        out.append(str(p_svg))
    plt.close(fig)
    return out


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Lifetime landscape scan and plotting workflow")
    p.add_argument("--input", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--mA-values", default="")
    p.add_argument("--tan-beta-values", default="")
    p.add_argument("--lambda6-values", default="")
    p.add_argument("--lambda7", type=float, default=0.0)
    p.add_argument("--sin-ba", type=float, default=1.0)
    p.add_argument("--lambda1-min", type=float, default=None)
    p.add_argument("--lambda1-max", type=float, default=None)
    p.add_argument("--lambda1-targets", default="")
    p.add_argument("--lambda1-tolerance", type=float, default=1e-3)
    p.add_argument("--mode", choices=["summary", "heatmap", "scatter", "slices", "all"], default="all")
    p.add_argument("--ctau-ymin", type=float, default=None)
    p.add_argument("--ctau-ymax", type=float, default=None)
    p.add_argument("--log-ctau", action="store_true", default=True)
    p.add_argument("--max-points-per-panel", type=int, default=50000)
    p.add_argument("--random-seed", type=int, default=42)
    p.add_argument("--style-preset", choices=["paper", "talk", "draft", "compact"], default="paper")
    p.add_argument("--save-svg", action="store_true")
    return p


def _set_style(style_name: str) -> dict:
    s = STYLE_PRESETS[style_name]
    plt.style.use("default")
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 12,
            "axes.labelsize": 13,
            "axes.titlesize": 13,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "axes.grid": True,
            "grid.alpha": s["grid_alpha"],
            "grid.linestyle": ":",
            "legend.frameon": True,
            "figure.figsize": s["figsize"],
        }
    )
    return s


def _group_summary(df: pd.DataFrame, group_cols: list[str], include_lambda1_bin: bool) -> pd.DataFrame:
    d = df.copy()
    if include_lambda1_bin:
        try:
            d["lambda1_bin"] = pd.qcut(d["lambda1"], q=8, duplicates="drop")
            d["lambda1_bin"] = d["lambda1_bin"].astype(str)
            group_cols = group_cols + ["lambda1_bin"]
        except Exception:
            pass

    rows = []
    for keys, g in d.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = {k: v for k, v in zip(group_cols, keys)}
        ct = g["ctau_mm"].to_numpy(dtype=float)
        tw = g["total_width"].to_numpy(dtype=float)
        bg = g["br_gaga_final"].to_numpy(dtype=float)
        bz = g["br_Zga"].to_numpy(dtype=float)
        bb = g["br_bb"].to_numpy(dtype=float)
        pa = g["pareto_gamma_bb"].to_numpy(dtype=float)
        l1 = g["lambda1"].to_numpy(dtype=float)

        row.update(
            {
                "n": int(len(g)),
                "ctau_min": _safe_min(ct),
                "ctau_q01": _quantile_or_none(ct, 0.01),
                "ctau_q05": _quantile_or_none(ct, 0.05),
                "ctau_q50": _quantile_or_none(ct, 0.50),
                "ctau_q95": _quantile_or_none(ct, 0.95),
                "ctau_q99": _quantile_or_none(ct, 0.99),
                "ctau_max": _safe_max(ct),
                "total_width_min": _safe_min(tw),
                "total_width_q05": _quantile_or_none(tw, 0.05),
                "total_width_q50": _quantile_or_none(tw, 0.50),
                "total_width_q95": _quantile_or_none(tw, 0.95),
                "total_width_max": _safe_max(tw),
                "br_gaga_q95": _quantile_or_none(bg, 0.95),
                "br_gaga_max": _safe_max(bg),
                "br_Zga_q95": _quantile_or_none(bz, 0.95),
                "br_Zga_max": _safe_max(bz),
                "br_bb_q05": _quantile_or_none(bb, 0.05),
                "br_bb_min": _safe_min(bb),
                "pareto_gamma_bb_q95": _quantile_or_none(pa, 0.95),
                "pareto_gamma_bb_max": _safe_max(pa),
                "lambda1_min": _safe_min(l1),
                "lambda1_q05": _quantile_or_none(l1, 0.05),
                "lambda1_q50": _quantile_or_none(l1, 0.50),
                "lambda1_q95": _quantile_or_none(l1, 0.95),
                "lambda1_max": _safe_max(l1),
            }
        )

        g_ct = g[np.isfinite(g["ctau_mm"].to_numpy(dtype=float))]
        if len(g_ct):
            i_min = g_ct["ctau_mm"].idxmin()
            i_max = g_ct["ctau_mm"].idxmax()
            row["m_phi_at_min_ctau"] = float(g_ct.loc[i_min, "m_phi"])
            row["m_phi_at_max_ctau"] = float(g_ct.loc[i_max, "m_phi"])
            row["lambda1_at_min_ctau"] = float(g_ct.loc[i_min, "lambda1"])
            row["lambda1_at_max_ctau"] = float(g_ct.loc[i_max, "lambda1"])
        else:
            row["m_phi_at_min_ctau"] = None
            row["m_phi_at_max_ctau"] = None
            row["lambda1_at_min_ctau"] = None
            row["lambda1_at_max_ctau"] = None

        rows.append(row)

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(group_cols).reset_index(drop=True)


def _plot_heatmap_block(df: pd.DataFrame, out_dir: Path, style: dict, save_svg: bool) -> list[str]:
    files = []
    metrics = [
        ("log10_ctau_q50", "Median log10(cτ [mm])", "heatmap_ctau_q50"),
        ("log10_ctau_q05", "Q05 log10(cτ [mm])", "heatmap_ctau_q05"),
        ("log10_ctau_q95", "Q95 log10(cτ [mm])", "heatmap_ctau_q95"),
        ("pareto_gamma_bb_max", "Max BR(γγ)/sqrt(BR(bb))", "heatmap_pareto_max"),
    ]

    gb = (
        df.groupby(["lambda6", "mA", "tan_beta"], dropna=False)
        .agg(
            log10_ctau_q50=("log10_ctau_mm", lambda x: np.quantile(x, 0.5)),
            log10_ctau_q05=("log10_ctau_mm", lambda x: np.quantile(x, 0.05)),
            log10_ctau_q95=("log10_ctau_mm", lambda x: np.quantile(x, 0.95)),
            pareto_gamma_bb_max=("pareto_gamma_bb", "max"),
        )
        .reset_index()
    )

    for lam6, d6 in gb.groupby("lambda6"):
        for col, cbar_label, stem_base in metrics:
            pvt = d6.pivot(index="mA", columns="tan_beta", values=col)
            if pvt.empty:
                continue
            xvals = np.array(sorted(pvt.columns.to_list()), dtype=float)
            yvals = np.array(sorted(pvt.index.to_list()), dtype=float)
            z = pvt.loc[yvals, xvals].to_numpy(dtype=float)

            fig, ax = plt.subplots(figsize=style["figsize"])
            im = ax.imshow(z, origin="lower", aspect="auto", cmap="viridis")
            ax.set_xticks(np.arange(len(xvals)))
            ax.set_xticklabels([f"{x:g}" for x in xvals], rotation=35, ha="right")
            ax.set_yticks(np.arange(len(yvals)))
            ax.set_yticklabels([f"{y:g}" for y in yvals])
            ax.set_xlabel(r"tan$\beta$")
            ax.set_ylabel(r"$m_A$ [GeV]")
            if np.nanmax(xvals) / max(np.nanmin(xvals), 1e-12) > 20:
                ax.set_xlabel(r"tan$\beta$ (log-spaced ticks)")
            ax.set_title(f"{cbar_label}; λ6={lam6:g}")
            cb = fig.colorbar(im, ax=ax)
            cb.set_label(cbar_label)

            stem = out_dir / f"{stem_base}_lambda6_{_fmt_val(lam6)}"
            files += _savefig(fig, stem, save_svg)
    return files


def _plot_lambda1_slices(df: pd.DataFrame, out_dir: Path, style: dict, args: argparse.Namespace) -> tuple[list[str], list[str]]:
    files, notes = [], []
    combos = [
        (300.0, 400000.0, 0.0013),
        (450.0, 250000.0, 0.0010),
    ]
    user_mA = _parse_float_list(args.mA_values)
    user_tb = _parse_float_list(args.tan_beta_values)
    user_l6 = _parse_float_list(args.lambda6_values)
    if user_mA and user_tb and user_l6:
        combos = [(a, b, c) for a in user_mA for b in user_tb for c in user_l6]

    for mA0, tb0, l60 in combos:
        d = df[
            (np.abs(df["mA"] - mA0) <= ATOL["mA"])
            & (np.abs(df["tan_beta"] - tb0) <= ATOL["tan_beta"])
            & (np.abs(df["lambda6"] - l60) <= ATOL["lambda6"])
        ].copy()
        if d.empty:
            notes.append(f"Empty slice: mA={mA0}, tan_beta={tb0}, lambda6={l60}")
            continue

        d = _subsample(d, args.max_points_per_panel, args.random_seed)
        fig, ax = plt.subplots(figsize=style["figsize"])
        sc = ax.scatter(
            d["lambda1"],
            d["ctau_mm"],
            c=d["m_phi"],
            s=style["point_size"],
            alpha=style["alpha"],
            cmap="plasma",
            linewidths=0,
        )
        ax.set_xlabel(r"$\lambda_1$")
        ax.set_ylabel(r"$c\tau$ [mm]")
        if args.log_ctau:
            ax.set_yscale("log")
        if args.ctau_ymin is not None or args.ctau_ymax is not None:
            ax.set_ylim(args.ctau_ymin, args.ctau_ymax)
        ax.set_title(f"Lifetime landscape slice: mA={mA0:g}, tanβ={tb0:g}, λ6={l60:g}")
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label(r"$m_\phi$ [GeV]")

        stem = out_dir / f"lambda1_ctau_slice_mA_{_fmt_val(mA0)}_tb_{_fmt_val(tb0)}_l6_{_fmt_val(l60)}"
        files += _savefig(fig, stem, args.save_svg)
    return files, notes


def _choose_lambda1_targets(df: pd.DataFrame, requested: list[float]) -> list[float]:
    if requested:
        return requested
    lam = np.sort(df["lambda1"].to_numpy(dtype=float))
    lam = lam[np.isfinite(lam)]
    if lam.size == 0:
        return []
    out = [float(np.quantile(lam, q)) for q in (0.2, 0.5, 0.8)]
    dedup = []
    for x in out:
        if not any(abs(x - y) <= 1e-10 for y in dedup):
            dedup.append(x)
    return dedup


def _plot_mphi_curves_fixed_lambda1(
    df: pd.DataFrame,
    out_dir: Path,
    style: dict,
    args: argparse.Namespace,
    lambda1_targets: list[float],
) -> tuple[list[str], dict]:
    files = []
    meta = {"lambda1_targets": lambda1_targets, "target_population": {}}

    for target in lambda1_targets:
        d = df[np.abs(df["lambda1"] - target) <= args.lambda1_tolerance].copy()
        meta["target_population"][f"{target:.8g}"] = int(len(d))
        if d.empty:
            continue

        for (mA0, l60), g in d.groupby(["mA", "lambda6"], dropna=False):
            fig, ax = plt.subplots(figsize=style["figsize"])
            for tb, gt in g.groupby("tan_beta", dropna=False):
                gt2 = gt.sort_values("m_phi")
                gt2 = _subsample(gt2, args.max_points_per_panel, args.random_seed)
                ax.plot(gt2["m_phi"], gt2["ctau_mm"], marker="o", ms=3.4, lw=1.0, alpha=0.85, label=f"tanβ={tb:g}")

            ax.set_xlabel(r"$m_\phi$ [GeV]")
            ax.set_ylabel(r"$c\tau$ [mm]")
            if args.log_ctau:
                ax.set_yscale("log")
            if args.ctau_ymin is not None or args.ctau_ymax is not None:
                ax.set_ylim(args.ctau_ymin, args.ctau_ymax)
            ax.set_title(f"m_phi lifetime curves at λ1≈{target:.4g}, mA={mA0:g}, λ6={l60:g}")
            ax.legend(ncol=2, fontsize=9)

            stem = out_dir / f"mphi_ctau_fixed_lambda1_{_fmt_val(target)}_mA_{_fmt_val(mA0)}_l6_{_fmt_val(l60)}"
            files += _savefig(fig, stem, args.save_svg)

    return files, meta


def _write_top_tables(df: pd.DataFrame, out_dir: Path) -> dict[str, str]:
    d = df.copy()
    d["br_gaga"] = d["br_gaga_final"]
    cols_out = [
        "m_phi",
        "mA",
        "tan_beta",
        "lambda6",
        "lambda7",
        "lambda1",
        "m12",
        "total_width",
        "ctau_mm",
        "br_gaga",
        "br_Zga",
        "br_bb",
        "pareto_gamma_bb",
    ]

    minf = d.sort_values("ctau_mm", ascending=True).head(200)
    maxf = d.sort_values("ctau_mm", ascending=False).head(200)
    parf = d.sort_values("pareto_gamma_bb", ascending=False).head(200)

    p1 = out_dir / "top_min_ctau.csv"
    p2 = out_dir / "top_max_ctau.csv"
    p3 = out_dir / "top_pareto_gamma_bb.csv"
    minf[cols_out].to_csv(p1, index=False)
    maxf[cols_out].to_csv(p2, index=False)
    parf[cols_out].to_csv(p3, index=False)
    return {"top_min_ctau": str(p1), "top_max_ctau": str(p2), "top_pareto_gamma_bb": str(p3)}


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    style = _set_style(args.style_preset)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    lf0 = pl.scan_parquet(args.input)
    schema = lf0.collect_schema().names()

    m_phi_col = _resolve_column(schema, "m_phi")
    mA_col = _resolve_column(schema, "mA")
    tb_col = _resolve_column(schema, "tan_beta")
    l6_col = _resolve_column(schema, "lambda6")
    l7_col = _resolve_column(schema, "lambda7")
    sba_col = _resolve_column(schema, "sin_ba")
    tw_col = _resolve_column(schema, "total_width")
    l1_col = _resolve_column(schema, "lambda1", aliases=["lam1", "lambda_1"])
    brg_col = _first_existing(schema, ["br_gaga"])
    wb_col = _resolve_column(schema, "width_bb")
    wz_col = _resolve_column(schema, "width_Zga")
    wg_col = _resolve_column(schema, "width_gaga")
    m12_col = _first_existing(schema, ["m12"])

    base_cols = [m_phi_col, mA_col, tb_col, l6_col, l7_col, sba_col, l1_col, tw_col, wb_col, wz_col, wg_col]
    if brg_col is not None:
        base_cols.append(brg_col)
    if m12_col is not None:
        base_cols.append(m12_col)

    lf = lf0.select(base_cols)

    mA_values = _parse_float_list(args.mA_values)
    tb_values = _parse_float_list(args.tan_beta_values)
    l6_values = _parse_float_list(args.lambda6_values)
    l1_targets_req = _parse_float_list(args.lambda1_targets)

    lf = _apply_value_filter(lf, mA_col, mA_values, ATOL["mA"])
    lf = _apply_value_filter(lf, tb_col, tb_values, ATOL["tan_beta"])
    lf = _apply_value_filter(lf, l6_col, l6_values, ATOL["lambda6"])

    lf = lf.filter(_float_tol_cut(l7_col, args.lambda7, ATOL["lambda7"]))
    lf = lf.filter(_float_tol_cut(sba_col, args.sin_ba, ATOL["sin_ba"]))

    if args.lambda1_min is not None:
        lf = lf.filter(pl.col(l1_col).cast(pl.Float64) >= float(args.lambda1_min))
    if args.lambda1_max is not None:
        lf = lf.filter(pl.col(l1_col).cast(pl.Float64) <= float(args.lambda1_max))

    brg_expr = pl.col(brg_col).cast(pl.Float64) if brg_col is not None else (pl.col(wg_col).cast(pl.Float64) / pl.col(tw_col).cast(pl.Float64))

    lf = (
        lf.with_columns(
            [
                pl.col(m_phi_col).cast(pl.Float64).alias("m_phi"),
                pl.col(mA_col).cast(pl.Float64).alias("mA"),
                pl.col(tb_col).cast(pl.Float64).alias("tan_beta"),
                pl.col(l6_col).cast(pl.Float64).alias("lambda6"),
                pl.col(l7_col).cast(pl.Float64).alias("lambda7"),
                pl.col(sba_col).cast(pl.Float64).alias("sin_ba"),
                pl.col(l1_col).cast(pl.Float64).alias("lambda1"),
                pl.col(tw_col).cast(pl.Float64).alias("total_width"),
                pl.col(wb_col).cast(pl.Float64).alias("width_bb"),
                pl.col(wz_col).cast(pl.Float64).alias("width_Zga"),
                pl.col(wg_col).cast(pl.Float64).alias("width_gaga"),
                brg_expr.alias("br_gaga"),
                (pl.col(m12_col).cast(pl.Float64).alias("m12") if m12_col is not None else pl.lit(None).cast(pl.Float64).alias("m12")),
            ]
        )
        .with_columns(
            [
                (pl.lit(HBAR_C_GEV_MM) / pl.col("total_width")).alias("ctau_mm"),
                (pl.col("width_bb") / pl.col("total_width")).alias("br_bb"),
                (pl.col("width_Zga") / pl.col("total_width")).alias("br_Zga"),
            ]
        )
        .with_columns(
            [
                pl.col("br_gaga").alias("br_gaga_final"),
                (pl.col("br_gaga") / (pl.col("br_bb") + pl.lit(1e-12)).sqrt()).alias("pareto_gamma_bb"),
                pl.col("ctau_mm").log10().alias("log10_ctau_mm"),
            ]
        )
        .filter(pl.col("total_width") > 0)
        .filter(pl.col("ctau_mm").is_finite() & (pl.col("ctau_mm") > 0))
        .filter(pl.col("br_bb").is_finite() & pl.col("br_Zga").is_finite() & pl.col("br_gaga_final").is_finite())
    )

    df_pl = _stream_collect(lf)
    if df_pl.height == 0:
        raise SystemExit("No rows remain after filters; nothing to plot.")
    df = df_pl.to_pandas()

    summary_meta = {
        "input": str(args.input),
        "output_dir": str(out_dir),
        "rows_after_filters": int(len(df)),
        "m12_label": "$m_{12}^{2}$" if m12_col is not None else "m12",
        "m12_semantics_note": "PARQUET_SCHEMA.md documents m12 as m12^2" if m12_col is not None else "m12 column not present; semantics uncertain",
        "empty_slices": [],
        "chosen_lambda1_targets": [],
    }

    run_summary = args.mode in ("summary", "all")
    run_heatmap = args.mode in ("heatmap", "all")
    run_scatter = args.mode in ("scatter", "all")
    run_slices = args.mode in ("slices", "all")

    created_figs: list[str] = []

    if run_summary:
        include_lambda1_bin = df["lambda1"].nunique() > 12
        gsum = _group_summary(df, ["mA", "tan_beta", "lambda6"], include_lambda1_bin=include_lambda1_bin)
        p_csv = out_dir / "ctau_landscape_summary.csv"
        p_json = out_dir / "ctau_landscape_summary.json"
        gsum.to_csv(p_csv, index=False)
        with open(p_json, "w", encoding="utf-8") as f:
            json.dump(
                {
                    "meta": summary_meta,
                    "summary_rows": gsum.to_dict(orient="records"),
                },
                f,
                indent=2,
            )

    if run_heatmap:
        created_figs += _plot_heatmap_block(df, out_dir, style, args.save_svg)

    if run_scatter:
        f_scatter, notes = _plot_lambda1_slices(df, out_dir, style, args)
        created_figs += f_scatter
        summary_meta["empty_slices"].extend(notes)

    if run_slices:
        targets = _choose_lambda1_targets(df, l1_targets_req)
        summary_meta["chosen_lambda1_targets"] = targets
        f2, target_meta = _plot_mphi_curves_fixed_lambda1(df, out_dir, style, args, targets)
        created_figs += f2
        summary_meta["target_population"] = target_meta["target_population"]

    top_paths = _write_top_tables(df, out_dir)

    # update JSON summary with post-run metadata
    if run_summary:
        p_json = out_dir / "ctau_landscape_summary.json"
        with open(p_json, "r", encoding="utf-8") as f:
            payload = json.load(f)
        payload["meta"] = {**payload.get("meta", {}), **summary_meta, "figures": created_figs, "top_tables": top_paths}
        with open(p_json, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)

    logger.info("Done. Rows after filters: %d", len(df))
    logger.info("Output dir: %s", out_dir)


if __name__ == "__main__":
    main()
