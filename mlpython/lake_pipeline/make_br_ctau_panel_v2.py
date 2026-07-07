#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import polars as pl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# shared constant (was 1.973269804e-13); see physics_conventions.py
from physics_conventions import HBAR_C_GEV_MM as HBARC_GEV_MM

CHANNELS = [
    ("br_gaga", "width_gaga", r"BR$(H\to\gamma\gamma)$", "o"),
    ("br_bb", "width_bb", r"BR$(H\to b\bar{b})$", "s"),
    ("br_tautau", "width_tautau", r"BR$(H\to\tau\tau)$", "^"),
    ("br_WW", "width_WW", r"BR$(H\to WW)$", "D"),
    ("br_ZZ", "width_ZZ", r"BR$(H\to ZZ)$", "v"),
]


def tag(x: float | str) -> str:
    s = f"{x:g}" if isinstance(x, float) else str(x)
    return s.replace("+", "").replace("-", "m").replace(".", "p")


def parse_only(s: str | None):
    if not s:
        return None
    vals = [float(x.strip()) for x in s.split(",")]
    if len(vals) != 3:
        raise SystemExit("--only must be: mA_target,tan_beta,lambda6")
    return tuple(vals)


def close(a: float, b: float, rel: float = 1e-9, abs_: float = 1e-12) -> bool:
    return math.isclose(float(a), float(b), rel_tol=rel, abs_tol=abs_)


def add_observables(df: pl.DataFrame) -> pl.DataFrame:
    exprs = []

    if "ctau_mm" not in df.columns and "total_width" in df.columns:
        exprs.append(
            pl.when(pl.col("total_width") > 0)
            .then(HBARC_GEV_MM / pl.col("total_width"))
            .otherwise(None)
            .alias("ctau_mm")
        )

    for br, width, _, _ in CHANNELS:
        if br not in df.columns and width in df.columns and "total_width" in df.columns:
            exprs.append(
                pl.when(pl.col("total_width") > 0)
                .then(pl.col(width) / pl.col("total_width"))
                .otherwise(None)
                .alias(br)
            )

    return df.with_columns(exprs) if exprs else df


def phys_filter(df: pl.DataFrame) -> pl.DataFrame:
    need = {"positivity_ok", "unitarity_ok", "perturbativity_ok"}
    if need <= set(df.columns):
        return df.filter(
            (pl.col("positivity_ok") == 1)
            & (pl.col("unitarity_ok") == 1)
            & (pl.col("perturbativity_ok") == 1)
        )
    return df


def fixed_values(df: pl.DataFrame) -> dict[str, float]:
    r = df.head(1).to_dicts()[0]
    return {
        "mA_target": float(r["mA_target"] if "mA_target" in r else r["mA"]),
        "tan_beta": float(r["tan_beta"]),
        "lambda6": float(r["lambda6"]),
        "lambda7": float(r.get("lambda7", 0.0)),
    }


def make_title(vals: dict[str, float], suffix: str) -> str:
    return (
        rf"$m_A=m_{{H^\pm}}={vals['mA_target']:g}$ GeV, "
        rf"$\tan\beta={vals['tan_beta']:g}$, "
        rf"$\lambda_6={vals['lambda6']:g}$, "
        rf"$\lambda_7={vals['lambda7']:g}$, "
        + suffix
    )


def plot_exact(df: pl.DataFrame, outbase: Path, vals: dict[str, float], lam_col: str, target: float) -> None:
    df = add_observables(phys_filter(df))
    df = df.filter((pl.col(lam_col) - target).abs() <= 1e-12)

    if df.height == 0:
        print(f"[skip] exact empty: {outbase}")
        return

    xcol = "mH_target" if "mH_target" in df.columns else "m_phi"
    pdf = df.sort(xcol).to_pandas()

    fig, (ax_br, ax_ctau) = plt.subplots(
        2, 1,
        figsize=(8.8, 6.5),
        sharex=True,
        gridspec_kw={"height_ratios": [2.1, 1.0]},
    )

    for br, _, label, marker in CHANNELS:
        if br in pdf.columns:
            ax_br.plot(
                pdf[xcol],
                pdf[br],
                marker=marker,
                markersize=3.5,
                linewidth=1.5,
                label=label,
            )

    ax_br.set_yscale("log")
    ax_br.set_ylabel("Branching ratio")
    ax_br.grid(True, which="both", alpha=0.25)
    ax_br.legend(fontsize=8, loc="best")
    ax_br.set_title(make_title(vals, rf"exact $\lambda_1={target:g}$"))

    ax_ctau.plot(
        pdf[xcol],
        pdf["ctau_mm"],
        marker="x",
        markersize=3.5,
        linewidth=1.4,
        label=r"$c\tau=\hbar c/\Gamma_\mathrm{tot}$",
    )
    ax_ctau.set_yscale("log")
    ax_ctau.set_ylabel(r"$c\tau$ [mm]")
    ax_ctau.set_xlabel(r"$m_H$ target / $m_\phi$ [GeV]")
    ax_ctau.grid(True, which="both", alpha=0.25)
    ax_ctau.legend(fontsize=8, loc="best")

    lam_min = float(pdf[lam_col].min())
    lam_max = float(pdf[lam_col].max())
    maxerr = max(abs(lam_min - target), abs(lam_max - target))

    fig.text(
        0.01,
        0.01,
        rf"$n={len(pdf)}$, min/max $\lambda_1$ = {lam_min:.16g} / {lam_max:.16g}, "
        rf"max$|\lambda_1-1|$ = {maxerr:.3e}",
        fontsize=8,
    )

    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(outbase.with_suffix(".png"), dpi=190)
    fig.savefig(outbase.with_suffix(".pdf"))
    plt.close(fig)
    print(f"[wrote] {outbase.with_suffix('.png')}")


def collect_near_from_full(
    full_parquet: Path,
    vals: dict[str, float],
    target: float,
    near_tol: float,
    max_rows: int,
) -> pl.DataFrame:
    lf = pl.scan_parquet(full_parquet)
    schema = set(lf.collect_schema().names())
    lam_col = "lam1" if "lam1" in schema else "computed_lam1"

    filters = [
        (pl.col("mA_target") - vals["mA_target"]).abs() <= 1e-7,
        (pl.col("tan_beta") - vals["tan_beta"]).abs() <= max(1e-7, abs(vals["tan_beta"]) * 1e-12),
        (pl.col("lambda6") - vals["lambda6"]).abs() <= max(1e-14, abs(vals["lambda6"]) * 1e-9),
        (pl.col("lambda7") - vals["lambda7"]).abs() <= 1e-14,
        (pl.col(lam_col) - target).abs() <= near_tol,
    ]

    lf = lf.filter(filters[0])
    for f in filters[1:]:
        lf = lf.filter(f)

    if {"positivity_ok", "unitarity_ok", "perturbativity_ok"} <= schema:
        lf = lf.filter(
            (pl.col("positivity_ok") == 1)
            & (pl.col("unitarity_ok") == 1)
            & (pl.col("perturbativity_ok") == 1)
        )

    df = lf.collect()
    df = add_observables(df)
    df = df.with_columns((pl.col(lam_col) - target).alias("delta_lam1"))

    if df.height > max_rows:
        xcol = "mH_target" if "mH_target" in df.columns else "m_phi"
        df = df.sort([xcol, "delta_lam1"])
        idx = [int(i * (df.height - 1) / (max_rows - 1)) for i in range(max_rows)]
        df = df.with_row_index("_i").filter(pl.col("_i").is_in(idx)).drop("_i")

    return df


def plot_near(df: pl.DataFrame, outbase: Path, vals: dict[str, float], target: float, near_tol: float) -> None:
    if df.height == 0:
        print(f"[skip] near empty: {outbase}")
        return

    xcol = "mH_target" if "mH_target" in df.columns else "m_phi"
    pdf = df.sort(xcol).to_pandas()

    vmax = max(abs(float(pdf["delta_lam1"].min())), abs(float(pdf["delta_lam1"].max())))
    if vmax == 0:
        vmax = near_tol
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)

    fig, (ax_br, ax_ctau) = plt.subplots(
        2, 1,
        figsize=(9.0, 6.8),
        sharex=True,
        gridspec_kw={"height_ratios": [2.1, 1.0]},
    )

    sc = None
    for br, _, label, marker in CHANNELS:
        if br in pdf.columns:
            sc = ax_br.scatter(
                pdf[xcol],
                pdf[br],
                c=pdf["delta_lam1"],
                cmap="coolwarm",
                norm=norm,
                marker=marker,
                s=18,
                alpha=0.82,
                label=label,
                linewidths=0.2,
            )

    ax_br.set_yscale("log")
    ax_br.set_ylabel("Branching ratio")
    ax_br.grid(True, which="both", alpha=0.25)
    ax_br.legend(fontsize=8, loc="best")
    ax_br.set_title(make_title(vals, rf"near $\lambda_1={target:g}$, $|\lambda_1-1|\leq {near_tol:g}$"))

    sc_ctau = ax_ctau.scatter(
        pdf[xcol],
        pdf["ctau_mm"],
        c=pdf["delta_lam1"],
        cmap="coolwarm",
        norm=norm,
        marker="x",
        s=18,
        alpha=0.82,
        label=r"$c\tau=\hbar c/\Gamma_\mathrm{tot}$",
    )
    if sc is None:
        sc = sc_ctau

    ax_ctau.set_yscale("log")
    ax_ctau.set_ylabel(r"$c\tau$ [mm]")
    ax_ctau.set_xlabel(r"$m_H$ target / $m_\phi$ [GeV]")
    ax_ctau.grid(True, which="both", alpha=0.25)
    ax_ctau.legend(fontsize=8, loc="best")

    cbar = fig.colorbar(sc, ax=[ax_br, ax_ctau], pad=0.02)
    cbar.set_label(r"$\lambda_1-1$")

    fig.text(
        0.01,
        0.01,
        rf"$n={len(pdf)}$, min/max $(\lambda_1-1)$ = "
        rf"{float(pdf['delta_lam1'].min()):.3e} / {float(pdf['delta_lam1'].max()):.3e}",
        fontsize=8,
    )

    fig.tight_layout(rect=(0, 0.035, 0.93, 1))
    fig.savefig(outbase.with_suffix(".png"), dpi=190)
    fig.savefig(outbase.with_suffix(".pdf"))
    plt.close(fig)
    print(f"[wrote] {outbase.with_suffix('.png')}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--slice-root", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--only", default=None, help="mA_target,tan_beta,lambda6")
    ap.add_argument("--lambda1-target", type=float, default=1.0)
    ap.add_argument("--near-tol", type=float, default=1e-4)
    ap.add_argument("--max-near-rows", type=int, default=6000)
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    only = parse_only(args.only)

    slice_paths = sorted(args.slice_root.glob("mA=*/lambda6=*/tan_beta=*/slice.parquet"))
    print(f"[info] found slice files: {len(slice_paths)}")

    made = 0

    for sp in slice_paths:
        df = pl.read_parquet(sp)
        if df.height == 0:
            continue

        vals = fixed_values(df)

        if only is not None:
            path_s = str(sp)
            path_match = (
                f"mA={only[0]:g}" in path_s
                and (f"tan_beta={only[1]:g}" in path_s or f"tan_beta={only[1]:.0f}" in path_s)
                and f"lambda6={only[2]:g}" in path_s
            )
            value_match = (
                close(vals["mA_target"], only[0], abs_=1e-5)
                and close(vals["tan_beta"], only[1], rel=1e-9, abs_=1e-3)
                and close(vals["lambda6"], only[2], rel=1e-6, abs_=1e-12)
            )
            print(
                f"[debug] slice={sp} "
                f"mA={vals['mA_target']:.17g}, tb={vals['tan_beta']:.17g}, l6={vals['lambda6']:.17g}, "
                f"path_match={path_match}, value_match={value_match}"
            )
            if not (path_match or value_match):
                continue

        lam_col = "lam1" if "lam1" in df.columns else "computed_lam1"
        base = f"mA{tag(vals['mA_target'])}_tb{tag(vals['tan_beta'])}_l6{tag(vals['lambda6'])}"

        exact_out = args.output_dir / f"br_ctau_panel_exact_lam1_{base}"
        plot_exact(df, exact_out, vals, lam_col, args.lambda1_target)
        made += 1

        near = collect_near_from_full(
            full_parquet=args.input,
            vals=vals,
            target=args.lambda1_target,
            near_tol=args.near_tol,
            max_rows=args.max_near_rows,
        )

        print(
            f"[info] near points for {base}: n={near.height}, "
            f"tol={args.near_tol}"
        )

        near_out = args.output_dir / f"br_ctau_panel_near_lam1_{base}_tol{tag(args.near_tol)}"
        plot_near(near, near_out, vals, args.lambda1_target, args.near_tol)

        audit_cols = [
            "mH_target", "m_phi", "mA_target", "mA", "tan_beta", "lambda6", "lambda7",
            lam_col, "delta_lam1", "total_width", "ctau_mm",
            "br_gaga", "br_bb", "br_tautau", "br_WW", "br_ZZ",
        ]
        audit_cols = [c for c in audit_cols if c in near.columns]
        near.select(audit_cols).write_csv(args.output_dir / f"near_lam1_points_{base}.csv")
        made += 1

    print(f"[done] made plot groups: {made}")
    print(f"[done] output dir: {args.output_dir}")

    if made == 0:
        raise SystemExit("No plots were made. Check --only values or slice-root.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
