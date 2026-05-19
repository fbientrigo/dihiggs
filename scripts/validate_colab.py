#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_COLAB_CSV = Path("colab_points_for_2hdmc_validation.csv")
DEFAULT_CALC_CSV = Path("colab_points_validated_2hdmc.csv")
DEFAULT_OUT_CSV = Path("colab_vs_2hdmc_merged_comparison.csv")

REQUIRED_COLAB_COLUMNS = [
    "point_id",
    "br_bb_colab", "br_gaga_colab", "br_Zga_colab",
    "width_bb_colab", "width_cc_colab", "width_tautau_colab",
    "width_gg_colab", "width_gaga_colab", "width_Zga_colab",
    "total_width_colab",
]

REQUIRED_2HDMC_COLUMNS = [
    "point_id",
    "br_bb_2hdmc", "br_gaga_2hdmc", "br_Zga_2hdmc",
    "width_bb", "width_tautau", "width_gg", "width_gaga", "width_Zga", "total_width",
]


def pick_col(df: pd.DataFrame, base: str) -> str:
    candidates = [base, f"{base}_x", f"{base}_y"]
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"Could not find any of: {candidates}")


def to_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def safe_rel_err(a: pd.Series, b: pd.Series, eps: float = 1e-30) -> pd.Series:
    a = to_num(a)
    b = to_num(b)
    denom = b.abs().clip(lower=eps)
    out = (a - b).abs() / denom
    out = out.where(b.abs() > eps, np.nan)
    return out


def safe_fraction(num: pd.Series, den: pd.Series, eps: float = 1e-30) -> pd.Series:
    num = to_num(num)
    den = to_num(den)
    out = num / den.clip(lower=eps)
    out = out.where(den.abs() > eps, np.nan)
    return out


def _validate_missing_point_id(df: pd.DataFrame, label: str) -> None:
    if "point_id" not in df.columns:
        raise ValueError(f"{label} CSV missing required column: point_id")


def _validate_duplicate_point_id(df: pd.DataFrame, label: str) -> None:
    dup = df["point_id"].duplicated(keep=False)
    if dup.any():
        n = int(dup.sum())
        examples = df.loc[dup, "point_id"].astype(str).head(5).tolist()
        raise ValueError(f"{label} CSV has duplicate point_id values (n={n}), examples={examples}")


def _validate_required_columns(df: pd.DataFrame, required: list[str], label: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{label} CSV missing required columns: {missing}")


def validate_input_frames(colab: pd.DataFrame, calc: pd.DataFrame) -> None:
    _validate_missing_point_id(colab, "Colab")
    _validate_missing_point_id(calc, "2HDMC")
    _validate_duplicate_point_id(colab, "Colab")
    _validate_duplicate_point_id(calc, "2HDMC")
    _validate_required_columns(colab, REQUIRED_COLAB_COLUMNS, "Colab")
    _validate_required_columns(calc, REQUIRED_2HDMC_COLUMNS, "2HDMC")


def build_merged_comparison(colab_csv: Path, calc_csv: Path) -> pd.DataFrame:
    colab = pd.read_csv(colab_csv)
    calc = pd.read_csv(calc_csv)
    validate_input_frames(colab, calc)

    colab_reduced = colab[REQUIRED_COLAB_COLUMNS].copy()
    df = calc.merge(colab_reduced, on="point_id", how="left", validate="one_to_one")

    bb_colab = pick_col(df, "br_bb_colab")
    gaga_colab = pick_col(df, "br_gaga_colab")
    zga_colab = pick_col(df, "br_Zga_colab")

    bb_2hdmc = pick_col(df, "br_bb_2hdmc")
    gaga_2hdmc = pick_col(df, "br_gaga_2hdmc")
    zga_2hdmc = pick_col(df, "br_Zga_2hdmc")

    wbb_colab = pick_col(df, "width_bb_colab")
    wcc_colab = pick_col(df, "width_cc_colab")
    wtautau_colab = pick_col(df, "width_tautau_colab")
    wgg_colab = pick_col(df, "width_gg_colab")
    wgaga_colab = pick_col(df, "width_gaga_colab")
    wzga_colab = pick_col(df, "width_Zga_colab")
    wtot_colab = pick_col(df, "total_width_colab")

    wbb_2hdmc = pick_col(df, "width_bb")
    wtautau_2hdmc = pick_col(df, "width_tautau")
    wgg_2hdmc = pick_col(df, "width_gg")
    wgaga_2hdmc = pick_col(df, "width_gaga")
    wzga_2hdmc = pick_col(df, "width_Zga")
    wtot_2hdmc = pick_col(df, "total_width")

    has_cc_2hdmc = "width_cc" in df.columns
    wcc_2hdmc = pick_col(df, "width_cc") if has_cc_2hdmc else None

    br_pairs = {
        "bb": (bb_colab, bb_2hdmc),
        "gaga": (gaga_colab, gaga_2hdmc),
        "Zga": (zga_colab, zga_2hdmc),
    }
    for ch, (c_colab, c_2hdmc) in br_pairs.items():
        df[f"abs_err_br_{ch}"] = (to_num(df[c_2hdmc]) - to_num(df[c_colab])).abs()
        df[f"rel_err_br_{ch}"] = safe_rel_err(df[c_colab], df[c_2hdmc])

    width_pairs = {
        "bb": (wbb_colab, wbb_2hdmc),
        "tautau": (wtautau_colab, wtautau_2hdmc),
        "gg": (wgg_colab, wgg_2hdmc),
        "gaga": (wgaga_colab, wgaga_2hdmc),
        "Zga": (wzga_colab, wzga_2hdmc),
        "total": (wtot_colab, wtot_2hdmc),
    }
    if has_cc_2hdmc:
        width_pairs["cc"] = (wcc_colab, wcc_2hdmc)

    for ch, (c_colab, c_2hdmc) in width_pairs.items():
        df[f"abs_err_width_{ch}"] = (to_num(df[c_2hdmc]) - to_num(df[c_colab])).abs()
        df[f"rel_err_width_{ch}"] = safe_rel_err(df[c_colab], df[c_2hdmc])

    subset_colab = (
        to_num(df[wbb_colab]) + to_num(df[wcc_colab]) + to_num(df[wtautau_colab]) +
        to_num(df[wgg_colab]) + to_num(df[wgaga_colab]) + to_num(df[wzga_colab])
    )
    subset_2hdmc = (
        to_num(df[wbb_2hdmc]) + to_num(df[wtautau_2hdmc]) + to_num(df[wgg_2hdmc]) +
        to_num(df[wgaga_2hdmc]) + to_num(df[wzga_2hdmc])
    )
    if has_cc_2hdmc:
        subset_2hdmc = subset_2hdmc + to_num(df[wcc_2hdmc])

    df["subset_width_colab"] = subset_colab
    df["subset_width_2hdmc"] = subset_2hdmc

    frac_defs = {
        "bb": (wbb_colab, wbb_2hdmc),
        "tautau": (wtautau_colab, wtautau_2hdmc),
        "gg": (wgg_colab, wgg_2hdmc),
        "gaga": (wgaga_colab, wgaga_2hdmc),
        "Zga": (wzga_colab, wzga_2hdmc),
    }
    if has_cc_2hdmc:
        frac_defs["cc"] = (wcc_colab, wcc_2hdmc)

    for ch, (c_colab, c_2hdmc) in frac_defs.items():
        df[f"f_{ch}_colab"] = safe_fraction(df[c_colab], df["subset_width_colab"])
        df[f"f_{ch}_2hdmc"] = safe_fraction(df[c_2hdmc], df["subset_width_2hdmc"])
        df[f"abs_err_f_{ch}"] = (df[f"f_{ch}_2hdmc"] - df[f"f_{ch}_colab"]).abs()
        df[f"rel_err_f_{ch}"] = safe_rel_err(df[f"f_{ch}_colab"], df[f"f_{ch}_2hdmc"])

    return df


def summarize_comparison(df: pd.DataFrame) -> pd.Series:
    summary = {
        "n_points": len(df),
        "n_set_param_ok": int(df["set_param_ok"].sum()) if "set_param_ok" in df.columns else None,
        "n_triple_ok": int(df["triple_ok"].sum()) if "triple_ok" in df.columns else None,
        "n_warning_flag": int(df["warning_flag"].sum()) if "warning_flag" in df.columns else None,
        "mae_br_bb": df["abs_err_br_bb"].mean() if "abs_err_br_bb" in df.columns else np.nan,
        "mae_br_gaga": df["abs_err_br_gaga"].mean() if "abs_err_br_gaga" in df.columns else np.nan,
        "mae_br_Zga": df["abs_err_br_Zga"].mean() if "abs_err_br_Zga" in df.columns else np.nan,
        "median_rel_br_bb": df["rel_err_br_bb"].median() if "rel_err_br_bb" in df.columns else np.nan,
        "median_rel_br_gaga": df["rel_err_br_gaga"].median() if "rel_err_br_gaga" in df.columns else np.nan,
        "median_rel_br_Zga": df["rel_err_br_Zga"].median() if "rel_err_br_Zga" in df.columns else np.nan,
        "mae_width_bb": df["abs_err_width_bb"].mean() if "abs_err_width_bb" in df.columns else np.nan,
        "mae_width_tautau": df["abs_err_width_tautau"].mean() if "abs_err_width_tautau" in df.columns else np.nan,
        "mae_width_gg": df["abs_err_width_gg"].mean() if "abs_err_width_gg" in df.columns else np.nan,
        "mae_width_gaga": df["abs_err_width_gaga"].mean() if "abs_err_width_gaga" in df.columns else np.nan,
        "mae_width_Zga": df["abs_err_width_Zga"].mean() if "abs_err_width_Zga" in df.columns else np.nan,
        "mae_width_total": df["abs_err_width_total"].mean() if "abs_err_width_total" in df.columns else np.nan,
        "median_rel_width_total": df["rel_err_width_total"].median() if "rel_err_width_total" in df.columns else np.nan,
    }
    return pd.Series(summary)


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge Colab and 2HDMC validation CSVs and compute comparison metrics.")
    parser.add_argument("--colab-csv", type=Path, default=DEFAULT_COLAB_CSV)
    parser.add_argument("--calc-csv", type=Path, default=DEFAULT_CALC_CSV)
    parser.add_argument("--out-csv", type=Path, default=DEFAULT_OUT_CSV)
    parser.add_argument("--top-n", type=int, default=10)
    return parser.parse_args()


def main() -> int:
    args = cli()
    df = build_merged_comparison(args.colab_csv, args.calc_csv)

    print(summarize_comparison(df))

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    print(f"\nWrote merged comparison to: {args.out_csv}")

    rank_col = "abs_err_f_bb" if "abs_err_f_bb" in df.columns else "abs_err_width_gaga"
    show_cols = [
        "point_id", "mH", "lambda1", "lambda6", "tan_beta", "triple_ok",
        "abs_error", "warning_flag", rank_col,
    ]
    show_cols = [c for c in show_cols if c in df.columns]
    print(f"\nTop {args.top_n} by {rank_col}:")
    print(df.sort_values(rank_col, ascending=False)[show_cols].head(args.top_n).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
