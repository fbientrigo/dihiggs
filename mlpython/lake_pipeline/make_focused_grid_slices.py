#!/usr/bin/env python3
"""Build deterministic focused-grid slices from selected combos.

For each co-occurring ``(mA_target, tan_beta, lambda6)`` combo in
``selected_values.json`` this writes one ``slice.parquet`` containing the exact
rows from the Parquet (no interpolation), with derived ``ctau_mm`` and BR
columns, plus a CSV/JSON manifest.

Never runs unbounded global family plotting; only writes to focused-grid output
paths.
"""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import polars as pl

HBARC_GEV_MM = 1.973269804e-13

PHYS_COLS = ["positivity_ok", "unitarity_ok", "perturbativity_ok"]
# BR channel name -> partial width column.
BR_CHANNELS = {
    "br_bb": "width_bb",
    "br_tautau": "width_tautau",
    "br_WW": "width_WW",
    "br_ZZ": "width_ZZ",
    "br_gaga": "width_gaga",
    "br_Zga": "width_Zga",
    "br_gg": "width_gg",
    "br_hh": "width_hh",
}
# Candidate charged-Higgs mass columns for the mA = mH+ check.
MHPLUS_CANDIDATES = ["mHp_calibrated", "mHp", "mHplus", "mH_plus", "mHpm"]
GROUP_KEY = [
    "mA_target", "tan_beta", "lambda6", "lambda7", "variation_idx", "mH_target",
]
MANIFEST_FIELDS = [
    "slice_path", "mA_target", "tan_beta", "lambda6", "lambda7",
    "lambda1_target", "lambda1_max_abs_error", "n_rows", "n_mH_target",
    "mH_min", "mH_max", "total_width_min", "total_width_max", "ctau_min",
    "ctau_max", "br_gaga_max", "br_bb_max", "br_Zga_max",
    "physical_filter_status", "physical_filter_missing_columns",
    "mAeqmHp_check_status", "mAeqmHp_max_abs_diff", "created_utc",
    "input_parquet", "selected_values_json",
]


def fmt(value: float) -> str:
    """Filesystem-friendly deterministic float label (e.g. 400, 1e-05)."""
    return f"{value:g}"


def br_and_ctau_exprs(schema_names: set[str]) -> tuple[list[pl.Expr], list[str]]:
    """Build ctau + BR expressions, skipping channels with missing widths."""
    exprs: list[pl.Expr] = []
    missing: list[str] = []
    tw = pl.col("total_width").cast(pl.Float64)
    good_tw = tw.is_not_null() & (tw > 0)
    exprs.append(
        pl.when(good_tw).then(HBARC_GEV_MM / tw).otherwise(None).alias("ctau_mm")
    )
    for br_col, width_col in BR_CHANNELS.items():
        if width_col not in schema_names:
            missing.append(width_col)
            continue
        exprs.append(
            pl.when(good_tw)
            .then(pl.col(width_col).cast(pl.Float64) / tw)
            .otherwise(None)
            .alias(br_col)
        )
    return exprs, missing


def physical_filter(lf: pl.LazyFrame, schema_names: set[str], apply: bool):
    present = [c for c in PHYS_COLS if c in schema_names]
    missing = [c for c in PHYS_COLS if c not in schema_names]
    if not apply:
        return lf, "not_applied_disabled", missing
    if missing:
        return lf, "not_applied_missing_columns", missing
    cond = pl.lit(True)
    for c in present:
        cond = cond & (pl.col(c).cast(pl.Float64) == 1)
    return lf.filter(cond), "applied", missing


def nearest_lambda1(df: pl.DataFrame, target: float) -> pl.DataFrame:
    """Keep the row with lambda1 nearest to target per group key."""
    df = df.with_columns(
        (pl.col("lam1").cast(pl.Float64) - target).abs().alias("lambda1_abs_error")
    )
    keys = [k for k in GROUP_KEY if k in df.columns]
    df = df.sort(keys + ["lambda1_abs_error", "lam1"])
    return df.group_by(keys, maintain_order=True).first()


def mAeqmHp_status(df: pl.DataFrame, schema_names: set[str]):
    col = next((c for c in MHPLUS_CANDIDATES if c in schema_names), None)
    if col is None or "mA" not in schema_names or df.height == 0:
        return "not_assessable_missing_columns", None
    diff = (
        df.select(
            (pl.col("mA").cast(pl.Float64) - pl.col(col).cast(pl.Float64))
            .abs()
            .max()
            .alias("d")
        )["d"].item()
    )
    if diff is None:
        return "not_assessable_null", None
    status = "ok" if diff <= 1e-6 else "mismatch"
    return status, float(diff)


def safe_stat(df: pl.DataFrame, col: str, agg: str):
    if col not in df.columns or df.height == 0:
        return None
    expr = getattr(pl.col(col), agg)()
    val = df.select(expr.alias("v"))["v"].item()
    return val


def build_combo_slice(
    lf_base: pl.LazyFrame,
    combo: dict,
    target: float,
    schema_names: set[str],
    br_exprs: list[pl.Expr],
) -> pl.DataFrame:
    lf = lf_base.filter(
        (pl.col("mA_target") == combo["mA_target"])
        & (pl.col("tan_beta") == combo["tan_beta"])
        & (pl.col("lambda6") == combo["lambda6"])
    )
    df = lf.with_columns(br_exprs).collect()
    if df.height == 0:
        return df
    return nearest_lambda1(df, target)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Path to silver_all.parquet")
    p.add_argument("--selected-values", required=True,
                   help="Path to selected_values.json")
    p.add_argument("--output-dir", required=True,
                   help="Base focused-grid output directory")
    p.add_argument("--apply-phys-filter", action="store_true",
                   help="Apply *_ok==1 filters when those columns are present")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    input_path = Path(args.input)
    selected_path = Path(args.selected_values)
    output_dir = Path(args.output_dir)

    selected = json.loads(selected_path.read_text())
    combos = selected["combos"]
    lambda7 = float(selected["lambda7"])
    target = float(selected["lambda1_target"])
    variation_idx = int(selected["variation_idx"])

    group_dir = output_dir / (
        f"lambda1={fmt(target)}_lambda7={fmt(lambda7)}_mAeqmHp"
    )
    group_dir.mkdir(parents=True, exist_ok=True)

    lf = pl.scan_parquet(input_path)
    schema_names = set(lf.collect_schema().names())

    base = lf.filter(
        (pl.col("lambda7") == lambda7) & (pl.col("variation_idx") == variation_idx)
    )
    base, phys_status, phys_missing = physical_filter(
        base, schema_names, args.apply_phys_filter
    )
    br_exprs, br_missing = br_and_ctau_exprs(schema_names)
    created = datetime.now(timezone.utc).isoformat()

    manifest_rows: list[dict] = []
    for combo in combos:
        df = build_combo_slice(base, combo, target, schema_names, br_exprs)
        slice_dir = (
            group_dir
            / f"mA={fmt(combo['mA_target'])}"
            / f"lambda6={fmt(combo['lambda6'])}"
            / f"tan_beta={fmt(combo['tan_beta'])}"
        )
        slice_dir.mkdir(parents=True, exist_ok=True)
        slice_path = slice_dir / "slice.parquet"
        df.write_parquet(slice_path)

        mhp_status, mhp_diff = mAeqmHp_status(df, schema_names)
        manifest_rows.append({
            "slice_path": str(slice_path),
            "mA_target": combo["mA_target"],
            "tan_beta": combo["tan_beta"],
            "lambda6": combo["lambda6"],
            "lambda7": combo["lambda7"],
            "lambda1_target": target,
            "lambda1_max_abs_error": safe_stat(df, "lambda1_abs_error", "max"),
            "n_rows": df.height,
            "n_mH_target": (df["mH_target"].n_unique()
                            if "mH_target" in df.columns else None),
            "mH_min": safe_stat(df, "mH_target", "min"),
            "mH_max": safe_stat(df, "mH_target", "max"),
            "total_width_min": safe_stat(df, "total_width", "min"),
            "total_width_max": safe_stat(df, "total_width", "max"),
            "ctau_min": safe_stat(df, "ctau_mm", "min"),
            "ctau_max": safe_stat(df, "ctau_mm", "max"),
            "br_gaga_max": safe_stat(df, "br_gaga", "max"),
            "br_bb_max": safe_stat(df, "br_bb", "max"),
            "br_Zga_max": safe_stat(df, "br_Zga", "max"),
            "physical_filter_status": phys_status,
            "physical_filter_missing_columns": ";".join(phys_missing),
            "mAeqmHp_check_status": mhp_status,
            "mAeqmHp_max_abs_diff": mhp_diff,
            "created_utc": created,
            "input_parquet": str(input_path),
            "selected_values_json": str(selected_path),
        })
        print(
            f"[slices] mA={fmt(combo['mA_target'])} "
            f"lambda6={fmt(combo['lambda6'])} tan_beta={fmt(combo['tan_beta'])}"
            f" -> {df.height} rows ({slice_path})"
        )

    manifest_rows.sort(key=lambda r: r["slice_path"])
    manifest_csv = group_dir / "manifest.csv"
    manifest_json = group_dir / "manifest.json"
    with manifest_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerows(manifest_rows)
    with manifest_json.open("w") as f:
        json.dump(
            {
                "br_channels_missing_columns": br_missing,
                "rows": manifest_rows,
            },
            f, indent=2, sort_keys=True, default=str,
        )
        f.write("\n")

    print(f"[slices] wrote {manifest_csv}")
    print(f"[slices] wrote {manifest_json}")
    if br_missing:
        print(f"[slices] missing width columns (BRs skipped): {br_missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
