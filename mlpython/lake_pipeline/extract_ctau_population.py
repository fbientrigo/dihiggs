#!/usr/bin/env python3
"""Extract the small/big ctau population from silver_all.parquet.

This is the one bounded, single-pass scan over the full 6,000,000-row
``silver_all.parquet`` permitted by project rules (no per-combo looping over
the unfiltered file). It produces three outputs in one go:

1. Per-combo ``slice.parquet`` files (+ manifest) for every (mA_target,
   tan_beta, lambda6, lambda7) combo that has any row with ``ctau_mm >=
   --boundary`` ("big" population). These are shaped exactly like the
   existing focused-grid slices so ``run_focused_grid_plots.py`` can be reused
   unmodified to make BR-vs-ctau plots for them.
2. A combined small+big sample parquet for parallel-coordinates plotting: all
   "big" rows (already bounded) plus a deterministic stride-sampled subset of
   the much larger "small" population.
3. A deterministic JSON summary of exact counts.

Never touches ``data/.../plots_chris/_partial_unbounded_family/`` and never
writes into the existing ``focused_grid/`` output tree.
"""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import polars as pl

# Shared physics constant (formerly a local copy of 1.973269804e-13). Sourced
# from mlpython/lake_pipeline/physics_conventions.py, which reads the
# ecosystem-wide conventions/physics_conventions.yaml. Kept as the local name
# HBARC_GEV_MM so the rest of this script is unchanged.
from physics_conventions import HBAR_C_GEV_MM as HBARC_GEV_MM

PHYS_COLS = ["positivity_ok", "unitarity_ok", "perturbativity_ok"]
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
COMBO_KEY = ["mA_target", "tan_beta", "lambda6", "lambda7"]
PARCOORDS_COLS = [
    "mA_target", "tan_beta", "lambda6", "m_phi", "lam1", "total_width", "br_gaga",
]
MANIFEST_FIELDS = ["slice_path", "mA_target", "tan_beta", "lambda6", "lambda7", "n_rows"]


def fmt(value: float) -> str:
    return f"{value:g}"


def _stream_collect(lf: pl.LazyFrame) -> pl.DataFrame:
    try:
        return lf.collect(engine="streaming")
    except TypeError:
        return lf.collect(streaming=True)
    except Exception:
        return lf.collect()


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


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Path to silver_all.parquet")
    p.add_argument("--output-dir", required=True,
                   help="Base output directory (e.g. .../plots_chris/ctau_population)")
    p.add_argument("--boundary", type=float, default=1.0,
                   help="ctau_mm threshold defining the 'big' population")
    p.add_argument("--small-sample-n", type=int, default=200_000,
                   help="Cap on small-population rows kept for parallel coords")
    p.add_argument("--seed", type=int, default=42,
                   help="Unused for stride sampling but recorded for reproducibility")
    p.add_argument("--apply-phys-filter", action="store_true",
                   help="Apply *_ok==1 filters when those columns are present")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    big_slices_dir = output_dir / "big_slices"
    big_slices_dir.mkdir(parents=True, exist_ok=True)

    lf = pl.scan_parquet(input_path)
    schema_names = set(lf.collect_schema().names())

    br_exprs, br_missing = br_and_ctau_exprs(schema_names)
    base = lf.with_columns(br_exprs)
    base, phys_status, phys_missing = physical_filter(
        base, schema_names, args.apply_phys_filter
    )

    created = datetime.now(timezone.utc).isoformat()

    # --- 1) Big population: bounded, materialize fully ---
    big_lf = base.filter(pl.col("ctau_mm") >= args.boundary)
    big_df = _stream_collect(big_lf)
    n_big = big_df.height
    print(f"[extract] big population (ctau_mm >= {args.boundary}): {n_big} rows")

    manifest_rows: list[dict] = []
    per_combo_counts: list[dict] = []
    if n_big > 0:
        for combo_vals, group_df in big_df.group_by(COMBO_KEY, maintain_order=True):
            combo = dict(zip(COMBO_KEY, combo_vals))
            slice_dir = (
                big_slices_dir
                / f"mA={fmt(combo['mA_target'])}"
                / f"lambda6={fmt(combo['lambda6'])}"
                / f"tan_beta={fmt(combo['tan_beta'])}"
            )
            slice_dir.mkdir(parents=True, exist_ok=True)
            slice_path = slice_dir / "slice.parquet"
            group_df.write_parquet(slice_path)
            manifest_rows.append({
                "slice_path": str(slice_path),
                "mA_target": combo["mA_target"],
                "tan_beta": combo["tan_beta"],
                "lambda6": combo["lambda6"],
                "lambda7": combo["lambda7"],
                "n_rows": group_df.height,
            })
            per_combo_counts.append({**combo, "n_rows": group_df.height})
            print(
                f"[extract] combo mA={fmt(combo['mA_target'])} "
                f"tan_beta={fmt(combo['tan_beta'])} lambda6={fmt(combo['lambda6'])} "
                f"-> {group_df.height} rows ({slice_path})"
            )

    manifest_rows.sort(key=lambda r: r["slice_path"])
    manifest_csv = big_slices_dir / "manifest.csv"
    manifest_json = big_slices_dir / "manifest.json"
    with manifest_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerows(manifest_rows)
    with manifest_json.open("w") as f:
        json.dump(manifest_rows, f, indent=2, sort_keys=True, default=str)
        f.write("\n")
    print(f"[extract] wrote {manifest_csv}")
    print(f"[extract] wrote {manifest_json}")

    # --- 2) Small population: lazy row-index stride sample, never fully materialized ---
    small_lf = base.filter(pl.col("ctau_mm") < args.boundary)
    n_small = int(
        _stream_collect(small_lf.select(pl.len().alias("n")))["n"].item()
    )
    print(f"[extract] small population (ctau_mm < {args.boundary}): {n_small} rows")

    stride = max(1, n_small // max(1, args.small_sample_n))
    small_sampled_lf = (
        small_lf.with_row_index("rid")
        .filter((pl.col("rid") % stride) == 0)
        .limit(args.small_sample_n)
        .select(PARCOORDS_COLS + ["ctau_mm"])
        .with_columns(pl.lit("small").alias("population_class"))
    )
    small_sample_df = _stream_collect(small_sampled_lf)
    n_small_sampled = small_sample_df.height
    print(f"[extract] small population sampled: {n_small_sampled} rows (stride={stride})")

    big_proj_df = (
        big_df.select(PARCOORDS_COLS + ["ctau_mm"])
        .with_columns(pl.lit("big").alias("population_class"))
    )

    combined = pl.concat([big_proj_df, small_sample_df])
    parcoords_path = output_dir / "population_for_parcoords.parquet"
    combined.write_parquet(parcoords_path)
    print(f"[extract] wrote {parcoords_path} ({combined.height} rows)")

    summary = {
        "input_parquet": str(input_path),
        "output_dir": str(output_dir),
        "boundary": args.boundary,
        "apply_phys_filter": bool(args.apply_phys_filter),
        "physical_filter_status": phys_status,
        "physical_filter_missing_columns": phys_missing,
        "br_channels_missing_columns": br_missing,
        "n_big_total": n_big,
        "n_small_total": n_small,
        "n_small_sampled": n_small_sampled,
        "small_sample_stride": stride,
        "small_sample_n_requested": args.small_sample_n,
        "seed": args.seed,
        "per_combo_counts": sorted(per_combo_counts, key=lambda r: r["mA_target"]),
        "manifest_csv": str(manifest_csv),
        "population_for_parcoords": str(parcoords_path),
        "created_utc": created,
    }
    summary_path = output_dir / "extraction_summary.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2, sort_keys=True, default=str)
        f.write("\n")
    print(f"[extract] wrote {summary_path}")

    if n_big + n_small != 6_000_000:
        print(
            f"[extract][WARN] big+small={n_big + n_small} != 6,000,000 total rows "
            "(input schema/row-count may have changed)."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
