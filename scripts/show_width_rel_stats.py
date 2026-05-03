#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

REL_COLS = {
    "bb": "rel_err_width_bb",
    "cc": "rel_err_width_cc",
    "tautau": "rel_err_width_tautau",
    "gg": "rel_err_width_gg",
    "gaga": "rel_err_width_gaga",
    "Zga": "rel_err_width_Zga",
    "total": "rel_err_width_total",
}

DEFAULT_MAX_POINT_COLS = [
    "point_id", "mH", "lambda1", "lambda6", "tan_beta", "triple_ok",
    "rel_err_width_bb", "rel_err_width_cc", "rel_err_width_tautau",
    "rel_err_width_gg", "rel_err_width_gaga", "rel_err_width_Zga", "rel_err_width_total",
]


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for channel, col in REL_COLS.items():
        if col not in df.columns:
            continue
        s = pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan)
        finite = s.dropna()
        rows.append({
            "channel": channel,
            "n_finite": int(finite.shape[0]),
            "mean_rel": float(finite.mean()) if len(finite) else np.nan,
            "median_rel": float(finite.median()) if len(finite) else np.nan,
            "std_rel": float(finite.std()) if len(finite) else np.nan,
            "min_rel": float(finite.min()) if len(finite) else np.nan,
            "max_rel": float(finite.max()) if len(finite) else np.nan,
        })
    return pd.DataFrame(rows)


def select_max_point(df: pd.DataFrame, channel: str) -> pd.Series | None:
    col = REL_COLS.get(channel)
    if col is None or col not in df.columns:
        return None
    s = pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan)
    if s.dropna().empty:
        return None
    return df.loc[s.idxmax()]


def parse_cols(csv_cols: list[str], spec: str | None) -> list[str]:
    if not spec:
        return [c for c in DEFAULT_MAX_POINT_COLS if c in csv_cols]
    cols = [c.strip() for c in spec.split(",") if c.strip()]
    return [c for c in cols if c in csv_cols]


def main() -> int:
    parser = argparse.ArgumentParser(description="Show terminal table with relative-width error stats from merged comparison CSV.")
    parser.add_argument("csv", nargs="?", default="colab_vs_2hdmc_merged_comparison.csv")
    parser.add_argument("--sort-by", default="mean_rel", choices=["channel", "n_finite", "mean_rel", "median_rel", "std_rel", "min_rel", "max_rel"])
    parser.add_argument("--ascending", action="store_true")
    parser.add_argument("--max-point", action="store_true")
    parser.add_argument("--max-point-cols", default=None, help="Comma-separated columns for --max-point row display")
    parser.add_argument("--out-csv", default=None, help="Optional output CSV path for summary table")
    args = parser.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise SystemExit(f"CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    if "triple_ok" in df.columns:
        df = df[df["triple_ok"] == True]

    summary = build_summary(df)
    if summary.empty:
        raise SystemExit("No relative-width columns were found in the CSV.")

    summary = summary.sort_values(args.sort_by, ascending=args.ascending).reset_index(drop=True)

    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 200)
    pd.set_option("display.float_format", lambda x: f"{x:.6e}")

    print("\nRelative width error stats\n")
    print(summary.to_string(index=False))

    if args.out_csv:
        out_csv = Path(args.out_csv)
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out_csv, index=False)
        print(f"\nWrote summary CSV: {out_csv}")

    if args.max_point:
        cols = parse_cols(list(df.columns), args.max_point_cols)
        print("\nPoints with maximum relative width errors\n")
        for channel in REL_COLS:
            max_point = select_max_point(df, channel)
            if max_point is not None:
                print(f"{channel}:")
                row = max_point[cols] if cols else max_point
                print(row.to_string())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
