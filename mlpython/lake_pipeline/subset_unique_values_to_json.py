#!/usr/bin/env python3
import argparse
import json
import math
import os
import sys
from datetime import datetime
from pathlib import Path

import polars as pl


def _resolve_col(schema_names, aliases):
    lower_map = {c.lower(): c for c in schema_names}
    for a in aliases:
        if a.lower() in lower_map:
            return lower_map[a.lower()]
    for c in schema_names:
        cl = c.lower()
        if any(tok.lower() in cl for tok in aliases):
            return c
    return None


def _is_float_dtype(dt) -> bool:
    return dt in (pl.Float32, pl.Float64)


def _json_safe_number(x):
    if x is None:
        return None
    if isinstance(x, bool):
        return x
    if isinstance(x, int):
        return x
    if isinstance(x, float):
        if math.isnan(x) or math.isinf(x):
            return None
        return x
    try:
        xf = float(x)
        if math.isnan(xf) or math.isinf(xf):
            return None
        return xf
    except Exception:
        return str(x)


def _append_log(msg: str, log_path: Path):
    print(msg)
    with log_path.open("a", encoding="utf-8") as f:
        f.write(msg + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Read a subset parquet and export unique values from one column into JSON."
    )
    parser.add_argument(
        "--parquet",
        type=str,
        default="temp_subspace.parquet",
        help="Path to the temporary subset parquet",
    )
    parser.add_argument(
        "--column",
        type=str,
        default="lambda6",
        help="Column to inspect. Default: lambda6",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="subset_unique_values",
        help="Directory where JSON and log will be written.",
    )
    parser.add_argument(
        "--out-json",
        type=str,
        default=None,
        help="Optional explicit JSON output path. Overrides --out-dir filename.",
    )
    parser.add_argument(
        "--sort",
        action="store_true",
        help="Sort unique values before writing JSON.",
    )
    parser.add_argument(
        "--with-counts",
        action="store_true",
        help="Also export occurrence counts for each unique value.",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON.",
    )
    args = parser.parse_args()

    parquet_path = Path(args.parquet)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    output_json = Path(args.out_json) if args.out_json else out_dir / f"unique_{args.column}.json"
    log_path = out_dir / "unique_export_log.txt"
    if log_path.exists():
        log_path.unlink()

    _append_log("============================================================", log_path)
    _append_log(" Unique values export from subset parquet", log_path)
    _append_log("============================================================", log_path)
    _append_log(f"[*] Timestamp : {datetime.now().isoformat(timespec='seconds')}", log_path)
    _append_log(f"[*] Parquet   : {parquet_path}", log_path)
    _append_log(f"[*] Out dir   : {out_dir.resolve()}", log_path)

    if not parquet_path.exists():
        _append_log(f"[!] Parquet not found: {parquet_path}", log_path)
        sys.exit(1)

    lf = pl.scan_parquet(parquet_path)
    schema = lf.collect_schema()
    schema_names = schema.names()
    _append_log(f"[*] Schema read successfully. Columns found: {len(schema_names)}", log_path)

    requested = args.column
    aliases = [requested]
    if requested.lower() in {"lambda6", "lambda_6", "lam6"}:
        aliases = ["lambda6", "lambda_6", "lam6"]
    elif requested.lower() in {"lambda7", "lambda_7", "lam7"}:
        aliases = ["lambda7", "lambda_7", "lam7"]

    target_col = _resolve_col(schema_names, aliases)
    if target_col is None:
        _append_log(f"[!] Could not resolve requested column '{requested}'.", log_path)
        _append_log(f"[*] Available columns: {', '.join(schema_names)}", log_path)
        sys.exit(1)

    dtype = schema[target_col]
    _append_log(f"[*] Requested column : {requested}", log_path)
    _append_log(f"[*] Resolved column  : {target_col}", log_path)
    _append_log(f"[*] Column dtype     : {dtype}", log_path)

    valid_expr = pl.col(target_col).is_not_null()
    nan_expr = pl.lit(False)
    if _is_float_dtype(dtype):
        nan_expr = pl.col(target_col).is_nan()
        valid_expr = valid_expr & ~nan_expr

    stats = lf.select(
        pl.len().alias("rows_total"),
        pl.col(target_col).is_null().sum().alias("null_count"),
        nan_expr.sum().alias("nan_count"),
        valid_expr.sum().alias("valid_count"),
        pl.col(target_col).filter(valid_expr).min().alias("min_value"),
        pl.col(target_col).filter(valid_expr).max().alias("max_value"),
    ).collect()

    rows_total = int(stats["rows_total"][0])
    null_count = int(stats["null_count"][0])
    nan_count = int(stats["nan_count"][0])
    valid_count = int(stats["valid_count"][0])
    min_value = _json_safe_number(stats["min_value"][0])
    max_value = _json_safe_number(stats["max_value"][0])

    _append_log(f"[*] Row count       : {rows_total}", log_path)
    _append_log(f"[*] Valid values    : {valid_count}", log_path)
    _append_log(f"[*] Null values     : {null_count}", log_path)
    _append_log(f"[*] NaN values      : {nan_count}", log_path)
    _append_log(f"[*] Min / Max       : {min_value} / {max_value}", log_path)
    _append_log(f"[*] Extracting unique values from '{target_col}'...", log_path)

    unique_lf = (
        lf.select(pl.col(target_col))
        .filter(valid_expr)
        .unique()
    )

    if args.sort:
        unique_lf = unique_lf.sort(target_col)

    try:
        unique_df = unique_lf.collect(streaming=True)
    except TypeError:
        unique_df = unique_lf.collect()

    unique_values = unique_df[target_col].to_list()
    unique_values_json = [_json_safe_number(v) for v in unique_values]
    unique_count = len(unique_values_json)

    _append_log(f"[+] Unique values found: {unique_count}", log_path)

    payload = {
        "metadata": {
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "parquet_path": str(parquet_path.resolve()),
            "requested_column": requested,
            "resolved_column": target_col,
            "dtype": str(dtype),
            "rows_total": rows_total,
            "valid_count": valid_count,
            "null_count": null_count,
            "nan_count": nan_count,
            "min_value": min_value,
            "max_value": max_value,
            "unique_count": unique_count,
            "sorted": bool(args.sort),
            "with_counts": bool(args.with_counts),
        },
        "unique_values": unique_values_json,
    }

    if args.with_counts:
        _append_log(f"[*] Computing per-value counts for '{target_col}'...", log_path)
        counts_lf = (
            lf.select(pl.col(target_col))
            .filter(valid_expr)
            .group_by(target_col)
            .len()
        )
        if args.sort:
            counts_lf = counts_lf.sort(target_col)
        try:
            counts_df = counts_lf.collect(streaming=True)
        except TypeError:
            counts_df = counts_lf.collect()

        value_counts = []
        for row in counts_df.iter_rows(named=True):
            value_counts.append(
                {
                    "value": _json_safe_number(row[target_col]),
                    "count": int(row["len"]),
                }
            )
        payload["value_counts"] = value_counts
        _append_log(f"[+] Exported counts for {len(value_counts)} unique values.", log_path)

    with output_json.open("w", encoding="utf-8") as f:
        if args.pretty:
            json.dump(payload, f, ensure_ascii=False, indent=2)
        else:
            json.dump(payload, f, ensure_ascii=False)

    _append_log(f"[+] JSON written to: {output_json.resolve()}", log_path)
    _append_log(f"[+] Log written to : {log_path.resolve()}", log_path)
    _append_log("[✓] Done.", log_path)


if __name__ == "__main__":
    main()
