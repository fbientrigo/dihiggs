#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import polars as pl

LOG_FILE: Optional[Path] = None


def append_log(msg: str = "") -> None:
    print(msg)
    if LOG_FILE is not None:
        with open(LOG_FILE, "a", encoding="utf-8") as f:
            f.write(msg + "\n")


NUMERIC_DTYPES = {
    pl.Int8, pl.Int16, pl.Int32, pl.Int64,
    pl.UInt8, pl.UInt16, pl.UInt32, pl.UInt64,
    pl.Float32, pl.Float64, pl.Decimal,
}

FLAG_ALIASES = {
    "positivity_ok": ["positivity_ok"],
    "unitarity_ok": ["unitarity_ok"],
    "perturbativity_ok": ["perturbativity_ok"],
}

KEY_ALIASES = {
    "m_phi": ["m_phi", "mphi"],
    "mA": ["mA", "ma"],
    "alpha": ["alpha"],
    "beta": ["beta"],
    "lambda6": ["lambda6", "lam6"],
    "lambda7": ["lambda7", "lam7"],
    "m12": ["m12"],
    "sin_ba": ["sin_ba", "sba"],
    "tan_beta": ["tan_beta", "tanb"],
    "width_bb": ["width_bb", "w_h2_bb"],
    "width_tautau": ["width_tautau", "w_h2_tautau"],
    "width_WW": ["width_WW", "w_h2_WW", "width_h2_WW"],
    "width_ZZ": ["width_ZZ", "w_h2_ZZ", "width_h2_ZZ"],
    "width_gaga": ["width_gaga", "w_h2_gaga", "width_hgaga"],
    "width_Zga": ["width_Zga", "w_h2_Zga"],
    "width_gg": ["width_gg", "w_h2_gg"],
    "width_hh": ["width_hh", "w_h2_hh"],
    "total_width": ["total_width", "w_total_h2", "total_decay_width", "decay_width"],
    "br_gaga": ["br_gaga", "branching_ratio_hgaga", "branching_ratio_h2_gaga"],
    "lam1": ["lam1", "lambda1"],
    "computed_lam1": ["computed_lam1"],
    "lam2": ["lam2", "lambda2"],
    "computed_lam2": ["computed_lam2"],
    "source_file": ["source_file", "csv_file", "source_csv", "input_file", "file_name", "filename"],
    "campaign": ["campaign", "scan_id", "batch_id", "run_id"],
}

PARAMETER_CANDIDATES = [
    "m_phi", "mA", "alpha", "beta", "lambda6", "lambda7", "m12", "sin_ba", "tan_beta"
]

WIDTH_CANDIDATES = [
    "width_bb", "width_tautau", "width_WW", "width_ZZ", "width_gaga", "width_Zga", "width_gg", "width_hh"
]

DEFAULT_PRIORITY_COLUMNS = [
    "positivity_ok", "unitarity_ok", "perturbativity_ok",
    "m_phi", "mA", "alpha", "beta", "lambda6", "lambda7", "m12", "sin_ba", "tan_beta",
    "width_bb", "width_tautau", "width_WW", "width_ZZ", "width_gaga", "width_Zga", "width_gg", "width_hh",
    "total_width", "br_gaga", "lam1", "computed_lam1", "lam2", "computed_lam2",
]


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def collect_lf(lf: pl.LazyFrame) -> pl.DataFrame:
    try:
        return lf.collect(engine="streaming")
    except TypeError:
        try:
            return lf.collect(streaming=True)
        except TypeError:
            return lf.collect()


def resolve_col(schema_names: List[str], aliases: List[str]) -> Optional[str]:
    lower_map = {c.lower(): c for c in schema_names}
    for a in aliases:
        if a.lower() in lower_map:
            return lower_map[a.lower()]
    for c in schema_names:
        c_low = c.lower()
        for a in aliases:
            if a.lower() in c_low:
                return c
    return None


def build_colmap(schema_names: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for logical, aliases in {**FLAG_ALIASES, **KEY_ALIASES}.items():
        col = resolve_col(schema_names, aliases)
        if col is not None:
            out[logical] = col
    return out


def is_float_dtype(dt: pl.DataType) -> bool:
    return dt in (pl.Float32, pl.Float64)


def is_numeric_dtype(dt: pl.DataType) -> bool:
    return dt in NUMERIC_DTYPES


def fmt_num(x) -> str:
    if x is None:
        return "None"
    if isinstance(x, float):
        if math.isnan(x):
            return "NaN"
        if math.isinf(x):
            return "+Inf" if x > 0 else "-Inf"
        ax = abs(x)
        if ax != 0 and (ax < 1e-4 or ax >= 1e4):
            return f"{x:.6e}"
        return f"{x:.12g}"
    return str(x)


def safe_rate(n: float, d: float) -> float:
    return float(n) / float(d) if d else 0.0


def fmt_pct(x: float) -> str:
    return f"{100.0 * x:.6f}%"


def flag_true_expr(col_name: str) -> pl.Expr:
    as_num = pl.col(col_name).cast(pl.Float64, strict=False)
    as_txt = (
        pl.col(col_name)
        .cast(pl.Utf8, strict=False)
        .str.strip_chars()
        .str.to_lowercase()
    )
    return (
        (as_num >= 0.5) |
        as_txt.is_in(["1", "1.0", "1.000000000000000", "true", "t", "yes", "y"])
    )


def pick_columns(schema: Dict[str, pl.DataType], colmap: Dict[str, str], mode: str, explicit_columns: str) -> List[str]:
    numeric_cols = [c for c, dt in schema.items() if is_numeric_dtype(dt)]
    if explicit_columns.strip().lower() != "all":
        requested = [x.strip() for x in explicit_columns.split(",") if x.strip()]
        return [c for c in requested if c in numeric_cols]

    if mode == "safe":
        picked = []
        for logical in DEFAULT_PRIORITY_COLUMNS:
            if logical in colmap and colmap[logical] in numeric_cols:
                picked.append(colmap[logical])
        return picked
    return numeric_cols


def numeric_column_stats_light(lf: pl.LazyFrame, col: str, include_quantiles: bool = False) -> Dict[str, object]:
    schema = lf.collect_schema()
    dt = schema[col]
    expr = pl.col(col)
    is_float = is_float_dtype(dt)

    finite_expr = expr.is_not_null()
    if is_float:
        finite_expr = finite_expr & ~expr.is_nan() & ~expr.is_infinite()

    select_exprs = [
        pl.len().alias("rows"),
        expr.is_null().sum().alias("nulls"),
        (expr.is_nan().sum().alias("nans") if is_float else pl.lit(0).alias("nans")),
        (expr.is_infinite().sum().alias("infs") if is_float else pl.lit(0).alias("infs")),
        (finite_expr & (expr < 0)).sum().alias("negative"),
        (finite_expr & (expr == 0)).sum().alias("zero"),
        (finite_expr & (expr > 0)).sum().alias("positive"),
        expr.filter(finite_expr).min().alias("min"),
        expr.filter(finite_expr).max().alias("max"),
        expr.filter(finite_expr).mean().alias("mean"),
        expr.filter(finite_expr).std().alias("std"),
    ]
    if include_quantiles:
        select_exprs.extend([
            expr.filter(finite_expr).quantile(0.01).alias("q01"),
            expr.filter(finite_expr).quantile(0.50).alias("q50"),
            expr.filter(finite_expr).quantile(0.99).alias("q99"),
        ])

    res = collect_lf(lf.select(select_exprs)).row(0, named=True)
    rows = int(res["rows"])
    valid = rows - int(res["nulls"]) - int(res["nans"]) - int(res["infs"])
    res["valid"] = valid
    res["valid_frac"] = safe_rate(valid, rows)
    res["null_frac"] = safe_rate(res["nulls"], rows)
    res["nan_frac"] = safe_rate(res["nans"], rows)
    res["inf_frac"] = safe_rate(res["infs"], rows)
    res["negative_frac"] = safe_rate(res["negative"], rows)
    res["zero_frac"] = safe_rate(res["zero"], rows)
    res["positive_frac"] = safe_rate(res["positive"], rows)
    res["dtype"] = str(dt)
    return res


def flag_column_stats(lf: pl.LazyFrame, col: str) -> Dict[str, object]:
    expr_true = flag_true_expr(col)
    expr_num = pl.col(col).cast(pl.Float64, strict=False)
    expr_txt = (
        pl.col(col)
        .cast(pl.Utf8, strict=False)
        .str.strip_chars()
        .str.to_lowercase()
    )
    known_boolish = expr_txt.is_in([
        "0", "0.0", "1", "1.0", "true", "false", "t", "f", "yes", "no", "y", "n", "none", "null", "nan"
    ]) | expr_num.is_not_null()

    res = collect_lf(lf.select(
        pl.len().alias("rows"),
        pl.col(col).is_null().sum().alias("nulls"),
        expr_true.sum().alias("true_like"),
        (~expr_true & pl.col(col).is_not_null()).sum().alias("not_true_like"),
        (~known_boolish & pl.col(col).is_not_null()).sum().alias("weird_tokens"),
        expr_num.min().alias("numeric_min"),
        expr_num.max().alias("numeric_max"),
    )).row(0, named=True)

    rows = int(res["rows"])
    res["true_frac"] = safe_rate(res["true_like"], rows)
    res["null_frac"] = safe_rate(res["nulls"], rows)
    res["weird_frac"] = safe_rate(res["weird_tokens"], rows)
    return res


def check_basic_consistency(lf: pl.LazyFrame, colmap: Dict[str, str], out_dir: Path, sample_rows: int) -> Tuple[List[str], List[Dict[str, object]]]:
    findings: List[str] = []
    summary_rows: List[Dict[str, object]] = []

    schema = lf.collect_schema()
    total_col = colmap.get("total_width")
    br_col = colmap.get("br_gaga")
    gaga_col = colmap.get("width_gaga")
    lam1 = colmap.get("lam1")
    clam1 = colmap.get("computed_lam1")
    lam2 = colmap.get("lam2")
    clam2 = colmap.get("computed_lam2")
    width_cols = [colmap[k] for k in WIDTH_CANDIDATES if k in colmap]

    if br_col:
        r = collect_lf(lf.select(
            pl.len().alias("rows"),
            ((pl.col(br_col) < 0) & pl.col(br_col).is_not_null()).sum().alias("neg_br"),
            ((pl.col(br_col) > 1) & pl.col(br_col).is_not_null()).sum().alias("gt1_br"),
            pl.col(br_col).is_null().sum().alias("null_br"),
            (pl.col(br_col).is_nan().sum().alias("nan_br") if schema[br_col] in (pl.Float32, pl.Float64) else pl.lit(0).alias("nan_br")),
        )).row(0, named=True)
        if r["neg_br"] or r["gt1_br"] or r["null_br"] or r["nan_br"]:
            findings.append(f"BR anomalies in '{br_col}': neg={r['neg_br']}, gt1={r['gt1_br']}, null={r['null_br']}, nan={r['nan_br']}")
        summary_rows.append({"check": "br_bounds", **r, "column": br_col})

    if total_col:
        r = collect_lf(lf.select(
            pl.len().alias("rows"),
            ((pl.col(total_col) <= 0) & pl.col(total_col).is_not_null()).sum().alias("non_positive_total_width"),
            pl.col(total_col).is_null().sum().alias("null_total_width"),
            (pl.col(total_col).is_nan().sum().alias("nan_total_width") if schema[total_col] in (pl.Float32, pl.Float64) else pl.lit(0).alias("nan_total_width")),
        )).row(0, named=True)
        if r["non_positive_total_width"] or r["null_total_width"] or r["nan_total_width"]:
            findings.append(f"Total width anomalies in '{total_col}': non_positive={r['non_positive_total_width']}, null={r['null_total_width']}, nan={r['nan_total_width']}")
        summary_rows.append({"check": "total_width_validity", **r, "column": total_col})

    for wc in width_cols:
        r = collect_lf(lf.select(
            pl.len().alias("rows"),
            ((pl.col(wc) < 0) & pl.col(wc).is_not_null()).sum().alias("negative_width"),
            pl.col(wc).is_null().sum().alias("null_width"),
            (pl.col(wc).is_nan().sum().alias("nan_width") if schema[wc] in (pl.Float32, pl.Float64) else pl.lit(0).alias("nan_width")),
        )).row(0, named=True)
        if r["negative_width"] or r["null_width"] or r["nan_width"]:
            findings.append(f"Width anomalies in '{wc}': neg={r['negative_width']}, null={r['null_width']}, nan={r['nan_width']}")
        summary_rows.append({"check": "width_validity", **r, "column": wc})

    if total_col and width_cols:
        expr_sum = None
        for c in width_cols:
            expr_sum = pl.col(c).fill_null(0) if expr_sum is None else (expr_sum + pl.col(c).fill_null(0))
        r = collect_lf(lf.select(
            pl.len().alias("rows"),
            ((expr_sum > pl.col(total_col)) & pl.col(total_col).is_not_null()).sum().alias("sum_partial_gt_total"),
        )).row(0, named=True)
        if r["sum_partial_gt_total"]:
            findings.append(f"Partial-width sum exceeds total_width in {r['sum_partial_gt_total']} rows.")
        summary_rows.append({"check": "partial_sum_vs_total", **r, "column": total_col})

    if total_col and gaga_col and br_col:
        tol = 1e-8
        ratio_expr = pl.when(pl.col(total_col) > 1e-30).then(pl.col(gaga_col) / pl.col(total_col)).otherwise(None)
        diff_expr = (pl.col(br_col) - ratio_expr).abs()
        r = collect_lf(lf.select(
            pl.len().alias("rows"),
            ((pl.col(total_col) > 1e-30) & diff_expr.is_not_null() & (diff_expr > tol)).sum().alias("br_ratio_mismatch"),
            diff_expr.max().alias("max_abs_diff"),
        )).row(0, named=True)
        if r["br_ratio_mismatch"]:
            findings.append(f"br_gaga != width_gaga/total_width above tol={tol} in {r['br_ratio_mismatch']} rows; max diff={fmt_num(r['max_abs_diff'])}.")
            mism = collect_lf(
                lf.filter((pl.col(total_col) > 1e-30) & diff_expr.is_not_null() & (diff_expr > tol))
                .select([c for c in [br_col, gaga_col, total_col] if c is not None] + [diff_expr.alias("abs_diff")])
                .head(sample_rows)
            )
            if mism.height > 0:
                mism.write_csv(out_dir / "samples_br_ratio_mismatch.csv")
        summary_rows.append({"check": "br_ratio_consistency", **r, "column": br_col})

    for a, b, name in [(lam1, clam1, "lam1_vs_computed_lam1"), (lam2, clam2, "lam2_vs_computed_lam2")]:
        if a and b:
            diff = (pl.col(a) - pl.col(b)).abs()
            r = collect_lf(lf.select(
                pl.len().alias("rows"),
                diff.max().alias("max_abs_diff"),
                (diff > 1e-8).sum().alias("mismatch_gt_1e8"),
            )).row(0, named=True)
            if r["mismatch_gt_1e8"]:
                findings.append(f"{name}: {r['mismatch_gt_1e8']} rows differ by more than 1e-8; max diff={fmt_num(r['max_abs_diff'])}.")
            summary_rows.append({"check": name, **r, "column": a})

    return findings, summary_rows


def suspicious_row_filters(colmap: Dict[str, str]) -> List[Tuple[str, pl.Expr]]:
    out: List[Tuple[str, pl.Expr]] = []
    total_col = colmap.get("total_width")
    br_col = colmap.get("br_gaga")
    gaga_col = colmap.get("width_gaga")
    width_cols = [colmap[k] for k in WIDTH_CANDIDATES if k in colmap]

    if total_col:
        out.append(("total_width_non_positive", pl.col(total_col).is_not_null() & (pl.col(total_col) <= 0)))
    if br_col:
        out.append(("br_negative", pl.col(br_col).is_not_null() & (pl.col(br_col) < 0)))
        out.append(("br_gt_one", pl.col(br_col).is_not_null() & (pl.col(br_col) > 1)))
    for wc in width_cols:
        out.append((f"{wc}_negative", pl.col(wc).is_not_null() & (pl.col(wc) < 0)))
    if total_col and gaga_col and br_col:
        ratio_expr = pl.when(pl.col(total_col) > 1e-30).then(pl.col(gaga_col) / pl.col(total_col)).otherwise(None)
        out.append(("br_ratio_mismatch", (pl.col(total_col) > 1e-30) & ((pl.col(br_col) - ratio_expr).abs() > 1e-8)))
    return out


def export_suspicious_samples(lf: pl.LazyFrame, colmap: Dict[str, str], out_dir: Path, sample_rows: int) -> List[str]:
    exported: List[str] = []
    cols_to_show = [
        c for logical, c in colmap.items()
        if logical in PARAMETER_CANDIDATES + WIDTH_CANDIDATES + ["total_width", "br_gaga", "lam1", "computed_lam1", "lam2", "computed_lam2"]
    ]
    for label, expr in suspicious_row_filters(colmap):
        out_path = out_dir / f"samples_{label}.csv"
        try:
            df = collect_lf(lf.filter(expr).select(cols_to_show).head(sample_rows))
            if df.height > 0:
                df.write_csv(out_path)
                exported.append(str(out_path))
        except Exception as e:
            append_log(f"[!] Could not export sample rows for {label}: {e}")
    return exported


def detect_group_col(colmap: Dict[str, str], explicit_group_by: Optional[str]) -> Optional[str]:
    if explicit_group_by:
        return explicit_group_by
    if "source_file" in colmap:
        return colmap["source_file"]
    if "campaign" in colmap:
        return colmap["campaign"]
    return None


def audit_groups_in_parquet(lf: pl.LazyFrame, group_col: str, colmap: Dict[str, str], out_dir: Path, top_k: int) -> Tuple[Optional[pl.DataFrame], List[str]]:
    findings: List[str] = []
    exprs: List[pl.Expr] = [pl.len().alias("rows")]

    for logical in ["positivity_ok", "unitarity_ok", "perturbativity_ok"]:
        if logical in colmap:
            exprs.append(flag_true_expr(colmap[logical]).sum().alias(f"{logical}_true"))

    if "total_width" in colmap:
        c = colmap["total_width"]
        exprs.extend([
            ((pl.col(c) <= 0) & pl.col(c).is_not_null()).sum().alias("non_positive_total_width"),
            pl.col(c).is_null().sum().alias("null_total_width"),
        ])

    if "br_gaga" in colmap:
        c = colmap["br_gaga"]
        exprs.extend([
            ((pl.col(c) < 0) & pl.col(c).is_not_null()).sum().alias("neg_br"),
            ((pl.col(c) > 1) & pl.col(c).is_not_null()).sum().alias("gt1_br"),
            pl.col(c).is_null().sum().alias("null_br"),
        ])

    agg = lf.group_by(group_col).agg(exprs)
    df = collect_lf(agg)
    if df.height == 0:
        return None, findings

    score_expr = None
    for name in df.columns:
        if name in {group_col, "rows"} or name.endswith("_true"):
            continue
        cur = pl.col(name).fill_null(0).cast(pl.Float64, strict=False)
        score_expr = cur if score_expr is None else (score_expr + cur)

    if score_expr is None:
        df = df.with_columns(pl.lit(0.0).alias("suspicion_score"))
    else:
        df = df.with_columns(score_expr.alias("suspicion_score"))

    df = df.sort(["suspicion_score", "rows"], descending=[True, True])
    df.head(top_k).write_parquet(out_dir / "group_audit_topk.parquet")

    top = df.head(top_k)
    for row in top.iter_rows(named=True):
        if row["suspicion_score"] > 0:
            findings.append(f"Group '{row[group_col]}' has suspicion_score={row['suspicion_score']} over {row['rows']} rows.")
    return top, findings


def write_summary_file(summary_path: Path, text: str) -> None:
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(text)
        if not text.endswith("\n"):
            f.write("\n")


def summarize_findings(findings: List[str], max_items: int = 25) -> str:
    if not findings:
        return "No strong anomalies were found by the implemented checks. That does not prove the lake is correct; it only means the current rules did not catch anything obvious."
    lines = ["Main findings:"]
    for i, item in enumerate(findings[:max_items], start=1):
        lines.append(f"{i}. {item}")
    if len(findings) > max_items:
        lines.append(f"... and {len(findings) - max_items} more findings.")
    return "\n".join(lines)


def main() -> None:
    global LOG_FILE

    parser = argparse.ArgumentParser(description="Memory-safe audit for a consolidated DiHiggs lake.")
    parser.add_argument("--parquet", required=True, help="Path to consolidated parquet/lake file.")
    parser.add_argument("--out-dir", default="lake_audit_logs_safe", help="Directory for logs and outputs.")
    parser.add_argument("--columns", default="all", help="Comma-separated numeric columns to audit, or 'all'.")
    parser.add_argument("--mode", choices=["safe", "deep"], default="safe", help="safe limits work to priority columns and disables heavy passes by default.")
    parser.add_argument("--include-quantiles", action="store_true", help="Enable q01/q50/q99 per column. More expensive.")
    parser.add_argument("--sample-rows", type=int, default=8, help="Rows to export per suspicious class.")
    parser.add_argument("--subset-br", type=float, default=None, help="Optional BR threshold for subset audit.")
    parser.add_argument("--subset-dw", type=float, default=None, help="Optional width threshold for subset audit.")
    parser.add_argument("--enable-group-audit", action="store_true", help="Enable provenance/campaign grouping inside parquet.")
    parser.add_argument("--group-by", default=None, help="Optional parquet column to use as provenance/campaign grouping.")
    parser.add_argument("--group-topk", type=int, default=100, help="Keep only top-K suspicious groups on disk.")
    parser.add_argument("--enable-duplicate-audit", action="store_true", help="Enable exact duplicate-like audit on parameter columns. Expensive.")
    args = parser.parse_args()

    parquet_path = Path(args.parquet).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    ensure_dir(out_dir)
    LOG_FILE = out_dir / "audit_log.txt"
    if LOG_FILE.exists():
        LOG_FILE.unlink()

    append_log("============================================================")
    append_log(" Memory-safe Audit of Consolidated DiHiggs Lake")
    append_log("============================================================")
    append_log(f"[*] Timestamp   : {datetime.now().isoformat(timespec='seconds')}")
    append_log(f"[*] Parquet     : {parquet_path}")
    append_log(f"[*] Output dir  : {out_dir}")
    append_log(f"[*] Mode        : {args.mode}")

    if not parquet_path.exists():
        append_log(f"[!] Parquet not found: {parquet_path}")
        raise SystemExit(1)

    lf = pl.scan_parquet(str(parquet_path))
    schema = lf.collect_schema()
    schema_names = schema.names()
    colmap = build_colmap(schema_names)

    append_log(f"[*] Schema read successfully. Columns found: {len(schema_names)}")
    append_log("[*] Resolved logical columns:")
    for logical in sorted(colmap):
        append_log(f"    - {logical:18s} -> {colmap[logical]}")

    row_count = collect_lf(lf.select(pl.len())).item()
    append_log(f"[*] Total rows in lake: {row_count}")

    cols_to_audit = pick_columns(schema, colmap, args.mode, args.columns)
    append_log(f"[*] Numeric columns selected for audit: {len(cols_to_audit)}")
    if args.mode == "safe" and args.columns.strip().lower() == "all":
        append_log("[*] Safe mode active: auditing only priority columns to reduce full-table rescans.")

    findings: List[str] = []
    column_records: List[Dict[str, object]] = []

    append_log("\n============================================================")
    append_log(" 1) Flag audit")
    append_log("============================================================")
    for logical in ["positivity_ok", "unitarity_ok", "perturbativity_ok"]:
        if logical not in colmap:
            append_log(f"[!] Flag column not found: {logical}")
            findings.append(f"Missing mandatory flag column: {logical}")
            continue
        stats = flag_column_stats(lf, colmap[logical])
        append_log(
            f"[*] {logical} -> rows={stats['rows']}, true={stats['true_like']} ({fmt_pct(stats['true_frac'])}), "
            f"null={stats['nulls']} ({fmt_pct(stats['null_frac'])}), weird_tokens={stats['weird_tokens']} ({fmt_pct(stats['weird_frac'])}), "
            f"numeric_min={fmt_num(stats['numeric_min'])}, numeric_max={fmt_num(stats['numeric_max'])}"
        )
        if stats["weird_tokens"]:
            findings.append(f"Flag column '{colmap[logical]}' contains {stats['weird_tokens']} weird/non-boolish tokens.")

    found_flags = [colmap[k] for k in ["positivity_ok", "unitarity_ok", "perturbativity_ok"] if k in colmap]
    if found_flags:
        combined = None
        for c in found_flags:
            expr = flag_true_expr(c)
            combined = expr if combined is None else (combined & expr)
        cstats = collect_lf(lf.select(pl.len().alias("rows"), combined.sum().alias("all_flags_true"))).row(0, named=True)
        append_log(
            f"[*] Combined theoretical mask: {cstats['all_flags_true']} / {cstats['rows']} rows pass all mandatory flags ({fmt_pct(safe_rate(cstats['all_flags_true'], cstats['rows']))})."
        )

    append_log("\n============================================================")
    append_log(" 2) Numeric column audit")
    append_log("============================================================")
    for i, c in enumerate(cols_to_audit, start=1):
        append_log(f"  -> [{i}/{len(cols_to_audit)}] Auditing column: {c}")
        stats = numeric_column_stats_light(lf, c, include_quantiles=args.include_quantiles)
        column_records.append({"column": c, **stats})
        head = (
            f"     dtype={stats['dtype']}, valid={stats['valid']} ({fmt_pct(stats['valid_frac'])}), "
            f"null={stats['nulls']} ({fmt_pct(stats['null_frac'])}), nan={stats['nans']} ({fmt_pct(stats['nan_frac'])}), "
            f"inf={stats['infs']} ({fmt_pct(stats['inf_frac'])}), neg={stats['negative']} ({fmt_pct(stats['negative_frac'])}), "
            f"zero={stats['zero']} ({fmt_pct(stats['zero_frac'])}), pos={stats['positive']} ({fmt_pct(stats['positive_frac'])})"
        )
        append_log(head)
        tail = f"     min={fmt_num(stats['min'])}, max={fmt_num(stats['max'])}, mean={fmt_num(stats['mean'])}, std={fmt_num(stats['std'])}"
        if args.include_quantiles:
            tail += f", q01={fmt_num(stats['q01'])}, median={fmt_num(stats['q50'])}, q99={fmt_num(stats['q99'])}"
        append_log(tail)

        lname = c.lower()
        if "width" in lname and stats["negative"]:
            findings.append(f"Width column '{c}' has {stats['negative']} negative values.")
        if lname.startswith("br") or "branching" in lname:
            if stats["negative"]:
                findings.append(f"BR-like column '{c}' has {stats['negative']} negative values.")
        if stats["nan_frac"] > 0 or stats["inf_frac"] > 0:
            findings.append(f"Column '{c}' contains NaN/Inf values (nan={stats['nans']}, inf={stats['infs']}).")

    if column_records:
        pl.DataFrame(column_records).write_parquet(out_dir / "column_audit.parquet")

    append_log("\n============================================================")
    append_log(" 3) Cross-column consistency checks")
    append_log("============================================================")
    cross_findings, summary_rows = check_basic_consistency(lf, colmap, out_dir, sample_rows=args.sample_rows)
    findings.extend(cross_findings)
    if summary_rows:
        pl.DataFrame(summary_rows).write_parquet(out_dir / "consistency_checks.parquet")
    if cross_findings:
        for item in cross_findings:
            append_log(f"[!] {item}")
    else:
        append_log("[*] No strong cross-column inconsistency was detected by the implemented rules.")

    append_log("\n============================================================")
    append_log(" 4) Suspicious row samples")
    append_log("============================================================")
    exported = export_suspicious_samples(lf, colmap, out_dir, sample_rows=args.sample_rows)
    if exported:
        for p in exported:
            append_log(f"[*] Exported suspicious sample rows -> {p}")
    else:
        append_log("[*] No suspicious sample files were exported by the current rules.")

    append_log("\n============================================================")
    append_log(" 5) Duplicate-like audit on parameter columns")
    append_log("============================================================")
    if args.enable_duplicate_audit:
        param_cols = [colmap[k] for k in PARAMETER_CANDIDATES if k in colmap]
        if param_cols:
            append_log("[*] Duplicate audit enabled. This may be memory-heavy on very large lakes.")
            dup = collect_lf(lf.select(
                pl.len().alias("rows"),
                pl.struct(param_cols).n_unique().alias("n_unique_param_points")
            )).row(0, named=True)
            n_dup = int(dup["rows"] - dup["n_unique_param_points"])
            append_log(f"[*] Parameter-point uniqueness over columns {param_cols}: unique={dup['n_unique_param_points']}, duplicates={n_dup}")
            if n_dup > 0:
                findings.append(f"Duplicate parameter points detected: {n_dup} repeated rows over parameter columns {param_cols}.")
        else:
            append_log("[!] Could not resolve enough parameter columns for duplicate-like audit.")
    else:
        append_log("[*] Duplicate audit skipped in safe configuration. Enable with --enable-duplicate-audit if you really need it.")

    append_log("\n============================================================")
    append_log(" 6) Group/campaign audit inside parquet")
    append_log("============================================================")
    if args.enable_group_audit:
        group_col = detect_group_col(colmap, args.group_by)
        if group_col and group_col in schema_names:
            append_log(f"[*] Using group column for provenance audit: {group_col}")
            group_df, group_findings = audit_groups_in_parquet(lf, group_col, colmap, out_dir, top_k=args.group_topk)
            findings.extend(group_findings)
            if group_df is not None and group_df.height > 0:
                append_log("[*] Top suspicious groups in parquet:")
                for row in group_df.head(10).iter_rows(named=True):
                    append_log(f"    - {row[group_col]} | rows={row['rows']} | suspicion_score={row['suspicion_score']}")
        else:
            append_log("[*] No provenance/group column detected inside the parquet.")
    else:
        append_log("[*] Group/campaign audit skipped by default to avoid materializing large grouped outputs.")

    append_log("\n============================================================")
    append_log(" 7) Optional subset audit")
    append_log("============================================================")
    if args.subset_br is not None or args.subset_dw is not None:
        subset_lf = lf
        subset_desc = []
        if args.subset_br is not None and "br_gaga" in colmap:
            subset_lf = subset_lf.filter(pl.col(colmap["br_gaga"]) > args.subset_br)
            subset_desc.append(f"{colmap['br_gaga']} > {args.subset_br}")
        if args.subset_dw is not None and "total_width" in colmap:
            subset_lf = subset_lf.filter(pl.col(colmap["total_width"]) < args.subset_dw)
            subset_desc.append(f"{colmap['total_width']} < {args.subset_dw}")

        sub_rows = collect_lf(subset_lf.select(pl.len())).item()
        append_log(f"[*] Subset defined by: {' AND '.join(subset_desc) if subset_desc else '(none)'}")
        append_log(f"[*] Subset rows: {sub_rows}")

        if sub_rows == 0:
            findings.append("Subset audit produced zero rows. Thresholds may be too restrictive or the data may be inconsistent with the expected region.")
        else:
            sub_records: List[Dict[str, object]] = []
            for c in cols_to_audit:
                stats = numeric_column_stats_light(subset_lf, c, include_quantiles=False)
                sub_records.append({"column": c, **stats})
                append_log(
                    f"    - SUBSET {c}: valid={stats['valid']} ({fmt_pct(stats['valid_frac'])}), null={stats['nulls']}, nan={stats['nans']}, inf={stats['infs']}, neg={stats['negative']}, zero={stats['zero']}, pos={stats['positive']}"
                )
            if sub_records:
                pl.DataFrame(sub_records).write_parquet(out_dir / "subset_column_audit.parquet")
            sub_findings, sub_summary_rows = check_basic_consistency(subset_lf, colmap, out_dir, sample_rows=args.sample_rows)
            findings.extend([f"[subset] {x}" for x in sub_findings])
            if sub_summary_rows:
                pl.DataFrame(sub_summary_rows).write_parquet(out_dir / "subset_consistency_checks.parquet")
    else:
        append_log("[*] Subset audit skipped because no subset thresholds were provided.")

    append_log("\n============================================================")
    append_log(" 8) Final summary")
    append_log("============================================================")

    findings = list(dict.fromkeys(findings))
    summary_lines = []
    summary_lines.append("Audit summary")
    summary_lines.append("=" * 60)
    summary_lines.append(f"Parquet: {parquet_path}")
    summary_lines.append(f"Rows: {row_count}")
    summary_lines.append(f"Mode: {args.mode}")
    summary_lines.append(f"Numeric columns audited: {len(cols_to_audit)}")
    summary_lines.append("")
    summary_lines.append(summarize_findings(findings))
    summary_lines.append("")
    summary_lines.append("What to inspect next")
    summary_lines.append("-" * 60)
    summary_lines.append("1. Start with samples_br_ratio_mismatch.csv if it exists.")
    summary_lines.append("2. Then inspect column_audit.parquet, especially columns with NaN/Inf/negative values.")
    summary_lines.append("3. Re-check rows where total_width <= 0 but any partial width is positive.")
    summary_lines.append("4. Only after that, enable group or duplicate audit if you still need campaign-level isolation.")
    summary_lines.append("")
    summary_lines.append("Memory notes")
    summary_lines.append("-" * 60)
    summary_lines.append("- Safe mode skips exact duplicate audit and group materialization by default.")
    summary_lines.append("- Quantiles are off by default because they force heavier passes over very large columns.")
    summary_lines.append("- If WSL is still unstable, audit only a short list of columns with --columns.")
    summary_lines.append("")
    summary_lines.append("Generated artifacts")
    summary_lines.append("-" * 60)
    summary_lines.append(f"- audit_log.txt: {out_dir / 'audit_log.txt'}")
    summary_lines.append(f"- audit_summary.txt: {out_dir / 'audit_summary.txt'}")
    if (out_dir / "column_audit.parquet").exists():
        summary_lines.append(f"- column_audit.parquet: {out_dir / 'column_audit.parquet'}")
    if (out_dir / "consistency_checks.parquet").exists():
        summary_lines.append(f"- consistency_checks.parquet: {out_dir / 'consistency_checks.parquet'}")
    if (out_dir / "group_audit_topk.parquet").exists():
        summary_lines.append(f"- group_audit_topk.parquet: {out_dir / 'group_audit_topk.parquet'}")
    if (out_dir / "subset_column_audit.parquet").exists():
        summary_lines.append(f"- subset_column_audit.parquet: {out_dir / 'subset_column_audit.parquet'}")
    if (out_dir / "subset_consistency_checks.parquet").exists():
        summary_lines.append(f"- subset_consistency_checks.parquet: {out_dir / 'subset_consistency_checks.parquet'}")
    for p in sorted(out_dir.glob("samples_*.csv")):
        summary_lines.append(f"- {p.name}: {p}")

    summary_txt = "\n".join(summary_lines)
    append_log(summary_txt)
    write_summary_file(out_dir / "audit_summary.txt", summary_txt)

    metadata = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "parquet": str(parquet_path),
        "row_count": int(row_count),
        "numeric_columns_audited": cols_to_audit,
        "resolved_columns": colmap,
        "findings": findings,
        "mode": args.mode,
    }
    with open(out_dir / "audit_metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)

    append_log("\n[✓] Audit finished.")


if __name__ == "__main__":
    main()
