#!/usr/bin/env python3
"""Report the real Parquet inventory for the focused-grid pipeline.

Reads ``silver_all.parquet`` with Polars lazy scanning and writes a
deterministic ``inventory.json`` / ``inventory.md`` pair describing the unique
values (and, for high-cardinality columns, summary statistics) of the grid
columns.

A key, campaign-specific fact this report surfaces explicitly is that the grid
is *diagonal*: ``mA_target`` is bound one-to-one to a single ``tan_beta`` and a
single ``lambda6``.  The co-occurrence table makes that visible so downstream
selection does not assume a full cross product.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import polars as pl

# Columns whose full unique value list is small enough to dump verbatim.
LOW_CARD_COLS = ["mA_target", "tan_beta", "lambda6", "lambda7"]
# Columns reported as count/min/max/sample only (never dump millions of values).
HIGH_CARD_COLS = ["lam1", "variation_idx", "mH_target"]
SAMPLE_N = 10


def summarize_column(lf: pl.LazyFrame, col: str, dump_values: bool) -> dict:
    n_unique = int(lf.select(pl.col(col).n_unique()).collect().item())
    cmin = lf.select(pl.col(col).min()).collect().item()
    cmax = lf.select(pl.col(col).max()).collect().item()
    info: dict = {
        "n_unique": n_unique,
        "min": cmin,
        "max": cmax,
    }
    if dump_values:
        values = lf.select(pl.col(col)).unique().collect()[col].sort().to_list()
        info["values"] = values
    else:
        sample = (
            lf.select(pl.col(col))
            .unique()
            .sort(col)
            .head(SAMPLE_N)
            .collect()[col]
            .to_list()
        )
        info["sample"] = sample
    return info


def cooccurrence_table(lf: pl.LazyFrame) -> list[dict]:
    cols = ["mA_target", "tan_beta", "lambda6", "lambda7"]
    counts = (
        lf.group_by(cols)
        .agg(pl.len().alias("n_rows"))
        .sort(cols)
        .collect()
    )
    return counts.to_dicts()


def build_inventory(input_path: Path) -> dict:
    lf = pl.scan_parquet(input_path)
    columns: dict = {}
    for col in LOW_CARD_COLS:
        columns[col] = summarize_column(lf, col, dump_values=True)
    for col in HIGH_CARD_COLS:
        columns[col] = summarize_column(lf, col, dump_values=False)
    total_rows = int(lf.select(pl.len()).collect().item())
    return {
        "input_parquet": str(input_path),
        "total_rows": total_rows,
        "columns": columns,
        "cooccurrence_mA_target_tan_beta_lambda6_lambda7": cooccurrence_table(lf),
    }


def render_markdown(inv: dict) -> str:
    lines: list[str] = []
    lines.append("# Focused-grid Parquet inventory")
    lines.append("")
    lines.append(f"- input: `{inv['input_parquet']}`")
    lines.append(f"- total_rows: {inv['total_rows']}")
    lines.append("")
    lines.append("## Column summaries")
    lines.append("")
    lines.append("| column | n_unique | min | max | values / sample |")
    lines.append("|---|---|---|---|---|")
    for col, info in inv["columns"].items():
        shown = info.get("values", info.get("sample"))
        suffix = "" if "values" in info else " (sample)"
        lines.append(
            f"| `{col}` | {info['n_unique']} | {info['min']} | {info['max']} | "
            f"{shown}{suffix} |"
        )
    lines.append("")
    lines.append("## Co-occurring (mA_target, tan_beta, lambda6, lambda7) combos")
    lines.append("")
    lines.append(
        "The grid is diagonal: each `mA_target` pairs with exactly one "
        "`tan_beta` and one `lambda6`."
    )
    lines.append("")
    lines.append("| mA_target | tan_beta | lambda6 | lambda7 | n_rows |")
    lines.append("|---|---|---|---|---|")
    for row in inv["cooccurrence_mA_target_tan_beta_lambda6_lambda7"]:
        lines.append(
            f"| {row['mA_target']} | {row['tan_beta']} | {row['lambda6']} | "
            f"{row['lambda7']} | {row['n_rows']} |"
        )
    lines.append("")
    return "\n".join(lines)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Path to silver_all.parquet")
    p.add_argument(
        "--output-dir",
        required=True,
        help="Directory for inventory.json and inventory.md",
    )
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    inv = build_inventory(input_path)

    json_path = output_dir / "inventory.json"
    md_path = output_dir / "inventory.md"
    with json_path.open("w") as f:
        json.dump(inv, f, indent=2, sort_keys=True)
        f.write("\n")
    md_path.write_text(render_markdown(inv))

    print(f"[inventory] total_rows={inv['total_rows']}")
    print(f"[inventory] wrote {json_path}")
    print(f"[inventory] wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
