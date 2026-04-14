#!/usr/bin/env python3

import csv
import math
import statistics
import sys
from pathlib import Path


def quantile(sorted_values: list[float], q: float) -> float:
    if not sorted_values:
        return float("nan")
    if len(sorted_values) == 1:
        return sorted_values[0]
    pos = q * (len(sorted_values) - 1)
    lower = math.floor(pos)
    upper = math.ceil(pos)
    if lower == upper:
        return sorted_values[lower]
    weight = pos - lower
    return sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight


def main() -> int:
    if len(sys.argv) not in {2, 3}:
        print("Usage: ./analyze_lam1_errors.py input_csv [output_md]", file=sys.stderr)
        return 1

    input_csv = Path(sys.argv[1])
    output_md = Path(sys.argv[2]) if len(sys.argv) == 3 else None

    rows: list[dict[str, str]] = []
    with input_csv.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row.get("abs_error") not in {None, ""}]

    if not rows:
        print("Input dataset is empty", file=sys.stderr)
        return 2

    abs_errors = [float(row["abs_error"]) for row in rows]
    warning_rows = [row for row in rows if int(row["warning_flag"]) == 1]
    sorted_errors = sorted(abs_errors)
    max_idx = max(range(len(rows)), key=lambda i: abs_errors[i])
    worst = rows[max_idx]

    lines = [
        f"# lam1 error analysis for `{input_csv.name}`",
        "",
        "## Dataset summary",
        "",
        f"- samples: `{len(rows)}`",
        f"- warnings (`abs_error > 1e-12`): `{len(warning_rows)}`",
        f"- warning rate: `{len(warning_rows) / len(rows):.6f}`",
        "",
        "## Absolute error statistics",
        "",
        f"- min: `{min(abs_errors):.17e}`",
        f"- median: `{statistics.median(abs_errors):.17e}`",
        f"- mean: `{statistics.fmean(abs_errors):.17e}`",
        f"- p90: `{quantile(sorted_errors, 0.90):.17e}`",
        f"- p99: `{quantile(sorted_errors, 0.99):.17e}`",
        f"- max: `{max(abs_errors):.17e}`",
        "",
        "## Worst sample",
        "",
        f"- sample_index: `{worst['sample_index']}`",
        f"- attempt_index: `{worst['attempt_index']}`",
        f"- lambda1_input: `{float(worst['lambda1_input']):.17e}`",
        f"- lambda1_recomputed: `{float(worst['lambda1_recomputed']):.17e}`",
        f"- abs_error: `{float(worst['abs_error']):.17e}`",
        f"- warning_flag: `{worst['warning_flag']}`",
        f"- point: `mh={float(worst['mh']):.8f}, mH={float(worst['mH']):.8f}, mA={float(worst['mA']):.8f}, mHp={float(worst['mHp']):.8f}, sba={float(worst['sba']):.8f}, lambda6={float(worst['lambda6']):.8f}, lambda7={float(worst['lambda7']):.8f}, m12_2={float(worst['m12_2']):.8f}, tan_beta={float(worst['tan_beta']):.8f}`",
        "",
        "## Interpretation",
        "",
        "- This dataset measures the round-trip error introduced by `set_param_phys_lam1(...)` after reconstructing `m12_2` and delegating through `set_param_phys(...)`.",
        "- `warning_flag=1` means the stored absolute difference exceeded the hard threshold `1e-12` used by `THDM::EPS`.",
        "- Compare the max and warning rate to decide whether your preferred scan region needs a looser warning threshold or further numerical stabilization.",
    ]

    report = "\n".join(lines) + "\n"
    if output_md is not None:
        output_md.write_text(report)
    print(report, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
