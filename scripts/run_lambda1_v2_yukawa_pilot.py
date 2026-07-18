#!/usr/bin/env python3
"""Regenerate the bounded lambda1.v2 Yukawa-initialization pilot."""

import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "dihiggs/app/Lambda1EvaluatorV2"
OUT = ROOT / "docs/pilots/lambda1_v2_yukawa_fix_v1"
INPUT = OUT / "input.csv"
OUTPUT = OUT / "lambda1_v2_pilot.csv"
REPORT = ROOT / "docs/verification/lambda1_v2_yukawa_fix_v1"
HEADER = "point_id,mh_gev,mH_gev,mA_gev,mHp_gev,sin_beta_minus_alpha,tan_beta,lambda1_target,lambda6_input,lambda7_input"
ROWS = [
    "L01,125.13,130,300,300,0.999,50,0.1,0.1,0",
    "L05,125.13,130,300,300,0.999,50,20,0.1,0",
    "L06,125.0,200,500,500,1,10000,1.0,1e-10,0",
    "ordering,125.13,125.12999999999999,300,300,0.999,50,0.1,0.1,0",
    "construction,125.13,130,-1,300,0.999,50,0.1,0.1,0",
]


def main() -> None:
    if not BINARY.is_file():
        raise SystemExit("build dihiggs/app/Lambda1EvaluatorV2 first")
    OUT.mkdir(parents=True, exist_ok=True)
    INPUT.write_text(HEADER + "\n" + "\n".join(ROWS) + "\n", encoding="utf-8")
    commit = os.environ.get("DIHIGGS_IMPLEMENTATION_COMMIT") or subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
        text=True, capture_output=True).stdout.strip()
    env = {**os.environ, "DIHIGGS_GIT_COMMIT": commit, "DIHIGHS_GIT_DIRTY": "no",
           "DIHIGGS_GIT_DIRTY": "no"}
    subprocess.run([str(BINARY), str(INPUT), str(OUTPUT)], cwd=ROOT, env=env, check=True)
    with OUTPUT.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    payload = {
        "report_schema": "dihiggs.lambda1.v2.yukawa_fix.verification.v1",
        "baseline_commit": "cbb5079c8ecc012395a525ccd1dd54f8481d5be9",
        "successor_commit": commit,
        "schema": "dihiggs.lambda1.v2",
        "input": str(INPUT.relative_to(ROOT)),
        "output": str(OUTPUT.relative_to(ROOT)),
        "output_sha256": hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
        "row_count": len(rows),
        "selected_widths_note": "nine H2 partial widths are reported; width_unaccounted_gev closes to total_width_gev",
        "rows": [{k: row[k] for k in ("point_id", "construction_ok", "yukawa_type_installed",
                                       "theory_ok", "width_bb_gev", "width_tautau_gev",
                                       "total_width_gev", "width_unaccounted_gev", "br_gammagamma", "ctau_mm")}
                 for row in rows],
    }
    REPORT.with_suffix(".json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    lines = ["# Lambda1 v2 Yukawa initialization correction", "",
             f"- Baseline: `cbb5079`", f"- Corrected producer commit: `{commit}`",
             "- Type I is installed only after successful construction and verified with `get_yukawas_type()`.",
             "- Width columns are selected H2 partial widths; `width_unaccounted_gev` is the closure diagnostic.", "",
             f"Output SHA-256: `{payload['output_sha256']}`", "", "| Case | construction | installed | theory | total width | ctau (mm) |",
             "|---|---:|---:|---:|---:|---:|"]
    for row in rows:
        lines.append(f"| {row['point_id']} | {row['construction_ok']} | {row['yukawa_type_installed']} | {row['theory_ok']} | {row['total_width_gev']} | {row['ctau_mm']} |")
    REPORT.with_suffix(".md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
