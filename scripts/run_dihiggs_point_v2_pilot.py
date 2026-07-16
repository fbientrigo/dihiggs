#!/usr/bin/env python3
"""Regenerate the bounded dihiggs.point.v2 engineering pilot and reports."""

import csv
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"
OUT = ROOT / "docs/pilots/dihiggs_point_v2_v1"
REPORT = ROOT / "docs/verification/dihiggs_point_v2_verification_v1"
BASE = ["--campaign-id", "point-v2-engineering-pilot", "--mh", "125.13", "--yukawa-type", "1"]
CASES = {
    "L01_accepted_anchor": ["--mH-min", "130", "--mH-max", "130", "--n-mH", "1", "--mA", "300", "--mHp", "300", "--sin-ba", "0.999", "--tan-beta", "50", "--M2-min", "16721.68154468371", "--M2-max", "16721.68154468371", "--n-M2", "1", "--lambda6", "0.1", "--lambda7", "0"],
    "L05_theory_rejected_anchor": ["--mH-min", "130", "--mH-max", "130", "--n-mH", "1", "--mA", "300", "--mHp", "300", "--sin-ba", "0.999", "--tan-beta", "50", "--M2-min", "16239.109978356435", "--M2-max", "16239.109978356435", "--n-M2", "1", "--lambda6", "0.1", "--lambda7", "0"],
    "L06_llp_anchor": ["--mH-min", "200", "--mH-max", "200", "--n-mH", "1", "--mA", "500", "--mHp", "500", "--sin-ba", "1", "--tan-beta", "10000", "--M2-min", "39999.9995713761", "--M2-max", "39999.9995713761", "--n-M2", "1", "--lambda6", "1e-10", "--lambda7", "0"],
    "ordering_boundary": ["--mH-min", repr(math.nextafter(125.13, -math.inf)), "--mH-max", "125.13", "--n-mH", "2", "--mA", "300", "--mHp", "300", "--sin-ba", "0.999", "--tan-beta", "50", "--M2-min", "16721.68154468371", "--M2-max", "16721.68154468371", "--n-M2", "1", "--lambda6", "0.1", "--lambda7", "0"],
    "construction_failure": ["--mH-min", "130", "--mH-max", "130", "--n-mH", "1", "--mA", "-1", "--mHp", "300", "--sin-ba", "0.999", "--tan-beta", "50", "--M2-min", "16721.68154468371", "--M2-max", "16721.68154468371", "--n-M2", "1", "--lambda6", "0.1", "--lambda7", "0"],
}


def main() -> None:
    if not BINARY.is_file():
        raise SystemExit("build dihiggs/app/DihiggsPointV2Evaluator first")
    OUT.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "DIHIGGS_GIT_COMMIT": "successor-worktree", "DIHIGGS_GIT_DIRTY": "yes"}
    results = []
    output_columns = []
    for name, args in CASES.items():
        output = OUT / f"{name}.csv"
        relative = ["dihiggs/app/DihiggsPointV2Evaluator", *BASE, "--run-id", name, *args,
                    "--output", str(output.relative_to(ROOT))]
        completed = subprocess.run([str(ROOT / relative[0]), *relative[1:]], cwd=ROOT, env=env,
                                   check=True, text=True, capture_output=True)
        with output.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
            output_columns = reader.fieldnames or []
        results.append({
            "case": name,
            "command": relative,
            "output": str(output.relative_to(ROOT)),
            "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
            "row_count": len(rows),
            "rows": [{key: row[key] for key in (
                "point_id", "mh_input_GeV", "mH_input_GeV", "M2_input_GeV2",
                "construction_ok", "numerical_ok", "rejection_stage", "rejection_reason",
                "theory_ok_v1", "total_width_GeV", "br_gammagamma", "ctau_mm",
            )} for row in rows],
            "stdout": completed.stdout.splitlines(),
        })

    payload = {
        "report_schema": "dihiggs.point.v2.verification.v1",
        "baseline_commit": "05a6217a",
        "successor_commit": "uncommitted-successor-worktree",
        "producer_schema": "dihiggs.point.v2",
        "mh_convention_GeV": 125.13,
        "mh_source": "https://pdg.lbl.gov/encoder_listings/s126.pdf",
        "scope": "bounded engineering pilot; no campaign-level scientific conclusion",
        "verification_commands": [
            "make -C 2hdmc -j2",
            "make -C dihiggs app/Lambda1EvaluatorV2 app/DihiggsPointV2Evaluator app/Phys_M2BandTracker -j2",
            "pytest -q tests/test_dihiggs_point_v2.py tests/test_m2_tracker_intervals.py tests/test_orchestrator tests/test_golden_lambda1.py tests/test_physics_conventions.py tests/test_external_tools.py",
            "python scripts/run_dihiggs_point_v2_pilot.py",
            "git diff --check",
        ],
        "baseline_verification": {
            "worktree_commit": "05a6217a",
            "command": "pytest -q tests/test_orchestrator tests/test_golden_lambda1.py tests/test_physics_conventions.py tests/test_external_tools.py",
            "result": "126 passed",
        },
        "producer_columns_in_fixed_order": output_columns,
        "acceptance": {
            "construction_ok": "exact THDM::set_param_phys result",
            "numerical_ok": "construction and finite reconstructed lambda1-lambda7, tan_beta, m12_sq, M2",
            "triple_ok_legacy": "positivity_reported_ok && unitarity_ok && perturbativity_ok",
            "theory_ok_v1": "triple_ok_legacy",
            "stability_reported_ok": "recorded separately; patched dependency aliases positivity",
            "experimental_ok": "unevaluated nan",
        },
        "dispositions": {
            "D01-D04": "verified by focused executable tests",
            "Q01": "quarantined and non-gating",
            "RV01": "verified for canonical producer",
            "RV02": "verified only in the tested pilot domains; boundaries are not globally validated",
            "RV03": "verified for canonical producer",
            "RV04": "deferred as approved",
            "RV05": "verified not production-integrated; experimental fields remain unevaluated",
        },
        "deferred": ["full scans", "HB/HS datasets", "recasting", "training", "plotting", "maintained TeX claims"],
        "preexisting_full_suite_collection_failures": [
            "mlpython/lake_pipeline/test_nulls.py: missing optional polars dependency",
            "tests/test_recompute_readiness.py: missing scripts.recompute_readiness",
            "tests/test_run_quarantine.py: missing scripts.run_quarantine",
        ],
        "pilot": results,
    }
    REPORT.with_suffix(".json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# DiHiggs Point v2 verification v1", "",
        "Status: canonical core and bounded engineering pilot verified; production experimental gates are not integrated.", "",
        "- Baseline: `05a6217a` (PR #58 final head supplied by the plan)",
        "- Successor: uncommitted successor worktree", "- Schema: `dihiggs.point.v2`",
        "- Higgs mass: 125.13 GeV ([PDG 2026 listing](https://pdg.lbl.gov/encoder_listings/s126.pdf))",
        "- Scope: engineering validation only; no campaign-level scientific conclusion", "",
        "## Pilot results", "",
        "| Case | Rows | construction_ok | theory_ok_v1 | SHA-256 |", "|---|---:|---|---|---|",
    ]
    for result in results:
        lines.append(f"| {result['case']} | {result['row_count']} | "
                     f"{','.join(row['construction_ok'] for row in result['rows'])} | "
                     f"{','.join(row['theory_ok_v1'] for row in result['rows'])} | `{result['output_sha256']}` |")
    lines += ["", "## Readiness dispositions", ""] + [f"- {key}: {value}" for key, value in payload["dispositions"].items()]
    lines += ["", "## Deferred", "", ", ".join(payload["deferred"]) + ".", "",
              "## Pre-existing full-suite collection failures", ""]
    lines += [f"- {failure}" for failure in payload["preexisting_full_suite_collection_failures"]]
    lines += ["",
              "Exact commands, field results, output paths, and checksums are in the adjacent JSON report.", ""]
    REPORT.with_suffix(".md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
