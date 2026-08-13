#!/usr/bin/env python3
"""Assemble results/high_mass_valid_points/manifest.json from the search
manifest fragment (written by high_mass_physical_point_search.py) plus
classification counts read back from high_mass_valid_points.csv, and the
HB/HS enrichment attempt outcome (recorded by hand in this script, since the
HB/HS pipeline in dihiggs_boundary was found to be schema-incompatible with
this campaign's frozen mh=125.13 convention -- see HIGH_MASS_VALID_POINT_REPORT.md
section on HB/HS enrichment for the full explanation).
"""
import csv
import json
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results/high_mass_valid_points"


def main():
    fragment = json.loads((RESULTS_DIR / "search_manifest_fragment.json").read_text())

    classification_counts = Counter()
    region_counts = Counter()
    with (RESULTS_DIR / "high_mass_valid_points.csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames:
            for row in reader:
                classification_counts[row["classification"]] += 1
                region_counts[row["mass_region"]] += 1

    with (RESULTS_DIR / "high_mass_attempted_points.csv").open(newline="", encoding="utf-8") as handle:
        attempted_total = sum(1 for _ in csv.DictReader(handle))

    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
                           text=True, capture_output=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, check=True,
                                 text=True, capture_output=True).stdout.strip())

    manifest = {
        "manifest_schema": "dihiggs.high_mass_h2.worker_a_manifest.v1",
        "mission": "Find theory-valid (and where possible HB/HS-valid) high-mass H2 2HDM points",
        "worktree": str(ROOT),
        "worktree_git_head": head,
        "worktree_git_dirty": dirty,
        "evaluator_binary": fragment["evaluator_binary"],
        "evaluator_producer_commit": fragment["producer_commit"],
        "evaluator_producer_dirty": fragment["producer_dirty"],
        "search_script": fragment["search_script"],
        "search_config_hash_sha256": fragment["config_hash_sha256"],
        "search_regions": fragment["regions"],
        "search_grid": fragment["search_grid"],
        "total_cli_invocations": fragment["total_cli_invocations"],
        "total_attempted_points": attempted_total,
        "total_accepted_theory_valid_points": fragment["total_accepted_points"],
        "search_wall_time_seconds": fragment["wall_time_seconds"],
        "classification_counts": dict(classification_counts),
        "accepted_points_by_mass_region": dict(region_counts),
        "hbhs_enrichment": {
            "attempted": True,
            "outcome": "HBHS_BLOCKED for all points",
            "reason": (
                "dihiggs_boundary's only tested, working HB/HS-linked producer "
                "(src/evaluate_point.cpp, wrapped by python/dhb/enrich.py + "
                "runner.HbhsRunner) hard-codes mh=125.09 GeV (constexpr kMh in "
                "evaluate_point.cpp; confirmed in docs/evaluate_point_contract.md "
                "and tests/test_golden_evaluate_point.py::test_g07_...). "
                "HIGH_MASS_H2_CONTRACT.md section 3 explicitly freezes mh=125.13 "
                "for this campaign and explicitly forbids mixing 125.09 results "
                "into any high-mass campaign artifact. evaluate_point also hard-codes "
                "sin(beta-alpha)=1.0 exactly (kSinBa) and only accepts a real M "
                "(so M2=M*M>=0), which does not cover every point in this "
                "campaign's search grid (sin_ba in {1.0,0.999,0.995}, M2 allowed "
                "negative). Bridging this would require either modifying "
                "dihiggs_boundary (explicitly out of scope, read-only access) or "
                "writing a new C++ program duplicating evaluate_point.cpp's "
                "2HDMC effective-coupling extraction (get_coupling_hdd/huu/hll "
                "for h1/h2/h3 normalized against a separately constructed SM-like "
                "reference model, plus charged-Higgs kappa_t/kappa_b) with "
                "mh=125.13 instead -- a nontrivial reimplementation, not a small "
                "conversion step, and out of this task's time-box. Per the "
                "mission's explicit allowance, this is reported as an unrelated "
                "operational/schema blocker rather than blocking the mission; "
                "every theory-valid point is classified HBHS_BLOCKED, not dropped."
            ),
            "files_read_before_concluding_blocked": [
                "dihiggs_boundary/python/dhb/runner.py",
                "dihiggs_boundary/python/dhb/enrich.py",
                "dihiggs_boundary/python/dhb/adapter.py",
                "dihiggs_boundary/python/dhb/schema.py",
                "dihiggs_boundary/src/evaluate_point.cpp",
                "dihiggs_boundary/docs/evaluate_point_contract.md",
                "dihiggs_boundary/tests/test_golden_evaluate_point.py",
                "dihiggs_boundary/README.md",
            ],
        },
        "outputs": {
            "high_mass_valid_points_csv": "high_mass_valid_points.csv",
            "high_mass_attempted_points_csv": "high_mass_attempted_points.csv",
            "rejection_summary_csv": "rejection_summary.csv",
            "mass_validity_map_png": "mass_validity_map.png",
            "mass_validity_map_pdf": "mass_validity_map.pdf",
            "report_md": "HIGH_MASS_VALID_POINT_REPORT.md",
            "checksums": "checksums.sha256",
        },
    }
    out = RESULTS_DIR / "manifest.json"
    out.write_text(json.dumps(manifest, indent=2, sort_keys=False) + "\n")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
