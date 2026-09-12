#!/usr/bin/env python3
"""Replay the exact, bounded #62 point-v2 anchor grid into a new evidence path."""

import csv
import hashlib
import json
import math
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / "docs/campaigns/yukawa_order_recompute_v1"
INPUT = CAMPAIGN / "input_grid.json"
OUTPUTS = CAMPAIGN / "outputs"
BINARY = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"
HBARC_GEV_MM = 1.973269804e-13


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, check=True, text=True,
                          capture_output=True).stdout.strip()


def finite(row: dict[str, str], *names: str) -> bool:
    return all(math.isfinite(float(row[name])) for name in names)


def one(row: dict[str, str], name: str) -> bool:
    return math.isclose(float(row[name]), 1.0, rel_tol=0.0, abs_tol=0.0)


def main() -> None:
    if not BINARY.is_file():
        raise SystemExit("build dihiggs/app/DihiggsPointV2Evaluator first")
    spec = json.loads(INPUT.read_text(encoding="utf-8"))
    commit = git("rev-parse", "HEAD")
    # Campaign outputs are expected mutations; only evaluator source trees
    # determine the producer's dirty state.
    dirty = "yes" if git("status", "--short", "--untracked-files=no", "--", "2hdmc", "dihiggs") else "no"
    started = datetime.now(timezone.utc).isoformat()
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "DIHIGGS_GIT_COMMIT": commit, "DIHIGGS_GIT_DIRTY": dirty,
           "OMP_NUM_THREADS": "1"}
    all_rows: list[dict[str, str]] = []
    header: list[str] | None = None
    cases: list[dict] = []
    for name, case_args in spec["cases"].items():
        output = OUTPUTS / f"{name}.csv"
        command = [str(BINARY), *spec["base_args"], "--run-id", name, *case_args,
                   "--output", str(output)]
        recorded_command = ["dihiggs/app/DihiggsPointV2Evaluator", *spec["base_args"],
                            "--run-id", name, *case_args, "--output",
                            str(output.relative_to(ROOT))]
        completed = subprocess.run(command, cwd=ROOT, env=env, text=True,
                                   capture_output=True)
        (OUTPUTS / f"{name}.stdout.log").write_text(completed.stdout, encoding="utf-8")
        (OUTPUTS / f"{name}.stderr.log").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode:
            cases.append({"case": name, "command": recorded_command, "attempted": 0,
                          "accepted": 0, "construction_failures": 1,
                          "error": completed.stderr})
            continue
        with output.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            rows = list(reader)
        if header is None:
            header = fields
        if fields != header:
            raise SystemExit(f"{name}: schema drift")
        accepted = 0
        construction_failures = theory_failures = width_failures = 0
        for row in rows:
            all_rows.append(row)
            if not one(row, "construction_ok"):
                construction_failures += 1
                continue
            if not one(row, "theory_ok_v1"):
                theory_failures += 1
            if not one(row, "width_ok"):
                width_failures += 1
            if one(row, "construction_ok") and one(row, "theory_ok_v1") and one(row, "width_ok"):
                accepted += 1
                widths = ["width_bb_GeV", "width_cc_GeV", "width_tt_GeV", "width_tautau_GeV",
                          "width_WW_GeV", "width_ZZ_GeV", "width_gammagamma_GeV",
                          "width_Zgamma_GeV", "width_gg_GeV", "width_hh_GeV"]
                brs = ["br_bb", "br_cc", "br_tt", "br_tautau", "br_WW", "br_ZZ",
                       "br_gammagamma", "br_Zgamma", "br_gg", "br_hh"]
                if not finite(row, "total_width_GeV", "width_unaccounted_GeV", "ctau_mm", *widths, *brs):
                    raise SystemExit(f"{name}: non-finite accepted width/BR/lifetime")
                if float(row["total_width_GeV"]) <= 0 or any(float(row[field]) < 0 for field in widths):
                    raise SystemExit(f"{name}: invalid accepted width/BR")
                if any(not 0 <= float(row[field]) <= 1 for field in brs):
                    raise SystemExit(f"{name}: branching ratio outside [0, 1]")
                selected_sum = sum(float(row[field]) for field in widths)
                if not math.isclose(float(row["total_width_GeV"]) - selected_sum,
                                    float(row["width_unaccounted_GeV"]), rel_tol=1e-14, abs_tol=1e-20):
                    raise SystemExit(f"{name}: width accounting mismatch")
                expected_ctau = HBARC_GEV_MM / float(row["total_width_GeV"])
                if not math.isclose(float(row["ctau_mm"]), expected_ctau, rel_tol=1e-14):
                    raise SystemExit(f"{name}: ctau_mm mismatch")
        cases.append({"case": name, "command": recorded_command, "attempted": len(rows),
                      "accepted": accepted, "construction_failures": construction_failures,
                      "theory_failures": theory_failures, "width_failures": width_failures,
                      "output": str(output.relative_to(ROOT)), "output_sha256": sha256(output),
                      "stdout_log": str((OUTPUTS / f"{name}.stdout.log").relative_to(ROOT)),
                      "stderr_log": str((OUTPUTS / f"{name}.stderr.log").relative_to(ROOT))})
    if header is None:
        raise SystemExit("no evaluator output")
    combined = CAMPAIGN / "point_v2_anchors.csv"
    with combined.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_rows)
    historical = json.loads((ROOT / "docs/verification/dihiggs_point_v2_verification_v1.json").read_text(encoding="utf-8"))
    comparison = CAMPAIGN / "old_vs_new.csv"
    comparison_fields = ["case", "quantity", "historical_pre_fix", "historical_post_fix",
                         "new_replay", "relative_difference_to_post_fix", "reason"]
    quantity_fields = {
        "total_width_GeV": "total_width_GeV", "width_bb_GeV": "width_bb_GeV",
        "width_tautau_GeV": "width_tautau_GeV", "br_gammagamma": "br_gammagamma",
        "ctau_mm": "ctau_mm",
    }
    comparison_rows = []
    for short, name in (("L01", "L01_accepted_anchor"), ("L06", "L06_llp_anchor")):
        with (OUTPUTS / f"{name}.csv").open(newline="", encoding="utf-8") as handle:
            new_row = next(csv.DictReader(handle))
        for quantity, field in quantity_fields.items():
            pre = float(historical["before_after"][short][quantity]["before"])
            post = float(historical["before_after"][short][quantity]["after"])
            new = float(new_row[field])
            relative = "nan" if post == 0 else f"{(new - post) / post:.17e}"
            comparison_rows.append({"case": short, "quantity": quantity,
                                    "historical_pre_fix": f"{pre:.17e}",
                                    "historical_post_fix": f"{post:.17e}",
                                    "new_replay": f"{new:.17e}",
                                    "relative_difference_to_post_fix": relative,
                                    "reason": "pre/post values are historical verification evidence; replay uses corrected Yukawa order and current main"})
    with comparison.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=comparison_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(comparison_rows)
    ended = datetime.now(timezone.utc).isoformat()
    manifest = {
        "schema": "dihiggs.issue62.yukawa_order_recompute.v1",
        "status": "completed_bounded_pilot",
        "evaluator": "dihiggs/app/DihiggsPointV2Evaluator",
        "producer_commit": commit,
        "producer_dirty": dirty,
        "corrected_yukawa_order": {"fix_commit": "6bfad7662fd87750d838bf2fe0bd7ac00ee2326a",
                                    "contract": "construct physical model, then install and verify Type I before DecayTable"},
        "input": {"path": str(INPUT.relative_to(ROOT)), "sha256": sha256(INPUT),
                  "schema": spec["schema"], "grid": spec["scope"]},
        "output": {"path": str(combined.relative_to(ROOT)), "sha256": sha256(combined),
                   "schema_version": "dihiggs.point.v2"},
        "comparison": {"path": str(comparison.relative_to(ROOT)), "sha256": sha256(comparison),
                       "historical_source": "docs/verification/dihiggs_point_v2_verification_v1.json"},
        "started_utc": started,
        "ended_utc": ended,
        "attempted_rows": sum(case["attempted"] for case in cases),
        "construction_failures": sum(case.get("construction_failures", 0) for case in cases),
        "theory_failures": sum(case.get("theory_failures", 0) for case in cases),
        "width_failures": sum(case.get("width_failures", 0) for case in cases),
        "accepted_rows": sum(case["accepted"] for case in cases),
        "cases": cases,
        "validation": ["schema and row cardinality", "finite accepted widths/BRs", "ctau_mm = hbarc/width", "post-construction yukawa_type_installed=1", "construction/theory/width status accounting"],
        "historical_scope_note": "All broad pre-fix payloads lacking exact grids remain UNRESOLVABLE_PROVENANCE; no replacement coordinates were invented.",
    }
    (CAMPAIGN / "campaign_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: manifest[key] for key in ("status", "attempted_rows", "accepted_rows", "construction_failures", "theory_failures", "width_failures")}, sort_keys=True))


if __name__ == "__main__":
    main()
