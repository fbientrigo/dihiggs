#!/usr/bin/env python3
"""Bounded 15-point H2 lifetime scan using the canonical dihiggs.lambda1.v2
evaluator (dihiggs/app/Lambda1EvaluatorV2). Grid and stopping policy are the
ones specified in benchmarks/FIRST_H2_RECAST_CANDIDATE.md ("Follow-up scan
specification"). Deterministic order, resumable, stops at the first clean
candidate (see classify() for the clean/provisional/rejected policy).

Usage: python3 benchmarks/run_first_h2_bounded_scan.py
"""
import csv
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "dihiggs/app/Lambda1EvaluatorV2"
CANONICAL_CSV = ROOT / "benchmarks/first_h2_bounded_scan.csv"
MANIFEST = ROOT / "benchmarks/first_h2_bounded_scan_manifest.json"

HEADER_IN = ("point_id,mh_gev,mH_gev,mA_gev,mHp_gev,sin_beta_minus_alpha,"
             "tan_beta,lambda1_target,lambda6_input,lambda7_input")

FIXED = dict(mh_gev=125.13, sin_beta_minus_alpha=1.0, lambda1_target=1.0,
             lambda6_input=1e-10, lambda7_input=0.0)

MASSES_GEV = [150.0, 200.0, 250.0]
# Priority order per FIRST_H2_RECAST_CANDIDATE.md: 3e5 and 1e5 first (bracket
# the cot^2(beta) lifetime-scaling estimate), then 1e6, then 3e4, then 1e4.
TAN_BETA_PRIORITY = [3.0e5, 1.0e5, 1.0e6, 3.0e4, 1.0e4]

CTAU_MIN_MM = 1.0
CTAU_MAX_MM = 1.0e4
BR_BB_MIN = 0.05
MAX_POINTS = 15


def build_grid():
    grid = []
    idx = 0
    for tan_beta in TAN_BETA_PRIORITY:
        for mH in MASSES_GEV:
            idx += 1
            point_id = f"H2scan_mH{int(mH)}_tb{int(tan_beta)}"
            grid.append({
                "priority_index": idx,
                "point_id": point_id,
                "mH_gev": mH,
                "mA_gev": mH + 300.0,
                "mHp_gev": mH + 300.0,
                "tan_beta": tan_beta,
            })
    assert len(grid) <= MAX_POINTS, "grid exceeds the bounded-scan point budget"
    return grid


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, check=True, text=True,
                           capture_output=True).stdout.strip()


def producer_state():
    commit = git("rev-parse", "HEAD")
    dirty_out = git("status", "--porcelain", "--", ".", ":!benchmarks")
    dirty = "yes" if dirty_out else "no"
    return commit, dirty


def classify(row):
    theory_ok = row["theory_ok"] == "1"
    width_ok = row["width_ok"] == "1"
    if not (theory_ok and width_ok):
        return "rejected", "theory_ok=0 or width_ok=0"
    ctau = float(row["ctau_mm"])
    br_bb = float(row["br_bb"])
    if not (CTAU_MIN_MM <= ctau <= CTAU_MAX_MM):
        return "rejected", f"ctau_mm={ctau!r} outside [{CTAU_MIN_MM},{CTAU_MAX_MM}] mm"
    if not (br_bb > BR_BB_MIN):
        return "rejected", f"br_bb={br_bb!r} <= {BR_BB_MIN}"
    if row["lambda1_residual_warning"] == "1":
        return "provisional", "lambda1_residual_warning=1"
    return "clean", "theory_ok=1, width_ok=1, ctau in window, br_bb>0.05, lambda1_residual_warning=0"


def load_existing():
    if not CANONICAL_CSV.is_file():
        return []
    with CANONICAL_CSV.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    if not BINARY.is_file():
        raise SystemExit("build dihiggs/app/Lambda1EvaluatorV2 first (see dihiggs/Makefile target)")

    commit, dirty = producer_state()
    if dirty == "yes":
        raise SystemExit("refusing to run: worktree has uncommitted changes outside benchmarks/")

    grid = build_grid()
    all_rows = load_existing()
    if all_rows:
        seen_states = {(r["evaluator_commit"], r["evaluator_dirty"]) for r in all_rows}
        if seen_states != {(commit, dirty)}:
            raise SystemExit(
                f"refusing to mix outputs: {CANONICAL_CSV} already contains rows from "
                f"producer state(s) {seen_states}, current state is ({commit!r}, {dirty!r}). "
                "Move or delete the existing file to rerun from scratch."
            )
    done_ids = {r["point_id"] for r in all_rows}
    fieldnames = list(all_rows[0].keys()) if all_rows else None

    clean_found = next((r for r in all_rows if r.get("classification") == "clean"), None)
    provisional_found = next((r for r in all_rows if r.get("classification") == "provisional"), None)

    env = {**os.environ, "DIHIGGS_GIT_COMMIT": commit, "DIHIGGS_GIT_DIRTY": dirty}
    scratch_in = ROOT / "benchmarks/.scan_point_input.csv.tmp"
    scratch_out = ROOT / "benchmarks/.scan_point_output.csv.tmp"

    evaluated_this_run = 0
    stop_reason = None

    if clean_found is not None:
        stop_reason = "clean_candidate_already_found_in_prior_run"
    else:
        for point in grid:
            if point["point_id"] in done_ids:
                continue
            row_csv = (f"{point['point_id']},{FIXED['mh_gev']},{point['mH_gev']},"
                       f"{point['mA_gev']},{point['mHp_gev']},{FIXED['sin_beta_minus_alpha']},"
                       f"{point['tan_beta']},{FIXED['lambda1_target']},{FIXED['lambda6_input']},"
                       f"{FIXED['lambda7_input']}")
            scratch_in.write_text(HEADER_IN + "\n" + row_csv + "\n", encoding="utf-8")
            subprocess.run([str(BINARY), str(scratch_in), str(scratch_out)],
                            cwd=ROOT, env=env, check=True, capture_output=True, text=True)
            with scratch_out.open(newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                if fieldnames is None:
                    fieldnames = reader.fieldnames + [
                        "scan_priority_index", "scan_mH_gev", "scan_tan_beta",
                        "classification", "classification_reason",
                    ]
                result_row = next(reader)

            classification, reason = classify(result_row)
            result_row["scan_priority_index"] = point["priority_index"]
            result_row["scan_mH_gev"] = point["mH_gev"]
            result_row["scan_tan_beta"] = point["tan_beta"]
            result_row["classification"] = classification
            result_row["classification_reason"] = reason
            all_rows.append(result_row)
            done_ids.add(point["point_id"])
            evaluated_this_run += 1

            # Rewrite the full canonical CSV after every point so an
            # interruption never loses more than the in-flight point.
            with CANONICAL_CSV.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for r in all_rows:
                    writer.writerow(r)

            if classification == "clean":
                clean_found = result_row
                stop_reason = "clean_candidate_found"
                break
            if classification == "provisional" and provisional_found is None:
                provisional_found = result_row
        else:
            stop_reason = "exhausted_grid"

    scratch_in.unlink(missing_ok=True)
    scratch_out.unlink(missing_ok=True)

    seen_states = {(r["evaluator_commit"], r["evaluator_dirty"]) for r in all_rows}
    if seen_states and seen_states != {(commit, dirty)}:
        raise SystemExit(f"post-run consistency check failed: mixed producer states {seen_states}")

    manifest = {
        "schema": "dihiggs.h2_bounded_scan.manifest.v1",
        "producer_commit": commit,
        "producer_dirty": dirty,
        "evaluator_binary": str(BINARY.relative_to(ROOT)),
        "evaluator_schema": "dihiggs.lambda1.v2",
        "fixed_parameters": FIXED,
        "mass_relation": "mA_gev = mHp_gev = mH_gev + 300",
        "grid_declared_in_priority_order": grid,
        "max_points": MAX_POINTS,
        "points_evaluated_total": len(all_rows),
        "points_evaluated_this_run": evaluated_this_run,
        "stop_reason": stop_reason,
        "classification_criteria": {
            "clean": "theory_ok=1 and width_ok=1 and 1<=ctau_mm<=1e4 and br_bb>0.05 and lambda1_residual_warning=0",
            "provisional": "theory_ok=1 and width_ok=1 and 1<=ctau_mm<=1e4 and br_bb>0.05 and lambda1_residual_warning=1",
            "rejected": "fails theory_ok, width_ok, the ctau window, or the br_bb threshold",
        },
        "clean_candidate_point_id": clean_found["point_id"] if clean_found else None,
        "provisional_candidate_point_id": provisional_found["point_id"] if provisional_found else None,
        "classification_counts": {
            "clean": sum(1 for r in all_rows if r.get("classification") == "clean"),
            "provisional": sum(1 for r in all_rows if r.get("classification") == "provisional"),
            "rejected": sum(1 for r in all_rows if r.get("classification") == "rejected"),
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "stop_reason": stop_reason,
        "points_evaluated_this_run": evaluated_this_run,
        "points_evaluated_total": len(all_rows),
        "clean_candidate_point_id": manifest["clean_candidate_point_id"],
        "provisional_candidate_point_id": manifest["provisional_candidate_point_id"],
    }, indent=2))


if __name__ == "__main__":
    main()
