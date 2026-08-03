#!/usr/bin/env python3
"""Run the bounded 15-point H2 scan with reproducible provenance."""
import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "dihiggs/app/Lambda1EvaluatorV2"
CSV_OUT = ROOT / "benchmarks/first_h2_bounded_scan.csv"
MANIFEST = ROOT / "benchmarks/first_h2_bounded_scan_manifest.json"
TMP_IN = ROOT / "benchmarks/.scan_point_input.csv.tmp"
TMP_OUT = ROOT / "benchmarks/.scan_point_output.csv.tmp"
GENERATED = (CSV_OUT, MANIFEST, TMP_IN, TMP_OUT)
PROVENANCE_FILES = (Path(__file__).resolve(), ROOT / "dihiggs/src/Lambda1EvaluatorV2.cpp", ROOT / "dihiggs/src/ReplaySafeOutput.cpp")

HEADER = "point_id,mh_gev,mH_gev,mA_gev,mHp_gev,sin_beta_minus_alpha,tan_beta,lambda1_target,lambda6_input,lambda7_input"
FIXED = {"mh_gev": 125.13, "sin_beta_minus_alpha": 1.0, "lambda1_target": 1.0, "lambda6_input": 1e-10, "lambda7_input": 0.0}
MASSES = [150.0, 200.0, 250.0]
TAN_BETAS = [3e5, 1e5, 1e6, 3e4, 1e4]


def git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, check=True, text=True, capture_output=True).stdout.strip()


def rel(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def producer_state() -> tuple[str, str]:
    commit = git("rev-parse", "HEAD")
    args = ["status", "--porcelain", "--", ".", *(f":(exclude){rel(path)}" for path in GENERATED)]
    dirty_output = git(*args)
    if dirty_output:
        raise SystemExit("refusing to run: uncommitted experiment-defining changes\n" + dirty_output)
    return commit, "no"


def grid() -> list[dict[str, object]]:
    rows = []
    for tan_beta in TAN_BETAS:
        for mass in MASSES:
            rows.append({"priority_index": len(rows) + 1, "point_id": f"H2scan_mH{int(mass)}_tb{int(tan_beta)}", "mH_gev": mass, "mA_gev": mass + 300.0, "mHp_gev": mass + 300.0, "tan_beta": tan_beta})
    assert len(rows) == 15
    return rows


def classify(row: dict[str, str]) -> tuple[str, str]:
    if row["theory_ok"] != "1" or row["width_ok"] != "1":
        return "rejected", "theory_ok=0 or width_ok=0"
    if not 1.0 <= float(row["ctau_mm"]) <= 1e4:
        return "rejected", "ctau outside [1,1e4] mm"
    if float(row["br_bb"]) <= 0.05:
        return "rejected", "br_bb<=0.05"
    if row["lambda1_residual_warning"] == "1":
        return "provisional", "lambda1_residual_warning=1"
    return "clean", "all scan gates pass"


def load_existing() -> list[dict[str, str]]:
    if not CSV_OUT.exists():
        return []
    with CSV_OUT.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(rows: list[dict[str, str]], fields: list[str]) -> None:
    with CSV_OUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if not BINARY.is_file():
        raise SystemExit("build dihiggs/app/Lambda1EvaluatorV2 first")
    commit, dirty = producer_state()
    hashes = {rel(path): sha256(path) for path in PROVENANCE_FILES}
    scan_grid = grid()
    rows = load_existing()
    if rows and {(row["evaluator_commit"], row["evaluator_dirty"]) for row in rows} != {(commit, dirty)}:
        raise SystemExit("existing CSV has different evaluator provenance; archive CSV/manifest and rerun from scratch")

    fields = list(rows[0].keys()) if rows else None
    done = {row["point_id"] for row in rows}
    clean = next((row for row in rows if row.get("classification") == "clean"), None)
    provisional = next((row for row in rows if row.get("classification") == "provisional"), None)
    evaluated = 0
    stop_reason = "clean_candidate_already_found_in_prior_run" if clean else None
    env = {**os.environ, "DIHIGGS_GIT_COMMIT": commit, "DIHIGGS_GIT_DIRTY": dirty}

    try:
        if not clean:
            for point in scan_grid:
                if point["point_id"] in done:
                    continue
                values = [point["point_id"], FIXED["mh_gev"], point["mH_gev"], point["mA_gev"], point["mHp_gev"], FIXED["sin_beta_minus_alpha"], point["tan_beta"], FIXED["lambda1_target"], FIXED["lambda6_input"], FIXED["lambda7_input"]]
                TMP_IN.write_text(HEADER + "\n" + ",".join(map(str, values)) + "\n", encoding="utf-8")
                subprocess.run([str(BINARY), str(TMP_IN), str(TMP_OUT)], cwd=ROOT, env=env, check=True, capture_output=True, text=True)
                with TMP_OUT.open(newline="", encoding="utf-8") as stream:
                    reader = csv.DictReader(stream)
                    if fields is None:
                        fields = list(reader.fieldnames or []) + ["scan_priority_index", "scan_mH_gev", "scan_tan_beta", "classification", "classification_reason"]
                    result = next(reader)
                result["scan_priority_index"], result["scan_mH_gev"], result["scan_tan_beta"] = point["priority_index"], point["mH_gev"], point["tan_beta"]
                result["classification"], result["classification_reason"] = classify(result)
                rows.append(result)
                done.add(point["point_id"])
                evaluated += 1
                write_rows(rows, fields)
                if result["classification"] == "clean":
                    clean, stop_reason = result, "clean_candidate_found"
                    break
                if result["classification"] == "provisional" and provisional is None:
                    provisional = result
            else:
                stop_reason = "exhausted_grid"
    finally:
        TMP_IN.unlink(missing_ok=True)
        TMP_OUT.unlink(missing_ok=True)

    manifest = {
        "schema": "dihiggs.h2_bounded_scan.manifest.v2",
        "provenance_status": "REPRODUCIBLE_FROM_CLEAN_COMMIT",
        "producer_commit": commit,
        "producer_dirty": dirty,
        "producer_files_sha256": hashes,
        "generated_paths_excluded_from_dirty_check": [rel(path) for path in GENERATED],
        "evaluator_binary": rel(BINARY),
        "evaluator_schema": "dihiggs.lambda1.v2",
        "fixed_parameters": FIXED,
        "mass_relation": "mA_gev = mHp_gev = mH_gev + 300",
        "grid_declared_in_priority_order": scan_grid,
        "points_evaluated_total": len(rows),
        "points_evaluated_this_run": evaluated,
        "stop_reason": stop_reason,
        "clean_candidate_point_id": clean["point_id"] if clean else None,
        "provisional_candidate_point_id": provisional["point_id"] if provisional else None,
        "classification_counts": {name: sum(row.get("classification") == name for row in rows) for name in ("clean", "provisional", "rejected")}
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stop_reason": stop_reason, "points": len(rows), "producer_commit": commit, "producer_files_sha256": hashes}, indent=2))


if __name__ == "__main__":
    main()
