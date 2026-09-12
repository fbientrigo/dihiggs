#!/usr/bin/env python3
"""Fresh, seed-informed Yukawa-order recomputation on the current point.v2 schema.

Historical points define seed locations only.  No historical value is used as a
numeric acceptance constraint; this run is evaluated entirely by the current
physical/theory/width contracts.
"""

import csv
import hashlib
import json
import math
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / "docs/campaigns/yukawa_order_recompute_v2"
INPUT = CAMPAIGN / "input_grid.json"
OUTPUTS = CAMPAIGN / "outputs"
BINARY = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"
HBARC_GEV_MM = 1.973269804e-13
WIDTHS = ("width_bb_GeV", "width_cc_GeV", "width_tt_GeV", "width_tautau_GeV",
          "width_WW_GeV", "width_ZZ_GeV", "width_gammagamma_GeV",
          "width_Zgamma_GeV", "width_gg_GeV", "width_hh_GeV")
BRS = tuple(field.replace("width_", "br_").replace("_GeV", "") for field in WIDTHS)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, check=True, text=True,
                          capture_output=True).stdout.strip()


def is_one(row: dict[str, str], name: str) -> bool:
    return math.isclose(float(row[name]), 1.0, rel_tol=0.0, abs_tol=0.0)


def validate(row: dict[str, str]) -> str:
    if not is_one(row, "construction_ok"):
        return "construction"
    if not math.isclose(float(row["yukawa_type_installed"]), 1.0, rel_tol=0.0, abs_tol=0.0):
        raise SystemExit(f"{row['point_id']}: corrected Type-I Yukawa was not installed")
    if not is_one(row, "theory_ok_v1"):
        return "theory"
    if not is_one(row, "width_ok"):
        return "width"
    numbers = ("total_width_GeV", "width_unaccounted_GeV", "ctau_mm", *WIDTHS, *BRS)
    if any(not math.isfinite(float(row[name])) for name in numbers):
        raise SystemExit(f"{row['point_id']}: non-finite accepted width/BR/lifetime")
    if float(row["total_width_GeV"]) <= 0 or any(float(row[name]) < 0 for name in WIDTHS):
        raise SystemExit(f"{row['point_id']}: invalid accepted width")
    if any(not 0 <= float(row[name]) <= 1 for name in BRS):
        raise SystemExit(f"{row['point_id']}: branching ratio outside [0, 1]")
    selected = sum(float(row[name]) for name in WIDTHS)
    if not math.isclose(float(row["total_width_GeV"]) - selected,
                        float(row["width_unaccounted_GeV"]), rel_tol=1e-14, abs_tol=1e-20):
        raise SystemExit(f"{row['point_id']}: width accounting mismatch")
    if not math.isclose(float(row["ctau_mm"]), HBARC_GEV_MM / float(row["total_width_GeV"]),
                        rel_tol=1e-14):
        raise SystemExit(f"{row['point_id']}: ctau_mm mismatch")
    return "accepted"


def main() -> None:
    if not BINARY.is_file():
        raise SystemExit("build dihiggs/app/DihiggsPointV2Evaluator first")
    spec = json.loads(INPUT.read_text(encoding="utf-8"))
    commit = git("rev-parse", "HEAD")
    dirty_paths = ("2hdmc", "dihiggs", str(INPUT.relative_to(ROOT)),
                   "scripts/run_yukawa_order_recompute_v2.py")
    dirty = "yes" if git("status", "--short", "--untracked-files=no", "--", *dirty_paths) else "no"
    started = datetime.now(timezone.utc).isoformat()
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "DIHIGGS_GIT_COMMIT": commit, "DIHIGGS_GIT_DIRTY": dirty,
           "OMP_NUM_THREADS": "1"}
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    cases = []
    for family in spec["families"]:
        name = family["name"]
        output = OUTPUTS / f"{name}.csv"
        args = ["--mH-min", str(family["mH_min"]), "--mH-max", str(family["mH_max"]),
                "--n-mH", str(family["n_mH"]), "--mA", str(family["mA"]),
                "--mHp", str(family["mHp"]), "--sin-ba", str(family["sin_ba"]),
                "--tan-beta", str(family["tan_beta"]), "--M2-min", str(family["M2_min"]),
                "--M2-max", str(family["M2_max"]), "--n-M2", str(family["n_M2"]),
                "--lambda6", str(family["lambda6"]), "--lambda7", str(family["lambda7"])]
        command = [str(BINARY), *spec["base_args"], "--run-id", name, *args,
                   "--output", str(output)]
        recorded = ["dihiggs/app/DihiggsPointV2Evaluator", *spec["base_args"],
                    "--run-id", name, *args, "--output", str(output.relative_to(ROOT))]
        result = subprocess.run(command, cwd=ROOT, env=env, text=True, capture_output=True)
        stdout = OUTPUTS / f"{name}.stdout.log"
        stderr = OUTPUTS / f"{name}.stderr.log"
        stdout.write_text(result.stdout, encoding="utf-8")
        stderr.write_text(result.stderr, encoding="utf-8")
        case = {"family": name, "seed": family["seed"], "command": recorded,
                "stdout_log": str(stdout.relative_to(ROOT)), "stderr_log": str(stderr.relative_to(ROOT))}
        if result.returncode:
            case.update(attempted=0, accepted=0, construction_failures=1, error=result.stderr)
            cases.append(case)
            continue
        with output.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            family_rows = list(reader)
        expected = int(family["n_mH"]) * int(family["n_M2"])
        if len(family_rows) != expected:
            raise SystemExit(f"{name}: expected {expected} rows, got {len(family_rows)}")
        if header is None:
            header = fields
        if fields != header:
            raise SystemExit(f"{name}: current-schema drift within campaign")
        counts = {"accepted": 0, "construction": 0, "theory": 0, "width": 0}
        for row in family_rows:
            outcome = validate(row)
            counts[outcome] += 1
            rows.append(row)
        case.update(attempted=len(family_rows), accepted=counts["accepted"],
                    construction_failures=counts["construction"], theory_failures=counts["theory"],
                    width_failures=counts["width"], output=str(output.relative_to(ROOT)),
                    output_sha256=sha256(output))
        cases.append(case)
    if header is None:
        raise SystemExit("no evaluator output")
    point_ids = [row["point_id"] for row in rows]
    if len(point_ids) != len(set(point_ids)):
        raise SystemExit("duplicate point_id in fresh campaign output")
    combined = CAMPAIGN / "fresh_point_v2.csv"
    with combined.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    ended = datetime.now(timezone.utc).isoformat()
    process_failures = sum("error" in case for case in cases)
    manifest = {
        "schema": "dihiggs.issue62.yukawa_order_recompute.v2",
        "status": "completed_fresh_seed_informed" if not process_failures else "incomplete_process_failure",
        "evaluator": "dihiggs/app/DihiggsPointV2Evaluator",
        "producer_commit": commit, "producer_dirty": dirty,
        "corrected_yukawa_order": {"fix_commit": "6bfad7662fd87750d838bf2fe0bd7ac00ee2326a",
                                    "contract": "construct physical model, then install and verify Type I before DecayTable"},
        "input": {"path": str(INPUT.relative_to(ROOT)), "sha256": sha256(INPUT),
                  "schema": spec["schema"], "grid": spec["scope"]},
        "output": {"path": str(combined.relative_to(ROOT)), "sha256": sha256(combined),
                   "schema_version": "dihiggs.point.v2"},
        "started_utc": started, "ended_utc": ended,
        "attempted_rows": sum(c.get("attempted", 0) for c in cases),
        "construction_failures": sum(c.get("construction_failures", 0) for c in cases),
        "theory_failures": sum(c.get("theory_failures", 0) for c in cases),
        "width_failures": sum(c.get("width_failures", 0) for c in cases),
        "process_failures": process_failures,
        "accepted_rows": sum(c.get("accepted", 0) for c in cases), "families": cases,
        "validation": ["current schema and row cardinality", "unique point IDs", "finite accepted widths/BRs",
                       "ctau_mm = hbarc/width", "width accounting", "post-construction Type-I installation"],
        "seed_policy": "Known historical points are seed locations only; no historical values, checksums, or column layouts constrain acceptance.",
        "historical_scope_note": "Broad pre-fix payloads without exact grids remain UNRESOLVABLE_PROVENANCE and were not regenerated from filenames.",
    }
    (CAMPAIGN / "campaign_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: manifest[k] for k in ("status", "attempted_rows", "accepted_rows", "construction_failures", "theory_failures", "width_failures")}, sort_keys=True))
    if process_failures:
        raise SystemExit("campaign incomplete: evaluator process failure(s), see manifest")


if __name__ == "__main__":
    main()
