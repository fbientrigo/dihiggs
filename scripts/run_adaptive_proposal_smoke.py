#!/usr/bin/env python3
"""Task 4 smoke fixture for the adaptive proposal-batch handoff (issue #72).

Three proposals, run through dihiggs.app.orchestrator.proposal_batch TWICE
with identical arguments, to demonstrate:

  Point A: exact frozen H2scan_mH150_tb300000 benchmark
           (benchmarks/H2scan_mH150_tb300000_production_coupling.json).
  Point B: a nearby constructed proposal (mH shifted 150 -> 160 GeV) that
           exercises arbitrary-point handling. NOT claimed as a new
           benchmark; its scientific status is whatever the evaluator
           returns, unmodified below.
  Point C: the repository's existing L05 theory-rejected anchor
           (tests/test_dihiggs_point_v2.py), which constructs but fails a
           theory predicate.

Writes docs/pilots/adaptive_proposal_v1/smoke_report.json: a plain-text,
reviewable report of both runs (inputs, per-attempt canonical status,
checksums, provenance, and the two-run replay comparison). The raw
attempts.csv/proposals.csv this script generates are NOT committed (all
*.csv in this repo are git-LFS filtered); their SHA-256 hashes and this
script's invocation are the durable, reviewable record instead.
"""
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dihiggs.app.orchestrator.io_utils import utc_now_iso  # noqa: E402
from dihiggs.app.orchestrator.proposal_batch import run_proposal_batch  # noqa: E402

EXECUTABLE = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"
SMOKE_DIR = ROOT / "docs/pilots/adaptive_proposal_v1"
REPORT_PATH = SMOKE_DIR / "smoke_report.json"

PROPOSALS = [
    {
        "proposal_id": "A_frozen_H2scan_mH150_tb300000",
        "mH_GeV": "150.0", "mA_GeV": "450.0", "mHp_GeV": "450.0",
        "M2_GeV2": "22499.999999500335", "tan_beta": "300000.0",
        "lambda6": "1e-10", "lambda7": "0.0",
        "sin_beta_minus_alpha": "1.0", "yukawa_type": "1",
    },
    {
        "proposal_id": "B_arbitrary_mH160_perturbation",
        "mH_GeV": "160.0", "mA_GeV": "450.0", "mHp_GeV": "450.0",
        "M2_GeV2": "22499.999999500335", "tan_beta": "300000.0",
        "lambda6": "1e-10", "lambda7": "0.0",
        "sin_beta_minus_alpha": "1.0", "yukawa_type": "1",
    },
    {
        "proposal_id": "C_theory_rejected_L05_anchor",
        "mH_GeV": "130.0", "mA_GeV": "300.0", "mHp_GeV": "300.0",
        "M2_GeV2": "16239.109978356435", "tan_beta": "50.0",
        "lambda6": "0.1", "lambda7": "0.0",
        "sin_beta_minus_alpha": "0.999", "yukawa_type": "1",
    },
]

FROZEN_ANCHOR = {
    "point_id": "point_c7afb83ab8127e47",
    "g_hH2H2_GeV": 63.5914252007596588,
    "total_width_GeV": 4.56118529862185007e-14,
    "ctau_mm": 4.32622152973311191,
    "br_bb": 0.756737485808578692,
}


def write_proposals_csv(path: Path) -> None:
    fieldnames = list(PROPOSALS[0].keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(PROPOSALS)


def run_once(run_name: str, workdir: Path) -> dict:
    proposals_csv = workdir / "input_proposals.csv"
    write_proposals_csv(proposals_csv)
    manifest = run_proposal_batch(
        executable=EXECUTABLE, proposals_csv=proposals_csv, outdir=workdir / "out",
        campaign_id="adaptive_proposal_smoke_v1", run_name=run_name, repo_root=ROOT,
    )
    attempts_csv = Path(manifest["attempts_csv"])
    with attempts_csv.open(newline="", encoding="utf-8") as handle:
        attempts = list(csv.DictReader(handle))
    return {"manifest": manifest, "attempts": attempts}


def summarize_attempt(row: dict) -> dict:
    return {
        "proposal_id": row["proposal_id"],
        "attempt_status": row["attempt_status"],
        "point_id": row.get("point_id") or None,
        "construction_ok": row.get("construction_ok") or None,
        "theory_ok_v1": row.get("theory_ok_v1") or None,
        "rejection_stage": row.get("rejection_stage") or None,
        "rejection_reason": row.get("rejection_reason") or None,
        "g_hH2H2_GeV": row.get("g_hH2H2_GeV") or None,
        "total_width_GeV": row.get("total_width_GeV") or None,
        "ctau_mm": row.get("ctau_mm") or None,
        "br_bb": row.get("br_bb") or None,
    }


def main() -> int:
    import tempfile

    ok = True
    runs = []
    with tempfile.TemporaryDirectory(prefix="adaptive_proposal_smoke_") as tmp:
        tmp_path = Path(tmp)
        for run_name in ("run1", "run2"):
            result = run_once(run_name, tmp_path / run_name)
            runs.append(result)

        run1_attempts, run2_attempts = runs[0]["attempts"], runs[1]["attempts"]

        replay = {
            "output_sha256_identical": runs[0]["manifest"]["output_sha256"] == runs[1]["manifest"]["output_sha256"],
            "point_ids_identical": [a["point_id"] for a in run1_attempts] == [b["point_id"] for b in run2_attempts],
            "canonical_fields_byte_identical": run1_attempts == run2_attempts,
        }
        ok = ok and all(replay.values())

        point_a = next(r for r in run1_attempts if r["proposal_id"] == "A_frozen_H2scan_mH150_tb300000")
        anchor_check = {
            "point_id_matches": point_a["point_id"] == FROZEN_ANCHOR["point_id"],
            "g_hH2H2_GeV_matches": abs(float(point_a["g_hH2H2_GeV"]) - FROZEN_ANCHOR["g_hH2H2_GeV"]) < 1e-10,
            "total_width_GeV_matches": abs(float(point_a["total_width_GeV"]) - FROZEN_ANCHOR["total_width_GeV"]) < 1e-24,
            "ctau_mm_matches": abs(float(point_a["ctau_mm"]) - FROZEN_ANCHOR["ctau_mm"]) < 1e-10,
            "br_bb_matches": abs(float(point_a["br_bb"]) - FROZEN_ANCHOR["br_bb"]) < 1e-12,
        }
        ok = ok and all(anchor_check.values())

        point_c = next(r for r in run1_attempts if r["proposal_id"] == "C_theory_rejected_L05_anchor")
        rejection_check = {
            "attempt_status_evaluated": point_c["attempt_status"] == "EVALUATED",
            "construction_ok": point_c["construction_ok"] == "1",
            "theory_rejected": point_c["theory_ok_v1"] == "0.00000000000000000e+00",
            "rejection_stage_present": bool(point_c["rejection_stage"]) and point_c["rejection_stage"] != "accepted",
        }
        ok = ok and all(rejection_check.values())

        row_counts_ok = len(run1_attempts) == len(PROPOSALS) == len(run2_attempts)
        ok = ok and row_counts_ok

        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
                                 text=True, capture_output=True).stdout.strip()

        report = {
            "report_schema": "dihiggs.adaptive_proposal_batch.smoke_report.v1",
            "generated_utc": utc_now_iso(),
            "repository_commit": commit,
            "regeneration_command": "python scripts/run_adaptive_proposal_smoke.py",
            "input_proposals": PROPOSALS,
            "frozen_anchor_reference": {
                "source": "benchmarks/H2scan_mH150_tb300000_production_coupling.json",
                **FROZEN_ANCHOR,
            },
            "runs": [
                {
                    "run_name": ("run1", "run2")[i],
                    "manifest": runs[i]["manifest"],
                    "attempts": [summarize_attempt(r) for r in runs[i]["attempts"]],
                }
                for i in range(2)
            ],
            "checks": {
                "row_counts_ok": row_counts_ok,
                "point_a_matches_frozen_anchor": anchor_check,
                "point_c_theory_rejected_and_preserved": rejection_check,
                "deterministic_replay": replay,
            },
            "all_checks_passed": ok,
        }

    SMOKE_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report["checks"], indent=2))
    print(f"all_checks_passed={ok}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
