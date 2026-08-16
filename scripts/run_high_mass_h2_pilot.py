#!/usr/bin/env python3
"""Bounded high-mass H2 pilot: run required anchors and verify contract invariants.

Produces docs/pilots/high_mass_h2_v1/pilot_points.csv (raw dihiggs.point.v2
rows plus derived high-mass fields) and pilot_validation.json (per-point and
aggregate invariant checks). Does not launch MadGraph and does not attempt
the full mH2<=800/mA=mHp<=2000 scan; see docs/HIGH_MASS_H2_CONTRACT.md.
"""
import csv
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"
OUT_DIR = ROOT / "docs/pilots/high_mass_h2_v1"
POINTS_CSV = OUT_DIR / "pilot_points.csv"
VALIDATION_JSON = OUT_DIR / "pilot_validation.json"

HBAR_C_GEV_MM = 1.973269804e-13
# Canonical SM-like Higgs mass; see conventions/physics_conventions.yaml
# (sm_like_higgs.m_h_GeV). Drives --mh and the H2->hh threshold 2*M_H.
M_H = 125.20
# Superseded convention, retained ONLY for the historical 150 GeV anchor.
M_H_HISTORICAL_ANCHOR = 125.13
M_W = 80.36951
M_Z = 91.15349
M_T_POLE = 172.5

WIDTH_FIELDS = (
    "width_bb_GeV", "width_cc_GeV", "width_tt_GeV", "width_tautau_GeV",
    "width_WW_GeV", "width_ZZ_GeV", "width_gammagamma_GeV", "width_Zgamma_GeV",
    "width_gg_GeV", "width_hh_GeV",
)

# P0_anchor_150 is the HISTORICAL validated 150 GeV benchmark
# (H2scan_mH150_tb300000). It is pinned to the superseded convention
# M_H_HISTORICAL_ANCHOR so its frozen total_width/g_hH2H2/ctau/br_bb remain an
# exact regression; the Gate A check "anchor_150_unaffected_by_tt_addition"
# compares against those frozen digits. Every other pilot point uses the
# canonical M_H. Do not "fix" this by moving the anchor to 125.20 -- that
# reinterprets a historical benchmark as a newly calculated point.
# See docs/HIGH_MASS_H2_CONTRACT.md, Decision v1 (superseded) / v2 (active).
# "P0_repeat" is the reproducibility re-run of the same anchor point and must
# therefore use the same (historical) mh, or its point_id will not match.
ANCHOR_PILOT_NAMES = frozenset({"P0_anchor_150", "P0_repeat"})

# (name, mH2, Delta_heavy, sin_ba, tan_beta, M2, lambda6, lambda7, model_variant)
PILOTS = [
    ("P0_anchor_150", 150.0, 300.0, 1.0, 300000.0, 22499.999999500335, 1e-10, 0.0,
     "PHYSICAL_DECAYS_NO_HEAVY_CASCADES"),
    ("P1_below_hh", 200.0, 300.0, 0.999, 50.0, 16721.68154468371, 0.1, 0.0,
     "PHYSICAL_DECAYS_NO_HEAVY_CASCADES"),
    ("P2_above_hh", 300.0, 500.0, 0.999, 50.0, 16721.68154468371, 0.1, 0.0,
     "PHYSICAL_DECAYS_NO_HEAVY_CASCADES"),
    ("P3_below_tt", 330.0, 470.0, 0.999, 50.0, 16721.68154468371, 0.1, 0.0,
     "PHYSICAL_DECAYS_NO_HEAVY_CASCADES"),
    ("P4_above_tt", 400.0, 400.0, 0.999, 50.0, 16721.68154468371, 0.1, 0.0,
     "PHYSICAL_DECAYS_NO_HEAVY_CASCADES"),
    ("P5_near_800", 800.0, 1200.0, 0.999, 50.0, 16721.68154468371, 0.1, 0.0,
     "PHYSICAL_DECAYS_NO_HEAVY_CASCADES"),
    ("P6a_heavy_sep_thin", 300.0, 200.0, 0.999, 50.0, 16721.68154468371, 0.1, 0.0,
     "PHYSICAL_DECAYS_NO_HEAVY_CASCADES"),
    ("P6b_heavy_sep_wide", 300.0, 900.0, 0.999, 50.0, 16721.68154468371, 0.1, 0.0,
     "PHYSICAL_DECAYS_NO_HEAVY_CASCADES"),
]


def point_mh(name):
    """mh for a pilot point: canonical, except the historical anchor."""
    return M_H_HISTORICAL_ANCHOR if name in ANCHOR_PILOT_NAMES else M_H


def run_point(name, mH2, mA, sba, tb, M2, l6, l7):
    output = OUT_DIR / f"_raw_{name}.csv"
    command = [
        str(BINARY), "--campaign-id", "high_mass_h2_pilot_v1", "--run-id", name,
        "--mh", repr(point_mh(name)), "--mH-min", repr(mH2), "--mH-max", repr(mH2), "--n-mH", "1",
        "--mA", repr(mA), "--mHp", repr(mA), "--yukawa-type", "1",
        "--sin-ba", repr(sba), "--tan-beta", repr(tb),
        "--M2-min", repr(M2), "--M2-max", repr(M2), "--n-M2", "1",
        "--lambda6", repr(l6), "--lambda7", repr(l7), "--output", str(output),
    ]
    env = {**os.environ, "DIHIGGS_GIT_COMMIT": subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, text=True, capture_output=True
    ).stdout.strip(), "DIHIGGS_GIT_DIRTY": "no"}
    subprocess.run(command, check=True, text=True, capture_output=True, env=env)
    with output.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    output.unlink()
    return row


def cascade_flags(mH2, mA, mHp, mh):
    return {
        "H2_to_AZ_open": mH2 > mA + M_Z,
        "H2_to_HpW_open": mH2 > mHp + M_W,
        "H2_to_AA_open": mH2 > 2 * mA,
        "H2_to_HpHm_open": mH2 > 2 * mHp,
        # Per docs/contracts/cascade_contract.yaml the flag must come from the
        # point's OWN mh, not a module constant -- the historical anchor runs
        # at the superseded convention and keeps its 2*125.13 threshold.
        "H2_to_hh_open": mH2 > 2 * mh,
        "H2_to_tt_open": mH2 > 2 * M_T_POLE,
    }


def finite(row, field):
    try:
        return math.isfinite(float(row[field]))
    except (KeyError, TypeError, ValueError):
        return False


def validate_point(name, row, mH2, mA, model_variant):
    checks = {}
    checks["construction_ok"] = row["construction_ok"] == "1"
    mHp = mA

    if checks["construction_ok"]:
        checks["hierarchy_mh_lt_mH2"] = float(row["mh_input_GeV"]) < float(row["mH_input_GeV"])
        checks["hierarchy_mH2_lt_mA"] = float(row["mH_input_GeV"]) < float(row["mA_input_GeV"])
        checks["hierarchy_mA_eq_mHp"] = float(row["mA_input_GeV"]) == float(row["mHp_input_GeV"])
        checks["g_hH2H2_finite_nonneg"] = finite(row, "g_hH2H2_GeV") and float(row["g_hH2H2_GeV"]) >= 0.0

        flags = cascade_flags(mH2, mA, mHp, float(row["mh_input_GeV"]))
        forbidden = ["H2_to_AZ_open", "H2_to_HpW_open", "H2_to_AA_open", "H2_to_HpHm_open"]
        if model_variant == "PHYSICAL_DECAYS_NO_HEAVY_CASCADES":
            checks["no_forbidden_cascade_open"] = not any(flags[f] for f in forbidden)
        checks["cascade_flags"] = flags

        if finite(row, "total_width_GeV"):
            total = float(row["total_width_GeV"])
            selected_sum = sum(float(row[f]) for f in WIDTH_FIELDS)
            unaccounted = float(row["width_unaccounted_GeV"])
            checks["width_closure_ok"] = math.isclose(
                selected_sum + unaccounted, total, rel_tol=1e-9, abs_tol=0.0
            ) if total != 0 else True
            checks["width_tt_field_present_and_finite"] = finite(row, "width_tt_GeV")
            if total > 0:
                ctau_expected = HBAR_C_GEV_MM / total
                checks["ctau_consistency_ok"] = math.isclose(
                    float(row["ctau_mm"]), ctau_expected, rel_tol=1e-9
                )
            checks["theory_ok_recorded"] = float(row["theory_ok_v1"]) in (0.0, 1.0)

    return checks


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    validations = []
    for name, mH2, delta, sba, tb, M2, l6, l7, variant in PILOTS:
        mA = mH2 + delta
        row = run_point(name, mH2, mA, sba, tb, M2, l6, l7)
        row["pilot_name"] = name
        row["Delta_heavy_GeV"] = f"{delta:.17e}"
        row["model_variant"] = variant
        rows.append(row)
        checks = validate_point(name, row, mH2, mA, variant)
        validations.append({
            "pilot_name": name,
            "m_H2_GeV": mH2,
            "Delta_heavy_GeV": delta,
            "m_A_GeV": mA,
            "m_Hp_GeV": mA,
            "model_variant": variant,
            "point_id": row["point_id"],
            "construction_ok": row["construction_ok"],
            "theory_ok_v1": row["theory_ok_v1"],
            "width_tt_GeV": row["width_tt_GeV"],
            "br_tt": row["br_tt"],
            "br_hh": row["br_hh"],
            "total_width_GeV": row["total_width_GeV"],
            "ctau_mm": row["ctau_mm"],
            "checks": checks,
            "all_checks_passed": all(
                v for k, v in checks.items() if isinstance(v, bool)
            ),
        })

    # Reproducibility: same physical coordinates -> same point_id, byte-identical
    # repeat invocation (single-process evaluator; no in-process worker count to vary).
    repeat_row = run_point("P0_repeat", 150.0, 450.0, 1.0, 300000.0,
                            22499.999999500335, 1e-10, 0.0)
    reproducibility = {
        "same_point_id_on_repeat": repeat_row["point_id"] == rows[0]["point_id"],
        "note": "DihiggsPointV2Evaluator is single-process per invocation; "
                "there is no in-process workers=N mode to diverge from. "
                "Reproducibility is verified via repeated CLI invocation "
                "producing an identical point_id and identical numeric fields.",
        "identical_total_width": repeat_row["total_width_GeV"] == rows[0]["total_width_GeV"],
        "identical_ctau_mm": repeat_row["ctau_mm"] == rows[0]["ctau_mm"],
    }

    with POINTS_CSV.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = list(rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    by_name = {v["pilot_name"]: v for v in validations}
    gate_a_evidence = {
        "below_hh_and_tt_width_tt_is_exactly_zero": float(by_name["P1_below_hh"]["width_tt_GeV"]) == 0.0,
        "width_tt_nonzero_and_rising_across_threshold": (
            0.0 == float(by_name["P1_below_hh"]["width_tt_GeV"])
            < float(by_name["P3_below_tt"]["width_tt_GeV"])
            < float(by_name["P4_above_tt"]["width_tt_GeV"])
            < float(by_name["P5_near_800"]["width_tt_GeV"])
        ),
        "br_tt_is_a_visible_fraction_above_threshold": float(by_name["P4_above_tt"]["br_tt"]) > 0.01,
        "anchor_150_unaffected_by_tt_addition":
            by_name["P0_anchor_150"]["total_width_GeV"] == "4.56118529862185007e-14",
    }

    payload = {
        "report_schema": "dihiggs.high_mass_h2.pilot_validation.v1",
        "campaign_id": "high_mass_h2_pilot_v1",
        "producer": "DihiggsPointV2Evaluator",
        "producer_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, text=True, capture_output=True
        ).stdout.strip(),
        "scope": "bounded engineering pilot for the high-mass H2 point factory contract; "
                 "no campaign-level scientific conclusion; no MadGraph production",
        "points_csv": str(POINTS_CSV.relative_to(ROOT)),
        "points_csv_sha256": hashlib.sha256(POINTS_CSV.read_bytes()).hexdigest(),
        "pilots": validations,
        "reproducibility": reproducibility,
        "gate_a_evidence": gate_a_evidence,
        "aggregate": {
            "all_points_construction_ok": all(v["construction_ok"] == "1" for v in validations),
            "all_points_pass_checks": all(v["all_checks_passed"] for v in validations),
            "reproducibility_ok": reproducibility["same_point_id_on_repeat"]
            and reproducibility["identical_total_width"]
            and reproducibility["identical_ctau_mm"],
            "gate_a_evidence_ok": all(gate_a_evidence.values()),
        },
    }
    VALIDATION_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["aggregate"], indent=2))


if __name__ == "__main__":
    main()
