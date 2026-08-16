"""Regression tests for the high-mass H2 point factory contract.

Covers: Gate A (explicit width_tt_GeV/br_tt), the cascade-state contract
(docs/contracts/cascade_contract.yaml), and the bounded pilot's own
invariants (docs/contracts/high_mass_point_schema.yaml). These tests
exercise the real DihiggsPointV2Evaluator binary; they do not launch
MadGraph and do not run the full mH2<=800/mA=mHp<=2000 scan.
"""
import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"
PILOT_SCRIPT = ROOT / "scripts/run_high_mass_h2_pilot.py"

# Canonical convention; mirrors scripts/run_high_mass_h2_pilot.py and
# conventions/physics_conventions.yaml (sm_like_higgs.m_h_GeV).
M_H = 125.20
M_W = 80.36951
M_Z = 91.15349
M_T_POLE = 172.5


@pytest.fixture(scope="module", autouse=True)
def require_binary():
    if not BINARY.is_file():
        pytest.skip("build dihiggs/app/DihiggsPointV2Evaluator first")


def run_point(tmp_path, name, *, mH2, mA, sba=0.999, tb=50.0,
              M2=16721.68154468371, l6=0.1, l7=0.0):
    output = tmp_path / f"{name}.csv"
    command = [
        str(BINARY), "--campaign-id", "contract-test", "--run-id", name,
        "--mh", repr(M_H), "--mH-min", repr(mH2), "--mH-max", repr(mH2), "--n-mH", "1",
        "--mA", repr(mA), "--mHp", repr(mA), "--yukawa-type", "1",
        "--sin-ba", repr(sba), "--tan-beta", repr(tb),
        "--M2-min", repr(M2), "--M2-max", repr(M2), "--n-M2", "1",
        "--lambda6", repr(l6), "--lambda7", repr(l7), "--output", str(output),
    ]
    env = {**os.environ, "DIHIGGS_GIT_COMMIT": "test-commit", "DIHIGGS_GIT_DIRTY": "no"}
    subprocess.run(command, check=True, text=True, capture_output=True, env=env)
    with output.open(newline="") as handle:
        return next(csv.DictReader(handle))


def cascade_flags(mH2, mA, mHp):
    return {
        "H2_to_AZ_open": mH2 > mA + M_Z,
        "H2_to_HpW_open": mH2 > mHp + M_W,
        "H2_to_AA_open": mH2 > 2 * mA,
        "H2_to_HpHm_open": mH2 > 2 * mHp,
        "H2_to_hh_open": mH2 > 2 * M_H,
        "H2_to_tt_open": mH2 > 2 * M_T_POLE,
    }


FORBIDDEN = ("H2_to_AZ_open", "H2_to_HpW_open", "H2_to_AA_open", "H2_to_HpHm_open")


def test_cascade_flags_forbidden_group_closed_for_valid_hierarchy():
    # mHp = mA > mH2 by construction; no forbidden channel can be open.
    for mH2, mA in ((150.0, 450.0), (300.0, 800.0), (800.0, 2000.0), (799.9, 800.1)):
        flags = cascade_flags(mH2, mA, mA)
        assert not any(flags[name] for name in FORBIDDEN), (mH2, mA, flags)


def test_cascade_flags_fail_closed_for_broken_hierarchy():
    # A deliberately broken hierarchy (mH2 far above mA=mHp) must trip at
    # least one forbidden flag so a campaign runner has something to fail
    # closed on.
    flags = cascade_flags(mH2=900.0, mA=400.0, mHp=400.0)
    assert flags["H2_to_AA_open"] is True
    assert flags["H2_to_HpHm_open"] is True
    assert any(flags[name] for name in FORBIDDEN)


def test_cascade_flags_physical_group_not_treated_as_forbidden():
    below = cascade_flags(mH2=200.0, mA=500.0, mHp=500.0)
    above = cascade_flags(mH2=800.0, mA=2000.0, mHp=2000.0)
    assert below["H2_to_hh_open"] is False and below["H2_to_tt_open"] is False
    assert above["H2_to_hh_open"] is True and above["H2_to_tt_open"] is True


def test_width_tt_field_present_for_every_selected_channel(tmp_path):
    row = run_point(tmp_path, "closure", mH2=400.0, mA=800.0)
    assert row["construction_ok"] == "1"
    assert math.isfinite(float(row["width_tt_GeV"]))
    assert math.isfinite(float(row["br_tt"]))


@pytest.mark.parametrize("mH2,mA", [(150.0, 450.0), (250.0, 800.0), (400.0, 800.0), (800.0, 2000.0)])
def test_width_closure_and_ctau_invariants_hold_across_mass_range(tmp_path, mH2, mA):
    row = run_point(tmp_path, f"closure-{mH2}", mH2=mH2, mA=mA)
    assert row["construction_ok"] == "1"
    total = float(row["total_width_GeV"])
    assert math.isfinite(total) and total > 0.0
    selected = sum(float(row[f]) for f in (
        "width_bb_GeV", "width_cc_GeV", "width_tt_GeV", "width_tautau_GeV",
        "width_WW_GeV", "width_ZZ_GeV", "width_gammagamma_GeV", "width_Zgamma_GeV",
        "width_gg_GeV", "width_hh_GeV",
    ))
    unaccounted = float(row["width_unaccounted_GeV"])
    assert selected + unaccounted == pytest.approx(total, rel=1e-9)
    expected_ctau = 1.973269804e-13 / total
    assert float(row["ctau_mm"]) == pytest.approx(expected_ctau, rel=1e-9)
    assert math.isfinite(float(row["g_hH2H2_GeV"])) and float(row["g_hH2H2_GeV"]) >= 0.0


def test_hierarchy_holds_for_all_pilot_anchors(tmp_path):
    for mH2, mA in ((150.0, 450.0), (300.0, 800.0), (800.0, 2000.0)):
        row = run_point(tmp_path, f"hier-{mH2}", mH2=mH2, mA=mA)
        assert float(row["mh_input_GeV"]) < float(row["mH_input_GeV"]) < float(row["mA_input_GeV"])
        assert float(row["mA_input_GeV"]) == float(row["mHp_input_GeV"])


def test_pilot_script_regenerates_and_all_checks_pass(tmp_path):
    completed = subprocess.run(
        [sys.executable, str(PILOT_SCRIPT)], cwd=ROOT, check=True, text=True, capture_output=True,
    )
    aggregate = json.loads(completed.stdout)
    assert aggregate["all_points_construction_ok"] is True
    assert aggregate["all_points_pass_checks"] is True
    assert aggregate["reproducibility_ok"] is True
    assert aggregate["gate_a_evidence_ok"] is True

    validation_path = ROOT / "docs/pilots/high_mass_h2_v1/pilot_validation.json"
    payload = json.loads(validation_path.read_text())
    assert len(payload["pilots"]) == 8
    assert payload["aggregate"] == aggregate
