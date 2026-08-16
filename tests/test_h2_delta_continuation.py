import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "dihiggs" / "app" / "DihiggsPointV2Evaluator"
sys.path.insert(0, str(ROOT / "scripts"))
import h2_delta_continuation as m  # noqa: E402


@pytest.fixture(scope="module", autouse=True)
def require_binary():
    if not BINARY.is_file():
        pytest.skip("build dihiggs/app/DihiggsPointV2Evaluator first")


def test_anchor_reproduces_frozen_record(tmp_path):
    ok, row, diffs = m.verify_anchor(tmp_path, open(tmp_path / "log.txt", "w"))
    assert ok
    assert row["construction_ok"] == "1"
    assert row["theory_ok_v1"] == m.VALID_TOKEN
    assert row["width_ok"] == m.VALID_TOKEN
    for field, reldiff in diffs.items():
        assert reldiff < 1e-9, f"{field} diverged from frozen record by {reldiff}"


def test_same_point_evaluated_twice_is_deterministic(tmp_path):
    row_a = m.evaluate_point(m.ANCHOR_MH2, m.ANCHOR_M2, 450.0, 450.0, tmp_path,
                              run_id="a")
    row_b = m.evaluate_point(m.ANCHOR_MH2, m.ANCHOR_M2, 450.0, 450.0, tmp_path,
                              run_id="b")
    numeric_fields = ["theory_ok_v1", "g_hH2H2_GeV", "total_width_GeV",
                       "ctau_mm", "br_bb", "M2_reconstructed_GeV2"]
    for field in numeric_fields:
        assert row_a[field] == row_b[field], field


def test_hierarchy_invariant_mA_greater_than_mH2(tmp_path):
    for mH2 in [150.0, 300.0, 600.0]:
        delta = 300.0
        mA = mH2 + delta
        mHp = mA
        assert mA > mH2
        assert mHp > mH2
        row = m.evaluate_point(mH2, m.ANCHOR_M2 * (mH2 / m.ANCHOR_MH2) ** 2,
                                mA, mHp, tmp_path, run_id=f"hier{mH2:g}")
        assert float(row["mA_input_GeV"]) > float(row["mH_input_GeV"])
        assert float(row["mHp_input_GeV"]) > float(row["mH_input_GeV"])


def test_candidate_selection_picks_nearest_to_predicted():
    rows = [
        {"M2_input_GeV2": "100.0"},
        {"M2_input_GeV2": "150.0"},
        {"M2_input_GeV2": "205.0"},
        {"M2_input_GeV2": "500.0"},
    ]
    best = m.pick_nearest(rows, target=200.0)
    assert best["M2_input_GeV2"] == "205.0"


def test_failed_theory_point_is_recorded_failed_not_retried(monkeypatch, tmp_path):
    """If the evaluator reports theory_ok_v1=0 everywhere in the search
    window, search_local_M2 must return None (a physics failure) without
    ever changing tan_beta or lambda6 -- both stay at the caller-supplied
    values on every evaluate_grid call it makes."""
    calls = []

    def fake_evaluate_grid(mH2, M2_lo, M2_hi, n_M2, mA, mHp, campaign_id,
                            run_id, out_path, tan_beta=m.TAN_BETA,
                            lambda6=m.LAMBDA6, **kw):
        calls.append((tan_beta, lambda6))
        return [{
            "M2_input_GeV2": f"{M2_lo:.6f}", "theory_ok_v1": "0.00000000000000000e+00",
            "construction_ok": "1", "rejection_stage": "perturbativity",
        }]

    monkeypatch.setattr(m, "evaluate_grid", fake_evaluate_grid)
    logf = open(tmp_path / "log.txt", "w")
    best, n_probes, n_rounds, moved, reason = m.search_local_M2(
        200.0, 500.0, 500.0, [{"mH2": 150.0, "M2": 22500.0}], tmp_path, logf,
        "test")
    assert best is None
    assert reason == "perturbativity"
    assert moved == "none"
    assert all(tb == m.TAN_BETA for tb, l6 in calls)
    assert all(l6 == m.LAMBDA6 for tb, l6 in calls)
