from __future__ import annotations

from autoresearch.harness.decision_layer import next_axis_selector, reachability_analysis


def test_reachability_basic_projection() -> None:
    out = reachability_analysis(
        path_trend_summary={"slope_dex_per_step": 0.145, "tan_beta_center": 124416, "lambda6_center": 0.0019683, "tan_beta_factor": 1.02, "lambda6_factor": 1.03},
        sensitivity_summary={"new_best_ctau_m": 5.390644846600711e-4},
        envelope={"tan_beta": [30000, 500000], "lambda6": [0.0003, 0.0060]},
    )
    assert out["best_ctau_m"] is not None
    assert out["improvement_factor_to_target"] > 1000
    assert out["estimated_steps_to_target"] is not None


def test_next_axis_selector_prefers_mphi_edge() -> None:
    tops = [{"m_phi": 200.0, "mA": 500.0} for _ in range(8)] + [{"m_phi": 230.0, "mA": 480.0} for _ in range(2)]
    rec = next_axis_selector(
        per_mA_summary=[],
        global_top_by_ctau=tops,
        path_trend_summary={"slope_dex_per_step": 0.145},
        triple_ok_rate=0.95,
        mphi_bounds=(200.0, 270.0),
        mA_bounds=(300.0, 500.0),
        target_ctau_m=1.0,
        best_ctau_m=5e-4,
    )
    assert rec["axis"] == "mphi"
    assert rec["action"] == "refine"


def test_next_axis_selector_backoff_on_low_ok_rate() -> None:
    rec = next_axis_selector(
        per_mA_summary=[],
        global_top_by_ctau=[],
        path_trend_summary={},
        triple_ok_rate=0.5,
        mphi_bounds=(200.0, 270.0),
        mA_bounds=(300.0, 500.0),
        target_ctau_m=1.0,
        best_ctau_m=5e-4,
    )
    assert rec["action"] == "backoff"
