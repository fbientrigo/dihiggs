from __future__ import annotations

from pathlib import Path

from scripts.gemini_contract_templates import materialize_template
from autoresearch.harness.bounded_adaptive_search import (
    attempt_isolated_summary,
    build_orchestrator_command_for_subrun,
    build_orchestrator_commands_for_subcampaign,
    build_subcampaign,
    calculate_raw_points,
    calculate_raw_points_multiplicative,
    canonicalize_physical_row,
    detect_ctau_metric_source,
    derive_step_accounting_v3,
    enforce_raw_point_rules,
    next_thread_attempts,
    proposal_inside_envelope,
    run_bounded_search,
    select_best_anchor_row,
    should_retry_due_to_gate,
    validate_budget,
    validate_envelope,
)


def _contract(tmp_path: Path) -> dict[str, object]:
    return {
        "search_envelope": {
            "sin_ba": [1.0, 1.0],
            "mA": [300, 300],
            "lambda7": [0.0, 0.0],
            "lambda6": [0.0025, 0.0030],
            "tan_beta": [50000, 60000],
            "mphi": [285, 295],
            "lambda1": [4.0, 4.8],
        },
        "budget": {
            "max_iterations": 2,
            "min_raw_points_per_subrun": 1000,
            "max_raw_points_per_subrun": 2500,
            "max_total_points": 5000,
            "max_failed_subcampaigns": 2,
            "preferred_min_n_lam1": 100,
            "production_target_n_lam1": 500,
            "min_n_mphi_when_varied": 10,
            "max_per_subrun_minutes": 1,
            "min_remaining_minutes_to_start": 0,
        },
        "runtime": {
            "exec_path": "/home/fabi/dihiggs/dihiggs/app/PhysScanWithFixings",
            "outdir": str(tmp_path / "out"),
            "lake_name": "dihiggs_lake",
            "campaign": "unit_campaign",
        },
        "runtime_thread_policy": {
            "initial_threads": 8,
            "min_threads": 1,
            "backoff_factor": 2,
            "max_retries_per_subcampaign": 4,
            "retry_on_gsl_abort": True,
        },
    }


def test_raw_point_calculation_includes_tanbeta_count() -> None:
    assert calculate_raw_points([50000, 60000], 10, 100) == 2000


def test_min_raw_points_enforced_nlam1_first() -> None:
    b = validate_budget(_contract(Path("/tmp"))["budget"])  # type: ignore[arg-type]
    n_mphi, n_lam1, pts, reason = enforce_raw_point_rules(tanbeta_values=[50000, 60000], n_mphi=5, n_lam1=10, mphi_varies=True, budget=b)
    assert reason is None
    assert n_lam1 >= 100
    assert pts >= 1000
    assert n_mphi >= 10


def test_n_lam1_minimum_preferred_rule() -> None:
    b = validate_budget(_contract(Path("/tmp"))["budget"])  # type: ignore[arg-type]
    n_mphi, n_lam1, _, _ = enforce_raw_point_rules(tanbeta_values=[1, 2], n_mphi=10, n_lam1=50, mphi_varies=True, budget=b)
    assert n_lam1 >= 100
    assert n_mphi >= 10


def test_subcampaign_inside_envelope() -> None:
    c = _contract(Path("/tmp"))
    env = validate_envelope(c["search_envelope"])  # type: ignore[arg-type]
    b = validate_budget(c["budget"])  # type: ignore[arg-type]
    sc = build_subcampaign(1, "box_partition", env, {"tan_beta": 55000, "mphi": 290, "lambda1": 4.4, "lambda6": 0.0027}, b)
    assert sc.inside_envelope
    assert proposal_inside_envelope(sc.variable_ranges, env)


def test_orchestrator_command_generation() -> None:
    c = _contract(Path("/tmp"))
    env = validate_envelope(c["search_envelope"])  # type: ignore[arg-type]
    b = validate_budget(c["budget"])  # type: ignore[arg-type]
    sc = build_subcampaign(1, "box_partition", env, {"tan_beta": 55000, "mphi": 290, "lambda1": 4.4, "lambda6": 0.0027}, b)
    cmd = build_orchestrator_command_for_subrun(sc, c["runtime"], 4, "rn1")  # type: ignore[arg-type]
    assert "dihiggs/app/orchestrate_scans.py" in cmd
    assert "--threads" in cmd and cmd[cmd.index("--threads") + 1] == "4"
    assert "--force" not in cmd


def test_backoff_retries_on_gsl_failure_report() -> None:
    stop = {"gates": {"failures": [{"code": "gate_failed_due_to_gsl_signature", "explanation": "sigabrt"}]}}
    assert should_retry_due_to_gate(stop, {"retry_on_gsl_abort": True}) is True


def test_backoff_stops_after_min_threads() -> None:
    plan = next_thread_attempts({"initial_threads": 8, "min_threads": 1, "backoff_factor": 2, "max_retries_per_subcampaign": 4})
    assert plan == [8, 4, 2, 1]


def test_plan_only_does_not_execute(tmp_path: Path) -> None:
    c = _contract(tmp_path)
    state = run_bounded_search(c, execute=False, plan_only=True, output_dir=tmp_path / "artifacts", max_iterations_override=1)
    assert len(state["iterations"]) == 1
    assert (tmp_path / "artifacts" / "attempts.jsonl").exists()


def test_max_total_points_enforced_across_iterations(tmp_path: Path) -> None:
    c = _contract(tmp_path)
    c["budget"]["max_total_points"] = 2000  # type: ignore[index]
    state = run_bounded_search(c, execute=False, plan_only=True, output_dir=tmp_path / "artifacts")
    assert len(state["iterations"]) <= 1


def test_accepted_attempt_updates_state_with_mocked_safe(monkeypatch, tmp_path: Path) -> None:
    from autoresearch.harness import bounded_adaptive_search as mod

    def fake_run_safe_campaign(*args, **kwargs):
        return {
            "campaign_state": {"summary": {"n_health_failed": 0}, "runs": [], "slices": []},
            "rerun_manifest": {"summary": {"n_runs_selected": 0, "n_slices_selected": 0}},
            "stop_report": {"gates": {"passed": True, "failures": []}, "production_validation": True},
        }

    monkeypatch.setattr(mod, "run_safe_campaign", fake_run_safe_campaign)

    c = _contract(tmp_path)
    state = run_bounded_search(c, execute=False, plan_only=False, output_dir=tmp_path / "artifacts", max_iterations_override=1)
    assert state["iterations"][0]["status"] == "accepted"


def test_failed_attempts_excluded_from_analysis(monkeypatch, tmp_path: Path) -> None:
    from autoresearch.harness import bounded_adaptive_search as mod

    calls = {"n": 0}

    def fake_run_safe_campaign(*args, **kwargs):
        calls["n"] += 1
        passed = calls["n"] == 2
        code = "gate_failed_due_to_gsl_signature" if not passed else "ok"
        failures = [] if passed else [{"code": code, "explanation": "gsl"}]
        return {
            "campaign_state": {"summary": {"n_health_failed": 0}, "runs": [], "slices": []},
            "rerun_manifest": {"summary": {"n_runs_selected": 0, "n_slices_selected": 0}},
            "stop_report": {"gates": {"passed": passed, "failures": failures}, "production_validation": passed},
        }

    monkeypatch.setattr(mod, "run_safe_campaign", fake_run_safe_campaign)
    c = _contract(tmp_path)
    state = run_bounded_search(c, execute=False, plan_only=False, output_dir=tmp_path / "artifacts", max_iterations_override=1)
    it = state["iterations"][0]
    assert len(it["attempts"]) >= 2
    assert it["status"] == "accepted"


def test_local_refine_inside_envelope() -> None:
    c = _contract(Path("/tmp"))
    env = validate_envelope(c["search_envelope"])  # type: ignore[arg-type]
    b = validate_budget(c["budget"])  # type: ignore[arg-type]
    sc = build_subcampaign(2, "local_refine", env, {"tan_beta": 56000, "mphi": 292, "lambda1": 4.6, "lambda6": 0.0028}, b)
    assert proposal_inside_envelope(sc.variable_ranges, env)


def test_local_refine_ranges_center_on_canonical_mphi_lambda1() -> None:
    c = _contract(Path("/tmp"))
    env = validate_envelope(c["search_envelope"])  # type: ignore[arg-type]
    b = validate_budget(c["budget"])  # type: ignore[arg-type]
    sc = build_subcampaign(
        2,
        "local_refine",
        env,
        {"tan_beta": 56000.0, "mphi": 292.0, "lambda1": 4.6, "lambda6": 0.0028},
        b,
        refine_cfg={"absolute_widths": {"tan_beta": 2000.0, "mphi": 2.0, "lambda1": 0.2, "lambda6": 0.0002}},
    )
    mphi_mid = (sc.variable_ranges["mphi"][0] + sc.variable_ranges["mphi"][1]) / 2.0
    lam1_mid = (sc.variable_ranges["lambda1"][0] + sc.variable_ranges["lambda1"][1]) / 2.0
    assert abs(mphi_mid - 292.0) < 1e-12
    assert abs(lam1_mid - 4.6) < 1e-12


def test_isolated_metrics_use_only_attempt_rows() -> None:
    rows = [
        {"positivity_ok": "1", "perturbativity_ok": "1", "unitarity_ok": "1", "br_gaga": "1e-4", "width_bb": "2e-4", "total_width": "1e-2", "lambda6": "0.0026", "tan_beta": "55000", "mH": "290", "lambda1": "4.2"},
        {"positivity_ok": "0", "perturbativity_ok": "1", "unitarity_ok": "1", "br_gaga": "9e-4", "width_bb": "9e-4", "total_width": "2e-2", "lambda6": "0.0027", "tan_beta": "56000", "mH": "291", "lambda1": "4.3"},
    ]
    iso = attempt_isolated_summary(rows)
    assert iso["isolated_n_rows"] == 2
    assert iso["isolated_n_triple_ok"] == 1


def test_best_row_anchor_uses_triple_ok_only() -> None:
    rows = [
        {"positivity_ok": "0", "perturbativity_ok": "1", "unitarity_ok": "1", "br_gaga": "9e-3", "width_bb": "9e-4", "total_width": "1e-2", "lambda6": "0.0029", "tan_beta": "59000", "mH": "294", "lambda1": "4.7"},
        {"positivity_ok": "1", "perturbativity_ok": "1", "unitarity_ok": "1", "br_gaga": "1e-4", "width_bb": "1e-4", "total_width": "1e-2", "lambda6": "0.0026", "tan_beta": "54000", "mH": "289", "lambda1": "4.2"},
    ]
    anchor = select_best_anchor_row(rows)
    assert anchor is not None
    assert anchor["lambda6"] == 0.0026
    assert anchor["br_gaga"] == 1e-4
    assert anchor["br_bb"] == 1e-2
    assert anchor["total_width"] == 1e-2


def test_anchor_normalization_accepts_mphi_lam1_aliases() -> None:
    row = {
        "lambda6": "0.0028",
        "tan_beta": "55000",
        "m_phi": "291.2",
        "lam1": "4.33",
        "br_gaga": "3e-4",
        "width_bb": "2e-4",
        "total_width": "1e-2",
    }
    canon = canonicalize_physical_row(row)
    assert canon["mphi"] == 291.2
    assert canon["lambda1"] == 4.33
    assert canon["br_bb"] == 0.02


def test_anchor_normalization_accepts_canonical_mphi_lambda1() -> None:
    row = {
        "lambda6": "0.0028",
        "tan_beta": "55000",
        "mphi": "292.0",
        "lambda1": "4.4",
        "br_gaga": "2e-4",
        "br_bb": "0.12",
        "total_width": "2e-2",
    }
    canon = canonicalize_physical_row(row)
    assert canon["mphi"] == 292.0
    assert canon["lambda1"] == 4.4
    assert canon["br_bb"] == 0.12


def test_select_best_anchor_row_maps_m_phi_and_lam1() -> None:
    rows = [
        {
            "positivity_ok": "1",
            "perturbativity_ok": "1",
            "unitarity_ok": "1",
            "lambda6": "0.0028",
            "tan_beta": "55000",
            "m_phi": "291.2",
            "lam1": "4.33",
            "br_gaga": "3e-4",
            "width_bb": "2e-4",
            "total_width": "1e-2",
        }
    ]
    anchor = select_best_anchor_row(rows)
    assert anchor is not None
    assert anchor["mphi"] == 291.2
    assert anchor["lambda1"] == 4.33


def test_missing_mphi_or_lambda1_prevents_local_refine(monkeypatch, tmp_path: Path) -> None:
    from autoresearch.harness import bounded_adaptive_search as mod

    def fake_rows(*args, **kwargs):
        rows = [
            {
                "positivity_ok": "1",
                "perturbativity_ok": "1",
                "unitarity_ok": "1",
                "lambda6": "0.0026",
                "tan_beta": "54000",
                "br_gaga": "1e-4",
                "width_bb": "1e-4",
                "total_width": "1e-2",
                "m_phi": "289.0",
            }
        ]
        return rows, ["dummy.csv"]

    def fake_run_safe_campaign(*args, **kwargs):
        return {
            "campaign_state": {"summary": {"n_health_failed": 0}, "runs": [], "slices": []},
            "rerun_manifest": {"summary": {"n_runs_selected": 0, "n_slices_selected": 0}},
            "stop_report": {"gates": {"passed": True, "failures": []}, "production_validation": True},
        }

    monkeypatch.setattr(mod, "_scan_rows_for_run", fake_rows)
    monkeypatch.setattr(mod, "run_safe_campaign", fake_run_safe_campaign)

    c = _contract(tmp_path)
    state = run_bounded_search(c, execute=False, plan_only=False, output_dir=tmp_path / "artifacts", max_iterations_override=2)
    it1 = state["iterations"][0]
    assert it1["decision"]["reason"] == "local_refine_unavailable"
    assert it1["decision"]["reason_code"] == "missing_anchor_fields"
    assert it1["decision"]["fallback_strategy"] == "coordinate_sweep"
    assert state["iterations"][1]["strategy"] == "coordinate_sweep"


def test_no_triple_ok_prevents_anchor() -> None:
    rows = [{"positivity_ok": "0", "perturbativity_ok": "1", "unitarity_ok": "1", "br_gaga": "1e-4", "width_bb": "1e-4"}]
    assert select_best_anchor_row(rows) is None


def test_local_refine_clipping_at_envelope_boundary() -> None:
    c = _contract(Path("/tmp"))
    env = validate_envelope(c["search_envelope"])  # type: ignore[arg-type]
    b = validate_budget(c["budget"])  # type: ignore[arg-type]
    sc = build_subcampaign(
        2,
        "local_refine",
        env,
        {"tan_beta": 60000.0, "mphi": 295.0, "lambda1": 4.8, "lambda6": 0.0030},
        b,
        refine_cfg={"absolute_widths": {"tan_beta": 20000.0, "mphi": 20.0, "lambda1": 2.0, "lambda6": 0.001}},
    )
    assert sc.variable_ranges["tan_beta"][1] <= env["tan_beta"][1]
    assert sc.variable_ranges["mphi"][1] <= env["mphi"][1]
    assert sc.variable_ranges["lambda1"][1] <= env["lambda1"][1]


def test_report_contains_isolated_and_cumulative(monkeypatch, tmp_path: Path) -> None:
    from autoresearch.harness import bounded_adaptive_search as mod

    def fake_run_safe_campaign(*args, **kwargs):
        return {
            "campaign_state": {"summary": {"n_health_failed": 0}, "runs": [], "slices": []},
            "rerun_manifest": {"summary": {"n_runs_selected": 0, "n_slices_selected": 0}},
            "stop_report": {"gates": {"passed": True, "failures": []}, "production_validation": True},
        }

    monkeypatch.setattr(mod, "run_safe_campaign", fake_run_safe_campaign)

    c = _contract(tmp_path)
    state = run_bounded_search(c, execute=False, plan_only=False, output_dir=tmp_path / "artifacts", max_iterations_override=1)
    summary = state["iterations"][0]["summary"]
    assert "isolated_summary" in summary
    assert "cumulative_summary" in summary


def test_multiplicative_path_generation_and_lambda1_fixed() -> None:
    c = _contract(Path("/tmp"))
    c["search_envelope"]["lambda1"] = [1.0, 1.0]  # type: ignore[index]
    c["search_envelope"]["tan_beta"] = [30000, 500000]  # type: ignore[index]
    c["search_envelope"]["lambda6"] = [0.0003, 0.0060]  # type: ignore[index]
    env = validate_envelope(c["search_envelope"])  # type: ignore[arg-type]
    b = validate_budget(c["budget"])  # type: ignore[arg-type]
    cfg = {"anchor": {"tan_beta": 60000, "lambda6": 0.0030, "mphi": [250, 330]}, "lambda1_fixed": 1.0, "tan_beta_factor": 1.2, "lambda6_factor": 0.9, "lambda6_local_multiplier_values": [0.97, 0.99, 1.0, 1.01], "tan_beta_local_multiplier_values": [0.96, 0.98, 1.0, 1.02, 1.04], "n_mphi": 50}
    sc0 = build_subcampaign(1, "multiplicative_ctau_path", env, {"tan_beta": 60000, "lambda6": 0.0030}, b, strategy_cfg=cfg)
    sc1 = build_subcampaign(2, "multiplicative_ctau_path", env, {"tan_beta": 60000, "lambda6": 0.0030}, b, strategy_cfg=cfg | {"path_step": 1})
    assert sc0.lambda1_fixed == 1.0
    assert sc0.variable_ranges["lambda1"] == (1.0, 1.0)
    assert abs(sc1.tan_beta_center - (60000 * 1.2)) < 1e-9
    assert abs(sc1.lambda6_center - (0.0030 * 0.9)) < 1e-12


def test_multiplicative_path_stops_on_envelope_escape(tmp_path: Path) -> None:
    c = _contract(tmp_path)
    c["search_envelope"]["lambda1"] = [1.0, 1.0]  # type: ignore[index]
    c["strategy"] = {"name": "multiplicative_ctau_path", "anchor": {"tan_beta": 60000, "lambda6": 0.0030, "mphi": [250, 330]}, "lambda1_fixed": 1.0, "tan_beta_factor": 10.0, "lambda6_factor": 0.9}  # type: ignore[index]
    state = run_bounded_search(c, execute=False, plan_only=True, output_dir=tmp_path / "artifacts", max_iterations_override=2)
    assert state["iterations"][-1]["status"] in {"stopped", "rejected"}


def test_multiplicative_raw_points_enforced_without_lambda1_variation() -> None:
    pts = calculate_raw_points_multiplicative([1, 2, 3, 4, 5], [1, 2, 3, 4], 50, 1)
    assert pts == 1000


def test_ctau_detection_explicit_and_proxy_labels() -> None:
    rows_ctau = [{"ctau": "1.2", "positivity_ok": "1", "unitarity_ok": "1", "perturbativity_ok": "1"}]
    rows_proxy = [{"total_width": "1e-3", "positivity_ok": "1", "unitarity_ok": "1", "perturbativity_ok": "1"}]
    assert detect_ctau_metric_source(rows_ctau) == "explicit_column"
    assert detect_ctau_metric_source(rows_proxy) == "derived_from_total_width"


def test_multiplicative_strategy_artifacts_include_path_factors(tmp_path: Path) -> None:
    c = _contract(tmp_path)
    c["search_envelope"]["lambda1"] = [1.0, 1.0]  # type: ignore[index]
    c["strategy"] = {"name": "multiplicative_ctau_path", "anchor": {"tan_beta": 60000, "lambda6": 0.0030, "mphi": [250, 330]}, "lambda1_fixed": 1.0, "tan_beta_factor": 1.2, "lambda6_factor": 0.9, "lambda6_local_multiplier_values": [0.97, 0.99, 1.0, 1.01], "tan_beta_local_multiplier_values": [0.96, 0.98, 1.0, 1.02, 1.04], "n_mphi": 50, "target_ctau_window_m": [1.0, 10.0], "objective_mode": "ctau_viability"}  # type: ignore[index]
    state = run_bounded_search(c, execute=False, plan_only=True, output_dir=tmp_path / "artifacts", max_iterations_override=1)
    it = state["iterations"][0]
    sub = it.get("subcampaign", {})
    assert sub.get("path_step") == 0
    assert sub.get("tan_beta_factor") == 1.2
    assert sub.get("lambda6_factor") == 0.9


def test_multiplicative_rejects_below_min_unless_operational_smoke() -> None:
    c = _contract(Path("/tmp"))
    c["search_envelope"]["lambda1"] = [1.0, 1.0]  # type: ignore[index]
    env = validate_envelope(c["search_envelope"])  # type: ignore[arg-type]
    b = validate_budget(c["budget"])  # type: ignore[arg-type]
    cfg = {
        "anchor": {"tan_beta": 60000, "lambda6": 0.0030, "mphi": [250, 330]},
        "lambda1_fixed": 1.0,
        "tan_beta_factor": 1.2,
        "lambda6_factor": 0.9,
        "lambda6_local_multiplier_values": [1.0],
        "tan_beta_local_multiplier_values": [1.0],
        "n_mphi": 50,
    }
    sc = build_subcampaign(1, "multiplicative_ctau_path", env, {"tan_beta": 60000, "lambda6": 0.0030}, b, strategy_cfg=cfg)
    assert sc.reason == "insufficient_points"

    b_smoke = dict(b)
    b_smoke["operational_smoke"] = True
    sc_smoke = build_subcampaign(1, "multiplicative_ctau_path", env, {"tan_beta": 60000, "lambda6": 0.0030}, b_smoke, strategy_cfg=cfg)
    assert sc_smoke.reason != "insufficient_points"
    assert sc_smoke.n_lam1 == 1


def test_multiplicative_command_generation_one_per_lambda6(tmp_path: Path) -> None:
    from autoresearch.harness import bounded_adaptive_search as mod

    def fake_run_safe_campaign(*args, **kwargs):
        return {
            "campaign_state": {"summary": {"n_health_failed": 0}, "runs": [], "slices": []},
            "rerun_manifest": {"summary": {"n_runs_selected": 0, "n_slices_selected": 0}},
            "stop_report": {"gates": {"passed": True, "failures": []}, "production_validation": True},
        }

    c = _contract(tmp_path)
    c["search_envelope"]["lambda1"] = [1.0, 1.0]  # type: ignore[index]
    c["search_envelope"]["tan_beta"] = [30000, 500000]  # type: ignore[index]
    c["search_envelope"]["lambda6"] = [0.0003, 0.0060]  # type: ignore[index]
    c["strategy"] = {
        "name": "multiplicative_ctau_path",
        "anchor": {"tan_beta": 60000, "lambda6": 0.0030, "mphi": [250, 330]},
        "lambda1_fixed": 1.0,
        "tan_beta_factor": 1.2,
        "lambda6_factor": 0.9,
        "lambda6_local_multiplier_values": [0.97, 0.99, 1.0, 1.01],
        "tan_beta_local_multiplier_values": [0.96, 0.98, 1.0, 1.02, 1.04],
        "n_mphi": 50,
    }  # type: ignore[index]
    mod.run_safe_campaign = fake_run_safe_campaign
    state = run_bounded_search(c, execute=False, plan_only=True, output_dir=tmp_path / "artifacts", max_iterations_override=1)
    assert state["iterations"][0]["subcampaign"]["estimated_raw_points"] == 1000

    cmds = (tmp_path / "artifacts" / "commands.jsonl").read_text().strip().splitlines()
    assert len(cmds) == 4
    import json
    for line in cmds:
        rec = json.loads(line)
        cmd = rec["command"]
        assert "--tanbeta" in cmd
        assert "--lam1-min" in cmd and cmd[cmd.index("--lam1-min") + 1] in {"1", "1.0"}
        assert "--lam1-max" in cmd and cmd[cmd.index("--lam1-max") + 1] in {"1", "1.0"}
        assert "--n-lam1" in cmd and cmd[cmd.index("--n-lam1") + 1] == "1"
        assert "--force" not in cmd


def test_ctau_summary_explicit_derived_and_unavailable() -> None:
    from autoresearch.harness.bounded_adaptive_search import _ctau_summary

    rows_explicit = [{"ctau": "2.0", "total_width": "1e-2", "positivity_ok": "1", "perturbativity_ok": "1", "unitarity_ok": "1"}]
    s1 = _ctau_summary(rows_explicit)
    assert s1["ctau_metric_source"] == "explicit_column"
    assert s1["ctau_m_q50"] == 2.0

    rows_proxy = [{"total_width": "1e-3", "positivity_ok": "1", "perturbativity_ok": "1", "unitarity_ok": "1"}]
    s2 = _ctau_summary(rows_proxy)
    assert s2["ctau_metric_source"] == "derived_from_total_width"
    assert s2["ctau_m_q50"] is not None

    rows_none = [{"br_gaga": "1e-4", "positivity_ok": "1", "perturbativity_ok": "1", "unitarity_ok": "1"}]
    s3 = _ctau_summary(rows_none)
    assert s3["ctau_metric_source"] == "unavailable"


def test_canonicalize_normalizes_alias_and_derives_br_bb() -> None:
    row = {"m_phi": "300", "lam1": "1.0", "lambda6": "0.003", "tan_beta": "60000", "br_gaga": "1e-4", "width_bb": "2e-4", "total_width": "1e-2"}
    c = canonicalize_physical_row(row)
    assert c["mphi"] == 300.0
    assert c["lambda1"] == 1.0
    assert abs(c["br_bb"] - 0.02) < 1e-12




def test_header_only_csv_counted_zero_data_rows_and_excluded(tmp_path: Path) -> None:
    good = tmp_path / "good.csv"
    good.write_text("a,b\n1,2\n3,4\n", encoding="utf-8")
    hdr = tmp_path / "header_only.csv"
    hdr.write_text("a,b\n", encoding="utf-8")

    rows = [
        {"positivity_ok": "1", "perturbativity_ok": "1", "unitarity_ok": "1", "total_width": "1e-12"},
        {"positivity_ok": "1", "perturbativity_ok": "1", "unitarity_ok": "1", "total_width": "2e-12"},
    ]
    acc = derive_step_accounting_v3(
        csv_paths=[str(good), str(hdr)],
        rows=rows,
        planned_points=21,
        attempted_points=21,
        failure_reason_codes=["gate_failed_due_to_header_only_csv"],
    )
    assert acc["raw_csv_rows"] == 2
    assert acc["accepted_csv_rows"] == 2
    assert acc["header_only_csv_count"] == 1
    assert acc["accepted_physics_points"] == 2
    assert acc["learning_points"] == 2


def test_accepted_physics_points_cannot_exceed_accepted_csv_rows(tmp_path: Path) -> None:
    one = tmp_path / "one.csv"
    one.write_text("a,b\n1,2\n", encoding="utf-8")
    rows = [
        {"positivity_ok": "1", "perturbativity_ok": "1", "unitarity_ok": "1", "total_width": "1e-12"},
        {"positivity_ok": "1", "perturbativity_ok": "1", "unitarity_ok": "1", "total_width": "2e-12"},
    ]
    try:
        derive_step_accounting_v3(
            csv_paths=[str(one)],
            rows=rows,
            planned_points=10,
            attempted_points=10,
            failure_reason_codes=[],
        )
        raised = False
    except ValueError as exc:
        raised = "learning_points exceeds accepted_physics_points" in str(exc)
    assert raised


def test_planned_points_do_not_become_executed_points(tmp_path: Path) -> None:
    one = tmp_path / "slice.csv"
    one.write_text("a,b\n1,2\n", encoding="utf-8")
    rows = [{"positivity_ok": "1", "perturbativity_ok": "1", "unitarity_ok": "1", "total_width": "1e-12"}]
    acc = derive_step_accounting_v3(
        csv_paths=[str(one)],
        rows=rows,
        planned_points=21000,
        attempted_points=21000,
        failure_reason_codes=[],
    )
    assert acc["planned_points"] == 21000
    assert acc["accepted_csv_rows"] == 1
    assert acc["accepted_physics_points"] == 1


def test_sensitivity_probe_contract_template_grid_probe_compatible() -> None:
    env = {
        "sin_ba": [1.0, 1.0], "lambda7": [0.0, 0.0], "lambda1": [1.0, 1.0],
        "lambda6": [0.0003, 0.0060], "tan_beta": [30000, 500000], "mphi": [160, 350], "mA": [300, 500],
    }
    tm = materialize_template(
        template_name="sensitivity_probe",
        overrides={
            "variable": "tan_beta",
            "relative_steps": [-0.10, 0.0, 0.10],
            "anchor": {"tan_beta": 126904.3, "lambda6": 0.0019, "mA": 500.0, "lambda1": 1.0, "mphi_min": 180.0, "mphi_max": 220.0},
        },
        envelope=env,
        min_points=1000,
        max_points=3000,
        iteration=1,
        max_total_points=9000,
        tan_beta_center_default=124416,
        lambda6_center_default=0.0019683,
    )
    sc = build_subcampaign(
        1,
        "template_grid_probe",
        {k: tuple(v) for k, v in tm.contract["search_envelope"].items()},
        anchor={"lambda6": 0.0019, "tan_beta": 126904.3, "mphi": 200.0, "lambda1": 1.0, "mA": 500.0, "lambda7": 0.0, "sin_ba": 1.0},
        budget=tm.contract["budget"],
        strategy_cfg=tm.contract["strategy"],
    )
    assert sc.strategy == "template_grid_probe"
    assert sc.n_lam1 == 1
    cmds = build_orchestrator_commands_for_subcampaign(sc, tm.contract["runtime"], 4, "r")
    assert len(cmds) == len(tm.contract["strategy"]["mA_values"]) * len(tm.contract["strategy"]["lambda6_values"])

def test_progress_files_written(monkeypatch, tmp_path: Path) -> None:
    from autoresearch.harness import bounded_adaptive_search as mod

    def fake_run_safe_campaign(*args, **kwargs):
        return {
            "campaign_state": {"summary": {"n_health_failed": 0}, "runs": [], "slices": []},
            "rerun_manifest": {"summary": {"n_runs_selected": 0, "n_slices_selected": 0}},
            "stop_report": {"gates": {"passed": True, "failures": []}, "production_validation": True},
        }

    monkeypatch.setattr(mod, "run_safe_campaign", fake_run_safe_campaign)
    c = _contract(tmp_path)
    outdir = tmp_path / "artifacts"
    run_bounded_search(c, execute=False, plan_only=False, output_dir=outdir, max_iterations_override=1)
    assert (outdir / "progress.jsonl").exists()
    assert (outdir / "progress.md").exists()
    text = (outdir / "progress.md").read_text(encoding="utf-8")
    assert "best_ctau_m" in text


def test_subcampaign_run_identity_suffix_propagates_to_run_names() -> None:
    c = _contract(Path("/tmp"))
    env = validate_envelope(c["search_envelope"])  # type: ignore[arg-type]
    b = validate_budget(c["budget"])  # type: ignore[arg-type]
    cfg = {
        "template_name": "ctau_mphi_mA_refine",
        "lambda1_fixed": 1.0,
        "n_lam1": 1,
        "mA_values": [300.0],
        "lambda6_center": 0.0027,
        "lambda6_local_multiplier_values": [1.0],
        "tan_beta_center": 55000.0,
        "tan_beta_local_multiplier_values": [1.0],
        "mphi_min": 285.0,
        "mphi_max": 295.0,
        "n_mphi": 1000,
        "run_identity_suffix": "hdeadbeef_cycleX_iterY",
    }
    sc = build_subcampaign(1, "template_grid_probe", env, {"tan_beta": 55000.0, "mphi": 290.0, "lambda1": 1.0, "lambda6": 0.0027, "mA": 300.0, "lambda7": 0.0, "sin_ba": 1.0}, b, strategy_cfg=cfg)
    assert "hdeadbeef" in sc.subcampaign_id
