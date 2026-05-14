from __future__ import annotations

import json
from pathlib import Path

from scripts.gemini_contract_templates import TEMPLATES, compute_template_points, materialize_template
from scripts.gemini_worker_loop import run_worker, validate_gemini_json
from autoresearch.harness.bounded_adaptive_search import build_subcampaign, build_orchestrator_commands_for_subcampaign


def _cfg() -> dict:
    return {
        "baseline_best_ctau_m": 5.390644846600711e-4,
        "target_ctau_m": 1.0781289693201422e-3,
        "max_iterations": 3,
        "max_walltime_minutes": 60,
        "max_points_total": 9000,
        "min_points_per_subrun": 1000,
        "max_points_per_subrun": 3000,
        "gemini_retry_on_invalid_json": 1,
        "allowed_actions": ["run_template", "summarize_only", "stop"],
        "forbidden_actions": ["edit_cpp", "edit_scoring", "edit_triple_ok", "delete_artifacts", "widen_envelope", "broad_scan"],
        "envelope": {
            "sin_ba": [1.0, 1.0], "lambda7": [0.0, 0.0], "lambda1": [1.0, 1.0], "lambda6": [0.0003, 0.0060],
            "tan_beta": [30000, 500000], "mphi": [160, 350], "mA": [300, 500]
        },
        "allow_direct_run_contract": False,
        "tan_beta_center": 124416,
        "lambda6_center": 0.0019683,
    }


def _run_template_obj() -> dict:
    return {
        "action": "run_template",
        "template_name": "ctau_mphi_mA_refine",
        "template_overrides": {"mphi_min": 180, "mphi_max": 220, "n_mphi": 40},
        "rationale": "test",
        "safety_checks": ["inside_envelope"],
        "stop_conditions": ["none"],
    }


def test_ctau_template_exact_1080_points() -> None:
    assert compute_template_points(TEMPLATES["ctau_mphi_mA_refine"]) == 1080


def test_templates_reject_bounds_and_forbidden() -> None:
    cfg = _cfg()
    # Above max must reject.
    try:
        materialize_template(
            template_name="control_box",
            overrides={"n_mphi": 5000},
            envelope=cfg["envelope"],
            min_points=1000,
            max_points=3000,
            iteration=1,
            max_total_points=9000,
            tan_beta_center_default=124416,
            lambda6_center_default=0.0019683,
        )
        assert False
    except ValueError as e:
        assert "above" in str(e)
    # Forbidden override must reject.
    try:
        materialize_template(
            template_name="control_box",
            overrides={"lambda1": [1, 1]},
            envelope=cfg["envelope"],
            min_points=1000,
            max_points=3000,
            iteration=1,
            max_total_points=9000,
            tan_beta_center_default=124416,
            lambda6_center_default=0.0019683,
        )
        assert False
    except ValueError as e:
        assert "forbidden" in str(e)


def test_low_point_template_auto_expands_to_min() -> None:
    cfg = _cfg()
    tm = materialize_template(
        template_name="control_box",
        overrides={"n_mphi": 10},
        envelope=cfg["envelope"],
        min_points=1000,
        max_points=3000,
        iteration=1,
        max_total_points=9000,
        tan_beta_center_default=124416,
        lambda6_center_default=0.0019683,
    )
    assert tm.real_points >= 1000


def test_template_grid_probe_dimensions_not_lost() -> None:
    cfg = _cfg()
    tm = materialize_template(
        template_name="ctau_mphi_mA_refine",
        overrides={},
        envelope=cfg["envelope"],
        min_points=1000,
        max_points=3000,
        iteration=1,
        max_total_points=9000,
        tan_beta_center_default=124416,
        lambda6_center_default=0.0019683,
    )
    assert tm.real_points == 1080
    sc = build_subcampaign(
        1,
        "template_grid_probe",
        {k: tuple(v) for k, v in tm.contract["search_envelope"].items()},
        anchor={"lambda6": 0.0019683, "tan_beta": 124416, "mphi": 200, "lambda1": 1.0, "mA": 450, "lambda7": 0.0, "sin_ba": 1.0},
        budget=tm.contract["budget"],
        strategy_cfg=tm.contract["strategy"],
    )
    assert sc.estimated_raw_points == 1080
    assert len(sc.mA_values or []) == 3
    assert len(sc.lambda6_local_values or []) == 3
    assert len(sc.tanbeta_values) == 3
    assert sc.n_mphi == 40
    assert sc.n_lam1 == 1
    assert sc.lambda1_fixed == 1.0


def test_template_grid_probe_command_count_per_ma_l6() -> None:
    cfg = _cfg()
    tm = materialize_template(
        template_name="ctau_mphi_mA_refine",
        overrides={},
        envelope=cfg["envelope"],
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
        anchor={"lambda6": 0.0019683, "tan_beta": 124416, "mphi": 200, "lambda1": 1.0, "mA": 450, "lambda7": 0.0, "sin_ba": 1.0},
        budget=tm.contract["budget"],
        strategy_cfg=tm.contract["strategy"],
    )
    cmds = build_orchestrator_commands_for_subcampaign(sc, tm.contract["runtime"], 4, "r")
    assert len(cmds) == 9


def test_wrapper_accepts_run_template_action_and_rejects_run_contract() -> None:
    v = validate_gemini_json(_run_template_obj(), _cfg(), {"executed_points_total": 0})
    assert v["valid"]
    bad = _run_template_obj(); bad["action"] = "run_contract"
    v2 = validate_gemini_json(bad, _cfg(), {"executed_points_total": 0})
    assert not v2["valid"]


def test_plancheck_points_equal_template_points(monkeypatch, tmp_path: Path) -> None:
    from scripts import gemini_worker_loop as mod

    cfg = _cfg(); cfg_path = tmp_path / "cfg.json"; cfg_path.write_text(json.dumps(cfg), encoding="utf-8")

    def fake_call(prompt: str):
        return 0, json.dumps(_run_template_obj()), ""

    def fake_run(cmd, cwd=None, text=None, capture_output=None):
        class R:
            returncode = 0
            stdout = ""
            stderr = ""
        if "--plan-only" in cmd:
            cpath = Path(cmd[cmd.index("--contract") + 1])
            contract = json.loads(cpath.read_text(encoding="utf-8"))
            pts = int(contract["strategy"]["expected_raw_points"])
            pdir = Path(cmd[cmd.index("--output-dir") + 1]); pdir.mkdir(parents=True, exist_ok=True)
            (pdir / "subcampaigns.jsonl").write_text(json.dumps({"estimated_raw_points": pts, "status": "planned"}) + "\n", encoding="utf-8")
            (pdir / "adaptive_state.json").write_text(json.dumps({"total_points_executed": pts}), encoding="utf-8")
        return R()

    monkeypatch.setattr(mod, "HARNESS_DIR", tmp_path / "h")
    monkeypatch.setattr(mod, "call_gemini", fake_call)
    monkeypatch.setattr(mod.subprocess, "run", fake_run)

    run_worker(cfg_path, dry_run=True, iterations_override=1)
    rep = json.loads((tmp_path / "h" / "latest_validation_report.json").read_text(encoding="utf-8"))
    assert rep["template_real_points"] == 1080
    assert rep["contract_estimated_points"] == 1080
    assert rep["planned_points"] == 1080
    assert "template_plan_mismatch" not in rep["plancheck_rejection_reasons"]


def test_template_plan_mismatch_blocks_execution(monkeypatch, tmp_path: Path) -> None:
    from scripts import gemini_worker_loop as mod

    cfg = _cfg(); cfg_path = tmp_path / "cfg.json"; cfg_path.write_text(json.dumps(cfg), encoding="utf-8")
    exec_calls = {"n": 0}

    def fake_call(prompt: str):
        return 0, json.dumps(_run_template_obj()), ""

    def fake_run(cmd, cwd=None, text=None, capture_output=None):
        class R:
            returncode = 0
            stdout = ""
            stderr = ""
        if "--execute" in cmd:
            exec_calls["n"] += 1
        if "--plan-only" in cmd:
            pdir = Path(cmd[cmd.index("--output-dir") + 1]); pdir.mkdir(parents=True, exist_ok=True)
            (pdir / "subcampaigns.jsonl").write_text(json.dumps({"estimated_raw_points": 1000, "status": "planned"}) + "\n", encoding="utf-8")
            (pdir / "adaptive_state.json").write_text(json.dumps({"total_points_executed": 1000}), encoding="utf-8")
        return R()

    monkeypatch.setattr(mod, "HARNESS_DIR", tmp_path / "h")
    monkeypatch.setattr(mod, "call_gemini", fake_call)
    monkeypatch.setattr(mod.subprocess, "run", fake_run)

    run_worker(cfg_path, dry_run=False, iterations_override=1)
    rep = json.loads((tmp_path / "h" / "latest_validation_report.json").read_text(encoding="utf-8"))
    assert "template_plan_mismatch" in rep["plancheck_rejection_reasons"]
    assert exec_calls["n"] == 0


def test_dry_run_no_execute(monkeypatch, tmp_path: Path) -> None:
    from scripts import gemini_worker_loop as mod

    cfg = _cfg(); cfg_path = tmp_path / "cfg.json"; cfg_path.write_text(json.dumps(cfg), encoding="utf-8")
    exec_calls = {"n": 0}

    def fake_call(prompt: str):
        return 0, json.dumps(_run_template_obj()), ""

    def fake_run(cmd, cwd=None, text=None, capture_output=None):
        class R:
            returncode = 0
            stdout = ""
            stderr = ""
        if "--execute" in cmd:
            exec_calls["n"] += 1
        if "--plan-only" in cmd:
            pdir = Path(cmd[cmd.index("--output-dir") + 1]); pdir.mkdir(parents=True, exist_ok=True)
            (pdir / "subcampaigns.jsonl").write_text(json.dumps({"estimated_raw_points": 1080, "status": "planned"}) + "\n", encoding="utf-8")
            (pdir / "adaptive_state.json").write_text(json.dumps({"total_points_executed": 1080}), encoding="utf-8")
        return R()

    monkeypatch.setattr(mod, "HARNESS_DIR", tmp_path / "h")
    monkeypatch.setattr(mod, "call_gemini", fake_call)
    monkeypatch.setattr(mod.subprocess, "run", fake_run)

    run_worker(cfg_path, dry_run=True, iterations_override=1)
    assert exec_calls["n"] == 0


def test_unique_run_identity_and_fresh_cycle_namespace(monkeypatch, tmp_path: Path) -> None:
    from scripts import gemini_worker_loop as mod

    cfg = _cfg(); cfg_path = tmp_path / "cfg.json"; cfg_path.write_text(json.dumps(cfg), encoding="utf-8")
    seen_contracts: list[dict] = []

    def fake_call(prompt: str):
        return 0, json.dumps(_run_template_obj()), ""

    def fake_run(cmd, cwd=None, text=None, capture_output=None):
        class R:
            returncode = 0
            stdout = ""
            stderr = ""
        if "--plan-only" in cmd:
            cpath = Path(cmd[cmd.index("--contract") + 1])
            contract = json.loads(cpath.read_text(encoding="utf-8"))
            seen_contracts.append(contract)
            pdir = Path(cmd[cmd.index("--output-dir") + 1]); pdir.mkdir(parents=True, exist_ok=True)
            pts = int(contract["strategy"]["expected_raw_points"])
            (pdir / "subcampaigns.jsonl").write_text(json.dumps({"estimated_raw_points": pts, "status": "planned"}) + "\n", encoding="utf-8")
            (pdir / "adaptive_state.json").write_text(json.dumps({"total_points_executed": pts}), encoding="utf-8")
        return R()

    monkeypatch.setattr(mod, "HARNESS_DIR", tmp_path / "h")
    monkeypatch.setattr(mod, "call_gemini", fake_call)
    monkeypatch.setattr(mod.subprocess, "run", fake_run)

    run_worker(cfg_path, dry_run=True, iterations_override=1)
    assert seen_contracts, "expected contract materialization"
    c = seen_contracts[0]
    suffix = c["strategy"].get("run_identity_suffix")
    assert isinstance(suffix, str) and suffix.startswith("h")
    assert c["runtime"]["campaign"].startswith("gemini_worker_iter01_cycle_")


def test_contamination_flag_when_isolated_rows_exceed_executed(monkeypatch, tmp_path: Path) -> None:
    from scripts import gemini_worker_loop as mod

    cfg = _cfg(); cfg_path = tmp_path / "cfg.json"; cfg_path.write_text(json.dumps(cfg), encoding="utf-8")

    def fake_call(prompt: str):
        return 0, json.dumps(_run_template_obj()), ""

    def fake_run(cmd, cwd=None, text=None, capture_output=None):
        class R:
            returncode = 0
            stdout = ""
            stderr = ""
        outdir = Path(cmd[cmd.index("--output-dir") + 1])
        outdir.mkdir(parents=True, exist_ok=True)
        if "--plan-only" in cmd:
            (outdir / "subcampaigns.jsonl").write_text(json.dumps({"estimated_raw_points": 1080, "status": "planned"}) + "\n", encoding="utf-8")
            (outdir / "adaptive_state.json").write_text(json.dumps({"total_points_executed": 1080}), encoding="utf-8")
        if "--execute" in cmd:
            (outdir / "adaptive_state.json").write_text(json.dumps({
                "total_points_executed": 1000,
                "iterations": [{"summary": {"isolated_summary": {"isolated_n_rows": 1200}}}]
            }), encoding="utf-8")
        return R()

    monkeypatch.setattr(mod, "HARNESS_DIR", tmp_path / "h")
    monkeypatch.setattr(mod, "call_gemini", fake_call)
    monkeypatch.setattr(mod.subprocess, "run", fake_run)

    run_worker(cfg_path, dry_run=False, iterations_override=1)
    rep = json.loads((tmp_path / "h" / "latest_validation_report.json").read_text(encoding="utf-8"))
    state = json.loads((tmp_path / "h" / "worker_state.json").read_text(encoding="utf-8"))
    assert rep["contaminated_summary"] is True
    assert state["status"] == "stopped_contamination_detected"
    assert state["contaminated_summary_count"] == 1
    assert state["last_iteration_summary_for_prompt"] is None


def test_sensitivity_probe_registered() -> None:
    assert "sensitivity_probe" in TEMPLATES


def test_sensitivity_probe_tan_beta_relative_steps() -> None:
    cfg = _cfg()
    tm = materialize_template(
        template_name="sensitivity_probe",
        overrides={
            "variable": "tan_beta",
            "relative_steps": [-0.10, 0.0, 0.10],
            "anchor": {
                "tan_beta": 126904.3,
                "lambda6": 0.0019,
                "mA": 500.0,
                "lambda1": 1.0,
                "mphi_min": 180.0,
                "mphi_max": 220.0,
            },
        },
        envelope=cfg["envelope"],
        min_points=1000,
        max_points=3000,
        iteration=1,
        max_total_points=9000,
        tan_beta_center_default=124416,
        lambda6_center_default=0.0019683,
    )
    vals = tm.contract["strategy"]["tan_beta_values"]
    assert len(vals) == 3
    assert vals[1] == 126904.3
    assert tm.contract["strategy"]["lambda6_values"] == [0.0019]
    assert tm.contract["strategy"]["mA_values"] == [500.0]
    assert tm.contract["strategy"]["n_lam1"] == 1
    assert tm.real_points >= 1000


def test_sensitivity_probe_lambda6_relative_steps() -> None:
    cfg = _cfg()
    tm = materialize_template(
        template_name="sensitivity_probe",
        overrides={
            "variable": "lambda6",
            "relative_steps": [-0.10, 0.0, 0.10],
            "anchor": {"tan_beta": 126904.3, "lambda6": 0.0019, "mA": 500.0, "lambda1": 1.0, "mphi_min": 180.0, "mphi_max": 220.0},
        },
        envelope=cfg["envelope"],
        min_points=1000,
        max_points=3000,
        iteration=1,
        max_total_points=9000,
        tan_beta_center_default=124416,
        lambda6_center_default=0.0019683,
    )
    vals = tm.contract["strategy"]["lambda6_values"]
    assert len(vals) == 3
    assert vals[1] == 0.0019
    assert tm.contract["strategy"]["tan_beta_values"] == [126904.3]
    assert tm.contract["strategy"]["n_lam1"] == 1


def test_sensitivity_probe_mA_relative_steps_respects_envelope() -> None:
    cfg = _cfg()
    tm = materialize_template(
        template_name="sensitivity_probe",
        overrides={
            "variable": "mA",
            "relative_steps": [-0.10, 0.0, 0.10],
            "anchor": {"tan_beta": 126904.3, "lambda6": 0.0019, "mA": 450.0, "lambda1": 1.0, "mphi_min": 180.0, "mphi_max": 220.0},
        },
        envelope=cfg["envelope"],
        min_points=1000,
        max_points=3000,
        iteration=1,
        max_total_points=9000,
        tan_beta_center_default=124416,
        lambda6_center_default=0.0019683,
    )
    mAs = tm.contract["strategy"]["mA_values"]
    assert min(mAs) >= cfg["envelope"]["mA"][0]
    assert max(mAs) <= cfg["envelope"]["mA"][1]


def test_sensitivity_probe_lambda1_rejected_by_default() -> None:
    cfg = _cfg()
    try:
        materialize_template(
            template_name="sensitivity_probe",
            overrides={"variable": "lambda1", "relative_steps": [-0.10, 0.0, 0.10], "anchor": {"lambda1": 1.0}},
            envelope=cfg["envelope"],
            min_points=1000,
            max_points=3000,
            iteration=1,
            max_total_points=9000,
            tan_beta_center_default=124416,
            lambda6_center_default=0.0019683,
        )
        assert False
    except ValueError as e:
        assert "lambda1 variation disabled" in str(e)
