from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import patch

import scripts.jev_run_once as jev_run_once
from scripts.gemini_contract_templates import DEFAULT_LAMBDA6_CENTER, DEFAULT_TAN_BETA_CENTER, materialize_template
from scripts.jev_policy import JevDecision


def _cfg() -> dict:
    return {
        "baseline_best_ctau_m": 5.390644846600711e-4,
        "target_ctau_m": 1.0781289693201422e-3,
        "max_iterations": 3,
        "max_walltime_minutes": 60,
        "max_points_total": 9000,
        "min_points_per_subrun": 1000,
        "max_points_per_subrun": 3000,
        "allowed_actions": ["run_template", "summarize_only", "stop"],
        "forbidden_actions": ["edit_cpp", "edit_scoring", "edit_triple_ok", "delete_artifacts", "widen_envelope", "broad_scan"],
        "envelope": {
            "sin_ba": [1.0, 1.0], "lambda7": [0.0, 0.0], "lambda1": [1.0, 1.0], "lambda6": [0.0003, 0.0060],
            "tan_beta": [30000, 500000], "mphi": [160, 350], "mA": [300, 500],
        },
    }


def test_materialize_template_without_explicit_centers_uses_canonical_defaults() -> None:
    # No tan_beta_center_default/lambda6_center_default passed: this is exactly the
    # call shape jev_worker_config.sample.json produces (it carries no center keys).
    cfg = _cfg()
    tm = materialize_template(
        template_name="ctau_mphi_mA_refine",
        overrides={},
        envelope=cfg["envelope"],
        min_points=cfg["min_points_per_subrun"],
        max_points=cfg["max_points_per_subrun"],
        iteration=1,
        max_total_points=cfg["max_points_total"],
    )
    assert tm.contract["strategy"]["tan_beta_center"] == DEFAULT_TAN_BETA_CENTER
    assert tm.contract["strategy"]["lambda6_center"] == DEFAULT_LAMBDA6_CENTER


def test_materialize_and_plancheck_runs_bounded_search_plan_only(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(jev_run_once, "HARNESS_DIR", tmp_path / "jev_harness")
    outcome = jev_run_once.materialize_and_plancheck(
        {"action": "run_template", "template_name": "control_box", "template_overrides": {}},
        _cfg(),
        iteration=1,
        cycle_id="cycle_test",
        iteration_id="iter001_test",
    )
    assert outcome["template_name"] == "control_box"
    assert outcome["plancheck_ok"] is True
    assert 1000 <= outcome["planned_points"] <= 3000
    assert Path(outcome["contract_path"]).exists()
    assert Path(outcome["plan_dir"]).exists()


def test_materialize_and_plancheck_rejects_forbidden_template() -> None:
    try:
        jev_run_once.materialize_and_plancheck(
            {"action": "run_template", "template_name": "edit_cpp", "template_overrides": {}},
            _cfg(),
            iteration=1,
            cycle_id="cycle_test",
            iteration_id="iter001_test",
        )
        assert False, "expected ValueError for unrecognized template"
    except ValueError as exc:
        assert "unknown template" in str(exc)


def _decision(choice: str) -> JevDecision:
    action = (
        {"action": "stop", "rationale": "r", "safety_checks": [], "stop_conditions": ["jev_selected_stop"]}
        if choice == "stop"
        else {"action": "run_template", "template_name": choice, "template_overrides": {}, "rationale": "r", "safety_checks": [], "stop_conditions": []}
    )
    return JevDecision(
        worker_action=action,
        request_payload={"state": {}, "questions": {}},
        raw_response={"model": "jev-test", "answers": {}, "usage": {"input_tokens": 1, "output_tokens": 1}},
    )


def test_main_writes_decision_event_with_materialized_outcome(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(jev_run_once, "HARNESS_DIR", tmp_path / "jev_harness")
    cfg_path = tmp_path / "cfg.json"
    cfg_path.write_text(json.dumps(_cfg()), encoding="utf-8")
    events_path = tmp_path / "events.jsonl"

    with patch("scripts.jev_run_once.choose_with_jev", return_value=_decision("control_box")):
        monkeypatch.setattr(sys, "argv", ["jev_run_once.py", "--config", str(cfg_path), "--events", str(events_path)])
        assert jev_run_once.main() == 0

    lines = events_path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1
    event = json.loads(lines[0])
    assert event["schema_version"] == "dihiggs_policy_decision_v1"
    assert event["selected_action"]["template_name"] == "control_box"
    assert event["outcome"]["status"] == "planned"
    assert event["outcome"]["plancheck_ok"] is True


def test_main_decision_fixture_never_calls_jev(tmp_path, monkeypatch) -> None:
    """--decision-fixture must go through response_to_worker_action on a saved raw
    response and never call choose_with_jev (i.e. never touch the network/proxy)."""
    monkeypatch.setattr(jev_run_once, "HARNESS_DIR", tmp_path / "jev_harness")
    cfg_path = tmp_path / "cfg.json"
    cfg_path.write_text(json.dumps(_cfg()), encoding="utf-8")
    events_path = tmp_path / "events.jsonl"

    fixture_path = tmp_path / "fixture.json"
    fixture_path.write_text(
        json.dumps(
            {
                "model": "jev-fixture-test",
                "answers": {
                    "next_operator": {
                        "type": "choice",
                        "choice": "control_box",
                        "confidence": 0.9,
                        "probabilities": {"control_box": 0.9, "stop": 0.1},
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    def _boom(*args, **kwargs):
        raise AssertionError("choose_with_jev must not be called on the --decision-fixture path")

    with patch("scripts.jev_run_once.choose_with_jev", side_effect=_boom):
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "jev_run_once.py",
                "--config",
                str(cfg_path),
                "--events",
                str(events_path),
                "--decision-fixture",
                str(fixture_path),
            ],
        )
        assert jev_run_once.main() == 0

    event = json.loads(events_path.read_text(encoding="utf-8").splitlines()[0])
    assert event["backend"] == "offline_fixture"
    assert event["decision_fixture_path"] == str(fixture_path)
    assert event["selected_action"]["template_name"] == "control_box"
    assert event["outcome"]["plancheck_ok"] is True


def test_main_stop_choice_skips_materialization(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(jev_run_once, "HARNESS_DIR", tmp_path / "jev_harness")
    cfg_path = tmp_path / "cfg.json"
    cfg_path.write_text(json.dumps(_cfg()), encoding="utf-8")
    events_path = tmp_path / "events.jsonl"

    with patch("scripts.jev_run_once.choose_with_jev", return_value=_decision("stop")):
        monkeypatch.setattr(sys, "argv", ["jev_run_once.py", "--config", str(cfg_path), "--events", str(events_path)])
        assert jev_run_once.main() == 0

    event = json.loads(events_path.read_text(encoding="utf-8").splitlines()[0])
    assert event["outcome"] == {"status": "stopped_by_jev"}
    assert not (tmp_path / "jev_harness" / "candidates").exists()
