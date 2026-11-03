from __future__ import annotations

import json
from unittest.mock import patch

from scripts.jev_policy import (
    JevError,
    build_decision_event,
    build_jev_request,
    build_scientific_state,
    call_jev,
    choose_with_gemini,
    choose_with_jev,
    response_to_worker_action,
)


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
            "sin_ba": [1.0, 1.0],
            "lambda7": [0.0, 0.0],
            "lambda1": [1.0, 1.0],
            "lambda6": [0.0003, 0.0060],
            "tan_beta": [30000, 500000],
            "mphi": [160, 350],
            "mA": [300, 500],
        },
    }


def _state() -> dict:
    return {
        "best_ctau_m": 6.0e-4,
        "executed_points_total": 2200,
        "iterations_completed": 1,
        "no_improvement_count": 0,
        "history": [{"iteration": 1, "result": "executed", "action": "run_template"}],
        "last_iteration_summary_for_prompt": {
            "iteration": 1,
            "executed_points": 1080,
            "isolated_summary_rows": 12,
            "planned_points": 1080,
            "template_name": "ctau_mphi_mA_refine",
        },
    }


def _response(choice: str = "mphi_lower_edge_probe") -> dict:
    probs = {
        "ctau_mphi_mA_refine": 0.10,
        "mphi_lower_edge_probe": 0.55,
        "mA_upper_edge_probe": 0.10,
        "sensitivity_probe": 0.10,
        "control_box": 0.05,
        "stop": 0.10,
    }
    return {
        "model": "jev-test",
        "answers": {
            "next_operator": {
                "type": "choice",
                "choice": choice,
                "confidence": 0.82,
                "probabilities": probs,
            },
            "local_saturation": {
                "type": "score",
                "score": 1.7,
                "confidence": 0.7,
                "legend": {"0": "none", "1": "weak", "2": "mixed", "3": "substantial", "4": "strong"},
                "probabilities": {"0": 0.1, "1": 0.3, "2": 0.5, "3": 0.1, "4": 0.0},
            },
            "exploitation_value": {"type": "score", "score": 2.4, "confidence": 0.7, "legend": {}, "probabilities": {}},
            "diversification_value": {"type": "score", "score": 3.1, "confidence": 0.8, "legend": {}, "probabilities": {}},
            "frontier_escalation": {"type": "noul", "noul": 0.18},
        },
        "usage": {"input_tokens": 900, "output_tokens": 40},
    }


def test_state_contains_budget_constraints_and_recent_outcome() -> None:
    state = build_scientific_state(_cfg(), _state())
    assert state["budget"]["remaining_points"] == 6800
    assert state["objective"]["best_ctau_m"] == 6.0e-4
    assert state["hard_constraints"]["search_envelope"]["mA"] == [300, 500]
    assert state["search_context"]["last_uncontaminated_iteration"]["planned_points"] == 1080


def test_request_uses_only_typed_systemone_questions() -> None:
    req = build_jev_request(_cfg(), _state())
    assert "model" not in req
    assert req["questions"]["next_operator"]["type"] == "choice"
    assert req["questions"]["local_saturation"]["type"] == "score"
    assert req["questions"]["frontier_escalation"]["type"] == "noul"
    assert set(req["questions"]["next_operator"]["criteria"]) == {
        "ctau_mphi_mA_refine",
        "mphi_lower_edge_probe",
        "mA_upper_edge_probe",
        "sensitivity_probe",
        "control_box",
        "stop",
    }


def test_choice_maps_to_existing_worker_contract_shape() -> None:
    action = response_to_worker_action(_response())
    assert action["action"] == "run_template"
    assert action["template_name"] == "mphi_lower_edge_probe"
    assert action["template_overrides"] == {}
    assert "plancheck_required_before_execution" in action["safety_checks"]


def test_stop_choice_maps_to_worker_stop() -> None:
    action = response_to_worker_action(_response("stop"))
    assert action["action"] == "stop"
    assert "template_name" not in action


def test_non_allowlisted_choice_is_rejected() -> None:
    bad = _response()
    bad["answers"]["next_operator"]["choice"] = "edit_cpp"
    try:
        response_to_worker_action(bad)
        assert False, "expected JevError"
    except JevError as exc:
        assert "non-allowlisted" in str(exc)


def test_choose_with_jev_is_deterministic_adapter_around_api_response() -> None:
    with patch("scripts.jev_policy.call_jev", return_value=_response()) as mocked:
        d = choose_with_jev(_cfg(), _state(), api_key="not-used")
    assert d.selected == "mphi_lower_edge_probe"
    assert d.raw_response["usage"]["input_tokens"] == 900
    mocked.assert_called_once()


def test_training_event_preserves_state_probabilities_and_usage() -> None:
    with patch("scripts.jev_policy.call_jev", return_value=_response()):
        d = choose_with_jev(_cfg(), _state(), api_key="not-used")
    event = build_decision_event(
        d,
        timestamp="2026-09-24T20:00:00Z",
        iteration=2,
        iteration_id="iter02_test",
        cycle_id="cycle_test",
    )
    assert event["schema_version"] == "dihiggs_policy_decision_v1"
    assert event["backend"] == "jev_via_vercel_ai_gateway"
    assert event["answers"]["next_operator"]["probabilities"]["mphi_lower_edge_probe"] == 0.55
    assert event["usage"]["input_tokens"] == 900
    json.dumps(event)  # serialization is part of the logging contract


def _proc(stdout: str, returncode: int = 0) -> object:
    class _Proc:
        pass

    p = _Proc()
    p.returncode = returncode
    p.stdout = stdout
    p.stderr = ""
    return p


def test_choose_with_gemini_parses_json_reply_through_same_allowlist_path() -> None:
    reply = json.dumps(_response(choice="control_box"))
    with patch("scripts.jev_policy.subprocess.run", return_value=_proc(reply)) as mocked:
        d = choose_with_gemini(_cfg(), _state())
    assert d.selected == "control_box"
    assert d.raw_response["usage"]["input_tokens"] == 900
    mocked.assert_called_once()
    args = mocked.call_args.args[0]
    assert args[0] == "agy" and args[1] == "-p"


def test_choose_with_gemini_tolerates_markdown_fence_around_json() -> None:
    reply = "```json\n" + json.dumps(_response(choice="stop")) + "\n```"
    with patch("scripts.jev_policy.subprocess.run", return_value=_proc(reply)):
        d = choose_with_gemini(_cfg(), _state())
    assert d.selected == "stop"


def test_choose_with_gemini_rejects_non_allowlisted_choice() -> None:
    reply = json.dumps(_response(choice="edit_cpp"))
    with patch("scripts.jev_policy.subprocess.run", return_value=_proc(reply)):
        try:
            choose_with_gemini(_cfg(), _state())
            assert False, "expected JevError for non-allowlisted choice"
        except JevError as exc:
            assert "non-allowlisted" in str(exc)


def test_choose_with_gemini_raises_on_non_json_output() -> None:
    with patch("scripts.jev_policy.subprocess.run", return_value=_proc("not json at all")):
        try:
            choose_with_gemini(_cfg(), _state())
            assert False, "expected JevError for non-JSON gemini output"
        except JevError as exc:
            assert "non-JSON" in str(exc)


def test_choose_with_gemini_raises_on_nonzero_exit() -> None:
    with patch("scripts.jev_policy.subprocess.run", return_value=_proc("", returncode=1)):
        try:
            choose_with_gemini(_cfg(), _state())
            assert False, "expected JevError for nonzero gemini exit"
        except JevError as exc:
            assert "gemini CLI exited 1" in str(exc)
