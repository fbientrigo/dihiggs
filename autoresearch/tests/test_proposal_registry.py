from __future__ import annotations

import copy

import pytest

from autoresearch.harness.dihiggs_physics_contract import get_default_physics_contract
from autoresearch.harness.proposal_registry import (
    append_event,
    apply_status_update,
    check_parent_links,
    is_valid_status_transition,
    load_events,
    make_proposal_id,
    normalize_proposal,
    replay_registry,
    validate_proposal,
)


def _contract() -> dict[str, object]:
    return get_default_physics_contract()


def _valid_proposal(proposal_id: str = "ar_20260503T220000Z_000001", parent_id: str | None = None) -> dict[str, object]:
    return {
        "schema_version": "1.0",
        "proposal_id": proposal_id,
        "parent_id": parent_id,
        "generation": 0,
        "operator": "seed",
        "bounds": {
            "lambda6": [1e-6, 0.5],
            "tan_beta": [2e3, 1e7],
            "lambda1": [0.2, 11.0],
            "mH": [130.0, 250.0],
        },
        "fixed": {
            "mh": 125.0,
            "sin_ba": 1.0,
            "mA": 300.0,
            "mHp": 300.0,
            "yukawa_type": 1,
            "lambda7": 0.0,
        },
        "resolution": {"lambda6": 8, "tan_beta": 12, "lambda1": 10, "mH": 16},
        "budget": {"requested_points": 500, "max_walltime_sec": 60.0},
        "status": "planned",
        "created_utc": "2026-05-03T22:00:00Z",
        "updated_utc": "2026-05-03T22:00:00Z",
        "run_dir": None,
        "metrics": {},
        "error": None,
    }


def test_make_proposal_id_deterministic_with_timestamp_and_index() -> None:
    assert make_proposal_id(timestamp="20260503T220000Z", index=1) == "ar_20260503T220000Z_000001"


def test_make_proposal_id_prefix_and_zero_padded_index() -> None:
    assert make_proposal_id(prefix="p2", timestamp="20260503T220000Z", index=42) == "p2_20260503T220000Z_000042"


def test_normalize_proposal_does_not_mutate_input() -> None:
    raw = _valid_proposal()
    raw["bounds"] = dict(raw["bounds"])
    raw["bounds"]["lambda6"] = ["1e-6", "0.5"]
    raw_copy = copy.deepcopy(raw)

    out = normalize_proposal(raw)
    assert raw == raw_copy
    assert out["bounds"]["lambda6"] == [1e-06, 0.5]


def test_normalize_proposal_normalizes_numeric_strings() -> None:
    raw = _valid_proposal()
    raw["generation"] = "1"
    raw["resolution"] = {"lambda6": "8", "tan_beta": "12"}
    raw["budget"] = {"requested_points": "500", "max_walltime_sec": "60.5"}
    raw["bounds"] = {
        "lambda6": ["1e-6", "0.5"],
        "tan_beta": ["2000", "1e7"],
        "lambda1": ["0.2", "11.0"],
        "mH": ["130", "250"],
    }

    out = normalize_proposal(raw)
    assert out["generation"] == 1
    assert out["resolution"]["lambda6"] == 8
    assert out["budget"]["requested_points"] == 500
    assert out["budget"]["max_walltime_sec"] == 60.5
    assert out["bounds"]["mH"] == [130.0, 250.0]


def test_normalize_proposal_invalid_region_fails_with_contract() -> None:
    raw = _valid_proposal()
    raw["bounds"] = dict(raw["bounds"])
    raw["bounds"]["mH"] = [125.0, 250.0]
    with pytest.raises(ValueError):
        normalize_proposal(raw, contract=_contract())


def test_validate_proposal_valid_minimal_passes() -> None:
    proposal = _valid_proposal()
    assert validate_proposal(proposal, contract=_contract()) == []


def test_validate_proposal_missing_required_fields_fail() -> None:
    proposal = _valid_proposal()
    del proposal["budget"]
    errors = validate_proposal(proposal)
    assert any("budget" in msg for msg in errors)


def test_validate_proposal_invalid_status_fails() -> None:
    proposal = _valid_proposal()
    proposal["status"] = "queued"
    errors = validate_proposal(proposal)
    assert any("status" in msg for msg in errors)


def test_validate_proposal_invalid_generation_fails() -> None:
    proposal = _valid_proposal()
    proposal["generation"] = -1
    errors = validate_proposal(proposal)
    assert any("generation" in msg for msg in errors)


def test_validate_proposal_empty_operator_fails() -> None:
    proposal = _valid_proposal()
    proposal["operator"] = ""
    errors = validate_proposal(proposal)
    assert any("operator" in msg for msg in errors)


def test_validate_proposal_invalid_bounds_fail_against_contract() -> None:
    proposal = _valid_proposal()
    proposal["bounds"] = dict(proposal["bounds"])
    proposal["bounds"]["lambda6"] = [1e-9, 0.5]
    errors = validate_proposal(proposal, contract=_contract())
    assert any("lambda6" in msg for msg in errors)


def test_validate_proposal_non_positive_resolution_fails() -> None:
    proposal = _valid_proposal()
    proposal["resolution"] = {"lambda6": 0}
    errors = validate_proposal(proposal)
    assert any("resolution.lambda6" in msg for msg in errors)


def test_validate_proposal_non_positive_requested_points_fails() -> None:
    proposal = _valid_proposal()
    proposal["budget"] = {"requested_points": 0}
    errors = validate_proposal(proposal)
    assert any("budget.requested_points" in msg for msg in errors)


def test_validate_proposal_non_positive_max_walltime_sec_fails() -> None:
    proposal = _valid_proposal()
    proposal["budget"] = {"requested_points": 1, "max_walltime_sec": 0}
    errors = validate_proposal(proposal)
    assert any("budget.max_walltime_sec" in msg for msg in errors)


def test_validate_proposal_error_messages_include_field_names() -> None:
    proposal = _valid_proposal()
    proposal["proposal_id"] = ""
    proposal["operator"] = ""
    errors = validate_proposal(proposal)
    joined = " | ".join(errors)
    assert "proposal_id" in joined
    assert "operator" in joined


def test_load_events_missing_file_returns_empty(tmp_path) -> None:
    path = tmp_path / "registry" / "events.jsonl"
    assert load_events(path) == []


def test_append_event_then_load_events_roundtrip(tmp_path) -> None:
    path = tmp_path / "registry" / "events.jsonl"
    event = {"event_type": "create_proposal", "proposal": _valid_proposal()}
    append_event(path, event)
    loaded = load_events(path)
    assert len(loaded) == 1
    assert loaded[0]["event_type"] == "create_proposal"


def test_load_events_ignores_blank_lines(tmp_path) -> None:
    path = tmp_path / "events.jsonl"
    path.write_text("\n\n{\"event_type\":\"noop\"}\n\n", encoding="utf-8")
    loaded = load_events(path)
    assert len(loaded) == 1


def test_load_events_invalid_json_raises_with_line_number(tmp_path) -> None:
    path = tmp_path / "events.jsonl"
    path.write_text("{\"event_type\":\"ok\"}\nnot-json\n", encoding="utf-8")
    with pytest.raises(ValueError, match="line 2"):
        load_events(path)


def test_load_events_non_dict_json_raises_with_line_number(tmp_path) -> None:
    path = tmp_path / "events.jsonl"
    path.write_text("[]\n", encoding="utf-8")
    with pytest.raises(ValueError, match="line 1"):
        load_events(path)


def test_replay_create_proposal_reconstructs_state() -> None:
    proposal = _valid_proposal()
    state = replay_registry([{"event_type": "create_proposal", "proposal": proposal}], contract=_contract())
    assert proposal["proposal_id"] in state


def test_replay_create_plus_update_modifies_state() -> None:
    proposal = _valid_proposal()
    events = [
        {"event_type": "create_proposal", "proposal": proposal},
        {
            "event_type": "update_proposal",
            "proposal_id": proposal["proposal_id"],
            "updates": {"status": "validated", "updated_utc": "2026-05-03T22:01:00Z"},
        },
    ]
    state = replay_registry(events, contract=_contract())
    assert state[proposal["proposal_id"]]["status"] == "validated"


def test_replay_duplicate_create_proposal_rejected() -> None:
    proposal = _valid_proposal()
    events = [
        {"event_type": "create_proposal", "proposal": proposal},
        {"event_type": "create_proposal", "proposal": proposal},
    ]
    with pytest.raises(ValueError, match="duplicate create_proposal"):
        replay_registry(events)


def test_replay_update_missing_proposal_id_rejected() -> None:
    proposal = _valid_proposal()
    events = [
        {"event_type": "create_proposal", "proposal": proposal},
        {"event_type": "update_proposal", "proposal_id": "missing", "updates": {"status": "validated"}},
    ]
    with pytest.raises(ValueError, match="missing proposal_id"):
        replay_registry(events)


def test_replay_final_proposal_validation_is_enforced() -> None:
    proposal = _valid_proposal()
    events = [
        {"event_type": "create_proposal", "proposal": proposal},
        {
            "event_type": "update_proposal",
            "proposal_id": proposal["proposal_id"],
            "updates": {"status": "bad-status"},
        },
    ]
    with pytest.raises(ValueError, match="invalid status transition"):
        replay_registry(events, contract=_contract())


def test_check_parent_links_reports_missing_parent() -> None:
    child = _valid_proposal(proposal_id="ar_20260503T220000Z_000002", parent_id="ar_missing")
    state = {child["proposal_id"]: child}
    errors = check_parent_links(state)
    assert any("missing parent_id" in msg for msg in errors)


def test_check_parent_links_valid_parent_passes() -> None:
    root = _valid_proposal("ar_20260503T220000Z_000001", None)
    child = _valid_proposal("ar_20260503T220000Z_000002", root["proposal_id"])
    child["generation"] = 1
    state = {root["proposal_id"]: root, child["proposal_id"]: child}
    assert check_parent_links(state) == []


def test_replay_event_order_is_preserved() -> None:
    proposal = _valid_proposal()
    events = [
        {"event_type": "create_proposal", "proposal": proposal},
        {
            "event_type": "update_proposal",
            "proposal_id": proposal["proposal_id"],
            "updates": {"status": "validated"},
        },
        {
            "event_type": "update_proposal",
            "proposal_id": proposal["proposal_id"],
            "updates": {"status": "submitted"},
        },
    ]
    state = replay_registry(events, contract=_contract())
    assert state[proposal["proposal_id"]]["status"] == "submitted"


def test_valid_forward_transitions_accepted() -> None:
    assert is_valid_status_transition("planned", "validated")
    assert is_valid_status_transition("running", "done")
    assert is_valid_status_transition("failed", "submitted")


def test_same_status_transition_accepted() -> None:
    assert is_valid_status_transition("running", "running")


def test_invalid_backward_transitions_rejected() -> None:
    assert not is_valid_status_transition("running", "planned")
    assert not is_valid_status_transition("validated", "planned")


def test_apply_status_update_does_not_mutate_input() -> None:
    proposal = _valid_proposal()
    original = copy.deepcopy(proposal)
    updated = apply_status_update(proposal, "validated")
    assert proposal == original
    assert updated["status"] == "validated"


def test_apply_status_update_can_attach_run_dir_metrics_or_error() -> None:
    proposal = _valid_proposal()
    validated = apply_status_update(proposal, "validated")
    submitted = apply_status_update(validated, "submitted", run_dir="/tmp/run", metrics={"x": 1}, error="")
    assert submitted["run_dir"] == "/tmp/run"
    assert submitted["metrics"] == {"x": 1}
    assert submitted["error"] == ""


def test_make_proposal_id_default_index_is_1_when_timestamp_reused() -> None:
    # Default index is intentionally stable but not globally unique for repeated calls.
    assert make_proposal_id(timestamp="20260503T220000Z") == "ar_20260503T220000Z_000001"


def test_replay_rejects_planned_to_done_direct_update() -> None:
    proposal = _valid_proposal()
    events = [
        {"event_type": "create_proposal", "proposal": proposal},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {"status": "done"}},
    ]
    with pytest.raises(ValueError, match="invalid status transition"):
        replay_registry(events)


def test_replay_rejects_done_to_planned_update() -> None:
    proposal = _valid_proposal()
    events = [
        {"event_type": "create_proposal", "proposal": proposal},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {"status": "validated"}},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {"status": "submitted"}},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {"status": "running"}},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {"status": "done"}},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {"status": "planned"}},
    ]
    with pytest.raises(ValueError, match="invalid status transition"):
        replay_registry(events)


def test_replay_accepts_full_forward_status_progression() -> None:
    proposal = _valid_proposal()
    events = [
        {"event_type": "create_proposal", "proposal": proposal},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {"status": "validated"}},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {"status": "submitted"}},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {"status": "running"}},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {"status": "done"}},
    ]
    state = replay_registry(events)
    assert state[proposal["proposal_id"]]["status"] == "done"


@pytest.mark.parametrize(
    "field,new_value",
    [("proposal_id", "other"), ("created_utc", "2026-01-01T00:00:00Z"), ("parent_id", "p2"), ("generation", 3)],
)
def test_replay_rejects_immutable_field_updates(field: str, new_value: object) -> None:
    proposal = _valid_proposal()
    events = [
        {"event_type": "create_proposal", "proposal": proposal},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {field: new_value}},
    ]
    with pytest.raises(ValueError, match="immutable field update"):
        replay_registry(events)


def test_replay_allows_bounds_update_while_planned() -> None:
    proposal = _valid_proposal()
    updated_bounds = dict(proposal["bounds"])
    updated_bounds["mH"] = [135.0, 240.0]
    events = [
        {"event_type": "create_proposal", "proposal": proposal},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {"bounds": updated_bounds}},
    ]
    state = replay_registry(events, contract=_contract())
    assert state[proposal["proposal_id"]]["bounds"]["mH"] == [135.0, 240.0]


def test_replay_rejects_bounds_update_after_validated() -> None:
    proposal = _valid_proposal()
    events = [
        {"event_type": "create_proposal", "proposal": proposal},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {"status": "validated"}},
        {
            "event_type": "update_proposal",
            "proposal_id": proposal["proposal_id"],
            "updates": {"bounds": {**proposal["bounds"], "mH": [140.0, 220.0]}},
        },
    ]
    with pytest.raises(ValueError, match="cannot update bounds"):
        replay_registry(events)


def test_replay_rejects_budget_update_after_validated() -> None:
    proposal = _valid_proposal()
    events = [
        {"event_type": "create_proposal", "proposal": proposal},
        {"event_type": "update_proposal", "proposal_id": proposal["proposal_id"], "updates": {"status": "validated"}},
        {
            "event_type": "update_proposal",
            "proposal_id": proposal["proposal_id"],
            "updates": {"budget": {"requested_points": 999, "max_walltime_sec": 120.0}},
        },
    ]
    with pytest.raises(ValueError, match="cannot update budget"):
        replay_registry(events)


def test_replay_rejects_create_key_payload_mismatch() -> None:
    proposal = _valid_proposal()
    events = [
        {
            "event_type": "create_proposal",
            "proposal_id": "other_id",
            "proposal": proposal,
        }
    ]
    with pytest.raises(ValueError, match="state key/payload proposal_id mismatch"):
        replay_registry(events)


def test_replay_rejects_update_key_payload_mismatch() -> None:
    proposal = _valid_proposal()
    events = [
        {"event_type": "create_proposal", "proposal": proposal},
        {
            "event_type": "update_proposal",
            "proposal_id": proposal["proposal_id"],
            "updates": {"proposal_id": "other_id"},
        },
    ]
    with pytest.raises(ValueError, match="immutable field update"):
        replay_registry(events)


def test_validate_proposal_run_dir_type_checks() -> None:
    for bad in ({"x": 1}, ["x"], 3):
        proposal = _valid_proposal()
        proposal["run_dir"] = bad
        errors = validate_proposal(proposal)
        assert any("run_dir" in msg for msg in errors)


def test_validate_proposal_run_dir_accepts_none_or_non_empty_string() -> None:
    p1 = _valid_proposal()
    p1["run_dir"] = None
    assert not any("run_dir" in msg for msg in validate_proposal(p1))
    p2 = _valid_proposal()
    p2["run_dir"] = "/tmp/run"
    assert not any("run_dir" in msg for msg in validate_proposal(p2))


def test_validate_proposal_rejects_empty_created_utc() -> None:
    proposal = _valid_proposal()
    proposal["created_utc"] = ""
    errors = validate_proposal(proposal)
    assert any("created_utc" in msg for msg in errors)


def test_validate_proposal_rejects_empty_updated_utc() -> None:
    proposal = _valid_proposal()
    proposal["updated_utc"] = ""
    errors = validate_proposal(proposal)
    assert any("updated_utc" in msg for msg in errors)


def test_check_parent_links_rejects_self_parent() -> None:
    proposal = _valid_proposal("ar_20260503T220000Z_000010", "ar_20260503T220000Z_000010")
    errors = check_parent_links({proposal["proposal_id"]: proposal})
    assert any("self-parent" in msg for msg in errors)


def test_check_parent_links_rejects_parent_cycle() -> None:
    a = _valid_proposal("a", "b")
    a["generation"] = 1
    b = _valid_proposal("b", "a")
    b["generation"] = 2
    errors = check_parent_links({"a": a, "b": b})
    assert any("parent cycle detected" in msg for msg in errors)


def test_check_parent_links_rejects_child_generation_not_greater_than_parent() -> None:
    root = _valid_proposal("root", None)
    root["generation"] = 1
    child = _valid_proposal("child", "root")
    child["generation"] = 1
    errors = check_parent_links({"root": root, "child": child})
    assert any("must be > parent" in msg for msg in errors)


def test_validate_proposal_accepts_fixed_numeric_strings_matching_contract() -> None:
    proposal = _valid_proposal()
    proposal["fixed"] = dict(proposal["fixed"])
    proposal["fixed"]["sin_ba"] = "1.0"
    proposal["fixed"]["mA"] = "300.0"
    proposal["fixed"]["mHp"] = "300"
    proposal["fixed"]["lambda7"] = "0.0"
    proposal["fixed"]["yukawa_type"] = "1"
    errors = validate_proposal(proposal, contract=_contract())
    assert not any("fixed." in msg for msg in errors)


def test_validate_proposal_rejects_non_numeric_fixed_mismatch_with_field_name() -> None:
    proposal = _valid_proposal()
    proposal["fixed"] = dict(proposal["fixed"])
    proposal["fixed"]["mA"] = "three hundred"
    errors = validate_proposal(proposal, contract=_contract())
    assert any("fixed.mA" in msg for msg in errors)
