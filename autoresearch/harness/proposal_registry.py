from __future__ import annotations

import copy
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from collections.abc import Iterable, Mapping

IMMUTABLE_FIELDS: tuple[str, ...] = (
    "schema_version",
    "proposal_id",
    "parent_id",
    "generation",
    "operator",
    "created_utc",
)

FROZEN_AFTER_PLANNED_FIELDS: tuple[str, ...] = (
    "bounds",
    "fixed",
    "resolution",
    "budget",
)

from autoresearch.harness.dihiggs_validators import validate_region

REQUIRED_FIELDS: tuple[str, ...] = (
    "schema_version",
    "proposal_id",
    "parent_id",
    "generation",
    "operator",
    "bounds",
    "fixed",
    "resolution",
    "budget",
    "status",
    "created_utc",
    "updated_utc",
    "run_dir",
    "metrics",
    "error",
)

SUPPORTED_STATUSES: tuple[str, ...] = (
    "planned",
    "validated",
    "submitted",
    "running",
    "done",
    "failed",
    "timeout",
    "skipped",
    "rejected",
)

_ALLOWED_TRANSITIONS: dict[str, set[str]] = {
    "planned": {"validated", "rejected"},
    "validated": {"submitted", "skipped"},
    "submitted": {"running"},
    "running": {"done", "failed", "timeout"},
    "failed": {"submitted"},
    "timeout": {"submitted"},
    "done": set(),
    "skipped": set(),
    "rejected": set(),
}


def _as_finite_float(value: object, field: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{field} must be numeric")
    if isinstance(value, (int, float, str)):
        try:
            out = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{field} must be numeric") from exc
        if not math.isfinite(out):
            raise ValueError(f"{field} must be finite")
        return out
    raise ValueError(f"{field} must be numeric")


def _as_int(value: object, field: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{field} must be integer")
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if not math.isfinite(value) or not value.is_integer():
            raise ValueError(f"{field} must be integer")
        return int(value)
    if isinstance(value, str):
        text = value.strip()
        if text == "":
            raise ValueError(f"{field} must be integer")
        try:
            parsed = float(text)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{field} must be integer") from exc
        if not math.isfinite(parsed) or not parsed.is_integer():
            raise ValueError(f"{field} must be integer")
        return int(parsed)
    raise ValueError(f"{field} must be integer")


def _parse_utc_timestamp(value: object, field: str) -> datetime:
    if not isinstance(value, str) or value.strip() == "":
        raise ValueError(f"{field} must be a non-empty string")
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError(f"{field} must be an ISO-8601 timestamp") from exc


def _normalize_for_compare(value: object) -> object:
    if isinstance(value, Mapping):
        return {str(k): _normalize_for_compare(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_normalize_for_compare(v) for v in value]
    return value


def _values_semantically_equal(left: object, right: object) -> bool:
    return _normalize_for_compare(left) == _normalize_for_compare(right)


def _validate_fixed_field_against_contract(fixed: Mapping[str, object], contract_fixed: Mapping[str, object], errors: list[str]) -> None:
    for field in ("sin_ba", "mA", "mHp", "lambda7"):
        if field in fixed and field in contract_fixed:
            try:
                lhs = _as_finite_float(fixed[field], f"fixed.{field}")
            except ValueError as exc:
                errors.append(str(exc))
                continue
            try:
                rhs = _as_finite_float(contract_fixed[field], f"contract.fixed.{field}")
            except ValueError:
                errors.append(f"contract.fixed.{field} must be numeric")
                continue
            if lhs != rhs:
                errors.append(f"fixed.{field} must match contract fixed value")

    field = "yukawa_type"
    if field in fixed and field in contract_fixed:
        try:
            lhs_i = _as_int(fixed[field], f"fixed.{field}")
        except ValueError as exc:
            errors.append(str(exc))
            lhs_i = None
        try:
            rhs_i = _as_int(contract_fixed[field], f"contract.fixed.{field}")
        except ValueError:
            errors.append(f"contract.fixed.{field} must be integer")
            rhs_i = None
        if lhs_i is not None and rhs_i is not None and lhs_i != rhs_i:
            errors.append(f"fixed.{field} must match contract fixed value")


def _normalize_bounds(raw_bounds: object) -> dict[str, object]:
    if not isinstance(raw_bounds, Mapping):
        raise ValueError("bounds must be a mapping")
    out: dict[str, object] = {}
    for key, raw_value in raw_bounds.items():
        if isinstance(raw_value, (list, tuple)) and len(raw_value) == 2:
            lo = _as_finite_float(raw_value[0], f"bounds.{key}[0]")
            hi = _as_finite_float(raw_value[1], f"bounds.{key}[1]")
            out[str(key)] = [lo, hi]
        else:
            out[str(key)] = raw_value
    return out


def make_proposal_id(prefix: str = "ar", *, timestamp: str | None = None, index: int | None = None) -> str:
    ts = timestamp if timestamp is not None else datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    idx = 1 if index is None else index
    if idx < 0:
        raise ValueError("index must be >= 0")
    return f"{prefix}_{ts}_{idx:06d}"


def normalize_proposal(raw: Mapping[str, object], contract: Mapping[str, object] | None = None) -> dict[str, object]:
    out: dict[str, object] = copy.deepcopy(dict(raw))

    if "bounds" in out:
        out["bounds"] = _normalize_bounds(out["bounds"])

    if "generation" in out:
        out["generation"] = _as_int(out["generation"], "generation")

    resolution = out.get("resolution")
    if isinstance(resolution, Mapping):
        res_out: dict[str, object] = dict(resolution)
        for key, value in list(res_out.items()):
            res_out[str(key)] = _as_int(value, f"resolution.{key}")
        out["resolution"] = res_out

    budget = out.get("budget")
    if isinstance(budget, Mapping):
        budget_out: dict[str, object] = dict(budget)
        if "requested_points" in budget_out:
            budget_out["requested_points"] = _as_int(budget_out["requested_points"], "budget.requested_points")
        if "max_walltime_sec" in budget_out:
            budget_out["max_walltime_sec"] = _as_finite_float(
                budget_out["max_walltime_sec"], "budget.max_walltime_sec"
            )
        out["budget"] = budget_out

    parent_id = out.get("parent_id")
    if parent_id is None:
        out["parent_id"] = None
    elif isinstance(parent_id, str):
        text = parent_id.strip()
        out["parent_id"] = text if text != "" else None
    else:
        raise ValueError("parent_id must be None or non-empty string")

    metrics = out.get("metrics")
    if metrics is None:
        out["metrics"] = {}
    elif isinstance(metrics, Mapping):
        out["metrics"] = dict(metrics)
    else:
        raise ValueError("metrics must be a mapping")

    if "error" in out:
        err = out["error"]
        if err is None or isinstance(err, str):
            out["error"] = err
        else:
            raise ValueError("error must be None or string")
    else:
        out["error"] = None

    if contract is not None and "bounds" in out:
        region_errors = validate_region(out["bounds"] if isinstance(out["bounds"], Mapping) else {}, contract)
        if region_errors:
            raise ValueError("invalid bounds for contract: " + "; ".join(region_errors))

    return out


def validate_proposal(proposal: Mapping[str, object], contract: Mapping[str, object] | None = None) -> list[str]:
    errors: list[str] = []

    for field in REQUIRED_FIELDS:
        if field not in proposal:
            errors.append(f"missing required field: {field}")

    proposal_id = proposal.get("proposal_id")
    if not isinstance(proposal_id, str) or proposal_id.strip() == "":
        errors.append("proposal_id must be a non-empty string")

    status = proposal.get("status")
    if not isinstance(status, str) or status not in SUPPORTED_STATUSES:
        errors.append("status must be one of supported statuses")

    generation = proposal.get("generation")
    try:
        gen_i = _as_int(generation, "generation")
        if gen_i < 0:
            errors.append("generation must be >= 0")
    except ValueError as exc:
        errors.append(str(exc))

    operator = proposal.get("operator")
    if not isinstance(operator, str) or operator.strip() == "":
        errors.append("operator must be a non-empty string")

    parent_id = proposal.get("parent_id")
    if parent_id is not None and (not isinstance(parent_id, str) or parent_id.strip() == ""):
        errors.append("parent_id must be None or non-empty string")

    bounds = proposal.get("bounds")
    if not isinstance(bounds, Mapping):
        errors.append("bounds must be a mapping")
    elif contract is not None:
        region_errors = validate_region(bounds, contract)
        errors.extend(region_errors)

    fixed = proposal.get("fixed")
    if not isinstance(fixed, Mapping):
        errors.append("fixed must be a mapping")
    elif contract is not None:
        contract_fixed = contract.get("fixed")
        if isinstance(contract_fixed, Mapping):
            _validate_fixed_field_against_contract(fixed, contract_fixed, errors)

    resolution = proposal.get("resolution")
    if not isinstance(resolution, Mapping):
        errors.append("resolution must be a mapping")
    else:
        for key, value in resolution.items():
            try:
                parsed = _as_int(value, f"resolution.{key}")
                if parsed <= 0:
                    errors.append(f"resolution.{key} must be positive")
            except ValueError as exc:
                errors.append(str(exc))

    budget = proposal.get("budget")
    if not isinstance(budget, Mapping):
        errors.append("budget must be a mapping")
    else:
        if "requested_points" not in budget:
            errors.append("budget.requested_points is required")
        else:
            try:
                rp = _as_int(budget["requested_points"], "budget.requested_points")
                if rp <= 0:
                    errors.append("budget.requested_points must be positive")
            except ValueError as exc:
                errors.append(str(exc))

        if "max_walltime_sec" in budget:
            try:
                wall = _as_finite_float(budget["max_walltime_sec"], "budget.max_walltime_sec")
                if wall <= 0:
                    errors.append("budget.max_walltime_sec must be positive")
            except ValueError as exc:
                errors.append(str(exc))

    metrics = proposal.get("metrics")
    if metrics is not None and not isinstance(metrics, Mapping):
        errors.append("metrics must be a mapping")

    run_dir = proposal.get("run_dir")
    if run_dir is not None and (not isinstance(run_dir, str) or run_dir.strip() == ""):
        errors.append("run_dir must be None or non-empty string")

    err = proposal.get("error")
    if err is not None and not isinstance(err, str):
        errors.append("error must be None or string")

    for ts_field in ("created_utc", "updated_utc"):
        try:
            _parse_utc_timestamp(proposal.get(ts_field), ts_field)
        except ValueError as exc:
            errors.append(str(exc))

    return errors


def append_event(path: str | Path, event: Mapping[str, object]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps(dict(event), sort_keys=True, separators=(",", ":")) + "\n"
    with target.open("a", encoding="utf-8") as handle:
        handle.write(line)
        handle.flush()
        os.fsync(handle.fileno())


def load_events(path: str | Path) -> list[dict[str, object]]:
    target = Path(path)
    if not target.exists():
        return []

    events: list[dict[str, object]] = []
    with target.open("r", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.strip()
            if line == "":
                continue
            try:
                parsed = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"invalid JSON at line {line_no}") from exc
            if not isinstance(parsed, dict):
                raise ValueError(f"non-dict JSON at line {line_no}")
            events.append(parsed)
    return events


def replay_registry(
    events: Iterable[Mapping[str, object]], contract: Mapping[str, object] | None = None
) -> dict[str, dict[str, object]]:
    state: dict[str, dict[str, object]] = {}

    for index, event in enumerate(events, start=1):
        event_type = event.get("event_type")
        if event_type == "create_proposal":
            proposal_raw = event.get("proposal")
            if not isinstance(proposal_raw, Mapping):
                raise ValueError(f"event {index}: create_proposal requires proposal mapping")
            proposal = normalize_proposal(proposal_raw, contract=None)
            proposal_id = proposal.get("proposal_id")
            if not isinstance(proposal_id, str) or proposal_id.strip() == "":
                raise ValueError(f"event {index}: proposal.proposal_id must be non-empty string")
            create_key = event.get("proposal_id")
            if create_key is not None:
                if not isinstance(create_key, str) or create_key.strip() == "":
                    raise ValueError(f"event {index}: create_proposal proposal_id must be non-empty string when provided")
                if create_key != proposal_id:
                    raise ValueError(
                        f"event {index}: state key/payload proposal_id mismatch for {create_key} vs {proposal_id}"
                    )
                proposal_id = create_key
            if proposal_id in state:
                raise ValueError(f"event {index}: duplicate create_proposal for {proposal_id}")
            if proposal_id != proposal.get("proposal_id"):
                raise ValueError(f"event {index}: state key/payload proposal_id mismatch for {proposal_id}")
            state[proposal_id] = proposal
        elif event_type == "update_proposal":
            proposal_id = event.get("proposal_id")
            if not isinstance(proposal_id, str) or proposal_id.strip() == "":
                raise ValueError(f"event {index}: update_proposal requires non-empty proposal_id")
            if proposal_id not in state:
                raise ValueError(f"event {index}: update for missing proposal_id {proposal_id}")
            updates = event.get("updates")
            if not isinstance(updates, Mapping):
                raise ValueError(f"event {index}: update_proposal requires updates mapping")

            current = state[proposal_id]
            for field in IMMUTABLE_FIELDS:
                if field in updates and not _values_semantically_equal(updates[field], current.get(field)):
                    raise ValueError(f"event {index}: proposal_id {proposal_id} attempted immutable field update: {field}")

            current_status_obj = current.get("status")
            if "status" in updates:
                new_status_obj = updates.get("status")
                if not isinstance(current_status_obj, str) or not isinstance(new_status_obj, str):
                    raise ValueError(f"event {index}: proposal_id {proposal_id} status must be string")
                if not is_valid_status_transition(current_status_obj, new_status_obj):
                    raise ValueError(
                        f"event {index}: proposal_id {proposal_id} invalid status transition "
                        f"{current_status_obj} -> {new_status_obj}"
                    )

            if isinstance(current_status_obj, str) and current_status_obj != "planned":
                for field in FROZEN_AFTER_PLANNED_FIELDS:
                    if field in updates and not _values_semantically_equal(updates[field], current.get(field)):
                        raise ValueError(
                            f"event {index}: proposal_id {proposal_id} cannot update {field} after status {current_status_obj}"
                        )

            merged = copy.deepcopy(current)
            for key, value in updates.items():
                merged[str(key)] = copy.deepcopy(value)
            normalized = normalize_proposal(merged, contract=None)
            payload_id = normalized.get("proposal_id")
            if payload_id != proposal_id:
                raise ValueError(
                    f"event {index}: state key/payload proposal_id mismatch for {proposal_id} vs {payload_id}"
                )
            state[proposal_id] = normalized
        else:
            raise ValueError(f"event {index}: unsupported event_type {event_type!r}")

    for proposal_id, proposal in state.items():
        proposal_errors = validate_proposal(proposal, contract=contract)
        if proposal_errors:
            raise ValueError(f"invalid final proposal {proposal_id}: " + "; ".join(proposal_errors))

    return state


def check_parent_links(state: Mapping[str, Mapping[str, object]]) -> list[str]:
    errors: list[str] = []

    for proposal_id, proposal in state.items():
        if not isinstance(proposal, Mapping):
            errors.append(f"proposal {proposal_id} has invalid payload")
            continue

        parent_id = proposal.get("parent_id")
        if parent_id is None:
            continue
        if not isinstance(parent_id, str) or parent_id.strip() == "":
            errors.append(f"proposal {proposal_id} has invalid parent_id")
            continue
        if parent_id == proposal_id:
            errors.append(f"proposal {proposal_id} has self-parent")
            continue
        if parent_id not in state:
            errors.append(f"proposal {proposal_id} references missing parent_id {parent_id}")
            continue

        child_gen = proposal.get("generation")
        parent_gen = state[parent_id].get("generation") if isinstance(state[parent_id], Mapping) else None
        try:
            child_gen_i = _as_int(child_gen, f"proposal {proposal_id} generation")
            parent_gen_i = _as_int(parent_gen, f"proposal {parent_id} generation")
            if child_gen_i <= parent_gen_i:
                errors.append(
                    f"proposal {proposal_id} generation {child_gen_i} must be > parent {parent_id} generation {parent_gen_i}"
                )
        except ValueError as exc:
            errors.append(str(exc))

    visiting: set[str] = set()
    visited: set[str] = set()

    def _dfs(node: str, path: list[str]) -> None:
        if node in visited:
            return
        if node in visiting:
            cycle_start = path.index(node) if node in path else 0
            cycle_path = path[cycle_start:] + [node]
            errors.append("parent cycle detected: " + " -> ".join(cycle_path))
            return

        visiting.add(node)
        proposal = state.get(node)
        if isinstance(proposal, Mapping):
            parent_id = proposal.get("parent_id")
            if isinstance(parent_id, str) and parent_id in state:
                _dfs(parent_id, path + [parent_id])
        visiting.remove(node)
        visited.add(node)

    for proposal_id in state.keys():
        _dfs(proposal_id, [proposal_id])

    return errors


def is_valid_status_transition(old: str, new: str) -> bool:
    if old not in SUPPORTED_STATUSES or new not in SUPPORTED_STATUSES:
        return False
    if old == new:
        return True
    return new in _ALLOWED_TRANSITIONS.get(old, set())


def apply_status_update(proposal: Mapping[str, object], new_status: str, **extra: object) -> dict[str, object]:
    old_status_obj = proposal.get("status")
    if not isinstance(old_status_obj, str):
        raise ValueError("proposal.status must be a string")
    if new_status not in SUPPORTED_STATUSES:
        raise ValueError(f"unsupported status: {new_status}")
    if not is_valid_status_transition(old_status_obj, new_status):
        raise ValueError(f"invalid status transition: {old_status_obj} -> {new_status}")

    out = copy.deepcopy(dict(proposal))
    out["status"] = new_status
    for key, value in extra.items():
        out[str(key)] = copy.deepcopy(value)
    return out
