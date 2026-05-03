from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import cast

SCHEMA_VERSION = 1
ROUND_DELTA_WINDOW = 5
METRIC_NAMES = ("yield", "coverage", "diversity", "composite")


def _as_int(value: object, default: int = 0) -> int:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    if isinstance(value, str):
        try:
            return int(value)
        except ValueError:
            return default
    return default


def _as_str(value: object, default: str = "") -> str:
    return value if isinstance(value, str) else default


def _as_dict(value: object) -> dict[str, object]:
    return cast(dict[str, object], value) if isinstance(value, dict) else {}


def _as_list(value: object) -> list[object]:
    return cast(list[object], value) if isinstance(value, list) else []


def _as_float(value: object, default: float = 0.0) -> float:
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return default
    return default


def classify_round(round_data: dict[str, object]) -> dict[str, object]:
    events_emitted = _as_int(round_data.get("events_emitted", 0))
    discoveries = _as_int(round_data.get("discoveries", 0))
    empty_round = events_emitted == 0
    return {
        "productive": discoveries > 0 and events_emitted > 0,
        "duplicate_only": discoveries > 0 and empty_round,
        "empty_round": empty_round,
        "subprocess_status": _as_str(round_data.get("subprocess_status", "unknown"), "unknown"),
    }


def _load_events(db: sqlite3.Connection) -> list[dict[str, object]]:
    cursor = db.execute(
        "SELECT id, campaign_id, event_type, utc, payload FROM events ORDER BY utc ASC, id ASC"
    )
    events: list[dict[str, object]] = []
    rows = cast(list[tuple[int, str, str, str, str]], cursor.fetchall())
    for row_id, campaign_id, event_type, utc, payload_text in rows:
        payload_obj = cast(object, json.loads(payload_text))
        payload = _as_dict(payload_obj)
        events.append(
            {
                "id": row_id,
                "campaign_id": campaign_id,
                "event_type": event_type,
                "utc": utc,
                "payload": payload,
            }
        )
    return events


def _round_payload(event: dict[str, object]) -> dict[str, object] | None:
    payload = event["payload"]
    if not isinstance(payload, dict):
        return None
    payload_dict = cast(dict[str, object], payload)
    required = {"round_index", "discoveries", "events_emitted", "subprocess_status"}
    if not required.issubset(payload_dict.keys()):
        return None
    if "arm_id" not in payload_dict and "selected_arm" not in payload_dict:
        return None
    return payload_dict


def _normalize_metrics(payload: dict[str, object]) -> dict[str, float]:
    metrics = _as_dict(payload.get("metrics", {}))
    return {name: _as_float(metrics.get(name, 0.0)) for name in METRIC_NAMES}


def _compute_deltas(round_payloads: list[dict[str, object]], current_index: int) -> dict[str, float | None]:
    deltas: dict[str, float | None] = {}
    prior_index = current_index - ROUND_DELTA_WINDOW
    current_metrics = _normalize_metrics(round_payloads[current_index])
    if prior_index < 0:
        for name in METRIC_NAMES:
            deltas[f"{name}_delta"] = None
        return deltas
    prior_metrics = _normalize_metrics(round_payloads[prior_index])
    for name in METRIC_NAMES:
        deltas[f"{name}_delta"] = current_metrics[name] - prior_metrics[name]
    return deltas


def project_round_snapshot(
    db: sqlite3.Connection,
    round_index: int,
    *,
    campaign_state: str | None = None,
    active_alerts: list[object] | None = None,
) -> dict[str, object]:
    events = _load_events(db)
    round_payloads = [payload for event in events if (payload := _round_payload(event)) is not None]
    round_payloads.sort(key=lambda payload: _as_int(payload["round_index"]))

    current_position = next(
        (index for index, payload in enumerate(round_payloads) if _as_int(payload["round_index"]) == round_index),
        None,
    )
    if current_position is None:
        raise ValueError(f"No normalized telemetry found for round_index={round_index}")

    round_payload = round_payloads[current_position]
    round_utc = _as_str(round_payload.get("utc", ""))
    matching_events = [event for event in events if _as_str(event["utc"]) <= round_utc] if round_utc else events
    last_activity_utc = round_payload.get("last_activity_utc")
    if last_activity_utc is None:
        for event in reversed(matching_events):
            if event["event_type"] == "ATTEMPT_EVALUATED":
                last_activity_utc = event["utc"]
                break
    if last_activity_utc is None:
        last_activity_utc = round_payload.get("utc")

    selected_arm = _as_str(round_payload.get("selected_arm", round_payload.get("arm_id", "")))
    classification = classify_round(round_payload)
    campaign_id = next((event["campaign_id"] for event in events if event.get("campaign_id")), round_payload.get("campaign_id", ""))
    raw_alerts = active_alerts if active_alerts is not None else _as_list(round_payload.get("active_alerts", []))
    snapshot = {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": _as_str(campaign_id),
        "campaign_state": campaign_state or _as_str(round_payload.get("campaign_state", "RUNNING"), "RUNNING"),
        "round_index": _as_int(round_payload["round_index"]),
        "selected_arm": selected_arm,
        "subprocess_status": _as_str(classification["subprocess_status"], "unknown"),
        "discoveries": _as_int(round_payload.get("discoveries", 0)),
        "events_emitted": _as_int(round_payload.get("events_emitted", 0)),
        "empty_round": classification["empty_round"],
        "duplicate_only": classification["duplicate_only"],
        "productive": classification["productive"],
        "rolling_metrics": _normalize_metrics(round_payload),
        "rolling_deltas": _compute_deltas(round_payloads, current_position),
        "last_activity_utc": _as_str(last_activity_utc),
        "active_alerts": raw_alerts,
    }
    return snapshot


def write_campaign_status(outdir: str | Path, snapshot: dict[str, object]) -> Path:
    target = Path(outdir) / "campaign_status.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    _ = target.write_text(json.dumps(snapshot, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    return target
