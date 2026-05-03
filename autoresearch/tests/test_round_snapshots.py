from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import cast

from autoresearch.harness.round_snapshots import (
    classify_round,
    project_round_snapshot,
    write_campaign_status,
)
from autoresearch.harness.telemetry_store import init_db


def _insert_event(
    db: sqlite3.Connection,
    *,
    campaign_id: str,
    event_type: str,
    utc: str,
    payload: dict[str, object],
    dedupe_key: str,
) -> None:
    _ = db.execute(
        """
        INSERT INTO events (
            source_file, line_offset, checksum, dedupe_key, campaign_id, event_type, utc, payload, ingested_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            "events.jsonl",
            1,
            dedupe_key,
            dedupe_key,
            campaign_id,
            event_type,
            utc,
            json.dumps(payload),
            utc,
        ),
    )
    db.commit()


def _round_payload(
    round_index: int,
    *,
    arm_id: str = "adaptive-smoke",
    discoveries: int,
    events_emitted: int,
    subprocess_status: str = "success",
    metrics: dict[str, float] | None = None,
    utc: str,
) -> dict[str, object]:
    return {
        "round_index": round_index,
        "arm_id": arm_id,
        "discoveries": discoveries,
        "events_emitted": events_emitted,
        "subprocess_status": subprocess_status,
        "metrics": metrics or {
            "yield": 0.0,
            "coverage": 0.0,
            "diversity": 0.0,
            "composite": 0.0,
        },
        "utc": utc,
    }


def _as_float(value: object) -> float:
    assert isinstance(value, (int, float))
    return float(value)


def test_classify_round_productive() -> None:
    result = classify_round({"discoveries": 2, "events_emitted": 2, "subprocess_status": "success"})
    assert result == {
        "productive": True,
        "duplicate_only": False,
        "empty_round": False,
        "subprocess_status": "success",
    }


def test_classify_round_duplicate_only_empty() -> None:
    result = classify_round({"discoveries": 3, "events_emitted": 0, "subprocess_status": "success"})
    assert result == {
        "productive": False,
        "duplicate_only": True,
        "empty_round": True,
        "subprocess_status": "success",
    }


def test_project_round_snapshot_includes_failed_subprocess_and_last_activity(tmp_path: Path) -> None:
    db = init_db(tmp_path / "campaign_monitor.db")
    _insert_event(
        db,
        campaign_id="campaign-1",
        event_type="ATTEMPT_EVALUATED",
        utc="2026-04-06T12:00:01Z",
        payload={"attempt_id": "a-1"},
        dedupe_key="attempt-1",
    )
    _insert_event(
        db,
        campaign_id="campaign-1",
        event_type="ROUND_COMPLETED",
        utc="2026-04-06T12:00:02Z",
        payload=_round_payload(
            0,
            discoveries=0,
            events_emitted=0,
            subprocess_status="fail",
            utc="2026-04-06T12:00:02Z",
        ),
        dedupe_key="round-0",
    )

    snapshot = project_round_snapshot(db, 0, campaign_state="FAILED")

    assert snapshot["campaign_id"] == "campaign-1"
    assert snapshot["campaign_state"] == "FAILED"
    assert snapshot["selected_arm"] == "adaptive-smoke"
    assert snapshot["subprocess_status"] == "fail"
    assert snapshot["empty_round"] is True
    assert snapshot["productive"] is False
    assert snapshot["last_activity_utc"] == "2026-04-06T12:00:01Z"
    assert snapshot["rolling_metrics"] == {
        "yield": 0.0,
        "coverage": 0.0,
        "diversity": 0.0,
        "composite": 0.0,
    }
    assert snapshot["rolling_deltas"] == {
        "yield_delta": None,
        "coverage_delta": None,
        "diversity_delta": None,
        "composite_delta": None,
    }


def test_project_round_snapshot_computes_round_six_rolling_deltas(tmp_path: Path) -> None:
    db = init_db(tmp_path / "campaign_monitor.db")
    metric_rows = [
        {"yield": 0.10, "coverage": 0.20, "diversity": 0.30, "composite": 0.40},
        {"yield": 0.11, "coverage": 0.22, "diversity": 0.33, "composite": 0.44},
        {"yield": 0.12, "coverage": 0.24, "diversity": 0.36, "composite": 0.48},
        {"yield": 0.13, "coverage": 0.26, "diversity": 0.39, "composite": 0.52},
        {"yield": 0.14, "coverage": 0.28, "diversity": 0.42, "composite": 0.56},
        {"yield": 0.25, "coverage": 0.35, "diversity": 0.55, "composite": 0.75},
    ]

    for round_index, metrics in enumerate(metric_rows):
        _insert_event(
            db,
            campaign_id="campaign-rolling",
            event_type="ROUND_COMPLETED",
            utc=f"2026-04-06T12:00:0{round_index}Z",
            payload=_round_payload(
                round_index,
                discoveries=round_index + 1,
                events_emitted=1,
                metrics=metrics,
                utc=f"2026-04-06T12:00:0{round_index}Z",
            ),
            dedupe_key=f"round-{round_index}",
        )

    snapshot = project_round_snapshot(db, 5)

    assert snapshot["productive"] is True
    assert snapshot["duplicate_only"] is False
    assert snapshot["rolling_metrics"] == metric_rows[-1]
    deltas = snapshot["rolling_deltas"]
    assert isinstance(deltas, dict)
    delta_map = cast(dict[str, object], deltas)
    assert abs(_as_float(delta_map["yield_delta"]) - 0.15) < 1e-12
    assert abs(_as_float(delta_map["coverage_delta"]) - 0.15) < 1e-12
    assert abs(_as_float(delta_map["diversity_delta"]) - 0.25000000000000006) < 1e-12
    assert abs(_as_float(delta_map["composite_delta"]) - 0.35) < 1e-12


def test_write_campaign_status_writes_latest_snapshot_file(tmp_path: Path) -> None:
    snapshot: dict[str, object] = {
        "schema_version": 1,
        "campaign_id": "campaign-1",
        "campaign_state": "RUNNING",
        "round_index": 2,
        "selected_arm": "adaptive-smoke",
        "subprocess_status": "success",
        "discoveries": 2,
        "events_emitted": 2,
        "empty_round": False,
        "duplicate_only": False,
        "productive": True,
        "rolling_metrics": {
            "yield": 0.42,
            "coverage": 0.67,
            "diversity": 1.23,
            "composite": 0.77,
        },
        "rolling_deltas": {
            "yield_delta": 0.01,
            "coverage_delta": 0.02,
            "diversity_delta": 0.05,
            "composite_delta": 0.03,
        },
        "last_activity_utc": "2026-04-06T12:34:56Z",
        "active_alerts": [],
    }

    target = write_campaign_status(tmp_path, snapshot)

    assert target == tmp_path / "campaign_status.json"
    assert json.loads(target.read_text(encoding="utf-8")) == snapshot
