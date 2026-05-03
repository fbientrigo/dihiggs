from __future__ import annotations
# pyright: reportArgumentType=false, reportIndexIssue=false, reportGeneralTypeIssues=false, reportMissingTypeArgument=false

import json
import tempfile
from pathlib import Path

import pytest

from ..harness.alerting import ALERT_SEVERITY_MAP, Alert, AlertEngine


def _minimal_state() -> dict[str, object]:
    return {
        "campaign_state": "RUNNING",
        "preflight_overall": "pass",
        "active_alerts": [],
    }


def test_alert_fingerprint_stable_for_same_type_and_context() -> None:
    alert1 = Alert(
        alert_type="preflight_blocked",
        severity="critical",
        message="Test",
        context={"key": "value"},
        first_seen="2026-04-06T00:00:00Z",
        last_seen="2026-04-06T00:00:00Z",
    )
    alert2 = Alert(
        alert_type="preflight_blocked",
        severity="critical",
        message="Different message",
        context={"key": "value"},
        first_seen="2026-04-06T01:00:00Z",
        last_seen="2026-04-06T01:00:00Z",
    )
    assert alert1.fingerprint == alert2.fingerprint


def test_alert_fingerprint_differs_for_different_context() -> None:
    alert1 = Alert(
        alert_type="preflight_blocked",
        severity="critical",
        message="Test",
        context={"key": "value1"},
        first_seen="2026-04-06T00:00:00Z",
        last_seen="2026-04-06T00:00:00Z",
    )
    alert2 = Alert(
        alert_type="preflight_blocked",
        severity="critical",
        message="Test",
        context={"key": "value2"},
        first_seen="2026-04-06T00:00:00Z",
        last_seen="2026-04-06T00:00:00Z",
    )
    assert alert1.fingerprint != alert2.fingerprint


def test_severity_map_complete() -> None:
    expected_types = {
        "preflight_blocked",
        "dataset_missing",
        "artifact_quarantine",
        "round_failed",
        "round_crashed",
        "round_timed_out",
        "round_stalled",
        "plateau_warning",
        "converged",
        "telemetry_corruption",
    }
    assert set(ALERT_SEVERITY_MAP.keys()) == expected_types


def test_preflight_blocked_alert_critical() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["preflight_overall"] = "fail"

    result = engine.evaluate(state)

    assert len(result["new_alerts"]) == 1
    alert = result["new_alerts"][0]
    assert alert["alert_type"] == "preflight_blocked"
    assert alert["severity"] == "critical"
    assert "blocked" in alert["message"].lower()


def test_dataset_missing_alert_warning() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["preflight_checks"] = [
        {"check": "datasets", "status": "warn", "message": "HB_DATASET not set"}
    ]

    result = engine.evaluate(state)

    assert len(result["new_alerts"]) == 1
    alert = result["new_alerts"][0]
    assert alert["alert_type"] == "dataset_missing"
    assert alert["severity"] == "warning"


def test_artifact_quarantine_alert_critical() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["quarantine_status"] = "quarantine_halt"
    state["quarantine_summary"] = {
        "decision": "halt",
        "corrupted_count": 1,
        "incomplete_count": 0,
        "reason": "corrupted_threshold_reached",
    }

    result = engine.evaluate(state)

    assert len(result["new_alerts"]) == 1
    alert = result["new_alerts"][0]
    assert alert["alert_type"] == "artifact_quarantine"
    assert alert["severity"] == "critical"
    assert "quarantine_halt" in alert["message"]


def test_round_failed_alert_warning() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["subprocess_status"] = "nonzero_exit"
    state["round_index"] = 5

    result = engine.evaluate(state)

    assert len(result["new_alerts"]) == 1
    alert = result["new_alerts"][0]
    assert alert["alert_type"] == "round_failed"
    assert alert["severity"] == "warning"
    assert "5" in alert["message"]


def test_round_crashed_alert_critical() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["subprocess_status"] = "crashed"
    state["round_index"] = 3

    result = engine.evaluate(state)

    assert len(result["new_alerts"]) == 1
    alert = result["new_alerts"][0]
    assert alert["alert_type"] == "round_crashed"
    assert alert["severity"] == "critical"


def test_round_timed_out_alert_warning() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["subprocess_status"] = "timeout"
    state["round_index"] = 2

    result = engine.evaluate(state)

    assert len(result["new_alerts"]) == 1
    alert = result["new_alerts"][0]
    assert alert["alert_type"] == "round_timed_out"
    assert alert["severity"] == "warning"


def test_round_stalled_alert_warning() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["consecutive_empty_rounds"] = 2
    state["max_empty_rounds"] = 3

    result = engine.evaluate(state)

    assert len(result["new_alerts"]) == 1
    alert = result["new_alerts"][0]
    assert alert["alert_type"] == "round_stalled"
    assert alert["severity"] == "warning"


def test_plateau_warning_alert_info() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["plateau_warning"] = True
    state["rolling_metrics"] = {"yield": 0.4, "coverage": 0.5}

    result = engine.evaluate(state)

    assert len(result["new_alerts"]) == 1
    alert = result["new_alerts"][0]
    assert alert["alert_type"] == "plateau_warning"
    assert alert["severity"] == "info"


def test_converged_alert_info() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["convergence_state"] = "converged_plateau"

    result = engine.evaluate(state)

    assert len(result["new_alerts"]) == 1
    alert = result["new_alerts"][0]
    assert alert["alert_type"] == "converged"
    assert alert["severity"] == "info"


def test_telemetry_corruption_alert_critical() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["telemetry_corruption"] = True

    result = engine.evaluate(state)

    assert len(result["new_alerts"]) == 1
    alert = result["new_alerts"][0]
    assert alert["alert_type"] == "telemetry_corruption"
    assert alert["severity"] == "critical"


def test_duplicate_alert_suppressed_within_cooldown() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["plateau_warning"] = True

    result1 = engine.evaluate(state)
    assert len(result1["new_alerts"]) == 1

    result2 = engine.evaluate(state)
    assert len(result2["new_alerts"]) == 0
    assert len(result2["active_alerts"]) == 1


def test_alert_resolved_when_condition_clears() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["plateau_warning"] = True

    result1 = engine.evaluate(state)
    assert len(result1["new_alerts"]) == 1
    assert len(result1["active_alerts"]) == 1

    state["plateau_warning"] = False
    result2 = engine.evaluate(state)
    assert len(result2["new_alerts"]) == 0
    assert len(result2["active_alerts"]) == 0
    assert len(result2["resolved_alerts"]) == 1
    assert result2["resolved_alerts"][0]["resolved_at"] is not None


def test_multiple_alert_types_simultaneously() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["preflight_overall"] = "fail"
    state["subprocess_status"] = "crashed"
    state["round_index"] = 2
    state["plateau_warning"] = True

    result = engine.evaluate(state)

    assert len(result["new_alerts"]) == 3
    alert_types = {a["alert_type"] for a in result["new_alerts"]}
    assert alert_types == {"preflight_blocked", "round_crashed", "plateau_warning"}


def test_active_alerts_restored_from_state() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    existing_alert = {
        "alert_type": "plateau_warning",
        "severity": "info",
        "message": "Test",
        "context": {"rolling_metrics": {"yield": 0.4}, "plateau_warning": True},
        "first_seen": "2026-04-06T00:00:00Z",
        "last_seen": "2026-04-06T00:00:00Z",
    }
    state = _minimal_state()
    state["active_alerts"] = [existing_alert]
    state["plateau_warning"] = True
    state["rolling_metrics"] = {"yield": 0.4}

    result = engine.evaluate(state)

    assert len(result["new_alerts"]) == 0
    assert len(result["active_alerts"]) == 1


def test_emit_alerts_creates_jsonl_and_stderr(capsys: pytest.CaptureFixture) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        engine = AlertEngine({"dedupe_window_sec": 900})
        outdir = Path(tmpdir)

        alerts = [
            {
                "alert_type": "preflight_blocked",
                "severity": "critical",
                "message": "Test critical alert",
                "context": {},
                "first_seen": "2026-04-06T00:00:00Z",
                "last_seen": "2026-04-06T00:00:00Z",
                "fingerprint": "abc123",
            }
        ]

        engine.emit_alerts(alerts, outdir)

        alerts_file = outdir / "alerts.jsonl"
        assert alerts_file.exists()

        with open(alerts_file) as f:
            lines = f.readlines()
        assert len(lines) == 1
        parsed = json.loads(lines[0])
        assert parsed["alert_type"] == "preflight_blocked"

        captured = capsys.readouterr()
        assert "[CRITICAL] preflight_blocked" in captured.err
        assert "Test critical alert" in captured.err


def test_acceptance_preflight_blocked_critical() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = {
        "campaign_state": "RUNNING",
        "preflight_overall": "fail",
        "active_alerts": [],
    }
    result = engine.evaluate(state)
    assert any(
        a["alert_type"] == "preflight_blocked" and a["severity"] == "critical"
        for a in result["new_alerts"]
    )


def test_acceptance_plateau_dedupe() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = {
        "campaign_state": "RUNNING",
        "plateau_warning": True,
        "active_alerts": [],
    }
    result1 = engine.evaluate(state)
    result2 = engine.evaluate(state)
    assert len(result1["new_alerts"]) == 1
    assert len(result2["new_alerts"]) == 0


def test_no_alerts_on_healthy_state() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = {
        "campaign_state": "RUNNING",
        "preflight_overall": "pass",
        "subprocess_status": "success",
        "active_alerts": [],
    }
    result = engine.evaluate(state)
    assert len(result["new_alerts"]) == 0
    assert len(result["active_alerts"]) == 0


def test_alert_context_preserved() -> None:
    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["subprocess_status"] = "timeout"
    state["round_index"] = 7

    result = engine.evaluate(state)

    alert = result["new_alerts"][0]
    assert alert["context"]["round_index"] == 7
    assert alert["context"]["subprocess_status"] == "timeout"


def test_alert_persistence_to_sqlite() -> None:
    import sqlite3
    from autoresearch.harness.telemetry_store import init_db

    engine = AlertEngine({"dedupe_window_sec": 900})
    state = _minimal_state()
    state["campaign_state"] = "RUNNING"
    state["plateau_warning"] = True
    state["rolling_metrics"] = {"yield": 0.4}

    result = engine.evaluate(state)
    assert len(result["new_alerts"]) == 1

    # Create temporary database
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "test.db"
        db = init_db(db_path)

        # Persist alerts to SQLite
        engine.persist_alerts_to_sqlite(db, "test-campaign")

        # Verify alert was persisted
        cursor = db.cursor()
        cursor.execute("SELECT COUNT(*) FROM alerts WHERE campaign_id = ? AND status = ?",
                       ("test-campaign", "active"))
        count = cursor.fetchone()[0]
        assert count == 1

        # Verify alert details
        cursor.execute("SELECT alert_type, severity, status FROM alerts WHERE campaign_id = ?",
                       ("test-campaign",))
        row = cursor.fetchone()
        assert row[0] == "plateau_warning"
        assert row[1] == "info"
        assert row[2] == "active"

        db.close()


def test_alert_load_from_sqlite() -> None:
    import sqlite3
    from autoresearch.harness.telemetry_store import init_db

    # Create database with pre-existing alert
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = Path(tmpdir) / "test.db"
        db = init_db(db_path)

        cursor = db.cursor()
        cursor.execute(
            """INSERT INTO alerts
            (campaign_id, fingerprint, alert_type, severity, message, context,
             first_seen, last_seen, resolved_at, status)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"""
,
            (
                "test-campaign",
                "abc123def456",
                "plateau_warning",
                "info",
                "Test plateau",
                json.dumps({"rolling_metrics": {"yield": 0.4}, "plateau_warning": True}),
                "2026-04-06T00:00:00Z",
                "2026-04-06T00:00:00Z",
                None,
                "active",
            ),
        )
        db.commit()

        # Load alerts from SQLite
        engine = AlertEngine({"dedupe_window_sec": 900})
        loaded_alerts = engine.load_alerts_from_sqlite(db, "test-campaign")

        assert len(loaded_alerts) == 1
        alert = loaded_alerts[0]
        assert alert["alert_type"] == "plateau_warning"
        assert alert["severity"] == "info"
        assert alert["message"] == "Test plateau"
        assert alert["fingerprint"] == "abc123def456"

        db.close()
