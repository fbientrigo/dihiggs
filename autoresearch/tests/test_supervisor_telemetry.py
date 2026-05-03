"""Unit tests for SQLite telemetry store and ingestion.

Test coverage:
- Empty replay produces zero rows
- Single valid event ingests correctly
- Double replay is idempotent (no duplicates)
- Malformed JSON line is quarantined, doesn't crash ingestion
- Corrupted checksum detected and logged
- WAL mode verification
- Multi-file replay
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest

from autoresearch.harness.autonomy_scheduler import FixedIntervalAutonomyScheduler
from autoresearch.harness.telemetry_ingest import ingest_jsonl_file, replay_sources
from autoresearch.harness.reconcile_watermark import get_watermark, update_watermark
from autoresearch.harness.telemetry_store import db_counts, init_db


def test_init_db_creates_wal_mode(tmp_path: Path) -> None:
    """Verify WAL mode is enabled on database initialization."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    cursor = conn.cursor()
    cursor.execute("PRAGMA journal_mode")
    mode = cursor.fetchone()[0]
    
    assert mode.upper() == "WAL"
    conn.close()


def test_init_db_creates_schema(tmp_path: Path) -> None:
    """Verify all required tables are created."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    cursor = conn.cursor()
    cursor.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
    )
    tables = [row[0] for row in cursor.fetchall()]
    
    assert "events" in tables
    assert "ingestion_errors" in tables
    assert "ingestion_log" in tables
    assert "schema_metadata" in tables
    
    conn.close()


def test_init_db_stores_schema_version(tmp_path: Path) -> None:
    """Verify schema version is recorded."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    cursor = conn.cursor()
    cursor.execute(
        "SELECT value FROM schema_metadata WHERE key = 'schema_version'"
    )
    version = cursor.fetchone()[0]
    
    assert version == "1"
    conn.close()


def test_db_counts_empty_database(tmp_path: Path) -> None:
    """Empty database should return zero counts."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    counts = db_counts(conn)
    
    assert counts["events"] == 0
    assert counts["ingestion_errors"] == 0
    assert counts["ingestion_log"] == 0
    
    conn.close()


def test_ingest_empty_replay_produces_zero_rows(tmp_path: Path) -> None:
    """Acceptance criterion: empty replay produces zero rows."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    # Replay empty list of sources
    results = replay_sources(conn, [])
    
    assert len(results) == 0
    
    counts = db_counts(conn)
    assert counts["events"] == 0
    assert counts["ingestion_errors"] == 0
    assert counts["ingestion_log"] == 0
    
    conn.close()


def test_ingest_single_valid_event(tmp_path: Path) -> None:
    """Single valid event should ingest correctly."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    # Create JSONL file with single event
    jsonl_path = tmp_path / "events.jsonl"
    event = {
        "schema_version": 1,
        "campaign_id": "test-campaign",
        "event_type": "ATTEMPT_EVALUATED",
        "utc": "2026-04-06T12:00:00Z",
        "payload": {
            "arm_id": "adaptive-v1",
            "attempt_id": "test123",
            "iter_index": 0,
            "attempt_index": 0,
            "cell_id": "tb=1000|bin=5",
            "eval_status": "success",
            "successes": 42.0,
            "trials": 100,
            "elapsed_sec": 0.5,
            "axes_binned": {"tb": 1000, "lam1_bin": 5},
            "axes_raw": {"tb": 1000, "lam1": -12.5},
        },
    }
    
    with jsonl_path.open("w", encoding="utf-8") as f:
        f.write(json.dumps(event) + "\n")
    
    # Ingest
    result = ingest_jsonl_file(conn, jsonl_path)
    
    assert result["status"] == "success"
    assert result["lines_read"] == 1
    assert result["lines_ingested"] == 1
    assert result["lines_skipped"] == 0
    assert result["lines_errored"] == 0
    
    counts = db_counts(conn)
    assert counts["events"] == 1
    assert counts["ingestion_errors"] == 0
    assert counts["ingestion_log"] == 1
    
    # Verify event content
    cursor = conn.cursor()
    cursor.execute("SELECT campaign_id, event_type, source_file FROM events")
    row = cursor.fetchone()
    
    assert row[0] == "test-campaign"
    assert row[1] == "ATTEMPT_EVALUATED"
    assert str(jsonl_path) in row[2]
    
    conn.close()


def test_ingest_double_replay_is_idempotent(tmp_path: Path) -> None:
    """Acceptance criterion: double replay produces no duplicates."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    # Create JSONL file with two events
    jsonl_path = tmp_path / "events.jsonl"
    events = [
        {
            "schema_version": 1,
            "campaign_id": "test-campaign",
            "event_type": "ATTEMPT_EVALUATED",
            "utc": "2026-04-06T12:00:00Z",
            "payload": {"attempt_id": "event001"},
        },
        {
            "schema_version": 1,
            "campaign_id": "test-campaign",
            "event_type": "ATTEMPT_EVALUATED",
            "utc": "2026-04-06T12:01:00Z",
            "payload": {"attempt_id": "event002"},
        },
    ]
    
    with jsonl_path.open("w", encoding="utf-8") as f:
        for event in events:
            f.write(json.dumps(event) + "\n")
    
    # First replay
    result1 = ingest_jsonl_file(conn, jsonl_path)
    assert result1["lines_ingested"] == 2
    
    counts1 = db_counts(conn)
    assert counts1["events"] == 2
    
    # Second replay (idempotency test)
    result2 = ingest_jsonl_file(conn, jsonl_path)
    assert result2["lines_ingested"] == 0  # No new rows
    assert result2["lines_skipped"] == 2   # Both skipped as duplicates
    
    counts2 = db_counts(conn)
    assert counts2["events"] == 2  # Row count stable
    assert counts2["ingestion_log"] == 2  # Two ingestion runs logged
    
    conn.close()


def test_ingest_malformed_json_quarantined(tmp_path: Path) -> None:
    """Malformed JSON line should be logged, not crash ingestion."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    # Create JSONL with malformed line
    jsonl_path = tmp_path / "events.jsonl"
    with jsonl_path.open("w", encoding="utf-8") as f:
        f.write('{"valid": "event", "event_type": "TEST", "campaign_id": "test"}\n')
        f.write("THIS IS NOT JSON\n")  # Malformed line
        f.write('{"valid": "event2", "event_type": "TEST", "campaign_id": "test"}\n')
    
    # Ingest
    result = ingest_jsonl_file(conn, jsonl_path)
    
    assert result["status"] == "success"
    assert result["lines_read"] == 3
    assert result["lines_ingested"] == 2  # Two valid events
    assert result["lines_errored"] == 1   # One malformed line
    
    counts = db_counts(conn)
    assert counts["events"] == 2
    assert counts["ingestion_errors"] == 1
    
    # Verify error was logged
    cursor = conn.cursor()
    cursor.execute("SELECT error_type, error_message FROM ingestion_errors")
    error_row = cursor.fetchone()
    
    assert error_row[0] == "json_decode_error"
    assert "Expecting" in error_row[1]  # JSONDecodeError message
    
    conn.close()


def test_ingest_missing_required_fields(tmp_path: Path) -> None:
    """Events missing required fields should be quarantined."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    # Create JSONL with incomplete event
    jsonl_path = tmp_path / "events.jsonl"
    incomplete_event = {
        "schema_version": 1,
        # Missing campaign_id and event_type
        "utc": "2026-04-06T12:00:00Z",
    }
    
    with jsonl_path.open("w", encoding="utf-8") as f:
        f.write(json.dumps(incomplete_event) + "\n")
    
    # Ingest (should not crash)
    result = ingest_jsonl_file(conn, jsonl_path)
    
    # Event should use fallback values for missing fields
    # (event_type="UNKNOWN", campaign_id="UNKNOWN")
    assert result["status"] == "success"
    assert result["lines_ingested"] == 1
    
    # Verify default values
    cursor = conn.cursor()
    cursor.execute("SELECT campaign_id, event_type FROM events")
    row = cursor.fetchone()
    
    assert row[0] == "UNKNOWN"
    assert row[1] == "UNKNOWN"
    
    conn.close()


def test_ingest_multiple_files(tmp_path: Path) -> None:
    """Replay multiple JSONL files in single operation."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    # Create two JSONL files with events format
    file1 = tmp_path / "events.jsonl"
    file2 = tmp_path / "other_events.jsonl"
    
    event1 = {
        "schema_version": 1,
        "campaign_id": "test",
        "event_type": "TYPE_A",
        "utc": "2026-04-06T12:00:00Z",
        "payload": {},
    }
    
    event2 = {
        "schema_version": 1,
        "campaign_id": "test",
        "event_type": "TYPE_B",
        "utc": "2026-04-06T12:01:00Z",
        "payload": {},
    }
    
    with file1.open("w", encoding="utf-8") as f:
        f.write(json.dumps(event1) + "\n")
    
    with file2.open("w", encoding="utf-8") as f:
        f.write(json.dumps(event2) + "\n")
    
    # Replay both files
    results = replay_sources(conn, [file1, file2])
    
    assert len(results) == 2
    assert results[0]["lines_ingested"] == 1
    assert results[1]["lines_ingested"] == 1
    
    counts = db_counts(conn)
    assert counts["events"] == 2
    assert counts["ingestion_log"] == 2
    
    conn.close()


def test_ingest_preserves_source_linkage(tmp_path: Path) -> None:
    """Verify source_file and line_offset are stored correctly."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    jsonl_path = tmp_path / "events.jsonl"
    events = [
        {"schema_version": 1, "campaign_id": "test", "event_type": "E1", "utc": "2026-04-06T12:00:00Z", "payload": {}},
        {"schema_version": 1, "campaign_id": "test", "event_type": "E2", "utc": "2026-04-06T12:00:00Z", "payload": {}},
        {"schema_version": 1, "campaign_id": "test", "event_type": "E3", "utc": "2026-04-06T12:00:00Z", "payload": {}},
    ]
    
    with jsonl_path.open("w", encoding="utf-8") as f:
        for event in events:
            f.write(json.dumps(event) + "\n")
    
    ingest_jsonl_file(conn, jsonl_path)
    
    cursor = conn.cursor()
    cursor.execute("SELECT source_file, line_offset, event_type FROM events ORDER BY line_offset")
    rows = cursor.fetchall()
    
    assert len(rows) == 3
    assert all(str(jsonl_path) in row[0] for row in rows)
    assert rows[0][1] == 1
    assert rows[1][1] == 2
    assert rows[2][1] == 3
    assert rows[0][2] == "E1"
    assert rows[1][2] == "E2"
    assert rows[2][2] == "E3"
    
    conn.close()


def test_ingest_checksum_stored(tmp_path: Path) -> None:
    """Verify checksum is computed and stored."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    jsonl_path = tmp_path / "events.jsonl"
    event = {
        "schema_version": 1,
        "campaign_id": "test",
        "event_type": "TEST",
        "utc": "2026-04-06T12:00:00Z",
        "payload": {},
    }
    
    with jsonl_path.open("w", encoding="utf-8") as f:
        f.write(json.dumps(event) + "\n")
    
    ingest_jsonl_file(conn, jsonl_path)
    
    cursor = conn.cursor()
    cursor.execute("SELECT checksum FROM events")
    checksum = cursor.fetchone()[0]
    
    # Checksum should be 64-char hex (SHA256)
    assert len(checksum) == 64
    assert all(c in "0123456789abcdef" for c in checksum)
    
    conn.close()


def test_ingest_file_not_found(tmp_path: Path) -> None:
    """Missing source file should return file_not_found status."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    missing_file = tmp_path / "nonexistent.jsonl"
    
    result = ingest_jsonl_file(conn, missing_file)
    
    assert result["status"] == "file_not_found"
    assert result["lines_read"] == 0
    assert result["lines_ingested"] == 0
    
    conn.close()


def test_ingest_empty_lines_skipped(tmp_path: Path) -> None:
    """Empty lines in JSONL should be skipped gracefully."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    jsonl_path = tmp_path / "events.jsonl"
    with jsonl_path.open("w", encoding="utf-8") as f:
        f.write("\n")
        f.write('{"schema_version": 1, "campaign_id": "test", "event_type": "TEST", "utc": "2026-04-06T12:00:00Z", "payload": {}}\n')
        f.write("\n")
        f.write("\n")
        f.write('{"schema_version": 1, "campaign_id": "test", "event_type": "TEST", "utc": "2026-04-06T12:00:00Z", "payload": {}}\n')
    
    result = ingest_jsonl_file(conn, jsonl_path)
    
    assert result["lines_read"] == 5
    assert result["lines_ingested"] == 2
    assert result["lines_skipped"] == 3  # Three empty lines
    assert result["lines_errored"] == 0
    
    conn.close()


def test_ingest_indexes_created(tmp_path: Path) -> None:
    """Verify indexes are created for query performance."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    cursor = conn.cursor()
    cursor.execute(
        "SELECT name FROM sqlite_master WHERE type='index' AND tbl_name='events'"
    )
    indexes = [row[0] for row in cursor.fetchall()]
    
    assert "idx_events_campaign" in indexes
    assert "idx_events_dedupe" in indexes
    
    conn.close()


def test_ingest_task_summary_jsonl(tmp_path: Path) -> None:
    """Ingest task_summary.jsonl records into task_summaries table."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    # Create task_summary.jsonl
    task_summary_path = tmp_path / "task_summary.jsonl"
    tasks = [
        {
            "campaign_id": "test-campaign",
            "task_id": "task-001",
            "task_type": "adaptive",
            "status": "completed",
            "payload": {"result": "success"},
        },
        {
            "campaign_id": "test-campaign",
            "task_id": "task-002",
            "task_type": "branch",
            "status": "failed",
            "payload": {"error": "timeout"},
        },
    ]
    
    with task_summary_path.open("w", encoding="utf-8") as f:
        for task in tasks:
            f.write(json.dumps(task) + "\n")
    
    # Ingest task summaries
    from autoresearch.harness.telemetry_ingest import ingest_task_summary_jsonl
    result = ingest_task_summary_jsonl(conn, task_summary_path)
    
    assert result["status"] == "success"
    assert result["lines_ingested"] == 2
    
    counts = db_counts(conn)
    assert counts["task_summaries"] == 2
    assert counts["events"] == 0  # Should not be in events table
    
    # Verify task records
    cursor = conn.cursor()
    cursor.execute("SELECT task_id, task_type, status FROM task_summaries ORDER BY task_id")
    rows = cursor.fetchall()
    
    assert len(rows) == 2
    assert rows[0][0] == "task-001"
    assert rows[0][1] == "adaptive"
    assert rows[0][2] == "completed"
    assert rows[1][0] == "task-002"
    assert rows[1][1] == "branch"
    assert rows[1][2] == "failed"
    
    conn.close()


def test_ingest_task_summary_double_replay_idempotent(tmp_path: Path) -> None:
    """Task summary replay should be idempotent."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    task_summary_path = tmp_path / "task_summary.jsonl"
    task = {
        "campaign_id": "test",
        "task_id": "task-001",
        "task_type": "test",
        "status": "done",
    }
    
    with task_summary_path.open("w", encoding="utf-8") as f:
        f.write(json.dumps(task) + "\n")
    
    from autoresearch.harness.telemetry_ingest import ingest_task_summary_jsonl
    
    # First replay
    result1 = ingest_task_summary_jsonl(conn, task_summary_path)
    assert result1["lines_ingested"] == 1
    
    counts1 = db_counts(conn)
    assert counts1["task_summaries"] == 1
    
    # Second replay (should not ingest duplicate)
    result2 = ingest_task_summary_jsonl(conn, task_summary_path)
    assert result2["lines_ingested"] == 0
    assert result2["lines_skipped"] == 1
    
    counts2 = db_counts(conn)
    assert counts2["task_summaries"] == 1  # Still one record
    
    conn.close()


def test_ingest_orchestrator_log(tmp_path: Path) -> None:
    """Ingest plain-text orchestrator.log entries."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    # Create orchestrator.log
    log_path = tmp_path / "orchestrator.log"
    log_lines = [
        "[2026-04-06 12:00:00] [INIT] Starting campaign",
        "[2026-04-06 12:00:01] [RUN] Processing arm-1",
        "[2026-04-06 12:00:02] [DONE] Completed 10 discoveries",
        "[2026-04-06 12:00:03] [ERROR] Connection timeout",
    ]
    
    with log_path.open("w", encoding="utf-8") as f:
        for line in log_lines:
            f.write(line + "\n")
    
    # Ingest orchestrator log
    from autoresearch.harness.telemetry_ingest import ingest_orchestrator_log
    result = ingest_orchestrator_log(conn, log_path, campaign_id="test-campaign")
    
    assert result["status"] == "success"
    assert result["lines_ingested"] == 4
    
    counts = db_counts(conn)
    assert counts["orchestrator_log_entries"] == 4
    
    # Verify log entries
    cursor = conn.cursor()
    cursor.execute(
        "SELECT log_level, log_message FROM orchestrator_log_entries ORDER BY line_offset"
    )
    rows = cursor.fetchall()
    
    assert len(rows) == 4
    assert rows[0][0] == "INIT"
    assert "Starting campaign" in rows[0][1]
    assert rows[1][0] == "RUN"
    assert rows[3][0] == "ERROR"
    
    conn.close()


def test_ingest_artifact_metadata(tmp_path: Path) -> None:
    """Ingest stdout/stderr artifact metadata."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    # Create mock artifact files
    stdout_file = tmp_path / "round_001.stdout"
    stderr_file = tmp_path / "round_001.stderr"
    stdout_file.write_text("Some output")
    stderr_file.write_text("Some error")
    
    # Ingest metadata
    from autoresearch.harness.telemetry_ingest import ingest_artifact_metadata
    result1 = ingest_artifact_metadata(
        conn,
        artifact_type="stdout",
        artifact_path=stdout_file,
        campaign_id="test-campaign",
        task_id="task-001",
        source_file=tmp_path / "events.jsonl",
        line_offset=5,
    )
    
    result2 = ingest_artifact_metadata(
        conn,
        artifact_type="stderr",
        artifact_path=stderr_file,
        campaign_id="test-campaign",
        task_id="task-001",
        source_file=tmp_path / "events.jsonl",
        line_offset=6,
    )
    
    assert result1["status"] == "success"
    assert result2["status"] == "success"
    
    counts = db_counts(conn)
    assert counts["artifact_metadata"] == 2
    
    # Verify metadata records
    cursor = conn.cursor()
    cursor.execute(
        "SELECT artifact_type, related_task_id FROM artifact_metadata ORDER BY artifact_type"
    )
    rows = cursor.fetchall()
    
    assert len(rows) == 2
    assert rows[0][0] == "stderr"
    assert rows[1][0] == "stdout"
    assert rows[0][1] == "task-001"
    
    conn.close()


def test_replay_sources_auto_detects_file_types(tmp_path: Path) -> None:
    """replay_sources should auto-detect and route different file types."""
    db_path = tmp_path / "telemetry.db"
    conn = init_db(db_path)
    
    # Create events.jsonl
    events_path = tmp_path / "events.jsonl"
    events_path.write_text(
        json.dumps({
            "campaign_id": "test",
            "event_type": "ATTEMPT_EVALUATED",
            "utc": "2026-04-06T12:00:00Z",
            "payload": {},
        }) + "\n"
    )
    
    # Create task_summary.jsonl
    task_summary_path = tmp_path / "task_summary.jsonl"
    task_summary_path.write_text(
        json.dumps({
            "campaign_id": "test",
            "task_id": "task-001",
            "task_type": "adaptive",
            "status": "done",
        }) + "\n"
    )
    
    # Create orchestrator.log
    log_path = tmp_path / "orchestrator.log"
    log_path.write_text("[2026-04-06 12:00:00] [INIT] Starting\n")
    
    # Replay all sources
    results = replay_sources(conn, [events_path, task_summary_path, log_path])
    
    assert len(results) == 3
    assert results[0]["status"] == "success"  # events.jsonl
    assert results[1]["status"] == "success"  # task_summary.jsonl
    assert results[2]["status"] == "success"  # orchestrator.log
    
    counts = db_counts(conn)
    assert counts["events"] == 1
    assert counts["task_summaries"] == 1
    assert counts["orchestrator_log_entries"] == 1

    conn.close()


def test_reconcile_incremental_replay_updates_only_delta_rows(tmp_path: Path) -> None:
    """Second reconcile after a new append should ingest only the delta rows."""
    db_path = tmp_path / "telemetry.db"
    events_path = tmp_path / "events.jsonl"
    first_event = {
        "schema_version": 1,
        "campaign_id": "test-campaign",
        "event_type": "ATTEMPT_EVALUATED",
        "utc": "2026-04-13T12:00:00Z",
        "payload": {"attempt_id": "attempt-001", "arm_id": "adaptive-v1"},
    }
    second_event = {
        "schema_version": 1,
        "campaign_id": "test-campaign",
        "event_type": "ATTEMPT_EVALUATED",
        "utc": "2026-04-13T12:05:00Z",
        "payload": {"attempt_id": "attempt-002", "arm_id": "adaptive-v1"},
    }
    events_path.write_text(json.dumps(first_event) + "\n", encoding="utf-8")

    conn = init_db(db_path)
    first_results = replay_sources(conn, [events_path])

    assert len(first_results) == 1
    assert first_results[0]["status"] == "success"
    assert first_results[0]["lines_ingested"] == 1
    assert first_results[0]["lines_skipped"] == 0

    cursor = conn.cursor()
    cursor.execute(
        "SELECT line_offset, checksum FROM events WHERE source_file = ? ORDER BY line_offset DESC LIMIT 1",
        (str(events_path),),
    )
    first_offset, first_checksum = cursor.fetchone()
    update_watermark(conn, events_path, first_offset, first_checksum, lambda: "2026-04-13T12:00:01Z")
    first_watermark = get_watermark(conn, events_path)

    assert first_watermark == {
        "last_line_offset": 1,
        "last_checksum": first_checksum,
        "updated_at": "2026-04-13T12:00:01Z",
    }
    conn.close()

    with events_path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(second_event) + "\n")

    conn = init_db(db_path)
    rehydrated_watermark = get_watermark(conn, events_path)
    assert rehydrated_watermark == first_watermark

    second_results = replay_sources(conn, [events_path])

    assert len(second_results) == 1
    assert second_results[0]["status"] == "success"
    assert second_results[0]["lines_ingested"] == 1
    assert second_results[0]["lines_skipped"] == 1

    cursor = conn.cursor()
    cursor.execute(
        "SELECT line_offset, checksum FROM events WHERE source_file = ? ORDER BY line_offset DESC LIMIT 1",
        (str(events_path),),
    )
    second_offset, second_checksum = cursor.fetchone()
    update_watermark(conn, events_path, second_offset, second_checksum, lambda: "2026-04-13T12:05:01Z")

    counts = db_counts(conn)
    assert counts["events"] == 2
    assert counts["ingestion_log"] == 2
    assert get_watermark(conn, events_path) == {
        "last_line_offset": 2,
        "last_checksum": second_checksum,
        "updated_at": "2026-04-13T12:05:01Z",
    }

    conn.close()


def test_reconcile_replay_noop_when_sources_are_unchanged(tmp_path: Path) -> None:
    """Second reconcile over unchanged sources should be a no-op replay."""
    db_path = tmp_path / "telemetry.db"
    events_path = tmp_path / "events.jsonl"
    event = {
        "schema_version": 1,
        "campaign_id": "test-campaign",
        "event_type": "ATTEMPT_EVALUATED",
        "utc": "2026-04-13T12:00:00Z",
        "payload": {"attempt_id": "attempt-001", "arm_id": "adaptive-v1"},
    }
    events_path.write_text(json.dumps(event) + "\n", encoding="utf-8")

    conn = init_db(db_path)
    first_results = replay_sources(conn, [events_path])

    assert len(first_results) == 1
    assert first_results[0]["lines_ingested"] == 1
    assert first_results[0]["lines_skipped"] == 0

    cursor = conn.cursor()
    cursor.execute(
        "SELECT line_offset, checksum FROM events WHERE source_file = ? ORDER BY line_offset DESC LIMIT 1",
        (str(events_path),),
    )
    line_offset, checksum = cursor.fetchone()
    update_watermark(conn, events_path, line_offset, checksum, lambda: "2026-04-13T12:00:01Z")
    conn.close()

    conn = init_db(db_path)
    assert get_watermark(conn, events_path) == {
        "last_line_offset": 1,
        "last_checksum": checksum,
        "updated_at": "2026-04-13T12:00:01Z",
    }

    second_results = replay_sources(conn, [events_path])

    assert len(second_results) == 1
    assert second_results[0]["status"] == "success"
    assert second_results[0]["lines_ingested"] == 0
    assert second_results[0]["lines_skipped"] == 1

    counts = db_counts(conn)
    assert counts["events"] == 1
    assert counts["ingestion_log"] == 2
    assert counts["reconcile_watermarks"] == 1

    conn.close()


def test_scan_datalake_happy_path(tmp_path: Path) -> None:
    """Scanner should enumerate campaign artifacts deterministically."""
    campaign_dir = tmp_path / "campaign-alpha"
    campaign_dir.mkdir()

    (campaign_dir / "events.jsonl").write_text(
        json.dumps(
            {
                "campaign_id": "campaign-alpha",
                "event_type": "ATTEMPT_EVALUATED",
                "utc": "2026-04-13T12:00:00Z",
                "payload": {"attempt_id": "attempt-001", "arm_id": "adaptive-v1"},
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (campaign_dir / "task_summary.jsonl").write_text(
        json.dumps(
            {
                "campaign_id": "campaign-alpha",
                "task_id": "task-001",
                "task_type": "adaptive",
                "status": "completed",
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (campaign_dir / "orchestrator.log").write_text("[2026-04-13 12:00:00] [INIT] start\n", encoding="utf-8")
    (campaign_dir / "campaign_status.json").write_text(
        json.dumps({"campaign_state": "RUNNING"}),
        encoding="utf-8",
    )

    from autoresearch.harness.datalake_manifest import scan_datalake

    report = scan_datalake(tmp_path)

    assert report["status"] == "success"
    assert report["campaign_count"] == 1
    assert report["manifest_row_count"] == 4
    assert report["parser_status_counts"] == {"success": 4}

    rows = report["manifest"]
    assert [row["artifact_name"] for row in rows] == [
        "events.jsonl",
        "task_summary.jsonl",
        "orchestrator.log",
        "campaign_status.json",
    ]
    assert all(row["campaign_dir"] == str(campaign_dir) for row in rows)
    assert all(row["source_mtime"].endswith("Z") for row in rows)
    assert all(row["records_parsed"] >= 1 for row in rows)


def test_scan_datalake_missing_artifact_quarantines_cleanly(tmp_path: Path) -> None:
    """Missing required artifacts should be flagged without crashing the scan."""
    campaign_dir = tmp_path / "campaign-beta"
    campaign_dir.mkdir()
    (campaign_dir / "task_summary.jsonl").write_text(
        json.dumps(
            {
                "campaign_id": "campaign-beta",
                "task_id": "task-002",
                "task_type": "branch",
                "status": "running",
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (campaign_dir / "campaign_state.json").write_text("{not valid json}\n", encoding="utf-8")

    from autoresearch.harness.datalake_manifest import scan_datalake

    report = scan_datalake(tmp_path)

    assert report["status"] == "success_with_quarantine"
    assert report["campaign_count"] == 1
    assert report["manifest_row_count"] == 3
    assert report["parser_status_counts"] == {
        "missing": 1,
        "quarantined_malformed": 1,
        "success": 1,
    }

    rows_by_name = {row["artifact_name"]: row for row in report["manifest"]}
    assert rows_by_name["events.jsonl"]["parser_status"] == "missing"
    assert rows_by_name["events.jsonl"]["source_mtime"] is None
    assert rows_by_name["task_summary.jsonl"]["parser_status"] == "success"
    assert rows_by_name["campaign_state.json"]["parser_status"] == "quarantined_malformed"


def test_task13_scheduler_persists_stage_telemetry_and_audit_trail(tmp_path: Path) -> None:
    conn = init_db(tmp_path / "telemetry.db")
    state_path = tmp_path / "scheduler_state.json"
    scheduler = FixedIntervalAutonomyScheduler(
        conn,
        config={
            "search": {
                "tb_values": [1000, 5000],
                "lam1": {"min": -20.0, "max": 20.0, "n_bins": 4},
                "mphi": {"min": 0.0, "max": 1000.0, "n_bins": 1},
            },
            "limits": {"max_new_run_dirs_per_round": 3},
            "metrics": {
                "weights": {"yield": 0.5, "coverage": 0.3, "diversity": 0.2},
                "multi_axis": {
                    "collapse_axes": ["tb", "lam1_bin"],
                    "coverage_axes": [
                        {"name": "tb", "kind": "categorical", "domain": [1000, 5000], "weight": 0.5},
                        {"name": "lam1_bin", "kind": "discrete", "domain_size": 4, "weight": 0.5},
                    ],
                    "diversity_axes": [
                        {"name": "tb", "weight": 0.5},
                        {"name": "lam1_bin", "weight": 0.5},
                    ],
                    "diversity_pairs": [{"axes": ["tb", "lam1_bin"], "weight": 1.0}],
                },
            },
            "autonomy_scheduler": {
                "interval_seconds": 900,
                "max_proposals": 2,
                "max_dispatches": 1,
                "dispatch_mode": "dry_run",
                "quarantine": {"corrupted_threshold": 1, "incomplete_threshold": 1},
            },
        },
        source_file=tmp_path / "scanner_rows.jsonl",
        state_path=state_path,
        utc_now_func=lambda: "2026-04-14T00:00:00+00:00",
        monotonic_func=iter([0.0, 0.1, 0.2]).__next__,
    )

    result = scheduler.tick(
        scanner_rows=[
            {
                "campaign_id": "camp-001",
                "arm_id": "adaptive-v1",
                "run_dir_fingerprint": "run-001",
                "attempt_id": "attempt-001",
                "artifact_status": "complete",
                "discovery_status": "new",
                "axes_binned": {"tb": 1000, "lam1_bin": 0},
            },
            {
                "campaign_id": "camp-001",
                "arm_id": "adaptive-v1",
                "run_dir_fingerprint": "run-002",
                "attempt_id": "attempt-002",
                "artifact_status": "complete",
                "discovery_status": "new",
                "axes_binned": {"tb": 1000, "lam1_bin": 1},
            },
        ]
    )

    snapshot = json.loads(state_path.read_text(encoding="utf-8"))
    stage_names = [entry["stage"] for entry in snapshot["stage_telemetry"]]

    assert result["status"] == "success"
    assert stage_names == ["reconcile", "proposal", "dispatch", "scheduler"]
    assert snapshot["quarantine"]["status"] == "clear"
    assert len(snapshot["proposal_audit"]) == 2
    assert snapshot["proposal_audit"][0]["decision"]["outcome"] == "dispatched"
    assert snapshot["proposal_audit"][0]["rationale"]["signals"]
    assert snapshot["proposal_audit"][1]["decision"]["outcome"] == "deferred"


def test_task13_corrupted_artifact_quarantine_halts_dispatch(tmp_path: Path) -> None:
    conn = init_db(tmp_path / "telemetry.db")
    state_path = tmp_path / "scheduler_state.json"
    scheduler = FixedIntervalAutonomyScheduler(
        conn,
        config={
            "search": {
                "tb_values": [1000],
                "lam1": {"min": -20.0, "max": 20.0, "n_bins": 2},
                "mphi": {"min": 0.0, "max": 1000.0, "n_bins": 1},
            },
            "limits": {"max_new_run_dirs_per_round": 2},
            "metrics": {
                "weights": {"yield": 0.5, "coverage": 0.3, "diversity": 0.2},
                "multi_axis": {
                    "collapse_axes": ["tb", "lam1_bin"],
                    "coverage_axes": [
                        {"name": "tb", "kind": "categorical", "domain": [1000], "weight": 0.5},
                        {"name": "lam1_bin", "kind": "discrete", "domain_size": 2, "weight": 0.5},
                    ],
                    "diversity_axes": [
                        {"name": "tb", "weight": 0.5},
                        {"name": "lam1_bin", "weight": 0.5},
                    ],
                    "diversity_pairs": [{"axes": ["tb", "lam1_bin"], "weight": 1.0}],
                },
            },
            "autonomy_scheduler": {
                "interval_seconds": 900,
                "max_proposals": 1,
                "max_dispatches": 1,
                "dispatch_mode": "dry_run",
                "quarantine": {"corrupted_threshold": 1, "incomplete_threshold": 1},
            },
        },
        source_file=tmp_path / "scanner_rows.jsonl",
        state_path=state_path,
        utc_now_func=lambda: "2026-04-14T00:10:00+00:00",
        monotonic_func=iter([0.0, 0.1, 0.2]).__next__,
    )

    result = scheduler.tick(
        scanner_rows=[
            {
                "campaign_id": "camp-001",
                "arm_id": "adaptive-v1",
                "run_dir_fingerprint": "run-001",
                "attempt_id": "attempt-001",
                "artifact_status": "complete",
                "discovery_status": "new",
                "axes_binned": {"tb": 1000, "lam1_bin": 0},
            },
            {"parser_status": "quarantined_malformed"},
        ]
    )

    snapshot = json.loads(state_path.read_text(encoding="utf-8"))

    assert result["status"] == "quarantine_halt"
    assert result["dispatch_count"] == 0
    assert snapshot["dispatch_result"]["reason"] == "quarantine_halt"
    assert snapshot["quarantine"]["status"] == "quarantine_halt"
    assert snapshot["quarantine"]["corrupted_count"] == 1
    assert snapshot["proposal_audit"][0]["decision"]["reason"] == "quarantine_halt"
