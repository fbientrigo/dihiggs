"""SQLite telemetry store in WAL mode - rebuildable from raw JSONL artifacts.

Schema version: 1

Tables:
- events: Normalized event records from events.jsonl
- ingestion_errors: Malformed/corrupted lines
- ingestion_log: Ingestion run metadata
- reconcile_watermarks: Track last processed offsets for incremental replay

Design principles:
- SQLite in WAL mode for concurrent read access
- Raw JSONL files are source of truth, DB is disposable
- Idempotent replay via composite dedupe keys
- Checksum validation for corruption detection
"""
from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

from autoresearch.harness.telemetry_ingest import replay_sources  # noqa: F401

SCHEMA_VERSION = 1


def init_db(db_path: str | Path) -> sqlite3.Connection:
    """Initialize SQLite database in WAL mode with telemetry schema.
    
    Args:
        db_path: Path to SQLite database file (created if missing)
    
    Returns:
        SQLite connection with WAL mode enabled
    
    Design:
        - WAL mode enables concurrent readers during ingestion
        - Single connection for MVP (no pooling needed)
        - Foreign keys enforced for referential integrity
    """
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    
    # Store schema version for future migrations
    conn.execute("""
        CREATE TABLE IF NOT EXISTS schema_metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        )
    """)
    
    conn.execute(
        "INSERT OR IGNORE INTO schema_metadata (key, value) VALUES (?, ?)",
        ("schema_version", str(SCHEMA_VERSION))
    )
    
    # Main events table with source linkage
    conn.execute("""
        CREATE TABLE IF NOT EXISTS events (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            source_file TEXT NOT NULL,
            line_offset INTEGER NOT NULL,
            checksum TEXT NOT NULL,
            dedupe_key TEXT UNIQUE NOT NULL,
            campaign_id TEXT NOT NULL,
            event_type TEXT NOT NULL,
            utc TEXT NOT NULL,
            payload JSON NOT NULL,
            ingested_at TEXT NOT NULL
        )
    """)
    
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_events_campaign "
        "ON events(campaign_id, event_type, utc)"
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_events_dedupe ON events(dedupe_key)"
    )
    
    # Ingestion errors (corruption/malformed JSON)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS ingestion_errors (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            source_file TEXT NOT NULL,
            line_offset INTEGER NOT NULL,
            error_type TEXT NOT NULL,
            error_message TEXT NOT NULL,
            raw_line TEXT,
            ingested_at TEXT NOT NULL
        )
    """)
    
    # Ingestion run log
    conn.execute("""
        CREATE TABLE IF NOT EXISTS ingestion_log (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            source_file TEXT NOT NULL,
            lines_read INTEGER NOT NULL,
            lines_ingested INTEGER NOT NULL,
            lines_skipped INTEGER NOT NULL,
            lines_errored INTEGER NOT NULL,
            started_at TEXT NOT NULL,
            completed_at TEXT NOT NULL
        )
    """)
    
    # Task summary records table (from task_summary.jsonl)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS task_summaries (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            source_file TEXT NOT NULL,
            line_offset INTEGER NOT NULL,
            checksum TEXT NOT NULL,
            dedupe_key TEXT UNIQUE NOT NULL,
            campaign_id TEXT NOT NULL,
            task_id TEXT NOT NULL,
            task_type TEXT,
            status TEXT,
            payload JSON NOT NULL,
            ingested_at TEXT NOT NULL
        )
    """)
    
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_task_summaries_campaign "
        "ON task_summaries(campaign_id, task_id, task_type)"
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_task_summaries_dedupe ON task_summaries(dedupe_key)"
    )
    
    # Orchestrator log entries (from orchestrator.log)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS orchestrator_log_entries (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            source_file TEXT NOT NULL,
            line_offset INTEGER NOT NULL,
            checksum TEXT NOT NULL,
            dedupe_key TEXT UNIQUE NOT NULL,
            campaign_id TEXT NOT NULL,
            log_level TEXT,
            log_timestamp TEXT,
            log_message TEXT NOT NULL,
            ingested_at TEXT NOT NULL
        )
    """)
    
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_orchestrator_log_campaign "
        "ON orchestrator_log_entries(campaign_id, log_level, log_timestamp)"
    )
    
    # Stdout/stderr metadata and source linkage
    conn.execute("""
        CREATE TABLE IF NOT EXISTS artifact_metadata (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            artifact_type TEXT NOT NULL,
            artifact_path TEXT NOT NULL,
            artifact_hash TEXT,
            related_campaign_id TEXT,
            related_task_id TEXT,
            source_file TEXT NOT NULL,
            line_offset INTEGER NOT NULL,
            dedupe_key TEXT UNIQUE NOT NULL,
            ingested_at TEXT NOT NULL
        )
    """)
    
    # Alerts state persistence table
    conn.execute("""
        CREATE TABLE IF NOT EXISTS alerts (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            campaign_id TEXT NOT NULL,
            fingerprint TEXT NOT NULL,
            alert_type TEXT NOT NULL,
            severity TEXT NOT NULL,
            message TEXT NOT NULL,
            context JSON NOT NULL,
            first_seen TEXT NOT NULL,
            last_seen TEXT NOT NULL,
            resolved_at TEXT,
            status TEXT NOT NULL DEFAULT 'active',
            UNIQUE(campaign_id, fingerprint)
        )
    """)
    
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_alerts_campaign_status "
        "ON alerts(campaign_id, status, last_seen)"
    )

    # Task 6 MVP registry/discovery/coverage persistence tables
    conn.execute("""
        CREATE TABLE IF NOT EXISTS run_registry (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            campaign_id TEXT NOT NULL,
            arm_id TEXT NOT NULL,
            run_dir_fingerprint TEXT NOT NULL,
            attempt_id TEXT NOT NULL,
            artifact_status TEXT NOT NULL,
            source_file TEXT NOT NULL,
            line_offset INTEGER NOT NULL,
            source_checksum TEXT NOT NULL,
            dedupe_key TEXT UNIQUE NOT NULL,
            payload JSON NOT NULL,
            updated_at TEXT NOT NULL,
            UNIQUE(campaign_id, arm_id, run_dir_fingerprint, attempt_id)
        )
    """)
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_run_registry_identity "
        "ON run_registry(campaign_id, arm_id, run_dir_fingerprint, attempt_id)"
    )

    conn.execute("""
        CREATE TABLE IF NOT EXISTS discovery_records (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            campaign_id TEXT NOT NULL,
            arm_id TEXT NOT NULL,
            run_dir_fingerprint TEXT NOT NULL,
            attempt_id TEXT NOT NULL,
            artifact_status TEXT NOT NULL,
            discovery_status TEXT NOT NULL,
            source_file TEXT NOT NULL,
            line_offset INTEGER NOT NULL,
            source_checksum TEXT NOT NULL,
            dedupe_key TEXT UNIQUE NOT NULL,
            payload JSON NOT NULL,
            updated_at TEXT NOT NULL,
            UNIQUE(campaign_id, arm_id, run_dir_fingerprint, attempt_id)
        )
    """)
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_discovery_records_identity "
        "ON discovery_records(campaign_id, arm_id, run_dir_fingerprint, attempt_id)"
    )

    conn.execute("""
        CREATE TABLE IF NOT EXISTS coverage_state (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            campaign_id TEXT NOT NULL,
            arm_id TEXT NOT NULL,
            run_dir_fingerprint TEXT NOT NULL,
            attempt_id TEXT NOT NULL,
            artifact_status TEXT NOT NULL,
            coverage_fraction REAL NOT NULL,
            axes_binned JSON NOT NULL,
            source_file TEXT NOT NULL,
            line_offset INTEGER NOT NULL,
            source_checksum TEXT NOT NULL,
            dedupe_key TEXT UNIQUE NOT NULL,
            payload JSON NOT NULL,
            updated_at TEXT NOT NULL,
            UNIQUE(campaign_id, arm_id, run_dir_fingerprint, attempt_id)
        )
    """)
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_coverage_state_identity "
        "ON coverage_state(campaign_id, arm_id, run_dir_fingerprint, attempt_id)"
    )

    conn.execute("""
        CREATE TABLE IF NOT EXISTS upsert_quarantine (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            source_file TEXT NOT NULL,
            line_offset INTEGER NOT NULL,
            source_checksum TEXT NOT NULL,
            dedupe_key TEXT UNIQUE NOT NULL,
            error_type TEXT NOT NULL,
            error_message TEXT NOT NULL,
            raw_payload JSON,
            quarantined_at TEXT NOT NULL
        )
    """)
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_upsert_quarantine_source "
        "ON upsert_quarantine(source_file, line_offset)"
    )

    # Reconciliation watermarks for incremental replay
    conn.execute("""
        CREATE TABLE IF NOT EXISTS reconcile_watermarks (
            source_file TEXT PRIMARY KEY,
            last_line_offset INTEGER NOT NULL,
            last_checksum TEXT NOT NULL,
            updated_at TEXT NOT NULL
        )
    """)
    
    conn.commit()
    return conn


def db_counts(conn: sqlite3.Connection) -> dict[str, Any]:
    """Get row counts for validation.
    
    Args:
        conn: SQLite connection
    
    Returns:
        Dictionary with row counts for each table
    """
    cursor = conn.cursor()
    
    def get_count(table: str) -> int:
        try:
            cursor.execute(f"SELECT COUNT(*) FROM {table}")
            return int(cursor.fetchone()[0])
        except sqlite3.OperationalError:
            return 0

    return {
        "events": get_count("events"),
        "ingestion_errors": get_count("ingestion_errors"),
        "ingestion_log": get_count("ingestion_log"),
        "task_summaries": get_count("task_summaries"),
        "orchestrator_log_entries": get_count("orchestrator_log_entries"),
        "artifact_metadata": get_count("artifact_metadata"),
        "alerts": get_count("alerts"),
        "run_registry": get_count("run_registry"),
        "discovery_records": get_count("discovery_records"),
        "coverage_state": get_count("coverage_state"),
        "upsert_quarantine": get_count("upsert_quarantine"),
        "reconcile_watermarks": get_count("reconcile_watermarks"),
    }
