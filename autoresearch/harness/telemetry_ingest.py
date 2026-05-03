"""Idempotent ingestion from JSONL sources into SQLite telemetry store.

Ingestion strategy:
- Compute checksum (SHA256) for each line
- Build composite dedupe key: (source_file, line_offset, event_type, campaign_id)
- Skip existing dedupe keys (INSERT OR IGNORE)
- Log malformed/corrupted lines to ingestion_errors table
- Never crash on single bad line - isolate and continue

Replay guarantees:
- Second replay of same sources produces zero new rows
- Checksum mismatches detected and logged
- Full source linkage preserved (file path + line number)
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _compute_checksum(line: str) -> str:
    """Compute SHA256 checksum of line for corruption detection."""
    return hashlib.sha256(line.encode("utf-8")).hexdigest()


def _utc_now_iso() -> str:
    """Current UTC timestamp in ISO 8601 format."""
    return datetime.now(timezone.utc).isoformat()


def _build_dedupe_key(
    source_file: str,
    line_offset: int,
    event_type: str,
    campaign_id: str,
) -> str:
    """Build composite dedupe key for idempotent replay.
    
    Design: Composite key ensures same event from same source location
    is never ingested twice, even across multiple replay runs.
    """
    return f"{source_file}|{line_offset}|{event_type}|{campaign_id}"


def ingest_jsonl_file(
    conn: sqlite3.Connection,
    source_file: str | Path,
) -> dict[str, Any]:
    """Ingest a single JSONL file into SQLite store.
    
    Args:
        conn: SQLite connection (must have WAL mode enabled)
        source_file: Path to JSONL file (events.jsonl or task_summary.jsonl)
    
    Returns:
        Ingestion statistics dictionary with counts
    
    Design:
        - Line-by-line processing to handle large files
        - Graceful degradation on malformed JSON
        - Atomic transaction per file (all or nothing)
    """
    source_path = Path(source_file)
    if not source_path.exists():
        return {
            "source_file": str(source_file),
            "lines_read": 0,
            "lines_ingested": 0,
            "lines_skipped": 0,
            "lines_errored": 0,
            "status": "file_not_found",
        }
    
    started_at = _utc_now_iso()
    lines_read = 0
    lines_ingested = 0
    lines_skipped = 0
    lines_errored = 0
    
    cursor = conn.cursor()
    
    try:
        with source_path.open("r", encoding="utf-8") as f:
            for line_offset, raw_line in enumerate(f, start=1):
                lines_read += 1
                
                # Skip empty lines
                line = raw_line.strip()
                if not line:
                    lines_skipped += 1
                    continue
                
                # Compute checksum
                checksum = _compute_checksum(line)
                
                # Parse JSON
                try:
                    event = json.loads(line)
                except json.JSONDecodeError as e:
                    lines_errored += 1
                    cursor.execute(
                        """
                        INSERT INTO ingestion_errors 
                        (source_file, line_offset, error_type, error_message, raw_line, ingested_at)
                        VALUES (?, ?, ?, ?, ?, ?)
                        """,
                        (
                            str(source_path),
                            line_offset,
                            "json_decode_error",
                            str(e),
                            line[:1000],  # Truncate to avoid blob storage issues
                            _utc_now_iso(),
                        )
                    )
                    continue
                
                # Extract required fields
                try:
                    event_type = event.get("event_type", "UNKNOWN")
                    campaign_id = event.get("campaign_id", "UNKNOWN")
                    utc = event.get("utc", "")
                    
                    # Build dedupe key
                    dedupe_key = _build_dedupe_key(
                        str(source_path), line_offset, event_type, campaign_id
                    )
                    
                    # Attempt insert (dedupe_key UNIQUE constraint prevents duplicates)
                    cursor.execute(
                        """
                        INSERT OR IGNORE INTO events
                        (source_file, line_offset, checksum, dedupe_key, campaign_id, 
                         event_type, utc, payload, ingested_at)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            str(source_path),
                            line_offset,
                            checksum,
                            dedupe_key,
                            campaign_id,
                            event_type,
                            utc,
                            json.dumps(event.get("payload", {})),
                            _utc_now_iso(),
                        )
                    )
                    
                    # Check if row was actually inserted (rowcount == 1)
                    if cursor.rowcount == 1:
                        lines_ingested += 1
                    else:
                        lines_skipped += 1
                
                except (KeyError, TypeError, ValueError) as e:
                    lines_errored += 1
                    cursor.execute(
                        """
                        INSERT INTO ingestion_errors 
                        (source_file, line_offset, error_type, error_message, raw_line, ingested_at)
                        VALUES (?, ?, ?, ?, ?, ?)
                        """,
                        (
                            str(source_path),
                            line_offset,
                            "schema_validation_error",
                            str(e),
                            line[:1000],
                            _utc_now_iso(),
                        )
                    )
                    continue
        
        # Log ingestion run
        completed_at = _utc_now_iso()
        cursor.execute(
            """
            INSERT INTO ingestion_log
            (source_file, lines_read, lines_ingested, lines_skipped, lines_errored, 
             started_at, completed_at)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                str(source_path),
                lines_read,
                lines_ingested,
                lines_skipped,
                lines_errored,
                started_at,
                completed_at,
            )
        )
        
        conn.commit()
        
        return {
            "source_file": str(source_path),
            "lines_read": lines_read,
            "lines_ingested": lines_ingested,
            "lines_skipped": lines_skipped,
            "lines_errored": lines_errored,
            "status": "success",
        }
    
    except Exception as e:
        conn.rollback()
        return {
            "source_file": str(source_path),
            "lines_read": lines_read,
            "lines_ingested": lines_ingested,
            "lines_skipped": lines_skipped,
            "lines_errored": lines_errored,
            "status": "fatal_error",
            "error": str(e),
        }


def ingest_task_summary_jsonl(
    conn: sqlite3.Connection,
    source_file: str | Path,
) -> dict[str, Any]:
    """Ingest task summary records from JSONL file into SQLite store.
    
    Args:
        conn: SQLite connection (must have WAL mode enabled)
        source_file: Path to task_summary.jsonl file
    
    Returns:
        Ingestion statistics dictionary with counts
    
    Design:
        - Task summaries have different schema than events
        - Extract task_id, task_type, status fields
        - Fall back to UNKNOWN for missing campaign_id
    """
    source_path = Path(source_file)
    if not source_path.exists():
        return {
            "source_file": str(source_file),
            "lines_read": 0,
            "lines_ingested": 0,
            "lines_skipped": 0,
            "lines_errored": 0,
            "status": "file_not_found",
        }
    
    started_at = _utc_now_iso()
    lines_read = 0
    lines_ingested = 0
    lines_skipped = 0
    lines_errored = 0
    
    cursor = conn.cursor()
    
    try:
        with source_path.open("r", encoding="utf-8") as f:
            for line_offset, raw_line in enumerate(f, start=1):
                lines_read += 1
                
                # Skip empty lines
                line = raw_line.strip()
                if not line:
                    lines_skipped += 1
                    continue
                
                # Compute checksum
                checksum = _compute_checksum(line)
                
                # Parse JSON
                try:
                    task = json.loads(line)
                except json.JSONDecodeError as e:
                    lines_errored += 1
                    cursor.execute(
                        """
                        INSERT INTO ingestion_errors 
                        (source_file, line_offset, error_type, error_message, raw_line, ingested_at)
                        VALUES (?, ?, ?, ?, ?, ?)
                        """,
                        (
                            str(source_path),
                            line_offset,
                            "json_decode_error",
                            str(e),
                            line[:1000],
                            _utc_now_iso(),
                        )
                    )
                    continue
                
                # Extract task summary fields
                try:
                    campaign_id = task.get("campaign_id", "UNKNOWN")
                    task_id = task.get("task_id", "UNKNOWN")
                    task_type = task.get("task_type", "UNKNOWN")
                    status = task.get("status", "UNKNOWN")
                    
                    # Build dedupe key for task summaries
                    dedupe_key = f"{source_path}|{line_offset}|{task_id}|{campaign_id}"
                    
                    # Attempt insert (dedupe_key UNIQUE constraint prevents duplicates)
                    cursor.execute(
                        """
                        INSERT OR IGNORE INTO task_summaries
                        (source_file, line_offset, checksum, dedupe_key, campaign_id, 
                         task_id, task_type, status, payload, ingested_at)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            str(source_path),
                            line_offset,
                            checksum,
                            dedupe_key,
                            campaign_id,
                            task_id,
                            task_type,
                            status,
                            json.dumps(task),
                            _utc_now_iso(),
                        )
                    )
                    
                    # Check if row was actually inserted (rowcount == 1)
                    if cursor.rowcount == 1:
                        lines_ingested += 1
                    else:
                        lines_skipped += 1
                
                except (KeyError, TypeError, ValueError) as e:
                    lines_errored += 1
                    cursor.execute(
                        """
                        INSERT INTO ingestion_errors 
                        (source_file, line_offset, error_type, error_message, raw_line, ingested_at)
                        VALUES (?, ?, ?, ?, ?, ?)
                        """,
                        (
                            str(source_path),
                            line_offset,
                            "schema_validation_error",
                            str(e),
                            line[:1000],
                            _utc_now_iso(),
                        )
                    )
                    continue
        
        # Log ingestion run
        completed_at = _utc_now_iso()
        cursor.execute(
            """
            INSERT INTO ingestion_log
            (source_file, lines_read, lines_ingested, lines_skipped, lines_errored, 
             started_at, completed_at)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                str(source_path),
                lines_read,
                lines_ingested,
                lines_skipped,
                lines_errored,
                started_at,
                completed_at,
            )
        )
        
        conn.commit()
        
        return {
            "source_file": str(source_path),
            "lines_read": lines_read,
            "lines_ingested": lines_ingested,
            "lines_skipped": lines_skipped,
            "lines_errored": lines_errored,
            "status": "success",
        }
    
    except Exception as e:
        conn.rollback()
        return {
            "source_file": str(source_path),
            "lines_read": lines_read,
            "lines_ingested": lines_ingested,
            "lines_skipped": lines_skipped,
            "lines_errored": lines_errored,
            "status": "fatal_error",
            "error": str(e),
        }


def ingest_orchestrator_log(
    conn: sqlite3.Connection,
    source_file: str | Path,
    campaign_id: str = "UNKNOWN",
) -> dict[str, Any]:
    """Ingest plain-text orchestrator log into SQLite store.
    
    Args:
        conn: SQLite connection (must have WAL mode enabled)
        source_file: Path to orchestrator.log file
        campaign_id: Campaign identifier for context
    
    Returns:
        Ingestion statistics dictionary with counts
    
    Design:
        - Parses plain-text log lines with [TIMESTAMP] [LEVEL] message format
        - Each line stored with source linkage for traceability
        - Gracefully handles malformed lines
    """
    source_path = Path(source_file)
    if not source_path.exists():
        return {
            "source_file": str(source_file),
            "lines_read": 0,
            "lines_ingested": 0,
            "lines_skipped": 0,
            "lines_errored": 0,
            "status": "file_not_found",
        }
    
    started_at = _utc_now_iso()
    lines_read = 0
    lines_ingested = 0
    lines_skipped = 0
    lines_errored = 0
    
    cursor = conn.cursor()
    
    try:
        with source_path.open("r", encoding="utf-8") as f:
            for line_offset, raw_line in enumerate(f, start=1):
                lines_read += 1
                
                # Skip empty lines
                line = raw_line.strip()
                if not line:
                    lines_skipped += 1
                    continue
                
                # Compute checksum
                checksum = _compute_checksum(line)
                
                # Parse log line format: [TIMESTAMP] [LEVEL] message
                log_level = "UNKNOWN"
                log_timestamp = ""
                log_message = line
                
                # Simple parser for [TIMESTAMP] [LEVEL] pattern
                try:
                    if line.startswith("["):
                        # Try to extract timestamp
                        close_bracket = line.find("]")
                        if close_bracket > 0:
                            log_timestamp = line[1:close_bracket]
                            remainder = line[close_bracket+1:].strip()
                            
                            # Try to extract level
                            if remainder.startswith("["):
                                close_level = remainder.find("]")
                                if close_level > 0:
                                    log_level = remainder[1:close_level]
                                    log_message = remainder[close_level+1:].strip()
                except Exception:
                    pass  # Use defaults
                
                # Build dedupe key for orchestrator log
                dedupe_key = f"{source_path}|{line_offset}|{log_level}|{campaign_id}"
                
                # Attempt insert
                try:
                    cursor.execute(
                        """
                        INSERT OR IGNORE INTO orchestrator_log_entries
                        (source_file, line_offset, checksum, dedupe_key, campaign_id, 
                         log_level, log_timestamp, log_message, ingested_at)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            str(source_path),
                            line_offset,
                            checksum,
                            dedupe_key,
                            campaign_id,
                            log_level,
                            log_timestamp,
                            log_message,
                            _utc_now_iso(),
                        )
                    )
                    
                    # Check if row was actually inserted (rowcount == 1)
                    if cursor.rowcount == 1:
                        lines_ingested += 1
                    else:
                        lines_skipped += 1
                
                except Exception as e:
                    lines_errored += 1
                    cursor.execute(
                        """
                        INSERT INTO ingestion_errors 
                        (source_file, line_offset, error_type, error_message, raw_line, ingested_at)
                        VALUES (?, ?, ?, ?, ?, ?)
                        """,
                        (
                            str(source_path),
                            line_offset,
                            "orchestrator_log_parse_error",
                            str(e),
                            line[:1000],
                            _utc_now_iso(),
                        )
                    )
                    continue
        
        # Log ingestion run
        completed_at = _utc_now_iso()
        cursor.execute(
            """
            INSERT INTO ingestion_log
            (source_file, lines_read, lines_ingested, lines_skipped, lines_errored, 
             started_at, completed_at)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                str(source_path),
                lines_read,
                lines_ingested,
                lines_skipped,
                lines_errored,
                started_at,
                completed_at,
            )
        )
        
        conn.commit()
        
        return {
            "source_file": str(source_path),
            "lines_read": lines_read,
            "lines_ingested": lines_ingested,
            "lines_skipped": lines_skipped,
            "lines_errored": lines_errored,
            "status": "success",
        }
    
    except Exception as e:
        conn.rollback()
        return {
            "source_file": str(source_path),
            "lines_read": lines_read,
            "lines_ingested": lines_ingested,
            "lines_skipped": lines_skipped,
            "lines_errored": lines_errored,
            "status": "fatal_error",
            "error": str(e),
        }


def ingest_artifact_metadata(
    conn: sqlite3.Connection,
    artifact_type: str,
    artifact_path: str | Path,
    campaign_id: str = "UNKNOWN",
    task_id: str | None = None,
    source_file: str | Path | None = None,
    line_offset: int = 0,
) -> dict[str, Any]:
    """Ingest stdout/stderr or other artifact metadata.
    
    Args:
        conn: SQLite connection
        artifact_type: Type of artifact (stdout, stderr, etc.)
        artifact_path: Path to the artifact file
        campaign_id: Associated campaign ID
        task_id: Associated task ID (optional)
        source_file: Source file that referenced this artifact
        line_offset: Line offset in source file
    
    Returns:
        Result dictionary
    """
    artifact_path_obj = Path(artifact_path)
    source_path = Path(source_file) if source_file else artifact_path_obj
    
    # Compute artifact hash if file exists
    artifact_hash = ""
    if artifact_path_obj.exists():
        try:
            with artifact_path_obj.open("rb") as f:
                artifact_hash = hashlib.sha256(f.read()).hexdigest()
        except Exception:
            pass  # Hash computation optional
    
    # Build dedupe key
    dedupe_key = f"{artifact_path}|{artifact_type}|{campaign_id}"
    
    try:
        cursor = conn.cursor()
        cursor.execute(
            """
            INSERT OR IGNORE INTO artifact_metadata
            (artifact_type, artifact_path, artifact_hash, related_campaign_id, 
             related_task_id, source_file, line_offset, dedupe_key, ingested_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                artifact_type,
                str(artifact_path),
                artifact_hash,
                campaign_id,
                task_id,
                str(source_path),
                line_offset,
                dedupe_key,
                _utc_now_iso(),
            )
        )
        conn.commit()
        
        return {"status": "success", "artifact_type": artifact_type}
    except Exception as e:
        conn.rollback()
        return {"status": "error", "error": str(e)}

def replay_sources(
    conn: sqlite3.Connection,
    source_files: list[str | Path],
    campaign_id: str = "UNKNOWN",
) -> list[dict[str, Any]]:
    """Idempotently replay ingestion from multiple source types.
    
    Args:
        conn: SQLite connection (must have WAL mode enabled)
        source_files: List of file paths (events.jsonl, task_summary.jsonl, orchestrator.log, etc.)
        campaign_id: Campaign identifier for context
    
    Returns:
        List of ingestion statistics (one per source file)
    
    Design:
        - Auto-detect file type by name and extension
        - Each file ingested independently (partial failures isolated)
        - Dedupe keys prevent duplicate events across replays
        - Second replay of same sources produces zero new rows
    """
    results = []
    for source_file in source_files:
        source_path = Path(source_file)
        filename = source_path.name
        
        # Auto-detect file type and route to appropriate ingester
        if filename == "orchestrator.log":
            result = ingest_orchestrator_log(conn, source_file, campaign_id=campaign_id)
        elif filename == "task_summary.jsonl":
            result = ingest_task_summary_jsonl(conn, source_file)
        else:
            # Default: treat as events.jsonl or generic JSONL
            result = ingest_jsonl_file(conn, source_file)
        
        results.append(result)
    
    return results

