"""Reconciliation watermark engine to track last processed offsets in telemetry sources.

Enables incremental replay by skipping already-processed lines in large JSONL artifacts.
Guarantees deterministic no-op behavior on second pass.
"""
from __future__ import annotations

import sqlite3
from collections.abc import Callable
from pathlib import Path
from typing import Any

def init_watermark_table(conn: sqlite3.Connection) -> None:
    """Initialize the reconciliation watermarks table."""
    conn.execute("""
        CREATE TABLE IF NOT EXISTS reconcile_watermarks (
            source_file TEXT PRIMARY KEY,
            last_line_offset INTEGER NOT NULL,
            last_checksum TEXT NOT NULL,
            updated_at TEXT NOT NULL
        )
    """)
    conn.commit()

def get_watermark(conn: sqlite3.Connection, source_file: str | Path) -> dict[str, Any] | None:
    """Get the last processed watermark for a source file."""
    cursor = conn.cursor()
    cursor.execute(
        "SELECT last_line_offset, last_checksum, updated_at FROM reconcile_watermarks WHERE source_file = ?",
        (str(source_file),)
    )
    row = cursor.fetchone()
    if row:
        return {
            "last_line_offset": row[0],
            "last_checksum": row[1],
            "updated_at": row[2]
        }
    return None

def update_watermark(
    conn: sqlite3.Connection, 
    source_file: str | Path, 
    line_offset: int, 
    checksum: str,
    utc_now_func: Callable[[], str],
) -> None:
    """Update or insert a watermark for a source file."""
    upsert_watermark(conn, source_file, line_offset, checksum, utc_now_func)
    conn.commit()


def upsert_watermark(
    conn: sqlite3.Connection,
    source_file: str | Path,
    line_offset: int,
    checksum: str,
    utc_now_func: Callable[[], str],
) -> None:
    """Update or insert a watermark without committing.

    This keeps watermark persistence available inside a broader transaction boundary
    where related domain upserts must commit atomically with the watermark delta.
    """
    conn.execute(
        """
        INSERT OR REPLACE INTO reconcile_watermarks (source_file, last_line_offset, last_checksum, updated_at)
        VALUES (?, ?, ?, ?)
        """,
        (str(source_file), line_offset, checksum, utc_now_func())
    )

class ReconcileStatus:
    """Summary of reconciliation activity."""
    def __init__(self, source_file: str):
        self.source_file = source_file
        self.new_count = 0
        self.unchanged_count = 0
        self.quarantined_count = 0
        self.error_count = 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "source_file": self.source_file,
            "new": self.new_count,
            "unchanged": self.unchanged_count,
            "quarantined": self.quarantined_count,
            "errors": self.error_count,
            "status": "success" if self.error_count == 0 else "partial_failure"
        }
