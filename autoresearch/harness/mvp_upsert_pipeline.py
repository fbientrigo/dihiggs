"""Task 6 atomic upsert pipeline for registry/discovery/coverage state."""
from __future__ import annotations

import hashlib
import json
import sqlite3
from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, cast

from autoresearch.harness.coverage_contract_bridge import extend_coverage_state
from autoresearch.harness.reconcile_watermark import get_watermark, upsert_watermark
from autoresearch.harness.run_identity_contract import make_run_identity_tuple

from autoresearch.harness.mvp_graph import (
    GraphReconcileCursor,
    MvpGraphState,
    apply_reconciled_graph_delta,
)

CANONICAL_IDENTITY_FIELDS = (
    "campaign_id",
    "arm_id",
    "run_dir_fingerprint",
    "attempt_id",
)


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _stable_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _row_checksum(row: object) -> str:
    return hashlib.sha256(_stable_json(row).encode("utf-8")).hexdigest()


def _identity_dedupe_key(record: Mapping[str, object]) -> str:
    identity = make_run_identity_tuple(
        str(record["campaign_id"]),
        str(record["arm_id"]),
        str(record["run_dir_fingerprint"]),
        str(record["attempt_id"]),
    )
    return "|".join(identity.as_tuple())


def _quarantine_dedupe_key(source_file: str, line_offset: int, checksum: str) -> str:
    return f"{source_file}|{line_offset}|{checksum}"


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return cast(Mapping[str, object], value)


def _mapping_or_empty(value: object) -> dict[str, object]:
    if isinstance(value, Mapping):
        return dict(cast(Mapping[str, object], value))
    return {}


def _artifact_status(record: Mapping[str, object]) -> str:
    explicit = record.get("artifact_status")
    if explicit not in (None, ""):
        return str(explicit)
    parser_status = str(record.get("parser_status", "")).strip()
    if parser_status and parser_status != "success":
        return "partial"
    return "complete"


def _discovery_status(record: Mapping[str, object], artifact_status: str) -> str:
    explicit = record.get("discovery_status")
    if explicit not in (None, ""):
        return str(explicit)
    parser_status = str(record.get("parser_status", "")).strip()
    if parser_status.startswith("quarantined"):
        return "quarantined"
    if parser_status == "missing" or artifact_status == "partial":
        return "partial"
    return "new"


def _normalize_scanner_row(raw_row: object) -> dict[str, object]:
    record = dict(_mapping(raw_row, "scanner row"))
    for field in CANONICAL_IDENTITY_FIELDS:
        value = record.get(field)
        if value in (None, ""):
            raise ValueError(f"missing required identity field: {field}")
        record[field] = str(value)

    record["artifact_status"] = _artifact_status(record)
    record["discovery_status"] = _discovery_status(record, str(record["artifact_status"]))
    record["axes_binned"] = _mapping_or_empty(record.get("axes_binned"))
    return record


def _fetch_existing_coverage_state(conn: sqlite3.Connection) -> list[dict[str, object]]:
    cursor = conn.cursor()
    cursor.execute(
        "SELECT payload FROM coverage_state ORDER BY source_file, line_offset, id"
    )
    rows: list[dict[str, object]] = []
    for payload_row in cursor.fetchall():
        payload_text = payload_row[0]
        if payload_text in (None, ""):
            continue
        rows.append(cast(dict[str, object], json.loads(str(payload_text))))
    return rows


def _float_value(value: object, *, field_name: str) -> float:
    if isinstance(value, (int, float)):
        return float(value)
    raise ValueError(f"{field_name} must be numeric")


def _upsert_quarantine_row(
    cursor: sqlite3.Cursor,
    *,
    source_file: str,
    line_offset: int,
    checksum: str,
    raw_row: object,
    error_type: str,
    error_message: str,
    utc_now_func: Callable[[], str],
) -> None:
    cursor.execute(
        """
        INSERT OR IGNORE INTO upsert_quarantine
        (source_file, line_offset, source_checksum, dedupe_key, error_type, error_message, raw_payload, quarantined_at)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            source_file,
            line_offset,
            checksum,
            _quarantine_dedupe_key(source_file, line_offset, checksum),
            error_type,
            error_message,
            _stable_json(raw_row),
            utc_now_func(),
        ),
    )


def _upsert_domain_row(
    cursor: sqlite3.Cursor,
    *,
    table: str,
    payload: Mapping[str, object],
    source_file: str,
    line_offset: int,
    checksum: str,
    utc_now_func: Callable[[], str],
) -> str:
    dedupe_key = _identity_dedupe_key(payload)
    payload_json = _stable_json(payload)
    cursor.execute(f"SELECT payload FROM {table} WHERE dedupe_key = ?", (dedupe_key,))
    row = cursor.fetchone()
    if row is None:
        if table == "run_registry":
            cursor.execute(
                """
                INSERT INTO run_registry
                (campaign_id, arm_id, run_dir_fingerprint, attempt_id, artifact_status,
                 source_file, line_offset, source_checksum, dedupe_key, payload, updated_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    payload["campaign_id"],
                    payload["arm_id"],
                    payload["run_dir_fingerprint"],
                    payload["attempt_id"],
                    payload["artifact_status"],
                    source_file,
                    line_offset,
                    checksum,
                    dedupe_key,
                    payload_json,
                    utc_now_func(),
                ),
            )
        elif table == "discovery_records":
            cursor.execute(
                """
                INSERT INTO discovery_records
                (campaign_id, arm_id, run_dir_fingerprint, attempt_id, artifact_status, discovery_status,
                 source_file, line_offset, source_checksum, dedupe_key, payload, updated_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    payload["campaign_id"],
                    payload["arm_id"],
                    payload["run_dir_fingerprint"],
                    payload["attempt_id"],
                    payload["artifact_status"],
                    payload["discovery_status"],
                    source_file,
                    line_offset,
                    checksum,
                    dedupe_key,
                    payload_json,
                    utc_now_func(),
                ),
            )
        elif table == "coverage_state":
            cursor.execute(
                """
                INSERT INTO coverage_state
                (campaign_id, arm_id, run_dir_fingerprint, attempt_id, artifact_status, coverage_fraction, axes_binned,
                 source_file, line_offset, source_checksum, dedupe_key, payload, updated_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    payload["campaign_id"],
                    payload["arm_id"],
                    payload["run_dir_fingerprint"],
                    payload["attempt_id"],
                    payload["artifact_status"],
                    _float_value(payload["coverage_fraction"], field_name="coverage_fraction"),
                    _stable_json(payload.get("axes_binned", {})),
                    source_file,
                    line_offset,
                    checksum,
                    dedupe_key,
                    payload_json,
                    utc_now_func(),
                ),
            )
        else:
            raise ValueError(f"unsupported table: {table}")
        return "insert"

    existing_payload = str(row[0])
    if existing_payload == payload_json:
        return "noop"

    if table == "run_registry":
        cursor.execute(
            """
            UPDATE run_registry
            SET artifact_status = ?, source_file = ?, line_offset = ?, source_checksum = ?, payload = ?, updated_at = ?
            WHERE dedupe_key = ?
            """,
            (payload["artifact_status"], source_file, line_offset, checksum, payload_json, utc_now_func(), dedupe_key),
        )
    elif table == "discovery_records":
        cursor.execute(
            """
            UPDATE discovery_records
            SET artifact_status = ?, discovery_status = ?, source_file = ?, line_offset = ?,
                source_checksum = ?, payload = ?, updated_at = ?
            WHERE dedupe_key = ?
            """,
            (
                payload["artifact_status"],
                payload["discovery_status"],
                source_file,
                line_offset,
                checksum,
                payload_json,
                utc_now_func(),
                dedupe_key,
            ),
        )
    elif table == "coverage_state":
        cursor.execute(
            """
            UPDATE coverage_state
            SET artifact_status = ?, coverage_fraction = ?, axes_binned = ?, source_file = ?, line_offset = ?,
                source_checksum = ?, payload = ?, updated_at = ?
            WHERE dedupe_key = ?
            """,
            (
                payload["artifact_status"],
                _float_value(payload["coverage_fraction"], field_name="coverage_fraction"),
                _stable_json(payload.get("axes_binned", {})),
                source_file,
                line_offset,
                checksum,
                payload_json,
                utc_now_func(),
                dedupe_key,
            ),
        )
    else:
        raise ValueError(f"unsupported table: {table}")
    return "update"


def execute_registry_discovery_coverage_upsert(
    conn: sqlite3.Connection,
    *,
    source_file: str | Path,
    scanner_rows: Sequence[object],
    config: Mapping[str, object],
    utc_now_func: Callable[[], str] = _utc_now_iso,
) -> dict[str, Any]:
    """Persist Task 6 registry, discovery, and coverage-state deltas atomically.

    Valid rows continue through the pipeline even when neighboring rows are malformed.
    Watermark persistence happens inside the same transaction boundary as all domain writes.
    """

    source_path = str(Path(source_file))
    watermark = get_watermark(conn, source_path)
    last_offset = int(watermark["last_line_offset"]) if watermark is not None else 0
    last_checksum = str(watermark["last_checksum"]) if watermark is not None else ""
    last_processed: tuple[int, str] | None = None

    existing_coverage_state = _fetch_existing_coverage_state(conn)
    valid_records: list[dict[str, object]] = []
    record_meta: list[tuple[int, str]] = []
    quarantined = 0
    skipped_by_watermark = 0

    for line_offset, raw_row in enumerate(scanner_rows, start=1):
        checksum = _row_checksum(raw_row)
        if line_offset < last_offset or (line_offset == last_offset and checksum == last_checksum):
            skipped_by_watermark += 1
            continue

        try:
            valid_records.append(_normalize_scanner_row(raw_row))
            record_meta.append((line_offset, checksum))
        except ValueError as exc:
            quarantined += 1
            conn.execute("BEGIN")
            try:
                _upsert_quarantine_row(
                    conn.cursor(),
                    source_file=source_path,
                    line_offset=line_offset,
                    checksum=checksum,
                    raw_row=raw_row,
                    error_type="schema_validation_error",
                    error_message=str(exc),
                    utc_now_func=utc_now_func,
                )
                upsert_watermark(conn, source_path, line_offset, checksum, utc_now_func)
                conn.commit()
            except Exception:
                conn.rollback()
                raise
            last_processed = (line_offset, checksum)

    coverage_rows = extend_coverage_state(existing_coverage_state, valid_records, config)
    ops = {
        "run_registry": {"insert": 0, "update": 0, "noop": 0},
        "discovery_records": {"insert": 0, "update": 0, "noop": 0},
        "coverage_state": {"insert": 0, "update": 0, "noop": 0},
    }

    if valid_records:
        conn.execute("BEGIN")
        try:
            cursor = conn.cursor()
            for record, coverage_row, (line_offset, checksum) in zip(valid_records, coverage_rows, record_meta, strict=True):
                run_payload = {
                    key: record[key]
                    for key in (*CANONICAL_IDENTITY_FIELDS, "artifact_status")
                }
                discovery_payload = {
                    key: record[key]
                    for key in (*CANONICAL_IDENTITY_FIELDS, "artifact_status", "discovery_status")
                }
                coverage_payload = dict(coverage_row)

                run_op = _upsert_domain_row(
                    cursor,
                    table="run_registry",
                    payload=run_payload,
                    source_file=source_path,
                    line_offset=line_offset,
                    checksum=checksum,
                    utc_now_func=utc_now_func,
                )
                discovery_op = _upsert_domain_row(
                    cursor,
                    table="discovery_records",
                    payload=discovery_payload,
                    source_file=source_path,
                    line_offset=line_offset,
                    checksum=checksum,
                    utc_now_func=utc_now_func,
                )
                coverage_op = _upsert_domain_row(
                    cursor,
                    table="coverage_state",
                    payload=coverage_payload,
                    source_file=source_path,
                    line_offset=line_offset,
                    checksum=checksum,
                    utc_now_func=utc_now_func,
                )
                ops["run_registry"][run_op] += 1
                ops["discovery_records"][discovery_op] += 1
                ops["coverage_state"][coverage_op] += 1
                last_processed = (line_offset, checksum)

            if last_processed is not None:
                upsert_watermark(conn, source_path, last_processed[0], last_processed[1], utc_now_func)
            conn.commit()
        except Exception:
            conn.rollback()
            raise

    return {
        "status": "success",
        "source_file": source_path,
        "rows_seen": len(scanner_rows),
        "rows_valid": len(valid_records),
        "rows_quarantined": quarantined,
        "rows_skipped_by_watermark": skipped_by_watermark,
        "run_registry": ops["run_registry"],
        "discovery_records": ops["discovery_records"],
        "coverage_state": ops["coverage_state"],
        "watermark": get_watermark(conn, source_path),
    }

def trigger_incremental_update(
    conn: sqlite3.Connection,
    graph: MvpGraphState,
    *,
    source_file: str | Path,
    scanner_rows: Sequence[object],
    config: Mapping[str, object],
    applied_cursor: GraphReconcileCursor | None = None,
    utc_now_func: Callable[[], str] = _utc_now_iso,
) -> dict[str, Any]:
    """Quick incremental update trigger for new run arrivals.

    Processes only new/changed artifacts, updates registry + coverage + graph,
    and returns a deterministic no-op summary when no changes are detected.
    """
    upsert_result = execute_registry_discovery_coverage_upsert(
        conn,
        source_file=source_file,
        scanner_rows=scanner_rows,
        config=config,
        utc_now_func=utc_now_func,
    )

    graph_result = apply_reconciled_graph_delta(
        graph,
        conn,
        upsert_result,
        applied_cursor=applied_cursor,
    )

    is_noop = (
        upsert_result["rows_valid"] == 0
        and upsert_result["rows_quarantined"] == 0
        and graph_result.summary.inserted == 0
        and graph_result.summary.updated == 0
    )

    return {
        "status": "success",
        "is_noop": is_noop,
        "upsert": upsert_result,
        "graph": {
            "summary": {
                "inserted": graph_result.summary.inserted,
                "updated": graph_result.summary.updated,
                "noops": graph_result.summary.noops,
                "node_count": graph_result.summary.node_count,
                "edge_count": graph_result.summary.edge_count,
            },
            "cursor": graph_result.cursor,
        },
        "utc": utc_now_func(),
    }

