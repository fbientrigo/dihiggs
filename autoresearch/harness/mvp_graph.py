from __future__ import annotations

import json
import sqlite3
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import cast

from autoresearch.harness.run_identity_contract import validate_mvp_contract_artifact

CANONICAL_IDENTITY_FIELDS = (
    "campaign_id",
    "arm_id",
    "run_dir_fingerprint",
    "attempt_id",
)
SUPPORTED_NODE_ENTITY_TYPES = ("run", "discovery", "coverage_state")
SUPPORTED_EDGE_ENTITY_TYPES = ("run_lineage", "discovery_provenance", "coverage_state_link")
SUPPORTED_GRAPH_ENTITY_TYPES = SUPPORTED_NODE_ENTITY_TYPES + SUPPORTED_EDGE_ENTITY_TYPES


class GraphContractError(ValueError):
    """Raised when MVP graph deltas fall outside the supported contract."""


@dataclass(frozen=True)
class GraphDelta:
    entity_type: str
    identity: tuple[str, ...]
    record: dict[str, object]


@dataclass(frozen=True)
class GraphApplyResult:
    entity_type: str
    identity: tuple[str, ...]
    operation: str


@dataclass(frozen=True)
class GraphApplySummary:
    inserted: int
    updated: int
    noops: int
    node_count: int
    edge_count: int


@dataclass(frozen=True)
class GraphReconcileCursor:
    source_file: str
    last_line_offset: int
    last_checksum: str


@dataclass(frozen=True)
class ReconciledGraphDeltaResult:
    deltas: tuple[GraphDelta, ...]
    summary: GraphApplySummary
    cursor: GraphReconcileCursor | None


@dataclass
class MvpGraphState:
    nodes: dict[tuple[str, ...], dict[str, object]] = field(default_factory=dict)
    edges: dict[tuple[str, ...], dict[str, object]] = field(default_factory=dict)


def _require_mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise GraphContractError(f"{name} must be a mapping")
    return cast(Mapping[str, object], value)


def _require_identity(record: Mapping[str, object], fields: Sequence[str], *, entity_type: str) -> tuple[str, ...]:
    identity: list[str] = []
    for field_name in fields:
        value = record.get(field_name)
        if value in (None, ""):
            raise GraphContractError(f"{entity_type} delta missing required identity field: {field_name}")
        identity.append(str(value))
    return tuple(identity)


def _supported_types_message() -> str:
    return ", ".join(SUPPORTED_GRAPH_ENTITY_TYPES)


def _coerce_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise GraphContractError(f"{field_name} must be an integer")
    return value


def _cursor_from_upsert_result(upsert_result: Mapping[str, object]) -> GraphReconcileCursor | None:
    source_file = upsert_result.get("source_file")
    if source_file in (None, ""):
        raise GraphContractError("upsert result missing required field: source_file")
    watermark = _require_mapping(upsert_result.get("watermark"), "upsert result watermark")
    last_line_offset = watermark.get("last_line_offset")
    last_checksum = watermark.get("last_checksum")
    if last_line_offset in (None, "") or last_checksum in (None, ""):
        return None
    return GraphReconcileCursor(
        source_file=str(source_file),
        last_line_offset=_coerce_int(last_line_offset, field_name="watermark.last_line_offset"),
        last_checksum=str(last_checksum),
    )


def _decode_payload(row: Sequence[object], *, table: str) -> dict[str, object]:
    payload_text = row[1]
    if not isinstance(payload_text, str) or payload_text == "":
        raise GraphContractError(f"{table} row payload must be non-empty JSON text")
    payload = cast(object, json.loads(payload_text))
    if not isinstance(payload, Mapping):
        raise GraphContractError(f"{table} row payload must decode to a mapping")
    record = dict(cast(Mapping[str, object], payload))
    record["source_file"] = str(row[2])
    record["line_offset"] = _coerce_int(row[3], field_name=f"{table}.line_offset")
    record["source_checksum"] = str(row[4])
    record["updated_at"] = str(row[5])
    record["source_table"] = table
    return record


def _fetch_domain_rows(
    conn: sqlite3.Connection,
    *,
    table: str,
    source_file: str,
    applied_cursor: GraphReconcileCursor | None,
    target_cursor: GraphReconcileCursor,
) -> list[dict[str, object]]:
    start_line_offset = 0
    if applied_cursor is not None:
        if applied_cursor.source_file != source_file:
            raise GraphContractError("graph cursor source_file does not match upsert result source_file")
        start_line_offset = applied_cursor.last_line_offset

    cursor = conn.cursor()
    _ = cursor.execute(
        f"""
        SELECT id, payload, source_file, line_offset, source_checksum, updated_at
        FROM {table}
        WHERE source_file = ? AND line_offset > ? AND line_offset <= ?
        ORDER BY line_offset, id
        """,
        (source_file, start_line_offset, target_cursor.last_line_offset),
    )
    fetched_rows = cast(list[Sequence[object]], cursor.fetchall())
    return [_decode_payload(row, table=table) for row in fetched_rows]


def _lineage_edge_record(source_run: Mapping[str, object], target_run: Mapping[str, object]) -> dict[str, object]:
    return {
        "campaign_id": str(target_run["campaign_id"]),
        "arm_id": str(target_run["arm_id"]),
        "source_attempt_id": str(source_run["attempt_id"]),
        "target_attempt_id": str(target_run["attempt_id"]),
        "edge_type": "run_lineage",
        "source_entity_type": "run",
        "target_entity_type": "run",
        "source_run_dir_fingerprint": str(source_run["run_dir_fingerprint"]),
        "target_run_dir_fingerprint": str(target_run["run_dir_fingerprint"]),
        "source_file": str(target_run["source_file"]),
        "line_offset": _coerce_int(target_run["line_offset"], field_name="run_lineage.line_offset"),
        "source_checksum": str(target_run["source_checksum"]),
        "updated_at": str(target_run["updated_at"]),
        "source_table": "run_registry",
    }


def _fetch_previous_run_row(
    conn: sqlite3.Connection,
    *,
    campaign_id: str,
    arm_id: str,
    current_line_offset: int,
) -> dict[str, object] | None:
    cursor = conn.cursor()
    _ = cursor.execute(
        """
        SELECT id, payload, source_file, line_offset, source_checksum, updated_at
        FROM run_registry
        WHERE campaign_id = ? AND arm_id = ? AND line_offset < ?
        ORDER BY line_offset DESC, id DESC
        LIMIT 1
        """,
        (campaign_id, arm_id, current_line_offset),
    )
    row = cast(Sequence[object] | None, cursor.fetchone())
    if row is None:
        return None
    return _decode_payload(row, table="run_registry")


def _run_lineage_deltas_from_rows(
    conn: sqlite3.Connection,
    run_rows: Sequence[Mapping[str, object]],
) -> list[GraphDelta]:
    deltas: list[GraphDelta] = []
    previous_by_scope: dict[tuple[str, str], Mapping[str, object]] = {}
    for run_row in run_rows:
        scope = (str(run_row["campaign_id"]), str(run_row["arm_id"]))
        source_run = previous_by_scope.get(scope)
        if source_run is None:
            source_run = _fetch_previous_run_row(
                conn,
                campaign_id=scope[0],
                arm_id=scope[1],
                current_line_offset=_coerce_int(run_row["line_offset"], field_name="run.line_offset"),
            )
        previous_by_scope[scope] = run_row
        if source_run is None:
            continue

        record = _lineage_edge_record(source_run, run_row)
        identity = _require_identity(
            record,
            ("campaign_id", "arm_id", "source_attempt_id", "target_attempt_id"),
            entity_type="run_lineage",
        )
        deltas.append(GraphDelta(entity_type="run_lineage", identity=identity, record=record))
    return deltas


def graph_deltas_from_reconciled_state(
    conn: sqlite3.Connection,
    upsert_result: Mapping[str, object],
    *,
    applied_cursor: GraphReconcileCursor | None = None,
) -> tuple[list[GraphDelta], GraphReconcileCursor | None]:
    target_cursor = _cursor_from_upsert_result(upsert_result)
    if target_cursor is None:
        return [], None

    source_file = target_cursor.source_file
    run_rows = _fetch_domain_rows(
        conn,
        table="run_registry",
        source_file=source_file,
        applied_cursor=applied_cursor,
        target_cursor=target_cursor,
    )
    discovery_rows = _fetch_domain_rows(
        conn,
        table="discovery_records",
        source_file=source_file,
        applied_cursor=applied_cursor,
        target_cursor=target_cursor,
    )
    coverage_rows = _fetch_domain_rows(
        conn,
        table="coverage_state",
        source_file=source_file,
        applied_cursor=applied_cursor,
        target_cursor=target_cursor,
    )

    deltas: list[GraphDelta] = []
    for run_row in run_rows:
        identity = _require_identity(run_row, CANONICAL_IDENTITY_FIELDS, entity_type="run")
        deltas.append(GraphDelta(entity_type="run", identity=identity, record={**run_row, "entity_type": "run"}))

    for discovery_row in discovery_rows:
        identity = _require_identity(discovery_row, CANONICAL_IDENTITY_FIELDS, entity_type="discovery")
        deltas.append(
            GraphDelta(entity_type="discovery", identity=identity, record={**discovery_row, "entity_type": "discovery"})
        )
        deltas.append(
            GraphDelta(
                entity_type="discovery_provenance",
                identity=identity,
                record={
                    **discovery_row,
                    "edge_type": "discovery_provenance",
                    "source_entity_type": "discovery",
                    "target_entity_type": "run",
                },
            )
        )

    for coverage_row in coverage_rows:
        identity = _require_identity(coverage_row, CANONICAL_IDENTITY_FIELDS, entity_type="coverage_state")
        deltas.append(
            GraphDelta(entity_type="coverage_state", identity=identity, record={**coverage_row, "entity_type": "coverage_state"})
        )
        deltas.append(
            GraphDelta(
                entity_type="coverage_state_link",
                identity=identity,
                record={
                    **coverage_row,
                    "edge_type": "coverage_state_link",
                    "source_entity_type": "coverage_state",
                    "target_entity_type": "run",
                },
            )
        )

    deltas.extend(_run_lineage_deltas_from_rows(conn, run_rows))
    return deltas, target_cursor


def _normalize_node_delta(delta: GraphDelta) -> tuple[tuple[str, ...], dict[str, object]]:
    record = dict(_require_mapping(delta.record, f"{delta.entity_type} delta record"))
    identity = _require_identity(record, CANONICAL_IDENTITY_FIELDS, entity_type=delta.entity_type)
    if delta.identity != identity:
        raise GraphContractError(
            f"{delta.entity_type} delta identity does not match canonical record identity"
        )
    _ = record.setdefault("entity_type", delta.entity_type)
    return (delta.entity_type, *identity), record


def _normalize_attempt_edge_delta(delta: GraphDelta) -> tuple[tuple[str, ...], dict[str, object]]:
    record = dict(_require_mapping(delta.record, f"{delta.entity_type} delta record"))
    identity = _require_identity(record, CANONICAL_IDENTITY_FIELDS, entity_type=delta.entity_type)
    if delta.identity != identity:
        raise GraphContractError(
            f"{delta.entity_type} delta identity does not match canonical record identity"
        )
    _ = record.setdefault("edge_type", delta.entity_type)
    source_entity = str(record.get("source_entity_type", delta.entity_type.removesuffix("_provenance").removesuffix("_link")))
    target_entity = str(record.get("target_entity_type", "run"))
    record["source_entity_type"] = source_entity
    record["target_entity_type"] = target_entity
    return (delta.entity_type, *identity), record


def _normalize_run_lineage_delta(delta: GraphDelta) -> tuple[tuple[str, ...], dict[str, object]]:
    record = dict(_require_mapping(delta.record, "run_lineage delta record"))
    identity = _require_identity(
        record,
        ("campaign_id", "arm_id", "source_attempt_id", "target_attempt_id"),
        entity_type=delta.entity_type,
    )
    if delta.identity != identity:
        raise GraphContractError("run_lineage delta identity does not match lineage record identity")
    edge_type = str(record.get("edge_type", ""))
    if edge_type != "run_lineage":
        raise GraphContractError(
            f"Unsupported graph edge type: {edge_type or '<missing>'}. Supported edge types: {', '.join(SUPPORTED_EDGE_ENTITY_TYPES)}"
        )
    return (delta.entity_type, *identity), record


def apply_graph_delta(graph: MvpGraphState, delta: GraphDelta) -> GraphApplyResult:
    if delta.entity_type not in SUPPORTED_GRAPH_ENTITY_TYPES:
        raise GraphContractError(
            f"Unsupported graph entity type: {delta.entity_type}. Supported graph entity types: {_supported_types_message()}"
        )

    if delta.entity_type in SUPPORTED_NODE_ENTITY_TYPES:
        key, record = _normalize_node_delta(delta)
        bucket = graph.nodes
    elif delta.entity_type == "run_lineage":
        key, record = _normalize_run_lineage_delta(delta)
        bucket = graph.edges
    else:
        key, record = _normalize_attempt_edge_delta(delta)
        bucket = graph.edges

    existing = bucket.get(key)
    if existing is None:
        bucket[key] = record
        return GraphApplyResult(entity_type=delta.entity_type, identity=delta.identity, operation="insert")
    if existing == record:
        return GraphApplyResult(entity_type=delta.entity_type, identity=delta.identity, operation="noop")

    bucket[key] = record
    return GraphApplyResult(entity_type=delta.entity_type, identity=delta.identity, operation="update")


def apply_graph_deltas(graph: MvpGraphState, deltas: Sequence[GraphDelta]) -> GraphApplySummary:
    inserted = 0
    updated = 0
    noops = 0
    for delta in deltas:
        result = apply_graph_delta(graph, delta)
        if result.operation == "insert":
            inserted += 1
        elif result.operation == "update":
            updated += 1
        else:
            noops += 1
    return GraphApplySummary(
        inserted=inserted,
        updated=updated,
        noops=noops,
        node_count=len(graph.nodes),
        edge_count=len(graph.edges),
    )


def graph_deltas_from_contract_artifact(artifact: Mapping[str, object]) -> list[GraphDelta]:
    validated = validate_mvp_contract_artifact(artifact)
    deltas: list[GraphDelta] = []

    for run_record in validated["run_registry"]:
        identity = _require_identity(run_record, CANONICAL_IDENTITY_FIELDS, entity_type="run")
        deltas.append(GraphDelta(entity_type="run", identity=identity, record={**run_record, "entity_type": "run"}))

    for discovery_record in validated["discovery_records"]:
        identity = _require_identity(discovery_record, CANONICAL_IDENTITY_FIELDS, entity_type="discovery")
        deltas.append(
            GraphDelta(
                entity_type="discovery",
                identity=identity,
                record={**discovery_record, "entity_type": "discovery"},
            )
        )
        deltas.append(
            GraphDelta(
                entity_type="discovery_provenance",
                identity=identity,
                record={
                    **{field: discovery_record[field] for field in CANONICAL_IDENTITY_FIELDS},
                    "edge_type": "discovery_provenance",
                    "source_entity_type": "discovery",
                    "target_entity_type": "run",
                },
            )
        )

    for coverage_record in validated["coverage_state"]:
        identity = _require_identity(coverage_record, CANONICAL_IDENTITY_FIELDS, entity_type="coverage_state")
        deltas.append(
            GraphDelta(
                entity_type="coverage_state",
                identity=identity,
                record={**coverage_record, "entity_type": "coverage_state"},
            )
        )
        deltas.append(
            GraphDelta(
                entity_type="coverage_state_link",
                identity=identity,
                record={
                    **{field: coverage_record[field] for field in CANONICAL_IDENTITY_FIELDS},
                    "edge_type": "coverage_state_link",
                    "source_entity_type": "coverage_state",
                    "target_entity_type": "run",
                },
            )
        )

    for lineage_record in validated["graph_edges"]:
        edge_type = str(lineage_record.get("edge_type", ""))
        if edge_type != "run_lineage":
            raise GraphContractError(
                f"Unsupported graph edge type: {edge_type or '<missing>'}. Supported edge types: {', '.join(SUPPORTED_EDGE_ENTITY_TYPES)}"
            )
        identity = _require_identity(
            lineage_record,
            ("campaign_id", "arm_id", "source_attempt_id", "target_attempt_id"),
            entity_type="run_lineage",
        )
        deltas.append(GraphDelta(entity_type="run_lineage", identity=identity, record=dict(lineage_record)))

    return deltas


def apply_reconciled_graph_artifact(
    graph: MvpGraphState,
    artifact: Mapping[str, object],
) -> GraphApplySummary:
    return apply_graph_deltas(graph, graph_deltas_from_contract_artifact(artifact))


def apply_reconciled_graph_delta(
    graph: MvpGraphState,
    conn: sqlite3.Connection,
    upsert_result: Mapping[str, object],
    *,
    applied_cursor: GraphReconcileCursor | None = None,
) -> ReconciledGraphDeltaResult:
    deltas, cursor = graph_deltas_from_reconciled_state(conn, upsert_result, applied_cursor=applied_cursor)
    summary = apply_graph_deltas(graph, deltas)
    return ReconciledGraphDeltaResult(deltas=tuple(deltas), summary=summary, cursor=cursor)
