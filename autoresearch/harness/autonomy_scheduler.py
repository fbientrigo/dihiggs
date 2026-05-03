from __future__ import annotations

import json
import sqlite3
import threading
import time
from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import cast

from autoresearch.harness.mvp_goal_proposer import propose_parameter_space_goals
from autoresearch.harness.mvp_graph import GraphReconcileCursor, MvpGraphState, apply_reconciled_graph_delta
from autoresearch.harness.mvp_upsert_pipeline import execute_registry_discovery_coverage_upsert
from autoresearch.harness.orchestrator_adapter import DryRunOrchestratorAdapter, OrchestratorAdapter


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _as_dict(value: object) -> dict[str, object]:
    return dict(cast(Mapping[str, object], value)) if isinstance(value, Mapping) else {}


def _as_int(value: object, default: int) -> int:
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


def _as_float(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return None
    return None


def _stable_json(value: object) -> str:
    return json.dumps(value, indent=2, sort_keys=True) + "\n"


def _scheduler_config(config: Mapping[str, object]) -> dict[str, object]:
    scheduler_cfg = _as_dict(config.get("autonomy_scheduler"))
    if not scheduler_cfg:
        scheduler_cfg = _as_dict(config.get("scheduler"))
    if not scheduler_cfg:
        scheduler_cfg = _as_dict(config.get("autonomy"))
    return scheduler_cfg


def _graph_cursor_payload(cursor: GraphReconcileCursor | None) -> dict[str, object] | None:
    if cursor is None:
        return None
    return {
        "source_file": cursor.source_file,
        "last_line_offset": cursor.last_line_offset,
        "last_checksum": cursor.last_checksum,
    }


def _graph_counts(graph: MvpGraphState) -> dict[str, int]:
    return {"node_count": len(graph.nodes), "edge_count": len(graph.edges)}


def _fetch_payload_rows(conn: sqlite3.Connection, table: str) -> list[dict[str, object]]:
    cursor = conn.cursor()
    _ = cursor.execute(f"SELECT payload FROM {table} ORDER BY source_file, line_offset, id")
    rows = cast(list[tuple[object]], cursor.fetchall())
    payload_rows: list[dict[str, object]] = []
    for row in rows:
        payload_rows.append(cast(dict[str, object], json.loads(str(row[0]))))
    return payload_rows


def _budget_config(config: Mapping[str, object]) -> dict[str, object]:
    scheduler_cfg = _scheduler_config(config)

    interval_seconds = _as_float(scheduler_cfg.get("interval_seconds"))
    if interval_seconds is None:
        interval_seconds = 900.0

    max_proposals = _as_int(scheduler_cfg.get("max_proposals"), 3)
    max_dispatches = _as_int(scheduler_cfg.get("max_dispatches"), 1)
    time_cap_seconds = _as_float(scheduler_cfg.get("time_cap_seconds"))
    dispatch_mode = str(scheduler_cfg.get("dispatch_mode", "dry_run") or "dry_run")

    return {
        "interval_seconds": max(0.0, interval_seconds),
        "max_proposals": max(0, max_proposals),
        "max_dispatches": max(0, max_dispatches),
        "time_cap_seconds": time_cap_seconds,
        "dispatch_mode": dispatch_mode,
    }


def _dispatch_payload(proposer_result: Mapping[str, object], dispatch_limit: int) -> dict[str, object]:
    proposals = proposer_result.get("proposals")
    if not isinstance(proposals, Sequence) or isinstance(proposals, (str, bytes)):
        bounded_proposals: list[object] = []
    else:
        bounded_proposals = list(proposals[:dispatch_limit])
    return {
        "status": proposer_result.get("status", "saturated"),
        "proposals": bounded_proposals,
    }


def _stage_event(
    *,
    tick_index: int,
    stage: str,
    status: str,
    started_at: str,
    finished_at: str,
    details: Mapping[str, object],
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "tick_index": tick_index,
        "stage": stage,
        "status": status,
        "started_at": started_at,
        "finished_at": finished_at,
        "details": dict(details),
    }


def _quarantine_config(config: Mapping[str, object]) -> dict[str, int]:
    quarantine_cfg = _as_dict(_scheduler_config(config).get("quarantine"))
    corrupted_threshold = _as_int(quarantine_cfg.get("corrupted_threshold", 1), 1)
    incomplete_threshold = _as_int(quarantine_cfg.get("incomplete_threshold", 1), 1)
    return {
        "corrupted_threshold": max(1, corrupted_threshold),
        "incomplete_threshold": max(1, incomplete_threshold),
    }


def _artifact_quarantine_summary(
    *,
    upsert_result: Mapping[str, object],
    scanner_rows: Sequence[object],
    config: Mapping[str, object],
) -> dict[str, object]:
    corrupted_count = _as_int(upsert_result.get("rows_quarantined", 0), 0)
    incomplete_count = 0
    affected_examples: list[str] = []
    for index, raw_row in enumerate(scanner_rows, start=1):
        if not isinstance(raw_row, Mapping):
            continue
        artifact_status = str(raw_row.get("artifact_status", "") or "")
        parser_status = str(raw_row.get("parser_status", "") or "")
        discovery_status = str(raw_row.get("discovery_status", "") or "")
        is_corrupted = discovery_status == "quarantined"
        is_incomplete = not is_corrupted and (
            artifact_status == "partial"
            or discovery_status == "partial"
            or parser_status == "missing"
        )
        if is_corrupted:
            corrupted_count += 1
        elif is_incomplete:
            incomplete_count += 1
        if (is_corrupted or is_incomplete) and len(affected_examples) < 5:
            affected_examples.append(f"line:{index}")

    thresholds = _quarantine_config(config)
    status = "clear"
    decision = "proceed"
    reason = "no_quarantine"
    if corrupted_count >= thresholds["corrupted_threshold"]:
        status = "quarantine_halt"
        decision = "halt"
        reason = "corrupted_threshold_reached"
    elif incomplete_count >= thresholds["incomplete_threshold"]:
        status = "quarantine_skip_dispatch"
        decision = "skip_dispatch"
        reason = "incomplete_threshold_reached"

    return {
        "status": status,
        "decision": decision,
        "reason": reason,
        "corrupted_count": corrupted_count,
        "incomplete_count": incomplete_count,
        "thresholds": thresholds,
        "affected_examples": affected_examples,
    }


def _proposal_audit_records(
    proposals: Sequence[Mapping[str, object]],
    dispatch_records: Sequence[Mapping[str, object]],
    *,
    quarantine: Mapping[str, object],
    dispatch_status: str,
) -> list[dict[str, object]]:
    records_by_goal_id = {
        str(record.get("goal_id", "")): dict(record)
        for record in dispatch_records
        if str(record.get("goal_id", ""))
    }
    quarantine_decision = str(quarantine.get("decision", "proceed"))
    quarantine_reason = str(quarantine.get("reason", "no_quarantine"))
    audit_records: list[dict[str, object]] = []
    for proposal in proposals:
        goal_id = str(proposal.get("goal_id", ""))
        dispatch_record = records_by_goal_id.get(goal_id)
        if dispatch_record is not None:
            outcome = "dispatched"
            outcome_reason = str(dispatch_record.get("status", dispatch_status or "success"))
        elif quarantine_decision == "halt":
            outcome = "skipped"
            outcome_reason = "quarantine_halt"
        elif quarantine_decision == "skip_dispatch":
            outcome = "skipped"
            outcome_reason = "quarantine_skip_dispatch"
        else:
            outcome = "deferred"
            outcome_reason = "dispatch_budget_exhausted"
        audit_records.append(
            {
                "goal_id": goal_id,
                "cell_id": str(proposal.get("cell_id", "")),
                "axes_binned": proposal.get("axes_binned", {}),
                "rationale": proposal.get("rationale", {}),
                "decision": {
                    "outcome": outcome,
                    "reason": outcome_reason,
                    "dispatch_status": dispatch_status,
                    "quarantine_reason": quarantine_reason,
                },
                "dispatch_record": dispatch_record,
            }
        )
    return audit_records


class FixedIntervalAutonomyScheduler:
    def __init__(
        self,
        conn: sqlite3.Connection,
        *,
        config: Mapping[str, object],
        source_file: str | Path,
        state_path: str | Path,
        adapter: OrchestratorAdapter | None = None,
        graph: MvpGraphState | None = None,
        utc_now_func: Callable[[], str] = _utc_now_iso,
        monotonic_func: Callable[[], float] = time.monotonic,
        sleep_func: Callable[[float], None] = time.sleep,
    ) -> None:
        self.conn: sqlite3.Connection = conn
        self.config: dict[str, object] = dict(config)
        self.source_file: Path = Path(source_file)
        self.state_path: Path = Path(state_path)
        self.graph: MvpGraphState = graph or MvpGraphState()
        self.adapter: OrchestratorAdapter = adapter or DryRunOrchestratorAdapter()
        self.utc_now_func: Callable[[], str] = utc_now_func
        self.monotonic_func: Callable[[], float] = monotonic_func
        self.sleep_func: Callable[[float], None] = sleep_func
        self._lock: threading.Lock = threading.Lock()
        self._tick_index: int = 0
        self._graph_cursor: GraphReconcileCursor | None = None
        self._scheduled_goals: list[dict[str, object]] = []
        self._last_snapshot: dict[str, object] = {
            "schema_version": 1,
            "tick_index": 0,
            "last_status": "idle",
            "graph_cursor": None,
            "graph_counts": _graph_counts(self.graph),
            "scheduled_goals": [],
            "dispatch_records": [],
        }

    def run(
        self,
        *,
        scanner_rows_provider: Callable[[], Sequence[object]],
        max_ticks: int | None = None,
    ) -> list[dict[str, object]]:
        results: list[dict[str, object]] = []
        interval_seconds = cast(float, _budget_config(self.config)["interval_seconds"])
        completed_ticks = 0
        while max_ticks is None or completed_ticks < max_ticks:
            results.append(self.tick(scanner_rows=scanner_rows_provider()))
            completed_ticks += 1
            if max_ticks is not None and completed_ticks >= max_ticks:
                break
            if interval_seconds > 0:
                self.sleep_func(interval_seconds)
        return results

    def tick(self, *, scanner_rows: Sequence[object]) -> dict[str, object]:
        requested_at = self.utc_now_func()
        if not self._lock.acquire(blocking=False):
            return {
                "status": "lock_conflict",
                "requested_at": requested_at,
                "tick_index": self._tick_index,
                "mutated": False,
            }

        budgets = _budget_config(self.config)
        started_at = requested_at
        start_monotonic = self.monotonic_func()
        tick_index = self._tick_index + 1
        status = "success"
        mutated = False
        dispatch_result: Mapping[str, object] = {"status": "skipped", "dispatch_count": 0, "records": []}
        stage_telemetry: list[dict[str, object]] = []
        try:
            upsert_result = cast(
                dict[str, object],
                execute_registry_discovery_coverage_upsert(
                self.conn,
                source_file=self.source_file,
                scanner_rows=scanner_rows,
                config=self.config,
                utc_now_func=self.utc_now_func,
                ),
            )
            domain_result = cast(dict[str, dict[str, int]], upsert_result)
            rows_quarantined: object = upsert_result.get("rows_quarantined", 0)
            mutated = any(
                domain_result[key][name] > 0
                for key in ("run_registry", "discovery_records", "coverage_state")
                for name in ("insert", "update")
            ) or _as_int(rows_quarantined, 0) > 0

            graph_result = apply_reconciled_graph_delta(
                self.graph,
                self.conn,
                upsert_result,
                applied_cursor=self._graph_cursor,
            )
            if graph_result.cursor is not None:
                self._graph_cursor = graph_result.cursor
            quarantine = _artifact_quarantine_summary(
                upsert_result=upsert_result,
                scanner_rows=scanner_rows,
                config=self.config,
            )
            stage_telemetry.append(
                _stage_event(
                    tick_index=tick_index,
                    stage="reconcile",
                    status=str(quarantine["status"]),
                    started_at=started_at,
                    finished_at=self.utc_now_func(),
                    details={
                        "rows_seen": upsert_result.get("rows_seen", 0),
                        "rows_valid": upsert_result.get("rows_valid", 0),
                        "rows_quarantined": upsert_result.get("rows_quarantined", 0),
                        "rows_skipped_by_watermark": upsert_result.get("rows_skipped_by_watermark", 0),
                        "graph_summary": {
                            "inserted": graph_result.summary.inserted,
                            "updated": graph_result.summary.updated,
                            "noops": graph_result.summary.noops,
                            "node_count": graph_result.summary.node_count,
                            "edge_count": graph_result.summary.edge_count,
                        },
                        "quarantine": quarantine,
                    },
                )
            )

            coverage_state = _fetch_payload_rows(self.conn, "coverage_state")
            discovery_records = _fetch_payload_rows(self.conn, "discovery_records")

            if self._time_cap_exceeded(start_monotonic, budgets):
                status = "time_cap_exceeded"
                proposer_result: dict[str, object] = {
                    "status": "saturated",
                    "max_new_runs": _as_int(budgets.get("max_proposals"), 0),
                    "coverage_fraction": coverage_state[-1]["coverage_fraction"] if coverage_state else 0.0,
                    "candidate_cells_considered": 0,
                    "suppressed_duplicates": {
                        "scheduled_duplicates": 0,
                        "covered_regions": 0,
                        "scheduled_regions": 0,
                    },
                    "fallback_policy": "time_cap_exceeded",
                    "proposals": [],
                }
            else:
                proposer_result = propose_parameter_space_goals(
                    coverage_state=coverage_state,
                    discovery_records=discovery_records,
                    graph=self.graph,
                    config=self.config,
                    scheduled_goals=self._scheduled_goals,
                    max_new_runs=_as_int(budgets.get("max_proposals"), 0),
                )
            proposal_stage_status = "success"
            if str(quarantine["decision"]) == "halt":
                proposal_stage_status = "quarantine_halt"
            elif str(quarantine["decision"]) == "skip_dispatch":
                proposal_stage_status = "quarantine_skip_dispatch"
            elif status == "time_cap_exceeded":
                proposal_stage_status = "time_cap_exceeded"
            elif proposer_result.get("status") == "saturated":
                proposal_stage_status = "saturated"
            stage_telemetry.append(
                _stage_event(
                    tick_index=tick_index,
                    stage="proposal",
                    status=proposal_stage_status,
                    started_at=started_at,
                    finished_at=self.utc_now_func(),
                    details={
                        "proposal_status": proposer_result.get("status", "unknown"),
                        "proposal_count": len(cast(list[object], proposer_result.get("proposals", []))),
                        "candidate_cells_considered": proposer_result.get("candidate_cells_considered", 0),
                        "suppressed_duplicates": proposer_result.get("suppressed_duplicates", {}),
                        "fallback_policy": proposer_result.get("fallback_policy"),
                        "quarantine": quarantine,
                    },
                )
            )

            dispatch_limit = _as_int(budgets.get("max_dispatches"), 0)
            dispatch_payload = _dispatch_payload(proposer_result, dispatch_limit)
            bounded_proposals = cast(list[dict[str, object]], dispatch_payload["proposals"])
            if str(quarantine["decision"]) == "halt":
                status = "quarantine_halt"
                dispatch_result = {
                    "status": "skipped",
                    "dispatch_count": 0,
                    "records": [],
                    "reason": "quarantine_halt",
                    "quarantine": quarantine,
                }
            elif str(quarantine["decision"]) == "skip_dispatch":
                status = "quarantine_skip_dispatch"
                dispatch_result = {
                    "status": "skipped",
                    "dispatch_count": 0,
                    "records": [],
                    "reason": "quarantine_skip_dispatch",
                    "quarantine": quarantine,
                }
            elif status != "time_cap_exceeded" and bounded_proposals:
                if str(budgets["dispatch_mode"]) == "dry_run":
                    dispatch_result = self.adapter.dispatch(dispatch_payload)
                else:
                    dispatch_result = self.adapter.dispatch(dispatch_payload)
                if dispatch_result.get("status") == "error":
                    status = "dispatch_error"
                else:
                    self._scheduled_goals.extend(bounded_proposals)
            elif proposer_result.get("status") == "saturated":
                status = "saturated"
            stage_telemetry.append(
                _stage_event(
                    tick_index=tick_index,
                    stage="dispatch",
                    status=str(dispatch_result.get("status", status)),
                    started_at=started_at,
                    finished_at=self.utc_now_func(),
                    details={
                        "dispatch_count": _as_int(dispatch_result.get("dispatch_count", 0), 0),
                        "dispatch_mode": budgets.get("dispatch_mode"),
                        "reason": dispatch_result.get("reason"),
                        "quarantine": quarantine,
                    },
                )
            )

            finished_at = self.utc_now_func()
            elapsed_seconds = self.monotonic_func() - start_monotonic
            dispatch_records = dispatch_result.get("records", [])
            proposal_count = len(cast(list[object], proposer_result.get("proposals", [])))
            proposal_audit = _proposal_audit_records(
                cast(list[Mapping[str, object]], proposer_result.get("proposals", [])),
                cast(list[Mapping[str, object]], dispatch_records),
                quarantine=quarantine,
                dispatch_status=str(dispatch_result.get("status", "skipped")),
            )
            stage_telemetry.append(
                _stage_event(
                    tick_index=tick_index,
                    stage="scheduler",
                    status=status,
                    started_at=started_at,
                    finished_at=finished_at,
                    details={
                        "mutated": mutated,
                        "elapsed_seconds": elapsed_seconds,
                        "proposal_count": proposal_count,
                        "dispatch_count": _as_int(dispatch_result.get("dispatch_count", 0), 0),
                        "quarantine": quarantine,
                    },
                )
            )
            snapshot = {
                "schema_version": 1,
                "tick_index": tick_index,
                "started_at": started_at,
                "finished_at": finished_at,
                "last_status": status,
                "source_file": str(self.source_file),
                "budgets": budgets,
                "elapsed_seconds": elapsed_seconds,
                "mutated": mutated,
                "upsert_result": upsert_result,
                "graph_summary": {
                    "inserted": graph_result.summary.inserted,
                    "updated": graph_result.summary.updated,
                    "noops": graph_result.summary.noops,
                    "node_count": graph_result.summary.node_count,
                    "edge_count": graph_result.summary.edge_count,
                },
                "graph_cursor": _graph_cursor_payload(self._graph_cursor),
                "graph_counts": _graph_counts(self.graph),
                "proposal_summary": {
                    "status": proposer_result.get("status"),
                    "proposal_count": proposal_count,
                    "candidate_cells_considered": proposer_result.get("candidate_cells_considered", 0),
                    "suppressed_duplicates": proposer_result.get("suppressed_duplicates", {}),
                },
                "quarantine": quarantine,
                "stage_telemetry": stage_telemetry,
                "proposals": proposer_result.get("proposals", []),
                "proposal_audit": proposal_audit,
                "dispatch_result": dispatch_result,
                "scheduled_goals": self._scheduled_goals,
                "dispatch_records": dispatch_records,
            }
            self._persist_snapshot(snapshot)
            self._last_snapshot = snapshot
            self._tick_index = tick_index
            return {
                "status": status,
                "tick_index": tick_index,
                "mutated": mutated,
                "upsert_result": upsert_result,
                "graph_summary": snapshot["graph_summary"],
                "proposal_count": proposal_count,
                "dispatch_count": _as_int(dispatch_result.get("dispatch_count", 0), 0),
                "proposals": snapshot["proposals"],
                "proposal_audit": proposal_audit,
                "dispatch_records": dispatch_records,
                "stage_telemetry": stage_telemetry,
                "quarantine": quarantine,
                "state_path": str(self.state_path),
                "started_at": started_at,
                "finished_at": finished_at,
            }
        finally:
            self._lock.release()

    def _persist_snapshot(self, snapshot: Mapping[str, object]) -> None:
        self.state_path.parent.mkdir(parents=True, exist_ok=True)
        _ = self.state_path.write_text(_stable_json(snapshot), encoding="utf-8")

    def _time_cap_exceeded(self, started_monotonic: float, budgets: Mapping[str, object]) -> bool:
        time_cap_seconds = _as_float(budgets.get("time_cap_seconds"))
        if time_cap_seconds is None:
            return False
        return (self.monotonic_func() - started_monotonic) >= time_cap_seconds
