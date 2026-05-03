from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from collections.abc import Mapping, Sequence
from typing import cast

from autoresearch.harness.mvp_goal_proposer import propose_parameter_space_goals
from autoresearch.harness.mvp_graph import MvpGraphState, apply_reconciled_graph_delta
from autoresearch.harness.mvp_upsert_pipeline import execute_registry_discovery_coverage_upsert
from autoresearch.harness.telemetry_store import init_db


def _basic_config() -> dict[str, object]:
    return {
        "search": {
            "tb_values": [1000, 5000],
            "lam1": {"min": -20.0, "max": 20.0, "n_bins": 4},
            "mphi": {"min": 0.0, "max": 1000.0, "n_bins": 1},
        },
        "limits": {
            "max_new_run_dirs_per_round": 3,
        },
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
            },
        },
    }


def _fetch_payload_rows(conn: sqlite3.Connection, table: str) -> list[dict[str, object]]:
    cursor = conn.cursor()
    _ = cursor.execute(f"SELECT payload FROM {table} ORDER BY line_offset, id")
    fetched_rows = cast(list[tuple[object]], cursor.fetchall())
    return [cast(dict[str, object], json.loads(str(row[0]))) for row in fetched_rows]


def _reconcile_state(
    tmp_path: Path,
    scanner_rows: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], MvpGraphState]:
    conn = init_db(tmp_path / "telemetry.db")
    source_file = tmp_path / "scanner_rows.jsonl"
    graph = MvpGraphState()
    upsert_result = execute_registry_discovery_coverage_upsert(
        conn,
        source_file=source_file,
        scanner_rows=list(scanner_rows),
        config=_basic_config(),
    )
    _ = apply_reconciled_graph_delta(graph, conn, upsert_result)
    coverage_state = _fetch_payload_rows(conn, "coverage_state")
    discovery_records = _fetch_payload_rows(conn, "discovery_records")
    return coverage_state, discovery_records, graph


def test_task8_goal_proposer_generates_bounded_gap_driven_goals(tmp_path: Path) -> None:
    scanner_rows = [
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
        {
            "campaign_id": "camp-001",
            "arm_id": "adaptive-v1",
            "run_dir_fingerprint": "run-003",
            "attempt_id": "attempt-003",
            "artifact_status": "complete",
            "discovery_status": "new",
            "axes_binned": {"tb": 5000, "lam1_bin": 0},
        },
    ]
    coverage_state, discovery_records, graph = _reconcile_state(tmp_path, scanner_rows)

    result = propose_parameter_space_goals(
        coverage_state=coverage_state,
        discovery_records=discovery_records,
        graph=graph,
        config=_basic_config(),
        scheduled_goals=[
            {"axes_binned": {"tb": 5000, "lam1_bin": 1}},
            {"axes_binned": {"tb": 5000, "lam1_bin": 1}},
        ],
        max_new_runs=2,
    )

    assert result["status"] == "success"
    assert result["max_new_runs"] == 2
    proposals = cast(list[dict[str, object]], result["proposals"])
    suppressed = cast(dict[str, int], result["suppressed_duplicates"])
    assert len(proposals) == 2
    assert suppressed["scheduled_duplicates"] == 1
    assert suppressed["scheduled_regions"] >= 1

    first = proposals[0]
    second = proposals[1]
    assert first["axes_binned"] == {"tb": 1000, "lam1_bin": 2}
    assert second["axes_binned"] not in (
        {"tb": 1000, "lam1_bin": 0},
        {"tb": 1000, "lam1_bin": 1},
        {"tb": 5000, "lam1_bin": 0},
        {"tb": 5000, "lam1_bin": 1},
    )

    for index, proposal in enumerate(proposals, start=1):
        rationale = cast(dict[str, object], proposal["rationale"])
        assert rationale["signals"]
        coverage_gap = cast(dict[str, object], rationale["coverage_gap"])
        budget = cast(dict[str, object], rationale["budget"])
        assert cast(float, coverage_gap["coverage_gain"]) > 0.0
        assert budget["max_new_runs"] == 2
        assert budget["selected_rank"] == index
        assert rationale["source_refs"]

    first_rationale = cast(dict[str, object], first["rationale"])
    first_frontier = cast(dict[str, object], first_rationale["discovery_frontier"])
    assert cast(int, first_frontier["frontier_neighbor_count"]) >= 1
    assert "tb=1000|lam1_bin=1" in cast(list[str], first_frontier["neighbor_cell_ids"])


def test_task8_goal_proposer_returns_zero_when_slice_is_saturated(tmp_path: Path) -> None:
    scanner_rows: list[dict[str, object]] = []
    attempt_index = 0
    for tb in (1000, 5000):
        for lam1_bin in range(4):
            attempt_index += 1
            scanner_rows.append(
                {
                    "campaign_id": "camp-001",
                    "arm_id": "adaptive-v1",
                    "run_dir_fingerprint": f"run-{attempt_index:03d}",
                    "attempt_id": f"attempt-{attempt_index:03d}",
                    "artifact_status": "complete",
                    "discovery_status": "new",
                    "axes_binned": {"tb": tb, "lam1_bin": lam1_bin},
                }
            )

    coverage_state, discovery_records, graph = _reconcile_state(tmp_path, scanner_rows)
    result = propose_parameter_space_goals(
        coverage_state=coverage_state,
        discovery_records=discovery_records,
        graph=graph,
        config=_basic_config(),
        max_new_runs=5,
    )

    assert result == {
        "status": "saturated",
        "max_new_runs": 5,
        "coverage_fraction": 1.0,
        "candidate_cells_considered": 0,
        "suppressed_duplicates": {
            "scheduled_duplicates": 0,
            "covered_regions": 8,
            "scheduled_regions": 0,
        },
        "fallback_policy": "coverage_saturated",
        "proposals": [],
    }
