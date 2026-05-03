from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from autoresearch.harness.mvp_graph import MvpGraphState
from autoresearch.harness.mvp_upsert_pipeline import trigger_incremental_update
from autoresearch.harness.telemetry_store import db_counts, init_db


@pytest.fixture
def basic_config() -> dict[str, Any]:
    return {
        "search": {
            "tb_values": [1000, 5000, 10000, 30000],
            "lam1": {"min": -20.0, "max": 20.0, "n_bins": 40},
            "mphi": {"min": 0.0, "max": 1000.0, "n_bins": 50},
        },
        "metrics": {
            "weights": {"yield": 0.5, "coverage": 0.3, "diversity": 0.2},
            "multi_axis": {
                "collapse_axes": ["tb", "lam1_bin"],
                "coverage_axes": [
                    {"name": "tb", "kind": "categorical", "domain": [1000, 5000, 10000, 30000], "weight": 0.5},
                    {"name": "lam1_bin", "kind": "discrete", "domain_size": 40, "weight": 0.5},
                ],
                "diversity_axes": [
                    {"name": "tb", "weight": 0.5},
                    {"name": "lam1_bin", "weight": 0.5},
                ],
                "diversity_pairs": [
                    {"axes": ["tb", "lam1_bin"], "weight": 1.0},
                ],
            },
        },
    }


def _scanner_rows() -> list[dict[str, object]]:
    return [
        {
            "campaign_id": "camp-001",
            "arm_id": "adaptive-v1",
            "run_dir_fingerprint": "run-001",
            "attempt_id": "attempt-001",
            "artifact_status": "complete",
            "discovery_status": "new",
            "axes_binned": {"tb": 1000, "lam1_bin": 0},
            "coverage_fraction": 0.1,
        },
        {
            "campaign_id": "camp-001",
            "arm_id": "adaptive-v1",
            "run_dir_fingerprint": "run-002",
            "attempt_id": "attempt-002",
            "artifact_status": "complete",
            "discovery_status": "new",
            "axes_binned": {"tb": 5000, "lam1_bin": 1},
            "coverage_fraction": 0.2,
        },
    ]


def test_task11_trigger_delta_path(tmp_path: Path, basic_config: dict[str, Any]) -> None:
    conn = init_db(tmp_path / "telemetry.db")
    graph = MvpGraphState()
    source_file = tmp_path / "scanner_rows.jsonl"

    # First trigger - should have deltas
    result = trigger_incremental_update(
        conn,
        graph,
        source_file=source_file,
        scanner_rows=_scanner_rows(),
        config=basic_config,
    )

    assert result["status"] == "success"
    assert result["is_noop"] is False
    assert result["upsert"]["rows_valid"] == 2
    assert result["graph"]["summary"]["inserted"] > 0
    
    counts = db_counts(conn)
    assert counts["run_registry"] == 2
    assert len(graph.nodes) > 0


def test_task11_trigger_noop_path(tmp_path: Path, basic_config: dict[str, Any]) -> None:
    conn = init_db(tmp_path / "telemetry.db")
    graph = MvpGraphState()
    source_file = tmp_path / "scanner_rows.jsonl"

    # First trigger to establish state
    trigger_incremental_update(
        conn,
        graph,
        source_file=source_file,
        scanner_rows=_scanner_rows(),
        config=basic_config,
    )
    
    # Capture cursor from first result for the second call
    applied_cursor = trigger_incremental_update(
        conn,
        graph,
        source_file=source_file,
        scanner_rows=_scanner_rows(),
        config=basic_config,
    )["graph"]["cursor"]

    # Second trigger with same data - should be noop
    result = trigger_incremental_update(
        conn,
        graph,
        source_file=source_file,
        scanner_rows=_scanner_rows(),
        config=basic_config,
        applied_cursor=applied_cursor,
    )

    assert result["status"] == "success"
    assert result["is_noop"] is True
    assert result["upsert"]["rows_valid"] == 0
    assert result["graph"]["summary"]["inserted"] == 0
    assert result["graph"]["summary"]["updated"] == 0
