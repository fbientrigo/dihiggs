from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from autoresearch.harness.mvp_upsert_pipeline import execute_registry_discovery_coverage_upsert
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
        },
        {
            "campaign_id": "camp-001",
            "arm_id": "adaptive-v1",
            "run_dir_fingerprint": "run-002",
            "attempt_id": "attempt-002",
            "artifact_status": "complete",
            "discovery_status": "new",
            "axes_binned": {"tb": 5000, "lam1_bin": 1},
        },
    ]


def test_task6_upsert_first_insert_success_path(tmp_path: Path, basic_config: dict[str, Any]) -> None:
    conn = init_db(tmp_path / "telemetry.db")

    result = execute_registry_discovery_coverage_upsert(
        conn,
        source_file=tmp_path / "scanner_rows.jsonl",
        scanner_rows=_scanner_rows(),
        config=basic_config,
    )

    counts = db_counts(conn)
    assert result["status"] == "success"
    assert result["rows_valid"] == 2
    assert result["rows_quarantined"] == 0
    assert result["run_registry"]["insert"] == 2
    assert result["discovery_records"]["insert"] == 2
    assert result["coverage_state"]["insert"] == 2
    assert counts["run_registry"] == 2
    assert counts["discovery_records"] == 2
    assert counts["coverage_state"] == 2
    assert result["watermark"]["last_line_offset"] == 2

    cursor = conn.cursor()
    cursor.execute("SELECT coverage_fraction FROM coverage_state ORDER BY line_offset")
    fractions = [row[0] for row in cursor.fetchall()]
    assert fractions[0] > 0.0
    assert fractions[-1] >= fractions[0]


def test_task6_upsert_rerun_duplicate_is_noop(tmp_path: Path, basic_config: dict[str, Any]) -> None:
    conn = init_db(tmp_path / "telemetry.db")
    source_file = tmp_path / "scanner_rows.jsonl"

    first = execute_registry_discovery_coverage_upsert(
        conn,
        source_file=source_file,
        scanner_rows=_scanner_rows(),
        config=basic_config,
    )
    second = execute_registry_discovery_coverage_upsert(
        conn,
        source_file=source_file,
        scanner_rows=_scanner_rows(),
        config=basic_config,
    )

    counts = db_counts(conn)
    assert first["run_registry"]["insert"] == 2
    assert second["rows_valid"] == 0
    assert second["rows_skipped_by_watermark"] == 2
    assert second["run_registry"] == {"insert": 0, "update": 0, "noop": 0}
    assert second["discovery_records"] == {"insert": 0, "update": 0, "noop": 0}
    assert second["coverage_state"] == {"insert": 0, "update": 0, "noop": 0}
    assert counts["run_registry"] == 2
    assert counts["discovery_records"] == 2
    assert counts["coverage_state"] == 2


def test_task6_upsert_quarantines_malformed_rows_while_valid_rows_persist(
    tmp_path: Path,
    basic_config: dict[str, Any],
) -> None:
    conn = init_db(tmp_path / "telemetry.db")
    rows: list[object] = [
        _scanner_rows()[0],
        {
            "campaign_id": "camp-001",
            "arm_id": "adaptive-v1",
            "run_dir_fingerprint": "run-bad",
            "artifact_status": "complete",
            "discovery_status": "new",
        },
        _scanner_rows()[1],
    ]

    result = execute_registry_discovery_coverage_upsert(
        conn,
        source_file=tmp_path / "scanner_rows.jsonl",
        scanner_rows=rows,
        config=basic_config,
    )

    counts = db_counts(conn)
    assert result["rows_valid"] == 2
    assert result["rows_quarantined"] == 1
    assert counts["run_registry"] == 2
    assert counts["discovery_records"] == 2
    assert counts["coverage_state"] == 2
    assert counts["upsert_quarantine"] == 1
    assert result["watermark"]["last_line_offset"] == 3

    cursor = conn.cursor()
    cursor.execute("SELECT error_type, error_message FROM upsert_quarantine")
    error_type, error_message = cursor.fetchone()
    assert error_type == "schema_validation_error"
    assert "attempt_id" in error_message
