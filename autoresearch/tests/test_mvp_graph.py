from __future__ import annotations
# pyright: reportMissingImports=false, reportUnknownVariableType=false, reportUnknownMemberType=false, reportAny=false

import copy
import json
from pathlib import Path

import pytest

from autoresearch.harness.mvp_graph import (
    GraphContractError,
    GraphDelta,
    GraphReconcileCursor,
    MvpGraphState,
    apply_graph_delta,
    apply_graph_deltas,
    apply_reconciled_graph_delta,
    apply_reconciled_graph_artifact,
    graph_deltas_from_reconciled_state,
    graph_deltas_from_contract_artifact,
)
from autoresearch.harness.mvp_upsert_pipeline import execute_registry_discovery_coverage_upsert
from autoresearch.harness.telemetry_store import init_db


def _artifact_fixture() -> dict[str, list[dict[str, object]]]:
    return {
        "run_registry": [
            {
                "campaign_id": "camp-001",
                "arm_id": "adaptive-v1",
                "run_dir_fingerprint": "runfp-001",
                "attempt_id": "attempt-001",
                "artifact_status": "complete",
            }
        ],
        "discovery_records": [
            {
                "campaign_id": "camp-001",
                "arm_id": "adaptive-v1",
                "run_dir_fingerprint": "runfp-001",
                "attempt_id": "attempt-001",
                "artifact_status": "complete",
                "discovery_status": "new",
            }
        ],
        "coverage_state": [
            {
                "campaign_id": "camp-001",
                "arm_id": "adaptive-v1",
                "run_dir_fingerprint": "runfp-001",
                "attempt_id": "attempt-001",
                "artifact_status": "complete",
                "coverage_fraction": 0.5,
            }
        ],
        "graph_edges": [
            {
                "campaign_id": "camp-001",
                "arm_id": "adaptive-v1",
                "source_attempt_id": "attempt-001",
                "target_attempt_id": "attempt-002",
                "edge_type": "run_lineage",
            }
        ],
    }


def test_graph_deltas_from_contract_artifact_stays_within_mvp_scope() -> None:
    deltas = graph_deltas_from_contract_artifact(_artifact_fixture())

    assert [delta.entity_type for delta in deltas] == [
        "run",
        "discovery",
        "discovery_provenance",
        "coverage_state",
        "coverage_state_link",
        "run_lineage",
    ]
    assert deltas[0].record["campaign_id"] == "camp-001"
    assert deltas[0].record["arm_id"] == "adaptive-v1"
    assert deltas[0].record["run_dir_fingerprint"] == "runfp-001"
    assert deltas[0].record["attempt_id"] == "attempt-001"


def test_apply_graph_delta_is_idempotent_for_duplicate_reapply() -> None:
    graph = MvpGraphState()
    delta = graph_deltas_from_contract_artifact(_artifact_fixture())[0]

    first = apply_graph_delta(graph, delta)
    second = apply_graph_delta(graph, delta)

    assert first.operation == "insert"
    assert second.operation == "noop"
    assert len(graph.nodes) == 1
    assert len(graph.edges) == 0


def test_apply_graph_delta_updates_existing_record_when_payload_changes() -> None:
    graph = MvpGraphState()
    delta = graph_deltas_from_contract_artifact(_artifact_fixture())[0]
    updated_delta = GraphDelta(
        entity_type=delta.entity_type,
        identity=delta.identity,
        record={**delta.record, "artifact_status": "partial"},
    )

    _ = apply_graph_delta(graph, delta)
    result = apply_graph_delta(graph, updated_delta)

    assert result.operation == "update"
    assert graph.nodes[("run", *delta.identity)]["artifact_status"] == "partial"


def test_apply_reconciled_graph_artifact_is_idempotent_on_reapply() -> None:
    graph = MvpGraphState()
    artifact = _artifact_fixture()

    first = apply_reconciled_graph_artifact(graph, artifact)
    second = apply_reconciled_graph_artifact(graph, artifact)

    assert first.inserted == 6
    assert first.updated == 0
    assert first.noops == 0
    assert first.node_count == 3
    assert first.edge_count == 3
    assert second.inserted == 0
    assert second.updated == 0
    assert second.noops == 6
    assert second.node_count == 3
    assert second.edge_count == 3


def test_apply_graph_delta_rejects_unsupported_entity_type() -> None:
    graph = MvpGraphState()
    delta = GraphDelta(
        entity_type="campaign",
        identity=("camp-001",),
        record={"campaign_id": "camp-001"},
    )

    with pytest.raises(GraphContractError, match="Unsupported graph entity type: campaign"):
        _ = apply_graph_delta(graph, delta)


def test_graph_deltas_from_contract_artifact_rejects_unsupported_edge_scope() -> None:
    artifact = _artifact_fixture()
    artifact["graph_edges"] = [
        {
            "campaign_id": "camp-001",
            "arm_id": "adaptive-v1",
            "source_attempt_id": "attempt-001",
            "target_attempt_id": "attempt-002",
            "edge_type": "follows",
        }
    ]

    with pytest.raises(GraphContractError, match="Unsupported graph edge type: follows"):
        _ = graph_deltas_from_contract_artifact(artifact)


def test_apply_graph_delta_rejects_identity_mismatch() -> None:
    delta = graph_deltas_from_contract_artifact(_artifact_fixture())[0]
    mismatched = GraphDelta(
        entity_type=delta.entity_type,
        identity=("camp-001", "adaptive-v1", "runfp-other", "attempt-001"),
        record=copy.deepcopy(delta.record),
    )

    with pytest.raises(GraphContractError, match="identity does not match canonical record identity"):
        _ = apply_graph_delta(MvpGraphState(), mismatched)


@pytest.fixture
def basic_config() -> dict[str, object]:
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


def _graph_fingerprint(graph: MvpGraphState) -> str:
    payload = {
        "nodes": {"|".join(key): graph.nodes[key] for key in sorted(graph.nodes)},
        "edges": {"|".join(key): graph.edges[key] for key in sorted(graph.edges)},
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def test_task7_graph_delta_executor_adds_only_new_graph_records(
    tmp_path: Path,
    basic_config: dict[str, object],
) -> None:
    conn = init_db(tmp_path / "telemetry.db")
    source_file = tmp_path / "scanner_rows.jsonl"
    graph = MvpGraphState()

    baseline_upsert = execute_registry_discovery_coverage_upsert(
        conn,
        source_file=source_file,
        scanner_rows=_scanner_rows(),
        config=basic_config,
    )
    baseline_apply = apply_reconciled_graph_delta(graph, conn, baseline_upsert)

    assert baseline_apply.summary.inserted == 11
    assert baseline_apply.summary.updated == 0
    assert baseline_apply.summary.noops == 0
    assert baseline_apply.summary.node_count == 6
    assert baseline_apply.summary.edge_count == 5
    assert baseline_apply.cursor == GraphReconcileCursor(
        source_file=str(source_file),
        last_line_offset=2,
        last_checksum=baseline_upsert["watermark"]["last_checksum"],
    )

    next_rows = [
        *_scanner_rows(),
        {
            "campaign_id": "camp-001",
            "arm_id": "adaptive-v1",
            "run_dir_fingerprint": "run-003",
            "attempt_id": "attempt-003",
            "artifact_status": "complete",
            "discovery_status": "new",
            "axes_binned": {"tb": 10000, "lam1_bin": 2},
        },
    ]
    incremental_upsert = execute_registry_discovery_coverage_upsert(
        conn,
        source_file=source_file,
        scanner_rows=next_rows,
        config=basic_config,
    )
    incremental_apply = apply_reconciled_graph_delta(
        graph,
        conn,
        incremental_upsert,
        applied_cursor=baseline_apply.cursor,
    )

    assert incremental_apply.summary.inserted == 6
    assert incremental_apply.summary.updated == 0
    assert incremental_apply.summary.noops == 0
    assert incremental_apply.summary.node_count == 9
    assert incremental_apply.summary.edge_count == 8
    assert [delta.entity_type for delta in incremental_apply.deltas] == [
        "run",
        "discovery",
        "discovery_provenance",
        "coverage_state",
        "coverage_state_link",
        "run_lineage",
    ]

    discovery_edge = graph.edges[("discovery_provenance", "camp-001", "adaptive-v1", "run-003", "attempt-003")]
    assert discovery_edge["source_file"] == str(source_file)
    assert discovery_edge["line_offset"] == 3
    assert discovery_edge["source_checksum"] == incremental_upsert["watermark"]["last_checksum"]
    lineage_edge = graph.edges[("run_lineage", "camp-001", "adaptive-v1", "attempt-002", "attempt-003")]
    assert lineage_edge["source_run_dir_fingerprint"] == "run-002"
    assert lineage_edge["target_run_dir_fingerprint"] == "run-003"
    assert lineage_edge["source_file"] == str(source_file)


def test_task7_duplicate_delta_replay_is_graph_noop(tmp_path: Path, basic_config: dict[str, object]) -> None:
    conn = init_db(tmp_path / "telemetry.db")
    source_file = tmp_path / "scanner_rows.jsonl"
    graph = MvpGraphState()

    upsert_result = execute_registry_discovery_coverage_upsert(
        conn,
        source_file=source_file,
        scanner_rows=_scanner_rows(),
        config=basic_config,
    )
    deltas, cursor = graph_deltas_from_reconciled_state(conn, upsert_result)
    first = apply_graph_deltas(graph, deltas)
    before_hash = _graph_fingerprint(graph)
    second = apply_graph_deltas(graph, deltas)
    after_hash = _graph_fingerprint(graph)

    assert cursor == GraphReconcileCursor(
        source_file=str(source_file),
        last_line_offset=2,
        last_checksum=upsert_result["watermark"]["last_checksum"],
    )
    assert first.inserted == 11
    assert second.inserted == 0
    assert second.updated == 0
    assert second.noops == len(deltas)
    assert first.node_count == second.node_count == 6
    assert first.edge_count == second.edge_count == 5
    assert before_hash == after_hash
