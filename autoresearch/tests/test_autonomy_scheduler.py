from __future__ import annotations

import json
import sqlite3
import threading
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import cast

from autoresearch.harness.autonomy_scheduler import FixedIntervalAutonomyScheduler
from autoresearch.harness.telemetry_store import init_db


class BlockingAdapter:
    def __init__(self) -> None:
        self.entered: threading.Event = threading.Event()
        self.release: threading.Event = threading.Event()

    def dispatch(self, payload: Mapping[str, object]) -> Mapping[str, object]:
        proposals = cast(list[dict[str, object]], payload.get("proposals", []))
        self.entered.set()
        _ = self.release.wait(timeout=5)
        return {
            "status": "success",
            "dispatch_count": len(proposals),
            "records": [
                {
                    "goal_id": str(proposal.get("goal_id")),
                    "status": "dry_run",
                    "dispatched_at": "2026-04-14T00:00:04+00:00",
                    "metadata": {"axes_binned": proposal.get("axes_binned")},
                }
                for proposal in proposals
            ],
        }


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
                "diversity_pairs": [
                    {"axes": ["tb", "lam1_bin"], "weight": 1.0},
                ],
            },
        },
        "autonomy_scheduler": {
            "interval_seconds": 900,
            "max_proposals": 2,
            "max_dispatches": 1,
            "dispatch_mode": "dry_run",
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
            "axes_binned": {"tb": 1000, "lam1_bin": 1},
        },
    ]


def _fixed_clock(value: str) -> Callable[[], str]:
    return lambda: value


def test_task10_scheduler_tick_happy_path_persists_bounded_snapshot(tmp_path: Path) -> None:
    conn = init_db(tmp_path / "telemetry.db")
    state_path = tmp_path / "scheduler_state.json"
    scheduler = FixedIntervalAutonomyScheduler(
        conn,
        config=_basic_config(),
        source_file=tmp_path / "scanner_rows.jsonl",
        state_path=state_path,
        utc_now_func=_fixed_clock("2026-04-14T00:00:00+00:00"),
        monotonic_func=iter([0.0, 0.1, 0.2]).__next__,
    )

    result = scheduler.tick(scanner_rows=_scanner_rows())

    snapshot = cast(dict[str, object], json.loads(state_path.read_text(encoding="utf-8")))
    assert result["status"] == "success"
    assert result["tick_index"] == 1
    assert result["proposal_count"] == 2
    assert result["dispatch_count"] == 1
    assert len(cast(list[object], result["dispatch_records"])) == 1
    assert snapshot["last_status"] == "success"
    assert snapshot["tick_index"] == 1
    assert snapshot["graph_counts"] == {"node_count": 6, "edge_count": 5}
    assert cast(dict[str, object], snapshot["graph_cursor"])["last_line_offset"] == 2
    assert cast(dict[str, object], snapshot["budgets"]) == {
        "dispatch_mode": "dry_run",
        "interval_seconds": 900.0,
        "max_dispatches": 1,
        "max_proposals": 2,
        "time_cap_seconds": None,
    }
    proposals = cast(list[dict[str, object]], snapshot["proposals"])
    assert len(proposals) == 2
    dispatch_result = cast(dict[str, object], snapshot["dispatch_result"])
    assert dispatch_result["dispatch_count"] == 1
    scheduled_goals = cast(list[dict[str, object]], snapshot["scheduled_goals"])
    assert len(scheduled_goals) == 1
    assert scheduled_goals[0]["goal_id"] == proposals[0]["goal_id"]


def test_task10_scheduler_rejects_overlapping_tick_without_concurrent_mutation(tmp_path: Path) -> None:
    db_path = tmp_path / "telemetry.db"
    bootstrap_conn = init_db(db_path)
    bootstrap_conn.close()
    conn = sqlite3.connect(str(db_path), check_same_thread=False)
    state_path = tmp_path / "scheduler_state.json"
    adapter = BlockingAdapter()
    scheduler = FixedIntervalAutonomyScheduler(
        conn,
        config=_basic_config(),
        source_file=tmp_path / "scanner_rows.jsonl",
        state_path=state_path,
        adapter=adapter,
        utc_now_func=_fixed_clock("2026-04-14T00:10:04+00:00"),
        monotonic_func=iter([0.0, 0.1, 0.2]).__next__,
    )

    first_result: dict[str, object] = {}

    def run_first_tick() -> None:
        first_result.update(scheduler.tick(scanner_rows=_scanner_rows()))

    worker = threading.Thread(target=run_first_tick)
    worker.start()
    assert adapter.entered.wait(timeout=5)

    conflict = scheduler.tick(scanner_rows=_scanner_rows())
    assert conflict == {
        "status": "lock_conflict",
        "requested_at": "2026-04-14T00:10:04+00:00",
        "tick_index": 0,
        "mutated": False,
    }
    assert not state_path.exists()

    adapter.release.set()
    worker.join(timeout=5)
    assert not worker.is_alive()
    assert first_result["status"] == "success"

    snapshot = cast(dict[str, object], json.loads(state_path.read_text(encoding="utf-8")))
    assert snapshot["tick_index"] == 1
    assert snapshot["last_status"] == "success"
