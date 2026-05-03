from __future__ import annotations
# pyright: reportMissingImports=false, reportUnknownVariableType=false, reportUnannotatedClassAttribute=false, reportIndexIssue=false, reportUnknownArgumentType=false, reportUnusedCallResult=false, reportUnknownMemberType=false, reportUnusedParameter=false, reportAny=false, reportArgumentType=false

import json
import os
import sys
import sqlite3
from pathlib import Path
from typing import cast

import pytest

from autoresearch.harness.autonomy_scheduler import FixedIntervalAutonomyScheduler
from autoresearch.harness.campaign_supervisor import CampaignSupervisor
from autoresearch.harness.dihiggs_adaptation import BanditState
from autoresearch.harness.mvp_graph import MvpGraphState
from autoresearch.harness.status_dashboard import print_cli_status
from autoresearch.harness.telemetry_store import db_counts, init_db


def _base_supervisor_config(tmp_path: Path) -> dict[str, object]:
    return {
        "campaign_id": "supervisor-e2e-smoke",
        "paths": {
            "repo_root": str(tmp_path),
            "outdir": str(tmp_path / "campaign"),
            "lake_name": "events.jsonl",
        },
        "runtime": {
            "python_exe": sys.executable,
            "threads": 1,
            "timeout_sec": 30,
            "max_empty_rounds": 2,
        },
        "dihiggs": {
            "phys_exec": str(tmp_path / "fake-phys"),
            "hb_dataset_env": "HB_DATASET",
            "hs_dataset_env": "HS_DATASET",
        },
        "supervisor": {
            "campaign_id": "single-host-smoke",
            "enable_autoscaling": True,
            "enable_alerting": True,
            "enable_convergence": True,
            "max_rounds": 2,
            "rebuild_telemetry": True,
        },
        "convergence": {
            "warmup_rounds": 1,
            "window_rounds": 1,
            "confirmation_required": 1,
            "thresholds": {
                "coverage_delta": 0.01,
                "diversity_delta": 0.01,
                "yield_delta": 0.01,
                "composite_delta": 0.01,
            },
        },
        "alerts": {
            "dedupe_window_sec": 900,
            "channels": ["stderr", "file"],
            "severity_levels": ["CRITICAL", "WARNING", "INFO"],
        },
        "autoscaling": {
            "min_threads": 1,
            "max_threads": 2,
            "cooldown_rounds": 1,
            "scale_up_factor": 1.5,
            "scale_down_factor": 0.5,
        },
        "search": {
            "tb_values": [50, 100],
            "lam1": {"min": -1.0, "max": 1.0, "n_bins": 4},
            "mphi": {"min": 130.0, "max": 130.0, "n_bins": 1},
        },
        "limits": {
            "max_new_run_dirs_per_round": 4,
            "max_new_run_dirs_per_arm_call": 2,
        },
        "metrics": {
            "floors": {"yield_norm": 0.01, "coverage_norm": 0.0, "diversity_norm": 0.0},
            "weights": {"yield": 0.4, "coverage": 0.3, "diversity": 0.3},
            "multi_axis": {
                "enabled": True,
                "collapse_axes": ["tb", "lam1_bin"],
                "coverage_axes": [
                    {"name": "tb", "kind": "categorical", "domain": [50, 100], "weight": 0.5},
                    {"name": "lam1_bin", "kind": "discrete", "domain_size": 4, "weight": 0.5},
                ],
                "diversity_axes": [
                    {"name": "tb", "weight": 0.5},
                    {"name": "lam1_bin", "weight": 0.5},
                ],
                "diversity_pairs": [],
            },
        },
        "adaptation": {
            "ucb_c": 1.414,
            "warm_start_each_arm": False,
        },
        "arms": [{"id": "adaptive-smoke"}],
    }


class FakeRunner:
    plans: list[dict[str, object]] = []

    def __init__(self, config: dict[str, object], outdir: str):
        self.config = config
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.event_log_path = self.outdir / str(config["paths"]["lake_name"])
        self.bandit_state = BanditState(
            arm_stats={"adaptive-smoke": {"pulls": 0, "total_reward": 0.0}},
            global_history=[],
            arm_histories={"adaptive-smoke": []},
        )
        self.calls = 0

    def run_adaptation_round(self) -> dict[str, object]:
        plan = type(self).plans[self.calls]
        self.calls += 1
        attempt = plan.get("attempt")
        if isinstance(attempt, dict):
            self.bandit_state.arm_stats["adaptive-smoke"]["pulls"] += 1
            self.bandit_state.global_history.append(attempt)
            self.bandit_state.arm_histories["adaptive-smoke"].append(attempt)
            event = {
                "schema_version": 1,
                "campaign_id": self.config["campaign_id"],
                "event_type": "ATTEMPT_EVALUATED",
                "utc": plan.get("utc", "2026-04-06T00:00:00Z"),
                "payload": {
                    "arm_id": "adaptive-smoke",
                    "attempt_id": f"attempt-{self.calls}",
                    "iter_index": self.calls - 1,
                    "attempt_index": 0,
                    "cell_id": f"tb={attempt['axes_binned']['tb']}|lam1_bin={attempt['axes_binned']['lam1_bin']}",
                    "eval_status": "success",
                    "successes": attempt.get("yield", 0.0),
                    "trials": 1,
                    "elapsed_sec": 0.1,
                    "axes_binned": attempt["axes_binned"],
                    "axes_raw": attempt["axes_raw"],
                },
            }
            with self.event_log_path.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(event) + "\n")
        return {
            "arm_id": "adaptive-smoke",
            "discoveries": plan.get("discoveries", 0),
            "events_emitted": plan.get("events_emitted", 0),
            "subprocess_status": plan.get("subprocess_status", "success"),
            "utc": plan.get("utc", "2026-04-06T00:00:00Z"),
        }


class StubConvergence:
    states: list[str] = []

    def __init__(self, config: dict[str, object]):
        self.calls = 0

    def evaluate(self, snapshot_history: list[dict[str, object]], arm_pulls: dict[str, object]) -> dict[str, object]:
        state = type(self).states[self.calls]
        self.calls += 1
        return {"state": state, "reason": state}


def _preflight_pass(config: dict[str, object]) -> dict[str, object]:
    return {"overall": "pass", "checks": []}


def _load_lam1_supervisor_config(tmp_path: Path) -> dict[str, object]:
    repo_root = Path(__file__).resolve().parents[2]
    config_path = repo_root / "autoresearch" / "configs" / "dihiggs_explorers_lam1.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    config["paths"]["repo_root"] = str(repo_root)
    config["paths"]["outdir"] = str(tmp_path / "campaign")
    config["runtime"]["python_exe"] = sys.executable
    config["supervisor"]["max_rounds"] = 1
    config["supervisor"]["rebuild_telemetry"] = True
    return config


def _task12_config() -> dict[str, object]:
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


def _task12_scanner_rows() -> list[dict[str, object]]:
    return [
        {
            "campaign_id": "camp-e2e-001",
            "arm_id": "adaptive-v1",
            "run_dir_fingerprint": "run-001",
            "attempt_id": "attempt-001",
            "artifact_status": "complete",
            "discovery_status": "new",
            "axes_binned": {"tb": 1000, "lam1_bin": 0},
        },
        {
            "campaign_id": "camp-e2e-001",
            "arm_id": "adaptive-v1",
            "run_dir_fingerprint": "run-002",
            "attempt_id": "attempt-002",
            "artifact_status": "complete",
            "discovery_status": "new",
            "axes_binned": {"tb": 1000, "lam1_bin": 1},
        },
        {
            "campaign_id": "camp-e2e-001",
            "arm_id": "adaptive-v1",
            "run_dir_fingerprint": "run-003",
            "attempt_id": "attempt-003",
            "artifact_status": "complete",
            "discovery_status": "new",
            "axes_binned": {"tb": 5000, "lam1_bin": 0},
        },
    ]


def _task12_fetch_payload_rows(conn: sqlite3.Connection, table: str) -> list[dict[str, object]]:
    cursor = conn.cursor()
    _ = cursor.execute(f"SELECT payload FROM {table} ORDER BY source_file, line_offset, id")
    rows = cast(list[tuple[object]], cursor.fetchall())
    return [cast(dict[str, object], json.loads(str(row[0]))) for row in rows]


def _task12_graph_payload(graph: MvpGraphState) -> dict[str, object]:
    return {
        "nodes": {"|".join(key): graph.nodes[key] for key in sorted(graph.nodes)},
        "edges": {"|".join(key): graph.edges[key] for key in sorted(graph.edges)},
    }


def _normalize_task12_value(value: object, *, replacements: dict[str, str]) -> object:
    if isinstance(value, dict):
        return {key: _normalize_task12_value(item, replacements=replacements) for key, item in value.items()}
    if isinstance(value, list):
        return [_normalize_task12_value(item, replacements=replacements) for item in value]
    if isinstance(value, str):
        normalized = value
        for original, replacement in replacements.items():
            normalized = normalized.replace(original, replacement)
        return normalized
    return value


def _task12_evidence_path(name: str) -> Path:
    return Path(__file__).resolve().parents[2] / ".sisyphus" / "evidence" / name


def _write_task12_json_evidence(name: str, payload: object) -> None:
    path = _task12_evidence_path(name)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_task12_text_evidence(name: str, content: str) -> None:
    path = _task12_evidence_path(name)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


class _FixedDryRunAdapter:
    def __init__(self, *, dispatched_at: str, timestamp: str) -> None:
        self._dispatched_at = dispatched_at
        self._timestamp = timestamp

    def dispatch(self, payload: dict[str, object]) -> dict[str, object]:
        proposals = cast(list[dict[str, object]], payload.get("proposals", []))
        records = [
            {
                "goal_id": str(proposal.get("goal_id")),
                "status": "dry_run",
                "dispatched_at": self._dispatched_at,
                "metadata": {
                    "axes_binned": proposal.get("axes_binned"),
                    "rationale": proposal.get("rationale"),
                },
            }
            for proposal in proposals
        ]
        return {
            "status": "success",
            "dispatch_count": len(records),
            "records": records,
            "timestamp": self._timestamp,
        }


def _run_task12_scheduler_flow(tmp_path: Path) -> tuple[dict[str, object], str]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    conn = init_db(tmp_path / "telemetry.db")
    state_path = tmp_path / "scheduler_state.json"
    source_file = tmp_path / "scanner_rows.jsonl"
    graph = MvpGraphState()
    scheduler = FixedIntervalAutonomyScheduler(
        conn,
        config=_task12_config(),
        source_file=source_file,
        state_path=state_path,
        graph=graph,
        adapter=_FixedDryRunAdapter(
            dispatched_at="2026-04-14T12:00:00+00:00",
            timestamp="2026-04-14T12:00:00+00:00",
        ),
        utc_now_func=lambda: "2026-04-14T12:00:00+00:00",
        monotonic_func=iter([0.0, 0.25]).__next__,
    )

    tick_result = scheduler.tick(scanner_rows=_task12_scanner_rows())
    snapshot = cast(dict[str, object], json.loads(state_path.read_text(encoding="utf-8")))
    coverage_state = _task12_fetch_payload_rows(conn, "coverage_state")
    discovery_records = _task12_fetch_payload_rows(conn, "discovery_records")
    run_registry = _task12_fetch_payload_rows(conn, "run_registry")
    counts = db_counts(conn)
    conn.close()

    first_proposal = cast(list[dict[str, object]], snapshot["proposals"])[0]
    first_rationale = cast(dict[str, object], first_proposal["rationale"])
    stage_markers = {
        "reconcile": "pass" if cast(dict[str, object], tick_result["upsert_result"])["rows_valid"] == 3 else "fail",
        "graph_update": "pass" if cast(dict[str, object], snapshot["graph_summary"])["inserted"] == 17 else "fail",
        "propose": "pass" if tick_result["proposal_count"] == 2 else "fail",
        "dry_run_dispatch": "pass" if tick_result["dispatch_count"] == 1 else "fail",
        "scheduler_tick": "pass" if tick_result["status"] == "success" else "fail",
    }
    summary = {
        "fixture": "task12-deterministic-e2e",
        "seed": 0,
        "clock": "2026-04-14T12:00:00+00:00",
        "stage_markers": stage_markers,
        "tick_result": tick_result,
        "snapshot": snapshot,
        "db_counts": counts,
        "run_registry": run_registry,
        "discovery_records": discovery_records,
        "coverage_state": coverage_state,
        "graph": _task12_graph_payload(graph),
        "assertions": {
            "first_proposal_goal_id": first_proposal["goal_id"],
            "first_proposal_axes": first_proposal["axes_binned"],
            "first_proposal_signals": first_rationale["signals"],
            "first_proposal_source_refs": first_rationale["source_refs"],
        },
    }
    normalized_summary = cast(
        dict[str, object],
        _normalize_task12_value(
            summary,
            replacements={
                str(source_file): "FIXTURE_SOURCE_FILE",
                str(state_path): "FIXTURE_STATE_PATH",
                str(tmp_path): "FIXTURE_TMPDIR",
            },
        ),
    )
    stable_summary = json.dumps(normalized_summary, indent=2, sort_keys=True) + "\n"
    return normalized_summary, stable_summary


def test_supervisor_smoke_creates_operator_artifacts(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    config = _base_supervisor_config(tmp_path)
    outdir = Path(config["paths"]["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "telemetry.db").write_text("stale telemetry", encoding="utf-8")

    FakeRunner.plans = [
        {
            "discoveries": 1,
            "events_emitted": 1,
            "utc": "2026-04-06T00:00:00Z",
            "attempt": {
                "yield": 0.6,
                "axes_binned": {"tb": 50, "lam1_bin": 1},
                "axes_raw": {"tb": 50, "lam1": -0.5},
            },
        },
        {
            "discoveries": 1,
            "events_emitted": 1,
            "utc": "2026-04-06T00:05:00Z",
            "attempt": {
                "yield": 0.7,
                "axes_binned": {"tb": 100, "lam1_bin": 2},
                "axes_raw": {"tb": 100, "lam1": 0.25},
            },
        },
    ]
    StubConvergence.states = ["plateau_warning", "plateau_warning"]

    supervisor = CampaignSupervisor(
        config,
        runner_factory=FakeRunner,
        preflight_runner=_preflight_pass,
        convergence_factory=StubConvergence,
    )

    terminal = supervisor.run()

    manifest_path = outdir / "campaign_manifest.json"
    state_path = outdir / "campaign_state.json"
    status_path = outdir / "campaign_status.json"
    alerts_path = outdir / "alerts.jsonl"
    supervisor_events_path = outdir / "supervisor_events.jsonl"

    assert terminal == "COMPLETED"
    assert manifest_path.exists()
    assert state_path.exists()
    assert status_path.exists()
    assert alerts_path.exists()
    assert supervisor_events_path.exists()

    state = json.loads(state_path.read_text(encoding="utf-8"))
    status = json.loads(status_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    alerts = [json.loads(line) for line in alerts_path.read_text(encoding="utf-8").splitlines()]
    supervisor_events = [json.loads(line) for line in supervisor_events_path.read_text(encoding="utf-8").splitlines()]

    conn = init_db(outdir / "telemetry.db")
    counts = db_counts(conn)
    conn.close()

    assert state["campaign_state"] == "COMPLETED"
    assert state["round_index"] == 1
    assert state["preflight_overall"] == "pass"
    assert state["convergence_state"] == "plateau_warning"
    assert state["plateau_warning"] is True
    assert len(state["active_alerts"]) == 1

    assert status["campaign_state"] == "RUNNING"
    assert status["round_index"] == 1
    assert status["selected_arm"] == "adaptive-smoke"
    assert len(status["active_alerts"]) == 1

    assert manifest["supervisor"]["enable_alerting"] is True
    assert manifest["alerts"]["channels"] == ["stderr", "file"]
    assert manifest["autoscaling"]["min_threads"] == 1
    assert counts["events"] == 4

    assert len(supervisor_events) == 2
    assert supervisor_events[-1]["event_type"] == "ROUND_COMPLETED"
    assert alerts[0]["alert_type"] == "plateau_warning"

    print_cli_status(outdir)
    captured = capsys.readouterr()
    assert "Campaign Status: RUNNING" in captured.out
    assert "Active Round:    1" in captured.out
    assert "Alerts:          1 active alerts" in captured.out
    assert "[INFO] plateau_warning" in captured.out
    assert "[INFO] plateau_warning" in captured.err


@pytest.mark.skipif(
    os.environ.get("DIHIGGS_E2E_TEST") != "1",
    reason="Supervisor E2E requires DIHIGGS_E2E_TEST=1 and real prerequisites",
)
def test_supervisor_smoke_gated_real_backend(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    config = _load_lam1_supervisor_config(tmp_path)
    outdir = Path(config["paths"]["outdir"])

    supervisor = CampaignSupervisor(config, outdir=outdir)
    terminal = supervisor.run()

    assert terminal in {"COMPLETED", "CONVERGED"}
    assert (outdir / "campaign_manifest.json").exists()
    assert (outdir / "campaign_state.json").exists()
    assert (outdir / "campaign_status.json").exists()
    assert (outdir / "supervisor_events.jsonl").exists()
    assert (outdir / "telemetry.db").exists()

    print_cli_status(outdir)
    captured = capsys.readouterr()
    assert "Campaign Status:" in captured.out
    assert "Metrics:" in captured.out


def test_task12_deterministic_e2e_dry_run_flow(tmp_path: Path) -> None:
    summary, _stable_summary = _run_task12_scheduler_flow(tmp_path)

    assert summary["stage_markers"] == {
        "reconcile": "pass",
        "graph_update": "pass",
        "propose": "pass",
        "dry_run_dispatch": "pass",
        "scheduler_tick": "pass",
    }
    tick_result = cast(dict[str, object], summary["tick_result"])
    snapshot = cast(dict[str, object], summary["snapshot"])
    graph_summary = cast(dict[str, object], snapshot["graph_summary"])
    dispatch_result = cast(dict[str, object], snapshot["dispatch_result"])
    proposals = cast(list[dict[str, object]], snapshot["proposals"])
    scheduled_goals = cast(list[dict[str, object]], snapshot["scheduled_goals"])
    first_proposal = proposals[0]
    first_rationale = cast(dict[str, object], first_proposal["rationale"])

    assert tick_result["status"] == "success"
    assert tick_result["tick_index"] == 1
    assert tick_result["mutated"] is True
    assert tick_result["proposal_count"] == 2
    assert tick_result["dispatch_count"] == 1
    assert graph_summary == {
        "inserted": 17,
        "updated": 0,
        "noops": 0,
        "node_count": 9,
        "edge_count": 8,
    }
    assert cast(dict[str, object], snapshot["graph_cursor"])["last_line_offset"] == 3
    assert first_proposal["goal_id"] == "goal::tb=1000|lam1_bin=2"
    assert first_proposal["axes_binned"] == {"tb": 1000, "lam1_bin": 2}
    assert first_rationale["signals"] == ["coverage_gap", "discovery_frontier"]
    assert cast(dict[str, object], first_rationale["discovery_frontier"])["neighbor_cell_ids"] == ["tb=1000|lam1_bin=1"]
    assert len(cast(list[object], first_rationale["source_refs"])) >= 3
    assert dispatch_result["status"] == "success"
    assert dispatch_result["dispatch_count"] == 1
    assert scheduled_goals == [proposals[0]]
    assert summary["db_counts"] == {
        "events": 0,
        "ingestion_errors": 0,
        "ingestion_log": 0,
        "task_summaries": 0,
        "orchestrator_log_entries": 0,
        "artifact_metadata": 0,
        "alerts": 0,
        "run_registry": 3,
        "discovery_records": 3,
        "coverage_state": 3,
        "upsert_quarantine": 0,
        "reconcile_watermarks": 1,
    }

    _write_task12_json_evidence("task-12-e2e-dryrun.json", summary)


def test_task12_deterministic_e2e_replay_is_byte_stable(tmp_path: Path) -> None:
    first_summary, first_stable = _run_task12_scheduler_flow(tmp_path / "run-one")
    second_summary, second_stable = _run_task12_scheduler_flow(tmp_path / "run-two")

    assert first_summary == second_summary
    assert first_stable == second_stable

    evidence = (
        "run_one==run_two\n"
        f"bytes={len(first_stable.encode('utf-8'))}\n"
        f"sha_preview={first_stable[:160]}\n"
    )
    _write_task12_text_evidence("task-12-e2e-determinism.txt", evidence)
