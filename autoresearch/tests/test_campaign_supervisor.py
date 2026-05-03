from __future__ import annotations
# pyright: reportUnusedImport=false, reportUnannotatedClassAttribute=false, reportIndexIssue=false, reportUnknownArgumentType=false, reportUnusedCallResult=false, reportUnknownMemberType=false, reportUnusedParameter=false, reportPrivateUsage=false, reportArgumentType=false, reportAny=false, reportUnknownVariableType=false

import json
from pathlib import Path

from autoresearch.harness.campaign_supervisor import CampaignSupervisor
from autoresearch.harness.dihiggs_adaptation import BanditState
from autoresearch.harness.telemetry_store import db_counts, init_db


def _base_config(tmp_path: Path) -> dict[str, object]:
    return {
        "campaign_id": "campaign-supervisor-test",
        "paths": {
            "repo_root": str(tmp_path),
            "outdir": str(tmp_path / "out"),
            "lake_name": "events.jsonl",
        },
        "runtime": {
            "threads": 1,
            "timeout_sec": 30,
            "max_empty_rounds": 2,
        },
        "limits": {
            "max_new_run_dirs_per_round": 5,
            "max_new_run_dirs_per_arm_call": 3,
        },
        "metrics": {
            "weights": {"yield": 0.4, "coverage": 0.3, "diversity": 0.3},
            "multi_axis": {"collapse_axes": ["tb", "lam1_bin"]},
        },
        "dihiggs": {
            "phys_exec": str(tmp_path / "fake-phys"),
            "hb_dataset_env": "HB_DATASET",
            "hs_dataset_env": "HS_DATASET",
        },
        "supervisor": {"max_rounds": 1},
        "arms": [{"id": "adaptive-smoke"}],
    }


class FakeRunner:
    plans: list[dict[str, object]] = []

    def __init__(self, config: dict[str, object], outdir: str):
        self.config = config
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.event_log_path = self.outdir / str(config["paths"]["lake_name"])
        self.bandit_state = BanditState(arm_stats={"adaptive-smoke": {"pulls": 0, "total_reward": 0.0}}, global_history=[], arm_histories={"adaptive-smoke": []})
        self.calls = 0
        self.seen_threads: list[int] = []

    def run_adaptation_round(self) -> dict[str, object]:
        plan = type(self).plans[self.calls]
        self.calls += 1
        self.seen_threads.append(int(self.config["runtime"]["threads"]))
        payload = plan.get("attempt")
        if isinstance(payload, dict):
            self.bandit_state.arm_stats.setdefault("adaptive-smoke", {"pulls": 0, "total_reward": 0.0})
            self.bandit_state.arm_stats["adaptive-smoke"]["pulls"] += 1
            self.bandit_state.global_history.append(payload)
            self.bandit_state.arm_histories.setdefault("adaptive-smoke", []).append(payload)
            with self.event_log_path.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps({
                    "schema_version": 1,
                    "campaign_id": self.config["campaign_id"],
                    "event_type": "ATTEMPT_EVALUATED",
                    "utc": plan.get("utc", "2026-04-06T00:00:00Z"),
                    "payload": {
                        "arm_id": "adaptive-smoke",
                        "attempt_id": f"attempt-{self.calls}",
                        "cell_id": "tb=1|lam1_bin=1",
                        "successes": payload.get("yield", 0.0),
                        "axes_binned": payload.get("axes_binned", {}),
                        "axes_raw": payload.get("axes_raw", {}),
                    },
                }) + "\n")
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
        self.config = config
        self.calls = 0

    def evaluate(self, snapshot_history: list[dict[str, object]], arm_pulls: dict[str, object]) -> dict[str, object]:
        state = type(self).states[self.calls]
        self.calls += 1
        return {"state": state, "reason": state}


class StubAutoscaling:
    actions: list[dict[str, object]] = []

    def __init__(self, config: dict[str, object]):
        self.config = config
        self.calls = 0

    def evaluate(self, state: dict[str, object]) -> dict[str, object]:
        result = type(self).actions[self.calls]
        self.calls += 1
        return result


def _preflight_pass(config: dict[str, object]) -> dict[str, object]:
    return {"overall": "pass", "checks": []}


def test_transition_acceptance_states() -> None:
    supervisor = CampaignSupervisor.for_test(
        {
            "campaign_id": "mvp",
            "runtime": {"threads": 1, "timeout_sec": 10, "max_empty_rounds": 2},
            "arms": [{"id": "adaptive-smoke"}],
        }
    )

    assert supervisor._transition("RUNNING", {"convergence_state": "converged_plateau"}) == "CONVERGED"
    assert supervisor._transition("INIT", {"preflight_overall": "fail"}) == "PREFLIGHT_BLOCKED"


def test_run_blocks_before_rounds_on_preflight_failure(tmp_path: Path) -> None:
    FakeRunner.plans = [{"discoveries": 1, "events_emitted": 1}]
    supervisor = CampaignSupervisor(
        _base_config(tmp_path),
        runner_factory=FakeRunner,
        preflight_runner=lambda config: {"overall": "fail", "checks": [{"check": "phys_exec", "status": "fail"}]},
    )

    terminal = supervisor.run()

    state = json.loads((Path(supervisor.outdir) / "campaign_state.json").read_text(encoding="utf-8"))
    assert terminal == "PREFLIGHT_BLOCKED"
    assert state["campaign_state"] == "PREFLIGHT_BLOCKED"
    assert supervisor.runner.calls == 0


def test_run_initializes_telemetry_and_persists_completed_state(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    FakeRunner.plans = [
        {
            "discoveries": 1,
            "events_emitted": 1,
            "utc": "2026-04-06T00:00:00Z",
            "attempt": {
                "yield": 0.5,
                "axes_binned": {"tb": 1, "lam1_bin": 1},
                "axes_raw": {"tb": 1, "lam1": 0.1},
            },
        }
    ]
    StubConvergence.states = ["running"]
    StubAutoscaling.actions = [{
        "action": "hold",
        "last_scaling_round": None,
        "next_threads": 1,
        "next_timeout_sec": 30,
        "next_max_new_run_dirs_per_round": 5,
        "next_max_new_run_dirs_per_arm_call": 3,
    }]
    supervisor = CampaignSupervisor(
        config,
        runner_factory=FakeRunner,
        preflight_runner=_preflight_pass,
        convergence_factory=StubConvergence,
        autoscaling_factory=StubAutoscaling,
    )

    terminal = supervisor.run()

    state = json.loads((Path(supervisor.outdir) / "campaign_state.json").read_text(encoding="utf-8"))
    status = json.loads((Path(supervisor.outdir) / "campaign_status.json").read_text(encoding="utf-8"))
    conn = init_db(Path(supervisor.outdir) / "telemetry.db")
    counts = db_counts(conn)
    conn.close()

    assert terminal == "COMPLETED"
    assert (Path(supervisor.outdir) / "campaign_status.html").exists()
    html_content = (Path(supervisor.outdir) / "campaign_status.html").read_text(encoding="utf-8")
    assert "COMPLETED" in html_content
    assert "last_successful_round" in status
    assert "manifest_compatibility" in status
    assert "last_scaling_action" in status
    assert "file_timestamps" in status
    assert state["campaign_state"] == "COMPLETED"
    assert state["round_index"] == 0
    assert (Path(supervisor.outdir) / "campaign_manifest.json").exists()
    assert counts["events"] == 2
    assert status["round_index"] == 0


def test_run_transitions_to_converged_between_rounds(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    config["supervisor"]["max_rounds"] = 3
    FakeRunner.plans = [
        {"discoveries": 1, "events_emitted": 1, "attempt": {"yield": 0.4, "axes_binned": {"tb": 1, "lam1_bin": 1}, "axes_raw": {"tb": 1, "lam1": 0.1}}},
        {"discoveries": 0, "events_emitted": 0},
    ]
    StubConvergence.states = ["running", "converged_plateau"]
    StubAutoscaling.actions = [
        {"action": "hold", "last_scaling_round": None, "next_threads": 1, "next_timeout_sec": 30, "next_max_new_run_dirs_per_round": 5, "next_max_new_run_dirs_per_arm_call": 3},
        {"action": "hold", "last_scaling_round": None, "next_threads": 1, "next_timeout_sec": 30, "next_max_new_run_dirs_per_round": 5, "next_max_new_run_dirs_per_arm_call": 3},
    ]
    supervisor = CampaignSupervisor(
        config,
        runner_factory=FakeRunner,
        preflight_runner=_preflight_pass,
        convergence_factory=StubConvergence,
        autoscaling_factory=StubAutoscaling,
    )

    terminal = supervisor.run()

    assert terminal == "CONVERGED"
    assert supervisor.runner.calls == 2


def test_autoscaling_applies_only_between_rounds(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    config["supervisor"]["max_rounds"] = 2
    FakeRunner.plans = [
        {"discoveries": 1, "events_emitted": 1, "attempt": {"yield": 0.4, "axes_binned": {"tb": 1, "lam1_bin": 1}, "axes_raw": {"tb": 1, "lam1": 0.1}}},
        {"discoveries": 1, "events_emitted": 1, "attempt": {"yield": 0.6, "axes_binned": {"tb": 2, "lam1_bin": 2}, "axes_raw": {"tb": 2, "lam1": 0.2}}},
    ]
    StubConvergence.states = ["running", "running"]
    StubAutoscaling.actions = [
        {"action": "scale_up", "last_scaling_round": 0, "next_threads": 2, "next_timeout_sec": 40, "next_max_new_run_dirs_per_round": 6, "next_max_new_run_dirs_per_arm_call": 4},
        {"action": "hold", "last_scaling_round": 0, "next_threads": 2, "next_timeout_sec": 40, "next_max_new_run_dirs_per_round": 6, "next_max_new_run_dirs_per_arm_call": 4},
    ]
    supervisor = CampaignSupervisor(
        config,
        runner_factory=FakeRunner,
        preflight_runner=_preflight_pass,
        convergence_factory=StubConvergence,
        autoscaling_factory=StubAutoscaling,
    )

    terminal = supervisor.run()

    assert terminal == "COMPLETED"
    assert supervisor.runner.seen_threads == [1, 2]


def test_supervisor_replays_expanded_raw_telemetry_sources(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    outdir = Path(config["paths"]["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)

    stdout_path = outdir / "task-001.stdout.txt"
    stderr_path = outdir / "task-001.stderr.txt"
    stdout_path.write_text("stdout artifact\n", encoding="utf-8")
    stderr_path.write_text("stderr artifact\n", encoding="utf-8")

    task_summary_record = {
        "campaign_id": config["campaign_id"],
        "task_id": "task-001",
        "task_type": "adaptive",
        "status": "completed",
        "payload": {
            "stdout_path": str(stdout_path),
            "stderr_path": str(stderr_path),
        },
    }
    (outdir / "task_summary.jsonl").write_text(json.dumps(task_summary_record) + "\n", encoding="utf-8")
    (outdir / "orchestrator.log").write_text(
        "[2026-04-06T00:00:00Z] [INFO] supervisor replay bootstrap\n",
        encoding="utf-8",
    )

    FakeRunner.plans = [
        {
            "discoveries": 1,
            "events_emitted": 1,
            "utc": "2026-04-06T00:00:00Z",
            "attempt": {
                "yield": 0.5,
                "axes_binned": {"tb": 1, "lam1_bin": 1},
                "axes_raw": {"tb": 1, "lam1": 0.1},
            },
        }
    ]
    StubConvergence.states = ["running"]
    StubAutoscaling.actions = [{
        "action": "hold",
        "last_scaling_round": None,
        "next_threads": 1,
        "next_timeout_sec": 30,
        "next_max_new_run_dirs_per_round": 5,
        "next_max_new_run_dirs_per_arm_call": 3,
    }]

    supervisor = CampaignSupervisor(
        config,
        runner_factory=FakeRunner,
        preflight_runner=_preflight_pass,
        convergence_factory=StubConvergence,
        autoscaling_factory=StubAutoscaling,
    )

    terminal = supervisor.run()

    conn = init_db(Path(supervisor.outdir) / "telemetry.db")
    counts = db_counts(conn)
    cursor = conn.cursor()
    cursor.execute("SELECT task_id, task_type, status FROM task_summaries")
    task_summary_row = cursor.fetchone()
    cursor.execute("SELECT log_level, log_message FROM orchestrator_log_entries")
    orchestrator_row = cursor.fetchone()
    cursor.execute("SELECT artifact_type, related_task_id FROM artifact_metadata ORDER BY artifact_type")
    artifact_rows = cursor.fetchall()
    conn.close()

    assert terminal == "COMPLETED"
    assert counts["task_summaries"] == 1
    assert counts["orchestrator_log_entries"] == 1
    assert counts["artifact_metadata"] == 2
    assert task_summary_row == ("task-001", "adaptive", "completed")
    assert orchestrator_row == ("INFO", "supervisor replay bootstrap")
    assert artifact_rows == [("stderr", "task-001"), ("stdout", "task-001")]
