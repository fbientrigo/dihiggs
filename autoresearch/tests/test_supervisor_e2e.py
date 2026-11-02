from __future__ import annotations
# pyright: reportMissingImports=false, reportUnknownVariableType=false, reportUnannotatedClassAttribute=false, reportIndexIssue=false, reportUnknownArgumentType=false, reportUnusedCallResult=false, reportUnknownMemberType=false, reportUnusedParameter=false, reportAny=false, reportArgumentType=false

import json
import os
import sys
from pathlib import Path

import pytest

from autoresearch.harness.campaign_supervisor import CampaignSupervisor
from autoresearch.harness.dihiggs_adaptation import BanditState
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




