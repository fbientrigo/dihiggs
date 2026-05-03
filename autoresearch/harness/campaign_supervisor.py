from __future__ import annotations
# pyright: reportUnusedImport=false, reportPrivateLocalImportUsage=false, reportUnknownArgumentType=false, reportExplicitAny=false, reportAny=false, reportUnannotatedClassAttribute=false, reportArgumentType=false, reportIndexIssue=false, reportUnusedCallResult=false, reportUnknownMemberType=false, reportUnusedParameter=false, reportPrivateUsage=false, reportUnknownVariableType=false

import json
import time
import tempfile
import tempfile
from collections.abc import Callable, Mapping
from contextlib import suppress
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from autoresearch.benchmarks.score import compute_composite, compute_coverage, compute_diversity
from autoresearch.harness.alerting import AlertEngine
from autoresearch.harness.autoscaling import AutoscalingPolicy
from autoresearch.harness.campaign_manifest import validate_manifest, write_manifest
from autoresearch.harness.convergence import ConvergenceDetector
from autoresearch.harness.dihiggs_preflight import run_all_preflight_checks
from autoresearch.harness.dihiggs_runner import DiHiggsRunner
from autoresearch.harness.round_snapshots import project_round_snapshot
from autoresearch.harness.status_dashboard import write_status_html, write_status_json
from autoresearch.harness.telemetry_ingest import ingest_artifact_metadata
from autoresearch.harness.telemetry_store import init_db, replay_sources

TERMINAL_STATES = {"COMPLETED", "CONVERGED", "FAILED", "PREFLIGHT_BLOCKED"}


def _as_dict(value: object) -> dict[str, object]:
    return dict(value) if isinstance(value, Mapping) else {}


def _as_int(value: object, default: int = 0) -> int:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    if isinstance(value, str):
        with suppress(ValueError):
            return int(value)
    return default


def _as_float(value: object, default: float = 0.0) -> float:
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        with suppress(ValueError):
            return float(value)
    return default


def _last_yield(history: list[dict[str, Any]]) -> float:
    if not history:
        return 0.0
    return _as_float(history[-1].get("yield", 0.0), 0.0)


@dataclass
class _Dependencies:
    runner_factory: Callable[[dict[str, object], str], Any]
    preflight_runner: Callable[[dict[str, object]], dict[str, object]]
    manifest_writer: Callable[[str | Path, Mapping[str, object]], Path]
    manifest_validator: Callable[[str, Mapping[str, object]], dict[str, str]]
    telemetry_init: Callable[[str | Path], Any]
    telemetry_replay: Callable[[Any, list[str | Path], str], list[dict[str, Any]]]
    convergence_factory: Callable[[Mapping[str, object]], Any]
    alert_engine_factory: Callable[[dict[str, object]], Any]
    autoscaling_factory: Callable[[Mapping[str, object]], Any]


class CampaignSupervisor:
    def __init__(
        self,
        config: dict[str, object],
        outdir: str | Path | None = None,
        *,
        runner_factory: Callable[[dict[str, object], str], Any] = DiHiggsRunner,
        preflight_runner: Callable[[dict[str, object]], dict[str, object]] = run_all_preflight_checks,
        manifest_writer: Callable[[str | Path, Mapping[str, object]], Path] = write_manifest,
        manifest_validator: Callable[[str, Mapping[str, object]], dict[str, str]] = validate_manifest,
        telemetry_init: Callable[[str | Path], Any] = init_db,
        telemetry_replay: Callable[[Any, list[str | Path]], list[dict[str, Any]]] = replay_sources,
        convergence_factory: Callable[[Mapping[str, object]], Any] = ConvergenceDetector,
        alert_engine_factory: Callable[[dict[str, object]], Any] = AlertEngine,
        autoscaling_factory: Callable[[Mapping[str, object]], Any] = AutoscalingPolicy,
    ):
        self.config = config
        paths = _as_dict(config.get("paths"))
        resolved_outdir = Path(outdir or paths.get("outdir") or tempfile.mkdtemp(prefix="campaign-supervisor-"))
        resolved_outdir.mkdir(parents=True, exist_ok=True)
        paths.setdefault("outdir", str(resolved_outdir))
        paths.setdefault("repo_root", str(Path.cwd()))
        paths.setdefault("lake_name", "events.jsonl")
        self.config["paths"] = paths

        self.outdir = resolved_outdir
        self.state_path = self.outdir / "campaign_state.json"
        self.telemetry_db_path = self.outdir / "telemetry.db"
        self.manifest_path = self.outdir / "campaign_manifest.json"
        self.events_path = self.outdir / str(paths.get("lake_name", "events.jsonl"))
        self.supervisor_events_path = self.outdir / "supervisor_events.jsonl"
        self.task_summary_path = self.outdir / "task_summary.jsonl"
        self.orchestrator_log_path = self.outdir / "orchestrator.log"

        self.dependencies = _Dependencies(
            runner_factory=runner_factory,
            preflight_runner=preflight_runner,
            manifest_writer=manifest_writer,
            manifest_validator=manifest_validator,
            telemetry_init=telemetry_init,
            telemetry_replay=telemetry_replay,
            convergence_factory=convergence_factory,
            alert_engine_factory=alert_engine_factory,
            autoscaling_factory=autoscaling_factory,
        )
        self.alert_engine = alert_engine_factory(_as_dict(self.config.get("alerts")))
        self.autoscaling_policy = autoscaling_factory(_as_dict(self.config.get("autoscaling")))
        self.convergence_detector = convergence_factory(self.config)
        self.runner = runner_factory(self.config, str(self.outdir))

    @classmethod
    def for_test(cls, config: dict[str, object]) -> "CampaignSupervisor":
        cfg = json.loads(json.dumps(config))
        paths = _as_dict(cfg.get("paths"))
        paths.setdefault("repo_root", str(Path.cwd()))
        paths.setdefault("outdir", tempfile.mkdtemp(prefix="campaign-supervisor-test-"))
        paths.setdefault("lake_name", "events.jsonl")
        cfg["paths"] = paths
        return cls(cfg, outdir=paths["outdir"])

    def _transition(self, current_state: str, signals: Mapping[str, object]) -> str:
        if current_state in TERMINAL_STATES:
            return current_state
        if signals.get("manifest_status") == "fail":
            return "FAILED"
        if signals.get("preflight_overall") == "fail":
            return "PREFLIGHT_BLOCKED"
        if signals.get("round_exception") is True:
            return "FAILED"
        convergence_state = signals.get("convergence_state")
        if isinstance(convergence_state, str) and convergence_state.startswith("converged"):
            return "CONVERGED"
        if signals.get("round_limit_reached") is True:
            return "COMPLETED"
        if signals.get("time_budget_reached") is True:
            return "COMPLETED"
        if current_state == "INIT":
            return "RUNNING"
        return current_state

    def run(self) -> str:
        state = self._load_state()
        current_state = str(state.get("campaign_state", "INIT"))
        if current_state in TERMINAL_STATES:
            return current_state

        # Initialize start time for time-budget tracking
        if "campaign_start_time" not in state:
            state["campaign_start_time"] = time.time()
            self._persist_state(state)

        manifest_status = self._ensure_manifest()
        state["manifest_status"] = manifest_status["status"]
        state["manifest_reason"] = manifest_status["reason"]
        state["manifest_path"] = str(self.manifest_path)
        current_state = self._transition(current_state, {"manifest_status": manifest_status["status"]})
        state["campaign_state"] = current_state
        self._persist_state(state)
        if current_state in TERMINAL_STATES:
            return current_state

        preflight = self.dependencies.preflight_runner(self.config)
        state["preflight_overall"] = preflight.get("overall", "fail")
        state["preflight_checks"] = preflight.get("checks", [])
        current_state = self._transition(current_state, {"preflight_overall": state["preflight_overall"]})
        state["campaign_state"] = current_state
        self._persist_state(state)
        if current_state in TERMINAL_STATES:
            return current_state

        replay_results = self._initialize_telemetry(rebuild=_as_dict(self.config.get("supervisor")).get("rebuild_telemetry") is True)
        state["telemetry_corruption"] = any(
            result.get("status") == "fatal_error" or _as_int(result.get("lines_errored", 0)) > 0
            for result in replay_results
        )
        self._persist_state(state)

        max_rounds = self._max_rounds()
        
        # Get time budget configuration
        supervisor_cfg = _as_dict(self.config.get("supervisor", {}))
        max_duration_hours = _as_float(supervisor_cfg.get("max_duration_hours"), None)
        stop_at_timestamp = _as_float(supervisor_cfg.get("stop_at_timestamp"), None)
        
        while current_state == "RUNNING":
            # Check time budget
            start_time = _as_float(state.get("campaign_start_time"), time.time())
            elapsed_hours = (time.time() - start_time) / 3600.0
            state["elapsed_hours"] = elapsed_hours
            
            time_budget_exceeded = False
            if max_duration_hours is not None and elapsed_hours >= max_duration_hours:
                time_budget_exceeded = True
                state["stop_reason"] = f"time_budget_exceeded (max_duration_hours={max_duration_hours}, elapsed={elapsed_hours:.2f})"
            elif stop_at_timestamp is not None and time.time() >= stop_at_timestamp:
                time_budget_exceeded = True
                state["stop_reason"] = f"time_budget_exceeded (stop_at_timestamp={stop_at_timestamp})"
            
            if time_budget_exceeded:
                current_state = self._transition(current_state, {"time_budget_reached": True})
                state["campaign_state"] = current_state
                self._persist_state(state)
                break
            
            # Check round limit
            completed_rounds = _as_int(state.get("round_index", -1), -1) + 1
            if max_rounds is not None and completed_rounds >= max_rounds:
                current_state = self._transition(current_state, {"round_limit_reached": True})
                state["campaign_state"] = current_state
                self._persist_state(state)
                break

            try:
                result = self.runner.run_adaptation_round()
            except Exception as exc:
                state["last_error"] = str(exc)
                current_state = self._transition(current_state, {"round_exception": True})
                state["campaign_state"] = current_state
                self._persist_state(state)
                break

            round_index = _as_int(state.get("round_index", -1), -1) + 1
            state["round_index"] = round_index
            state["last_round_result"] = result
            state["consecutive_empty_rounds"] = (
                _as_int(state.get("consecutive_empty_rounds", 0)) + 1
                if _as_int(result.get("events_emitted", 0)) == 0
                else 0
            )
            metrics = self._compute_metrics()
            self._append_supervisor_round_event(round_index, result, metrics, current_state)
            replay_results = self._replay_telemetry()
            state["telemetry_corruption"] = any(
                result_item.get("status") == "fatal_error" or _as_int(result_item.get("lines_errored", 0)) > 0
                for result_item in replay_results
            )

            snapshot = self._project_snapshot(round_index, current_state, state)
            snapshot_history = state.setdefault("snapshot_history", [])
            assert isinstance(snapshot_history, list)
            snapshot_history.append(snapshot)

            convergence = self.convergence_detector.evaluate(
                snapshot_history,
                {arm_id: _as_dict(stats).get("pulls", 0) for arm_id, stats in self.runner.bandit_state.arm_stats.items()},
            )
            state["convergence_state"] = convergence.get("state", "running")
            state["convergence_reason"] = convergence.get("reason", "")
            state["plateau_warning"] = convergence.get("state") == "plateau_warning"

            scaling = self.autoscaling_policy.evaluate(self._autoscaling_state(state, snapshot))
            state["last_scaling_action"] = scaling
            self._apply_autoscaling(scaling)

            alerts = self.alert_engine.evaluate(self._alert_state(state, snapshot))
            state["active_alerts"] = alerts.get("active_alerts", [])
            state["resolved_alerts"] = alerts.get("resolved_alerts", [])
            self.alert_engine.emit_alerts(alerts.get("new_alerts", []), self.outdir)
            self._persist_alerts_to_sqlite()

            current_state = self._transition(
                current_state,
                {
                    "convergence_state": state["convergence_state"],
                    "round_limit_reached": max_rounds is not None and (round_index + 1) >= max_rounds,
                },
            )
            state["campaign_state"] = current_state
            self._persist_state(state)

        return current_state

    def _ensure_manifest(self) -> dict[str, str]:
        if self.manifest_path.exists():
            return self.dependencies.manifest_validator(str(self.outdir), self.config)
        self.dependencies.manifest_writer(self.outdir, self.config)
        return {"status": "pass", "reason": "Manifest created"}

    def _initialize_telemetry(self, *, rebuild: bool) -> list[dict[str, Any]]:
        if rebuild and self.telemetry_db_path.exists():
            self.telemetry_db_path.unlink()
        self._ensure_supervisor_sources()
        return self._replay_telemetry()

    def _replay_telemetry(self) -> list[dict[str, Any]]:
        conn = self.dependencies.telemetry_init(self.telemetry_db_path)
        try:
            campaign_id = str(self.config.get("campaign_id", "UNKNOWN"))
            results = self.dependencies.telemetry_replay(conn, self._telemetry_source_paths(), campaign_id)
            results.extend(self._ingest_artifact_metadata_sources(conn, campaign_id))
            return results
        finally:
            conn.close()

    def _telemetry_source_paths(self) -> list[Path]:
        return [
            self.events_path,
            self.supervisor_events_path,
            self.task_summary_path,
            self.orchestrator_log_path,
        ]

    def _ingest_artifact_metadata_sources(
        self,
        conn: Any,
        campaign_id: str,
    ) -> list[dict[str, Any]]:
        if not self.task_summary_path.exists():
            return []

        results: list[dict[str, Any]] = []
        with self.task_summary_path.open("r", encoding="utf-8") as handle:
            for line_offset, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                with suppress(json.JSONDecodeError):
                    record = json.loads(line)
                    if not isinstance(record, dict):
                        continue
                    task_id_obj = record.get("task_id")
                    task_id = str(task_id_obj) if task_id_obj is not None else None
                    payload = record.get("payload")
                    payload_dict = payload if isinstance(payload, dict) else {}
                    for artifact_type, key in (("stdout", "stdout_path"), ("stderr", "stderr_path")):
                        path_obj = payload_dict.get(key, record.get(key))
                        if isinstance(path_obj, str) and path_obj:
                            results.append(
                                ingest_artifact_metadata(
                                    conn,
                                    artifact_type=artifact_type,
                                    artifact_path=path_obj,
                                    campaign_id=campaign_id,
                                    task_id=task_id,
                                    source_file=self.task_summary_path,
                                    line_offset=line_offset,
                                )
                            )
        return results

    def _project_snapshot(self, round_index: int, campaign_state: str, state: dict[str, object]) -> dict[str, object]:
        conn = self.dependencies.telemetry_init(self.telemetry_db_path)
        try:
            snapshot = project_round_snapshot(
                conn,
                round_index,
                campaign_state=campaign_state,
                active_alerts=state.get("active_alerts") if isinstance(state.get("active_alerts"), list) else None,
            )
        finally:
            conn.close()
        snapshot["last_scaling_action"] = state.get("last_scaling_action", {})
        snapshot["manifest_compatibility"] = state.get("manifest_status", "unknown")

        history = state.get("snapshot_history", [])
        last_success = -1
        if isinstance(history, list):
            for snap in reversed(history):
                if isinstance(snap, dict) and snap.get("subprocess_status") == "success":
                    last_success = _as_int(snap.get("round_index", -1), -1)
                    break
        if snapshot.get("subprocess_status") == "success" and round_index > last_success:
            last_success = round_index
        snapshot["last_successful_round"] = last_success if last_success >= 0 else None

        file_ts = {}
        run_dirs_path = self.outdir / "run_dirs"
        if run_dirs_path.exists():
            import datetime
            run_dirs = sorted([d for d in run_dirs_path.iterdir() if d.is_dir()], key=lambda p: p.stat().st_mtime, reverse=True)
            if run_dirs:
                latest_dir = run_dirs[0]
                for f in latest_dir.rglob("*"):
                    if f.is_file() and f.suffix in {".csv", ".json", ".jsonl"}:
                        mtime = datetime.datetime.fromtimestamp(f.stat().st_mtime, datetime.timezone.utc).isoformat().replace("+00:00", "Z")
                        file_ts[f"{latest_dir.name}/{f.name}"] = mtime
        snapshot["file_timestamps"] = file_ts

        write_status_json(self.outdir, snapshot)
        write_status_html(self.outdir, snapshot)
        return snapshot


    def _persist_alerts_to_sqlite(self) -> None:
        """Persist current alert state to SQLite telemetry database."""
        conn = self.dependencies.telemetry_init(self.telemetry_db_path)
        try:
            campaign_id = str(self.config.get("campaign_id", "unknown"))
            self.alert_engine.persist_alerts_to_sqlite(conn, campaign_id)
        finally:
            conn.close()
    def _compute_metrics(self) -> dict[str, float]:
        history = list(self.runner.bandit_state.global_history)
        coverage = compute_coverage(history, self.config)
        diversity = compute_diversity(history, self.config)
        yield_val = _last_yield(history)
        composite = compute_composite(yield_val, coverage, diversity, self.config)
        return {
            "yield": yield_val,
            "coverage": coverage,
            "diversity": diversity,
            "composite": composite,
        }

    def _append_supervisor_round_event(
        self,
        round_index: int,
        result: Mapping[str, object],
        metrics: Mapping[str, float],
        campaign_state: str,
    ) -> None:
        self._ensure_supervisor_sources()
        event = {
            "schema_version": 1,
            "campaign_id": self.config.get("campaign_id", ""),
            "event_type": "ROUND_COMPLETED",
            "utc": result.get("utc", self._utc_now()),
            "payload": {
                "round_index": round_index,
                "arm_id": result.get("arm_id", ""),
                "selected_arm": result.get("arm_id", ""),
                "discoveries": _as_int(result.get("discoveries", 0)),
                "events_emitted": _as_int(result.get("events_emitted", 0)),
                "subprocess_status": result.get("subprocess_status", "unknown"),
                "metrics": dict(metrics),
                "campaign_state": campaign_state,
                "active_alerts": [],
            },
        }
        with self.supervisor_events_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(event) + "\n")

    def _autoscaling_state(self, state: Mapping[str, object], snapshot: Mapping[str, object]) -> dict[str, object]:
        runtime = _as_dict(self.config.get("runtime"))
        limits = _as_dict(self.config.get("limits"))
        rolling_deltas = _as_dict(snapshot.get("rolling_deltas"))
        return {
            "round_index": _as_int(state.get("round_index", 0)),
            "threads": _as_int(runtime.get("threads", 1), 1),
            "timeout_sec": _as_int(runtime.get("timeout_sec", 1), 1),
            "max_new_run_dirs_per_round": _as_int(limits.get("max_new_run_dirs_per_round", 1), 1),
            "max_new_run_dirs_per_arm_call": _as_int(limits.get("max_new_run_dirs_per_arm_call", 1), 1),
            "empty_rounds": _as_int(state.get("consecutive_empty_rounds", 0)),
            "plateau_warning": state.get("plateau_warning", False),
            "coverage_gain": _as_float(rolling_deltas.get("coverage_delta", 0.0), 0.0),
            "diversity_gain": _as_float(rolling_deltas.get("diversity_delta", 0.0), 0.0),
            "host_pressure": "low",
            "last_scaling_round": _as_dict(state.get("last_scaling_action")).get("last_scaling_round"),
        }

    def _apply_autoscaling(self, scaling: Mapping[str, object]) -> None:
        if scaling.get("action") not in {"scale_up", "scale_down"}:
            return
        runtime = _as_dict(self.config.get("runtime"))
        limits = _as_dict(self.config.get("limits"))
        runtime["threads"] = _as_int(scaling.get("next_threads", runtime.get("threads", 1)), 1)
        runtime["timeout_sec"] = _as_int(scaling.get("next_timeout_sec", runtime.get("timeout_sec", 1)), 1)
        limits["max_new_run_dirs_per_round"] = _as_int(
            scaling.get("next_max_new_run_dirs_per_round", limits.get("max_new_run_dirs_per_round", 1)),
            1,
        )
        limits["max_new_run_dirs_per_arm_call"] = _as_int(
            scaling.get("next_max_new_run_dirs_per_arm_call", limits.get("max_new_run_dirs_per_arm_call", 1)),
            1,
        )
        self.config["runtime"] = runtime
        self.config["limits"] = limits

    def _alert_state(self, state: Mapping[str, object], snapshot: Mapping[str, object]) -> dict[str, object]:
        return {
            "campaign_state": state.get("campaign_state", "RUNNING"),
            "preflight_overall": state.get("preflight_overall", "unknown"),
            "preflight_checks": state.get("preflight_checks", []),
            "subprocess_status": snapshot.get("subprocess_status", "unknown"),
            "round_index": state.get("round_index", -1),
            "consecutive_empty_rounds": state.get("consecutive_empty_rounds", 0),
            "max_empty_rounds": _as_int(_as_dict(self.config.get("runtime")).get("max_empty_rounds", 3), 3),
            "plateau_warning": state.get("plateau_warning", False),
            "convergence_state": state.get("convergence_state", "running"),
            "telemetry_corruption": state.get("telemetry_corruption", False),
            "quarantine_status": state.get("quarantine_status", snapshot.get("quarantine_status", "clear")),
            "quarantine_summary": state.get("quarantine_summary", snapshot.get("quarantine", {})),
            "active_alerts": state.get("active_alerts", []),
            "rolling_metrics": snapshot.get("rolling_metrics", {}),
        }

    def _max_rounds(self) -> int | None:
        supervisor = _as_dict(self.config.get("supervisor"))
        runtime = _as_dict(self.config.get("runtime"))
        raw = supervisor.get("max_rounds", runtime.get("max_rounds"))
        if raw is None:
            return None
        value = _as_int(raw, 0)
        return value if value > 0 else None

    def _load_state(self) -> dict[str, object]:
        if not self.state_path.exists():
            return {
                "campaign_id": self.config.get("campaign_id", ""),
                "campaign_state": "INIT",
                "round_index": -1,
                "consecutive_empty_rounds": 0,
                "snapshot_history": [],
                "active_alerts": [],
            }
        return json.loads(self.state_path.read_text(encoding="utf-8"))

    def _persist_state(self, state: Mapping[str, object]) -> None:
        payload = dict(state)
        payload["campaign_id"] = self.config.get("campaign_id", "")
        self.state_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    def _ensure_supervisor_sources(self) -> None:
        self.outdir.mkdir(parents=True, exist_ok=True)
        if not self.events_path.exists():
            self.events_path.touch()
        if not self.supervisor_events_path.exists():
            self.supervisor_events_path.touch()

    @staticmethod
    def _utc_now() -> str:
        import datetime as _datetime

        return _datetime.datetime.now(_datetime.timezone.utc).isoformat().replace("+00:00", "Z")


__all__ = ["CampaignSupervisor", "TERMINAL_STATES"]
