from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import subprocess
import sys
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from autoresearch.harness.dihiggs_axis_contract import DIHIGGS_EXPLORERS_MODE, encode_cell_id, substitute_placeholders
from autoresearch.harness.dihiggs_parsers import (
    AdaptiveProposalInventory,
    parse_adaptive_checkpoint_inventory,
    parse_branch_checkpoint_inventory,
)
from autoresearch.harness.dihiggs_preflight import run_physscanwithfixings_preflight
from autoresearch.harness.evaluator import evaluate_run_dir
from autoresearch.harness.models import ArtifactRef, EvalResult
from autoresearch.harness.scoring import fingerprint_artifacts
from autoresearch.harness.state import append_event, replay_events
from autoresearch.benchmarks.score import compute_composite, compute_metrics


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class DiscoveredRun:
    arm_id: str
    explorer: str
    run_dir: str | None
    discovered_at_utc: str
    source_ref: str


@dataclass(frozen=True)
class _SliceAttempt:
    discovered_run: DiscoveredRun
    iter_index: int
    attempt_index: int
    tb: int
    lam1_bin: int | None
    lam1_raw: float | None
    successes: int | None
    trials: int | None
    slice_key: str


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _sha256_text(*parts: object) -> str:
    joined = "".join(str(part) for part in parts)
    return hashlib.sha256(joined.encode("utf-8")).hexdigest()


def _coerce_int(value: object, *, default: int = 0) -> int:
    if isinstance(value, bool):
        return default
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    if isinstance(value, str):
        try:
            return int(value.strip())
        except ValueError:
            return default
    return default


def _coerce_float(value: object) -> float | None:
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            return None
    return None


def _coerce_bool(value: object, *, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"1", "true", "yes", "y", "on"}:
            return True
        if lowered in {"0", "false", "no", "n", "off"}:
            return False
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return bool(value)
    return default


def _failure_result(status: str, *, elapsed_sec: float = 0.0, details: str | None = None) -> EvalResult:
    artifact_ref = ArtifactRef(run_dir="", manifest_path=None, jsonl_path=None, csv_paths=[])
    return EvalResult(
        status=status,
        successes=0,
        trials=0,
        elapsed_sec=float(elapsed_sec),
        event_counts={"done": 0, "skip": 0, "dry_run": 0, "fail": 0, "crash": 0},
        reuse_detected=False,
        artifact_ref=artifact_ref,
        details=details,
    )


class DiHiggsRunner:
    def __init__(self, config: dict[str, object], outdir: str):
        self.config: dict[str, object] = config
        self.outdir: str = outdir
        self.known_attempt_ids: set[str] = set()
        self.known_run_dirs: set[str] = set()
        self.seen_fingerprints: set[str] = set()
        self._round_counter: dict[str, int] = {}
        self.state_file_path: Path = self._state_file_path()
        self._bandit_state: dict[str, object] = {"total_rounds": 0, "arms": {}}

        paths_cfg = self._mapping(self.config.get("paths"), "paths")
        self.event_log_path: Path = Path(outdir) / str(paths_cfg.get("lake_name", "events.jsonl"))

        self._ensure_outdir_layout()
        self._load_existing_state()
        self._load_mode_state()
        self._ensure_arm_state_entries()
        self._save_mode_state()

    @staticmethod
    def _mapping(value: object, name: str) -> Mapping[str, object]:
        if value is None:
            return {}
        if not isinstance(value, Mapping):
            raise ValueError(f"config.{name} must be a mapping")
        return value

    @staticmethod
    def make_attempt_id(campaign_id: str, arm_id: str, fingerprint: str, slice_key: str) -> str:
        return _sha256_text(campaign_id, arm_id, fingerprint, slice_key)

    def _ensure_outdir_layout(self) -> None:
        Path(self.outdir).mkdir(parents=True, exist_ok=True)
        self.event_log_path.parent.mkdir(parents=True, exist_ok=True)
        self.state_file_path.parent.mkdir(parents=True, exist_ok=True)

    def _state_file_path(self) -> Path:
        campaign_id = str(self.config.get("campaign_id", "")).strip()
        if not campaign_id:
            campaign_id = "unknown-campaign"
        return Path("autoresearch") / "opencode_logs" / campaign_id / "state" / "dihiggs_mode_state.json"

    def _arm_ids(self) -> list[str]:
        arms = self.config.get("arms")
        if not isinstance(arms, list):
            raise ValueError("config.arms must be a list")
        arm_ids: list[str] = []
        for arm in arms:
            if not isinstance(arm, Mapping):
                raise ValueError("Each config.arms entry must be a mapping")
            arm_id = arm.get("id")
            if not isinstance(arm_id, str) or not arm_id.strip():
                raise ValueError("Each config.arms entry must define a non-empty id")
            arm_ids.append(arm_id.strip())
        if not arm_ids:
            raise ValueError("config.arms must contain at least one arm")
        return arm_ids

    def _adaptation_cfg(self) -> Mapping[str, object]:
        return self._mapping(self.config.get("adaptation"), "adaptation")

    def _load_mode_state(self) -> None:
        if not self.state_file_path.exists():
            return
        try:
            raw = json.loads(self.state_file_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            logger.warning("Failed to read dihiggs mode state from %s; reinitializing", self.state_file_path)
            return
        if not isinstance(raw, Mapping):
            return
        raw_arms = raw.get("arms")
        arms_state: dict[str, dict[str, float | int]] = {}
        if isinstance(raw_arms, Mapping):
            for arm_id, arm_state in raw_arms.items():
                if not isinstance(arm_id, str) or not isinstance(arm_state, Mapping):
                    continue
                arms_state[arm_id] = {
                    "n": max(0, _coerce_int(arm_state.get("n"), default=0)),
                    "mean_reward": max(0.0, min(1.0, float(_coerce_float(arm_state.get("mean_reward")) or 0.0))),
                }
        self._bandit_state = {
            "total_rounds": max(0, _coerce_int(raw.get("total_rounds"), default=0)),
            "arms": arms_state,
        }

    def _ensure_arm_state_entries(self) -> None:
        arms_state_obj = self._bandit_state.setdefault("arms", {})
        if not isinstance(arms_state_obj, dict):
            arms_state_obj = {}
            self._bandit_state["arms"] = arms_state_obj
        for arm_id in self._arm_ids():
            if arm_id not in arms_state_obj or not isinstance(arms_state_obj[arm_id], dict):
                arms_state_obj[arm_id] = {"n": 0, "mean_reward": 0.0}
                continue
            arm_state = arms_state_obj[arm_id]
            arm_state["n"] = max(0, _coerce_int(arm_state.get("n"), default=0))
            arm_state["mean_reward"] = max(0.0, min(1.0, float(_coerce_float(arm_state.get("mean_reward")) or 0.0)))

    def _save_mode_state(self) -> None:
        self._ensure_arm_state_entries()
        arms_state = self._bandit_state.get("arms")
        if not isinstance(arms_state, dict):
            arms_state = {}
        serialized_arms: dict[str, dict[str, float | int]] = {}
        for arm_id in self._arm_ids():
            arm_state = arms_state.get(arm_id)
            if not isinstance(arm_state, Mapping):
                arm_state = {}
            serialized_arms[arm_id] = {
                "n": max(0, _coerce_int(arm_state.get("n"), default=0)),
                "mean_reward": max(0.0, min(1.0, float(_coerce_float(arm_state.get("mean_reward")) or 0.0))),
            }
        state_payload = {
            "schema_version": 1,
            "campaign_id": str(self.config.get("campaign_id", "")),
            "total_rounds": max(0, _coerce_int(self._bandit_state.get("total_rounds"), default=0)),
            "arms": serialized_arms,
        }
        self.state_file_path.write_text(json.dumps(state_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    def _compute_campaign_score(self) -> float:
        metrics = compute_metrics(str(self.event_log_path), self.config)
        composite = compute_composite(metrics)
        if isinstance(composite, Mapping):
            score = _coerce_float(composite.get("score")) or 0.0
        else:
            score = _coerce_float(composite) or 0.0
        return max(0.0, min(1.0, float(score)))

    def _select_arm_ucb1(self) -> str:
        self._ensure_arm_state_entries()
        adaptation_cfg = self._adaptation_cfg()
        warm_start_each_arm = _coerce_bool(adaptation_cfg.get("warm_start_each_arm"), default=True)
        arms_state = self._bandit_state.get("arms")
        if not isinstance(arms_state, Mapping):
            raise ValueError("bandit arm state must be a mapping")
        arm_ids = self._arm_ids()
        if warm_start_each_arm:
            for arm_id in arm_ids:
                arm_state = arms_state.get(arm_id)
                pulls = _coerce_int(arm_state.get("n"), default=0) if isinstance(arm_state, Mapping) else 0
                if pulls == 0:
                    return arm_id

        total_rounds = max(1, _coerce_int(self._bandit_state.get("total_rounds"), default=0))
        ucb_c = float(_coerce_float(adaptation_cfg.get("ucb_c")) or math.sqrt(2.0))
        best_arm = arm_ids[0]
        best_score = -math.inf
        for arm_id in arm_ids:
            arm_state = arms_state.get(arm_id)
            pulls = _coerce_int(arm_state.get("n"), default=0) if isinstance(arm_state, Mapping) else 0
            if pulls <= 0:
                return arm_id
            mean_reward = float(_coerce_float(arm_state.get("mean_reward")) or 0.0) if isinstance(arm_state, Mapping) else 0.0
            ucb_score = mean_reward + ucb_c * math.sqrt(math.log(total_rounds) / pulls)
            if ucb_score > best_score:
                best_score = ucb_score
                best_arm = arm_id
        return best_arm

    def _record_arm_reward(self, arm_id: str, reward: float) -> None:
        self._ensure_arm_state_entries()
        arms_state = self._bandit_state.get("arms")
        if not isinstance(arms_state, dict):
            raise ValueError("bandit arm state must be a dict")
        arm_state = arms_state.get(arm_id)
        if not isinstance(arm_state, dict):
            arm_state = {"n": 0, "mean_reward": 0.0}
            arms_state[arm_id] = arm_state
        pulls_before = max(0, _coerce_int(arm_state.get("n"), default=0))
        mean_before = max(0.0, min(1.0, float(_coerce_float(arm_state.get("mean_reward")) or 0.0)))
        pulls_after = pulls_before + 1
        bounded_reward = max(0.0, min(1.0, float(reward)))
        mean_after = ((mean_before * pulls_before) + bounded_reward) / pulls_after
        arm_state["n"] = pulls_after
        arm_state["mean_reward"] = mean_after
        self._bandit_state["total_rounds"] = max(0, _coerce_int(self._bandit_state.get("total_rounds"), default=0)) + 1
        self._save_mode_state()

    def _load_existing_state(self) -> None:
        if not self.event_log_path.exists():
            return
        for event in replay_events(self.event_log_path):
            payload = event.get("payload")
            if not isinstance(payload, Mapping):
                continue
            attempt_id = payload.get("attempt_id")
            if isinstance(attempt_id, str) and attempt_id:
                self.known_attempt_ids.add(attempt_id)
            run_dir = payload.get("run_dir")
            if isinstance(run_dir, str) and run_dir:
                self.known_run_dirs.add(run_dir)
            fingerprint = payload.get("fingerprint_sha256")
            if isinstance(fingerprint, str) and fingerprint:
                self.seen_fingerprints.add(fingerprint)

    def _event(self, event_type: str, payload: Mapping[str, object]) -> dict[str, object]:
        return {
            "schema_version": 1,
            "campaign_id": str(self.config.get("campaign_id", "")),
            "event_type": event_type,
            "utc": _utc_now_iso(),
            "payload": dict(payload),
        }

    def _append_event(self, event_type: str, payload: Mapping[str, object]) -> None:
        append_event(self.event_log_path, self._event(event_type, payload))

    def _first_arm_id(self) -> str:
        arms = self.config.get("arms")
        if not isinstance(arms, list) or not arms:
            raise ValueError("config.arms must contain at least one arm")
        arm = arms[0]
        if not isinstance(arm, Mapping):
            raise ValueError("config.arms[0] must be a mapping")
        arm_id = arm.get("id")
        if not isinstance(arm_id, str) or not arm_id.strip():
            raise ValueError("config.arms[0].id must be a non-empty string")
        return arm_id.strip()

    def _resolve_arm(self, arm_id: str) -> tuple[dict[str, object], dict[str, object]]:
        arms = self.config.get("arms")
        if not isinstance(arms, list):
            raise ValueError("config.arms must be a list")
        arm_cfg: dict[str, object] | None = None
        for arm in arms:
            if isinstance(arm, dict) and arm.get("id") == arm_id:
                arm_cfg = arm
                break
        if arm_cfg is None:
            raise ValueError(f"Unknown arm id: {arm_id}")

        paths = self._mapping(self.config.get("paths"), "paths")
        runtime = self._mapping(self.config.get("runtime"), "runtime")
        dihiggs = self._mapping(self.config.get("dihiggs"), "dihiggs")
        search = self._mapping(self.config.get("search"), "search")
        lam1_cfg_obj = search.get("lam1")
        lam1_cfg = lam1_cfg_obj if isinstance(lam1_cfg_obj, Mapping) else {}
        tb_values_obj = search.get("tb_values")
        tb_values = tb_values_obj if isinstance(tb_values_obj, list) else []
        first_tb = tb_values[0] if tb_values else 0

        runtime_ctx: dict[str, object] = {
            "sys_executable": sys.executable,
            "python": sys.executable,
            "repo_root": paths.get("repo_root", ""),
            "outdir": self.outdir,
            "checkpoint_root": str(Path(self.outdir) / "checkpoints"),
            "campaign": self.config.get("campaign_id", ""),
            "threads": runtime.get("threads", 1),
            "tb_values": ",".join(str(v) for v in tb_values),
            "tb": first_tb,
            "track_id": f"{arm_id}-track-0001",
            "iter": self._round_counter.get(arm_id, 0),
            "lam1_min": lam1_cfg.get("min", 0.0) if isinstance(lam1_cfg, Mapping) else 0.0,
            "lam1_max": lam1_cfg.get("max", 0.0) if isinstance(lam1_cfg, Mapping) else 0.0,
            "lam1_n_bins": lam1_cfg.get("n_bins", 1) if isinstance(lam1_cfg, Mapping) else 1,
            "hb_dataset": os.environ.get(str(dihiggs.get("hb_dataset_env", "HB_DATASET")), ""),
            "hs_dataset": os.environ.get(str(dihiggs.get("hs_dataset_env", "HS_DATASET")), ""),
        }

        resolved = substitute_placeholders(self.config, runtime_ctx)
        resolved_arms = resolved.get("arms")
        if not isinstance(resolved_arms, list):
            raise ValueError("resolved config.arms must be a list")
        for arm in resolved_arms:
            if isinstance(arm, dict) and arm.get("id") == arm_id:
                return arm, runtime_ctx
        raise ValueError(f"resolved arm config missing for {arm_id}")

    def _run_subprocess(self, arm_id: str, arm_cfg: Mapping[str, object]) -> str:
        cmd_raw = arm_cfg.get("cmd", [])
        if not isinstance(cmd_raw, list):
            raise ValueError(f"arm[{arm_id}].cmd must be a list")
        cmd = [str(item) for item in cmd_raw]

        env_raw = arm_cfg.get("env", {})
        if not isinstance(env_raw, Mapping):
            raise ValueError(f"arm[{arm_id}].env must be a mapping")
        env = {**os.environ, **{str(key): str(value) for key, value in env_raw.items()}}

        runtime_cfg = self._mapping(self.config.get("runtime"), "runtime")
        paths_cfg = self._mapping(self.config.get("paths"), "paths")
        timeout_sec = _coerce_int(arm_cfg.get("timeout_sec"), default=_coerce_int(runtime_cfg.get("timeout_sec"), default=3600))
        cwd = str(paths_cfg.get("repo_root", os.getcwd()))

        round_num = self._round_counter.get(arm_id, 0)
        stdout_path = Path(self.outdir) / f"{arm_id}_round_{round_num}.stdout.txt"
        stderr_path = Path(self.outdir) / f"{arm_id}_round_{round_num}.stderr.txt"

        try:
            result = subprocess.run(
                cmd,
                env=env,
                capture_output=True,
                timeout=timeout_sec,
                cwd=cwd,
                check=False,
            )
            _ = stdout_path.write_text(result.stdout.decode("utf-8", errors="replace"), encoding="utf-8")
            _ = stderr_path.write_text(result.stderr.decode("utf-8", errors="replace"), encoding="utf-8")
            if result.returncode != 0:
                return "nonzero_exit"
            return "success"
        except subprocess.TimeoutExpired as exc:
            stdout = exc.stdout if exc.stdout is not None else b""
            stderr = exc.stderr if exc.stderr is not None else b""
            stdout_text = stdout if isinstance(stdout, str) else stdout.decode("utf-8", errors="replace")
            stderr_text = stderr if isinstance(stderr, str) else stderr.decode("utf-8", errors="replace")
            _ = stdout_path.write_text(stdout_text, encoding="utf-8")
            _ = stderr_path.write_text(stderr_text, encoding="utf-8")
            return "timeout"

    def _checkpoint_root_for_arm(self, runtime_ctx: Mapping[str, object], arm_id: str) -> str:
        return str(Path(str(runtime_ctx["checkpoint_root"])) / arm_id)

    def _iter_slice_keys(self, proposal: AdaptiveProposalInventory) -> Iterable[tuple[int, int, str]]:
        ordered_keys: list[str] = []
        for key in sorted(proposal.successes_by_tb):
            if key not in ordered_keys:
                ordered_keys.append(key)
        for key in sorted(proposal.trials_by_tb):
            if key not in ordered_keys:
                ordered_keys.append(key)
        if not ordered_keys:
            fallback_tb = proposal.bin_index if proposal.bin_index is not None else 0
            ordered_keys.append(str(fallback_tb))
        for key in ordered_keys:
            yield (_coerce_int(key, default=0), proposal.iter_index, key)

    def _adaptive_attempts(self, arm_id: str, checkpoint_root: str) -> list[_SliceAttempt]:
        attempts: list[_SliceAttempt] = []
        for proposal in parse_adaptive_checkpoint_inventory(checkpoint_root):
            if proposal.status != "DONE":
                continue
            discovered_run = DiscoveredRun(
                arm_id=arm_id,
                explorer="adaptive",
                run_dir=proposal.run_dir,
                discovered_at_utc=_utc_now_iso(),
                source_ref=f"{checkpoint_root}/iter_{proposal.iter_index:04d}/adaptive_state.json#{proposal.proposal_id}",
            )
            for tb_value, iter_index, key in self._iter_slice_keys(proposal):
                attempts.append(
                    _SliceAttempt(
                        discovered_run=discovered_run,
                        iter_index=iter_index,
                        attempt_index=proposal.proposal_index,
                        tb=tb_value,
                        lam1_bin=None,
                        lam1_raw=proposal.lam1_min,
                        successes=proposal.successes_by_tb.get(key),
                        trials=proposal.trials_by_tb.get(key),
                        slice_key=f"tb={tb_value}",
                    )
                )
        return attempts

    def _branch_attempts(self, arm_id: str, checkpoint_root: str) -> list[_SliceAttempt]:
        attempts: list[_SliceAttempt] = []
        for event in parse_branch_checkpoint_inventory(checkpoint_root):
            if event.status and event.status != "DONE":
                continue
            tb_value = _coerce_int(event.tanbeta, default=0)
            discovered_run = DiscoveredRun(
                arm_id=arm_id,
                explorer="branch",
                run_dir=event.run_dir,
                discovered_at_utc=_utc_now_iso(),
                source_ref=f"{checkpoint_root}/{event.track_id}/events.jsonl#{event.line_index}",
            )
            attempts.append(
                _SliceAttempt(
                    discovered_run=discovered_run,
                    iter_index=0,
                    attempt_index=event.line_index,
                    tb=tb_value,
                    lam1_bin=event.lam1_bin,
                    lam1_raw=event.lam1_raw,
                    successes=None,
                    trials=None,
                    slice_key=f"tb={tb_value}",
                )
            )
        return attempts

    def _slice_attempts(self, arm_id: str, explorer: str, checkpoint_root: str) -> list[_SliceAttempt]:
        if explorer == "adaptive":
            attempts = self._adaptive_attempts(arm_id, checkpoint_root)
        elif explorer == "branch":
            attempts = self._branch_attempts(arm_id, checkpoint_root)
        else:
            raise ValueError(f"Unsupported explorer type: {explorer}")

        search_cfg = self._mapping(self.config.get("search"), "search")
        lam1_cfg = search_cfg.get("lam1")
        if not isinstance(lam1_cfg, Mapping):
            lam1_cfg = {}
        lam1_min_cfg = _coerce_float(lam1_cfg.get("min"))
        lam1_max_cfg = _coerce_float(lam1_cfg.get("max"))
        n_bins = max(1, _coerce_int(lam1_cfg.get("n_bins"), default=1))

        resolved: list[_SliceAttempt] = []
        for attempt in attempts:
            lam1_bin = attempt.lam1_bin
            if lam1_bin is None and attempt.lam1_raw is not None and lam1_min_cfg is not None and lam1_max_cfg is not None:
                span = lam1_max_cfg - lam1_min_cfg
                if span <= 0:
                    lam1_bin = 0
                elif attempt.lam1_raw <= lam1_min_cfg:
                    lam1_bin = 0
                elif attempt.lam1_raw >= lam1_max_cfg:
                    lam1_bin = n_bins - 1
                else:
                    rel = (attempt.lam1_raw - lam1_min_cfg) / span
                    lam1_bin = min(n_bins - 1, max(0, int(rel * n_bins)))
            resolved.append(
                _SliceAttempt(
                    discovered_run=attempt.discovered_run,
                    iter_index=attempt.iter_index,
                    attempt_index=attempt.attempt_index,
                    tb=attempt.tb,
                    lam1_bin=lam1_bin,
                    lam1_raw=attempt.lam1_raw,
                    successes=attempt.successes,
                    trials=attempt.trials,
                    slice_key=attempt.slice_key,
                )
            )
        return resolved

    def _evaluate_discovered_run(self, discovered_run: DiscoveredRun) -> tuple[EvalResult, ArtifactRef | None, str | None, bool]:
        run_dir = discovered_run.run_dir
        if not run_dir:
            failure = _failure_result(
                "CRASH",
                details=f"checkpoint entry missing run_dir: {discovered_run.source_ref}",
            )
            return failure, None, None, False

        eval_result = evaluate_run_dir(run_dir, timed_out=False)
        artifact_ref = eval_result.artifact_ref
        if artifact_ref is None:
            artifact_ref = ArtifactRef(run_dir=run_dir, manifest_path=None, jsonl_path=None, csv_paths=[])
            eval_result.artifact_ref = artifact_ref
        fingerprint: str | None = None
        is_duplicate = False
        if eval_result.status != "BROKEN_ARTIFACT":
            fingerprint, _ = fingerprint_artifacts(artifact_ref)
            artifact_ref.fingerprint_sha256 = fingerprint
            is_duplicate = fingerprint in self.seen_fingerprints
            self.seen_fingerprints.add(fingerprint)
        return eval_result, artifact_ref, fingerprint, is_duplicate

    def _terminal_event_type(self, status: str) -> str:
        if status in {"TIMEOUT", "CRASH", "FAIL", "BROKEN_ARTIFACT"}:
            return "ATTEMPT_FAILED"
        return "ATTEMPT_SUCCEEDED"

    def _emit_attempt_events(
        self,
        *,
        arm_id: str,
        explorer: str,
        slice_attempt: _SliceAttempt,
        eval_result: EvalResult,
        artifact_ref: ArtifactRef | None,
        fingerprint: str | None,
        is_duplicate: bool,
    ) -> bool:
        lam1_bin = slice_attempt.lam1_bin if slice_attempt.lam1_bin is not None else 0
        axes_binned = {"tb": int(slice_attempt.tb), "lam1_bin": int(lam1_bin)}
        axes_raw: dict[str, object] = {"tb": int(slice_attempt.tb)}
        if slice_attempt.lam1_raw is not None:
            axes_raw["lam1"] = float(slice_attempt.lam1_raw)
        cell_id = encode_cell_id(axes_binned)
        stable_fingerprint = fingerprint or _sha256_text(slice_attempt.discovered_run.source_ref, slice_attempt.slice_key)
        attempt_id = self.make_attempt_id(str(self.config.get("campaign_id", "")), arm_id, stable_fingerprint, slice_attempt.slice_key)
        if attempt_id in self.known_attempt_ids:
            return False

        payload_base: dict[str, object] = {
            "attempt_id": attempt_id,
            "arm_id": arm_id,
            "explorer": explorer,
            "mode": DIHIGGS_EXPLORERS_MODE,
            "iter_index": int(slice_attempt.iter_index),
            "attempt_index": int(slice_attempt.attempt_index),
            "cell_id": cell_id,
            "slice_key": slice_attempt.slice_key,
            "axes_binned": axes_binned,
            "axes_raw": axes_raw,
            "source_ref": slice_attempt.discovered_run.source_ref,
            "discovered_at_utc": slice_attempt.discovered_run.discovered_at_utc,
            "is_duplicate": bool(is_duplicate),
        }
        if slice_attempt.discovered_run.run_dir:
            payload_base["run_dir"] = slice_attempt.discovered_run.run_dir
        if fingerprint:
            payload_base["fingerprint_sha256"] = fingerprint
        if artifact_ref is not None and artifact_ref.run_dir:
            payload_base["run_dir"] = artifact_ref.run_dir

        self._append_event("ATTEMPT_CREATED", payload_base)
        self._append_event("ATTEMPT_STARTED", payload_base)

        emitted_status = eval_result.status
        if is_duplicate:
            emitted_status = "SKIP_REUSED"

        successes = 0 if is_duplicate else int(eval_result.successes)
        if slice_attempt.successes is not None:
            successes = int(slice_attempt.successes)
        trials = 0 if is_duplicate else int(eval_result.trials)
        if slice_attempt.trials is not None:
            trials = int(slice_attempt.trials)
        elif slice_attempt.successes is not None and slice_attempt.successes > 0:
            trials = max(trials, int(slice_attempt.successes))
        if is_duplicate:
            successes = 0
            trials = 0

        evaluated_payload = {
            **payload_base,
            "eval_status": emitted_status,
            "successes": successes,
            "trials": trials,
            "elapsed_sec": float(eval_result.elapsed_sec),
            "reuse_detected": bool(eval_result.reuse_detected),
        }
        if eval_result.details:
            evaluated_payload["details"] = eval_result.details
        self._append_event("ATTEMPT_EVALUATED", evaluated_payload)

        terminal_payload = {
            **payload_base,
            "final_status": eval_result.status,
        }
        if eval_result.details:
            terminal_payload["details"] = eval_result.details
        self._append_event(self._terminal_event_type(eval_result.status), terminal_payload)

        self.known_attempt_ids.add(attempt_id)
        if slice_attempt.discovered_run.run_dir:
            self.known_run_dirs.add(slice_attempt.discovered_run.run_dir)
        return True

    def run_single_round(
        self,
        arm_id: str | None = None,
        update_bandit: bool = False,
        *,
        preflight_override: Mapping[str, object] | None = None,
    ) -> dict[str, object]:
        del update_bandit
        if arm_id is None:
            arm_id = self._first_arm_id()

        preflight = (
            dict(preflight_override)
            if isinstance(preflight_override, Mapping)
            else run_physscanwithfixings_preflight(self.config)
        )
        if preflight.get("overall") == "fail":
            return {
                "arm_id": arm_id,
                "discoveries": 0,
                "events_emitted": 0,
                "preflight": preflight,
                "subprocess_status": "skipped_preflight_fail",
            }

        arm_cfg, runtime_ctx = self._resolve_arm(arm_id)
        subprocess_status = self._run_subprocess(arm_id, arm_cfg)

        explorer = str(arm_cfg.get("explorer", ""))
        checkpoint_root = self._checkpoint_root_for_arm(runtime_ctx, arm_id)
        attempts = self._slice_attempts(arm_id, explorer, checkpoint_root)

        limits = self.config.get("limits")
        max_new = len(attempts)
        if isinstance(limits, Mapping):
            max_new = _coerce_int(limits.get("max_new_run_dirs_per_arm_call"), default=len(attempts))
            round_cap = _coerce_int(limits.get("max_new_run_dirs_per_round"), default=max_new)
            if round_cap >= 0:
                max_new = min(max_new, round_cap)
        if max_new >= 0:
            attempts = attempts[:max_new]

        events_emitted = 0
        seen_run_sources: set[str] = set()
        run_eval_cache: dict[str, tuple[EvalResult, ArtifactRef | None, str | None, bool]] = {}

        for slice_attempt in attempts:
            source_key = slice_attempt.discovered_run.source_ref
            if source_key not in run_eval_cache:
                run_eval_cache[source_key] = self._evaluate_discovered_run(slice_attempt.discovered_run)
            eval_result, artifact_ref, fingerprint, is_duplicate = run_eval_cache[source_key]
            if self._emit_attempt_events(
                arm_id=arm_id,
                explorer=explorer,
                slice_attempt=slice_attempt,
                eval_result=eval_result,
                artifact_ref=artifact_ref,
                fingerprint=fingerprint,
                is_duplicate=is_duplicate,
            ):
                events_emitted += 4
                seen_run_sources.add(source_key)

        self._round_counter[arm_id] = self._round_counter.get(arm_id, 0) + 1
        return {
            "arm_id": arm_id,
            "discoveries": len(seen_run_sources),
            "events_emitted": events_emitted,
            "preflight": preflight,
            "subprocess_status": subprocess_status,
        }

    def run_adaptation_round(self) -> dict[str, object]:
        arm_id = self._select_arm_ucb1()
        score_before = self._compute_campaign_score()
        round_summary = self.run_single_round(arm_id)
        score_after = self._compute_campaign_score()
        delta = max(-1.0, min(1.0, score_after - score_before))
        reward = max(0.0, min(1.0, (delta + 1.0) / 2.0))
        self._record_arm_reward(arm_id, reward)
        return {
            **round_summary,
            "score_before": score_before,
            "score_after": score_after,
            "score_delta": delta,
            "reward": reward,
            "bandit_state_path": str(self.state_file_path),
        }

    def run(
        self,
        max_rounds: int | None = None,
        *,
        preflight_override: Mapping[str, object] | None = None,
    ) -> dict[str, object]:
        runtime_cfg = self._mapping(self.config.get("runtime"), "runtime")
        total_rounds = max_rounds if max_rounds is not None else _coerce_int(runtime_cfg.get("max_rounds"), default=1)
        total_rounds = max(0, total_rounds)
        rounds: list[dict[str, object]] = []
        for _ in range(total_rounds):
            arm_id = self._select_arm_ucb1()
            score_before = self._compute_campaign_score()
            round_summary = self.run_single_round(arm_id, preflight_override=preflight_override)
            score_after = self._compute_campaign_score()
            delta = max(-1.0, min(1.0, score_after - score_before))
            reward = max(0.0, min(1.0, (delta + 1.0) / 2.0))
            self._record_arm_reward(arm_id, reward)
            rounds.append(
                {
                    **round_summary,
                    "score_before": score_before,
                    "score_after": score_after,
                    "score_delta": delta,
                    "reward": reward,
                }
            )
        return {
            "rounds": rounds,
            "total_rounds": total_rounds,
            "bandit_state_path": str(self.state_file_path),
        }
