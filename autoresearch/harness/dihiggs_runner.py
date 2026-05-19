from __future__ import annotations
# pyright: reportImplicitRelativeImport=false, reportUnannotatedClassAttribute=false, reportUnknownVariableType=false, reportUnknownArgumentType=false, reportUnknownMemberType=false, reportAttributeAccessIssue=false, reportArgumentType=false, reportUnusedCallResult=false, reportAny=false, reportUnnecessaryIsInstance=false

import hashlib
import json
import logging
import os
import subprocess
import sys
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path

from autoresearch.benchmarks.score import (
    compute_composite,
    compute_coverage,
    compute_diversity,
    compute_metrics,
)
from autoresearch.harness.dihiggs_adaptation import (
    BanditState,
    compute_arm_reward,
    select_arm_ucb1,
    update_state,
)
from autoresearch.harness.dihiggs_axis_contract import (
    encode_cell_id,
    substitute_placeholders,
)
from autoresearch.harness.dihiggs_parsers import (
    parse_adaptive_checkpoint,
    parse_branch_checkpoint,
)
from autoresearch.harness.dihiggs_preflight import run_all_preflight_checks


logger = logging.getLogger(__name__)


class DiHiggsRunner:
    def __init__(self, config: dict[str, object], outdir: str):
        self.config = config
        self.outdir = outdir
        self.known_run_dirs: set[str] = set()
        self.known_attempt_ids: set[str] = set()
        paths_cfg = self._mapping(self.config.get("paths"), "paths")
        self.event_log_path = Path(outdir) / str(paths_cfg.get("lake_name", "events.jsonl"))
        self._round_counter: dict[str, int] = {}
        self.bandit_state = BanditState(arm_stats={}, global_history=[], arm_histories={})

        self._ensure_outdir_layout()
        self._load_existing_state()

    def _ensure_outdir_layout(self) -> None:
        Path(self.outdir).mkdir(parents=True, exist_ok=True)
        self.event_log_path.parent.mkdir(parents=True, exist_ok=True)
        if not self.event_log_path.exists():
            self.event_log_path.touch()

    def _load_existing_state(self) -> None:
        """
        Reconstruct known_attempt_ids and BanditState from existing event log.
        
        Reads the event log (if it exists) and:
        1. Populates known_attempt_ids to prevent duplicate event emission
        2. Reconstructs BanditState by replaying ATTEMPT_EVALUATED events
        
        This enables safe restart without duplicating work or corrupting state.
        """
        if not self.event_log_path.exists():
            logger.debug(f"No existing event log found at {self.event_log_path}")
            return
        
        logger.info(f"Loading existing state from {self.event_log_path}")
        events_loaded = 0
        
        try:
            with self.event_log_path.open("r", encoding="utf-8") as f:
                for line_num, line in enumerate(f, start=1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    try:
                        event = json.loads(line)
                    except json.JSONDecodeError as exc:
                        logger.warning(f"Failed to parse event at line {line_num}: {exc}")
                        continue
                    
                    if not isinstance(event, dict):
                        logger.warning(f"Event at line {line_num} is not a JSON object")
                        continue
                    
                    # Only process ATTEMPT_EVALUATED events
                    if event.get("event_type") != "ATTEMPT_EVALUATED":
                        continue
                    
                    payload = event.get("payload")
                    if not isinstance(payload, dict):
                        logger.warning(f"Event at line {line_num} has invalid payload")
                        continue
                    
                    # Extract attempt_id for deduplication
                    attempt_id = payload.get("attempt_id")
                    if attempt_id:
                        self.known_attempt_ids.add(str(attempt_id))
                    
                    # Reconstruct BanditState from event payload
                    arm_id = payload.get("arm_id")
                    if not arm_id:
                        logger.debug(f"Event at line {line_num} has no arm_id, skipping BanditState update")
                        continue
                    
                    # Extract scores from payload
                    axes_binned = payload.get("axes_binned", {})
                    axes_raw = payload.get("axes_raw", {})
                    cell_id = payload.get("cell_id", "")
                    yield_val = float(payload.get("successes", 0.0))
                    
                    # Reconstruct attempt_data for update_state
                    attempt_data = {
                        "cell_id": cell_id,
                        "axes_binned": axes_binned,
                        "axes_raw": axes_raw,
                        "yield": yield_val,
                        "coverage": 0.0,  # Will be recomputed by compute_arm_reward
                        "diversity": 0.0,  # Will be recomputed by compute_arm_reward
                    }
                    
                    # Compute reward and update state
                    reward = compute_arm_reward([attempt_data], self.config)
                    update_state(self.bandit_state, str(arm_id), reward, attempt_data)
                    
                    events_loaded += 1
        
        except OSError as exc:
            logger.error(f"Failed to read event log {self.event_log_path}: {exc}")
            return
        
        logger.info(f"Loaded {events_loaded} events, {len(self.known_attempt_ids)} unique attempts, "
                    f"reconstructed BanditState: {len(self.bandit_state.arm_stats)} arms")


    def _resolve_arm(self, arm_id: str) -> tuple[dict[str, object], dict[str, object], int]:
        arms = self.config.get("arms")
        if not isinstance(arms, list):
            raise ValueError("config.arms must be a list")

        arm_index = next(
            (idx for idx, arm in enumerate(arms) if isinstance(arm, Mapping) and arm.get("id") == arm_id),
            None,
        )
        if arm_index is None:
            raise ValueError(f"Unknown arm id: {arm_id}")

        paths = self._mapping(self.config.get("paths"), "paths")
        runtime = self._mapping(self.config.get("runtime"), "runtime")
        dihiggs = self._mapping(self.config.get("dihiggs"), "dihiggs")
        search = self._mapping(self.config.get("search"), "search")
        lam1_cfg = search.get("lam1", {}) if isinstance(search, Mapping) else {}
        tb_values = search.get("tb_values", []) if isinstance(search, Mapping) else []
        first_tb = tb_values[0] if isinstance(tb_values, list) and tb_values else 0

        checkpoint_root = str(Path(self.outdir) / "checkpoints")
        runtime_ctx: dict[str, object] = {
            "sys_executable": sys.executable,
            "python": sys.executable,
            "repo_root": paths.get("repo_root", ""),
            "outdir": self.outdir,
            "checkpoint_root": checkpoint_root,
            "campaign": self.config.get("campaign_id", ""),
            "threads": runtime.get("threads", 1),
            "tb_values": ",".join(str(v) for v in tb_values) if isinstance(tb_values, list) else "",
            "tb": first_tb,
            "track_id": f"{arm_id}-track-0001",
            "iter": self._round_counter.get(arm_id, 0),
            "lam1_min": lam1_cfg.get("min", 0.0),
            "lam1_max": lam1_cfg.get("max", 0.0),
            "lam1_n_bins": lam1_cfg.get("n_bins", 1),
            "hb_dataset": os.environ.get(str(dihiggs.get("hb_dataset_env", "HB_DATASET")), ""),
            "hs_dataset": os.environ.get(str(dihiggs.get("hs_dataset_env", "HS_DATASET")), ""),
        }

        resolved = substitute_placeholders(self.config, runtime_ctx)
        resolved_arms = resolved.get("arms")
        if not isinstance(resolved_arms, list):
            raise ValueError("Resolved config.arms must be a list")
        arm_cfg_raw = resolved_arms[arm_index]
        if not isinstance(arm_cfg_raw, dict):
            raise ValueError("Resolved arm config must be a dictionary")

        return arm_cfg_raw, runtime_ctx, arm_index

    @staticmethod
    def _mapping(value: object, name: str) -> Mapping[str, object]:
        if value is None:
            return {}
        if not isinstance(value, Mapping):
            raise ValueError(f"config.{name} must be a mapping")
        return value

    @staticmethod
    def make_attempt_id(campaign_id: str, arm_id: str, cell_id: str, run_dir: str) -> str:
        fingerprint = f"{campaign_id}|{arm_id}|{cell_id}|{run_dir}"
        return hashlib.sha256(fingerprint.encode("utf-8")).hexdigest()[:16]

    def run_single_round(self, arm_id: str, update_bandit: bool = False) -> dict[str, object]:
        preflight = run_all_preflight_checks(self.config)
        if preflight.get("overall") == "fail":
            logger.error("Preflight failed; skipping round for arm %s", arm_id)
            return {
                "discoveries": 0,
                "events_emitted": 0,
                "preflight": preflight,
                "subprocess_status": "skipped_preflight_fail",
            }

        arm_cfg, runtime_ctx, _ = self._resolve_arm(arm_id)

        cmd_raw = arm_cfg.get("cmd", [])
        if not isinstance(cmd_raw, list):
            raise ValueError(f"arm[{arm_id}].cmd must be a list")
        cmd = [str(item) for item in cmd_raw]

        env_raw = arm_cfg.get("env", {})
        if not isinstance(env_raw, Mapping):
            raise ValueError(f"arm[{arm_id}].env must be a mapping")
        env = {**os.environ, **{str(k): str(v) for k, v in env_raw.items()}}

        runtime_cfg = self._mapping(self.config.get("runtime"), "runtime")
        paths_cfg = self._mapping(self.config.get("paths"), "paths")

        timeout_sec = int(arm_cfg.get("timeout_sec", runtime_cfg.get("timeout_sec", 3600)))
        cwd = str(paths_cfg.get("repo_root", os.getcwd()))

        round_num = self._round_counter.get(arm_id, 0)
        stdout_path = Path(self.outdir) / f"{arm_id}_round_{round_num}.stdout.txt"
        stderr_path = Path(self.outdir) / f"{arm_id}_round_{round_num}.stderr.txt"

        subprocess_status = "success"
        try:
            result = subprocess.run(
                cmd,
                env=env,
                capture_output=True,
                timeout=timeout_sec,
                cwd=cwd,
                check=False,
            )
            stdout_path.write_text(result.stdout.decode("utf-8", errors="replace"), encoding="utf-8")
            stderr_path.write_text(result.stderr.decode("utf-8", errors="replace"), encoding="utf-8")
            if result.returncode != 0:
                subprocess_status = "nonzero_exit"
        except subprocess.TimeoutExpired as exc:
            subprocess_status = "timeout"
            logger.error("Timeout while executing arm %s after %ss", arm_id, timeout_sec)

            stdout = exc.stdout if exc.stdout is not None else b""
            stderr = exc.stderr if exc.stderr is not None else b""
            if isinstance(stdout, str):
                stdout_text = stdout
            else:
                stdout_text = stdout.decode("utf-8", errors="replace")
            if isinstance(stderr, str):
                stderr_text = stderr
            else:
                stderr_text = stderr.decode("utf-8", errors="replace")

            stdout_path.write_text(stdout_text, encoding="utf-8")
            stderr_path.write_text(stderr_text, encoding="utf-8")

        checkpoint_root = str(Path(str(runtime_ctx["checkpoint_root"])) / arm_id)
        explorer = str(arm_cfg.get("explorer", ""))
        if explorer == "adaptive":
            discoveries = parse_adaptive_checkpoint(checkpoint_root, self.config, self.known_run_dirs)
        elif explorer == "branch":
            discoveries = parse_branch_checkpoint(checkpoint_root, self.config, self.known_run_dirs)
        else:
            raise ValueError(f"Unsupported explorer type: {explorer}")

        limits = self.config.get("limits", {})
        max_new = int(limits.get("max_new_run_dirs_per_arm_call", len(discoveries))) if isinstance(limits, Mapping) else len(discoveries)
        discoveries = discoveries[:max_new]

        events_emitted = 0
        for discovery in discoveries:
            axes_binned = {"tb": int(getattr(discovery, "tb")), "lam1_bin": int(getattr(discovery, "lam1_bin"))}
            axes_raw = {
                "tb": int(getattr(discovery, "tb")),
                "lam1": float(getattr(discovery, "lam1_raw", 0.0)),
            }
            cell_id = encode_cell_id(axes_binned)
            attempt_id = self.make_attempt_id(
                str(self.config.get("campaign_id", "")), arm_id, cell_id, str(getattr(discovery, "run_dir"))
            )
            scores = self.evaluate_run_dir(str(getattr(discovery, "run_dir")), axes_binned)
            was_emitted = self.emit_attempt_event(
                attempt_id=attempt_id,
                arm_id=arm_id,
                cell_id=cell_id,
                eval_status="success",
                scores=scores,
                axes_binned=axes_binned,
                axes_raw=axes_raw,
            )
            if was_emitted:
                self.known_run_dirs.add(str(getattr(discovery, "run_dir")))
                
                if update_bandit:
                    attempt_data = {
                        "cell_id": cell_id,
                        "axes_binned": axes_binned,
                        "axes_raw": axes_raw,
                        "yield": scores["yield"],
                        "coverage": scores["coverage"],
                        "diversity": scores["diversity"],
                    }
                    reward = compute_arm_reward([attempt_data], self.config)
                    update_state(self.bandit_state, arm_id, reward, attempt_data)
                
                events_emitted += 1

        self._round_counter[arm_id] = round_num + 1
        return {
            "discoveries": len(discoveries),
            "events_emitted": events_emitted,
            "preflight": preflight,
            "subprocess_status": subprocess_status,
        }

    def evaluate_run_dir(self, run_dir: str, axes_binned: Mapping[str, object]) -> dict[str, float]:
        metrics_result = compute_metrics(run_dir, axes_binned, self.config)
        yield_val = metrics_result["yield_val"]
        
        coverage = compute_coverage(self.bandit_state.global_history, self.config)
        diversity = compute_diversity(self.bandit_state.global_history, self.config)
        composite = compute_composite(yield_val, coverage, diversity, self.config)
        
        return {
            "yield": yield_val,
            "coverage": coverage,
            "diversity": diversity,
            "composite": composite,
        }

    def emit_attempt_event(
        self,
        attempt_id: str,
        arm_id: str,
        cell_id: str,
        eval_status: str,
        scores: Mapping[str, float],
        axes_binned: Mapping[str, object],
        axes_raw: Mapping[str, object],
    ) -> bool:
        if attempt_id in self.known_attempt_ids:
            logger.info(f"Skipping duplicate attempt {attempt_id}")
            return False
        
        event = {
            "schema_version": 1,
            "campaign_id": self.config.get("campaign_id"),
            "event_type": "ATTEMPT_EVALUATED",
            "utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
            "payload": {
                "arm_id": arm_id,
                "attempt_id": attempt_id,
                "iter_index": 0,
                "attempt_index": 0,
                "cell_id": cell_id,
                "eval_status": eval_status,
                "successes": float(scores.get("yield", 0.0)),
                "trials": 1,
                "elapsed_sec": 0.0,
                "axes_binned": dict(axes_binned),
                "axes_raw": dict(axes_raw),
            },
        }
        with self.event_log_path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(event) + "\n")
        
        self.known_attempt_ids.add(attempt_id)
        return True

    def run_adaptation_round(self) -> dict[str, object]:
        arm_id = select_arm_ucb1(self.bandit_state, self.config)
        
        result = self.run_single_round(arm_id, update_bandit=True)
        
        return {
            "arm_id": arm_id,
            "discoveries": result["discoveries"],
            "events_emitted": result["events_emitted"],
            "preflight": result.get("preflight"),
            "subprocess_status": result.get("subprocess_status"),
        }

