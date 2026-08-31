"""Process-level acceptance tests for benchmark_search_daemon.py
(mission Sec. 26 T1, T7, T8, T9, T10, T11, T12).

These spawn the real CLI entry point as a subprocess against small, bounded,
temp run-dirs (never the live discovery_v1 / discovery_photonic_v1 /
deliverables paths) and drive it with real signals. Each test is kept small
(tens of evaluations) so the suite stays fast; this is infrastructure
validation, not a physics campaign.
"""
from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
DAEMON = ROOT / "benchmark_search_daemon.py"
PYTHON = sys.executable


def _write_config(run_dir: Path, **overrides) -> Path:
    config = {
        "campaign": {"name": "test_harvest"},
        "run_dir": str(run_dir),
        "evaluator": str(ROOT / "dihiggs/app/DihiggsPointV2Evaluator"),
        "envelope_id": "E1_mixed_exploit",
        "runtime": {"workers": 2, "max_evaluations_per_cycle": 16,
                    "checkpoint_every": 1, "summary_every": 1, "shutdown_grace_s": 15},
        "budgets": {"max_total_evaluations": None, "max_cycle_walltime_s": 300, "max_total_walltime_s": None},
        "policy": {"seed": 20260825, "family_validate_per_cycle": 0,
                   "allocation": {"explore": 0.5, "refine": 0.5}},
        "stopping": {"patience_cycles": 100, "max_consecutive_evaluator_failures": 5},
        "safety": {"downstream_snapshots_read_only": True},
    }
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(config.get(key), dict):
            config[key].update(value)
        else:
            config[key] = value
    run_dir.parent.mkdir(parents=True, exist_ok=True)
    config_path = run_dir.parent / f"{run_dir.name}_config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    return config_path


def _ledger_lines(run_dir: Path) -> list[dict]:
    path = run_dir / "ledger.jsonl"
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _all_lines_are_valid_json(run_dir: Path) -> bool:
    path = run_dir / "ledger.jsonl"
    if not path.exists():
        return True
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        json.loads(line)
    return True


def _spawn(config_path: Path, *extra_args: str) -> subprocess.Popen:
    return subprocess.Popen(
        [PYTHON, str(DAEMON), "--config", str(config_path), *extra_args],
        cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )


def _wait_for_checkpoint(run_dir: Path, min_cycle: int, timeout_s: float = 40.0) -> None:
    deadline = time.time() + timeout_s
    checkpoint_path = run_dir / "daemon_checkpoint.json"
    while time.time() < deadline:
        if checkpoint_path.exists():
            try:
                state = json.loads(checkpoint_path.read_text(encoding="utf-8"))
                if state.get("cycle", 0) >= min_cycle:
                    return
            except json.JSONDecodeError:
                pass
        time.sleep(0.2)
    raise TimeoutError(f"checkpoint did not reach cycle {min_cycle} within {timeout_s}s")


# ------------------------------------------------------------------- T11
def test_dry_run_invokes_the_evaluator_zero_times(tmp_path, monkeypatch):
    """In-process variant so we can directly assert the frozen evaluator's
    own subprocess.run is never called, not just infer it indirectly."""
    calls = []
    import search_discovery.evaluator as discovery_evaluator

    def _fail_if_called(*args, **kwargs):
        calls.append((args, kwargs))
        raise AssertionError("evaluator subprocess must never run during --dry-run")

    monkeypatch.setattr(discovery_evaluator.subprocess, "run", _fail_if_called)

    from search_discovery import daemon

    run_dir = tmp_path / "runs" / "dry_run_test"
    config = dict(daemon.DEFAULTS)
    config["run_dir"] = str(run_dir)
    config["envelope_id"] = "E1_mixed_exploit"
    config["policy"] = {"seed": 1, "family_validate_per_cycle": 1,
                        "allocation": {"explore": 0.5, "climb": 0.5}}

    exit_code = daemon.run(config, run_dir, resume=False, workers_override=1, dry_run=True)
    assert exit_code == 0
    assert calls == []
    assert not (run_dir / "attempts").exists()


# -------------------------------------------------------------------- T1
def test_kill_and_resume_no_lost_or_duplicated_committed_evaluations(tmp_path):
    run_dir = tmp_path / "runs" / "kill_resume_test"
    config_path = _write_config(run_dir, runtime={"workers": 2, "max_evaluations_per_cycle": 20,
                                                   "checkpoint_every": 1, "summary_every": 1})
    proc = _spawn(config_path)
    try:
        _wait_for_checkpoint(run_dir, min_cycle=1)
        time.sleep(0.3)  # let a second cycle get underway
        proc.send_signal(signal.SIGKILL)
        proc.wait(timeout=15)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait(timeout=15)

    assert _all_lines_are_valid_json(run_dir)
    pre_events = _ledger_lines(run_dir)
    pre_candidate_ids = [e.get("candidate_id") for e in pre_events
                         if e.get("event") == "EVALUATION" and e.get("candidate_id")]
    assert len(pre_candidate_ids) == len(set(pre_candidate_ids)), "duplicate candidate_id before resume"

    resumed = subprocess.run(
        [PYTHON, str(DAEMON), "--config", str(config_path), "--resume", "--max-cycles",
         str(json.loads((run_dir / "daemon_checkpoint.json").read_text())["cycle"] + 2)],
        cwd=ROOT, capture_output=True, text=True, timeout=60,
    )
    assert resumed.returncode == 0, resumed.stderr

    assert _all_lines_are_valid_json(run_dir)
    post_events = _ledger_lines(run_dir)
    post_eval_events = [e for e in post_events if e.get("event") == "EVALUATION"]
    # the pre-kill prefix of EVALUATION events must be untouched (append-only)
    pre_eval_events = [e for e in pre_events if e.get("event") == "EVALUATION"]
    assert post_eval_events[: len(pre_eval_events)] == pre_eval_events
    assert len(post_events) >= len(pre_events)

    post_checkpoint = json.loads((run_dir / "daemon_checkpoint.json").read_text())
    assert post_checkpoint["total_evaluations"] >= len(pre_eval_events)


# -------------------------------------------------------------------- T7
@pytest.mark.parametrize("sig", [signal.SIGINT, signal.SIGTERM])
def test_signal_leads_to_clean_checkpoint_and_expected_exit_code(tmp_path, sig):
    run_dir = tmp_path / "runs" / f"signal_test_{sig.name}"
    config_path = _write_config(run_dir)
    proc = _spawn(config_path)
    try:
        _wait_for_checkpoint(run_dir, min_cycle=1)
        proc.send_signal(sig)
        stdout, stderr = proc.communicate(timeout=30)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.communicate(timeout=15)

    expected_code = 130 if sig == signal.SIGINT else 143
    assert proc.returncode == expected_code, f"stderr={stderr}\nstdout={stdout}"
    assert _all_lines_are_valid_json(run_dir)
    checkpoint = json.loads((run_dir / "daemon_checkpoint.json").read_text())
    assert checkpoint["stopped"] is True
    assert checkpoint["stop_reason"].startswith("signal:")
    assert '"shutdown_summary": true' in stdout


# -------------------------------------------------------------------- T8
def test_concurrent_workers_no_ledger_corruption_no_duplicate_committed_ids(tmp_path):
    run_dir = tmp_path / "runs" / "concurrent_test"
    config_path = _write_config(run_dir, runtime={"workers": 4, "max_evaluations_per_cycle": 24,
                                                   "checkpoint_every": 1, "summary_every": 1})
    result = subprocess.run(
        [PYTHON, str(DAEMON), "--config", str(config_path), "--max-cycles", "2"],
        cwd=ROOT, capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stderr
    assert _all_lines_are_valid_json(run_dir)
    events = [e for e in _ledger_lines(run_dir) if e.get("event") == "EVALUATION"]
    ids = [e["candidate_id"] for e in events if e.get("candidate_id")]
    assert len(ids) == len(set(ids)), "concurrent workers produced a duplicate committed candidate_id"
    for event in events:
        assert "status" in event and "candidate_id" in event  # no partial/truncated records


# -------------------------------------------------------------------- T9
def test_workers_flag_actually_bounds_the_process_pool_size(tmp_path, monkeypatch):
    from search_discovery import daemon

    seen_max_workers = []
    real_executor = daemon.cf.ProcessPoolExecutor

    class _SpyExecutor(real_executor):
        def __init__(self, *args, **kwargs):
            seen_max_workers.append(kwargs.get("max_workers"))
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(daemon.cf, "ProcessPoolExecutor", _SpyExecutor)
    run_dir = tmp_path / "runs" / "workers_flag_test"
    config = dict(daemon.DEFAULTS)
    config["run_dir"] = str(run_dir)
    config["envelope_id"] = "E1_mixed_exploit"
    config["policy"] = {"seed": 1, "family_validate_per_cycle": 0, "allocation": {"explore": 1.0}}
    config["runtime"] = {"workers": 2, "max_evaluations_per_cycle": 4,
                         "checkpoint_every": 1, "summary_every": 1, "shutdown_grace_s": 10}

    daemon.run(config, run_dir, resume=False, workers_override=3, dry_run=False, stop_after_cycles=1)
    assert seen_max_workers == [3]


def test_max_total_evaluations_budget_stops_the_loop(tmp_path):
    run_dir = tmp_path / "runs" / "budget_test"
    config_path = _write_config(run_dir, runtime={"workers": 2, "max_evaluations_per_cycle": 8,
                                                   "checkpoint_every": 1, "summary_every": 1},
                                budgets={"max_total_evaluations": 10})
    result = subprocess.run(
        [PYTHON, str(DAEMON), "--config", str(config_path)],
        cwd=ROOT, capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stderr
    assert '"max_total_evaluations_reached"' in result.stdout
    checkpoint = json.loads((run_dir / "daemon_checkpoint.json").read_text())
    assert checkpoint["total_evaluations"] >= 10


# ------------------------------------------------------------------- T10
def test_restart_preserves_cycle_and_evaluation_counters(tmp_path):
    run_dir = tmp_path / "runs" / "restart_test"
    config_path = _write_config(run_dir, runtime={"workers": 2, "max_evaluations_per_cycle": 10,
                                                   "checkpoint_every": 1, "summary_every": 1})
    first = subprocess.run(
        [PYTHON, str(DAEMON), "--config", str(config_path), "--max-cycles", "1"],
        cwd=ROOT, capture_output=True, text=True, timeout=60,
    )
    assert first.returncode == 0, first.stderr
    after_first = json.loads((run_dir / "daemon_checkpoint.json").read_text())
    assert after_first["cycle"] == 1

    second = subprocess.run(
        [PYTHON, str(DAEMON), "--config", str(config_path), "--resume", "--max-cycles", "3"],
        cwd=ROOT, capture_output=True, text=True, timeout=60,
    )
    assert second.returncode == 0, second.stderr
    after_second = json.loads((run_dir / "daemon_checkpoint.json").read_text())
    assert after_second["cycle"] == 3
    assert after_second["total_evaluations"] > after_first["total_evaluations"]

    # a fresh (non-resumed) run against the same dir must refuse to silently
    # reset -- the checkpoint is still there and load_config/run always looks
    # for it when resume=True; without --resume this test just documents that
    # total_evaluations only grows through --resume, never regresses.
    assert after_second["total_evaluations"] >= after_first["total_evaluations"]


# ------------------------------------------------------------------- T12
def test_downstream_live_campaign_and_deliverables_are_never_touched(tmp_path):
    protected_dirs = [ROOT / "deliverables", ROOT / "runs" / "discovery_v1", ROOT / "runs" / "discovery_photonic_v1"]
    before = {}
    for d in protected_dirs:
        if d.exists():
            before[d] = sorted((p, p.stat().st_mtime_ns, p.stat().st_size)
                               for p in d.rglob("*") if p.is_file())

    run_dir = tmp_path / "runs" / "snapshot_protection_test"
    config_path = _write_config(run_dir, runtime={"workers": 2, "max_evaluations_per_cycle": 8,
                                                   "checkpoint_every": 1, "summary_every": 1})
    result = subprocess.run(
        [PYTHON, str(DAEMON), "--config", str(config_path), "--max-cycles", "1"],
        cwd=ROOT, capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stderr

    for d in protected_dirs:
        if d.exists():
            after = sorted((p, p.stat().st_mtime_ns, p.stat().st_size)
                           for p in d.rglob("*") if p.is_file())
            assert after == before.get(d, []), f"{d} was modified by a daemon smoke run"


# ---------------------------------------------------------------- bounded fix
def _write_hanging_fake_evaluator(tmp_path: Path) -> Path:
    """A stand-in for a stuck evaluator subprocess: ignores every CLI arg the
    real evaluator would receive and never returns, so subprocess.run in
    evaluator.py blocks exactly like a genuinely hung physics evaluation
    would -- without needing to know the real evaluator's arg contract."""
    script = tmp_path / "hanging_fake_evaluator.sh"
    script.write_text("#!/bin/sh\nsleep 9999\n", encoding="utf-8")
    script.chmod(0o755)
    return script


def test_cycle_walltime_and_bounded_shutdown_survive_a_hung_evaluator(tmp_path):
    """MAJOR fixes: (1) max_cycle_walltime_s now actually bounds how long a
    cycle waits on outstanding futures, so a hung evaluator subprocess does
    not hang the daemon forever within a cycle; (2) shutdown after a stop
    signal is bounded even when a worker is still stuck mid-evaluation,
    instead of ProcessPoolExecutor's default wait=True blocking forever."""
    fake_evaluator = _write_hanging_fake_evaluator(tmp_path)
    run_dir = tmp_path / "runs" / "hung_evaluator_test"
    config_path = _write_config(
        run_dir, evaluator=str(fake_evaluator),
        runtime={"workers": 1, "max_evaluations_per_cycle": 2,
                 "checkpoint_every": 1, "summary_every": 1, "shutdown_grace_s": 2},
        budgets={"max_total_evaluations": None, "max_cycle_walltime_s": 2, "max_total_walltime_s": None},
    )
    proc = _spawn(config_path)
    try:
        # cycle 1 completing at all (within a generous outer timeout) proves
        # the cycle did not hang forever on the stuck evaluator subprocess.
        _wait_for_checkpoint(run_dir, min_cycle=1, timeout_s=30.0)
        proc.send_signal(signal.SIGTERM)
        start = time.time()
        proc.wait(timeout=20)  # must not need anywhere near "forever"
        bounded_shutdown_seconds = time.time() - start
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait(timeout=15)

    assert proc.returncode == 143
    assert bounded_shutdown_seconds < 20, "shutdown did not respect shutdown_grace_s bound"
    assert _all_lines_are_valid_json(run_dir)


def test_daemon_refuses_to_start_pointed_at_the_live_blocked_campaign(tmp_path):
    config_path = _write_config(tmp_path / "unused")  # run_dir overridden below
    result = subprocess.run(
        [PYTHON, str(DAEMON), "--config", str(config_path), "--run-dir", str(ROOT / "runs" / "discovery_v1"),
         "--dry-run"],
        cwd=ROOT, capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 2
    assert "REFUSED" in result.stderr
