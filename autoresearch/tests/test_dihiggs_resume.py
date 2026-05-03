from __future__ import annotations
# pyright: reportImplicitRelativeImport=false, reportAny=false

import json
import subprocess
from pathlib import Path
from unittest.mock import patch

from autoresearch.harness.dihiggs_parsers import AdaptiveDiscovery
from autoresearch.harness.dihiggs_runner import DiHiggsRunner


def _base_config(tmp_path: Path) -> dict[str, object]:
    return {
        "campaign_id": "campaign-test",
        "paths": {
            "repo_root": str(tmp_path),
            "lake_name": "events.jsonl",
        },
        "runtime": {
            "threads": 4,
            "timeout_sec": 60,
        },
        "dihiggs": {
            "hb_dataset_env": "HB_DATASET",
            "hs_dataset_env": "HS_DATASET",
        },
        "search": {
            "tb_values": [10000],
            "lam1": {"min": -20.0, "max": 20.0, "n_bins": 40},
        },
        "limits": {
            "max_new_run_dirs_per_arm_call": 10,
        },
        "metrics": {
            "weights": {"yield": 0.5, "coverage": 0.3, "diversity": 0.2},
            "multi_axis": {
                "collapse_axes": ["tb", "lam1_bin"],
            },
        },
        "adaptation": {
            "ucb1_exploration_constant": 2.0,
        },
        "arms": [
            {
                "id": "adaptive-v1",
                "explorer": "adaptive",
                "timeout_sec": 30,
                "cmd": ["{python}", "tool.py", "--iter", "{iter}"],
                "env": {"OMP_NUM_THREADS": "{threads}"},
            }
        ],
    }


def _create_mock_results_csv(run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    csv_path = run_dir / "results.csv"
    csv_path.write_text("success,total_events\ntrue,100\nfalse,100\n", encoding="utf-8")


def test_load_existing_state_empty_log(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    assert len(runner.known_attempt_ids) == 0


def test_load_existing_state_populates_known_attempts(tmp_path: Path) -> None:
    event_log = tmp_path / "out" / "events.jsonl"
    event_log.parent.mkdir(parents=True, exist_ok=True)
    
    event = {
        "schema_version": 1,
        "campaign_id": "campaign-test",
        "event_type": "ATTEMPT_EVALUATED",
        "utc": "2026-04-03T00:00:00Z",
        "payload": {
            "attempt_id": "abc123def456",
            "iter_index": 0,
            "attempt_index": 0,
            "cell_id": "tb=10000|bin=7",
            "eval_status": "success",
            "successes": 0.5,
            "trials": 1,
            "elapsed_sec": 0.0,
            "axes_binned": {"tb": 10000, "lam1_bin": 7},
            "axes_raw": {"tb": 10000, "lam1": -12.5},
        },
    }
    
    with event_log.open("w", encoding="utf-8") as f:
        f.write(json.dumps(event) + "\n")
    
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    assert "abc123def456" in runner.known_attempt_ids
    assert len(runner.known_attempt_ids) == 1


def test_load_existing_state_multiple_events(tmp_path: Path) -> None:
    event_log = tmp_path / "out" / "events.jsonl"
    event_log.parent.mkdir(parents=True, exist_ok=True)
    
    events = [
        {
            "schema_version": 1,
            "campaign_id": "campaign-test",
            "event_type": "ATTEMPT_EVALUATED",
            "utc": "2026-04-03T00:00:00Z",
            "payload": {
                "attempt_id": f"attempt{i:03d}",
                "iter_index": 0,
                "attempt_index": 0,
                "cell_id": "tb=10000|bin=7",
                "eval_status": "success",
                "successes": 0.5,
                "trials": 1,
                "elapsed_sec": 0.0,
                "axes_binned": {"tb": 10000, "lam1_bin": 7},
                "axes_raw": {"tb": 10000, "lam1": -12.5},
            },
        }
        for i in range(5)
    ]
    
    with event_log.open("w", encoding="utf-8") as f:
        for event in events:
            f.write(json.dumps(event) + "\n")
    
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    assert len(runner.known_attempt_ids) == 5
    assert "attempt000" in runner.known_attempt_ids
    assert "attempt004" in runner.known_attempt_ids


def test_load_existing_state_malformed_json(tmp_path: Path) -> None:
    event_log = tmp_path / "out" / "events.jsonl"
    event_log.parent.mkdir(parents=True, exist_ok=True)
    
    event_log.write_text(
        '{"schema_version": 1, "event_type": "ATTEMPT_EVALUATED", "payload": {"attempt_id": "valid001"}}\n'
        'THIS IS INVALID JSON\n'
        '{"schema_version": 1, "event_type": "ATTEMPT_EVALUATED", "payload": {"attempt_id": "valid002"}}\n',
        encoding="utf-8"
    )
    
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    assert len(runner.known_attempt_ids) == 2
    assert "valid001" in runner.known_attempt_ids
    assert "valid002" in runner.known_attempt_ids


def test_load_existing_state_empty_lines(tmp_path: Path) -> None:
    event_log = tmp_path / "out" / "events.jsonl"
    event_log.parent.mkdir(parents=True, exist_ok=True)
    
    event_log.write_text(
        '\n'
        '{"schema_version": 1, "event_type": "ATTEMPT_EVALUATED", "payload": {"attempt_id": "valid001"}}\n'
        '\n'
        '\n'
        '{"schema_version": 1, "event_type": "ATTEMPT_EVALUATED", "payload": {"attempt_id": "valid002"}}\n'
        '\n',
        encoding="utf-8"
    )
    
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    assert len(runner.known_attempt_ids) == 2


def test_load_existing_state_ignores_non_evaluated_events(tmp_path: Path) -> None:
    event_log = tmp_path / "out" / "events.jsonl"
    event_log.parent.mkdir(parents=True, exist_ok=True)
    
    events_text = (
        '{"schema_version": 1, "event_type": "ATTEMPT_STARTED", "payload": {"attempt_id": "started001"}}\n'
        '{"schema_version": 1, "event_type": "ATTEMPT_EVALUATED", "payload": {"attempt_id": "evaluated001"}}\n'
        '{"schema_version": 1, "event_type": "ATTEMPT_FAILED", "payload": {"attempt_id": "failed001"}}\n'
        '{"schema_version": 1, "event_type": "ATTEMPT_EVALUATED", "payload": {"attempt_id": "evaluated002"}}\n'
    )
    
    event_log.write_text(events_text, encoding="utf-8")
    
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    assert len(runner.known_attempt_ids) == 2
    assert "evaluated001" in runner.known_attempt_ids
    assert "evaluated002" in runner.known_attempt_ids
    assert "started001" not in runner.known_attempt_ids
    assert "failed001" not in runner.known_attempt_ids


def test_emit_skips_duplicate_attempt_id(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    runner.known_attempt_ids.add("duplicate123")
    
    runner.emit_attempt_event(
        attempt_id="duplicate123",
        arm_id="adaptive-v1",
        cell_id="tb=10000|bin=7",
        eval_status="success",
        scores={"yield": 0.5},
        axes_binned={"tb": 10000, "lam1_bin": 7},
        axes_raw={"tb": 10000, "lam1": -12.5},
    )
    
    event_log = tmp_path / "out" / "events.jsonl"
    assert not event_log.exists() or event_log.stat().st_size == 0


def test_emit_adds_new_attempt_to_known_set(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    assert "new123" not in runner.known_attempt_ids
    
    runner.emit_attempt_event(
        attempt_id="new123",
        arm_id="adaptive-v1",
        cell_id="tb=10000|bin=7",
        eval_status="success",
        scores={"yield": 0.5},
        axes_binned={"tb": 10000, "lam1_bin": 7},
        axes_raw={"tb": 10000, "lam1": -12.5},
    )
    
    assert "new123" in runner.known_attempt_ids
    
    event_log = tmp_path / "out" / "events.jsonl"
    lines = event_log.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1
    event = json.loads(lines[0])
    assert event["payload"]["attempt_id"] == "new123"


def test_restart_after_crash_no_duplicates(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    
    runner1 = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    runner1.emit_attempt_event(
        attempt_id="stable123",
        arm_id="adaptive-v1",
        cell_id="tb=10000|bin=7",
        eval_status="success",
        scores={"yield": 0.5},
        axes_binned={"tb": 10000, "lam1_bin": 7},
        axes_raw={"tb": 10000, "lam1": -12.5},
    )
    
    del runner1
    
    runner2 = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    assert "stable123" in runner2.known_attempt_ids
    
    runner2.emit_attempt_event(
        attempt_id="stable123",
        arm_id="adaptive-v1",
        cell_id="tb=10000|bin=7",
        eval_status="success",
        scores={"yield": 0.5},
        axes_binned={"tb": 10000, "lam1_bin": 7},
        axes_raw={"tb": 10000, "lam1": -12.5},
    )
    
    event_log = tmp_path / "out" / "events.jsonl"
    lines = event_log.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1


def test_restart_preserves_known_run_dirs(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    
    runner1 = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    runner1.known_run_dirs.add("/runs/0001")
    runner1.known_run_dirs.add("/runs/0002")
    
    runner1.emit_attempt_event(
        attempt_id="run0001",
        arm_id="adaptive-v1",
        cell_id="tb=10000|bin=7",
        eval_status="success",
        scores={"yield": 0.5},
        axes_binned={"tb": 10000, "lam1_bin": 7},
        axes_raw={"tb": 10000, "lam1": -12.5},
    )
    
    del runner1
    
    runner2 = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    assert "run0001" in runner2.known_attempt_ids


def test_checkpoint_stability_integration(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    checkpoint_root = tmp_path / "out" / "checkpoints" / "adaptive-v1"
    iter_dir = checkpoint_root / "iter_0000"
    iter_dir.mkdir(parents=True, exist_ok=True)
    
    state_file = iter_dir / "adaptive_state.json"
    state_data = {
        "proposals": [
            {
                "run_dir": str(tmp_path / "runs" / "0001"),
                "bin_index": 10000,
                "lam1_min": -12.5,
            }
        ]
    }
    state_file.write_text(json.dumps(state_data), encoding="utf-8")
    
    _create_mock_results_csv(tmp_path / "runs" / "0001")
    
    fake_discoveries = [
        AdaptiveDiscovery(
            run_dir=str(tmp_path / "runs" / "0001"),
            tb=10000,
            lam1_raw=-12.5,
            lam1_bin=7,
            iter_index=0,
            proposal_index=0,
        )
    ]
    
    with (
        patch("autoresearch.harness.dihiggs_runner.run_all_preflight_checks", return_value={"overall": "pass", "checks": []}),
        patch("autoresearch.harness.dihiggs_runner.parse_adaptive_checkpoint", return_value=fake_discoveries),
        patch("autoresearch.harness.dihiggs_runner.subprocess.run") as run_mock,
    ):
        run_mock.return_value = subprocess.CompletedProcess(
            args=["python", "tool.py"],
            returncode=0,
            stdout=b"",
            stderr=b"",
        )
        
        runner.run_single_round("adaptive-v1")


def test_subprocess_logging_integration(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    fake_discoveries = [
        AdaptiveDiscovery(
            run_dir=str(tmp_path / "runs" / "0001"),
            tb=10000,
            lam1_raw=-12.5,
            lam1_bin=7,
            iter_index=0,
            proposal_index=0,
        )
    ]
    
    _create_mock_results_csv(tmp_path / "runs" / "0001")
    
    with (
        patch("autoresearch.harness.dihiggs_runner.run_all_preflight_checks", return_value={"overall": "pass", "checks": []}),
        patch("autoresearch.harness.dihiggs_runner.parse_adaptive_checkpoint", return_value=fake_discoveries),
        patch("autoresearch.harness.dihiggs_runner.subprocess.run") as run_mock,
    ):
        run_mock.return_value = subprocess.CompletedProcess(
            args=["python", "tool.py"],
            returncode=0,
            stdout=b"subprocess stdout output",
            stderr=b"subprocess stderr output",
        )
        
        runner.run_single_round("adaptive-v1")
    
    stdout_log = tmp_path / "out" / "adaptive-v1_round_0.stdout.txt"
    stderr_log = tmp_path / "out" / "adaptive-v1_round_0.stderr.txt"
    
    assert stdout_log.exists()
    assert stderr_log.exists()
    assert "subprocess stdout output" in stdout_log.read_text(encoding="utf-8")
    assert "subprocess stderr output" in stderr_log.read_text(encoding="utf-8")


def test_full_resume_cycle(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    
    runner1 = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    _create_mock_results_csv(tmp_path / "runs" / "0001")
    _create_mock_results_csv(tmp_path / "runs" / "0002")
    
    fake_discoveries_1 = [
        AdaptiveDiscovery(
            run_dir=str(tmp_path / "runs" / "0001"),
            tb=10000,
            lam1_raw=-12.5,
            lam1_bin=7,
            iter_index=0,
            proposal_index=0,
        )
    ]
    
    with (
        patch("autoresearch.harness.dihiggs_runner.run_all_preflight_checks", return_value={"overall": "pass", "checks": []}),
        patch("autoresearch.harness.dihiggs_runner.parse_adaptive_checkpoint", return_value=fake_discoveries_1),
        patch("autoresearch.harness.dihiggs_runner.subprocess.run") as run_mock,
    ):
        run_mock.return_value = subprocess.CompletedProcess(
            args=["python", "tool.py"],
            returncode=0,
            stdout=b"",
            stderr=b"",
        )
        
        result1 = runner1.run_single_round("adaptive-v1")
    
    assert result1["discoveries"] == 1
    assert result1["events_emitted"] == 1
    
    event_log = tmp_path / "out" / "events.jsonl"
    events_after_first = event_log.read_text(encoding="utf-8").splitlines()
    assert len(events_after_first) == 1
    
    del runner1
    
    runner2 = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    fake_discoveries_2 = [
        AdaptiveDiscovery(
            run_dir=str(tmp_path / "runs" / "0001"),
            tb=10000,
            lam1_raw=-12.5,
            lam1_bin=7,
            iter_index=0,
            proposal_index=0,
        ),
        AdaptiveDiscovery(
            run_dir=str(tmp_path / "runs" / "0002"),
            tb=10000,
            lam1_raw=-15.0,
            lam1_bin=5,
            iter_index=0,
            proposal_index=1,
        ),
    ]
    
    with (
        patch("autoresearch.harness.dihiggs_runner.run_all_preflight_checks", return_value={"overall": "pass", "checks": []}),
        patch("autoresearch.harness.dihiggs_runner.parse_adaptive_checkpoint", return_value=fake_discoveries_2),
        patch("autoresearch.harness.dihiggs_runner.subprocess.run") as run_mock,
    ):
        run_mock.return_value = subprocess.CompletedProcess(
            args=["python", "tool.py"],
            returncode=0,
            stdout=b"",
            stderr=b"",
        )
        
        result2 = runner2.run_single_round("adaptive-v1")
    
    assert result2["discoveries"] == 2
    
    events_after_restart = event_log.read_text(encoding="utf-8").splitlines()
    
    attempt_ids = []
    for line in events_after_restart:
        event = json.loads(line)
        attempt_ids.append(event["payload"]["attempt_id"])
    
    assert len(attempt_ids) == len(set(attempt_ids)), "Found duplicate attempt_ids"
    assert len(events_after_restart) == 2
