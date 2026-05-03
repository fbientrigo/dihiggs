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


def test_run_single_round_subprocess_parse_emit(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))

    fake_discoveries = [
        AdaptiveDiscovery(
            run_dir="/runs/0001",
            tb=10000,
            lam1_raw=-12.5,
            lam1_bin=7,
            iter_index=1,
            proposal_index=2,
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
            stdout=b"hello stdout",
            stderr=b"hello stderr",
        )

        summary = runner.run_single_round("adaptive-v1")

    assert summary["discoveries"] == 1
    assert summary["events_emitted"] == 1
    assert summary["subprocess_status"] == "success"

    assert run_mock.call_count == 1
    kwargs = run_mock.call_args.kwargs
    assert kwargs["capture_output"] is True
    assert kwargs["timeout"] == 30
    assert kwargs["cwd"] == str(tmp_path)
    assert kwargs["env"]["OMP_NUM_THREADS"] == "4"

    stdout_log = (tmp_path / "out" / "adaptive-v1_round_0.stdout.txt").read_text(encoding="utf-8")
    stderr_log = (tmp_path / "out" / "adaptive-v1_round_0.stderr.txt").read_text(encoding="utf-8")
    assert "hello stdout" in stdout_log
    assert "hello stderr" in stderr_log

    events = (tmp_path / "out" / "events.jsonl").read_text(encoding="utf-8").splitlines()
    assert len(events) == 1
    event = json.loads(events[0])
    assert event["event_type"] == "ATTEMPT_EVALUATED"
    assert event["payload"]["eval_status"] == "success"
    assert event["payload"]["axes_binned"] == {"tb": 10000, "lam1_bin": 7}


def test_attempt_id_stability() -> None:
    first = DiHiggsRunner.make_attempt_id("c1", "a1", "tb=10000|bin=7", "/runs/1")
    second = DiHiggsRunner.make_attempt_id("c1", "a1", "tb=10000|bin=7", "/runs/1")
    different_cell = DiHiggsRunner.make_attempt_id("c1", "a1", "tb=10000|bin=8", "/runs/1")
    different_run = DiHiggsRunner.make_attempt_id("c1", "a1", "tb=10000|bin=7", "/runs/2")

    assert first == second
    assert first != different_cell
    assert first != different_run
    assert len(first) == 16


def test_subprocess_timeout_handled_and_logs_written(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))

    timeout_exc = subprocess.TimeoutExpired(
        cmd=["python", "tool.py"],
        timeout=30,
        output=b"partial stdout",
        stderr=b"partial stderr",
    )

    with (
        patch("autoresearch.harness.dihiggs_runner.run_all_preflight_checks", return_value={"overall": "pass", "checks": []}),
        patch("autoresearch.harness.dihiggs_runner.parse_adaptive_checkpoint", return_value=[]),
        patch("autoresearch.harness.dihiggs_runner.subprocess.run", side_effect=timeout_exc),
    ):
        summary = runner.run_single_round("adaptive-v1")

    assert summary["discoveries"] == 0
    assert summary["events_emitted"] == 0
    assert summary["subprocess_status"] == "timeout"

    stdout_log = (tmp_path / "out" / "adaptive-v1_round_0.stdout.txt").read_text(encoding="utf-8")
    stderr_log = (tmp_path / "out" / "adaptive-v1_round_0.stderr.txt").read_text(encoding="utf-8")
    assert "partial stdout" in stdout_log
    assert "partial stderr" in stderr_log


def test_preflight_fail_skips_subprocess_and_parsing(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))

    with (
        patch("autoresearch.harness.dihiggs_runner.run_all_preflight_checks", return_value={"overall": "fail", "checks": []}),
        patch("autoresearch.harness.dihiggs_runner.parse_adaptive_checkpoint") as parse_mock,
        patch("autoresearch.harness.dihiggs_runner.subprocess.run") as run_mock,
    ):
        summary = runner.run_single_round("adaptive-v1")

    assert summary["subprocess_status"] == "skipped_preflight_fail"
    assert summary["discoveries"] == 0
    assert summary["events_emitted"] == 0
    run_mock.assert_not_called()
    parse_mock.assert_not_called()


def test_emit_attempt_event_schema_and_jsonl_append(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))

    runner.emit_attempt_event(
        attempt_id="abc123",
        arm_id="test-arm-1",
        cell_id="tb=10000|bin=7",
        eval_status="success",
        scores={"yield": 1.0, "coverage": 0.5, "diversity": 0.3, "composite": 0.6},
        axes_binned={"tb": 10000, "lam1_bin": 7},
        axes_raw={"tb": 10000, "lam1": -12.5},
    )

    lines = (tmp_path / "out" / "events.jsonl").read_text(encoding="utf-8").splitlines()
    assert len(lines) == 1

    event = json.loads(lines[0])
    assert event["schema_version"] == 1
    assert event["campaign_id"] == "campaign-test"
    assert event["event_type"] == "ATTEMPT_EVALUATED"

    payload = event["payload"]
    assert payload["attempt_id"] == "abc123"
    assert payload["iter_index"] == 0
    assert payload["attempt_index"] == 0
    assert payload["cell_id"] == "tb=10000|bin=7"
    assert payload["eval_status"] == "success"
    assert payload["successes"] == 1.0
    assert payload["trials"] == 1
    assert payload["elapsed_sec"] == 0.0
    assert payload["axes_binned"] == {"tb": 10000, "lam1_bin": 7}
    assert payload["axes_raw"] == {"tb": 10000, "lam1": -12.5}


def test_evaluate_run_dir_returns_mvp_scores(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))
    
    # Create a mock run_dir with results.csv
    run_dir = tmp_path / "runs" / "0001"
    run_dir.mkdir(parents=True)
    csv_path = run_dir / "results.csv"
    csv_path.write_text("success,total_events\nTrue,1000\n")

    scores = runner.evaluate_run_dir(str(run_dir), {"tb": 10000, "lam1_bin": 7})
    
    # With empty history, coverage and diversity should be 0
    # Yield should be 1.0 (1 success / 1 trial)
    assert scores["yield"] == 1.0
    assert scores["coverage"] >= 0.0
    assert scores["diversity"] >= 0.0
    assert scores["composite"] >= 0.0
