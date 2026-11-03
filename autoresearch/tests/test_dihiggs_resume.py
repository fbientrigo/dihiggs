from __future__ import annotations
# pyright: reportImplicitRelativeImport=false, reportAny=false

import json
from pathlib import Path

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
            {"id": "adaptive-v1"}
        ],
    }


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
