from __future__ import annotations
# pyright: reportImplicitRelativeImport=false, reportAny=false

import json
from pathlib import Path
from unittest.mock import patch

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
            {"id": "adaptive-v1"}
        ],
    }


def test_attempt_id_stability() -> None:
    first = DiHiggsRunner.make_attempt_id("c1", "a1", "tb=10000|bin=7", "/runs/1")
    second = DiHiggsRunner.make_attempt_id("c1", "a1", "tb=10000|bin=7", "/runs/1")
    different_cell = DiHiggsRunner.make_attempt_id("c1", "a1", "tb=10000|bin=8", "/runs/1")
    different_run = DiHiggsRunner.make_attempt_id("c1", "a1", "tb=10000|bin=7", "/runs/2")

    assert first == second
    assert first != different_cell
    assert first != different_run
    assert len(first) == 16


def test_preflight_fail_skips_subprocess_and_parsing(tmp_path: Path) -> None:
    config = _base_config(tmp_path)
    runner = DiHiggsRunner(config=config, outdir=str(tmp_path / "out"))

    with (
        patch("autoresearch.harness.dihiggs_runner.run_all_preflight_checks", return_value={"overall": "fail", "checks": []}),
        patch("autoresearch.harness.dihiggs_runner.subprocess.run") as run_mock,
    ):
        summary = runner.run_single_round("adaptive-v1")

    assert summary["subprocess_status"] == "skipped_preflight_fail"
    assert summary["discoveries"] == 0
    assert summary["events_emitted"] == 0
    run_mock.assert_not_called()


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


def test_arm_without_action_is_refused(tmp_path: Path) -> None:
    runner = DiHiggsRunner(config=_base_config(tmp_path), outdir=str(tmp_path / "out"))
    with patch("autoresearch.harness.dihiggs_runner.run_all_preflight_checks", return_value={"overall": "pass", "checks": []}):
        try:
            runner.run_single_round("adaptive-v1")
        except ValueError as exc:
            assert "no 'action'" in str(exc)
        else:
            raise AssertionError("a legacy cmd arm must be refused")
