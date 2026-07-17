from __future__ import annotations

import json
import shutil
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, cast
from unittest.mock import patch

import pytest

pytestmark = pytest.mark.skip(reason="autoresearch is frozen outside the canonical evaluation core")

from autoresearch.harness.dihiggs_runner import DiHiggsRunner


def _base_config(tmp_path: Path, *, campaign_id: str) -> dict[str, object]:
    return {
        "campaign_id": campaign_id,
        "paths": {
            "repo_root": str(tmp_path),
            "lake_name": "events.jsonl",
        },
        "runtime": {
            "threads": 1,
            "timeout_sec": 30,
            "max_rounds": 3,
        },
        "search": {
            "tb_values": [1000, 5000],
            "lam1": {"min": -20.0, "max": 20.0, "n_bins": 4},
        },
        "metrics": {
            "weights": {"yield_norm": 1.0, "coverage_norm": 0.0, "diversity_norm": 0.0},
            "floors": {"yield_norm": 0.0, "coverage_norm": 0.0, "diversity_norm": 0.0},
            "multi_axis": {
                "coverage_axes": [
                    {"name": "tb", "kind": "categorical", "domain": [1000, 5000], "weight": 1.0},
                ],
                "diversity_axes": [{"name": "tb", "weight": 1.0}],
            },
        },
        "adaptation": {
            "ucb_c": 1.0,
            "warm_start_each_arm": True,
        },
        "arms": [
            {"id": "arm-a", "explorer": "adaptive", "cmd": ["python", "a.py"]},
            {"id": "arm-b", "explorer": "branch", "cmd": ["python", "b.py"]},
        ],
    }


class DiHiggsUcbSelectionTests(unittest.TestCase):
    tempdir: TemporaryDirectory[str] | None = None
    tmp_path: Path = Path(".")
    campaign_id: str = ""
    state_root: Path = Path(".")

    def setUp(self) -> None:
        self.tempdir = TemporaryDirectory()
        self.tmp_path = Path(self.tempdir.name)
        self.campaign_id = f"ucb-test-{self._testMethodName}"
        self.state_root = Path("autoresearch") / "opencode_logs" / self.campaign_id
        if self.state_root.exists():
            shutil.rmtree(self.state_root)

    def tearDown(self) -> None:
        if self.tempdir is not None:
            self.tempdir.cleanup()
        if self.state_root.exists():
            shutil.rmtree(self.state_root)

    def test_warm_start_pulls_each_arm_once_before_ucb(self) -> None:
        runner = DiHiggsRunner(_base_config(self.tmp_path, campaign_id=self.campaign_id), outdir=str(self.tmp_path / "out"))
        summaries: list[dict[str, object]] = [
            {"arm_id": "arm-a", "discoveries": 0, "events_emitted": 0, "preflight": {}, "subprocess_status": "success"},
            {"arm_id": "arm-b", "discoveries": 0, "events_emitted": 0, "preflight": {}, "subprocess_status": "success"},
        ]
        scores = [0.0, 0.4, 0.4, 0.7]

        with (
            patch.object(DiHiggsRunner, "run_single_round", side_effect=summaries),
            patch.object(DiHiggsRunner, "_compute_campaign_score", side_effect=scores),
        ):
            first = runner.run_adaptation_round()
            second = runner.run_adaptation_round()

        self.assertEqual(first["arm_id"], "arm-a")
        self.assertEqual(second["arm_id"], "arm-b")

        state = cast(dict[str, Any], json.loads(runner.state_file_path.read_text(encoding="utf-8")))
        self.assertEqual(state["arms"]["arm-a"]["n"], 1)
        self.assertEqual(state["arms"]["arm-b"]["n"], 1)

    def test_ucb_prefers_arm_with_higher_bonus_after_warm_start(self) -> None:
        config = _base_config(self.tmp_path, campaign_id=self.campaign_id)
        config["adaptation"] = {"ucb_c": 1.5, "warm_start_each_arm": False}
        runner = DiHiggsRunner(config, outdir=str(self.tmp_path / "out"))
        runner._bandit_state = {
            "total_rounds": 11,
            "arms": {
                "arm-a": {"n": 10, "mean_reward": 0.8},
                "arm-b": {"n": 1, "mean_reward": 0.5},
            },
        }
        runner._save_mode_state()

        with (
            patch.object(DiHiggsRunner, "run_single_round", return_value={"arm_id": "arm-b", "discoveries": 0, "events_emitted": 0, "preflight": {}, "subprocess_status": "success"}),
            patch.object(DiHiggsRunner, "_compute_campaign_score", side_effect=[0.2, 0.6]),
        ):
            summary = runner.run_adaptation_round()

        self.assertEqual(summary["arm_id"], "arm-b")
        self.assertAlmostEqual(cast(float, summary["score_delta"]), 0.4)
        self.assertAlmostEqual(cast(float, summary["reward"]), 0.7)

        state = cast(dict[str, Any], json.loads(runner.state_file_path.read_text(encoding="utf-8")))
        self.assertEqual(state["total_rounds"], 12)
        self.assertEqual(state["arms"]["arm-b"]["n"], 2)
        self.assertAlmostEqual(cast(float, state["arms"]["arm-b"]["mean_reward"]), 0.6)

    def test_negative_score_delta_maps_to_bounded_reward(self) -> None:
        config = _base_config(self.tmp_path, campaign_id=self.campaign_id)
        config["adaptation"] = {"ucb_c": 0.0, "warm_start_each_arm": False}
        runner = DiHiggsRunner(config, outdir=str(self.tmp_path / "out"))
        runner._bandit_state = {
            "total_rounds": 3,
            "arms": {
                "arm-a": {"n": 2, "mean_reward": 0.75},
                "arm-b": {"n": 2, "mean_reward": 0.25},
            },
        }
        runner._save_mode_state()

        with (
            patch.object(DiHiggsRunner, "run_single_round", return_value={"arm_id": "arm-a", "discoveries": 0, "events_emitted": 0, "preflight": {}, "subprocess_status": "success"}),
            patch.object(DiHiggsRunner, "_compute_campaign_score", side_effect=[0.9, 0.3]),
        ):
            summary = runner.run_adaptation_round()

        self.assertEqual(summary["arm_id"], "arm-a")
        self.assertAlmostEqual(cast(float, summary["score_delta"]), -0.6)
        self.assertAlmostEqual(cast(float, summary["reward"]), 0.2)

        state = cast(dict[str, Any], json.loads(runner.state_file_path.read_text(encoding="utf-8")))
        self.assertEqual(state["total_rounds"], 4)
        self.assertEqual(state["arms"]["arm-a"]["n"], 3)
        self.assertAlmostEqual(cast(float, state["arms"]["arm-a"]["mean_reward"]), (0.75 * 2 + 0.2) / 3)

    def test_run_uses_max_rounds_and_persists_mean_reward(self) -> None:
        runner = DiHiggsRunner(_base_config(self.tmp_path, campaign_id=self.campaign_id), outdir=str(self.tmp_path / "out"))
        round_summaries: list[dict[str, object]] = [
            {"arm_id": "arm-a", "discoveries": 1, "events_emitted": 4, "preflight": {}, "subprocess_status": "success"},
            {"arm_id": "arm-b", "discoveries": 1, "events_emitted": 4, "preflight": {}, "subprocess_status": "success"},
            {"arm_id": "arm-a", "discoveries": 1, "events_emitted": 4, "preflight": {}, "subprocess_status": "success"},
        ]
        scores = [0.1, 0.5, 0.5, 0.7, 0.7, 1.0]

        with (
            patch.object(DiHiggsRunner, "run_single_round", side_effect=round_summaries),
            patch.object(DiHiggsRunner, "_compute_campaign_score", side_effect=scores),
        ):
            result = runner.run(max_rounds=3)

        self.assertEqual(result["total_rounds"], 3)
        rounds = cast(list[dict[str, object]], result["rounds"])
        self.assertEqual([round_summary["arm_id"] for round_summary in rounds], ["arm-a", "arm-b", "arm-a"])

        state = cast(dict[str, Any], json.loads(runner.state_file_path.read_text(encoding="utf-8")))
        self.assertEqual(state["total_rounds"], 3)
        self.assertEqual(state["arms"]["arm-a"]["n"], 2)
        self.assertEqual(state["arms"]["arm-b"]["n"], 1)
        self.assertAlmostEqual(cast(float, state["arms"]["arm-a"]["mean_reward"]), 0.675)
        self.assertAlmostEqual(cast(float, state["arms"]["arm-b"]["mean_reward"]), 0.6)


if __name__ == "__main__":
    _ = unittest.main()
