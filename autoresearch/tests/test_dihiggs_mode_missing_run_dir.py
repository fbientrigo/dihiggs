from __future__ import annotations

import json
import os
import stat
import subprocess
import tempfile
import time
import unittest
import uuid
from pathlib import Path
from unittest.mock import patch

from autoresearch.harness.dihiggs_axis_contract import DIHIGGS_EXPLORERS_MODE
from autoresearch.harness.scheduler import run_campaign


def _base_config(repo_root: Path, outdir: Path) -> dict[str, object]:
    return {
        "campaign_id": f"mode-missing-run-dir-campaign-{uuid.uuid4().hex[:8]}",
        "paths": {
            "repo_root": str(repo_root),
            "outdir": str(outdir),
            "lake_name": "events.jsonl",
        },
        "runtime": {
            "threads": 2,
            "timeout_sec": 60,
        },
        "dihiggs": {
            "phys_exec": "{repo_root}/dihiggs/app/PhysScanWithFixings",
            "hb_dataset_env": "HB_DATASET",
            "hs_dataset_env": "HS_DATASET",
        },
        "search": {
            "tb_values": [1000],
            "lam1": {"min": -20.0, "max": 20.0, "n_bins": 8},
        },
        "limits": {
            "max_new_run_dirs_per_round": 5,
            "max_new_run_dirs_per_arm_call": 5,
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


class DihiggsModeMissingRunDirTest(unittest.TestCase):
    tmpdir: tempfile.TemporaryDirectory[str] | None = None
    repo_root: Path = Path(".")
    outdir: Path = Path(".")

    def setUp(self) -> None:
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.repo_root = Path(self.tmpdir.name)
        self.outdir = self.repo_root / "out"
        (self.repo_root / "hb").mkdir(parents=True, exist_ok=True)
        (self.repo_root / "hs").mkdir(parents=True, exist_ok=True)
        binary_path = self.repo_root / "dihiggs" / "app" / "PhysScanWithFixings"
        binary_path.parent.mkdir(parents=True, exist_ok=True)
        binary_path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        _ = binary_path.chmod(binary_path.stat().st_mode | stat.S_IXUSR)

    def test_missing_run_dir_emits_crash_and_broken_artifact(self) -> None:
        config = _base_config(self.repo_root, self.outdir)

        checkpoint_root = self.outdir / "checkpoints" / "adaptive-v1" / "iter_0000"
        checkpoint_root.mkdir(parents=True, exist_ok=True)
        broken_run_dir = self.repo_root / "runs" / "broken"
        broken_run_dir.mkdir(parents=True, exist_ok=True)
        state_file = checkpoint_root / "adaptive_state.json"
        state_file.write_text(
            json.dumps(
                {
                    "proposals": [
                        {
                            "proposal_id": "missing-run-dir",
                            "status": "DONE",
                            "bin_index": 0,
                            "lam1_min": -8.0,
                            "lam1_max": -7.0,
                            "successes_by_tb": {"1000": 0},
                            "trials_by_tb": {"1000": 0},
                        },
                        {
                            "proposal_id": "broken-artifact",
                            "status": "DONE",
                            "bin_index": 1,
                            "lam1_min": -6.0,
                            "lam1_max": -5.0,
                            "run_dir": str(broken_run_dir),
                            "successes_by_tb": {"1000": 0},
                            "trials_by_tb": {"1000": 0},
                        },
                    ]
                }
            )
            + "\n",
            encoding="utf-8",
        )
        old_mtime = time.time() - 5.0
        os.utime(state_file, (old_mtime, old_mtime))

        with patch("autoresearch.harness.dihiggs_preflight.find_library", return_value="libgsl.so"), patch(
            "autoresearch.harness.dihiggs_preflight.shutil.which", return_value="/usr/bin/gfortran"
        ), patch.dict(
            os.environ,
            {"HB_DATASET": str(self.repo_root / "hb"), "HS_DATASET": str(self.repo_root / "hs")},
            clear=False,
        ), patch("autoresearch.harness.dihiggs_runner.subprocess.run") as run_mock:
            run_mock.return_value = subprocess.CompletedProcess(
                args=["python", "tool.py"],
                returncode=0,
                stdout=b"ok",
                stderr=b"",
            )
            summary = run_campaign(
                config=config,
                mode=DIHIGGS_EXPLORERS_MODE,
                iters=1,
                attempts_per_iter=1,
            )

        self.assertEqual("pass", summary["status"])

        events_path = Path(str(summary["events_path"]))
        rows = [json.loads(line) for line in events_path.read_text(encoding="utf-8").splitlines() if line.strip()]
        evaluated = [row["payload"] for row in rows if row["event_type"] == "ATTEMPT_EVALUATED"]
        self.assertEqual(2, len(evaluated))

        crash_payload = next(payload for payload in evaluated if payload.get("eval_status") == "CRASH")
        self.assertEqual(0, crash_payload["successes"])
        self.assertEqual(0, crash_payload["trials"])
        self.assertIn("source_ref", crash_payload)

        broken_payload = next(payload for payload in evaluated if payload.get("eval_status") == "BROKEN_ARTIFACT")
        self.assertEqual(0, broken_payload["successes"])
        self.assertEqual(0, broken_payload["trials"])
        self.assertEqual(str(broken_run_dir), broken_payload["run_dir"])


if __name__ == "__main__":
    _ = unittest.main()
