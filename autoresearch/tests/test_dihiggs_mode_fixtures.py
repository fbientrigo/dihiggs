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
        "campaign_id": f"mode-fixture-campaign-{uuid.uuid4().hex[:8]}",
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
            "tb_values": [1000, 2000, 3000],
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
            },
            {
                "id": "branch-v1",
                "explorer": "branch",
                "timeout_sec": 30,
                "cmd": ["{python}", "branch.py"],
                "env": {"OMP_NUM_THREADS": "{threads}"},
            },
        ],
    }


def _write_run_artifacts(run_dir: Path, *, valid_rows: int = 2) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest = run_dir / "run_manifest.json"
    manifest.write_text(json.dumps({"run_dir": str(run_dir)}) + "\n", encoding="utf-8")
    csv_path = run_dir / "results.csv"
    csv_path.write_text(
        "positivity_ok,unitarity_ok,perturbativity_ok\n" + "\n".join(["1,1,1"] * valid_rows) + "\n",
        encoding="utf-8",
    )
    summary_path = run_dir / "task_summary.jsonl"
    records = [
        {
            "event": "done",
            "triple_ok_points": valid_rows,
            "attempts": valid_rows,
            "elapsed_sec": 1.25,
            "output_csv": str(csv_path),
        }
    ]
    summary_path.write_text("".join(json.dumps(row) + "\n" for row in records), encoding="utf-8")


class DihiggsModeFixturesTest(unittest.TestCase):
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

    def test_mode_runner_writes_fixture_events_and_duplicate_skip_reused(self) -> None:
        config = _base_config(self.repo_root, self.outdir)

        checkpoint_root = self.outdir / "checkpoints" / "adaptive-v1" / "iter_0000"
        run_dir = self.repo_root / "runs" / "shared"
        _write_run_artifacts(run_dir)
        checkpoint_root.mkdir(parents=True, exist_ok=True)
        state_file = checkpoint_root / "adaptive_state.json"
        state_file.write_text(
            json.dumps(
                {
                    "proposals": [
                        {
                            "proposal_id": "p0",
                            "status": "DONE",
                            "bin_index": 0,
                            "lam1_min": -12.5,
                            "lam1_max": -11.5,
                            "run_dir": str(run_dir),
                            "successes_by_tb": {"1000": 2, "2000": 1},
                            "trials_by_tb": {"1000": 2, "2000": 1},
                        },
                        {
                            "proposal_id": "p1",
                            "status": "DONE",
                            "bin_index": 1,
                            "lam1_min": -10.0,
                            "lam1_max": -9.0,
                            "run_dir": str(run_dir),
                            "successes_by_tb": {"3000": 1},
                            "trials_by_tb": {"3000": 1},
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
        self.assertEqual("adaptive-v1", summary["arm_id"])
        self.assertEqual(2, summary["discoveries"])

        events_path = Path(str(summary["events_path"]))
        rows = [json.loads(line) for line in events_path.read_text(encoding="utf-8").splitlines() if line.strip()]
        evaluated = [row for row in rows if row["event_type"] == "ATTEMPT_EVALUATED"]
        self.assertEqual(3, len(evaluated))

        first = evaluated[0]["payload"]
        self.assertEqual({"tb": 1000, "lam1_bin": 1}, first["axes_binned"])
        self.assertEqual(-12.5, first["axes_raw"]["lam1"])
        self.assertEqual("DONE", first["eval_status"])
        self.assertIn("source_ref", first)
        self.assertIn("fingerprint_sha256", first)

        duplicate = evaluated[-1]["payload"]
        self.assertTrue(duplicate["is_duplicate"])
        self.assertEqual("SKIP_REUSED", duplicate["eval_status"])
        self.assertEqual(0, duplicate["successes"])
        self.assertEqual(0, duplicate["trials"])

        attempt_ids = [row["payload"]["attempt_id"] for row in evaluated]
        self.assertEqual(len(attempt_ids), len(set(attempt_ids)))

        summary_json = json.loads((self.outdir / "reports" / "summary.json").read_text(encoding="utf-8"))
        self.assertEqual(3, summary_json["attempts"])


if __name__ == "__main__":
    _ = unittest.main()
