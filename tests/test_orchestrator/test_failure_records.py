"""
test_failure_records.py — Verify that a fake failing executable
produces correct fail records with captured stdout/stderr.

When the C++ binary exits non-zero:
- TaskResult.status == "fail"
- task_summary.jsonl has status="fail" with returncode != 0
- scan_meta.json has event="fail"
- stdout.log and stderr.log are created in the task directory
- stderr.log contains the error message printed by the binary
- stdout.log contains any stdout the binary printed
- The fail record does NOT qualify for checkpoint reuse
  (event != "done" in scan_meta.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dihiggs.app.orchestrator.engines.lambda1 import Lambda1Engine
from dihiggs.app.orchestrator.runner import ScanRunner
from tests.test_orchestrator.conftest import (
    make_fixed_lam1,
    make_grid,
)


class TestFailureRecords:

    def _make_runner(self, tmp_outdir, fake_failure_exe, tanbeta_list=None):
        return ScanRunner(
            engine=Lambda1Engine(),
            executable=fake_failure_exe,
            grid=make_grid(),
            fixed_base=make_fixed_lam1(tan_beta=0.0),
            tanbeta_list=tanbeta_list or [10000.0],
            campaign="test_fail",
            outdir=tmp_outdir,
            lake_name="lake",
            run_name="run_fail",
            dry_run=False,
        )

    def test_task_result_status_fail(self, tmp_outdir, fake_failure_exe):
        runner = self._make_runner(tmp_outdir, fake_failure_exe)
        results = runner.run()
        assert len(results) == 1
        assert results[0].status == "fail"

    def test_returncode_nonzero(self, tmp_outdir, fake_failure_exe):
        """returncode must be non-zero (1 from fake_failure_exe)."""
        runner = self._make_runner(tmp_outdir, fake_failure_exe)
        results = runner.run()
        assert results[0].returncode is not None
        assert results[0].returncode != 0

    def test_stderr_log_created(self, tmp_outdir, fake_failure_exe):
        """stderr.log must be created in the task directory."""
        runner = self._make_runner(tmp_outdir, fake_failure_exe)
        results = runner.run()

        stderr_path = results[0].stderr_path
        assert stderr_path is not None
        assert Path(stderr_path).exists(), f"stderr.log not found: {stderr_path}"

    def test_stderr_contains_error_message(self, tmp_outdir, fake_failure_exe):
        """stderr.log must contain the error message from fake_failure_exe."""
        runner = self._make_runner(tmp_outdir, fake_failure_exe)
        results = runner.run()

        stderr_content = Path(results[0].stderr_path).read_text(encoding="utf-8")
        # fake_failure_exe prints "Fatal: parameter out of range"
        assert "Fatal" in stderr_content or "parameter" in stderr_content, (
            f"stderr.log unexpected content: {stderr_content!r}"
        )

    def test_stdout_log_created(self, tmp_outdir, fake_failure_exe):
        """stdout.log must be created even when the binary fails."""
        runner = self._make_runner(tmp_outdir, fake_failure_exe)
        results = runner.run()

        stdout_path = results[0].stdout_path
        assert stdout_path is not None
        assert Path(stdout_path).exists(), f"stdout.log not found: {stdout_path}"

    def test_task_summary_has_fail_event(self, tmp_outdir, fake_failure_exe):
        runner = self._make_runner(tmp_outdir, fake_failure_exe)
        runner.run()

        jsonl = list(tmp_outdir.rglob("task_summary.jsonl"))
        assert len(jsonl) == 1

        records = [json.loads(l) for l in jsonl[0].read_text().splitlines() if l.strip()]
        fail_records = [r for r in records if r.get("status") == "fail"]
        assert len(fail_records) == 1

        rec = fail_records[0]
        assert rec["returncode"] != 0
        assert "stderr_path" in rec
        assert "stdout_path" in rec

    def test_scan_meta_has_fail_event(self, tmp_outdir, fake_failure_exe):
        runner = self._make_runner(tmp_outdir, fake_failure_exe)
        runner.run()

        metas = list(tmp_outdir.rglob("scan_meta.json"))
        assert len(metas) == 1

        meta = json.loads(metas[0].read_text())
        assert meta["event"] == "fail", f"Expected event=fail, got: {meta['event']}"

    def test_fail_not_eligible_for_resume(self, tmp_outdir, fake_failure_exe):
        """
        A fail record must NOT be treated as a done checkpoint.
        scan_meta.json event must be 'fail', not 'done'.
        """
        runner = self._make_runner(tmp_outdir, fake_failure_exe)
        runner.run()

        metas = list(tmp_outdir.rglob("scan_meta.json"))
        meta = json.loads(metas[0].read_text())
        assert meta["event"] != "done", (
            "Failed tasks must not be eligible for resume checkpoint reuse"
        )

    def test_manifest_summary_records_failure(self, tmp_outdir, fake_failure_exe):
        """run_manifest.json summary must count the task as failed."""
        runner = self._make_runner(
            tmp_outdir, fake_failure_exe, tanbeta_list=[10000.0, 20000.0]
        )
        runner.run()

        manifest = json.loads(
            next(tmp_outdir.rglob("run_manifest.json")).read_text()
        )
        summary = manifest["summary"]
        assert summary["tasks_failed"] == 2
        assert summary["tasks_ok"] == 0

    def test_mixed_tasks(self, tmp_outdir, fake_success_exe, fake_failure_exe):
        """
        Run one task with success and one with failure in separate runners.
        Verify each runner's summary is correct.
        (We can't easily mix in a single runner without more complex setup.)
        """
        # success runner: 1 task
        runner_ok = ScanRunner(
            engine=Lambda1Engine(),
            executable=fake_success_exe,
            grid=make_grid(),
            fixed_base=make_fixed_lam1(tan_beta=0.0),
            tanbeta_list=[10000.0],
            campaign="test_mixed_ok",
            outdir=tmp_outdir / "ok",
            lake_name="lake",
            run_name="run_ok",
        )
        results_ok = runner_ok.run()
        assert results_ok[0].status == "done"

        # failure runner: 1 task
        runner_fail = ScanRunner(
            engine=Lambda1Engine(),
            executable=fake_failure_exe,
            grid=make_grid(),
            fixed_base=make_fixed_lam1(tan_beta=0.0),
            tanbeta_list=[20000.0],
            campaign="test_mixed_fail",
            outdir=tmp_outdir / "fail",
            lake_name="lake",
            run_name="run_fail",
        )
        results_fail = runner_fail.run()
        assert results_fail[0].status == "fail"
