"""
test_success_records.py — Verify that a fake successful executable
produces correct done records.

When the C++ binary exits 0 and prints Total Attempts + TRIPLE_OK_POINTS:
- TaskResult.status == "done"
- task_summary.jsonl has event="done"
- scan_meta.json has event="done"
- total_attempts is parsed from stdout
- triple_ok_points is parsed from stdout
- Output CSV exists (created by the fake executable)
- Grid signature is stored correctly in the record
- axis_metadata appears in scan_meta.json
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dihiggs.app.orchestrator.engines.lambda1 import Lambda1Engine
from dihiggs.app.orchestrator.engines.m2 import M2Engine
from dihiggs.app.orchestrator.runner import ScanRunner
from tests.test_orchestrator.conftest import (
    make_fixed_lam1,
    make_fixed_m2,
    make_grid,
)


class TestSuccessRecordsLambda1:

    def _make_runner(self, tmp_outdir, fake_success_exe, tanbeta_list=None):
        return ScanRunner(
            engine=Lambda1Engine(),
            executable=fake_success_exe,
            grid=make_grid(),
            fixed_base=make_fixed_lam1(tan_beta=0.0),
            tanbeta_list=tanbeta_list or [10000.0],
            campaign="test_success",
            outdir=tmp_outdir,
            lake_name="lake",
            run_name="run_ok",
            dry_run=False,
        )

    def test_task_result_status_done(self, tmp_outdir, fake_success_exe):
        runner = self._make_runner(tmp_outdir, fake_success_exe)
        results = runner.run()
        assert len(results) == 1
        assert results[0].status == "done"

    def test_total_attempts_parsed(self, tmp_outdir, fake_success_exe):
        """total_attempts must be 100 (from fake_success_exe stdout)."""
        runner = self._make_runner(tmp_outdir, fake_success_exe)
        results = runner.run()
        assert results[0].total_attempts == 100

    def test_triple_ok_parsed(self, tmp_outdir, fake_success_exe):
        """triple_ok_points must be '42' (from fake_success_exe stdout)."""
        runner = self._make_runner(tmp_outdir, fake_success_exe)
        results = runner.run()
        assert results[0].triple_ok_points == "42"

    def test_output_csv_created(self, tmp_outdir, fake_success_exe):
        """The fake executable creates the CSV; it must exist after the run."""
        runner = self._make_runner(tmp_outdir, fake_success_exe)
        results = runner.run()
        assert results[0].output_csv.exists(), (
            f"Output CSV not found: {results[0].output_csv}"
        )

    def test_task_summary_has_done_event(self, tmp_outdir, fake_success_exe):
        runner = self._make_runner(tmp_outdir, fake_success_exe)
        runner.run()

        jsonl = list(tmp_outdir.rglob("task_summary.jsonl"))
        assert len(jsonl) == 1

        records = [json.loads(l) for l in jsonl[0].read_text().splitlines() if l.strip()]
        done_records = [r for r in records if r.get("status") == "done"]
        assert len(done_records) == 1

        rec = done_records[0]
        assert rec["total_attempts"] == 100
        assert rec["triple_ok_points"] == "42"
        assert "grid_signature" in rec
        assert rec["engine_name"] == "lambda1"

    def test_scan_meta_has_done_event(self, tmp_outdir, fake_success_exe):
        runner = self._make_runner(tmp_outdir, fake_success_exe)
        runner.run()

        metas = list(tmp_outdir.rglob("scan_meta.json"))
        assert len(metas) == 1

        meta = json.loads(metas[0].read_text())
        assert meta["event"] == "done"
        assert "axis_metadata" in meta
        assert meta["axis_metadata"]["scan_axis"] == "lambda1"
        assert meta["total_attempts"] == 100
        assert meta["triple_ok_points"] == "42"

    def test_grid_signature_in_task_result(self, tmp_outdir, fake_success_exe):
        runner = self._make_runner(tmp_outdir, fake_success_exe)
        results = runner.run()
        assert results[0].grid_sig
        assert len(results[0].grid_sig) == 16

    def test_multiple_tanbeta_all_done(self, tmp_outdir, fake_success_exe):
        """All tan(beta) tasks should complete with status=done."""
        runner = self._make_runner(
            tmp_outdir, fake_success_exe,
            tanbeta_list=[10000.0, 20000.0, 50000.0],
        )
        results = runner.run()
        assert len(results) == 3
        assert all(r.status == "done" for r in results)

    def test_manifest_summary_updated(self, tmp_outdir, fake_success_exe):
        """run_manifest.json summary block must be filled after run."""
        runner = self._make_runner(
            tmp_outdir, fake_success_exe,
            tanbeta_list=[10000.0, 20000.0],
        )
        runner.run()

        manifest = json.loads(
            next(tmp_outdir.rglob("run_manifest.json")).read_text()
        )
        summary = manifest["summary"]
        assert summary["tasks_ok"] == 2
        assert summary["tasks_failed"] == 0
        assert summary["total_attempts_parsed"] == 200  # 100 per task


class TestSuccessRecordsM2:

    def test_m2_success_records(self, tmp_outdir, fake_success_exe):
        """M2 engine success records must include M2 axis metadata."""
        runner = ScanRunner(
            engine=M2Engine(),
            executable=fake_success_exe,
            grid=make_grid(axis_min=0.0, axis_max=500_000.0, n_axis=50),
            fixed_base=make_fixed_m2(tan_beta=0.0, lambda1=1.0),
            tanbeta_list=[50.0],
            campaign="test_m2_success",
            outdir=tmp_outdir,
            lake_name="lake",
            run_name="run_ok",
            dry_run=False,
        )
        results = runner.run()
        assert results[0].status == "done"

        metas = list(tmp_outdir.rglob("scan_meta.json"))
        assert len(metas) == 1
        meta = json.loads(metas[0].read_text())
        assert meta["axis_metadata"]["scan_axis"] == "M2"
        assert "GeV^2" in meta["axis_metadata"]["axis_units"]
