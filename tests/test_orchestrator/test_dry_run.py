"""
test_dry_run.py — Verify dry-run mode writes manifests without executing C++.

Invariants:
- run_manifest.json is created in the run directory.
- task_summary.jsonl is created with one "dry_run" event per tan(beta).
- scan_meta.json is created per task with event="dry_run".
- No subprocess is spawned (the fake executable is never called).
- Commands are included in the task_summary.jsonl records.
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


class TestDryRunLambda1:
    """Dry-run with Lambda1Engine."""

    def test_manifest_written(self, tmp_outdir, fake_success_exe):
        """run_manifest.json must be created even in dry-run."""
        runner = ScanRunner(
            engine=Lambda1Engine(),
            executable=fake_success_exe,  # won't be called
            grid=make_grid(),
            fixed_base=make_fixed_lam1(tan_beta=0.0),
            tanbeta_list=[10000.0, 20000.0],
            campaign="test_dry",
            outdir=tmp_outdir,
            lake_name="lake",
            run_name="run_test",
            dry_run=True,
        )
        results = runner.run()

        # Find run_manifest.json
        manifests = list(tmp_outdir.rglob("run_manifest.json"))
        assert len(manifests) == 1, f"Expected 1 manifest, found: {manifests}"

        manifest = json.loads(manifests[0].read_text())
        assert manifest["runtime"]["dry_run"] is True
        assert manifest["runtime"]["engine_name"] == "lambda1"
        # axis_metadata must be present
        assert "axis_metadata" in manifest
        assert manifest["axis_metadata"]["scan_axis"] == "lambda1"

    def test_task_summary_written_with_dry_run_events(self, tmp_outdir, fake_success_exe):
        """task_summary.jsonl must have one dry_run event per tan(beta)."""
        tanbeta_list = [10000.0, 20000.0, 50000.0]
        runner = ScanRunner(
            engine=Lambda1Engine(),
            executable=fake_success_exe,
            grid=make_grid(),
            fixed_base=make_fixed_lam1(tan_beta=0.0),
            tanbeta_list=tanbeta_list,
            campaign="test_dry",
            outdir=tmp_outdir,
            lake_name="lake",
            run_name="run_test",
            dry_run=True,
        )
        runner.run()

        jsonl_files = list(tmp_outdir.rglob("task_summary.jsonl"))
        assert len(jsonl_files) == 1

        records = [json.loads(line) for line in jsonl_files[0].read_text().splitlines() if line.strip()]
        assert len(records) == len(tanbeta_list), (
            f"Expected {len(tanbeta_list)} records, got {len(records)}"
        )
        for rec in records:
            assert rec["status"] == "dry_run", f"Expected dry_run, got: {rec['status']}"
            assert "cmd" in rec or "cmd" in rec.get("extra", {}), (
                "Dry-run records must include the command"
            )
            assert "grid_signature" in rec
            assert "tanbeta" in rec

    def test_scan_meta_written_per_task(self, tmp_outdir, fake_success_exe):
        """Each task must have a scan_meta.json with event=dry_run."""
        runner = ScanRunner(
            engine=Lambda1Engine(),
            executable=fake_success_exe,
            grid=make_grid(),
            fixed_base=make_fixed_lam1(tan_beta=0.0),
            tanbeta_list=[10000.0],
            campaign="test_dry",
            outdir=tmp_outdir,
            lake_name="lake",
            run_name="run_test",
            dry_run=True,
        )
        runner.run()

        metas = list(tmp_outdir.rglob("scan_meta.json"))
        assert len(metas) == 1

        meta = json.loads(metas[0].read_text())
        assert meta["event"] == "dry_run"
        assert "axis_metadata" in meta
        assert "engine_name" in meta
        assert "grid_signature" in meta

    def test_no_csv_created_in_dry_run(self, tmp_outdir, fake_success_exe):
        """No CSV output files must be created in dry-run mode."""
        runner = ScanRunner(
            engine=Lambda1Engine(),
            executable=fake_success_exe,
            grid=make_grid(),
            fixed_base=make_fixed_lam1(tan_beta=0.0),
            tanbeta_list=[10000.0, 20000.0],
            campaign="test_dry",
            outdir=tmp_outdir,
            lake_name="lake",
            run_name="run_test",
            dry_run=True,
        )
        runner.run()

        csvs = list(tmp_outdir.rglob("*.csv"))
        assert csvs == [], f"No CSVs should be created in dry-run; found: {csvs}"

    def test_results_have_dry_run_status(self, tmp_outdir, fake_success_exe):
        """ScanRunner.run() must return TaskResults with status=dry_run."""
        runner = ScanRunner(
            engine=Lambda1Engine(),
            executable=fake_success_exe,
            grid=make_grid(),
            fixed_base=make_fixed_lam1(tan_beta=0.0),
            tanbeta_list=[10000.0, 20000.0],
            campaign="test_dry",
            outdir=tmp_outdir,
            lake_name="lake",
            run_name="run_test",
            dry_run=True,
        )
        results = runner.run()
        assert all(r.status == "dry_run" for r in results)


class TestDryRunM2:
    """Dry-run with M2Engine."""

    def test_manifest_has_m2_engine(self, tmp_outdir, fake_success_exe):
        """Manifest must record engine_name=m2 and axis_metadata.scan_axis=M2."""
        runner = ScanRunner(
            engine=M2Engine(),
            executable=fake_success_exe,
            grid=make_grid(axis_min=0.0, axis_max=500_000.0, n_axis=50),
            fixed_base=make_fixed_m2(tan_beta=0.0, lambda1=1.0),
            tanbeta_list=[50.0],
            campaign="test_m2_dry",
            outdir=tmp_outdir,
            lake_name="lake",
            run_name="run_test",
            dry_run=True,
        )
        runner.run()

        manifests = list(tmp_outdir.rglob("run_manifest.json"))
        assert len(manifests) == 1

        manifest = json.loads(manifests[0].read_text())
        assert manifest["runtime"]["engine_name"] == "m2"
        assert manifest["axis_metadata"]["scan_axis"] == "M2"
        assert "GeV^2" in manifest["axis_metadata"]["axis_units"]
