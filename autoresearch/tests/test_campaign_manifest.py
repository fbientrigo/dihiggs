from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from ..harness.campaign_manifest import (
    validate_manifest,
    write_manifest,
)


def _base_config() -> dict[str, object]:
    """Minimal valid config for testing."""
    return {
        "campaign_id": "test-campaign-001",
        "paths": {
            "repo_root": "/home/fabi/dihiggs",
            "outdir": "/tmp/test-outdir",
        },
        "dihiggs": {
            "phys_exec": "/home/fabi/dihiggs/dihiggs/app/PhysLam1Scan",
            "hb_dataset_env": "HB_DATASET",
            "hs_dataset_env": "HS_DATASET",
        },
        "search": {
            "tb_values": [50, 100, 500],
            "lam1": {"min": -1.0, "max": 1.0, "n_bins": 20},
        },
        "runtime": {
            "timeout_sec": 900,
            "max_empty_rounds": 2,
        },
        "limits": {
            "max_new_run_dirs_per_round": 20,
        },
        "supervisor": {
            "campaign_id": "dihiggs-lam1-exploration",
        },
        "convergence": {
            "warmup_rounds": 8,
        },
        "alerts": {
            "dedupe_window_sec": 900,
        },
        "autoscaling": {
            "min_threads": 1,
            "max_threads": 4,
        },
    }


def test_write_manifest_creates_file_with_all_required_fields() -> None:
    """Test manifest creation with all required fields."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = _base_config()
        outdir = Path(tmpdir)

        # Mock environment variables for dataset paths
        with patch.dict(os.environ, {"HB_DATASET": "/datasets/HB", "HS_DATASET": "/datasets/HS"}):
            write_manifest(str(outdir), config)

        manifest_path = outdir / "campaign_manifest.json"
        assert manifest_path.exists()

        with open(manifest_path) as f:
            manifest = json.load(f)

        # Verify all required fields exist
        assert "campaign_id" in manifest
        assert "config_hash" in manifest
        assert "runtime_ceilings" in manifest
        assert "datasets" in manifest
        assert "executable_path" in manifest
        assert "metric_schema_version" in manifest
        assert "telemetry_schema_version" in manifest
        assert "created_utc" in manifest

        # Verify specific values
        assert manifest["campaign_id"] == "test-campaign-001"
        assert manifest["metric_schema_version"] == 1
        assert manifest["telemetry_schema_version"] == 1
        assert manifest["executable_path"] == "/home/fabi/dihiggs/dihiggs/app/PhysLam1Scan"


def test_write_manifest_returns_created_path() -> None:
    """Test that write_manifest returns the created Path object."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = _base_config()
        outdir = Path(tmpdir)
        
        with patch.dict(os.environ, {"HB_DATASET": "/datasets/HB", "HS_DATASET": "/datasets/HS"}):
            returned_path = write_manifest(str(outdir), config)
        
        # Assert return value is a Path object
        assert isinstance(returned_path, Path)
        # Assert it points to campaign_manifest.json
        assert returned_path.name == "campaign_manifest.json"
        # Assert the file actually exists
        assert returned_path.exists()
        # Assert it can be opened and parsed as JSON
        with open(returned_path) as f:
            manifest = json.load(f)
        assert manifest["campaign_id"] == "test-campaign-001"


def test_write_manifest_includes_dataset_env_names_and_resolved_paths() -> None:
    """Test manifest stores both env var names and resolved paths."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = _base_config()
        outdir = Path(tmpdir)

        with patch.dict(os.environ, {"HB_DATASET": "/path/to/hb", "HS_DATASET": "/path/to/hs"}):
            write_manifest(str(outdir), config)

        manifest_path = outdir / "campaign_manifest.json"
        with open(manifest_path) as f:
            manifest = json.load(f)

        datasets = manifest["datasets"]
        assert "HB_DATASET" in datasets
        assert "HS_DATASET" in datasets
        assert datasets["HB_DATASET"]["env_name"] == "HB_DATASET"
        assert datasets["HB_DATASET"]["resolved_path"] == "/path/to/hb"
        assert datasets["HS_DATASET"]["env_name"] == "HS_DATASET"
        assert datasets["HS_DATASET"]["resolved_path"] == "/path/to/hs"


def test_write_manifest_includes_runtime_ceilings() -> None:
    """Test manifest includes runtime ceilings from config."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = _base_config()
        config["runtime"]["timeout_sec"] = 1200
        config["runtime"]["max_empty_rounds"] = 5
        outdir = Path(tmpdir)

        with patch.dict(os.environ, {"HB_DATASET": "/datasets/HB", "HS_DATASET": "/datasets/HS"}):
            write_manifest(str(outdir), config)

        manifest_path = outdir / "campaign_manifest.json"
        with open(manifest_path) as f:
            manifest = json.load(f)

        ceilings = manifest["runtime_ceilings"]
        assert ceilings["timeout_sec"] == 1200
        assert ceilings["max_empty_rounds"] == 5


def test_write_manifest_includes_supervisor_config_subset() -> None:
    """Test manifest stores supervisor, convergence, alerts, autoscaling config."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = _base_config()
        config["supervisor"] = {"campaign_id": "dihiggs-lam1-exploration", "enable_autoscaling": True}
        config["convergence"] = {"warmup_rounds": 8, "window_rounds": 5}
        config["alerts"] = {"dedupe_window_sec": 900, "channels": ["stderr"]}
        config["autoscaling"] = {"min_threads": 1, "max_threads": 4}
        outdir = Path(tmpdir)

        with patch.dict(os.environ, {"HB_DATASET": "/datasets/HB", "HS_DATASET": "/datasets/HS"}):
            write_manifest(str(outdir), config)

        manifest_path = outdir / "campaign_manifest.json"
        with open(manifest_path) as f:
            manifest = json.load(f)

        assert "supervisor" in manifest
        assert manifest["supervisor"]["campaign_id"] == "dihiggs-lam1-exploration"
        assert "convergence" in manifest
        assert manifest["convergence"]["warmup_rounds"] == 8
        assert "alerts" in manifest
        assert manifest["alerts"]["dedupe_window_sec"] == 900
        assert "autoscaling" in manifest
        assert manifest["autoscaling"]["min_threads"] == 1


def test_write_manifest_computes_stable_config_hash() -> None:
    """Test config hash is stable (sorted keys) and deterministic."""
    with tempfile.TemporaryDirectory() as tmpdir1, tempfile.TemporaryDirectory() as tmpdir2:
        config = _base_config()

        with patch.dict(os.environ, {"HB_DATASET": "/datasets/HB", "HS_DATASET": "/datasets/HS"}):
            write_manifest(tmpdir1, config)
            write_manifest(tmpdir2, config)

        with open(Path(tmpdir1) / "campaign_manifest.json") as f1:
            manifest1 = json.load(f1)
        with open(Path(tmpdir2) / "campaign_manifest.json") as f2:
            manifest2 = json.load(f2)

        # Config hash should be identical for identical configs
        assert manifest1["config_hash"] == manifest2["config_hash"]
        assert len(manifest1["config_hash"]) == 64  # SHA256 hex digest


def test_validate_manifest_pass_for_matching_config() -> None:
    """Test validation passes when manifest matches current config."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = _base_config()
        outdir = Path(tmpdir)

        with patch.dict(os.environ, {"HB_DATASET": "/datasets/HB", "HS_DATASET": "/datasets/HS"}):
            write_manifest(str(outdir), config)
            result = validate_manifest(str(outdir), config)

        assert result["status"] == "pass"
        assert "manifest matches" in result["reason"].lower()


def test_validate_manifest_fail_when_campaign_id_mismatch() -> None:
    """Test validation fails when campaign_id differs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = _base_config()
        outdir = Path(tmpdir)

        with patch.dict(os.environ, {"HB_DATASET": "/datasets/HB", "HS_DATASET": "/datasets/HS"}):
            write_manifest(str(outdir), config)

            # Change campaign_id
            config["campaign_id"] = "different-campaign-002"
            result = validate_manifest(str(outdir), config)

        assert result["status"] == "fail"
        assert "campaign_id" in result["reason"].lower()


def test_validate_manifest_fail_when_config_hash_mismatch() -> None:
    """Test validation fails when config content changes."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = _base_config()
        outdir = Path(tmpdir)

        with patch.dict(os.environ, {"HB_DATASET": "/datasets/HB", "HS_DATASET": "/datasets/HS"}):
            write_manifest(str(outdir), config)

            # Change config content
            config["runtime"]["timeout_sec"] = 2000
            result = validate_manifest(str(outdir), config)

        assert result["status"] == "fail"
        assert "config_hash" in result["reason"].lower()


def test_validate_manifest_fail_when_metric_schema_version_mismatch() -> None:
    """Test validation fails when metric schema version differs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = _base_config()
        outdir = Path(tmpdir)

        with patch.dict(os.environ, {"HB_DATASET": "/datasets/HB", "HS_DATASET": "/datasets/HS"}):
            write_manifest(str(outdir), config)

            # Manually edit manifest to change metric schema version
            manifest_path = outdir / "campaign_manifest.json"
            with open(manifest_path) as f:
                manifest = json.load(f)
            manifest["metric_schema_version"] = 2
            with open(manifest_path, "w") as f:
                json.dump(manifest, f)

            result = validate_manifest(str(outdir), config)

        assert result["status"] == "fail"
        assert "metric_schema_version" in result["reason"].lower()


def test_validate_manifest_fail_when_telemetry_schema_version_mismatch() -> None:
    """Test validation fails when telemetry schema version differs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = _base_config()
        outdir = Path(tmpdir)

        with patch.dict(os.environ, {"HB_DATASET": "/datasets/HB", "HS_DATASET": "/datasets/HS"}):
            write_manifest(str(outdir), config)

            # Manually edit manifest to change telemetry schema version
            manifest_path = outdir / "campaign_manifest.json"
            with open(manifest_path) as f:
                manifest = json.load(f)
            manifest["telemetry_schema_version"] = 2
            with open(manifest_path, "w") as f:
                json.dump(manifest, f)

            result = validate_manifest(str(outdir), config)

        assert result["status"] == "fail"
        assert "telemetry_schema_version" in result["reason"].lower()


def test_validate_manifest_fail_when_manifest_missing() -> None:
    """Test validation fails when manifest file doesn't exist."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = _base_config()
        result = validate_manifest(tmpdir, config)

        assert result["status"] == "fail"
        assert "not found" in result["reason"].lower()


def test_validate_manifest_fail_on_corrupted_json() -> None:
    """Test validation fails gracefully on corrupted manifest JSON."""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = _base_config()
        outdir = Path(tmpdir)

        # Write corrupted JSON
        manifest_path = outdir / "campaign_manifest.json"
        with open(manifest_path, "w") as f:
            f.write("{invalid json content")

        result = validate_manifest(str(outdir), config)

        assert result["status"] == "fail"
        assert "json" in result["reason"].lower() or "parse" in result["reason"].lower()
