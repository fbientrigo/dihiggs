"""Campaign manifest creation and validation for immutable run reproducibility.

This module implements the campaign manifest system that binds configuration,
code version, dataset paths, and schema versions into an immutable contract.
Manifest mismatches block campaign resume to ensure reproducibility.

Schema Version: 1 (initial release)

Manifest Fields:
- campaign_id: Unique identifier from config
- config_hash: SHA256 of sorted JSON config subset
- runtime_ceilings: Timeout, max_empty_rounds limits
- datasets: Env var names + resolved paths for HB/HS datasets
- executable_path: PhysLam1Scan or PhysScanWithFixings binary path
- metric_schema_version: Version 1 (yield/coverage/diversity/composite)
- telemetry_schema_version: Version 1 (events.jsonl schema)
- supervisor: supervisor config block
- convergence: convergence config block
- alerts: alerts config block
- autoscaling: autoscaling config block
- created_utc: ISO 8601 timestamp of manifest creation

Validation Contract:
- campaign_id mismatch: FAIL (different campaign)
- config_hash mismatch: FAIL (config changed)
- metric_schema_version mismatch: FAIL (incompatible metrics)
- telemetry_schema_version mismatch: FAIL (incompatible events)

Validation returns:
- {"status": "pass", "reason": "..."} on success
- {"status": "fail", "reason": "..."} on mismatch
"""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping


def _compute_config_hash(config: Mapping[str, object]) -> str:
    """Compute stable SHA256 hash of config subset for compatibility checking.
    
    Includes only the fields that affect campaign reproducibility:
    - campaign_id
    - runtime (timeout, max_empty_rounds)
    - limits
    - supervisor
    - convergence
    - alerts
    - autoscaling
    
    Args:
        config: Full configuration dictionary
        
    Returns:
        64-character hex string (SHA256 digest)
    """
    # Extract reproducibility-critical subset
    subset = {
        "campaign_id": config.get("campaign_id"),
        "runtime": config.get("runtime"),
        "limits": config.get("limits"),
        "supervisor": config.get("supervisor"),
        "convergence": config.get("convergence"),
        "alerts": config.get("alerts"),
        "autoscaling": config.get("autoscaling"),
    }
    
    # Serialize with sorted keys for stability
    canonical_json = json.dumps(subset, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical_json.encode("utf-8")).hexdigest()


def _resolve_dataset_paths(config: Mapping[str, object]) -> dict[str, dict[str, str]]:
    """Resolve dataset environment variables to their current paths.
    
    Args:
        config: Configuration dictionary with dihiggs.hb_dataset_env and dihiggs.hs_dataset_env
        
    Returns:
        Dictionary mapping env var name to {env_name, resolved_path}
        
    Example:
        {
            "HB_DATASET": {
                "env_name": "HB_DATASET",
                "resolved_path": "/datasets/HBDataset"
            },
            "HS_DATASET": {
                "env_name": "HS_DATASET",
                "resolved_path": "/datasets/HSDataset"
            }
        }
    """
    dihiggs_cfg = config.get("dihiggs", {})
    if not isinstance(dihiggs_cfg, dict):
        dihiggs_cfg = {}
    
    hb_env_name = dihiggs_cfg.get("hb_dataset_env", "HB_DATASET")
    hs_env_name = dihiggs_cfg.get("hs_dataset_env", "HS_DATASET")
    
    return {
        str(hb_env_name): {
            "env_name": str(hb_env_name),
            "resolved_path": os.environ.get(str(hb_env_name), ""),
        },
        str(hs_env_name): {
            "env_name": str(hs_env_name),
            "resolved_path": os.environ.get(str(hs_env_name), ""),
        },
    }


def _resolve_phys_executable(config: Mapping[str, object]) -> str:
    """Resolve physics executable path with {repo_root} substitution.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        Resolved absolute path to PhysLam1Scan or PhysScanWithFixings
    """
    dihiggs_cfg = config.get("dihiggs", {})
    if not isinstance(dihiggs_cfg, dict):
        return ""
    
    phys_exec = dihiggs_cfg.get("phys_exec", "")
    if not isinstance(phys_exec, str):
        return ""
    
    # Substitute {repo_root} if present
    if "{repo_root}" in phys_exec:
        paths_cfg = config.get("paths", {})
        if isinstance(paths_cfg, dict):
            repo_root = paths_cfg.get("repo_root", "")
            if isinstance(repo_root, str):
                phys_exec = phys_exec.format(repo_root=repo_root)
    
    return phys_exec


def write_manifest(outdir: str | Path, config: Mapping[str, object]) -> Path:
    """Create campaign_manifest.json in the campaign output directory.
    
    This function is called on first campaign start to bind the configuration,
    executable paths, dataset paths, and schema versions into an immutable manifest.
    
    Args:
        outdir: Campaign output directory path
        config: Full configuration dictionary
        
    Returns:
        Path to the created manifest file
        
    Side Effects:
        Creates {outdir}/campaign_manifest.json with all required fields
    """
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    
    manifest = {
        "campaign_id": config.get("campaign_id", ""),
        "config_hash": _compute_config_hash(config),
        "runtime_ceilings": config.get("runtime", {}),
        "datasets": _resolve_dataset_paths(config),
        "executable_path": _resolve_phys_executable(config),
        "metric_schema_version": 1,
        "telemetry_schema_version": 1,
        "supervisor": config.get("supervisor", {}),
        "convergence": config.get("convergence", {}),
        "alerts": config.get("alerts", {}),
        "autoscaling": config.get("autoscaling", {}),
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    
    manifest_path = outdir_path / "campaign_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    
    return manifest_path


def validate_manifest(outdir: str, config: Mapping[str, object]) -> dict[str, str]:
    """Validate campaign manifest against current configuration.
    
    Checks for compatibility mismatches that would break reproducibility:
    - campaign_id must match exactly
    - config_hash must match exactly
    - metric_schema_version must match (version 1)
    - telemetry_schema_version must match (version 1)
    
    Args:
        outdir: Campaign output directory path
        config: Current configuration dictionary
        
    Returns:
        Dictionary with keys:
        - status: "pass" or "fail"
        - reason: Human-readable explanation
        
    Examples:
        {"status": "pass", "reason": "Manifest matches current config"}
        {"status": "fail", "reason": "campaign_id mismatch: expected 'A', got 'B'"}
    """
    manifest_path = Path(outdir) / "campaign_manifest.json"
    
    # Check if manifest exists
    if not manifest_path.exists():
        return {
            "status": "fail",
            "reason": f"Manifest not found at {manifest_path}",
        }
    
    # Load existing manifest
    try:
        with open(manifest_path) as f:
            manifest = json.load(f)
    except json.JSONDecodeError as exc:
        return {
            "status": "fail",
            "reason": f"Failed to parse manifest JSON: {exc}",
        }
    except Exception as exc:
        return {
            "status": "fail",
            "reason": f"Failed to read manifest: {exc}",
        }
    
    # Check campaign_id
    expected_campaign_id = config.get("campaign_id", "")
    actual_campaign_id = manifest.get("campaign_id", "")
    if expected_campaign_id != actual_campaign_id:
        return {
            "status": "fail",
            "reason": f"campaign_id mismatch: expected '{expected_campaign_id}', got '{actual_campaign_id}'",
        }
    
    # Check config_hash
    expected_hash = _compute_config_hash(config)
    actual_hash = manifest.get("config_hash", "")
    if expected_hash != actual_hash:
        return {
            "status": "fail",
            "reason": f"config_hash mismatch: config has changed (expected {expected_hash}, got {actual_hash})",
        }
    
    # Check metric_schema_version
    expected_metric_version = 1
    actual_metric_version = manifest.get("metric_schema_version")
    if expected_metric_version != actual_metric_version:
        return {
            "status": "fail",
            "reason": f"metric_schema_version mismatch: expected {expected_metric_version}, got {actual_metric_version}",
        }
    
    # Check telemetry_schema_version
    expected_telemetry_version = 1
    actual_telemetry_version = manifest.get("telemetry_schema_version")
    if expected_telemetry_version != actual_telemetry_version:
        return {
            "status": "fail",
            "reason": f"telemetry_schema_version mismatch: expected {expected_telemetry_version}, got {actual_telemetry_version}",
        }
    
    # All checks passed
    return {
        "status": "pass",
        "reason": "Manifest matches current config (campaign_id, config_hash, schema versions all valid)",
    }
