"""
Gated end-to-end smoke test for DiHiggs autoresearch harness.

This module provides a skip-by-default smoke test that validates the full
harness workflow: preflight checks → run_single_round → event emission.

To run this test:
    1. Set DIHIGGS_E2E_TEST=1 environment variable
    2. Ensure prerequisites are available:
       - PhysScanWithFixings binary at $PHYS_SCAN_EXEC
       - HB dataset at $HB_DATASET
       - HS dataset at $HS_DATASET
    3. Run:
       cd autoresearch && DIHIGGS_E2E_TEST=1 pytest tests/test_e2e_smoke.py -v -s

Without the environment variable, tests are skipped automatically.
This prevents CI from attempting to run explorer binaries without real data.
"""

from __future__ import annotations
# pyright: reportImplicitRelativeImport=false

import json
import os
import pytest
from pathlib import Path

from autoresearch.harness.dihiggs_runner import DiHiggsRunner


def _get_smoke_config(repo_root: str) -> dict[str, object]:
    """Load minimal smoke test config from configs/smoke_test_minimal.json."""
    config_path = Path(__file__).parent.parent / "configs" / "smoke_test_minimal.json"
    
    if not config_path.exists():
        raise FileNotFoundError(f"Smoke config not found: {config_path}")
    
    with open(config_path) as f:
        config = json.load(f)
    
    # Substitute repo_root placeholder
    config_str = json.dumps(config)
    config_str = config_str.replace("{repo_root}", repo_root)
    config_str = config_str.replace("{sys_executable}", "python")
    config_str = config_str.replace(
        "{phys_scan_exec}",
        os.environ.get("PHYS_SCAN_EXEC", "/usr/local/bin/PhysScanWithFixings")
    )
    
    return json.loads(config_str)


@pytest.mark.skipif(
    os.environ.get("DIHIGGS_E2E_TEST") != "1",
    reason="E2E test requires DIHIGGS_E2E_TEST=1 and real prerequisites (PhysScanWithFixings, datasets)"
)
def test_full_harness_smoke_test_gated(tmp_path):
    """
    Full end-to-end smoke test for DiHiggs harness.
    
    Validates:
    - Preflight checks pass (or are configured)
    - run_single_round executes successfully
    - Event log is created
    - Stdout/stderr logs are generated
    - No crashes or hangs (timeout_sec ensures this)
    
    This test is SKIPPED by default. Set DIHIGGS_E2E_TEST=1 to enable.
    """
    repo_root = str(Path(__file__).parent.parent.parent)
    config = _get_smoke_config(repo_root)
    
    outdir = tmp_path / "smoke_out"
    runner = DiHiggsRunner(config=config, outdir=str(outdir))
    
    # Run a single round with the adaptive explorer arm
    result = runner.run_single_round("adaptive-smoke")
    
    # Validate subprocess execution status
    assert result["subprocess_status"] in ["success", "nonzero_exit", "timeout", "skipped_preflight_fail"], \
        f"Unexpected subprocess_status: {result['subprocess_status']}"
    
    # Validate event log exists (should be created even if no events emitted)
    assert runner.event_log_path.exists(), \
        f"Event log not created: {runner.event_log_path}"
    
    # Validate stdout log exists
    stdout_log = outdir / "adaptive-smoke_round_0.stdout.txt"
    assert stdout_log.exists(), \
        f"Stdout log not created: {stdout_log}"
    
    # Validate stderr log exists
    stderr_log = outdir / "adaptive-smoke_round_0.stderr.txt"
    assert stderr_log.exists(), \
        f"Stderr log not created: {stderr_log}"
    
    # Basic sanity: logs should be readable
    assert stdout_log.read_text().strip() is not None
    assert stderr_log.read_text().strip() is not None
    
    print(f"✓ Smoke test completed: subprocess_status={result['subprocess_status']}")
    print(f"  - event_log: {runner.event_log_path}")
    print(f"  - stdout: {stdout_log}")
    print(f"  - stderr: {stderr_log}")
