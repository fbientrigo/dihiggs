from __future__ import annotations

from autoresearch.harness.orchestrator_adapter import (
    DryRunOrchestratorAdapter,
    build_orchestrator_command,
    normalize_runtime_config,
    validate_executable_proposal,
    validate_proposer_payload,
)


def _proposal() -> dict[str, object]:
    return {
        "proposal_id": "p1",
        "operator": "autoresearch",
        "bounds": {
            "lambda6": [0.0027, 0.0027],
            "tan_beta": [50000, 60000],
            "lambda1": [1.0, 1.0],
            "mH": [180.0, 220.0],
        },
        "scan_values": {"tan_beta": [50000, 55000, 60000]},
        "fixed": {"sin_ba": 1.0, "lambda7": 0.0, "mA": 450.0, "mHp": 450.0, "yukawa_type": 1},
        "resolution": {"mH": 40, "lambda1": 1},
    }


def test_validate_proposer_payload() -> None:
    payload = {"status": "success", "proposals": [{"goal_id": "goal::1", "axes_binned": {"tb": 1000, "lam1_bin": 5}, "rationale": {"confidence": "high"}}]}
    assert validate_proposer_payload(payload) == []
    assert "Missing required key: 'status'" in validate_proposer_payload({"proposals": []})
    assert "Invalid status: unknown" in validate_proposer_payload({"status": "unknown", "proposals": []})
    assert "Missing required key: 'proposals'" in validate_proposer_payload({"status": "success"})


def test_dry_run_adapter_dispatch_success() -> None:
    adapter = DryRunOrchestratorAdapter()
    payload = {"status": "success", "proposals": [{"goal_id": "goal::tb=1000|bin=1", "axes_binned": {"tb": 1000, "lam1_bin": 1}, "rationale": {"confidence": "high"}}]}
    result = adapter.dispatch(payload)
    assert result["status"] == "success"
    assert result["dispatch_count"] == 1
    assert result["records"][0]["status"] == "dry_run"


def test_normalize_runtime_config_defaults_threads_to_1() -> None:
    cfg = normalize_runtime_config(None)
    assert cfg["cpp_omp_threads"] == 1
    assert cfg["threads"] == 1


def test_validate_executable_proposal_accepts_minimal_contract() -> None:
    assert validate_executable_proposal(_proposal()) == []


def test_build_orchestrator_command_includes_required_args_and_no_force_by_default() -> None:
    cmd = build_orchestrator_command(
        _proposal(),
        {
            "exec_path": "/home/fabi/dihiggs/dihiggs/app/PhysScanWithFixings",
            "outdir": "/home/fabi/dihiggs/scripts/out",
            "lake_name": "dihiggs_lake",
            "campaign": "unit_campaign",
            "cpp_omp_threads": 4,
        },
        dry_run=False,
        force=False,
    )
    required = {
        "--exec",
        "--campaign",
        "--run-name",
        "--sin-ba",
        "--lambda6",
        "--lambda7",
        "--mA",
        "--tanbeta",
        "--mphi-min",
        "--mphi-max",
        "--n-mphi",
        "--lam1-min",
        "--lam1-max",
        "--n-lam1",
        "--threads",
    }
    assert required.issubset(set(cmd))
    assert "dihiggs/app/orchestrate_scans.py" in cmd
    assert "--threads" in cmd and cmd[cmd.index("--threads") + 1] == "4"
    assert "--force" not in cmd


def test_build_orchestrator_command_lambda1_fixed_contract() -> None:
    cmd = build_orchestrator_command(_proposal(), {"cpp_omp_threads": 1}, dry_run=False, force=False)
    assert cmd[cmd.index("--lam1-min") + 1] in {"1.0", "1"}
    assert cmd[cmd.index("--lam1-max") + 1] in {"1.0", "1"}
    assert cmd[cmd.index("--n-lam1") + 1] == "1"
