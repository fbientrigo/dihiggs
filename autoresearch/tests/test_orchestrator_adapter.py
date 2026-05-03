import json
import pytest
from autoresearch.harness.orchestrator_adapter import DryRunOrchestratorAdapter, validate_proposer_payload

def test_validate_proposer_payload():
    # Valid payload
    payload = {
        "status": "success",
        "proposals": [
            {
                "goal_id": "goal::1",
                "axes_binned": {"tb": 1000, "lam1_bin": 5},
                "rationale": {"confidence": "high"}
            }
        ]
    }
    assert validate_proposer_payload(payload) == []

    # Missing status
    assert "Missing required key: 'status'" in validate_proposer_payload({"proposals": []})
    
    # Invalid status
    assert "Invalid status: unknown" in validate_proposer_payload({"status": "unknown", "proposals": []})
    
    # Missing proposals
    assert "Missing required key: 'proposals'" in validate_proposer_payload({"status": "success"})
    
    # Invalid proposals type
    assert "'proposals' must be a sequence" in validate_proposer_payload({"status": "success", "proposals": "not a list"})
    
    # Invalid proposal elements
    invalid_proposal_payload = {
        "status": "success",
        "proposals": [
            {"goal_id": "missing_axes"},
            {"axes_binned": {}},
            "not a mapping"
        ]
    }
    errors = validate_proposer_payload(invalid_proposal_payload)
    assert any("Proposal at index 0 missing 'axes_binned'" in e for e in errors)
    assert any("Proposal at index 1 missing 'goal_id'" in e for e in errors)
    assert any("Proposal at index 2 is not a mapping" in e for e in errors)

def test_dry_run_adapter_dispatch_success():
    adapter = DryRunOrchestratorAdapter()
    payload = {
        "status": "success",
        "proposals": [
            {
                "goal_id": "goal::tb=1000|bin=1",
                "axes_binned": {"tb": 1000, "lam1_bin": 1},
                "rationale": {"confidence": "high"}
            },
            {
                "goal_id": "goal::tb=5000|bin=2",
                "axes_binned": {"tb": 5000, "lam1_bin": 2},
                "rationale": {"confidence": "medium"}
            },
            {
                "goal_id": "goal::tb=10000|bin=3",
                "axes_binned": {"tb": 10000, "lam1_bin": 3},
                "rationale": {"confidence": "low"}
            }
        ]
    }
    
    result = adapter.dispatch(payload)
    
    assert result["status"] == "success"
    assert result["dispatch_count"] == 3
    assert len(result["records"]) == 3
    for record in result["records"]:
        assert record["status"] == "dry_run"
        assert "goal_id" in record
        assert "dispatched_at" in record
        assert "axes_binned" in record["metadata"]

def test_dry_run_adapter_invalid_payload():
    adapter = DryRunOrchestratorAdapter()
    payload = {"status": "success"} # Missing proposals
    
    result = adapter.dispatch(payload)
    
    assert result["status"] == "error"
    assert "errors" in result
    assert len(result["errors"]) > 0
