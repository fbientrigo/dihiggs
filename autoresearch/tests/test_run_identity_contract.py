from __future__ import annotations

import copy
from typing import cast

import pytest

from autoresearch.harness.run_identity_contract import (
    ContractValidationError,
    make_run_identity_tuple,
    validate_mvp_contract_artifact,
)


def _contract_fixture() -> dict[str, list[dict[str, object]]]:
    return {
        "run_registry": [
            {
                "campaign_id": "camp-001",
                "arm_id": "adaptive-v1",
                "run_dir_fingerprint": "runfp-001",
                "attempt_id": "attempt-001",
                "artifact_status": "complete",
            }
        ],
        "discovery_records": [
            {
                "campaign_id": "camp-001",
                "arm_id": "adaptive-v1",
                "run_dir_fingerprint": "runfp-001",
                "attempt_id": "attempt-001",
                "artifact_status": "complete",
                "discovery_status": "new",
            }
        ],
        "coverage_state": [
            {
                "campaign_id": "camp-001",
                "arm_id": "adaptive-v1",
                "run_dir_fingerprint": "runfp-001",
                "attempt_id": "attempt-001",
                "artifact_status": "complete",
                "coverage_fraction": 0.5,
            }
        ],
        "graph_edges": [
            {
                "campaign_id": "camp-001",
                "arm_id": "adaptive-v1",
                "source_attempt_id": "attempt-001",
                "target_attempt_id": "attempt-002",
                "edge_type": "follows",
            }
        ],
    }


def test_make_run_identity_tuple_is_canonical() -> None:
    identity = make_run_identity_tuple("camp-001", "adaptive-v1", "runfp-001", "attempt-001")
    assert identity.as_tuple() == ("camp-001", "adaptive-v1", "runfp-001", "attempt-001")


def test_validate_contract_artifact_accepts_complete_fixture() -> None:
    artifact = _contract_fixture()
    validated = validate_mvp_contract_artifact(artifact)

    assert set(validated) == {"run_registry", "discovery_records", "coverage_state", "graph_edges"}
    assert validated["run_registry"][0]["attempt_id"] == "attempt-001"


def test_validate_contract_artifact_rejects_duplicate_identity_tuple() -> None:
    artifact = _contract_fixture()
    run_registry = cast(list[dict[str, object]], artifact["run_registry"])
    artifact["run_registry"] = run_registry + [copy.deepcopy(run_registry[0])]

    with pytest.raises(ContractValidationError, match="duplicate identity tuple.*run_registry"):
        validate_mvp_contract_artifact(artifact)


def test_validate_contract_artifact_rejects_partial_run_missing_attempt_identity() -> None:
    artifact = _contract_fixture()
    discovery_records = cast(list[dict[str, object]], artifact["discovery_records"])
    partial_record = copy.deepcopy(discovery_records[0])
    partial_record["artifact_status"] = "partial"
    partial_record.pop("attempt_id")
    artifact["discovery_records"] = [partial_record]

    with pytest.raises(ContractValidationError, match="partial.*attempt_id"):
        validate_mvp_contract_artifact(artifact)


def test_validate_contract_artifact_rejects_malformed_domain_payload() -> None:
    artifact = _contract_fixture()
    artifact["graph_edges"] = [
        {
            "campaign_id": "camp-001",
            "arm_id": "adaptive-v1",
            "source_attempt_id": "attempt-001",
            "target_attempt_id": "attempt-002",
        }
    ]

    with pytest.raises(ContractValidationError, match="missing required field: edge_type"):
        validate_mvp_contract_artifact(artifact)
