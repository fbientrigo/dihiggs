from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class OrchestratorAdapter(Protocol):

    def dispatch(self, payload: Mapping[str, Any]) -> Mapping[str, Any]:
        ...


@dataclass(frozen=True)
class DispatchRecord:
    goal_id: str
    dispatched_at: str
    metadata: Mapping[str, Any] = field(default_factory=dict)
    metadata: Mapping[str, Any] = field(default_factory=dict)
    status: str = "dry_run"


def validate_proposer_payload(payload: Mapping[str, Any]) -> list[str]:
    errors = []
    if "status" not in payload:
        errors.append("Missing required key: 'status'")
    
    status = payload.get("status")
    if status not in ("success", "saturated"):
        errors.append(f"Invalid status: {status}")
        
    if "proposals" not in payload:
        errors.append("Missing required key: 'proposals'")
    else:
        proposals = payload.get("proposals")
        if not isinstance(proposals, Sequence) or isinstance(proposals, (str, bytes)):
            errors.append("'proposals' must be a sequence")
        else:
            for i, p in enumerate(proposals):
                if not isinstance(p, Mapping):
                    errors.append(f"Proposal at index {i} is not a mapping")
                    continue
                if "goal_id" not in p:
                    errors.append(f"Proposal at index {i} missing 'goal_id'")
                if "axes_binned" not in p:
                    errors.append(f"Proposal at index {i} missing 'axes_binned'")
                    
    return errors


class DryRunOrchestratorAdapter:

    def dispatch(self, payload: Mapping[str, Any]) -> Mapping[str, Any]:
        errors = validate_proposer_payload(payload)
        if errors:
            return {
                "status": "error",
                "errors": errors,
                "timestamp": datetime.now(timezone.utc).isoformat(),
            }

        proposals = payload.get("proposals", [])
        dispatch_records = []
        
        for proposal in proposals:
            record = DispatchRecord(
                goal_id=str(proposal.get("goal_id")),
                status="dry_run",
                dispatched_at=datetime.now(timezone.utc).isoformat(),
                metadata={
                    "axes_binned": proposal.get("axes_binned"),
                    "rationale": proposal.get("rationale"),
                }
            )
            dispatch_records.append({
                "goal_id": record.goal_id,
                "status": record.status,
                "dispatched_at": record.dispatched_at,
                "metadata": record.metadata
            })

        return {
            "status": "success",
            "dispatch_count": len(dispatch_records),
            "records": dispatch_records,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
