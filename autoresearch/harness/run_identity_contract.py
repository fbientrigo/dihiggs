"""MVP run identity and state-contract validation.

This module defines the canonical run identity tuple used by the MVP contract
artifact:

    (campaign_id, arm_id, run_dir_fingerprint, attempt_id)

The artifact is intentionally small and strict. It validates the four MVP
domains used by the harness state model:

- run_registry
- discovery_records
- coverage_state
- graph_edges

Validation is fail-closed:
- missing required fields raise deterministic errors
- duplicate identity tuples raise deterministic errors
- partial records are only accepted when they still satisfy the identity
  contract; otherwise they fail with an explicit partial-run message
"""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping, Sequence


class ContractValidationError(ValueError):
    """Raised when an MVP contract artifact is malformed."""


@dataclass(frozen=True)
class RunIdentityTuple:
    campaign_id: str
    arm_id: str
    run_dir_fingerprint: str
    attempt_id: str

    def as_tuple(self) -> tuple[str, str, str, str]:
        return (
            self.campaign_id,
            self.arm_id,
            self.run_dir_fingerprint,
            self.attempt_id,
        )


@dataclass(frozen=True)
class DomainSchema:
    name: str
    required_fields: tuple[str, ...]
    unique_key_fields: tuple[str, ...]
    identity_fields: tuple[str, ...]


MVP_CONTRACT_VERSION = 1

MVP_DOMAIN_SCHEMAS: dict[str, DomainSchema] = {
    "run_registry": DomainSchema(
        name="run_registry",
        required_fields=(
            "campaign_id",
            "arm_id",
            "run_dir_fingerprint",
            "attempt_id",
            "artifact_status",
        ),
        unique_key_fields=(
            "campaign_id",
            "arm_id",
            "run_dir_fingerprint",
            "attempt_id",
        ),
        identity_fields=(
            "campaign_id",
            "arm_id",
            "run_dir_fingerprint",
            "attempt_id",
        ),
    ),
    "discovery_records": DomainSchema(
        name="discovery_records",
        required_fields=(
            "campaign_id",
            "arm_id",
            "run_dir_fingerprint",
            "attempt_id",
            "discovery_status",
        ),
        unique_key_fields=(
            "campaign_id",
            "arm_id",
            "run_dir_fingerprint",
            "attempt_id",
        ),
        identity_fields=(
            "campaign_id",
            "arm_id",
            "run_dir_fingerprint",
            "attempt_id",
        ),
    ),
    "coverage_state": DomainSchema(
        name="coverage_state",
        required_fields=(
            "campaign_id",
            "arm_id",
            "run_dir_fingerprint",
            "attempt_id",
            "coverage_fraction",
        ),
        unique_key_fields=(
            "campaign_id",
            "arm_id",
            "run_dir_fingerprint",
            "attempt_id",
        ),
        identity_fields=(
            "campaign_id",
            "arm_id",
            "run_dir_fingerprint",
            "attempt_id",
        ),
    ),
    "graph_edges": DomainSchema(
        name="graph_edges",
        required_fields=(
            "campaign_id",
            "arm_id",
            "source_attempt_id",
            "target_attempt_id",
            "edge_type",
        ),
        unique_key_fields=(
            "campaign_id",
            "arm_id",
            "source_attempt_id",
            "target_attempt_id",
            "edge_type",
        ),
        identity_fields=(
            "campaign_id",
            "arm_id",
            "source_attempt_id",
            "target_attempt_id",
            "edge_type",
        ),
    ),
}


def _require_mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ContractValidationError(f"{name} must be a mapping")
    return value


def _require_sequence(value: object, name: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise ContractValidationError(f"{name} must be a list of records")
    return value


def _identity_tuple(entry: Mapping[str, object], schema: DomainSchema) -> tuple[str, ...]:
    parts: list[str] = []
    for field in schema.identity_fields:
        value = entry.get(field)
        if value is None or value == "":
            raise ContractValidationError(
                f"{schema.name} missing required identity field: {field}"
            )
        parts.append(str(value))
    return tuple(parts)


def _validate_entry(schema: DomainSchema, entry: Mapping[str, object], index: int) -> tuple[str, ...]:
    partial_state = entry.get("artifact_status") == "partial" or entry.get("partial_run") is True

    for field in schema.required_fields:
        if field not in entry or entry.get(field) in (None, ""):
            if partial_state:
                raise ContractValidationError(
                    f"{schema.name}[{index}] is partial and missing required field: {field}"
                )
            raise ContractValidationError(
                f"{schema.name}[{index}] missing required field: {field}"
            )

    return _identity_tuple(entry, schema)


def validate_mvp_contract_artifact(artifact: Mapping[str, object]) -> dict[str, list[dict[str, object]]]:
    """Validate the strict MVP contract artifact.

    Returns a normalized copy of the four domain lists when validation passes.
    """
    artifact_map = _require_mapping(artifact, "artifact")

    normalized: dict[str, list[dict[str, object]]] = {}
    for domain_name, schema in MVP_DOMAIN_SCHEMAS.items():
        domain_value = artifact_map.get(domain_name)
        domain_records = _require_sequence(domain_value, domain_name)

        seen: set[tuple[str, ...]] = set()
        normalized_records: list[dict[str, object]] = []
        for index, raw_entry in enumerate(domain_records):
            entry = _require_mapping(raw_entry, f"{domain_name}[{index}]")
            entry_dict = dict(entry)
            identity = _validate_entry(schema, entry_dict, index)
            if identity in seen:
                raise ContractValidationError(
                    f"duplicate identity tuple in {domain_name}: {identity}"
                )
            seen.add(identity)
            normalized_records.append(entry_dict)

        normalized[domain_name] = normalized_records

    return normalized


def make_run_identity_tuple(
    campaign_id: str,
    arm_id: str,
    run_dir_fingerprint: str,
    attempt_id: str,
) -> RunIdentityTuple:
    return RunIdentityTuple(
        campaign_id=campaign_id,
        arm_id=arm_id,
        run_dir_fingerprint=run_dir_fingerprint,
        attempt_id=attempt_id,
    )
