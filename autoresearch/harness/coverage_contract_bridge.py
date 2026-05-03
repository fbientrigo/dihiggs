"""Bridge reconciled MVP artifacts into ``compute_coverage``-compatible history.

The benchmark coverage math already exists in :mod:`autoresearch.benchmarks.score`.
This module does not redefine that math. It only adapts reconciled records into the
history/state shapes expected by ``compute_coverage`` and downstream MVP coverage
state upserts.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import cast

from autoresearch.benchmarks.score import compute_coverage


def _mapping_or_empty(value: object) -> Mapping[str, object]:
    if isinstance(value, Mapping):
        return cast(Mapping[str, object], value)
    return {}


def _artifact_status(record: Mapping[str, object]) -> str:
    status = record.get("artifact_status")
    if status in (None, ""):
        return "complete"
    return str(status)


def _extract_axes_binned(record: Mapping[str, object]) -> dict[str, object]:
    direct_axes = _mapping_or_empty(record.get("axes_binned"))
    if direct_axes:
        return dict(direct_axes)

    payload = _mapping_or_empty(record.get("payload"))
    payload_axes = _mapping_or_empty(payload.get("axes_binned"))
    if payload_axes:
        return dict(payload_axes)

    return {}


def _identity_value(record: Mapping[str, object], field: str) -> str | None:
    value = record.get(field)
    if value not in (None, ""):
        return str(value)

    payload = _mapping_or_empty(record.get("payload"))
    payload_value = payload.get(field)
    if payload_value in (None, ""):
        return None
    return str(payload_value)


def coverage_history_from_reconciled_records(
    records: Sequence[Mapping[str, object]],
) -> list[dict[str, dict[str, object]]]:
    """Convert reconciled records into ``compute_coverage`` history entries.

    Partial records and records with no usable ``axes_binned`` payload are treated as
    deterministic no-ops so replay and empty-history paths stay safe.
    """

    history: list[dict[str, dict[str, object]]] = []
    for record in records:
        if not all(
            (
                _identity_value(record, "campaign_id"),
                _identity_value(record, "arm_id"),
                _identity_value(record, "run_dir_fingerprint"),
                _identity_value(record, "attempt_id"),
            )
        ):
            continue

        if _artifact_status(record) == "partial":
            continue

        axes_binned = _extract_axes_binned(record)
        if not axes_binned:
            continue

        history.append({"axes_binned": axes_binned})
    return history


def coverage_state_from_reconciled_records(
    records: Sequence[Mapping[str, object]],
    config: Mapping[str, object],
) -> list[dict[str, object]]:
    """Build cumulative coverage-state rows from reconciled records.

    Each output row keeps the canonical Task 1 identity fields, the original artifact
    status, the reconciled ``axes_binned`` payload, and the cumulative coverage value
    produced by the existing benchmark math.
    """

    return extend_coverage_state([], records, config)


def extend_coverage_state(
    existing_state: Sequence[Mapping[str, object]],
    records: Sequence[Mapping[str, object]],
    config: Mapping[str, object],
) -> list[dict[str, object]]:
    """Build cumulative coverage-state rows on top of persisted coverage history."""

    state_rows: list[dict[str, object]] = []
    history = coverage_history_from_state_records(existing_state)

    for record in records:
        campaign_id = _identity_value(record, "campaign_id")
        arm_id = _identity_value(record, "arm_id")
        run_dir_fingerprint = _identity_value(record, "run_dir_fingerprint")
        attempt_id = _identity_value(record, "attempt_id")
        if not all((campaign_id, arm_id, run_dir_fingerprint, attempt_id)):
            continue

        axes_binned = _extract_axes_binned(record)
        artifact_status = _artifact_status(record)
        if artifact_status != "partial" and axes_binned:
            history.append({"axes_binned": axes_binned})

        state_rows.append(
            {
                "campaign_id": campaign_id,
                "arm_id": arm_id,
                "run_dir_fingerprint": run_dir_fingerprint,
                "attempt_id": attempt_id,
                "artifact_status": artifact_status,
                "axes_binned": axes_binned,
                "coverage_fraction": compute_coverage(history, config),
            }
        )

    return state_rows


def coverage_history_from_state_records(
    coverage_state: Sequence[Mapping[str, object]],
) -> list[dict[str, dict[str, object]]]:
    """Recover ``compute_coverage`` history from persisted coverage-state rows."""

    return coverage_history_from_reconciled_records(coverage_state)
