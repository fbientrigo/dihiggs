from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class ClassificationRules:
    """Runtime interface for a future scientific decision.

    The v1 repository contract intentionally supplies no rules. Tests may
    provide explicit rules to exercise the deterministic mechanics without
    promoting those test values into the scientific contract.
    """

    mixed: dict[str, tuple[float, float]]
    photonic: dict[str, tuple[float, float]]


def _matches(values: dict[str, float], rules: dict[str, tuple[float, float]]) -> bool:
    return all(name in values and low <= values[name] <= high for name, (low, high) in rules.items())


def classify_observables(values: dict[str, Any], rules: ClassificationRules | None = None) -> dict[str, Any]:
    if rules is None:
        return {
            "family": None,
            "status": "BLOCKED_SCIENTIFIC_DECISION",
            "reason": "no_authoritative_mixed_photonic_thresholds_or_active_12520_fixtures",
        }
    numeric: dict[str, float] = {}
    for name, value in values.items():
        if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)):
            numeric[name] = float(value)
    mixed = _matches(numeric, rules.mixed)
    photonic = _matches(numeric, rules.photonic)
    if mixed and photonic:
        return {"family": None, "status": "AMBIGUOUS", "reason": "overlapping_classification_rules"}
    if not mixed and not photonic:
        return {"family": None, "status": "UNCLASSIFIED", "reason": "no_family_rule_matched"}
    return {"family": "mixed" if mixed else "photonic", "status": "CLASSIFIED", "reason": "deterministic_rule_match"}

