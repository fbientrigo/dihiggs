from __future__ import annotations

import math
from typing import Any


PREFERRED = (500.0, 1000.0)
WEIGHTS = {
    "family_topology_match": 0.35,
    "anchor_lifetime_proximity": 0.35,
    "cross_mass_family_consistency": 0.25,
    "diversity_archive_contribution": 0.05,
}


def lifetime_proximity(ctau_mm: float | None) -> float:
    if ctau_mm is None or not math.isfinite(ctau_mm) or ctau_mm <= 0:
        return 0.0
    if PREFERRED[0] <= ctau_mm <= PREFERRED[1]:
        return 1.0
    return math.exp(-abs(math.log(ctau_mm / 700.0)) / math.log(2.0))


def score_components(*, valid: bool, classification: dict[str, Any], ctau_mm: float | None,
                     cross_mass_consistency: float = 0.0, diversity: float = 0.0) -> dict[str, float]:
    components = {
        "validity_gate": 1.0 if valid else 0.0,
        "family_topology_match": 1.0 if classification.get("status") == "CLASSIFIED" else 0.0,
        "anchor_lifetime_proximity": lifetime_proximity(ctau_mm),
        "cross_mass_family_consistency": max(0.0, min(1.0, float(cross_mass_consistency))),
        "diversity_archive_contribution": max(0.0, min(1.0, float(diversity))),
    }
    if not valid:
        components["total_score"] = 0.0
    else:
        components["total_score"] = sum(WEIGHTS[name] * components[name] for name in WEIGHTS)
    return components

