"""Frozen discovery objective. Workers optimize this and may not redefine it."""
from __future__ import annotations

import math
from typing import Any

TARGET_CTAU_MM = (500.0, 1000.0)
TARGET_CTAU_CENTRE_MM = 700.0
INVALID_PENALTY = 100.0

OBJECTIVE_ID = "discovery_objective_v1"


def lifetime_distance(ctau_mm: float | None) -> float:
    """0.0 inside the target interval; log2-octave distance outside it."""
    if ctau_mm is None or not math.isfinite(ctau_mm) or ctau_mm <= 0.0:
        return INVALID_PENALTY
    if TARGET_CTAU_MM[0] <= ctau_mm <= TARGET_CTAU_MM[1]:
        return 0.0
    return abs(math.log(ctau_mm / TARGET_CTAU_CENTRE_MM)) / math.log(2.0)


def family_distance(target_family: str, provisional: dict[str, Any]) -> float:
    """0.0 when the provisional label already matches the requested family."""
    pf = provisional.get("photon_fraction")
    ff = provisional.get("fermionic_fraction")
    if pf is None or ff is None:
        return INVALID_PENALTY
    if target_family == "photonic":
        return max(0.0, 0.5 - float(pf)) / 0.5
    if target_family == "mixed":
        return max(0.0, 0.5 - float(ff)) / 0.5
    raise ValueError(f"unknown_target_family:{target_family}")


def max_photon_objective(result: dict[str, Any], require_window: bool = False) -> dict[str, Any]:
    """Maximise photon fraction SUBJECT TO deterministic validity.

    Answers the quantitative question the target-directed objective cannot: how
    photon-rich can a theory-VALID point get? The lifetime term is dropped (or
    optional) so the optimizer spends its whole budget on photon fraction instead
    of satisfying the cheap lifetime constraint first. Lower is better.
    """
    if result.get("status") != "TERMINATED":
        return {"objective_id": "max_photon_valid_v1", "value": INVALID_PENALTY * 3,
                "reason": "evaluation_failed", "feasible": False}
    if not result.get("validity_gate"):
        failed = sorted(name for name, ok in (result.get("gates") or {}).items() if not ok)
        return {"objective_id": "max_photon_valid_v1", "value": INVALID_PENALTY * 2,
                "reason": "deterministically_invalid:" + ",".join(failed), "feasible": False}
    pf = (result.get("provisional_family") or {}).get("photon_fraction")
    if pf is None:
        return {"objective_id": "max_photon_valid_v1", "value": INVALID_PENALTY,
                "reason": "missing_branchings", "feasible": False}
    value = 1.0 - float(pf)
    extra = lifetime_distance(result.get("ctau_mm")) if require_window else 0.0
    return {"objective_id": "max_photon_valid_v1", "value": value + extra,
            "photon_fraction": float(pf), "lifetime_distance": extra, "feasible": True,
            "is_photonic": float(pf) >= 0.5}


def objective(target_family: str, result: dict[str, Any]) -> dict[str, Any]:
    """Lower is better. Deterministic invalidity is a HARD constraint."""
    if result.get("status") != "TERMINATED":
        return {"objective_id": OBJECTIVE_ID, "value": INVALID_PENALTY * 3,
                "reason": "evaluation_failed", "feasible": False}
    if not result.get("validity_gate"):
        failed = sorted(name for name, ok in (result.get("gates") or {}).items() if not ok)
        return {"objective_id": OBJECTIVE_ID, "value": INVALID_PENALTY * 2,
                "reason": "deterministically_invalid:" + ",".join(failed), "feasible": False}
    ld = lifetime_distance(result.get("ctau_mm"))
    fd = family_distance(target_family, result.get("provisional_family") or {})
    return {"objective_id": OBJECTIVE_ID, "value": ld + fd, "lifetime_distance": ld,
            "family_distance": fd, "feasible": True,
            "in_target_interval": ld == 0.0, "family_matched": fd == 0.0}
