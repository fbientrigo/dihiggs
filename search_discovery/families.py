"""PROVISIONAL family assignment on raw deterministic observables.

IMPORTANT SCOPE NOTE
--------------------
The frozen contract (dihiggs.llp_benchmark_search.v1) still carries
classification.status = BLOCKED_SCIENTIFIC_DECISION: no scientific owner has
supplied authoritative numerical mixed/photonic predicates or active-125.20
fixtures. Nothing in this module changes that, and nothing here is written back
into the frozen contract or the frozen fail-closed archive.

What this module provides is an explicitly PROVISIONAL, deterministic, purely
observable-driven operational split, so that the discovery campaign can proceed
while the categorical question is unresolved -- as the campaign brief directs.
Every record it produces is stamped provisional=True and carries the predicate
that produced it, so a scientific owner can re-derive or overrule the labels
without re-running the search.

The split is not invented from the descriptive x~1500 / x~10000 labels. It is
measured: the frozen evaluator produces a sharply bimodal branching structure,
    photon-dominated : br_gammagamma + br_Zgamma ~ 0.685 + 0.315 ~ 1.00
    fermion-dominated: br_bb + br_gg + br_tautau ~ 0.675 + 0.222 + 0.071 ~ 0.97
with essentially nothing between the modes. The 0.5 cut sits in the empty gap.
"""
from __future__ import annotations

from typing import Any

PHOTONIC_FRACTION_CUT = 0.5
MIXED_FRACTION_CUT = 0.5
PREDICATE_ID = "provisional_photon_fraction_v1"


def photon_fraction(observables: dict[str, Any]) -> float | None:
    gg = observables.get("br_gammagamma")
    zg = observables.get("br_Zgamma")
    if gg is None or zg is None:
        return None
    return float(gg) + float(zg)


def fermionic_fraction(observables: dict[str, Any]) -> float | None:
    parts = [observables.get(name) for name in ("br_bb", "br_gg", "br_tautau", "br_cc")]
    if any(part is None for part in parts):
        return None
    return float(sum(float(part) for part in parts))


def provisional_family(observables: dict[str, Any], validity_gate: bool) -> dict[str, Any]:
    """Deterministic provisional label. Never a scientific verdict."""
    pf = photon_fraction(observables)
    ff = fermionic_fraction(observables)
    record: dict[str, Any] = {
        "predicate_id": PREDICATE_ID,
        "provisional": True,
        "authoritative": False,
        "frozen_contract_status": "BLOCKED_SCIENTIFIC_DECISION",
        "photon_fraction": pf,
        "fermionic_fraction": ff,
        "validity_gate": bool(validity_gate),
    }
    if pf is None or ff is None:
        record.update({"family": None, "status": "UNDETERMINED",
                       "reason": "missing_branching_observables"})
        return record
    if pf >= PHOTONIC_FRACTION_CUT:
        record.update({"family": "photonic", "status": "PROVISIONALLY_CLASSIFIED",
                       "reason": f"photon_fraction_{pf:.4f}_ge_{PHOTONIC_FRACTION_CUT}"})
    elif ff >= MIXED_FRACTION_CUT:
        record.update({"family": "mixed", "status": "PROVISIONALLY_CLASSIFIED",
                       "reason": f"fermionic_fraction_{ff:.4f}_ge_{MIXED_FRACTION_CUT}"})
    else:
        record.update({"family": None, "status": "UNCLASSIFIED",
                       "reason": f"photon_fraction_{pf:.4f}_fermionic_fraction_{ff:.4f}_both_below_cut"})
    return record
