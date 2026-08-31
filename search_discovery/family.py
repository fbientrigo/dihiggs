"""Cross-mass family validation at mH = 150, 200, 250 GeV.

The continuation is the frozen one: Q, S, X and tan_beta are held fixed and the
other parameters are derived, so X = lambda6 * tan_beta is identical across all
three members by construction.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any

from .bounds import Envelope, check_global, BoundsError
from .evaluator import DiscoveryEvaluator, normalize
from .families import provisional_family
from .objective import lifetime_distance

VALIDATION_MASSES = (150.0, 200.0, 250.0)
ANCHOR_MASS = 200.0


def family_members(anchor_coords: dict[str, float]) -> list[dict[str, float]]:
    if abs(float(anchor_coords["mH_GeV"]) - ANCHOR_MASS) > 1e-9:
        raise ValueError("family_anchor_mass_must_be_200_GeV")
    s = anchor_coords["mA_GeV"] ** 2 - anchor_coords["mH_GeV"] ** 2
    out = []
    for mass in VALIDATION_MASSES:
        out.append({
            "mH_GeV": mass,
            "mA_GeV": math.sqrt(mass * mass + s),
            "tan_beta": anchor_coords["tan_beta"],
            "X": anchor_coords["X"],
            "Q": anchor_coords["Q"],
        })
    return out


def validate_family(anchor_coords: dict[str, float], evaluator: DiscoveryEvaluator,
                    run_dir: Path, envelope: Envelope, worker_id: str,
                    strategy_id: str, ledger=None) -> dict[str, Any]:
    members = family_members(anchor_coords)
    results = []
    for coords in members:
        try:
            check_global(**coords)
        except BoundsError as error:
            results.append({"status": "FAILED", "failure_stage": "family_continuation",
                            "failure_reason": str(error), "coordinates": coords})
            continue
        proposal = {
            "proposal_id": f"famval_{worker_id}_{int(coords['mH_GeV'])}_"
                           + hashlib.sha256(repr(sorted(coords.items())).encode()).hexdigest()[:10],
            "strategy_id": strategy_id, "worker_id": worker_id, "generation": 0,
            "parent_ids": [], "random_seed": None,
            "rationale": "frozen Q,S continuation cross-mass family validation",
            "coordinates": coords,
        }
        # the continuation may legitimately leave the *adaptive* envelope; global
        # physics bounds still apply and are checked above
        result = evaluator.evaluate(proposal, run_dir / "family_attempts" /
                                    proposal["proposal_id"], envelope=None)
        results.append(result)
        if ledger is not None:
            ledger.append({"event": "FAMILY_MEMBER_EVALUATION",
                           "lifecycle": result.get("status", "FAILED"), **result})

    terminated = [r for r in results if r.get("status") == "TERMINATED"]
    cross_mass_valid = len(terminated) == 3 and all(r.get("validity_gate") for r in terminated)
    xs = {round(float(r["coordinates"]["X"]), 12) for r in terminated} if terminated else set()
    same_x = len(xs) == 1
    labels = [(r.get("provisional_family") or {}).get("family") for r in terminated]
    consistent_family = len(set(labels)) == 1 and labels and labels[0] is not None
    anchor = next((r for r in terminated
                   if abs(r["coordinates"]["mH_GeV"] - ANCHOR_MASS) < 1e-9), None)
    ctau_200 = anchor.get("ctau_mm") if anchor else None
    ld = lifetime_distance(ctau_200)
    score = 0.0
    if cross_mass_valid and same_x:
        score = (0.35 * (1.0 if consistent_family else 0.0)
                 + 0.35 * (1.0 / (1.0 + ld))
                 + 0.30 * (1.0 if labels and labels[0] is not None else 0.0))
    ids = "|".join(str(r.get("candidate_id")) for r in results)
    return {
        "schema_version": "dihiggs.llp_discovery.family.v1",
        "family_id": "family_" + hashlib.sha256(ids.encode()).hexdigest()[:16],
        "family": labels[0] if consistent_family else None,
        "family_labels_by_mass": {str(int(r["coordinates"]["mH_GeV"])):
                                  (r.get("provisional_family") or {}).get("family") for r in terminated},
        "provisional_classification": True,
        "cross_mass_valid": cross_mass_valid,
        "same_X": same_x,
        "X": anchor_coords["X"], "Q": anchor_coords["Q"],
        "S": anchor_coords["mA_GeV"] ** 2 - anchor_coords["mH_GeV"] ** 2,
        "tan_beta": anchor_coords["tan_beta"],
        "ctau_200_mm": ctau_200,
        "ctau_by_mass_mm": {str(int(r["coordinates"]["mH_GeV"])): r.get("ctau_mm") for r in terminated},
        "anchor_lifetime_distance": ld,
        "in_target_interval": ld == 0.0,
        "consistent_family_across_masses": consistent_family,
        "score": score,
        "anchor": {"coordinates": anchor_coords},
        "members": results,
        "worker_id": worker_id, "strategy_id": strategy_id,
    }
