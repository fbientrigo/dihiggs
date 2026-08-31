"""Constrained continuation: ascend photon fraction along the FEASIBLE manifold.

Unlike search_discovery.climb (single-point ascent at the mH=200 anchor only),
this module evaluates the FULL 150/200/250 family at every step, tracks
per-mass photon fraction, and supports branching when two structurally
different directions both improve. It never accepts an intermediate point
that violates any hard condition -- there is no "point toward a better
answer through invalid territory".

Free coordinates: mA_GeV, X, Q. mH is pinned to the 200 GeV anchor by the
campaign contract. tan_beta is not a free search direction; it is always
RE-SOLVED deterministically (by the frozen evaluator, confirmed at each
iterate) so ctau_200 stays in [500, 1000] mm.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any

from .bounds import Envelope, check_global, BoundsError
from .evaluator import DiscoveryEvaluator
from .family import family_members

TARGET = (500.0, 1000.0)
CENTRE = 700.0
MIN_S = 501.0
FREE_COORDS = ("mA_GeV", "X", "Q")


def _pf(observables: dict[str, Any]) -> float | None:
    gg = observables.get("br_gammagamma")
    zg = observables.get("br_Zgamma")
    if gg is None or zg is None:
        return None
    return float(gg) + float(zg)


def _eval_member(coords: dict[str, float], evaluator: DiscoveryEvaluator, ledger, run_dir: Path,
                  worker_id: str, tag: str) -> dict[str, Any]:
    pid = f"{tag}_" + hashlib.sha256(repr(sorted(coords.items())).encode()).hexdigest()[:12]
    proposal = {"proposal_id": pid, "strategy_id": "constrained_photon_continuation",
                "worker_id": worker_id, "generation": 0, "parent_ids": [],
                "random_seed": None, "rationale": "feasible-manifold continuation step",
                "coordinates": dict(coords)}
    result = evaluator.evaluate(proposal, run_dir / "attempts" / pid, None)
    ledger.append({"event": "EVALUATION", "lifecycle": result.get("status", "FAILED"), **result})
    return result


def evaluate_family_full(anchor_coords: dict[str, float], evaluator: DiscoveryEvaluator, ledger,
                          run_dir: Path, worker_id: str, tag: str) -> dict[str, Any]:
    """Evaluate all three mass members. Never partial: reports what happened at each."""
    coords = dict(anchor_coords)
    coords["mH_GeV"] = 200.0
    try:
        members_coords = family_members(coords)
    except ValueError as error:
        return {"ok": False, "reason": f"family_derivation_failed:{error}", "anchor_coordinates": coords}
    per_mass: dict[int, dict[str, Any]] = {}
    for member in members_coords:
        mass = int(member["mH_GeV"])
        try:
            check_global(**member)
        except BoundsError as error:
            per_mass[mass] = {"status": "FAILED", "failure_stage": "global_bounds",
                              "failure_reason": str(error)}
            continue
        result = _eval_member(member, evaluator, ledger, run_dir, worker_id, f"{tag}_m{mass}")
        per_mass[mass] = result
    terminated = {m: r for m, r in per_mass.items() if r.get("status") == "TERMINATED"}
    all_valid = len(terminated) == 3 and all(r.get("validity_gate") for r in terminated.values())
    xs = {round(float(r["coordinates"]["X"]), 12) for r in terminated.values()} if terminated else set()
    same_x = len(xs) == 1
    s_val = coords["mA_GeV"] ** 2 - 200.0 ** 2
    ctau = {m: r.get("ctau_mm") for m, r in terminated.items()}
    pf = {m: _pf(r.get("observables") or {}) for m, r in terminated.items()}
    theory = {m: bool(r.get("validity_gate")) for m, r in terminated.items()}
    q_eff = {}
    for m, r in terminated.items():
        p = r.get("parameters") or {}
        mH2, M2, tb = p.get("mH_GeV", 0.0) ** 2, p.get("M2_GeV2"), p.get("tan_beta")
        q_eff[m] = (mH2 - M2) * tb * tb if M2 is not None and tb else None
    ctau200 = ctau.get(200)
    in_window = ctau200 is not None and TARGET[0] <= ctau200 <= TARGET[1]
    ok = all_valid and same_x and s_val > MIN_S and in_window
    pf_values = [pf[m] for m in (150, 200, 250) if pf.get(m) is not None]
    return {
        "ok": ok, "all_valid": all_valid, "same_X": same_x, "S_GeV2": s_val,
        "in_target_window": in_window,
        "anchor_coordinates": coords,
        "X": coords["X"], "tan_beta": coords["tan_beta"], "mA_GeV": coords["mA_GeV"],
        "Q_requested": coords["Q"], "lambda6": coords["X"] / coords["tan_beta"],
        "M2_anchor_GeV2": (terminated.get(200, {}).get("parameters") or {}).get("M2_GeV2"),
        "Q_effective_200": q_eff.get(200), "Q_effective_by_mass": q_eff,
        "ctau150": ctau.get(150), "ctau200": ctau.get(200), "ctau250": ctau.get(250),
        "pf150": pf.get(150), "pf200": pf.get(200), "pf250": pf.get(250),
        "min_pf": min(pf_values) if len(pf_values) == 3 else None,
        "mean_pf": sum(pf_values) / len(pf_values) if len(pf_values) == 3 else None,
        "theory150": theory.get(150), "theory200": theory.get(200), "theory250": theory.get(250),
        "candidate_id_200": (terminated.get(200) or {}).get("candidate_id"),
        "per_mass_results": per_mass,
    }


def _solve_lifetime_family(anchor_coords: dict[str, float], evaluator: DiscoveryEvaluator, ledger,
                            run_dir: Path, worker_id: str, tag: str, max_iter: int = 5) -> dict[str, Any]:
    """Re-solve tan_beta (only) so ctau_200 lands in window; every iterate is family-confirmed."""
    coords = dict(anchor_coords)
    fam = None
    for i in range(max_iter):
        fam = evaluate_family_full(coords, evaluator, ledger, run_dir, worker_id, f"{tag}s{i}")
        ctau200 = fam.get("ctau200")
        if fam.get("in_target_window") and fam.get("all_valid"):
            return fam
        if not ctau200 or ctau200 <= 0:
            return fam
        coords["tan_beta"] = coords["tan_beta"] * math.sqrt(CENTRE / ctau200)
        try:
            check_global(mH_GeV=200.0, mA_GeV=coords["mA_GeV"], tan_beta=coords["tan_beta"],
                        X=coords["X"], Q=coords["Q"])
        except BoundsError:
            return fam
    return fam


def _feasible(fam: dict[str, Any]) -> bool:
    return bool(fam and fam.get("ok"))


def ascend_lineage(seed_coords: dict[str, float], evaluator: DiscoveryEvaluator, ledger,
                    run_dir: Path, envelope: Envelope, worker_id: str, lineage_id: str,
                    max_generations: int = 30, step0: float = 0.10, min_step: float = 0.006,
                    patience: int = 6, allow_branch: bool = True) -> dict[str, Any]:
    """Trust-region ascent on pf200, subject to full-family feasibility at every step.

    Branching: at most ONE extra branch per lineage (deliberately cost-bounded --
    this keeps a "cheap worker" call to a few dozen evaluator calls instead of an
    unbounded tree). A branch is spawned only when two DIFFERENT coordinates both
    give an improving, feasible step in the same generation.
    """
    from .helpers import to_unit_vector, from_unit_vector, COORD_ORDER

    coords = dict(seed_coords)
    coords["mH_GeV"] = 200.0
    seed_fam = _solve_lifetime_family(coords, evaluator, ledger, run_dir, worker_id, f"{lineage_id}_seed")
    if not _feasible(seed_fam):
        return {"ok": False, "lineage_id": lineage_id, "reason": "seed_infeasible", "seed": seed_fam}
    coords = dict(seed_fam["anchor_coordinates"])

    def trace_row(fam, direction):
        return {"pf150": fam["pf150"], "pf200": fam["pf200"], "pf250": fam["pf250"],
                "min_pf": fam["min_pf"], "mean_pf": fam["mean_pf"],
                "ctau150": fam["ctau150"], "ctau200": fam["ctau200"], "ctau250": fam["ctau250"],
                "theory150": fam["theory150"], "theory200": fam["theory200"], "theory250": fam["theory250"],
                "X": fam["X"], "tan_beta": fam["tan_beta"], "lambda6": fam["lambda6"],
                "mA_GeV": fam["mA_GeV"], "S_GeV2": fam["S_GeV2"],
                "Q_requested": fam["Q_requested"], "Q_effective_200": fam["Q_effective_200"],
                "M2_anchor_GeV2": fam["M2_anchor_GeV2"],
                "candidate_id": fam["candidate_id_200"], "coordinates": dict(fam["anchor_coordinates"]),
                "direction": direction}

    def new_state(tag, coords, fam, step):
        return {"branch_id": tag, "coords": dict(coords), "pf200": fam["pf200"], "fam": fam,
                "step": step, "streak": 0, "blocked": set(), "done": False,
                "no_improve_run": 0, "generations": 0,
                "trace": [{"generation": 0, **trace_row(fam, "seed")}]}

    branches = [new_state(f"{lineage_id}", coords, seed_fam, step0)]
    branch_budget = 1 if allow_branch else 0
    generation = 0
    total_evals = 0
    while generation < max_generations and any(not b["done"] for b in branches):
        generation += 1
        for branch in list(branches):
            if branch["done"]:
                continue
            base_unit = to_unit_vector(branch["coords"], envelope)
            candidates = []
            for index, name in enumerate(COORD_ORDER):
                if name not in FREE_COORDS:
                    continue
                for sign in (+1.0, -1.0):
                    key = f"{name}{'+' if sign > 0 else '-'}"
                    if key in branch["blocked"]:
                        continue
                    trial_unit = list(base_unit)
                    trial_unit[index] = min(max(trial_unit[index] + sign * branch["step"], 0.0), 1.0)
                    trial = from_unit_vector(trial_unit, envelope)
                    trial["mH_GeV"] = 200.0
                    trial["tan_beta"] = branch["coords"]["tan_beta"]
                    if trial["mA_GeV"] ** 2 - 200.0 ** 2 <= MIN_S:
                        branch["blocked"].add(key)
                        continue
                    try:
                        check_global(**trial)
                    except BoundsError:
                        branch["blocked"].add(key)
                        continue
                    fam = _solve_lifetime_family(trial, evaluator, ledger, run_dir, worker_id,
                                                 f"{branch['branch_id']}_g{generation}_{key}")
                    total_evals += 1
                    feasible = _feasible(fam)
                    if not feasible:
                        continue
                    improved = fam["pf200"] > branch["pf200"] * (1.0 + 1e-9)
                    candidates.append({"key": key, "coord": name, "fam": fam, "improved": improved})
            improving = [c for c in candidates if c["improved"]]
            if not improving:
                branch["step"] *= 0.5
                branch["no_improve_run"] += 1
                branch["streak"] = 0
                if branch["step"] < min_step or branch["no_improve_run"] >= patience:
                    branch["done"] = True
                continue
            improving.sort(key=lambda c: -c["fam"]["pf200"])
            best = improving[0]
            distinct = [c for c in improving[1:] if c["coord"] != best["coord"]]
            branch["coords"] = dict(best["fam"]["anchor_coordinates"])
            branch["pf200"] = best["fam"]["pf200"]
            branch["fam"] = best["fam"]
            branch["no_improve_run"] = 0
            branch["streak"] += 1
            if branch["streak"] >= 3:
                branch["step"] = min(branch["step"] * 1.4, 0.35)
                branch["streak"] = 0
            branch["generations"] = generation
            branch["trace"].append({"generation": generation, **trace_row(best["fam"], best["key"])})
            if distinct and branch_budget > 0:
                branch_budget -= 1
                second = distinct[0]
                new_id = f"{lineage_id}_branch{2 - branch_budget}"
                nb = new_state(new_id, second["fam"]["anchor_coordinates"], second["fam"], branch["step"])
                nb["trace"][0]["direction"] = f"branch_from_{branch['branch_id']}_g{generation}_{second['key']}"
                branches.append(nb)
    results = []
    for branch in branches:
        results.append({"branch_id": branch["branch_id"], "generations": branch["generations"],
                        "final_pf200": branch["pf200"], "final_min_pf": branch["fam"]["min_pf"],
                        "final_mean_pf": branch["fam"]["mean_pf"], "final_ctau200": branch["fam"]["ctau200"],
                        "final_candidate_id": branch["fam"]["candidate_id_200"],
                        "final_coordinates": branch["coords"], "trace": branch["trace"]})
    best_branch = max(results, key=lambda r: r["final_pf200"])
    return {"ok": True, "lineage_id": lineage_id, "seed_pf200": seed_fam["pf200"],
            "seed_candidate_id": seed_fam["candidate_id_200"],
            "seed_coordinates": seed_fam["anchor_coordinates"],
            "branches": results, "best_branch_id": best_branch["branch_id"],
            "best_pf200": best_branch["final_pf200"],
            "enrichment_factor": (best_branch["final_pf200"] / seed_fam["pf200"]
                                  if seed_fam["pf200"] else None),
            "total_family_evaluations": total_evals}
