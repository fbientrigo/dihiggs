"""Constrained ascent on photon fraction along the FEASIBLE manifold.

Starts from an already-validated mixed benchmark and takes small steps in
(mA_GeV, X, Q) only. After every step tan_beta is re-solved so that ctau_200
stays inside the target window, then the frozen evaluator confirms the point.
A step is accepted only if ALL campaign conditions still hold:

    theory_ok_v1 at the mH=200 anchor,
    ctau_200 in [500, 1000] mm,
    S = mA^2 - mH^2 > 501 GeV^2  (so the 250 GeV member exists),

and it is kept only if the photon fraction increased. No condition is ever
relaxed to admit a step; when nothing improves the step size halves.

The anchor mass is pinned at 200 GeV by the campaign contract, so mH is not a
free direction here.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from .bounds import Envelope, check_global, BoundsError
from .evaluator import DiscoveryEvaluator
from .families import photon_fraction as _pf_of

TARGET = (500.0, 1000.0)
CENTRE = 700.0
MIN_S = 501.0
FREE = ("mA_GeV", "X", "Q")


def _eval(coords, evaluator, ledger, run_dir, worker_id, tag):
    import hashlib
    pid = f"{tag}_" + hashlib.sha256(repr(sorted(coords.items())).encode()).hexdigest()[:12]
    proposal = {"proposal_id": pid, "strategy_id": "constrained_photon_ascent",
                "worker_id": worker_id, "generation": 0, "parent_ids": [],
                "random_seed": None, "rationale": "feasible-manifold ascent step",
                "coordinates": dict(coords)}
    result = evaluator.evaluate(proposal, run_dir / "attempts" / pid, None)
    ledger.append({"event": "EVALUATION", "lifecycle": result.get("status", "FAILED"), **result})
    return result


def _solve_lifetime(coords, evaluator, ledger, run_dir, worker_id, tag, max_iter=4):
    """Re-solve tan_beta so ctau_200 lands in the window. Evaluator-confirmed."""
    c = dict(coords)
    last = None
    for i in range(max_iter):
        r = _eval(c, evaluator, ledger, run_dir, worker_id, f"{tag}s{i}")
        last = r
        if r.get("status") != "TERMINATED":
            return c, r
        ctau = r.get("ctau_mm")
        if ctau and TARGET[0] <= ctau <= TARGET[1] and r.get("validity_gate"):
            return c, r
        if not ctau or ctau <= 0:
            return c, r
        c["tan_beta"] = c["tan_beta"] * math.sqrt(CENTRE / ctau)
        try:
            check_global(**c)
        except BoundsError:
            return c, r
    return c, last


def _feasible(coords, result):
    if result is None or result.get("status") != "TERMINATED":
        return False
    if not result.get("validity_gate"):
        return False
    ctau = result.get("ctau_mm")
    if not ctau or not (TARGET[0] <= ctau <= TARGET[1]):
        return False
    s = coords["mA_GeV"] ** 2 - coords["mH_GeV"] ** 2
    return s > MIN_S


def ascend(start_coords: dict[str, float], evaluator: DiscoveryEvaluator, ledger,
           run_dir: Path, envelope: Envelope, worker_id: str, max_steps: int = 40,
           step0: float = 0.08, min_step: float = 0.004) -> dict[str, Any]:
    from .helpers import to_unit_vector, from_unit_vector, COORD_ORDER
    coords = dict(start_coords)
    coords["mH_GeV"] = 200.0
    coords, result = _solve_lifetime(coords, evaluator, ledger, run_dir, worker_id, "seed")
    if not _feasible(coords, result):
        return {"ok": False, "reason": "start_point_not_feasible",
                "start": coords,
                "start_valid": bool(result.get("validity_gate")) if result else None,
                "start_ctau": result.get("ctau_mm") if result else None}
    best_pf = _pf_of(result.get("observables") or {}) or 0.0
    trace = [{"step": 0, "photon_fraction": best_pf, "ctau_200_mm": result.get("ctau_mm"),
              "candidate_id": result.get("candidate_id"), "coordinates": dict(coords),
              "accepted_direction": "start"}]
    step = step0
    evaluations = 0
    for iteration in range(1, max_steps + 1):
        base_unit = to_unit_vector(coords, envelope)
        improved = None
        for index, name in enumerate(COORD_ORDER):
            if name not in FREE:
                continue
            for sign in (+1.0, -1.0):
                trial_unit = list(base_unit)
                trial_unit[index] = min(max(trial_unit[index] + sign * step, 0.0), 1.0)
                trial = from_unit_vector(trial_unit, envelope)
                trial["mH_GeV"] = 200.0
                trial["tan_beta"] = coords["tan_beta"]
                if trial["mA_GeV"] ** 2 - 200.0 ** 2 <= MIN_S:
                    continue
                try:
                    check_global(**trial)
                except BoundsError:
                    continue
                trial, r = _solve_lifetime(trial, evaluator, ledger, run_dir, worker_id,
                                           f"i{iteration}_{name}{int(sign)}")
                evaluations += 1
                if not _feasible(trial, r):
                    continue
                pf = _pf_of(r.get("observables") or {}) or 0.0
                if pf > best_pf * (1.0 + 1e-9) and (improved is None or pf > improved["pf"]):
                    improved = {"pf": pf, "coords": dict(trial), "result": r,
                                "direction": f"{name}{'+' if sign > 0 else '-'}"}
        if improved is None:
            step *= 0.5
            if step < min_step:
                break
            continue
        coords = improved["coords"]
        best_pf = improved["pf"]
        trace.append({"step": iteration, "photon_fraction": best_pf,
                      "ctau_200_mm": improved["result"].get("ctau_mm"),
                      "candidate_id": improved["result"].get("candidate_id"),
                      "coordinates": dict(coords), "accepted_direction": improved["direction"],
                      "step_size": step})
    return {"ok": True, "final_coordinates": coords, "final_photon_fraction": best_pf,
            "start_photon_fraction": trace[0]["photon_fraction"],
            "enrichment_factor": best_pf / trace[0]["photon_fraction"] if trace[0]["photon_fraction"] else None,
            "accepted_steps": len(trace) - 1, "evaluations": evaluations,
            "final_candidate_id": trace[-1]["candidate_id"], "trace": trace}
