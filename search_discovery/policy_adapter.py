"""Unified proposal-policy interface over the five existing strategies.

None of `helpers.py`, `climb.py`, or `continuation.py` are modified here --
this module only imports and calls them. Two shapes of strategy exist because
two shapes already exist in the frozen code:

  * "batch"          -- propose(context, budget) returns a list of proposals;
                        the daemon core evaluates them (worker pool) and
                        performs every Ledger.append itself, then calls
                        tell(...) so stateful strategies (CMA-ES) can update.
                        Covers explore / refine / optimize.

  * "self_contained" -- run(context, budget, evaluator, ledger) owns its own
                        ask-evaluate-tell loop internally (this is how
                        climb.ascend and continuation.ascend_lineage were
                        already written: each step re-solves and re-evaluates
                        before deciding the next step, so a fixed batch of
                        proposals cannot be handed off in advance). The
                        daemon runs these synchronously in its own process
                        (never in the worker pool), so Ledger.append is still
                        only ever called from one process at a time here --
                        the same single-writer discipline as the batch path,
                        just enforced by "run synchronously" instead of
                        "collect results in the main process".

A strategy never touches the ledger directly in the batch path, and the
self-contained path's ledger writes are the pre-existing frozen behavior of
climb.py/continuation.py, unmodified.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .bounds import Envelope
from .helpers import CMAES, COORD_ORDER, deduplicate as _dedupe_helper
from .helpers import perturb as _perturb
from .helpers import random_candidates as _random_candidates
from .helpers import sobol_candidates as _sobol_candidates
from .helpers import to_unit_vector
from .objective import max_photon_objective as _max_photon_objective
from .objective import objective as _objective


@dataclass
class ProposeContext:
    run_dir: Path
    envelope: Envelope
    worker_id: str
    generation: int
    rng_seed: int
    best_by_family: dict[str, Any] = field(default_factory=dict)


@dataclass
class ProposeResult:
    proposals: list[dict[str, Any]]
    state: dict[str, Any]


def _pick_parent(context: ProposeContext) -> dict[str, Any] | None:
    for family in ("mixed", "photonic"):
        best = context.best_by_family.get(family)
        if best and best.get("coordinates"):
            return {"candidate_id": best.get("candidate_id"), "coordinates": best["coordinates"]}
    return None


class BatchStrategy:
    kind = "batch"
    name = "base_batch"

    def propose(self, context: ProposeContext, budget: int, state: dict[str, Any]) -> ProposeResult:
        raise NotImplementedError

    def tell(self, proposals: list[dict[str, Any]], results: list[dict[str, Any]],
             state: dict[str, Any]) -> dict[str, Any]:
        return state


class ExploreStrategy(BatchStrategy):
    """Strategy A: broad/diversity coverage of the active envelope."""

    def __init__(self, mode: str = "sobol"):
        self.mode = mode
        self.name = "A_sobol_broad" if mode == "sobol" else "A_random_broad"

    def propose(self, context: ProposeContext, budget: int, state: dict[str, Any]) -> ProposeResult:
        if budget <= 0:
            return ProposeResult([], state)
        seed = int(state.get("seed", context.rng_seed))
        maker = _sobol_candidates if self.mode == "sobol" else _random_candidates
        proposals = maker(budget, seed, context.envelope, context.worker_id, context.generation)
        return ProposeResult(proposals, {"seed": seed + 1})


class RefineStrategy(BatchStrategy):
    """Strategy C: local mutation around the best known candidate per family."""

    name = "C_local_mutation"

    def __init__(self, radius: float = 0.05):
        self.default_radius = radius

    def propose(self, context: ProposeContext, budget: int, state: dict[str, Any]) -> ProposeResult:
        if budget <= 0:
            return ProposeResult([], state)
        parent = _pick_parent(context)
        if parent is None:
            return ProposeResult([], state)
        seed = int(state.get("seed", context.rng_seed))
        radius = float(state.get("radius", self.default_radius))
        proposals = _perturb(parent["coordinates"], budget, seed, radius,
                              context.envelope, context.worker_id, parent["candidate_id"],
                              context.generation)
        return ProposeResult(proposals, {"seed": seed + 1, "radius": radius})


class OptimizeStrategy(BatchStrategy):
    """Strategy B: deterministic (mu/mu_w, lambda)-CMA-ES, one generation per cycle.

    State (xmean/sigma/C/B/D/pc/ps/generation) is serialised to plain lists so
    it round-trips through the JSON checkpoint; the pending ask-vectors are
    carried in state between propose() and tell() within the same cycle.
    """

    def __init__(self, target: str = "mixed", objective_mode: str = "target",
                 require_window: bool = False, sigma0: float = 0.3):
        self.target = target
        self.objective_mode = objective_mode
        self.require_window = require_window
        self.sigma0 = sigma0
        self.name = f"B_cmaes_{target}"

    def _restore(self, state: dict[str, Any], context: ProposeContext | None) -> CMAES:
        if "xmean" in state:
            es = CMAES(state["xmean"], state["sigma"], state["seed"], state.get("popsize"))
            es.C = np.array(state["C"]); es.B = np.array(state["B"]); es.D = np.array(state["D"])
            es.pc = np.array(state["pc"]); es.ps = np.array(state["ps"])
            es.generation = int(state.get("generation", 0))
            # CMAES.__init__ always seeds a fresh generator from `seed`; restore
            # the *actual* evolved bit-generator state on top of that so a
            # restored optimizer does not replay the same random sequence it
            # already consumed in a prior cycle (see rng_state in _serialize).
            if "rng_state" in state:
                es.rng.bit_generator.state = state["rng_state"]
            return es
        assert context is not None
        rng = np.random.default_rng(context.rng_seed)
        x0 = rng.random(len(COORD_ORDER))
        return CMAES(x0, self.sigma0, context.rng_seed, None)

    @staticmethod
    def _serialize(es: CMAES) -> dict[str, Any]:
        return {
            "xmean": es.xmean.tolist(), "sigma": float(es.sigma), "seed": es.seed,
            "popsize": es.lam, "C": es.C.tolist(), "B": es.B.tolist(), "D": es.D.tolist(),
            "pc": es.pc.tolist(), "ps": es.ps.tolist(), "generation": es.generation,
            "rng_state": es.rng.bit_generator.state,
        }

    def propose(self, context: ProposeContext, budget: int, state: dict[str, Any]) -> ProposeResult:
        if budget <= 0:
            return ProposeResult([], state)
        es = self._restore(state, context)
        proposals = es.proposals(context.envelope, context.worker_id, es.seed)[:max(budget, 1)]
        vectors = [to_unit_vector({k: p["coordinates"][k] for k in COORD_ORDER}, context.envelope)
                   for p in proposals]
        new_state = self._serialize(es)  # captured *after* es.proposals()/ask() advanced the RNG
        new_state["_pending_vectors"] = {p["proposal_id"]: v.tolist() for p, v in zip(proposals, vectors)}
        return ProposeResult(proposals, new_state)

    def tell(self, proposals: list[dict[str, Any]], results: list[dict[str, Any]],
             state: dict[str, Any]) -> dict[str, Any]:
        """Update CMA-ES from whatever (proposal, result) pairs survived dedup.

        `proposals`/`results` may be shorter than the full asked population --
        dedup can drop some, and a cycle-walltime timeout can leave some
        unresolved (daemon.py). Pair by `proposal_id` (stable across dedup)
        rather than requiring the historical exact-length match, which used
        to silently no-op the optimizer's *learning* step on any partial drop.
        Progress is preserved either way because propose()'s rng_state is
        captured after ask() already advanced it -- the next propose() will
        not repeat the same population even when tell() can't run a real
        CMA-ES update this cycle.
        """
        pending = state.get("_pending_vectors") or {}
        base_state = {k: v for k, v in state.items() if k != "_pending_vectors"}
        if not pending or "xmean" not in state:
            return base_state
        vectors, fitnesses = [], []
        for proposal, result in zip(proposals, results):
            vector = pending.get(proposal.get("proposal_id"))
            if vector is None:
                continue
            score = (_max_photon_objective(result, self.require_window)
                      if self.objective_mode == "max_photon_valid" else _objective(self.target, result))
            vectors.append(np.array(vector))
            fitnesses.append(score["value"])
        es = self._restore(state, None)
        if len(vectors) < es.mu:
            # Not enough surviving pairs for a valid (mu/mu_w) update this
            # cycle -- skip the CMA-ES step but keep the already-advanced RNG
            # state from propose() so the next ask() still yields new points.
            return base_state
        es.tell(vectors, fitnesses)
        return self._serialize(es)


class SelfContainedStrategy:
    """run(context, budget, evaluator, ledger, state) -> report dict.

    The report's "state" key (if present) is persisted back into
    daemon checkpoint policy_state for this strategy name. These strategies
    are deterministic ask-evaluate-tell loops (climb.ascend /
    continuation.ascend_lineage): re-running from the exact same seed
    candidate_id reproduces the exact same trace and therefore contributes
    zero new evidence while re-spending the whole physics budget on
    already-recorded duplicate candidates. Each strategy therefore tracks the
    seed candidate_ids it has already fully ascended from and skips a repeat
    ascent from an unchanged seed -- it only spends budget again once a
    *different* (better) seed becomes available, e.g. because refine/optimize
    found a new best-by-family candidate.
    """
    kind = "self_contained"
    name = "base_self_contained"

    def run(self, context: ProposeContext, budget: int, evaluator, ledger,
            state: dict[str, Any]) -> dict[str, Any]:
        raise NotImplementedError


class ClimbStrategy(SelfContainedStrategy):
    """Strategy: constrained single-point photon-fraction ascent (search_discovery.climb)."""

    name = "climb"

    def run(self, context: ProposeContext, budget: int, evaluator, ledger,
            state: dict[str, Any]) -> dict[str, Any]:
        parent = _pick_parent(context)
        if parent is None:
            return {"ok": False, "reason": "no_seed_candidate_available", "evaluations": 0, "state": state}
        exhausted = set(state.get("exhausted_seed_candidate_ids", []))
        if parent["candidate_id"] in exhausted:
            return {"ok": True, "skipped": True, "reason": "already_ascended_from_this_seed",
                    "evaluations": 0, "state": state}
        from .climb import ascend
        out = ascend(dict(parent["coordinates"]), evaluator, ledger, context.run_dir,
                     context.envelope, context.worker_id, max_steps=max(1, min(budget, 40)))
        exhausted.add(parent["candidate_id"])
        new_state = {"exhausted_seed_candidate_ids": sorted(exhausted)[-64:]}
        return {"ok": bool(out.get("ok")), "evaluations": _climb_eval_count(out),
                "final_photon_fraction": out.get("final_photon_fraction"), "state": new_state}


def _climb_eval_count(out: dict[str, Any]) -> int:
    if not out.get("ok"):
        return 1
    return int(out.get("evaluations", 0)) + 1


class ContinueLineageStrategy(SelfContainedStrategy):
    """Strategy: multi-mass constrained lineage continuation (search_discovery.continuation)."""

    name = "continue_lineage"

    def run(self, context: ProposeContext, budget: int, evaluator, ledger,
            state: dict[str, Any]) -> dict[str, Any]:
        parent = _pick_parent(context)
        if parent is None:
            return {"ok": False, "reason": "no_seed_candidate_available", "evaluations": 0, "state": state}
        exhausted = set(state.get("exhausted_seed_candidate_ids", []))
        if parent["candidate_id"] in exhausted:
            return {"ok": True, "skipped": True, "reason": "already_ascended_from_this_seed",
                    "evaluations": 0, "state": state}
        from .continuation import ascend_lineage
        lineage_id = f"daemon_lineage_g{context.generation}"
        out = ascend_lineage(dict(parent["coordinates"]), evaluator, ledger, context.run_dir,
                             context.envelope, context.worker_id, lineage_id,
                             max_generations=max(1, min(budget, 30)))
        exhausted.add(parent["candidate_id"])
        new_state = {"exhausted_seed_candidate_ids": sorted(exhausted)[-64:]}
        return {"ok": bool(out.get("ok")),
                "evaluations": int(out.get("total_family_evaluations", 0)) + 3,
                "best_pf200": out.get("best_pf200"), "state": new_state}


BATCH_REGISTRY: dict[str, BatchStrategy] = {
    "explore": ExploreStrategy("sobol"),
    "explore_random": ExploreStrategy("random"),
    "refine": RefineStrategy(),
    "optimize_mixed": OptimizeStrategy(target="mixed"),
    "optimize_photonic": OptimizeStrategy(target="photonic"),
}
SELF_CONTAINED_REGISTRY: dict[str, SelfContainedStrategy] = {
    "climb": ClimbStrategy(),
    "continue_lineage": ContinueLineageStrategy(),
}
ALL_STRATEGY_NAMES = tuple(BATCH_REGISTRY) + tuple(SELF_CONTAINED_REGISTRY)


def resolve_allocation(allocation: dict[str, float]) -> list[tuple[str, float]]:
    if not allocation:
        allocation = {"explore": 1.0}
    unknown = sorted(set(allocation) - set(ALL_STRATEGY_NAMES))
    if unknown:
        raise ValueError(f"unknown_policy_strategy:{','.join(unknown)}")
    total = sum(float(w) for w in allocation.values())
    if total <= 0:
        raise ValueError("policy_allocation_must_sum_positive")
    return [(name, float(weight) / total) for name, weight in allocation.items()]


def split_budget(allocation: list[tuple[str, float]], total_budget: int) -> dict[str, int]:
    """Largest-remainder apportionment so the shares sum exactly to total_budget."""
    if total_budget <= 0 or not allocation:
        return {name: 0 for name, _ in allocation}
    raw = {name: weight * total_budget for name, weight in allocation}
    floors = {name: int(value) for name, value in raw.items()}
    remainder = total_budget - sum(floors.values())
    ranked = sorted(raw, key=lambda name: raw[name] - floors[name], reverse=True)
    for name in ranked[:remainder]:
        floors[name] += 1
    return floors
