"""Deterministic proposal mechanisms. No LLM chooses coordinates here.

Search space is (mH_GeV, mA_GeV, tan_beta, X, Q). tan_beta and X are handled in
log space; Q uses a signed log ("symlog") map because it spans both signs and
many decades; mH and mA are linear.
"""
from __future__ import annotations

import hashlib
import math
import random
from typing import Any, Iterable

import numpy as np

from .bounds import Envelope

LOG_COORDS = ("tan_beta",)
SYMLOG_COORDS = ("Q",)
LINEAR_COORDS = ("mH_GeV", "mA_GeV")
COORD_ORDER = ("mH_GeV", "mA_GeV", "tan_beta", "X", "Q")
Q_LINTHRESH = 1.0
X_FLOOR = 1e-6


# ---------------------------------------------------------------- transforms
def _to_unit(name: str, value: float, envelope: Envelope) -> float:
    low, high = getattr(envelope, name)
    if name in LOG_COORDS:
        low, high = math.log(low), math.log(high)
        value = math.log(value)
    elif name == "X":
        low = math.log(max(low, X_FLOOR)); high = math.log(max(high, X_FLOOR * 10))
        value = math.log(max(value, X_FLOOR))
    elif name in SYMLOG_COORDS:
        low, high = _symlog(low), _symlog(high)
        value = _symlog(value)
    if high == low:
        return 0.0
    return (value - low) / (high - low)


def _from_unit(name: str, unit: float, envelope: Envelope) -> float:
    unit = min(max(unit, 0.0), 1.0)
    low, high = getattr(envelope, name)
    if name in LOG_COORDS:
        return math.exp(math.log(low) + unit * (math.log(high) - math.log(low)))
    if name == "X":
        lo = math.log(max(low, X_FLOOR)); hi = math.log(max(high, X_FLOOR * 10))
        return math.exp(lo + unit * (hi - lo))
    if name in SYMLOG_COORDS:
        lo, hi = _symlog(low), _symlog(high)
        return _symexp(lo + unit * (hi - lo))
    return low + unit * (high - low)


def _symlog(value: float) -> float:
    return math.copysign(math.log1p(abs(value) / Q_LINTHRESH), value)


def _symexp(value: float) -> float:
    return math.copysign(Q_LINTHRESH * math.expm1(abs(value)), value)


def to_unit_vector(coords: dict[str, float], envelope: Envelope) -> np.ndarray:
    return np.array([_to_unit(name, coords[name], envelope) for name in COORD_ORDER], dtype=float)


def from_unit_vector(vector: Iterable[float], envelope: Envelope) -> dict[str, float]:
    return {name: _from_unit(name, float(u), envelope) for name, u in zip(COORD_ORDER, vector)}


def _repair(coords: dict[str, float], envelope: Envelope) -> dict[str, float]:
    """Enforce the structural mA > mH requirement without leaving the envelope."""
    out = dict(coords)
    lo, hi = envelope.mA_GeV
    minimum = max(lo, out["mH_GeV"] + 1.0)
    if minimum > hi:
        out["mH_GeV"] = min(out["mH_GeV"], hi - 1.0)
        minimum = max(lo, out["mH_GeV"] + 1.0)
    out["mA_GeV"] = min(max(out["mA_GeV"], minimum), hi)
    return out


def _proposal(index: int, seed: int, worker_id: str, strategy: str, coords: dict[str, float],
              generation: int = 0, parents: list[str] | None = None, rationale: str = "") -> dict[str, Any]:
    tag = hashlib.sha256(f"{strategy}:{seed}:{worker_id}:{index}".encode()).hexdigest()[:8]
    return {
        "proposal_id": f"{strategy}_{worker_id}_{seed}_{generation:03d}_{index:04d}_{tag}",
        "strategy_id": strategy, "worker_id": worker_id, "generation": generation,
        "parent_ids": parents or [], "random_seed": seed,
        "rationale": rationale or "deterministic proposal helper",
        "coordinates": coords,
    }


# ------------------------------------------------------- Strategy A: sampling
def sobol_candidates(count: int, seed: int, envelope: Envelope, worker_id: str,
                     generation: int = 0) -> list[dict[str, Any]]:
    """Scrambled Sobol' quasi-random coverage of the active envelope."""
    from scipy.stats import qmc
    engine = qmc.Sobol(d=len(COORD_ORDER), scramble=True, seed=seed)
    points = engine.random(count)
    out = []
    for index, row in enumerate(points):
        coords = _repair(from_unit_vector(row, envelope), envelope)
        out.append(_proposal(index, seed, worker_id, "A_sobol_broad", coords, generation,
                             rationale="scrambled Sobol quasi-random envelope coverage"))
    return out


def random_candidates(count: int, seed: int, envelope: Envelope, worker_id: str,
                      generation: int = 0) -> list[dict[str, Any]]:
    """Independent log/symlog-uniform sampling (distinct stream from Sobol)."""
    rng = random.Random(seed)
    out = []
    for index in range(count):
        coords = _repair(from_unit_vector([rng.random() for _ in COORD_ORDER], envelope), envelope)
        out.append(_proposal(index, seed, worker_id, "A_random_broad", coords, generation,
                             rationale="log-uniform independent envelope coverage"))
    return out


# ------------------------------------------------- Strategy C: local mutation
def perturb(parent_coords: dict[str, float], count: int, seed: int, radius: float,
            envelope: Envelope, worker_id: str, parent_id: str,
            generation: int = 0) -> list[dict[str, Any]]:
    """Bounded Gaussian mutation in unit space; radius is the unit-space sigma."""
    rng = np.random.default_rng(seed)
    base = to_unit_vector(parent_coords, envelope)
    out = []
    for index in range(count):
        vector = base + rng.normal(0.0, radius, size=len(COORD_ORDER))
        coords = _repair(from_unit_vector(np.clip(vector, 0.0, 1.0), envelope), envelope)
        out.append(_proposal(index, seed, worker_id, "C_local_mutation", coords, generation,
                             parents=[parent_id],
                             rationale=f"bounded local mutation radius={radius}"))
    return out


# ------------------------------------------- Strategy B: CMA-ES (deterministic)
class CMAES:
    """Compact (mu/mu_w, lambda)-CMA-ES over the unit hypercube.

    Deterministic given (x0, sigma0, seed). State is JSON-serialisable so a
    worker can run one generation per invocation and resume.
    """

    def __init__(self, x0: Iterable[float], sigma0: float, seed: int, popsize: int | None = None):
        self.n = len(COORD_ORDER)
        self.xmean = np.clip(np.array(list(x0), dtype=float), 0.0, 1.0)
        self.sigma = float(sigma0)
        self.seed = int(seed)
        self.rng = np.random.default_rng(seed)
        self.lam = int(popsize or (4 + int(3 * math.log(self.n))))
        self.mu = self.lam // 2
        weights = np.log(self.mu + 0.5) - np.log(np.arange(1, self.mu + 1))
        self.weights = weights / weights.sum()
        self.mueff = 1.0 / np.sum(self.weights ** 2)
        n = self.n
        self.cc = (4 + self.mueff / n) / (n + 4 + 2 * self.mueff / n)
        self.cs = (self.mueff + 2) / (n + self.mueff + 5)
        self.c1 = 2 / ((n + 1.3) ** 2 + self.mueff)
        self.cmu = min(1 - self.c1, 2 * (self.mueff - 2 + 1 / self.mueff) / ((n + 2) ** 2 + self.mueff))
        self.damps = 1 + 2 * max(0.0, math.sqrt((self.mueff - 1) / (n + 1)) - 1) + self.cs
        self.pc = np.zeros(n)
        self.ps = np.zeros(n)
        self.B = np.eye(n)
        self.D = np.ones(n)
        self.C = np.eye(n)
        self.chiN = math.sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n ** 2))
        self.generation = 0

    def ask(self) -> list[np.ndarray]:
        self._population = [
            self.xmean + self.sigma * (self.B @ (self.D * self.rng.standard_normal(self.n)))
            for _ in range(self.lam)
        ]
        return self._population

    def tell(self, solutions: list[np.ndarray], fitnesses: list[float]) -> None:
        order = np.argsort(np.array(fitnesses, dtype=float))
        selected = np.array([solutions[i] for i in order[: self.mu]])
        old_mean = self.xmean.copy()
        self.xmean = self.weights @ selected
        n = self.n
        invsqrtC = self.B @ np.diag(1.0 / self.D) @ self.B.T
        self.ps = (1 - self.cs) * self.ps + math.sqrt(self.cs * (2 - self.cs) * self.mueff) * (
            invsqrtC @ (self.xmean - old_mean) / self.sigma)
        self.generation += 1
        hsig = (np.linalg.norm(self.ps) /
                math.sqrt(1 - (1 - self.cs) ** (2 * self.generation)) / self.chiN) < (1.4 + 2 / (n + 1))
        self.pc = (1 - self.cc) * self.pc + (hsig * math.sqrt(self.cc * (2 - self.cc) * self.mueff)
                                             * (self.xmean - old_mean) / self.sigma)
        artmp = (selected - old_mean) / self.sigma
        self.C = ((1 - self.c1 - self.cmu) * self.C
                  + self.c1 * (np.outer(self.pc, self.pc) + (not hsig) * self.cc * (2 - self.cc) * self.C)
                  + self.cmu * artmp.T @ np.diag(self.weights) @ artmp)
        self.sigma *= math.exp((self.cs / self.damps) * (np.linalg.norm(self.ps) / self.chiN - 1))
        self.sigma = float(min(max(self.sigma, 1e-9), 1.0))
        self.C = np.triu(self.C) + np.triu(self.C, 1).T
        eigenvalues, self.B = np.linalg.eigh(self.C)
        self.D = np.sqrt(np.maximum(eigenvalues, 1e-20))

    def proposals(self, envelope: Envelope, worker_id: str, seed: int) -> list[dict[str, Any]]:
        out = []
        for index, vector in enumerate(self.ask()):
            coords = _repair(from_unit_vector(np.clip(vector, 0.0, 1.0), envelope), envelope)
            out.append(_proposal(index, seed, worker_id, "B_cmaes", coords, self.generation,
                                 rationale=f"CMA-ES gen={self.generation} sigma={self.sigma:.4g}"))
        return out


def deduplicate(proposals: Iterable[dict[str, Any]], seen_candidate_ids: set[str],
                envelope: Envelope) -> tuple[list[dict[str, Any]], int]:
    """Drop proposals whose physical candidate_id was already evaluated."""
    from .evaluator import normalize, DiscoveryError
    seen = set(seen_candidate_ids)
    kept, dropped = [], 0
    for proposal in proposals:
        try:
            normalized = normalize(proposal, envelope)
        except DiscoveryError:
            kept.append(proposal)
            continue
        if normalized["candidate_id"] in seen:
            dropped += 1
            continue
        seen.add(normalized["candidate_id"])
        kept.append(proposal)
    return kept, dropped
