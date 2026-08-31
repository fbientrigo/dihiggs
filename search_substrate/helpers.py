from __future__ import annotations

import math
import random
from typing import Any, Iterable

from .contract import BOUNDS, normalize_proposal


def _proposal(index: int, seed: int, worker_id: str, values: dict[str, float], strategy: str, parents: list[str] | None = None) -> dict[str, Any]:
    return {
        "proposal_id": f"{strategy}_{seed}_{index:04d}", "strategy_id": strategy, "worker_id": worker_id,
        "generation": 0, "parent_ids": parents or [], "random_seed": seed,
        "rationale": "non-LLM deterministic proposal helper", "parameters": values,
    }


def random_candidates(count: int, seed: int, worker_id: str = "helper") -> list[dict[str, Any]]:
    rng = random.Random(seed)
    result = []
    for index in range(count):
        tb = math.exp(rng.uniform(math.log(30000.0), math.log(3000000.0)))
        x_low = max(200.0, 0.0001 * tb)
        x_high = min(100000.0, 0.05 * tb)
        x = math.exp(rng.uniform(math.log(x_low), math.log(x_high)))
        mH = rng.uniform(150.0, 250.0)
        s = rng.uniform(180000.0, 220000.0)
        mA = math.sqrt(mH * mH + s)
        q = rng.uniform(-50000.0, 50000.0)
        m2 = mH * mH - q / (tb * tb)
        result.append(_proposal(index, seed, worker_id, {
            "mH_GeV": mH, "mA_GeV": mA, "M2_GeV2": min(max(m2, 20000.0), 63000.0),
            "tan_beta": tb, "lambda6": x / tb,
        }, "random_log_uniform"))
    return result


def perturb(parent: dict[str, Any], count: int, seed: int, radius: float, worker_id: str = "helper") -> list[dict[str, Any]]:
    normalized = normalize_proposal(parent)
    base = normalized["parameters"]
    rng = random.Random(seed)
    result = []
    for index in range(count):
        values = dict(base)
        for name in ("mH_GeV", "mA_GeV", "M2_GeV2"):
            low, high = BOUNDS[name]
            values[name] = min(max(values[name] + rng.uniform(-radius, radius) * (high - low), low), high)
        values["tan_beta"] = min(max(base["tan_beta"] * math.exp(rng.uniform(-radius, radius)), BOUNDS["tan_beta"][0]), BOUNDS["tan_beta"][1])
        x = base["tan_beta"] * base["lambda6"] * math.exp(rng.uniform(-radius, radius))
        values["lambda6"] = min(max(x / values["tan_beta"], BOUNDS["lambda6"][0]), BOUNDS["lambda6"][1])
        result.append(_proposal(index, seed, worker_id, values, "local_perturbation", [normalized["candidate_id"]]))
    return result


def deduplicate(candidates: Iterable[dict[str, Any]], evaluated_ids: set[str]) -> list[dict[str, Any]]:
    seen = set(evaluated_ids)
    result = []
    for candidate in candidates:
        normalized = normalize_proposal(candidate)
        if normalized["candidate_id"] not in seen:
            seen.add(normalized["candidate_id"])
            result.append({
                "proposal_id": normalized["proposal_id"], "strategy_id": normalized["strategy_id"],
                "worker_id": normalized["worker_id"], "generation": normalized["generation"],
                "parent_ids": normalized["parent_ids"], "random_seed": normalized["random_seed"],
                "rationale": normalized["rationale"], "parameters": normalized["parameters"],
            })
    return result
