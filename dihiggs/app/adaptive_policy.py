from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from typing import TypedDict


class Lam1Bin(TypedDict):
    id: str
    index: int
    lam1_lo: float
    lam1_hi: float


def make_lam1_bins(lam1_min: float, lam1_max: float, n_bins: int) -> list[Lam1Bin]:
    if n_bins <= 0:
        raise ValueError("n_bins must be > 0")
    if not (lam1_max > lam1_min):
        raise ValueError("lam1_max must be > lam1_min")

    width = lam1_max - lam1_min
    edges = [lam1_min + width * (i / n_bins) for i in range(n_bins + 1)]
    edges[-1] = lam1_max

    bins: list[Lam1Bin] = []
    for i in range(n_bins):
        bins.append(
            {
                "id": f"lam1_bin_{i:03d}",
                "index": i,
                "lam1_lo": float(edges[i]),
                "lam1_hi": float(edges[i + 1]),
            }
        )
    return bins


def beta_binomial_posterior_mean(
    successes: int, trials: int, alpha: float, beta: float
) -> float:
    if trials < 0:
        raise ValueError("trials must be >= 0")
    if successes < 0:
        raise ValueError("successes must be >= 0")
    if alpha < 0.0:
        raise ValueError("alpha must be >= 0")
    if beta < 0.0:
        raise ValueError("beta must be >= 0")

    k = min(successes, trials)
    denom = alpha + beta + trials
    if denom <= 0.0:
        return 0.0
    return (alpha + k) / denom


def allocate_budget(
    bins: Iterable[Mapping[str, object]],
    successes_by_bin: Mapping[str, int],
    trials_by_bin: Mapping[str, int],
    total_budget: int,
    floor_points: int,
    alpha: float = 1.0,
    beta: float = 1.0,
) -> dict[str, int]:
    bins_list = list(bins)
    if not bins_list:
        if total_budget == 0:
            return {}
        raise ValueError("bins must be non-empty")

    if floor_points < 0:
        raise ValueError("floor_points must be >= 0")
    if total_budget < 0:
        raise ValueError("total_budget must be >= 0")

    bin_ids: list[str] = []
    for b in bins_list:
        bid = b.get("id")
        if not isinstance(bid, str) or not bid:
            raise ValueError("each bin must have a non-empty string 'id'")
        bin_ids.append(bid)

    if len(set(bin_ids)) != len(bin_ids):
        raise ValueError("bin ids must be unique")

    n_bins = len(bin_ids)
    floor_total = n_bins * floor_points
    if total_budget < floor_total:
        raise ValueError("total_budget is smaller than n_bins*floor_points")

    budgets: dict[str, int] = {bid: floor_points for bid in bin_ids}
    remaining = total_budget - floor_total
    if remaining == 0:
        return budgets

    scores: list[float] = []
    for bid in bin_ids:
        s_raw = successes_by_bin.get(bid, 0)
        t_raw = trials_by_bin.get(bid, 0)
        try:
            s = int(s_raw)
            t = int(t_raw)
        except (TypeError, ValueError):
            s = 0
            t = 0

        score = beta_binomial_posterior_mean(successes=s, trials=t, alpha=alpha, beta=beta)
        scores.append(float(score))

    total_score = sum(scores)
    if total_score <= 0.0:
        scores = [1.0 for _ in scores]
        total_score = float(n_bins)

    quotas = [remaining * (s / total_score) for s in scores]
    floor_extras = [int(math.floor(q)) for q in quotas]
    used = sum(floor_extras)
    leftovers = remaining - used

    extras = list(floor_extras)
    if leftovers:
        remainders = [q - fe for q, fe in zip(quotas, floor_extras)]
        order = sorted(
            range(n_bins),
            key=lambda i: (-remainders[i], bin_ids[i]),
        )
        for i in order[:leftovers]:
            extras[i] += 1

    for bid, extra in zip(bin_ids, extras):
        budgets[bid] += int(extra)

    return budgets
