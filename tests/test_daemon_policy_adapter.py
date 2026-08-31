from __future__ import annotations

import numpy as np
import pytest

from search_discovery import policy_adapter as pa
from search_discovery.envelopes import ENVELOPES


ENVELOPE = ENVELOPES["E1_mixed_exploit"]


def _context(best_by_family=None):
    return pa.ProposeContext(run_dir="unused", envelope=ENVELOPE, worker_id="test",
                             generation=0, rng_seed=42, best_by_family=best_by_family or {})


def test_resolve_allocation_normalizes_weights():
    allocation = pa.resolve_allocation({"explore": 1.0, "refine": 3.0})
    weights = dict(allocation)
    assert weights["explore"] == pytest.approx(0.25)
    assert weights["refine"] == pytest.approx(0.75)


def test_resolve_allocation_rejects_unknown_strategy():
    with pytest.raises(ValueError, match="unknown_policy_strategy"):
        pa.resolve_allocation({"not_a_real_strategy": 1.0})


def test_resolve_allocation_rejects_nonpositive_total():
    with pytest.raises(ValueError, match="must_sum_positive"):
        pa.resolve_allocation({"explore": 0.0})


def test_split_budget_sums_exactly_to_total_via_largest_remainder():
    allocation = pa.resolve_allocation({"explore": 1.0, "refine": 1.0, "optimize_mixed": 1.0})
    for total in range(0, 23):
        shares = pa.split_budget(allocation, total)
        assert sum(shares.values()) == total
        assert all(v >= 0 for v in shares.values())


def test_explore_strategy_proposes_within_envelope_and_advances_seed():
    strategy = pa.ExploreStrategy("sobol")
    result = strategy.propose(_context(), 8, {})
    assert len(result.proposals) == 8
    for proposal in result.proposals:
        c = proposal["coordinates"]
        for name in ("mH_GeV", "mA_GeV", "tan_beta", "X", "Q"):
            low, high = getattr(ENVELOPE, name)
            assert low <= c[name] <= high
    assert result.state["seed"] != 42  # advanced so a resumed cycle doesn't repeat the stream


def test_explore_strategy_zero_budget_proposes_nothing():
    strategy = pa.ExploreStrategy("random")
    result = strategy.propose(_context(), 0, {})
    assert result.proposals == []


def test_refine_strategy_without_a_seed_candidate_proposes_nothing():
    strategy = pa.RefineStrategy()
    result = strategy.propose(_context(best_by_family={}), 8, {})
    assert result.proposals == []


def test_refine_strategy_perturbs_around_best_known_candidate():
    seed_coords = {"mH_GeV": 200.0, "mA_GeV": 400.0, "tan_beta": 5.0e6, "X": 0.1, "Q": 1000.0}
    best = {"mixed": {"candidate_id": "candidate_seed", "coordinates": seed_coords}}
    strategy = pa.RefineStrategy()
    result = strategy.propose(_context(best_by_family=best), 6, {})
    assert len(result.proposals) == 6
    for proposal in result.proposals:
        assert proposal["parent_ids"] == ["candidate_seed"]


def test_optimize_strategy_cma_es_state_round_trips_through_plain_json_types():
    strategy = pa.OptimizeStrategy(target="mixed")
    first = strategy.propose(_context(), 8, {})
    assert len(first.proposals) > 0
    # state must be JSON-serialisable (as the daemon checkpoint requires)
    import json
    json.dumps(first.state)
    # a fake evaluation batch, just enough shape for tell() to compute an objective
    fake_results = [{"status": "TERMINATED", "validity_gate": False, "gates": {"positivity": False}}
                    for _ in first.proposals]
    told_state = strategy.tell(first.proposals, fake_results, first.state)
    assert "_pending_vectors" not in told_state
    json.dumps(told_state)
    # resuming from persisted state must not crash and must keep the CMA-ES generation moving
    second = strategy.propose(_context(), 8, told_state)
    assert second.state["generation"] >= first.state["generation"]


def test_optimize_strategy_tell_is_noop_safe_without_matching_pending_vectors():
    strategy = pa.OptimizeStrategy(target="mixed")
    result = strategy.tell([], [], {"some": "state", "_pending_vectors": [[0.1] * 5]})
    assert "_pending_vectors" not in result


class _StubLedger:
    def append(self, event):
        raise AssertionError("ledger.append must not be called on the skip path")


def test_climb_strategy_returns_false_without_a_seed_candidate():
    strategy = pa.ClimbStrategy()
    out = strategy.run(_context(best_by_family={}), 10, evaluator=None, ledger=_StubLedger(), state={})
    assert out["ok"] is False
    assert out["evaluations"] == 0


def test_climb_strategy_skips_repeat_ascent_from_an_already_exhausted_seed():
    seed_coords = {"mH_GeV": 200.0, "mA_GeV": 400.0, "tan_beta": 5.0e6, "X": 0.1, "Q": 1000.0}
    best = {"mixed": {"candidate_id": "candidate_seed", "coordinates": seed_coords}}
    strategy = pa.ClimbStrategy()
    state = {"exhausted_seed_candidate_ids": ["candidate_seed"]}
    out = strategy.run(_context(best_by_family=best), 10, evaluator=None, ledger=_StubLedger(), state=state)
    assert out["ok"] is True
    assert out["skipped"] is True
    assert out["evaluations"] == 0


def test_continue_lineage_strategy_skips_repeat_ascent_from_an_already_exhausted_seed():
    seed_coords = {"mH_GeV": 200.0, "mA_GeV": 400.0, "tan_beta": 5.0e6, "X": 0.1, "Q": 1000.0}
    best = {"mixed": {"candidate_id": "candidate_seed", "coordinates": seed_coords}}
    strategy = pa.ContinueLineageStrategy()
    state = {"exhausted_seed_candidate_ids": ["candidate_seed"]}
    out = strategy.run(_context(best_by_family=best), 10, evaluator=None, ledger=_StubLedger(), state=state)
    assert out["ok"] is True
    assert out["skipped"] is True
    assert out["evaluations"] == 0
