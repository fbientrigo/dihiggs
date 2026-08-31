"""Regression tests for the daemon <-> self-contained-strategy seen_ids sync
gap (bounded fix loop, MAJOR finding): climb/continue_lineage call
Ledger.append internally, bypassing the daemon's own batch dedup path, so the
live `seen_ids` set never learned about candidates they had already
committed -- a later cycle's batch proposal could then silently re-evaluate
the same candidate_id as if it were new (real symptom observed in the
delivered smoke run: duplicate ledger rows while `total_duplicates_dropped`
reported 0). No physics is evaluated in this file.
"""
from __future__ import annotations

from pathlib import Path

from search_substrate.ledger import Ledger

from search_discovery import daemon
from search_discovery import policy_adapter as pa
from search_discovery.archive import DiscoveryArchive
from search_discovery.envelopes import ENVELOPES

ENVELOPE = ENVELOPES["E1_mixed_exploit"]


def test_new_ledger_candidate_ids_since_reads_only_the_new_tail(tmp_path):
    ledger_path = tmp_path / "ledger.jsonl"
    ledger = Ledger(tmp_path)
    ledger.append({"event": "EVALUATION", "candidate_id": "before_1"})
    offset = ledger_path.stat().st_size
    ledger.append({"event": "EVALUATION", "candidate_id": "after_1"})
    ledger.append({"event": "EVALUATION", "candidate_id": "after_2"})
    ledger.append({"event": "FAMILY_VALIDATION", "candidate_id": "after_3_not_evaluation"})

    ids = daemon._new_ledger_candidate_ids_since(ledger_path, offset)
    assert ids == {"after_1", "after_2"}


def test_new_ledger_candidate_ids_since_missing_file_returns_empty(tmp_path):
    assert daemon._new_ledger_candidate_ids_since(tmp_path / "missing.jsonl", 0) == set()


class _FakeSelfContainedStrategy:
    """Mimics climb/continue_lineage's shape without their own exhaustion
    tracking, to isolate testing the daemon's seen_ids merge specifically."""
    kind = "self_contained"
    name = "fake_sc"

    def __init__(self, candidate_id: str):
        self.candidate_id = candidate_id

    def run(self, context, budget, evaluator, ledger, state):
        ledger.append({"event": "EVALUATION", "status": "TERMINATED",
                       "candidate_id": self.candidate_id, "coordinates": {}})
        return {"ok": True, "evaluations": 1, "state": {}}


def test_run_cycle_merges_self_contained_strategy_ledger_writes_into_seen_ids(tmp_path, monkeypatch):
    run_dir = tmp_path / "runs" / "seen_ids_sync_test"
    run_dir.mkdir(parents=True)
    fake_strategy = _FakeSelfContainedStrategy("dup_candidate_1")
    monkeypatch.setitem(pa.SELF_CONTAINED_REGISTRY, "fake_sc", fake_strategy)
    monkeypatch.setattr(pa, "ALL_STRATEGY_NAMES", pa.ALL_STRATEGY_NAMES + ("fake_sc",))

    config = dict(daemon.DEFAULTS)
    config["run_dir"] = str(run_dir)
    config["envelope_id"] = "E1_mixed_exploit"
    config["policy"] = {"seed": 1, "family_validate_per_cycle": 0, "allocation": {"fake_sc": 1.0}}
    config["runtime"] = {"workers": 1, "max_evaluations_per_cycle": 1,
                         "checkpoint_every": 1, "summary_every": 1, "shutdown_grace_s": 5}
    config["budgets"] = {"max_total_evaluations": None, "max_cycle_walltime_s": 30, "max_total_walltime_s": None}

    ledger = Ledger(run_dir)
    archive = DiscoveryArchive(run_dir)
    state = {"cycle": 0, "policy_state": {}, "total_evaluations": 0, "total_proposed": 0,
             "total_duplicates_dropped": 0, "consecutive_no_improvement": 0,
             "consecutive_evaluator_failures": 0, "families_validated": 0}
    seen_ids: set[str] = set()
    stop_flag = {"stop": False, "signal": None}

    daemon._run_cycle(config, run_dir, ENVELOPE, ledger, archive,
                      Path("/nonexistent/evaluator"), pool=None, state=state,
                      seen_ids=seen_ids, stop_flag=stop_flag)

    assert "dup_candidate_1" in seen_ids, (
        "self-contained strategy committed a candidate to the ledger but the "
        "daemon's live dedup set never learned about it"
    )
