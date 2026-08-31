from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

import pytest

from search_substrate.archive import Archive
from search_substrate.classification import ClassificationRules, classify_observables
from search_substrate.contract import ContractError, normalize_proposal
from search_substrate.evaluator import CanonicalEvaluator, _validity
from search_substrate.helpers import random_candidates
from search_substrate.ledger import Ledger
from search_substrate.score import score_components


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"


def proposal(**updates):
    values = {"mH_GeV": 200.0, "mA_GeV": 469.041575982343, "M2_GeV2": 39999.9999995, "tan_beta": 300000.0, "lambda6": 0.001}
    values.update(updates)
    return {"proposal_id": "p0", "strategy_id": "test", "worker_id": "test", "generation": 0, "parent_ids": [], "random_seed": 7, "parameters": values}


def test_contract_aliases_units_bounds_and_hashes():
    a = normalize_proposal(proposal())
    b = normalize_proposal({**proposal(), "parameters": {"mH": 200.0, "mA": 469.041575982343, "M2": 39999.9999995, "tb": 300000.0, "l6": 0.001}})
    assert a["candidate_id"] == b["candidate_id"]
    assert a["derived"]["X"] == 300.0
    with pytest.raises(ContractError, match="unknown_or_forbidden"):
        normalize_proposal({**proposal(), "parameters": {**proposal()["parameters"], "mH_mm": 200.0}})
    with pytest.raises(ContractError, match="unit_or_type"):
        normalize_proposal({**proposal(), "parameters": {**proposal()["parameters"], "mH_GeV": "200 GeV"}})
    with pytest.raises(ContractError, match="out_of_range:X"):
        normalize_proposal({**proposal(), "parameters": {**proposal()["parameters"], "lambda6": 0.05, "tan_beta": 3000000.0}})


def test_classifier_interface_has_synthetic_fixtures_but_contract_is_blocked():
    rules = ClassificationRules(
        mixed={"br_gammagamma": (0.20, 0.40), "br_bb": (0.20, 0.35)},
        photonic={"br_gammagamma": (0.70, 1.0), "br_bb": (0.0, 0.05)},
    )
    assert classify_observables({"br_gammagamma": 0.30, "br_bb": 0.25}, rules)["family"] == "mixed"
    assert classify_observables({"br_gammagamma": 0.80, "br_bb": 0.02}, rules)["family"] == "photonic"
    assert classify_observables({"br_gammagamma": 0.50, "br_bb": 0.10}, rules)["status"] == "UNCLASSIFIED"
    assert classify_observables({"br_gammagamma": 0.30, "br_bb": 0.25})["status"] == "BLOCKED_SCIENTIFIC_DECISION"
    contract = json.loads((ROOT / "docs/contracts/llp_benchmark_search_v1.json").read_text())
    assert contract["classification"]["status"] == "BLOCKED_SCIENTIFIC_DECISION"
    assert contract["classification"]["rules"] is None


def test_hard_invalidity_cannot_be_compensated():
    components = score_components(valid=False, classification={"status": "CLASSIFIED"}, ctau_mm=700.0, cross_mass_consistency=1.0, diversity=1.0)
    assert components["validity_gate"] == 0.0
    assert components["total_score"] == 0.0


def _family_record(family_id: str, family: str, score: float, offset: float = 0.0, mass: float = 200.0):
    anchor = normalize_proposal(proposal(mH_GeV=mass, mA_GeV=469.041575982343 + offset, M2_GeV2=39999.9999995 + offset * 400.0, lambda6=0.005))
    return {"family_id": family_id, "family": family, "status": "PROMOTABLE", "validity_gate": True, "total_score": score, "anchor": anchor, "score_components": {"total_score": score}}


def test_archive_top5_duplicate_replacement_and_separate_families(tmp_path):
    archive = Archive(tmp_path)
    for index in range(6):
        archive.consider(_family_record(f"m{index}", "mixed", 0.1 + index * 0.1, offset=index * 10.0, mass=150.0 + index * 20.0))
    assert len(archive.load()["mixed"]) == 5
    assert "m0" not in {item["family_id"] for item in archive.load()["mixed"]}
    assert archive.consider(_family_record("m5", "mixed", 0.99, offset=50.0, mass=250.0))["promoted"]
    assert archive.consider(_family_record("m5", "mixed", 0.8, offset=50.0, mass=250.0))["reason"] == "existing_family_is_better_or_equal"
    assert archive.consider(_family_record("p0", "photonic", 0.8))["promoted"]
    assert len(archive.load()["photonic"]) == 1


def test_ledger_is_append_only_and_replay_counts_are_stable(tmp_path):
    first = Ledger(tmp_path / "one")
    first.append({"event": "PROPOSAL", "lifecycle": "PROPOSED", "candidate_id": "c"})
    first.append({"event": "EVALUATION_STARTED", "lifecycle": "EVALUATING", "candidate_id": "c"})
    first.append({"event": "EVALUATION_TERMINATED", "lifecycle": "TERMINATED", "candidate_id": "c"})
    events = list(first.events())
    second = Ledger(tmp_path / "two")
    for event in events:
        second.append({key: value for key, value in event.items() if key not in {"event_id", "timestamp_unix"}})
    assert first.summary()["event_counts"] == second.summary()["event_counts"]
    assert len(list(first.events())) == 3


def test_deterministic_helpers_and_mutation_guards():
    left = random_candidates(4, 123)
    right = random_candidates(4, 123)
    assert [normalize_proposal(item)["candidate_id"] for item in left] == [normalize_proposal(item)["candidate_id"] for item in right]
    with pytest.raises(ContractError):
        normalize_proposal({**proposal(), "archive": {}})
    with pytest.raises(ContractError):
        normalize_proposal({**proposal(), "parameters": {**proposal()["parameters"], "family": "mixed"}})


@pytest.mark.skipif(not BINARY.is_file(), reason="canonical evaluator not built")
def test_canonical_wrapper_uses_v2_and_corrected_yukawa_route(tmp_path):
    result = CanonicalEvaluator(ROOT, BINARY).evaluate(proposal(), tmp_path / "attempt")
    row = result["canonical_evaluator_row"]
    assert result["status"] == "TERMINATED"
    assert row["schema_version"] == "dihiggs.point.v2"
    assert row["producer"] == "DihiggsPointV2Evaluator"
    assert row["yukawa_type"] == "1"
    assert row["yukawa_type_installed"] == "1.00000000000000000e+00"
    assert row["evaluator_api"].startswith("THDM::set_param_phys")
    assert result["ctau_consistency_ok"] is True
    assert result["classification"]["status"] == "BLOCKED_SCIENTIFIC_DECISION"


def test_theory_and_width_failures_remain_data():
    row = {
        "construction_ok": "1", "numerical_ok": "1.00000000000000000e+00",
        "positivity_reported_ok": "0.00000000000000000e+00", "unitarity_ok": "1.00000000000000000e+00",
        "perturbativity_ok": "0.00000000000000000e+00", "theory_ok_v1": "0.00000000000000000e+00",
        "width_ok": "0.00000000000000000e+00", "total_width_GeV": "nan",
    }
    valid, gates = _validity(row, False)
    assert valid is False
    assert gates["theory"] is False and gates["width"] is False and gates["lifetime"] is False
