"""Pins the discovery layer's equivalence to the frozen substrate."""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from search_substrate.contract import normalize_proposal as frozen_normalize
from search_substrate.evaluator import CanonicalEvaluator, HBAR_C_GEV_MM
from search_discovery.bounds import GLOBAL_PHYSICS_BOUNDS, Envelope, BoundsError, derive
from search_discovery.envelopes import ENVELOPES
from search_discovery.evaluator import DiscoveryEvaluator, DiscoveryError, normalize
from search_discovery.families import provisional_family
from search_discovery.helpers import COORD_ORDER, from_unit_vector, to_unit_vector, sobol_candidates
from search_discovery.objective import objective, lifetime_distance

ROOT = Path(__file__).resolve().parents[1]
EXECUTABLE = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"

# a point inside BOTH the frozen v1 contract envelope and the discovery bounds
IN_V1 = {"mH_GeV": 200.0, "mA_GeV": 500.0, "M2_GeV2": 40000.0,
         "tan_beta": 100000.0, "lambda6": 0.005}   # X = 500, inside frozen v1 [200, 1e5]


def _discovery_proposal(params):
    Q = (params["mH_GeV"] ** 2 - params["M2_GeV2"]) * params["tan_beta"] ** 2
    X = params["tan_beta"] * params["lambda6"]
    return {"proposal_id": "equiv_test", "coordinates": {
        "mH_GeV": params["mH_GeV"], "mA_GeV": params["mA_GeV"],
        "tan_beta": params["tan_beta"], "X": X, "Q": Q}}


def test_candidate_identity_matches_frozen_substrate():
    """Same physical point reached either way must carry the same candidate_id."""
    frozen = frozen_normalize({"proposal_id": "p", "parameters": dict(IN_V1)})
    discovery = normalize(_discovery_proposal(IN_V1))
    assert discovery["candidate_id"] == frozen["candidate_id"]
    assert discovery["physical_identity"] == frozen["physical_identity"]
    assert discovery["canonical_physical_json"] == frozen["canonical_physical_json"]
    assert discovery["fixed_model_settings"] == frozen["fixed_model_settings"]


def test_derived_parameters_match_frozen_substrate():
    discovery = normalize(_discovery_proposal(IN_V1))
    for name, value in IN_V1.items():
        assert discovery["parameters"][name] == pytest.approx(value, rel=1e-12, abs=1e-12)
    assert discovery["derived"]["mHp_GeV"] == IN_V1["mA_GeV"]


def test_evaluator_command_matches_frozen_evaluator():
    """argv handed to the physics binary must be byte-identical."""
    frozen = frozen_normalize({"proposal_id": "equiv_test", "parameters": dict(IN_V1)})
    discovery = normalize(_discovery_proposal(IN_V1))
    out = Path("/tmp/equiv.csv")
    frozen_ev = CanonicalEvaluator(ROOT, EXECUTABLE)
    disc_ev = DiscoveryEvaluator(ROOT, EXECUTABLE)
    p = frozen["parameters"]
    expected = [
        str(frozen_ev.executable), "--campaign-id", "llp_benchmark_search_v1",
        "--run-id", frozen["proposal_id"], "--mh", str(frozen["fixed_model_settings"]["m_h_GeV"]),
        "--mH-min", str(p["mH_GeV"]), "--mH-max", str(p["mH_GeV"]), "--n-mH", "1",
        "--mA", str(p["mA_GeV"]), "--mHp", str(frozen["derived"]["mHp_GeV"]),
        "--yukawa-type", str(frozen["fixed_model_settings"]["yukawa_type"]),
        "--sin-ba", str(frozen["fixed_model_settings"]["sin_beta_minus_alpha"]),
        "--tan-beta", str(p["tan_beta"]), "--M2-min", str(p["M2_GeV2"]),
        "--M2-max", str(p["M2_GeV2"]), "--n-M2", "1", "--lambda6", str(p["lambda6"]),
        "--lambda7", str(frozen["fixed_model_settings"]["lambda7"]), "--output", str(out),
    ]
    assert disc_ev.build_command(discovery, out) == expected


def test_frozen_physics_constants_are_imported_not_redefined():
    import search_discovery.evaluator as de
    assert de.HBAR_C_GEV_MM is HBAR_C_GEV_MM
    assert de.HBAR_C_GEV_MM == 1.973269804e-13


def test_envelopes_are_inside_global_physics_bounds():
    for envelope in ENVELOPES.values():
        for name in ("mH_GeV", "mA_GeV", "tan_beta", "X", "Q"):
            low, high = getattr(envelope, name)
            glow, ghigh = GLOBAL_PHYSICS_BOUNDS[name]
            assert glow <= low <= high <= ghigh


def test_global_bounds_can_express_every_frozen_v1_contract_point():
    """The discovery layer must be able to reach any point the frozen contract allows."""
    from search_substrate.contract import BOUNDS as V1
    assert GLOBAL_PHYSICS_BOUNDS["X"][1] >= 1.0e5          # frozen v1 X ceiling
    assert GLOBAL_PHYSICS_BOUNDS["mH_GeV"] == V1["mH_GeV"]
    assert GLOBAL_PHYSICS_BOUNDS["tan_beta"][0] <= V1["tan_beta"][0]
    assert GLOBAL_PHYSICS_BOUNDS["tan_beta"][1] >= V1["tan_beta"][1]
    assert GLOBAL_PHYSICS_BOUNDS["mA_GeV"][1] >= V1["mA_GeV"][1]


def test_envelope_cannot_escape_global_bounds():
    with pytest.raises(BoundsError):
        Envelope("bad", "r", (150.0, 250.0), (155.0, 900.0),
                 (1.0e4, 1.0e12), (0.0, 0.3), (-1e7, 1e7))


def test_forbidden_coordinates_are_rejected():
    for bad in ("m_h_GeV", "lambda7", "sin_beta_minus_alpha", "family", "score", "M2_GeV2", "lambda6"):
        proposal = _discovery_proposal(IN_V1)
        proposal["coordinates"][bad] = 1.0
        with pytest.raises(DiscoveryError, match="forbidden_coordinate"):
            normalize(proposal)


def test_unit_transform_roundtrip():
    envelope = ENVELOPES["E0"]
    coords = {"mH_GeV": 200.0, "mA_GeV": 500.0, "tan_beta": 4.54e6, "X": 0.10, "Q": 1234.0}
    back = from_unit_vector(to_unit_vector(coords, envelope), envelope)
    for name in COORD_ORDER:
        assert back[name] == pytest.approx(coords[name], rel=1e-9)


def test_sobol_proposals_stay_in_envelope_and_respect_hierarchy():
    envelope = ENVELOPES["E0"]
    for proposal in sobol_candidates(64, 7, envelope, "test"):
        c = proposal["coordinates"]
        for name in COORD_ORDER:
            low, high = getattr(envelope, name)
            assert low <= c[name] <= high
        assert c["mA_GeV"] > c["mH_GeV"]
        normalize(proposal, envelope)


def test_provisional_family_is_never_authoritative():
    record = provisional_family({"br_gammagamma": 0.68, "br_Zgamma": 0.31, "br_bb": 0.0,
                                 "br_gg": 0.0, "br_tautau": 0.0, "br_cc": 0.0}, True)
    assert record["family"] == "photonic"
    assert record["provisional"] is True and record["authoritative"] is False
    assert record["frozen_contract_status"] == "BLOCKED_SCIENTIFIC_DECISION"


def test_objective_treats_invalidity_as_hard_constraint():
    invalid = {"status": "TERMINATED", "validity_gate": False,
               "gates": {"positivity": False, "unitarity": True}, "ctau_mm": 700.0}
    valid = {"status": "TERMINATED", "validity_gate": True, "gates": {},
             "ctau_mm": 700.0,
             "provisional_family": {"photon_fraction": 0.0, "fermionic_fraction": 0.97}}
    assert objective("mixed", invalid)["feasible"] is False
    assert objective("mixed", invalid)["value"] > objective("mixed", valid)["value"]
    assert objective("mixed", valid)["value"] == 0.0


def test_lifetime_distance_zero_inside_target_interval():
    assert lifetime_distance(500.0) == 0.0
    assert lifetime_distance(1000.0) == 0.0
    assert lifetime_distance(700.0) == 0.0
    assert lifetime_distance(350.0) == pytest.approx(1.0, rel=1e-9)
    assert lifetime_distance(None) >= 100.0
