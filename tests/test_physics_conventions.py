import importlib
import os
import sys
from decimal import Decimal

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAKE_PIPELINE = os.path.join(REPO_ROOT, "mlpython", "lake_pipeline")
CONVENTIONS = os.path.join(REPO_ROOT, "conventions", "physics_conventions.yaml")
HUMAN_CONTRACT = os.path.join(REPO_ROOT, "docs", "PROJECT_PHYSICS_AUTHORITY.md")
HIGH_MASS_CONTRACT = os.path.join(REPO_ROOT, "docs", "HIGH_MASS_H2_CONTRACT.md")
HIGH_MASS_TEMPLATE = os.path.join(
    REPO_ROOT, "docs", "contracts", "high_mass_campaign_template.yaml"
)
HIGH_MASS_SCHEMA = os.path.join(
    REPO_ROOT, "docs", "contracts", "high_mass_point_schema.yaml"
)
CASCADE_CONTRACT = os.path.join(REPO_ROOT, "docs", "contracts", "cascade_contract.yaml")

sys.path.insert(0, LAKE_PIPELINE)
physics_conventions = importlib.import_module("physics_conventions")


def _load_conventions():
    yaml = pytest.importorskip("yaml")
    with open(CONVENTIONS) as fh:
        return yaml.safe_load(fh)


def test_conventions_file_present():
    assert os.path.exists(CONVENTIONS)


def test_pinned_values():
    assert physics_conventions.HBAR_C_GEV_MM == 1.973269804e-13
    assert physics_conventions.C_MM_PER_NS == 299.792458


def test_loaded_matches_pinned_fallback():
    assert physics_conventions.HBAR_C_GEV_MM == physics_conventions._HBAR_C_GEV_MM_PINNED
    assert physics_conventions.C_MM_PER_NS == physics_conventions._C_MM_PER_NS_PINNED


def test_loaded_matches_conventions_yaml():
    conv = _load_conventions()
    constants = conv["constants"]
    assert physics_conventions.HBAR_C_GEV_MM == float(constants["hbar_c_gev_mm"])
    assert physics_conventions.C_MM_PER_NS == float(constants["c_mm_per_ns"])


def test_v3_schema_and_single_authority():
    conv = _load_conventions()
    assert conv["schema_version"] == "physics_conventions_v3"
    assert conv["authority"]["status"] == "active"
    assert conv["authority"]["repository"] == "fbientrigo/dihiggs"
    assert conv["authority"]["path"] == "conventions/physics_conventions.yaml"


def test_active_and_historical_higgs_mass_conventions():
    masses = _load_conventions()["mass_conventions"]
    active = masses["active"]
    assert active["id"] == "mh_125p13_pdg_2026"
    assert active["value_GeV"] == "125.13"
    assert isinstance(active["value_GeV"], str)
    assert Decimal(active["value_GeV"]) == Decimal("125.13")
    assert active["status"] == "active_new_production"
    assert active["authority_kind"] == "external_measurement"
    assert active["external_measurement_claim"] is True
    assert active["external_source"]["central_value_GeV"] == "125.13"
    assert active["external_source"]["uncertainty_GeV"] == "0.11"

    historical = {item["id"]: item for item in masses["historical"]}
    assert set(historical) == {
        "mh_125p20_project_pre_v3",
        "mh_125p13_high_mass_v1",
        "mh_125p09_lambda1_boundary_legacy",
        "mh_125p0_ufo_autoresearch_legacy",
    }
    assert historical["mh_125p20_project_pre_v3"]["value_GeV"] == "125.20"
    assert historical["mh_125p20_project_pre_v3"]["status"] == (
        "superseded_pending_consumer_migration"
    )
    assert all(item["new_production_allowed"] is False for item in historical.values())
    assert all(item["preserve_for_replay"] is True for item in historical.values())


def test_state_identity_and_pdg_map():
    states = _load_conventions()["states"]
    assert states["h"]["pdg_id"] == 25
    assert states["phi"]["pdg_id"] == 35
    assert states["phi"]["accepted_aliases"] == ["H2"]
    assert states["phi"]["forbidden_unqualified_aliases"] == ["H"]
    assert states["A"]["pdg_id"] == 36
    assert states["H_plus"]["pdg_id"] == 37
    assert states["H_plus"]["antiparticle_pdg_id"] == -37


def test_type_i_basis_and_z2_breaking_semantics():
    model = _load_conventions()["model"]
    assert model["family"] == "general_2HDM"
    assert model["cp_conserving"] is True
    assert model["scalar_potential_basis"] == "generic_Phi1_Phi2"
    yukawa = model["yukawa"]
    assert yukawa["assignment"] == "Type-I"
    assert yukawa["fermion_coupling_doublet"] == "Phi2"
    assert yukawa["tan_beta"]["definition"] == "v2/v1"
    assert yukawa["tan_beta"]["basis_dependent"] is True
    assert yukawa["tan_beta"]["meaningful_only_in_fixed_yukawa_basis"] is True
    z2 = model["z2_breaking"]
    assert z2["soft_coefficients"] == ["m12_sq"]
    assert z2["hard_coefficients"] == ["lambda6", "lambda7"]


def test_campaign_restrictions_are_explicit_choices():
    choices = _load_conventions()["campaign_choices"]
    assert choices["must_be_explicit_in_new_manifest"] is True
    assert choices["authority_supplies_no_default"] is True
    assert choices["exact_alignment"]["status"] == "campaign_choice"
    assert choices["lambda7_zero"]["status"] == "campaign_choice"
    assert choices["neutral_charged_degeneracy"]["status"] == "campaign_choice"


def test_mass_squared_parameters_are_distinct():
    params = _load_conventions()["parameters"]
    assert params["M2"]["definition"] == "m12_sq/(sin(beta)*cos(beta))"
    assert params["m12_sq"]["definition"] == "M2*sin(beta)*cos(beta)"
    assert params["m12_sq"]["z2_breaking"] == "soft"
    assert "Phi2^dagger*Phi2" in params["m22_sq"]["definition"]
    for name, others in {
        "M2": {"m12_sq", "m22_sq"},
        "m12_sq": {"M2", "m22_sq"},
        "m22_sq": {"M2", "m12_sq"},
    }.items():
        assert set(params[name]["distinct_from"]) == others


def test_three_trilinear_prescriptions_remain_distinct():
    prescriptions = _load_conventions()["trilinear_prescriptions"]
    assert set(prescriptions) == {
        "physical_2HDM",
        "effective_mass_only",
        "effective_m22_shifted",
    }
    assert prescriptions["physical_2HDM"]["category"] == "model_derived"
    assert prescriptions["physical_2HDM"]["signed_export_status"] == "unresolved"
    assert prescriptions["physical_2HDM"]["sign_sensitive_use"] == "fail_closed"
    assert prescriptions["effective_mass_only"]["formula"] == "8*m_phi^2/v"
    assert prescriptions["effective_m22_shifted"]["formula"] == (
        "8*(m_phi^2-m22_sq)/v"
    )
    assert prescriptions["effective_mass_only"]["not_a_general_2HDM_identity"] is True
    assert prescriptions["effective_m22_shifted"]["not_a_general_2HDM_identity"] is True


def test_stage_ownership_and_migration_fail_closed():
    conv = _load_conventions()
    owners = conv["stage_ownership"]
    assert set(owners) == {
        "dihiggs",
        "dihiggs_ufo",
        "dihiggs_hep_cross",
        "dihiggs_llp_recast",
        "dihiggs_boundary",
    }
    assert "model_points" in owners["dihiggs"]["owns"]
    assert "UFO_model_implementation" in owners["dihiggs_ufo"]["owns"]
    assert "cross_sections" in owners["dihiggs_hep_cross"]["owns"]
    assert "efficiency" in owners["dihiggs_llp_recast"]["owns"]
    assert "statistical_limits" in owners["dihiggs_boundary"]["owns"]

    copy = conv["migration"]["downstream_copy"]
    assert copy["authority"] is False
    assert copy["checksum_algorithm"] == "sha256"
    assert copy["mismatch_action"] == "reject"
    assert set(copy["required_sidecar_fields"]) == {
        "source_repository",
        "source_commit",
        "source_path",
        "source_sha256",
        "schema_version",
        "copied_utc",
    }
    assert conv["failure_policy"]["default"] == "reject"


def test_pending_125p20_migration_inventory_is_actionable():
    migration = _load_conventions()["migration"]
    pending = migration["pending_125p20_uses"]
    assert pending["required_before_new_production"] is True
    occurrences = {(item["repository"], item["path"]): item for item in pending["occurrences"]}
    assert ("fbientrigo/dihiggs_hep_cross", "src/llp_recast/constants.py") in occurrences
    assert ("fbientrigo/dihiggs_hep_cross", "state/CURRENT_CAMPAIGN.yaml") in occurrences
    assert all(item["action"] for item in occurrences.values())
    collision = migration["excluded_numeric_collision"]
    assert collision["repository"] == "fbientrigo/dihiggs_llp_recast"
    assert collision["path"].endswith("geometry.csv")


def test_human_contract_agrees_with_machine_contract():
    assert os.path.exists(HUMAN_CONTRACT)
    with open(HUMAN_CONTRACT) as fh:
        human = fh.read()
    required_literals = [
        "physics_conventions_v3",
        'm_h = "125.13" GeV',
        "mh_125p13_pdg_2026",
        "mh_125p20_project_pre_v3",
        "mh_125p13_high_mass_v1",
        "mh_125p09_lambda1_boundary_legacy",
        "mh_125p0_ufo_autoresearch_legacy",
        "`phi`",
        "`H2`",
        "`m12_sq` is the soft",
        "Nonzero\n`lambda6` or `lambda7` is hard",
        "M2 sin(beta) cos(beta)",
        "8 m_phi^2 / v",
        "8 (m_phi^2 - m22_sq) / v",
        "Signed export is deliberately `unresolved`",
        "source_commit",
        "source_sha256",
        "fail closed",
    ]
    for literal in required_literals:
        assert literal in human


def test_high_mass_contract_and_templates_consume_v3_authority():
    yaml = pytest.importorskip("yaml")
    with open(HIGH_MASS_TEMPLATE) as fh:
        template = yaml.safe_load(fh)
    with open(HIGH_MASS_SCHEMA) as fh:
        schema = yaml.safe_load(fh)
    with open(CASCADE_CONTRACT) as fh:
        cascade = yaml.safe_load(fh)
    with open(HIGH_MASS_CONTRACT) as fh:
        high_mass = fh.read()

    convention = template["frozen_assumptions"]["physics_convention"]
    active = _load_conventions()["mass_conventions"]["active"]
    assert convention["convention_id"] == active["id"]
    assert convention["m_h_GeV_text"] == active["value_GeV"]
    assert convention["schema_version"] == "physics_conventions_v3"
    assert len(convention["source_commit"]) == 40
    assert len(convention["source_sha256"]) == 64
    assert "physics_convention" in template["runner_contract"]["manifest_required_fields"]
    assert set(template["campaign_choices"]["must_be_explicit_per_manifest"]) == {
        "sin_beta_minus_alpha",
        "lambda7",
        "mA_equals_mHp",
    }

    fields = {field["name"]: field for field in schema["fields"]}
    assert fields["physics_convention"]["required"] is True
    assert set(fields["physics_convention"]["required_fields"]) == {
        "convention_id",
        "m_h_GeV_text",
        "schema_version",
        "source_repository",
        "source_commit",
        "source_path",
        "source_sha256",
    }
    assert "Z2-softly-broken" not in high_mass
    assert "nonzero `lambda7` is hard Z2 breaking" in high_mass
    assert "Campaign choice" in high_mass
    assert cascade["physics_convention"] == convention
    assert cascade["constants_GeV"]["m_h"] == float(convention["m_h_GeV_text"])


def test_new_production_rejects_authority_bytes_not_in_source_commit(tmp_path, monkeypatch):
    from dihiggs.app.orchestrator import physics_authority

    altered = tmp_path / "physics_conventions.yaml"
    altered.write_text(
        physics_authority._CONVENTIONS_PATH.read_text(encoding="utf-8").replace(
            'value_GeV: "125.13"', 'value_GeV: "125.14"', 1
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(physics_authority, "_CONVENTIONS_PATH", altered)
    physics_authority._load_contract.cache_clear()
    try:
        with pytest.raises(RuntimeError, match="bytes do not match source_commit"):
            physics_authority.convention_provenance(
                source_commit=physics_authority._authority_commit(),
                requested_mass_gev=125.14,
            )
    finally:
        physics_authority._load_contract.cache_clear()


def test_ctau_helper():
    w = 2e-13
    assert physics_conventions.ctau_mm_from_width_gev(w) == (
        physics_conventions.HBAR_C_GEV_MM / w
    )
