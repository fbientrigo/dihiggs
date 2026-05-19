from __future__ import annotations

import copy

from autoresearch.harness.dihiggs_physics_contract import (
    get_default_physics_contract,
    normalize_contract,
    validate_contract,
)


def test_default_contract_validates_cleanly() -> None:
    contract = get_default_physics_contract()
    assert validate_contract(contract) == []


def test_default_contract_contains_friendly_bounds() -> None:
    contract = get_default_physics_contract()
    vs = contract["variable_space"]
    assert isinstance(vs, dict)
    assert vs["lambda6"] == {"min": 1e-9, "max": 1.0}
    assert vs["tan_beta"] == {"min": 1e3, "max": 1e8}
    assert vs["lambda1"] == {"min": 1e-1, "max": 12.0}


def test_default_fixed_params_match_expected() -> None:
    contract = get_default_physics_contract()
    fixed = contract["fixed"]
    assert isinstance(fixed, dict)
    assert fixed["sin_ba"] == 1.0
    assert fixed["lambda7"] == 0.0
    assert fixed["mA"] == 300.0
    assert fixed["mHp"] == 300.0
    assert fixed["yukawa_type"] == 1


def test_normalize_contract_does_not_mutate_input() -> None:
    raw = get_default_physics_contract()
    raw["fixed"] = dict(raw["fixed"])
    raw["fixed"]["mA"] = "300.0"
    raw_copy = copy.deepcopy(raw)

    out = normalize_contract(raw)
    assert raw == raw_copy
    assert isinstance(out["fixed"], dict)
    assert out["fixed"]["mA"] == 300.0


def test_invalid_contract_detects_inverted_ranges() -> None:
    contract = get_default_physics_contract()
    contract["variable_space"] = dict(contract["variable_space"])
    contract["variable_space"]["lambda6"] = {"min": 2.0, "max": 1.0}

    errors = validate_contract(contract)
    assert any("variable_space.lambda6" in msg for msg in errors)


def test_invalid_contract_detects_missing_required_sections() -> None:
    contract = get_default_physics_contract()
    del contract["runtime_guards"]

    errors = validate_contract(contract)
    assert any("runtime_guards" in msg for msg in errors)


def test_normalize_contract_does_not_mutate_nested_structures() -> None:
    raw = get_default_physics_contract()
    raw["variable_space"] = dict(raw["variable_space"])
    raw["variable_space"]["lambda6"] = {"min": "1e-9", "max": "1.0"}
    raw_copy = copy.deepcopy(raw)

    out = normalize_contract(raw)
    assert raw == raw_copy
    assert out["variable_space"]["lambda6"] == {"min": 1e-9, "max": 1.0}


def test_validate_contract_rejects_nan_and_inf_values() -> None:
    contract_nan = get_default_physics_contract()
    contract_nan["variable_space"] = dict(contract_nan["variable_space"])
    contract_nan["variable_space"]["lambda1"] = {"min": float("nan"), "max": 12.0}
    errors_nan = validate_contract(contract_nan)
    assert any("variable_space.lambda1.min" in msg for msg in errors_nan)

    contract_inf = get_default_physics_contract()
    contract_inf["fixed"] = dict(contract_inf["fixed"])
    contract_inf["fixed"]["mA"] = float("inf")
    errors_inf = validate_contract(contract_inf)
    assert any("fixed.mA" in msg for msg in errors_inf)


def test_validate_contract_rejects_invalid_section_types() -> None:
    contract = get_default_physics_contract()
    contract["schema_version"] = 1
    contract["physics_contract"] = []
    contract["runtime_guards"] = "guards"

    errors = validate_contract(contract)
    assert any("schema_version" in msg for msg in errors)
    assert any("physics_contract" in msg for msg in errors)
    assert any("runtime_guards" in msg for msg in errors)
