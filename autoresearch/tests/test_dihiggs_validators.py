from __future__ import annotations

import math

from autoresearch.harness.dihiggs_physics_contract import get_default_physics_contract
from autoresearch.harness.dihiggs_validators import (
    compute_triple_ok,
    derive_br_bb,
    validate_point,
    validate_region,
)


def _contract() -> dict[str, object]:
    return get_default_physics_contract()


def _valid_region() -> dict[str, object]:
    return {
        "lambda6": [1e-6, 0.5],
        "tan_beta": [2e3, 1e7],
        "lambda1": [0.2, 11.0],
        "mH": [130.0, 250.0],
    }


def _valid_point() -> dict[str, object]:
    return {
        "lambda6": 1e-4,
        "lambda7": 0.0,
        "tan_beta": 1e4,
        "lambda1": 1.5,
        "mH": 200.0,
        "sin_ba": 1.0,
        "mA": 300.0,
        "mHp": 300.0,
        "yukawa_type": 1,
    }


def test_validate_region_accepts_valid_friendly_region() -> None:
    assert validate_region(_valid_region(), _contract()) == []


def test_validate_region_rejects_lambda6_strict_edges() -> None:
    contract = _contract()

    low = _valid_region()
    low["lambda6"] = [1e-9, 0.5]
    errors_low = validate_region(low, contract)
    assert any("lambda6" in msg for msg in errors_low)

    high = _valid_region()
    high["lambda6"] = [1e-6, 1.0]
    errors_high = validate_region(high, contract)
    assert any("lambda6" in msg for msg in errors_high)


def test_validate_region_rejects_tan_beta_outside_range() -> None:
    region = _valid_region()
    region["tan_beta"] = [999.0, 2000.0]
    errors = validate_region(region, _contract())
    assert any("tan_beta" in msg for msg in errors)


def test_validate_region_rejects_lambda1_outside_range() -> None:
    region = _valid_region()
    region["lambda1"] = [0.09, 5.0]
    errors = validate_region(region, _contract())
    assert any("lambda1" in msg for msg in errors)


def test_validate_region_rejects_mh_ordering_violations() -> None:
    region_low = _valid_region()
    region_low["mH"] = [125.0, 250.0]
    errors_low = validate_region(region_low, _contract())
    assert any("mH" in msg for msg in errors_low)

    region_high = _valid_region()
    region_high["mH"] = [130.0, 300.0]
    errors_high = validate_region(region_high, _contract())
    assert any("mH" in msg for msg in errors_high)


def test_validate_region_rejects_mh_ordering_violations_via_aliases() -> None:
    # The mh < mH < mA hierarchy must be enforced regardless of how the
    # heavy scalar is spelled in the region bounds.
    for alias in ("mH", "m_phi", "mphi"):
        region_low = _valid_region()
        del region_low["mH"]
        region_low[alias] = [125.0, 250.0]
        errors_low = validate_region(region_low, _contract())
        assert any("mH" in msg for msg in errors_low), alias

        region_high = _valid_region()
        del region_high["mH"]
        region_high[alias] = [130.0, 300.0]
        errors_high = validate_region(region_high, _contract())
        assert any("mH" in msg for msg in errors_high), alias


def test_validate_region_accepts_valid_heavy_scalar_aliases() -> None:
    for alias in ("mH", "m_phi", "mphi"):
        region = _valid_region()
        del region["mH"]
        region[alias] = [130.0, 250.0]
        assert validate_region(region, _contract()) == [], alias


def test_validate_region_rejects_lambda7_not_zero_if_supplied() -> None:
    region = _valid_region()
    region["lambda7"] = 1e-3
    errors = validate_region(region, _contract())
    assert any("lambda7" in msg for msg in errors)


def test_validate_region_accepts_nested_bounds_shape() -> None:
    region = {"bounds": _valid_region()}
    assert validate_region(region, _contract()) == []


def test_validate_region_error_messages_include_field_names() -> None:
    region = _valid_region()
    region["lambda6"] = [1e-9, 1.0]
    errors = validate_region(region, _contract())
    joined = " | ".join(errors)
    assert "lambda6" in joined


def test_validate_point_accepts_valid_point() -> None:
    assert validate_point(_valid_point(), _contract()) == []


def test_validate_point_accepts_aliases() -> None:
    point = {
        "lambda_6": 1e-4,
        "lambda_7": 0.0,
        "tanbeta": 1e4,
        "lam1": 0.5,
        "mphi": 150.0,
    }
    assert validate_point(point, _contract()) == []


def test_validate_point_rejects_mH_le_mh() -> None:
    point = _valid_point()
    point["mH"] = 125.0
    errors = validate_point(point, _contract())
    assert any("mH" in msg for msg in errors)


def test_validate_point_rejects_mH_ge_mA() -> None:
    point = _valid_point()
    point["mH"] = 300.0
    errors = validate_point(point, _contract())
    assert any("mH" in msg for msg in errors)


def test_validate_point_rejects_lambda7_not_zero() -> None:
    point = _valid_point()
    point["lambda7"] = 0.1
    errors = validate_point(point, _contract())
    assert any("lambda7" in msg for msg in errors)


def test_validate_point_rejects_fixed_param_mismatch() -> None:
    point = _valid_point()
    point["sin_ba"] = 0.9
    errors = validate_point(point, _contract())
    assert any("sin_ba" in msg for msg in errors)


def test_validate_point_rejects_nan_and_inf_values() -> None:
    point_nan = _valid_point()
    point_nan["lambda1"] = float("nan")
    errors_nan = validate_point(point_nan, _contract())
    assert any("lambda1" in msg for msg in errors_nan)

    point_inf = _valid_point()
    point_inf["lambda6"] = float("inf")
    errors_inf = validate_point(point_inf, _contract())
    assert any("lambda6" in msg for msg in errors_inf)


def test_compute_triple_ok_true_only_when_all_true() -> None:
    assert compute_triple_ok({"positivity_ok": 1, "unitarity_ok": True, "perturbativity_ok": 1.0})


def test_compute_triple_ok_handles_string_flags() -> None:
    assert compute_triple_ok({"positivity_ok": "1", "unitarity_ok": "1.0", "perturbativity_ok": "1"})
    assert compute_triple_ok({"positivity_ok": "1.000000e+00", "unitarity_ok": "1e0", "perturbativity_ok": "1.0E+00"})


def test_compute_triple_ok_false_when_missing_or_false() -> None:
    assert not compute_triple_ok({"positivity_ok": 1, "unitarity_ok": 1})
    assert not compute_triple_ok({"positivity_ok": 1, "unitarity_ok": 0, "perturbativity_ok": 1})


def test_derive_br_bb_returns_existing_valid_br_bb() -> None:
    assert derive_br_bb({"br_bb": 0.42, "width_bb": 1.0, "total_width": 2.0}) == 0.42


def test_derive_br_bb_derives_from_widths() -> None:
    val = derive_br_bb({"width_bb": 3.0, "total_width": 12.0})
    assert val == 0.25


def test_derive_br_bb_returns_none_for_bad_cases() -> None:
    assert derive_br_bb({"width_bb": 1.0, "total_width": 0.0}) is None
    assert derive_br_bb({"width_bb": 1.0}) is None
    assert derive_br_bb({"width_bb": "bad", "total_width": 2.0}) is None
    assert derive_br_bb({"br_bb": float("nan")}) is None
    assert derive_br_bb({"br_bb": float("inf")}) is None


def test_validate_point_rejects_conflicting_aliases() -> None:
    point = _valid_point()
    point["lambda_6"] = 1e-2
    errors = validate_point(point, _contract())
    assert any("conflicting aliases for lambda6" in msg for msg in errors)


def test_validate_point_rejects_conflicting_fixed_lambda7_aliases() -> None:
    point = _valid_point()
    point["lambda_7"] = 0.2
    errors = validate_point(point, _contract())
    assert any("conflicting aliases for lambda7" in msg for msg in errors)


def test_validate_region_rejects_conflicting_fixed_lambda7_aliases() -> None:
    region = _valid_region()
    region["lambda7"] = 0.0
    region["lambda_7"] = 0.1
    errors = validate_region(region, _contract())
    assert any("conflicting aliases for lambda7" in msg for msg in errors)


def test_validate_region_rejects_conflicting_aliases_for_scanned_axes() -> None:
    region = _valid_region()
    region["lambda_6"] = [1e-3, 2e-3]
    errors = validate_region(region, _contract())
    assert any("conflicting aliases for lambda6" in msg for msg in errors)


def test_validate_point_rejects_fixed_mA_and_mHp_mismatch() -> None:
    point = _valid_point()
    point["mA"] = 310.0
    point["mHp"] = 310.0
    errors = validate_point(point, _contract())
    assert any("mA must match contract fixed value" in msg for msg in errors)
    assert any("mHp must match contract fixed value" in msg for msg in errors)


def test_validate_point_rejects_all_strict_boundaries() -> None:
    contract = _contract()
    for field, value in (("lambda6", 1e-9), ("tan_beta", 1e3), ("lambda1", 1e-1), ("mH", 126.0)):
        point = _valid_point()
        point[field] = value
        errors = validate_point(point, contract)
        assert any(field in msg for msg in errors)

    for field, value in (("lambda6", 1.0), ("tan_beta", 1e8), ("lambda1", 12.0), ("mH", 299.0)):
        point = _valid_point()
        point[field] = value
        errors = validate_point(point, contract)
        assert any(field in msg for msg in errors)


def test_validate_region_rejects_all_strict_boundaries() -> None:
    contract = _contract()
    for field, lo, hi in (
        ("lambda6", 1e-9, 0.2),
        ("tan_beta", 1e3, 2e3),
        ("lambda1", 1e-1, 2.0),
        ("mH", 126.0, 200.0),
    ):
        region = _valid_region()
        region[field] = [lo, hi]
        errors = validate_region(region, contract)
        assert any(field in msg for msg in errors)

    for field, lo, hi in (
        ("lambda6", 1e-6, 1.0),
        ("tan_beta", 2e3, 1e8),
        ("lambda1", 0.2, 12.0),
        ("mH", 130.0, 299.0),
    ):
        region = _valid_region()
        region[field] = [lo, hi]
        errors = validate_region(region, contract)
        assert any(field in msg for msg in errors)


def test_validate_region_rejects_mH_equal_to_mh_or_mA() -> None:
    region_lo = _valid_region()
    region_lo["mH"] = [125.0, 250.0]
    errors_lo = validate_region(region_lo, _contract())
    assert any("mH" in msg for msg in errors_lo)

    region_hi = _valid_region()
    region_hi["mH"] = [130.0, 300.0]
    errors_hi = validate_region(region_hi, _contract())
    assert any("mH" in msg for msg in errors_hi)


def test_derive_br_bb_invalid_direct_value_falls_back_to_valid_widths() -> None:
    val = derive_br_bb({"br_bb": "bad", "width_bb": 2.0, "total_width": 8.0})
    assert val == 0.25
