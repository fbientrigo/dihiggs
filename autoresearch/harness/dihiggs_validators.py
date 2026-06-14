from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any, cast


ALIASES: dict[str, tuple[str, ...]] = {
    "lambda6": ("lambda6", "lambda_6"),
    "lambda7": ("lambda7", "lambda_7"),
    "tan_beta": ("tan_beta", "tanbeta"),
    "lambda1": ("lambda1", "lam1"),
    "mH": ("mH", "m_phi", "mphi"),
}


def _get_first(mapping: Mapping[str, object], keys: tuple[str, ...]) -> object | None:
    for key in keys:
        if key in mapping:
            return mapping[key]
    return None


def _to_finite_float(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float, str)):
        try:
            out = float(value)
        except (TypeError, ValueError):
            return None
        if math.isfinite(out):
            return out
    return None


def _to_int_like(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if not math.isfinite(value) or not value.is_integer():
            return None
        return int(value)
    if isinstance(value, str):
        try:
            parsed = float(value)
        except (TypeError, ValueError):
            return None
        if not math.isfinite(parsed) or not parsed.is_integer():
            return None
        return int(parsed)
    return None


def _find_conflicting_aliases(mapping: Mapping[str, object], canonical: str, keys: tuple[str, ...]) -> str | None:
    present = {k: mapping[k] for k in keys if k in mapping}
    if len(present) <= 1:
        return None
    parsed_values: list[float] = []
    for value in present.values():
        parsed = _to_finite_float(value)
        if parsed is None:
            return f"conflicting aliases for {canonical}: {', '.join(present.keys())}"
        parsed_values.append(parsed)
    first = parsed_values[0]
    if any(v != first for v in parsed_values[1:]):
        return f"conflicting aliases for {canonical}: {', '.join(present.keys())}"
    return None


def _strict_for(contract: Mapping[str, object], field: str) -> bool:
    pc = contract.get("physics_contract")
    if not isinstance(pc, Mapping):
        return False
    strict = pc.get("strict_bounds")
    if not isinstance(strict, Mapping):
        return False
    val = strict.get(field)
    return bool(val) if isinstance(val, bool) else False


def _get_contract_range(contract: Mapping[str, object], field: str) -> tuple[float, float] | None:
    vs = contract.get("variable_space")
    if not isinstance(vs, Mapping):
        return None
    axis = vs.get(field)
    if not isinstance(axis, Mapping):
        return None
    lo = _to_finite_float(axis.get("min"))
    hi = _to_finite_float(axis.get("max"))
    if lo is None or hi is None:
        return None
    return lo, hi


def _get_fixed(contract: Mapping[str, object], field: str) -> object | None:
    fixed = contract.get("fixed")
    if not isinstance(fixed, Mapping):
        return None
    return fixed.get(field)


def validate_region(region: Mapping[str, object], contract: Mapping[str, object]) -> list[str]:
    errors: list[str] = []
    bounds_obj: object = region.get("bounds", region)
    if not isinstance(bounds_obj, Mapping):
        return ["bounds must be a mapping"]
    bounds = cast(Mapping[str, object], bounds_obj)

    for canon in ("lambda6", "tan_beta", "lambda1", "mH", "lambda7"):
        conflict = _find_conflicting_aliases(bounds, canon, ALIASES[canon])
        if conflict is not None:
            errors.append(conflict)

    for field in ("lambda6", "tan_beta", "lambda1", "mH"):
        raw_pair = _get_first(bounds, ALIASES[field])
        if not isinstance(raw_pair, (list, tuple)) or len(raw_pair) != 2:
            errors.append(f"{field} must be [lo, hi]")
            continue

        lo = _to_finite_float(raw_pair[0])
        hi = _to_finite_float(raw_pair[1])
        if lo is None or hi is None:
            errors.append(f"{field} bounds must be finite numeric values")
            continue
        if not (lo < hi):
            errors.append(f"{field} must satisfy lo < hi")
            continue

        crange = _get_contract_range(contract, field)
        if crange is None:
            errors.append(f"contract missing variable_space.{field}")
            continue
        c_lo, c_hi = crange
        strict = _strict_for(contract, field)

        if strict:
            if lo <= c_lo:
                errors.append(f"{field} lower bound must be > contract min")
            if hi >= c_hi:
                errors.append(f"{field} upper bound must be < contract max")
        else:
            if lo < c_lo:
                errors.append(f"{field} lower bound below contract min")
            if hi > c_hi:
                errors.append(f"{field} upper bound above contract max")

    fixed_mh = _to_finite_float(_get_fixed(contract, "mh"))
    fixed_mA = _to_finite_float(_get_fixed(contract, "mA"))
    # Honor the same heavy-scalar aliases (mH / m_phi / mphi) used by the
    # ordinary bounds validation above so the mh < mH < mA hierarchy is
    # enforced regardless of how the heavy scalar is spelled.
    raw_mh = _get_first(bounds, ALIASES["mH"])
    if isinstance(raw_mh, (list, tuple)) and len(raw_mh) == 2 and fixed_mh is not None and fixed_mA is not None:
        m_lo = _to_finite_float(raw_mh[0])
        m_hi = _to_finite_float(raw_mh[1])
        if m_lo is not None and not (fixed_mh < m_lo):
            errors.append("mH must satisfy mh < mH lower bound")
        if m_hi is not None and not (m_hi < fixed_mA):
            errors.append("mH must satisfy mH upper bound < mA")

    lam7_raw = _get_first(bounds, ALIASES["lambda7"])
    if lam7_raw is not None:
        lam7 = _to_finite_float(lam7_raw)
        if lam7 is None or lam7 != 0.0:
            errors.append("lambda7 must be 0.0")

    for fixed_field in ("sin_ba", "mA", "mHp", "yukawa_type"):
        if fixed_field in bounds:
            expected = _get_fixed(contract, fixed_field)
            got = bounds[fixed_field]
            if fixed_field == "yukawa_type":
                exp_i = _to_int_like(expected)
                got_i = _to_int_like(got)
                if exp_i is None or got_i is None or exp_i != got_i:
                    errors.append(f"{fixed_field} must match contract fixed value")
            else:
                exp_f = _to_finite_float(expected)
                got_f = _to_finite_float(got)
                if exp_f is None or got_f is None or got_f != exp_f:
                    errors.append(f"{fixed_field} must match contract fixed value")

    return errors


def validate_point(point: Mapping[str, object], contract: Mapping[str, object]) -> list[str]:
    errors: list[str] = []

    for canon in ("lambda6", "tan_beta", "lambda1", "mH", "lambda7"):
        conflict = _find_conflicting_aliases(point, canon, ALIASES[canon])
        if conflict is not None:
            errors.append(conflict)

    values: dict[str, float] = {}
    for canon in ("lambda6", "tan_beta", "lambda1", "mH"):
        raw = _get_first(point, ALIASES[canon])
        if raw is None:
            errors.append(f"{canon} is required")
            continue
        num = _to_finite_float(raw)
        if num is None:
            errors.append(f"{canon} must be finite numeric")
            continue
        values[canon] = num

    for canon in ("lambda6", "tan_beta", "lambda1", "mH"):
        if canon not in values:
            continue
        crange = _get_contract_range(contract, canon)
        if crange is None:
            errors.append(f"contract missing variable_space.{canon}")
            continue
        c_lo, c_hi = crange
        val = values[canon]
        strict = _strict_for(contract, canon)
        if strict:
            if not (c_lo < val < c_hi):
                errors.append(f"{canon} must be inside strict range ({c_lo}, {c_hi})")
        else:
            if not (c_lo <= val <= c_hi):
                errors.append(f"{canon} must be inside range [{c_lo}, {c_hi}]")

    lam7_raw = _get_first(point, ALIASES["lambda7"])
    if lam7_raw is not None:
        lam7 = _to_finite_float(lam7_raw)
        if lam7 is None or lam7 != 0.0:
            errors.append("lambda7 must be 0.0")

    fixed_mh = _to_finite_float(_get_fixed(contract, "mh"))
    fixed_mA = _to_finite_float(_get_fixed(contract, "mA"))
    if "mH" in values and fixed_mh is not None and fixed_mA is not None:
        if not (fixed_mh < values["mH"]):
            errors.append("mH must satisfy mh < mH")
        if not (values["mH"] < fixed_mA):
            errors.append("mH must satisfy mH < mA")

    for fixed_field in ("sin_ba", "mA", "mHp", "yukawa_type"):
        if fixed_field not in point:
            continue
        expected = _get_fixed(contract, fixed_field)
        got = point[fixed_field]
        if fixed_field == "yukawa_type":
            exp_i = _to_int_like(expected)
            got_i = _to_int_like(got)
            if exp_i is None or got_i is None or exp_i != got_i:
                errors.append(f"{fixed_field} must match contract fixed value")
        else:
            exp_f = _to_finite_float(expected)
            got_f = _to_finite_float(got)
            if exp_f is None or got_f is None or got_f != exp_f:
                errors.append(f"{fixed_field} must match contract fixed value")

    return errors


def _flag_is_true(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return float(value) == 1.0
    if isinstance(value, str):
        text = value.strip()
        if text == "":
            return False
        try:
            return float(text) == 1.0
        except Exception:
            return False
    return False


def compute_triple_ok(point: Mapping[str, object]) -> bool:
    for field in ("positivity_ok", "unitarity_ok", "perturbativity_ok"):
        if field not in point:
            return False
        if not _flag_is_true(point[field]):
            return False
    return True


def derive_br_bb(point: Mapping[str, object]) -> float | None:
    br_bb_raw = point.get("br_bb")
    br_bb = _to_finite_float(br_bb_raw)
    if br_bb is not None:
        return br_bb

    width_bb = _to_finite_float(point.get("width_bb"))
    total_width = _to_finite_float(point.get("total_width"))
    if width_bb is None or total_width is None:
        return None
    if total_width <= 0:
        return None
    return width_bb / total_width
