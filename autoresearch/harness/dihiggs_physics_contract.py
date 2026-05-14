from __future__ import annotations

import copy
import math
from collections.abc import Mapping
from typing import Any, cast


def get_default_physics_contract() -> dict[str, object]:
    return {
        "schema_version": "1.0",
        "physics_contract": {
            "model": "2HDM",
            "strict_bounds": {
                "lambda6": True,
                "tan_beta": True,
                "lambda1": True,
                "mH": True,
            },
        },
        "variable_space": {
            "lambda6": {"min": 1e-9, "max": 1.0},
            "tan_beta": {"min": 1e3, "max": 1e8},
            "lambda1": {"min": 1e-1, "max": 12.0},
            # Operational margin: keep mH.min above fixed mh=125.0 by ~1 GeV.
            # Physics intent is strict mh < mH < mA; we encode 126.0 as
            # conservative scanning floor (instead of 125.1) to avoid edge noise.
            "mH": {"min": 126.0, "max": 299.0},
        },
        "fixed": {
            "mh": 125.0,
            "sin_ba": 1.0,
            "mA": 300.0,
            "mHp": 300.0,
            "yukawa_type": 1,
            "lambda7": 0.0,
        },
        "runtime_guards": {
            "ordering": ["mh < mH < mA"],
            "triple_ok_required": True,
        },
    }


def _to_float(value: object, field: str) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{field} must be numeric")
    if isinstance(value, (int, float, str)):
        out = float(value)
        if not math.isfinite(out):
            raise ValueError(f"{field} must be finite")
        return out
    raise TypeError(f"{field} must be numeric")


def _to_int(value: object, field: str) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{field} must be integer")
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if not value.is_integer():
            raise ValueError(f"{field} must be an integer value")
        return int(value)
    if isinstance(value, str):
        if value.strip() == "":
            raise ValueError(f"{field} must be integer")
        parsed = float(value)
        if not math.isfinite(parsed) or not parsed.is_integer():
            raise ValueError(f"{field} must be an integer value")
        return int(parsed)
    raise TypeError(f"{field} must be integer")


def normalize_contract(raw: Mapping[str, object]) -> dict[str, object]:
    out = cast(dict[str, object], copy.deepcopy(dict(raw)))

    variable_space_obj = out.get("variable_space")
    if isinstance(variable_space_obj, Mapping):
        variable_space = cast(dict[str, object], variable_space_obj)
        for axis in ("lambda6", "tan_beta", "lambda1", "mH"):
            axis_obj = variable_space.get(axis)
            if isinstance(axis_obj, Mapping):
                axis_map = cast(dict[str, object], axis_obj)
                if "min" in axis_map:
                    axis_map["min"] = _to_float(axis_map["min"], f"variable_space.{axis}.min")
                if "max" in axis_map:
                    axis_map["max"] = _to_float(axis_map["max"], f"variable_space.{axis}.max")

    fixed_obj = out.get("fixed")
    if isinstance(fixed_obj, Mapping):
        fixed = cast(dict[str, object], fixed_obj)
        float_fields = ("mh", "sin_ba", "mA", "mHp", "lambda7")
        for field in float_fields:
            if field in fixed:
                fixed[field] = _to_float(fixed[field], f"fixed.{field}")
        if "yukawa_type" in fixed:
            fixed["yukawa_type"] = _to_int(fixed["yukawa_type"], "fixed.yukawa_type")

    return out


def validate_contract(contract: Mapping[str, object]) -> list[str]:
    errors: list[str] = []

    for section in ("schema_version", "physics_contract", "variable_space", "fixed", "runtime_guards"):
        if section not in contract:
            errors.append(f"missing required section: {section}")

    schema_version = contract.get("schema_version")
    if schema_version is not None and not isinstance(schema_version, str):
        errors.append("schema_version must be a string")

    physics_contract_obj = contract.get("physics_contract")
    if physics_contract_obj is not None and not isinstance(physics_contract_obj, Mapping):
        errors.append("physics_contract must be a mapping")

    runtime_guards_obj = contract.get("runtime_guards")
    if runtime_guards_obj is not None and not isinstance(runtime_guards_obj, Mapping):
        errors.append("runtime_guards must be a mapping")

    variable_space_obj = contract.get("variable_space")
    if not isinstance(variable_space_obj, Mapping):
        errors.append("variable_space must be a mapping")
    else:
        variable_space = cast(Mapping[str, object], variable_space_obj)
        for axis in ("lambda6", "tan_beta", "lambda1", "mH"):
            axis_obj = variable_space.get(axis)
            if not isinstance(axis_obj, Mapping):
                errors.append(f"variable_space.{axis} must be a mapping")
                continue
            axis_map = cast(Mapping[str, object], axis_obj)
            if "min" not in axis_map or "max" not in axis_map:
                errors.append(f"variable_space.{axis} must define min and max")
                continue
            try:
                lo = _to_float(axis_map["min"], f"variable_space.{axis}.min")
                hi = _to_float(axis_map["max"], f"variable_space.{axis}.max")
            except (TypeError, ValueError) as exc:
                errors.append(str(exc))
                continue
            if hi <= lo:
                errors.append(f"variable_space.{axis} has invalid range: min >= max")

    fixed_obj = contract.get("fixed")
    if not isinstance(fixed_obj, Mapping):
        errors.append("fixed must be a mapping")
    else:
        fixed = cast(Mapping[str, object], fixed_obj)
        for field in ("mh", "sin_ba", "mA", "mHp", "yukawa_type", "lambda7"):
            if field not in fixed:
                errors.append(f"fixed.{field} is required")

        for field in ("mh", "sin_ba", "mA", "mHp", "lambda7"):
            if field in fixed:
                try:
                    _ = _to_float(fixed[field], f"fixed.{field}")
                except (TypeError, ValueError) as exc:
                    errors.append(str(exc))

        if "yukawa_type" in fixed:
            try:
                _ = _to_int(fixed["yukawa_type"], "fixed.yukawa_type")
            except (TypeError, ValueError) as exc:
                errors.append(str(exc))

    if isinstance(variable_space_obj, Mapping) and isinstance(fixed_obj, Mapping):
        try:
            m_h = _to_float(cast(Mapping[str, object], fixed_obj)["mh"], "fixed.mh")
            m_a = _to_float(cast(Mapping[str, object], fixed_obj)["mA"], "fixed.mA")
            m_h_bounds = cast(Mapping[str, object], cast(Mapping[str, object], variable_space_obj).get("mH", {}))
            m_h_lo = _to_float(m_h_bounds["min"], "variable_space.mH.min")
            m_h_hi = _to_float(m_h_bounds["max"], "variable_space.mH.max")
            if not (m_h < m_h_lo and m_h_hi < m_a):
                errors.append("mH bounds must satisfy fixed.mh < variable_space.mH.min and variable_space.mH.max < fixed.mA")
        except Exception:
            pass

    return errors
