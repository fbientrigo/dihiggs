from __future__ import annotations

import copy
import math
import re
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from typing import cast


PLACEHOLDER_PATTERN = re.compile(r"\{([a-zA-Z_][a-zA-Z0-9_]*)\}")


@dataclass(frozen=True)
class AxisContract:
    physical_axes: tuple[str, ...] = ("tb", "lam1_bin", "mphi_bin")
    meta_axes: tuple[str, ...] = ("explorer", "strategy")
    required_cell_axes: tuple[str, ...] = ("tb", "lam1_bin")

    @property
    def all_axes(self) -> tuple[str, ...]:
        return self.physical_axes + self.meta_axes


def _coerce_int_if_possible(value: str) -> str | int:
    if re.fullmatch(r"-?\d+", value):
        try:
            return int(value)
        except ValueError:
            return value
    return value


def encode_cell_id(axes_binned: Mapping[str, object]) -> str:
    if "tb" not in axes_binned:
        raise ValueError("Missing required axis: tb")
    if "lam1_bin" not in axes_binned:
        raise ValueError("Missing required axis: lam1_bin")

    parts = [f"tb={axes_binned['tb']}", f"bin={axes_binned['lam1_bin']}"]

    for key in sorted(k for k in axes_binned.keys() if k not in {"tb", "lam1_bin"}):
        value = axes_binned[key]
        if value is None:
            continue
        parts.append(f"{key}={value}")

    return "|".join(parts)


def decode_cell_id(cell_id: str) -> dict[str, str | int]:
    if not cell_id:
        raise ValueError("cell_id must be a non-empty string")

    out: dict[str, str | int] = {}
    for segment in cell_id.split("|"):
        if "=" not in segment:
            raise ValueError(f"Invalid segment in cell_id: {segment!r}")
        key, raw_value = segment.split("=", 1)
        if not key:
            raise ValueError(f"Invalid key in cell_id segment: {segment!r}")

        parsed_key = "lam1_bin" if key == "bin" else key
        out[parsed_key] = _coerce_int_if_possible(raw_value)

    if "tb" not in out or "lam1_bin" not in out:
        raise ValueError("cell_id must contain both tb and bin/lam1_bin")

    return out


def _linear_bin(raw_value: float, bin_cfg: Mapping[str, object], axis_name: str) -> int:
    def to_float(value: object, field: str) -> float:
        if isinstance(value, (int, float, str)):
            return float(value)
        raise TypeError(f"{axis_name}.{field} must be numeric")

    def to_int(value: object, field: str) -> int:
        if isinstance(value, (int, float, str)):
            return int(value)
        raise TypeError(f"{axis_name}.{field} must be integer-like")

    try:
        lo = to_float(bin_cfg["min"], "min")
        hi = to_float(bin_cfg["max"], "max")
        n_bins = to_int(bin_cfg["n_bins"], "n_bins")
    except KeyError as exc:
        raise KeyError(f"Missing {axis_name} binning config key: {exc.args[0]}") from exc

    if n_bins <= 0:
        raise ValueError(f"{axis_name}.n_bins must be > 0")
    if hi <= lo:
        raise ValueError(f"{axis_name}.max must be greater than {axis_name}.min")

    if raw_value <= lo:
        return 0
    if raw_value >= hi:
        return n_bins - 1

    width = (hi - lo) / n_bins
    idx = math.floor((raw_value - lo) / width)
    return max(0, min(n_bins - 1, idx))


def _expect_mapping(value: object, name: str) -> Mapping[str, object]:
    if isinstance(value, Mapping):
        return cast(Mapping[str, object], value)
    raise TypeError(f"Expected mapping for {name}")


def bin_lam1(lam1_raw: float, config: Mapping[str, object]) -> int:
    search = _expect_mapping(config["search"], "search")
    lam1_cfg = _expect_mapping(search["lam1"], "search.lam1")
    return _linear_bin(float(lam1_raw), lam1_cfg, "lam1")


def bin_mphi(mphi_raw: float, config: Mapping[str, object]) -> int:
    search = _expect_mapping(config["search"], "search")
    mphi_cfg = _expect_mapping(search["mphi"], "search.mphi")
    return _linear_bin(float(mphi_raw), mphi_cfg, "mphi")


def substitute_placeholders(
    config: Mapping[str, object], runtime_ctx: Mapping[str, object]
) -> dict[str, object]:
    src = copy.deepcopy(dict(config))

    ctx: dict[str, object] = dict(runtime_ctx)
    _ = ctx.setdefault("sys_executable", sys.executable)
    _ = ctx.setdefault("python", ctx.get("sys_executable"))

    resolving: set[str] = set()

    def resolve_ctx_key(name: str) -> object:
        if name not in ctx:
            raise ValueError(f"Missing placeholder: {name}")
        if name in resolving:
            raise ValueError(f"Circular placeholder reference detected: {name}")

        value = ctx[name]
        if isinstance(value, str):
            resolving.add(name)
            value = resolve_string(value)
            resolving.remove(name)
            ctx[name] = value
        return value

    def replace_match(match: re.Match[str]) -> str:
        placeholder = match.group(1)
        resolved = resolve_ctx_key(placeholder)
        return str(resolved)

    def resolve_string(value: str) -> str:
        current = value
        for _ in range(20):
            replaced = PLACEHOLDER_PATTERN.sub(replace_match, current)
            if replaced == current:
                return replaced
            current = replaced
        raise ValueError(f"Exceeded placeholder resolution depth for value: {value!r}")

    def walk(value: object) -> object:
        if isinstance(value, dict):
            raw_dict = cast(dict[object, object], value)
            out_dict: dict[str, object] = {}
            for k, v in raw_dict.items():
                out_dict[str(k)] = walk(v)
            return out_dict
        if isinstance(value, list):
            raw_list = cast(list[object], value)
            return [walk(v) for v in raw_list]
        if isinstance(value, str):
            return resolve_string(value)
        return value

    resolved = walk(src)
    if not isinstance(resolved, dict):
        raise TypeError("Resolved config must be a dictionary")
    return cast(dict[str, object], resolved)
