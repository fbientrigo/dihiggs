from __future__ import annotations

import math
from collections.abc import Mapping


TARGET_KNOBS = (
    "threads",
    "timeout_sec",
    "max_new_run_dirs_per_round",
    "max_new_run_dirs_per_arm_call",
)


def _as_bool(value: object) -> bool:
    return value if isinstance(value, bool) else False


def _as_int(value: object, default: int = 0) -> int:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    if isinstance(value, str):
        try:
            return int(value)
        except ValueError:
            return default
    return default


def _as_float(value: object, default: float = 0.0) -> float:
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return default
    return default


def _as_str(value: object, default: str = "") -> str:
    return value if isinstance(value, str) else default


def _bound_key(knob: str, side: str) -> str:
    return f"{side}_{knob}"


class AutoscalingPolicy:
    def __init__(self, config: Mapping[str, object] | None = None):
        self.config: dict[str, object] = dict(config or {})
        self._last_action_round: int | None = None

    def evaluate(self, state: Mapping[str, object]) -> dict[str, object]:
        round_index = _as_int(state.get("round_index", 0))
        cause = self._decide_cause(state)

        cooldown_round = self._cooldown_round(state)
        cooldown_rounds = _as_int(self.config.get("cooldown_rounds", 3), 3)
        cooldown_active = (
            cooldown_round is not None
            and cooldown_rounds > 0
            and round_index - cooldown_round <= cooldown_rounds
        )
        if cooldown_active and cause != "hold":
            assert cooldown_round is not None
            blocked_for = cooldown_rounds - (round_index - cooldown_round) + 1
            return self._build_result(
                action="hold",
                cause=f"cooldown_active:{max(blocked_for, 0)}_rounds_remaining",
                state=state,
                round_index=round_index,
                cooldown_active=True,
            )

        action = cause
        if action == "scale_up":
            result = self._build_result(
                action="scale_up",
                cause="healthy_rounds_with_positive_coverage_and_diversity_gain",
                state=state,
                round_index=round_index,
                cooldown_active=False,
            )
        elif action == "scale_down":
            result = self._build_result(
                action="scale_down",
                cause=self._scale_down_cause(state),
                state=state,
                round_index=round_index,
                cooldown_active=False,
            )
        else:
            result = self._build_result(
                action="hold",
                cause="no_scaling_signal",
                state=state,
                round_index=round_index,
                cooldown_active=False,
            )

        if result["action"] in {"scale_up", "scale_down"}:
            self._last_action_round = round_index
        return result

    def _decide_cause(self, state: Mapping[str, object]) -> str:
        plateau_warning = _as_bool(state.get("plateau_warning", False))
        empty_rounds = _as_int(state.get("empty_rounds", 0))
        coverage_gain = _as_float(state.get("coverage_gain", 0.0))
        diversity_gain = _as_float(state.get("diversity_gain", 0.0))
        host_pressure = _as_str(state.get("host_pressure", "unknown"), "unknown").lower()
        empty_round_threshold = _as_int(self.config.get("empty_rounds_scale_down_threshold", 2), 2)

        if plateau_warning or empty_rounds >= empty_round_threshold or host_pressure in {"high", "critical"}:
            return "scale_down"

        healthy_growth = (
            host_pressure == "low"
            and empty_rounds == 0
            and not plateau_warning
            and coverage_gain > 0.0
            and diversity_gain > 0.0
        )
        if healthy_growth:
            return "scale_up"
        return "hold"

    def _scale_down_cause(self, state: Mapping[str, object]) -> str:
        reasons: list[str] = []
        if _as_bool(state.get("plateau_warning", False)):
            reasons.append("plateau_warning")
        if _as_int(state.get("empty_rounds", 0)) >= _as_int(self.config.get("empty_rounds_scale_down_threshold", 2), 2):
            reasons.append("repeated_empty_rounds")
        host_pressure = _as_str(state.get("host_pressure", "unknown"), "unknown").lower()
        if host_pressure in {"high", "critical"}:
            reasons.append(f"host_pressure:{host_pressure}")
        return ", ".join(reasons) if reasons else "protective_scale_down"

    def _cooldown_round(self, state: Mapping[str, object]) -> int | None:
        state_round = state.get("last_scaling_round")
        if state_round is None:
            return self._last_action_round
        parsed = _as_int(state_round)
        if self._last_action_round is None:
            return parsed
        return max(parsed, self._last_action_round)

    def _build_result(
        self,
        *,
        action: str,
        cause: str,
        state: Mapping[str, object],
        round_index: int,
        cooldown_active: bool,
    ) -> dict[str, object]:
        current = {knob: max(_as_int(state.get(knob, 0)), 1) for knob in TARGET_KNOBS}
        next_values = dict(current)

        if action == "scale_up":
            factor = max(_as_float(self.config.get("scale_up_factor", 1.5), 1.5), 1.0)
            next_values = {knob: self._scale_value(knob, current[knob], factor, direction="up") for knob in TARGET_KNOBS}
        elif action == "scale_down":
            factor = min(_as_float(self.config.get("scale_down_factor", 0.67), 0.67), 1.0)
            next_values = {knob: self._scale_value(knob, current[knob], factor, direction="down") for knob in TARGET_KNOBS}

        next_values["max_new_run_dirs_per_arm_call"] = min(
            next_values["max_new_run_dirs_per_arm_call"],
            next_values["max_new_run_dirs_per_round"],
        )

        changes = {
            knob: {"before": current[knob], "after": next_values[knob]}
            for knob in TARGET_KNOBS
        }
        return {
            "round_index": round_index,
            "action": action,
            "cause": cause,
            "cooldown_active": cooldown_active,
            "last_scaling_round": self._last_action_round if action == "hold" else round_index,
            "next_threads": next_values["threads"],
            "next_timeout_sec": next_values["timeout_sec"],
            "next_max_new_run_dirs_per_round": next_values["max_new_run_dirs_per_round"],
            "next_max_new_run_dirs_per_arm_call": next_values["max_new_run_dirs_per_arm_call"],
            "changes": changes,
            "signals": {
                "empty_rounds": _as_int(state.get("empty_rounds", 0)),
                "plateau_warning": _as_bool(state.get("plateau_warning", False)),
                "coverage_gain": _as_float(state.get("coverage_gain", 0.0)),
                "diversity_gain": _as_float(state.get("diversity_gain", 0.0)),
                "host_pressure": _as_str(state.get("host_pressure", "unknown"), "unknown"),
            },
        }

    def _scale_value(self, knob: str, current: int, factor: float, *, direction: str) -> int:
        minimum = self._min_bound(knob)
        maximum = self._max_bound(knob, current)
        if direction == "up":
            raw = math.ceil(current * factor)
            if raw <= current and current < maximum:
                raw = current + 1
        else:
            raw = math.floor(current * factor)
            if raw >= current and current > minimum:
                raw = current - 1
        bounded = min(max(raw, minimum), maximum)
        return bounded

    def _min_bound(self, knob: str) -> int:
        if knob == "threads":
            default = 1
        else:
            default = 1
        return max(_as_int(self.config.get(_bound_key(knob, "min"), default), default), 1)

    def _max_bound(self, knob: str, current: int) -> int:
        if knob == "threads":
            default = max(current, 1)
        elif knob == "timeout_sec":
            default = max(current, 7200)
        else:
            default = max(current, 100)
        return max(_as_int(self.config.get(_bound_key(knob, "max"), default), default), self._min_bound(knob))


__all__ = ["AutoscalingPolicy", "TARGET_KNOBS"]
