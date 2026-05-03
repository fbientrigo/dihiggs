from __future__ import annotations
# pyright: reportMissingImports=false, reportUnknownVariableType=false, reportUnknownMemberType=false, reportUnknownArgumentType=false

from autoresearch.harness.autoscaling import AutoscalingPolicy


def _state(**overrides: object) -> dict[str, object]:
    state: dict[str, object] = {
        "round_index": 6,
        "threads": 1,
        "timeout_sec": 900,
        "max_new_run_dirs_per_round": 20,
        "max_new_run_dirs_per_arm_call": 10,
        "empty_rounds": 0,
        "plateau_warning": False,
        "coverage_gain": 0.04,
        "diversity_gain": 0.10,
        "host_pressure": "low",
    }
    state.update(overrides)
    return state


def test_healthy_round_scales_up_then_holds_during_cooldown() -> None:
    policy = AutoscalingPolicy(
        {
            "min_threads": 1,
            "max_threads": 4,
            "cooldown_rounds": 3,
            "scale_up_factor": 1.5,
            "max_timeout_sec": 1800,
            "max_max_new_run_dirs_per_round": 30,
            "max_max_new_run_dirs_per_arm_call": 15,
        }
    )

    first = policy.evaluate(_state())

    assert first["action"] == "scale_up"
    assert first["cause"] == "healthy_rounds_with_positive_coverage_and_diversity_gain"
    assert first["next_threads"] == 2
    assert first["next_timeout_sec"] == 1350
    assert first["next_max_new_run_dirs_per_round"] == 30
    assert first["next_max_new_run_dirs_per_arm_call"] == 15
    assert first["changes"]["threads"] == {"before": 1, "after": 2}

    second = policy.evaluate(_state(round_index=7, threads=2, timeout_sec=1350, max_new_run_dirs_per_round=30, max_new_run_dirs_per_arm_call=15))

    assert second["action"] == "hold"
    assert second["cooldown_active"] is True
    assert second["cause"].startswith("cooldown_active:")
    assert second["next_threads"] == 2
    assert second["next_timeout_sec"] == 1350


def test_plateau_and_high_pressure_scale_down() -> None:
    policy = AutoscalingPolicy(
        {
            "min_threads": 1,
            "max_threads": 4,
            "cooldown_rounds": 3,
            "scale_down_factor": 0.67,
            "min_timeout_sec": 600,
            "min_max_new_run_dirs_per_round": 8,
            "min_max_new_run_dirs_per_arm_call": 4,
        }
    )

    result = policy.evaluate(
        _state(
            threads=3,
            empty_rounds=2,
            plateau_warning=True,
            coverage_gain=0.0,
            diversity_gain=0.0,
            host_pressure="high",
        )
    )

    assert result["action"] == "scale_down"
    assert result["cause"] == "plateau_warning, repeated_empty_rounds, host_pressure:high"
    assert result["next_threads"] == 2
    assert result["next_timeout_sec"] == 603
    assert result["next_max_new_run_dirs_per_round"] == 13
    assert result["next_max_new_run_dirs_per_arm_call"] == 6


def test_scale_up_respects_bounds_and_arm_budget_never_exceeds_round_budget() -> None:
    policy = AutoscalingPolicy(
        {
            "min_threads": 1,
            "max_threads": 4,
            "cooldown_rounds": 3,
            "scale_up_factor": 2.0,
            "max_timeout_sec": 1600,
            "max_max_new_run_dirs_per_round": 32,
            "max_max_new_run_dirs_per_arm_call": 40,
        }
    )

    result = policy.evaluate(
        _state(
            threads=4,
            timeout_sec=1500,
            max_new_run_dirs_per_round=25,
            max_new_run_dirs_per_arm_call=24,
        )
    )

    assert result["action"] == "scale_up"
    assert result["next_threads"] == 4
    assert result["next_timeout_sec"] == 1600
    assert result["next_max_new_run_dirs_per_round"] == 32
    assert result["next_max_new_run_dirs_per_arm_call"] == 32


def test_cooldown_suppresses_scale_down_signal() -> None:
    policy = AutoscalingPolicy(
        {
            "min_threads": 1,
            "max_threads": 4,
            "cooldown_rounds": 3,
        }
    )

    first = policy.evaluate(_state(round_index=2, threads=2, timeout_sec=1000, max_new_run_dirs_per_round=16, max_new_run_dirs_per_arm_call=8))
    second = policy.evaluate(
        _state(
            round_index=4,
            threads=first["next_threads"],
            timeout_sec=first["next_timeout_sec"],
            max_new_run_dirs_per_round=first["next_max_new_run_dirs_per_round"],
            max_new_run_dirs_per_arm_call=first["next_max_new_run_dirs_per_arm_call"],
            empty_rounds=2,
            plateau_warning=True,
            coverage_gain=0.0,
            diversity_gain=0.0,
            host_pressure="high",
        )
    )

    assert first["action"] == "scale_up"
    assert second["action"] == "hold"
    assert second["cooldown_active"] is True
    assert second["next_threads"] == first["next_threads"]
