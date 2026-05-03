from __future__ import annotations

# pyright: reportMissingImports=false, reportUnknownVariableType=false, reportUnknownMemberType=false, reportUnknownArgumentType=false

from autoresearch.harness.convergence import ConvergenceDetector


def _base_detector_config() -> dict[str, object]:
    return {
        "arms": [{"id": "adaptive-a"}],
        "convergence": {
            "warmup_rounds": 8,
            "window_rounds": 5,
            "confirmation_required": 2,
        },
    }


def _snapshot(
    round_index: int,
    *,
    arm_id: str = "adaptive-a",
    events_emitted: int = 1,
    discoveries: int = 1,
    yield_value: float,
    coverage: float,
    diversity: float,
    composite: float,
) -> dict[str, object]:
    return {
        "round_index": round_index,
        "arm_id": arm_id,
        "events_emitted": events_emitted,
        "discoveries": discoveries,
        "yield": yield_value,
        "coverage": coverage,
        "diversity": diversity,
        "composite": composite,
    }


def test_convergence_warmup_incomplete_before_required_rounds() -> None:
    detector = ConvergenceDetector(_base_detector_config())
    history = [
        _snapshot(i, yield_value=0.1 + i * 0.01, coverage=0.2 + i * 0.02, diversity=1.0 + i * 0.1, composite=0.3 + i * 0.03)
        for i in range(4)
    ]

    result = detector.evaluate(history, {"adaptive-a": 4})

    assert result["state"] == "warmup_incomplete"
    assert result["warmup_rounds"] == 8
    assert "4/8 rounds observed" in str(result["reason"])


def test_convergence_running_after_warmup_without_confirmation() -> None:
    detector = ConvergenceDetector(_base_detector_config())
    history = [
        _snapshot(i, yield_value=0.10 + i * 0.02, coverage=0.20 + i * 0.03, diversity=1.0 + i * 0.08, composite=0.30 + i * 0.05)
        for i in range(8)
    ]

    result = detector.evaluate(history, {"adaptive-a": 8})

    assert result["state"] == "running"
    assert result["plateau_confirmations"] == 0
    assert result["empty_confirmations"] == 0


def test_convergence_plateau_warning_on_first_confirmation() -> None:
    detector = ConvergenceDetector(_base_detector_config())
    history = [
        _snapshot(0, yield_value=0.10, coverage=0.20, diversity=1.00, composite=0.30),
        _snapshot(1, yield_value=0.12, coverage=0.24, diversity=1.08, composite=0.35),
        _snapshot(2, yield_value=0.14, coverage=0.28, diversity=1.16, composite=0.41),
        _snapshot(3, yield_value=0.16, coverage=0.32, diversity=1.24, composite=0.47),
        _snapshot(4, yield_value=0.161, coverage=0.321, diversity=1.25, composite=0.475),
        _snapshot(5, yield_value=0.162, coverage=0.322, diversity=1.26, composite=0.478),
        _snapshot(6, yield_value=0.163, coverage=0.323, diversity=1.27, composite=0.479),
        _snapshot(7, yield_value=0.164, coverage=0.324, diversity=1.28, composite=0.4795),
        _snapshot(8, yield_value=0.165, coverage=0.325, diversity=1.29, composite=0.48),
    ]

    result = detector.evaluate(history, {"adaptive-a": 9})

    assert result["state"] == "plateau_warning"
    assert result["plateau_confirmations"] == 1
    assert "coverage<=0.01" in str(result["reason"])
    assert "diversity<=0.05" in str(result["reason"])


def test_convergence_plateau_after_second_confirmation() -> None:
    detector = ConvergenceDetector(_base_detector_config())
    history = [
        _snapshot(0, yield_value=0.10, coverage=0.20, diversity=1.00, composite=0.30),
        _snapshot(1, yield_value=0.12, coverage=0.24, diversity=1.08, composite=0.35),
        _snapshot(2, yield_value=0.14, coverage=0.28, diversity=1.16, composite=0.41),
        _snapshot(3, yield_value=0.16, coverage=0.32, diversity=1.24, composite=0.47),
        _snapshot(4, yield_value=0.161, coverage=0.321, diversity=1.25, composite=0.475),
        _snapshot(5, yield_value=0.162, coverage=0.322, diversity=1.26, composite=0.478),
        _snapshot(6, yield_value=0.163, coverage=0.323, diversity=1.27, composite=0.479),
        _snapshot(7, yield_value=0.164, coverage=0.324, diversity=1.28, composite=0.4795),
        _snapshot(8, yield_value=0.165, coverage=0.325, diversity=1.29, composite=0.48),
        _snapshot(9, yield_value=0.165, coverage=0.325, diversity=1.29, composite=0.48),
    ]

    result = detector.evaluate(history, {"adaptive-a": 10})

    assert result["state"] == "converged_plateau"
    assert result["plateau_confirmations"] == 2
    assert "5-round trailing window" in str(result["reason"])


def test_convergence_empty_exhaustion_after_second_confirmation() -> None:
    detector = ConvergenceDetector(_base_detector_config())
    history = [
        _snapshot(i, yield_value=0.10 + i * 0.05, coverage=0.20 + i * 0.05, diversity=1.00 + i * 0.20, composite=0.30 + i * 0.07)
        for i in range(8)
    ]
    history.extend(
        [
            _snapshot(8, events_emitted=0, discoveries=2, yield_value=0.90, coverage=0.90, diversity=2.50, composite=0.95),
            _snapshot(9, events_emitted=0, discoveries=1, yield_value=0.90, coverage=0.90, diversity=2.50, composite=0.95),
        ]
    )

    result = detector.evaluate(history, {"adaptive-a": 10})

    assert result["state"] == "converged_empty_exhaustion"
    assert result["empty_confirmations"] == 2
    assert "discoveries=1 kept for diagnostics" in str(result["reason"])


def test_convergence_requires_every_arm_seen_before_converging() -> None:
    detector = ConvergenceDetector(
        {
            "arms": [{"id": "adaptive-a"}, {"id": "branch-b"}],
            "convergence": {
                "warmup_rounds": 8,
                "window_rounds": 5,
                "confirmation_required": 2,
            },
        }
    )
    history = [
        _snapshot(i, yield_value=0.50, coverage=0.50, diversity=1.20, composite=0.71)
        for i in range(8)
    ]

    result = detector.evaluate(history, {"adaptive-a": 8, "branch-b": 0})

    assert result["state"] == "warmup_incomplete"
    assert "missing pulls for branch-b" in str(result["reason"])
