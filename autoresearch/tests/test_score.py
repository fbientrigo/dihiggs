"""
Comprehensive tests for autoresearch.benchmarks.score module.

Tests cover:
- compute_metrics: valid CSV, missing CSV, malformed CSV, empty CSV
- compute_coverage: single-axis, multi-axis, empty history, full coverage
- compute_diversity: uniform distribution, skewed distribution, empty history, pairwise
- compute_composite: various weight combinations
- End-to-end: full pipeline from history to composite score
"""

from __future__ import annotations

import csv
import os
from collections.abc import Mapping
from typing import Any

import pytest

from autoresearch.benchmarks.score import (
    compute_composite,
    compute_coverage,
    compute_diversity,
    compute_metrics,
)


@pytest.fixture
def basic_config() -> dict[str, Any]:
    return {
        "search": {
            "tb_values": [1000, 5000, 10000, 30000],
            "lam1": {"min": -20.0, "max": 20.0, "n_bins": 40},
            "mphi": {"min": 0.0, "max": 1000.0, "n_bins": 50},
        },
        "metrics": {
            "weights": {"yield": 0.5, "coverage": 0.3, "diversity": 0.2},
            "multi_axis": {
                "collapse_axes": ["tb", "lam1_bin"],
                "coverage_axes": [
                    {"name": "tb", "kind": "categorical", "domain": [1000, 5000, 10000, 30000], "weight": 0.5},
                    {"name": "lam1_bin", "kind": "discrete", "domain_size": 40, "weight": 0.5},
                ],
                "diversity_axes": [
                    {"name": "tb", "weight": 0.5},
                    {"name": "lam1_bin", "weight": 0.5},
                ],
                "diversity_pairs": [
                    {"axes": ["tb", "lam1_bin"], "weight": 1.0},
                ],
            },
        },
    }


def test_compute_metrics_valid_csv(tmp_path: Any, basic_config: dict[str, Any]) -> None:
    run_dir = tmp_path / "run1"
    run_dir.mkdir()
    csv_path = run_dir / "results.csv"
    
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["success", "total_events"])
        writer.writerow(["True", "1000"])
        writer.writerow(["False", "1000"])
        writer.writerow(["True", "1500"])
        writer.writerow(["True", "2000"])
    
    result = compute_metrics(str(run_dir), {}, basic_config)
    
    assert result["successes"] == 3
    assert result["trials"] == 4
    assert result["yield_val"] == 0.75


def test_compute_metrics_missing_csv(tmp_path: Any, basic_config: dict[str, Any]) -> None:
    run_dir = tmp_path / "run_missing"
    run_dir.mkdir()
    
    result = compute_metrics(str(run_dir), {}, basic_config)
    
    assert result["yield_val"] == 0.0
    assert result["successes"] == 0
    assert result["trials"] == 0


def test_compute_metrics_malformed_csv(tmp_path: Any, basic_config: dict[str, Any]) -> None:
    run_dir = tmp_path / "run_malformed"
    run_dir.mkdir()
    csv_path = run_dir / "results.csv"
    
    with open(csv_path, "w") as f:
        f.write("success,total_events\n")
        f.write("not_a_boolean,not_a_number\n")
    
    result = compute_metrics(str(run_dir), {}, basic_config)
    
    assert result["yield_val"] == 0.0
    assert result["successes"] == 0
    assert result["trials"] == 1


def test_compute_metrics_empty_csv(tmp_path: Any, basic_config: dict[str, Any]) -> None:
    run_dir = tmp_path / "run_empty"
    run_dir.mkdir()
    csv_path = run_dir / "results.csv"
    
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["success", "total_events"])
    
    result = compute_metrics(str(run_dir), {}, basic_config)
    
    assert result["yield_val"] == 0.0
    assert result["successes"] == 0
    assert result["trials"] == 0


def test_compute_metrics_various_success_formats(tmp_path: Any, basic_config: dict[str, Any]) -> None:
    run_dir = tmp_path / "run_formats"
    run_dir.mkdir()
    csv_path = run_dir / "results.csv"
    
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["success", "total_events"])
        writer.writerow(["true", "1000"])
        writer.writerow(["TRUE", "1000"])
        writer.writerow(["1", "1000"])
        writer.writerow(["yes", "1000"])
        writer.writerow(["False", "1000"])
        writer.writerow(["0", "1000"])
    
    result = compute_metrics(str(run_dir), {}, basic_config)
    
    assert result["successes"] == 4
    assert result["trials"] == 6


def test_compute_coverage_single_axis(basic_config: dict[str, Any]) -> None:
    history = [
        {"axes_binned": {"lam1_bin": 0}},
        {"axes_binned": {"lam1_bin": 1}},
        {"axes_binned": {"lam1_bin": 0}},
    ]
    
    config_single = {
        **basic_config,
        "metrics": {
            **basic_config["metrics"],
            "multi_axis": {
                "coverage_axes": [
                    {"name": "lam1_bin", "kind": "discrete", "domain_size": 40, "weight": 1.0},
                ],
            },
        },
    }
    
    coverage = compute_coverage(history, config_single)
    
    assert coverage == 2 / 40


def test_compute_coverage_multi_axis(basic_config: dict[str, Any]) -> None:
    history = [
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 5000, "lam1_bin": 1}},
        {"axes_binned": {"tb": 1000, "lam1_bin": 2}},
    ]
    
    coverage = compute_coverage(history, basic_config)
    
    tb_coverage = 2 / 4
    lam1_coverage = 3 / 40
    expected = 0.5 * tb_coverage + 0.5 * lam1_coverage
    
    assert coverage == pytest.approx(expected)


def test_compute_coverage_empty_history(basic_config: dict[str, Any]) -> None:
    history: list[dict[str, Any]] = []
    
    coverage = compute_coverage(history, basic_config)
    
    assert coverage == 0.0


def test_compute_coverage_full_coverage_single_axis(basic_config: dict[str, Any]) -> None:
    history = [{"axes_binned": {"tb": tb}} for tb in [1000, 5000, 10000, 30000]]
    
    config_single_axis = {
        **basic_config,
        "metrics": {
            **basic_config["metrics"],
            "multi_axis": {
                "coverage_axes": [
                    {"name": "tb", "kind": "categorical", "domain": [1000, 5000, 10000, 30000], "weight": 1.0},
                ],
            },
        },
    }
    
    coverage = compute_coverage(history, config_single_axis)
    
    assert coverage == 1.0


def test_compute_coverage_weighted_axes(basic_config: dict[str, Any]) -> None:
    history = [
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 5000, "lam1_bin": 1}},
    ]
    
    config_weighted = {
        **basic_config,
        "metrics": {
            **basic_config["metrics"],
            "multi_axis": {
                "coverage_axes": [
                    {"name": "tb", "kind": "categorical", "domain": [1000, 5000, 10000, 30000], "weight": 0.8},
                    {"name": "lam1_bin", "kind": "discrete", "domain_size": 40, "weight": 0.2},
                ],
            },
        },
    }
    
    coverage = compute_coverage(history, config_weighted)
    
    tb_coverage = 2 / 4
    lam1_coverage = 2 / 40
    expected = (0.8 * tb_coverage + 0.2 * lam1_coverage) / (0.8 + 0.2)
    
    assert coverage == pytest.approx(expected)


def test_compute_diversity_uniform_distribution(basic_config: dict[str, Any]) -> None:
    history = [
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 5000, "lam1_bin": 1}},
        {"axes_binned": {"tb": 10000, "lam1_bin": 2}},
        {"axes_binned": {"tb": 30000, "lam1_bin": 3}},
    ]
    
    diversity = compute_diversity(history, basic_config)
    
    assert diversity > 0.0


def test_compute_diversity_skewed_distribution(basic_config: dict[str, Any]) -> None:
    history = [
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 5000, "lam1_bin": 1}},
    ]
    
    diversity = compute_diversity(history, basic_config)
    
    assert diversity > 0.0


def test_compute_diversity_empty_history(basic_config: dict[str, Any]) -> None:
    history: list[dict[str, Any]] = []
    
    diversity = compute_diversity(history, basic_config)
    
    assert diversity == 0.0


def test_compute_diversity_single_value(basic_config: dict[str, Any]) -> None:
    history = [
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
    ]
    
    diversity = compute_diversity(history, basic_config)
    
    assert diversity == 0.0


def test_compute_diversity_with_pairwise(basic_config: dict[str, Any]) -> None:
    history = [
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 1000, "lam1_bin": 1}},
        {"axes_binned": {"tb": 5000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 5000, "lam1_bin": 1}},
    ]
    
    diversity = compute_diversity(history, basic_config)
    
    assert diversity > 0.0


def test_compute_diversity_uniform_higher_than_skewed(basic_config: dict[str, Any]) -> None:
    history_uniform = [
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 5000, "lam1_bin": 1}},
        {"axes_binned": {"tb": 10000, "lam1_bin": 2}},
        {"axes_binned": {"tb": 30000, "lam1_bin": 3}},
    ]
    
    history_skewed = [
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 5000, "lam1_bin": 1}},
    ]
    
    diversity_uniform = compute_diversity(history_uniform, basic_config)
    diversity_skewed = compute_diversity(history_skewed, basic_config)
    
    assert diversity_uniform > diversity_skewed


def test_compute_composite_basic(basic_config: dict[str, Any]) -> None:
    composite = compute_composite(1.0, 0.5, 0.3, basic_config)
    
    expected = 0.5 * 1.0 + 0.3 * 0.5 + 0.2 * 0.3
    assert composite == pytest.approx(expected)


def test_compute_composite_zero_weights() -> None:
    config = {
        "metrics": {
            "weights": {"yield": 0.0, "coverage": 0.0, "diversity": 0.0},
        },
    }
    
    composite = compute_composite(1.0, 1.0, 1.0, config)
    
    assert composite == 0.0


def test_compute_composite_various_weights() -> None:
    config = {
        "metrics": {
            "weights": {"yield": 0.7, "coverage": 0.2, "diversity": 0.1},
        },
    }
    
    composite = compute_composite(0.8, 0.6, 0.4, config)
    
    expected = 0.7 * 0.8 + 0.2 * 0.6 + 0.1 * 0.4
    assert composite == pytest.approx(expected)


def test_end_to_end_pipeline(tmp_path: Any, basic_config: dict[str, Any]) -> None:
    run_dir = tmp_path / "run_e2e"
    run_dir.mkdir()
    csv_path = run_dir / "results.csv"
    
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["success", "total_events"])
        writer.writerow(["True", "1000"])
        writer.writerow(["True", "1500"])
    
    history = [
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 5000, "lam1_bin": 1}},
        {"axes_binned": {"tb": 1000, "lam1_bin": 2}},
    ]
    
    metrics_result = compute_metrics(str(run_dir), {}, basic_config)
    coverage = compute_coverage(history, basic_config)
    diversity = compute_diversity(history, basic_config)
    composite = compute_composite(
        metrics_result["yield_val"],
        coverage,
        diversity,
        basic_config,
    )
    
    assert metrics_result["yield_val"] == 1.0
    assert coverage > 0.0
    assert diversity > 0.0
    assert composite > 0.0


def test_compute_coverage_fallback_to_collapse_axes() -> None:
    history = [
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 5000, "lam1_bin": 1}},
    ]
    
    config_no_coverage_axes = {
        "search": {
            "tb_values": [1000, 5000, 10000, 30000],
            "lam1": {"min": -20.0, "max": 20.0, "n_bins": 40},
        },
        "metrics": {
            "multi_axis": {
                "collapse_axes": ["tb", "lam1_bin"],
            },
        },
    }
    
    coverage = compute_coverage(history, config_no_coverage_axes)
    
    assert coverage > 0.0


def test_compute_diversity_fallback_to_collapse_axes() -> None:
    history = [
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}},
        {"axes_binned": {"tb": 5000, "lam1_bin": 1}},
    ]
    
    config_no_diversity_axes = {
        "search": {
            "tb_values": [1000, 5000, 10000, 30000],
            "lam1": {"min": -20.0, "max": 20.0, "n_bins": 40},
        },
        "metrics": {
            "multi_axis": {
                "collapse_axes": ["tb", "lam1_bin"],
            },
        },
    }
    
    diversity = compute_diversity(history, config_no_diversity_axes)
    
    assert diversity > 0.0
