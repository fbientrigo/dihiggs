"""
Multi-axis coverage and diversity metrics for DiHiggs autoresearch campaigns.

This module provides sparse axis tracking with per-axis bins and selected 
pairwise interactions for computing coverage, diversity, and composite scores.

Key functions:
- compute_metrics: Extract yield metrics from run_dir CSV
- compute_coverage: Per-axis bin coverage with weighted mean
- compute_diversity: Shannon entropy over per-axis distributions
- compute_composite: Weighted sum of yield, coverage, diversity
"""

from __future__ import annotations

import csv
import logging
import math
import os
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Any, TypedDict

logger = logging.getLogger(__name__)


class MetricsResult(TypedDict):
    """Result structure for compute_metrics."""
    yield_val: float
    successes: int
    trials: int


def compute_metrics(
    run_dir: str, 
    axes_binned: Mapping[str, Any], 
    config: Mapping[str, Any]
) -> MetricsResult:
    """
    Extract yield metrics from a run directory's results.csv.
    
    CSV format expected:
        success,total_events
        True,1000
        False,1000
        True,1500
    
    Args:
        run_dir: Path to run directory containing results.csv
        axes_binned: Binned axis values (unused, reserved for future use)
        config: Configuration dict (unused, reserved for future use)
    
    Returns:
        Dict with keys: yield_val, successes, trials
        
    Graceful degradation:
        - Missing file → yield=0.0, successes=0, trials=0
        - Malformed CSV → yield=0.0, successes=0, trials=0
    """
    csv_path = os.path.join(run_dir, "results.csv")
    
    if not os.path.exists(csv_path):
        logger.warning(f"Missing results.csv in {run_dir}, returning zero metrics")
        return {"yield_val": 0.0, "successes": 0, "trials": 0}
    
    try:
        with open(csv_path, "r") as f:
            reader = csv.DictReader(f)
            successes = 0
            trials = 0
            
            for row in reader:
                trials += 1
                success_str = row.get("success", "").strip().lower()
                if success_str in ("true", "1", "yes"):
                    successes += 1
            
            if trials == 0:
                logger.warning(f"Empty results.csv in {run_dir}")
                return {"yield_val": 0.0, "successes": 0, "trials": 0}
            
            yield_val = successes / trials
            return {"yield_val": yield_val, "successes": successes, "trials": trials}
    
    except (IOError, csv.Error, KeyError) as exc:
        logger.warning(f"Failed to parse results.csv in {run_dir}: {exc}")
        return {"yield_val": 0.0, "successes": 0, "trials": 0}


def _get_axis_domain_size(axis_name: str, config: Mapping[str, Any]) -> int:
    """
    Get domain size for an axis from config.
    
    For categorical axes (tb), returns number of values.
    For discrete axes (lam1_bin, mphi_bin), returns n_bins.
    """
    # Check coverage_axes first (preferred source)
    multi_axis = config.get("metrics", {}).get("multi_axis", {})
    coverage_axes = multi_axis.get("coverage_axes", [])
    
    for axis_spec in coverage_axes:
        if axis_spec.get("name") == axis_name:
            kind = axis_spec.get("kind")
            if kind == "categorical":
                domain = axis_spec.get("domain", [])
                return len(domain)
            elif kind == "discrete":
                return int(axis_spec.get("domain_size", 1))
    
    # Fallback: check search config
    search = config.get("search", {})
    
    if axis_name == "tb":
        tb_values = search.get("tb_values", [])
        return len(tb_values) if tb_values else 1
    
    if axis_name == "lam1_bin":
        lam1_cfg = search.get("lam1", {})
        return int(lam1_cfg.get("n_bins", 1))
    
    if axis_name == "mphi_bin":
        mphi_cfg = search.get("mphi", {})
        return int(mphi_cfg.get("n_bins", 1))
    
    # Default for unknown axes
    logger.warning(f"Unknown axis {axis_name}, defaulting domain_size=1")
    return 1


def compute_coverage(
    history: Sequence[Mapping[str, Any]], 
    config: Mapping[str, Any]
) -> float:
    """
    Compute per-axis bin coverage with weighted mean.
    
    Coverage formula per axis:
        bins_visited / total_bins
    
    Overall coverage:
        weighted_mean(axis_coverages) using coverage_axes weights
        or simple mean if weights not specified
    
    Args:
        history: List of attempt dicts with "axes_binned" key
        config: Config dict with metrics.multi_axis.coverage_axes
    
    Returns:
        Coverage score in [0.0, 1.0]
    """
    if not history:
        return 0.0
    
    multi_axis = config.get("metrics", {}).get("multi_axis", {})
    coverage_axes_specs = multi_axis.get("coverage_axes", [])
    
    # Extract axis names and weights
    axis_names: list[str] = []
    axis_weights: list[float] = []
    
    for spec in coverage_axes_specs:
        name = spec.get("name")
        if name:
            axis_names.append(name)
            axis_weights.append(float(spec.get("weight", 1.0)))
    
    # Fallback: use collapse_axes if coverage_axes not specified
    if not axis_names:
        axis_names = list(multi_axis.get("collapse_axes", ["tb", "lam1_bin"]))
        axis_weights = [1.0] * len(axis_names)
    
    if not axis_names:
        logger.warning("No axes configured for coverage, returning 0.0")
        return 0.0
    
    # Collect visited bins per axis
    visited_per_axis: dict[str, set[Any]] = {name: set() for name in axis_names}
    
    for attempt in history:
        axes_binned = attempt.get("axes_binned", {})
        for axis_name in axis_names:
            value = axes_binned.get(axis_name)
            if value is not None:
                visited_per_axis[axis_name].add(value)
    
    # Compute coverage per axis
    axis_coverages: list[float] = []
    valid_weights: list[float] = []
    
    for axis_name, weight in zip(axis_names, axis_weights):
        domain_size = _get_axis_domain_size(axis_name, config)
        visited_count = len(visited_per_axis[axis_name])
        
        if domain_size > 0:
            coverage = visited_count / domain_size
            axis_coverages.append(coverage)
            valid_weights.append(weight)
    
    if not axis_coverages:
        return 0.0
    
    # Weighted mean
    total_weight = sum(valid_weights)
    if total_weight == 0:
        return 0.0
    
    weighted_sum = sum(cov * w for cov, w in zip(axis_coverages, valid_weights))
    return weighted_sum / total_weight


def _compute_entropy(distribution: Counter[Any]) -> float:
    """
    Compute Shannon entropy of a distribution.
    
    H = -sum(p * log2(p)) for all p > 0
    
    Args:
        distribution: Counter of observed values
    
    Returns:
        Entropy in bits (0.0 for single value, higher for more uniform)
    """
    if not distribution:
        return 0.0
    
    total = sum(distribution.values())
    if total == 0:
        return 0.0
    
    entropy = 0.0
    for count in distribution.values():
        if count > 0:
            p = count / total
            entropy -= p * math.log2(p)
    
    return entropy


def compute_diversity(
    history: Sequence[Mapping[str, Any]], 
    config: Mapping[str, Any]
) -> float:
    """
    Compute Shannon entropy over per-axis distributions with weighted mean.
    
    Diversity formula per axis:
        -sum(p * log2(p)) where p = freq(value) / total_attempts
    
    Overall diversity:
        weighted_mean(axis_entropies) using diversity_axes weights
        or simple mean if weights not specified
    
    Also supports pairwise interactions if configured:
        diversity_pairs: [{"axes": ["tb", "lam1_bin"], "weight": 0.5}, ...]
    
    Args:
        history: List of attempt dicts with "axes_binned" key
        config: Config dict with metrics.multi_axis.diversity_axes
    
    Returns:
        Diversity score (unbounded, higher = more diverse)
    """
    if not history:
        return 0.0
    
    multi_axis = config.get("metrics", {}).get("multi_axis", {})
    diversity_axes_specs = multi_axis.get("diversity_axes", [])
    diversity_pairs = multi_axis.get("diversity_pairs", [])
    
    # Extract axis names and weights
    axis_names: list[str] = []
    axis_weights: list[float] = []
    
    for spec in diversity_axes_specs:
        name = spec.get("name")
        if name:
            axis_names.append(name)
            axis_weights.append(float(spec.get("weight", 1.0)))
    
    # Fallback: use collapse_axes if diversity_axes not specified
    if not axis_names:
        axis_names = list(multi_axis.get("collapse_axes", ["tb", "lam1_bin"]))
        axis_weights = [1.0] * len(axis_names)
    
    if not axis_names and not diversity_pairs:
        logger.warning("No axes configured for diversity, returning 0.0")
        return 0.0
    
    # Collect distributions per axis
    dist_per_axis: dict[str, Counter[Any]] = {name: Counter() for name in axis_names}
    
    for attempt in history:
        axes_binned = attempt.get("axes_binned", {})
        for axis_name in axis_names:
            value = axes_binned.get(axis_name)
            if value is not None:
                dist_per_axis[axis_name][value] += 1
    
    # Compute entropy per axis
    entropies: list[float] = []
    weights: list[float] = []
    
    for axis_name, weight in zip(axis_names, axis_weights):
        entropy = _compute_entropy(dist_per_axis[axis_name])
        entropies.append(entropy)
        weights.append(weight)
    
    # Compute pairwise interaction entropies
    for pair_spec in diversity_pairs:
        pair_axes = pair_spec.get("axes", [])
        pair_weight = float(pair_spec.get("weight", 1.0))
        
        if len(pair_axes) != 2:
            logger.warning(f"Skipping invalid diversity_pair: {pair_spec}")
            continue
        
        axis1, axis2 = pair_axes
        pair_dist: Counter[tuple[Any, Any]] = Counter()
        
        for attempt in history:
            axes_binned = attempt.get("axes_binned", {})
            val1 = axes_binned.get(axis1)
            val2 = axes_binned.get(axis2)
            if val1 is not None and val2 is not None:
                pair_dist[(val1, val2)] += 1
        
        pair_entropy = _compute_entropy(pair_dist)
        entropies.append(pair_entropy)
        weights.append(pair_weight)
    
    if not entropies:
        return 0.0
    
    # Weighted mean
    total_weight = sum(weights)
    if total_weight == 0:
        return 0.0
    
    weighted_sum = sum(e * w for e, w in zip(entropies, weights))
    return weighted_sum / total_weight


def compute_composite(
    yield_val: float,
    coverage: float,
    diversity: float,
    config: Mapping[str, Any]
) -> float:
    """
    Compute weighted composite score from yield, coverage, and diversity.
    
    Formula:
        w_yield * yield + w_coverage * coverage + w_diversity * diversity
    
    Args:
        yield_val: Yield metric (typically in [0.0, 1.0])
        coverage: Coverage metric (in [0.0, 1.0])
        diversity: Diversity metric (unbounded, typically [0.0, 10.0])
        config: Config dict with metrics.weights
    
    Returns:
        Composite score (unbounded, depends on weights and diversity scale)
    """
    weights = config.get("metrics", {}).get("weights", {})
    
    w_yield = float(weights.get("yield", 1.0))
    w_coverage = float(weights.get("coverage", 1.0))
    w_diversity = float(weights.get("diversity", 1.0))
    
    composite = (
        w_yield * yield_val +
        w_coverage * coverage +
        w_diversity * diversity
    )
    
    return composite
