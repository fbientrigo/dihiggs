"""
Bandit adaptation logic for DiHiggs autoresearch campaigns.

This module implements UCB1 arm selection, composite score rewards, and 
history tracking for multi-armed bandit-based configuration adaptation.

Key components:
- BanditState: Tracks arm statistics and attempt histories
- select_arm_ucb1: UCB1 formula with exploration constant
- update_state: Update arm stats and append to history
- compute_arm_reward: Composite score from yield + coverage + diversity
"""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from autoresearch.benchmarks.score import (
    compute_composite,
    compute_coverage,
    compute_diversity,
    compute_metrics,
)

logger = logging.getLogger(__name__)


@dataclass
class BanditState:
    """
    Bandit state tracking for multi-armed configuration search.
    
    Attributes:
        arm_stats: Per-arm statistics
            - pulls: int (number of times arm was selected)
            - total_reward: float (cumulative reward sum)
        global_history: All attempts across all arms with axes_binned
            Each entry: {"axes_binned": {...}, "yield": float, ...}
        arm_histories: Per-arm attempt histories for debugging/analysis
    """
    arm_stats: dict[str, dict[str, Any]] = field(default_factory=dict)
    global_history: list[dict[str, Any]] = field(default_factory=list)
    arm_histories: dict[str, list[dict[str, Any]]] = field(default_factory=dict)


def select_arm_ucb1(state: BanditState, config: Mapping[str, Any]) -> str:
    """
    Select arm using UCB1 algorithm.
    
    UCB1 formula:
        mean_reward + c * sqrt(log(total_pulls) / arm_pulls)
    
    Cold-start handling:
        Arms with 0 pulls get infinity score (always explored first)
    
    Args:
        state: Current bandit state with arm statistics
        config: Configuration dict with adaptation.ucb1_exploration_constant
    
    Returns:
        arm_id: Selected arm identifier
    
    Raises:
        ValueError: If no arms configured or all arms missing from state
    """
    arms_cfg = config.get("arms", [])
    if not isinstance(arms_cfg, list) or not arms_cfg:
        raise ValueError("config.arms must be a non-empty list")
    
    arm_ids = [
        str(arm.get("id")) 
        for arm in arms_cfg 
        if isinstance(arm, Mapping) and arm.get("id")
    ]
    
    if not arm_ids:
        raise ValueError("No valid arms found in config")
    
    adaptation_cfg = config.get("adaptation", {})
    if not isinstance(adaptation_cfg, Mapping):
        adaptation_cfg = {}
    
    exploration_constant = float(adaptation_cfg.get("ucb1_exploration_constant", 1.41))
    
    # Initialize arm_stats for any missing arms
    for arm_id in arm_ids:
        if arm_id not in state.arm_stats:
            state.arm_stats[arm_id] = {"pulls": 0, "total_reward": 0.0}
        if arm_id not in state.arm_histories:
            state.arm_histories[arm_id] = []
    
    # Compute total pulls across all arms
    total_pulls = sum(
        int(stats.get("pulls", 0)) 
        for stats in state.arm_stats.values()
    )
    
    # UCB1 scores
    best_arm = None
    best_score = -math.inf
    
    for arm_id in arm_ids:
        stats = state.arm_stats[arm_id]
        pulls = int(stats.get("pulls", 0))
        total_reward = float(stats.get("total_reward", 0.0))
        
        # Cold-start: 0 pulls → infinity score
        if pulls == 0:
            ucb1_score = math.inf
        else:
            mean_reward = total_reward / pulls
            if total_pulls > 0:
                exploration_bonus = exploration_constant * math.sqrt(math.log(total_pulls) / pulls)
            else:
                exploration_bonus = 0.0
            ucb1_score = mean_reward + exploration_bonus
        
        logger.debug(
            f"Arm {arm_id}: pulls={pulls}, mean_reward={total_reward/pulls if pulls > 0 else 0.0:.4f}, ucb1={ucb1_score}"
        )
        
        if ucb1_score > best_score:
            best_score = ucb1_score
            best_arm = arm_id
    
    if best_arm is None:
        raise ValueError("Failed to select arm (should not happen)")
    
    logger.info(f"UCB1 selected arm: {best_arm} (score={best_score:.4f})")
    return best_arm


def update_state(
    state: BanditState,
    arm_id: str,
    reward: float,
    attempt_data: Mapping[str, Any]
) -> None:
    """
    Update bandit state after arm execution.
    
    Updates:
        - Increment arm pull count
        - Add reward to arm total
        - Append attempt_data to global_history and arm_histories
    
    Args:
        state: Bandit state to update (modified in-place)
        arm_id: Arm that was executed
        reward: Composite score from compute_arm_reward
        attempt_data: Attempt record with axes_binned, yield, etc.
    """
    if arm_id not in state.arm_stats:
        state.arm_stats[arm_id] = {"pulls": 0, "total_reward": 0.0}
    if arm_id not in state.arm_histories:
        state.arm_histories[arm_id] = []
    
    stats = state.arm_stats[arm_id]
    stats["pulls"] = int(stats.get("pulls", 0)) + 1
    stats["total_reward"] = float(stats.get("total_reward", 0.0)) + reward
    
    # Append to histories
    state.global_history.append(dict(attempt_data))
    state.arm_histories[arm_id].append(dict(attempt_data))
    
    logger.debug(
        f"Updated arm {arm_id}: pulls={stats['pulls']}, "
        f"total_reward={stats['total_reward']:.4f}, "
        f"mean_reward={stats['total_reward']/stats['pulls']:.4f}"
    )


def compute_arm_reward(
    history: Sequence[Mapping[str, Any]],
    config: Mapping[str, Any]
) -> float:
    """
    Compute composite reward for an arm using latest attempt and global history.
    
    Reward components:
        - Yield: from latest attempt's metrics (immediate performance)
        - Coverage: from global_history (exploration breadth)
        - Diversity: from global_history (exploration balance)
    
    Formula:
        composite = w_yield * yield + w_coverage * coverage + w_diversity * diversity
    
    Args:
        history: Per-arm attempt history (used for latest attempt yield)
        config: Configuration dict with metrics.weights
    
    Returns:
        Composite reward score
    
    Note:
        If history is empty, returns 0.0
        Coverage/diversity use global history (passed via config in integration)
    """
    if not history:
        logger.warning("compute_arm_reward called with empty history, returning 0.0")
        return 0.0
    
    # Extract latest attempt for yield
    latest = history[-1]
    yield_val = float(latest.get("yield", 0.0))
    
    # Compute coverage and diversity from history
    coverage = compute_coverage(history, config)
    diversity = compute_diversity(history, config)
    
    composite = compute_composite(yield_val, coverage, diversity, config)
    
    logger.debug(
        f"Reward computed: yield={yield_val:.4f}, coverage={coverage:.4f}, "
        f"diversity={diversity:.4f}, composite={composite:.4f}"
    )
    
    return composite
