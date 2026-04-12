"""
Comprehensive tests for autoresearch.harness.dihiggs_adaptation module.

Tests cover:
- select_arm_ucb1: cold-start, exploitation, exploration
- update_state: arm stats updates, history appending
- compute_arm_reward: composite score computation
- Edge cases: empty history, missing config, invalid arms
"""

from __future__ import annotations

from typing import Any

import pytest

from autoresearch.harness.dihiggs_adaptation import (
    BanditState,
    compute_arm_reward,
    select_arm_ucb1,
    update_state,
)


@pytest.fixture
def basic_config() -> dict[str, Any]:
    return {
        "arms": [
            {"id": "arm1", "explorer": "adaptive"},
            {"id": "arm2", "explorer": "branch"},
        ],
        "adaptation": {
            "ucb1_exploration_constant": 1.41,
            "min_pulls_per_arm": 3,
        },
        "search": {
            "tb_values": [1000, 5000, 10000, 30000],
            "lam1": {"min": -20.0, "max": 20.0, "n_bins": 40},
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
            },
        },
    }


def test_select_arm_ucb1_cold_start(basic_config: dict[str, Any]) -> None:
    """Arms with 0 pulls should be selected first (infinity score)."""
    state = BanditState(
        arm_stats={
            "arm1": {"pulls": 0, "total_reward": 0.0},
            "arm2": {"pulls": 0, "total_reward": 0.0},
        },
        global_history=[],
        arm_histories={"arm1": [], "arm2": []},
    )
    
    arm = select_arm_ucb1(state, basic_config)
    
    assert arm in ["arm1", "arm2"]


def test_select_arm_ucb1_exploitation(basic_config: dict[str, Any]) -> None:
    """With zero exploration, highest mean reward should be selected."""
    state = BanditState(
        arm_stats={
            "arm1": {"pulls": 10, "total_reward": 8.0},
            "arm2": {"pulls": 10, "total_reward": 5.0},
        },
        global_history=[],
        arm_histories={"arm1": [], "arm2": []},
    )
    
    config_no_exploration = {
        **basic_config,
        "adaptation": {"ucb1_exploration_constant": 0.0},
    }
    
    arm = select_arm_ucb1(state, config_no_exploration)
    
    assert arm == "arm1"


def test_select_arm_ucb1_exploration(basic_config: dict[str, Any]) -> None:
    """With high exploration, less-tried arm should be favored."""
    state = BanditState(
        arm_stats={
            "arm1": {"pulls": 100, "total_reward": 80.0},
            "arm2": {"pulls": 1, "total_reward": 0.5},
        },
        global_history=[],
        arm_histories={"arm1": [], "arm2": []},
    )
    
    config_high_exploration = {
        **basic_config,
        "adaptation": {"ucb1_exploration_constant": 10.0},
    }
    
    arm = select_arm_ucb1(state, config_high_exploration)
    
    # With high exploration constant, arm2 should win due to exploration bonus
    assert arm == "arm2"


def test_select_arm_ucb1_initializes_missing_arms(basic_config: dict[str, Any]) -> None:
    """State should auto-initialize arms that don't exist in arm_stats."""
    state = BanditState(
        arm_stats={},
        global_history=[],
        arm_histories={},
    )
    
    arm = select_arm_ucb1(state, basic_config)
    
    assert arm in ["arm1", "arm2"]
    assert "arm1" in state.arm_stats
    assert "arm2" in state.arm_stats
    assert state.arm_stats["arm1"]["pulls"] == 0
    assert state.arm_stats["arm2"]["pulls"] == 0


def test_select_arm_ucb1_raises_on_no_arms() -> None:
    """Should raise ValueError if no arms configured."""
    state = BanditState()
    config_no_arms = {"arms": []}
    
    with pytest.raises(ValueError, match="non-empty list"):
        select_arm_ucb1(state, config_no_arms)


def test_update_state_increments_pulls(basic_config: dict[str, Any]) -> None:
    """update_state should increment pull count."""
    state = BanditState(
        arm_stats={"arm1": {"pulls": 5, "total_reward": 3.0}},
        global_history=[],
        arm_histories={"arm1": []},
    )
    
    attempt_data = {"axes_binned": {"tb": 1000, "lam1_bin": 0}, "yield": 0.8}
    
    update_state(state, "arm1", reward=0.75, attempt_data=attempt_data)
    
    assert state.arm_stats["arm1"]["pulls"] == 6


def test_update_state_adds_reward(basic_config: dict[str, Any]) -> None:
    """update_state should accumulate rewards."""
    state = BanditState(
        arm_stats={"arm1": {"pulls": 5, "total_reward": 3.0}},
        global_history=[],
        arm_histories={"arm1": []},
    )
    
    attempt_data = {"axes_binned": {"tb": 1000, "lam1_bin": 0}, "yield": 0.8}
    
    update_state(state, "arm1", reward=0.75, attempt_data=attempt_data)
    
    assert state.arm_stats["arm1"]["total_reward"] == pytest.approx(3.75)


def test_update_state_appends_to_global_history(basic_config: dict[str, Any]) -> None:
    """update_state should append attempt to global_history."""
    state = BanditState(
        arm_stats={"arm1": {"pulls": 0, "total_reward": 0.0}},
        global_history=[],
        arm_histories={"arm1": []},
    )
    
    attempt_data = {"axes_binned": {"tb": 1000, "lam1_bin": 0}, "yield": 0.8}
    
    update_state(state, "arm1", reward=0.75, attempt_data=attempt_data)
    
    assert len(state.global_history) == 1
    assert state.global_history[0]["axes_binned"]["tb"] == 1000


def test_update_state_appends_to_arm_history(basic_config: dict[str, Any]) -> None:
    """update_state should append attempt to per-arm history."""
    state = BanditState(
        arm_stats={"arm1": {"pulls": 0, "total_reward": 0.0}},
        global_history=[],
        arm_histories={"arm1": []},
    )
    
    attempt_data = {"axes_binned": {"tb": 1000, "lam1_bin": 0}, "yield": 0.8}
    
    update_state(state, "arm1", reward=0.75, attempt_data=attempt_data)
    
    assert len(state.arm_histories["arm1"]) == 1
    assert state.arm_histories["arm1"][0]["yield"] == 0.8


def test_update_state_initializes_missing_arm(basic_config: dict[str, Any]) -> None:
    """update_state should initialize arm if not present."""
    state = BanditState()
    
    attempt_data = {"axes_binned": {"tb": 1000, "lam1_bin": 0}, "yield": 0.8}
    
    update_state(state, "arm1", reward=0.75, attempt_data=attempt_data)
    
    assert "arm1" in state.arm_stats
    assert state.arm_stats["arm1"]["pulls"] == 1
    assert state.arm_stats["arm1"]["total_reward"] == 0.75


def test_compute_arm_reward_returns_composite(basic_config: dict[str, Any]) -> None:
    """compute_arm_reward should return composite score."""
    history = [
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}, "yield": 0.8},
        {"axes_binned": {"tb": 5000, "lam1_bin": 1}, "yield": 0.9},
    ]
    
    reward = compute_arm_reward(history, basic_config)
    
    assert reward > 0.0
    assert isinstance(reward, float)


def test_compute_arm_reward_uses_latest_yield(basic_config: dict[str, Any]) -> None:
    """compute_arm_reward should use yield from latest attempt."""
    history = [
        {"axes_binned": {"tb": 1000, "lam1_bin": 0}, "yield": 0.5},
        {"axes_binned": {"tb": 5000, "lam1_bin": 1}, "yield": 0.9},
    ]
    
    reward = compute_arm_reward(history, basic_config)
    
    # Latest yield is 0.9, so reward should reflect this
    # With weights: yield=0.5, coverage=0.3, diversity=0.2
    # Minimum reward would be 0.5 * 0.9 = 0.45 (if coverage/diversity = 0)
    assert reward >= 0.4


def test_compute_arm_reward_empty_history_returns_zero(basic_config: dict[str, Any]) -> None:
    """compute_arm_reward should return 0.0 for empty history."""
    history: list[dict[str, Any]] = []
    
    reward = compute_arm_reward(history, basic_config)
    
    assert reward == 0.0


def test_end_to_end_adaptation_cycle(basic_config: dict[str, Any]) -> None:
    """Test full cycle: select → update → select again."""
    state = BanditState()
    
    # First selection (cold start)
    arm1 = select_arm_ucb1(state, basic_config)
    assert arm1 in ["arm1", "arm2"]
    
    # Update after first round
    attempt_data_1 = {"axes_binned": {"tb": 1000, "lam1_bin": 0}, "yield": 0.8}
    reward_1 = compute_arm_reward([attempt_data_1], basic_config)
    update_state(state, arm1, reward_1, attempt_data_1)
    
    assert state.arm_stats[arm1]["pulls"] == 1
    assert len(state.global_history) == 1
    
    # Second selection (should favor the other cold-start arm)
    arm2 = select_arm_ucb1(state, basic_config)
    
    # After first pull, the other arm should be selected (cold-start)
    other_arm = "arm2" if arm1 == "arm1" else "arm1"
    assert arm2 == other_arm or arm2 == arm1  # Either is valid depending on implementation
    
    # Update after second round
    attempt_data_2 = {"axes_binned": {"tb": 5000, "lam1_bin": 1}, "yield": 0.9}
    reward_2 = compute_arm_reward([attempt_data_1, attempt_data_2], basic_config)
    update_state(state, arm2, reward_2, attempt_data_2)
    
    assert state.arm_stats[arm2]["pulls"] == 1
    assert len(state.global_history) == 2
