"""
Tests for adaptive checkpoint parser.

Verifies parse_adaptive_checkpoint correctly extracts discoveries from checkpoint
files, handles edge cases, and filters known run_dirs.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import pytest

from autoresearch.harness.dihiggs_parsers import parse_adaptive_checkpoint


@pytest.fixture
def sample_config():
    """Config with lam1 binning parameters."""
    return {
        "search": {
            "lam1": {
                "min": -20.0,
                "max": 20.0,
                "n_bins": 10,
            }
        }
    }


@pytest.fixture
def temp_checkpoint_root(tmp_path):
    """Create temporary checkpoint directory structure."""
    return tmp_path / "checkpoints"


def create_checkpoint(
    checkpoint_root: Path,
    iter_index: int,
    proposals: list[dict],
    metadata: dict | None = None,
) -> Path:
    """Helper to create checkpoint file with given proposals."""
    iter_dir = checkpoint_root / f"iter_{iter_index:04d}"
    iter_dir.mkdir(parents=True, exist_ok=True)

    state = {
        "proposals": proposals,
        "metadata": metadata or {"iter_index": iter_index},
    }

    state_file = iter_dir / "adaptive_state.json"
    with open(state_file, "w") as f:
        json.dump(state, f)

    return state_file


def test_parse_valid_checkpoint(temp_checkpoint_root, sample_config):
    """Parse checkpoint with valid proposals returns discoveries."""
    proposals = [
        {
            "run_dir": "/path/to/run1",
            "bin_index": 5,
            "lam1_min": -10.0,
        },
        {
            "run_dir": "/path/to/run2",
            "bin_index": 8,
            "lam1_min": 5.0,
        },
    ]

    state_file = create_checkpoint(temp_checkpoint_root, 1, proposals)
    
    # Ensure file is stable (modified >2s ago)
    old_time = time.time() - 10.0
    os.utime(state_file, (old_time, old_time))

    discoveries = parse_adaptive_checkpoint(
        str(temp_checkpoint_root), sample_config, set()
    )

    assert len(discoveries) == 2

    assert discoveries[0].run_dir == "/path/to/run1"
    assert discoveries[0].tb == 5
    assert discoveries[0].lam1_raw == -10.0
    assert discoveries[0].lam1_bin == 2
    assert discoveries[0].iter_index == 1
    assert discoveries[0].proposal_index == 5

    assert discoveries[1].run_dir == "/path/to/run2"
    assert discoveries[1].tb == 8
    assert discoveries[1].lam1_raw == 5.0
    assert discoveries[1].lam1_bin == 6
    assert discoveries[1].iter_index == 1
    assert discoveries[1].proposal_index == 8


def test_skip_known_run_dirs(temp_checkpoint_root, sample_config):
    """Filter out already-known run_dirs."""
    proposals = [
        {"run_dir": "/path/to/run1", "bin_index": 5, "lam1_min": -10.0},
        {"run_dir": "/path/to/run2", "bin_index": 8, "lam1_min": 5.0},
    ]

    state_file = create_checkpoint(temp_checkpoint_root, 1, proposals)
    old_time = time.time() - 10.0
    os.utime(state_file, (old_time, old_time))

    known_run_dirs = {"/path/to/run1"}

    discoveries = parse_adaptive_checkpoint(
        str(temp_checkpoint_root), sample_config, known_run_dirs
    )

    assert len(discoveries) == 1
    assert discoveries[0].run_dir == "/path/to/run2"


def test_skip_unstable_checkpoint(temp_checkpoint_root, sample_config):
    """Skip checkpoint modified within last 2 seconds."""
    proposals = [
        {"run_dir": "/path/to/run1", "bin_index": 5, "lam1_min": -10.0},
    ]

    create_checkpoint(temp_checkpoint_root, 1, proposals)
    # File is fresh (just created), should be skipped

    discoveries = parse_adaptive_checkpoint(
        str(temp_checkpoint_root), sample_config, set()
    )

    assert len(discoveries) == 0


def test_handle_missing_checkpoint_root(sample_config):
    """Missing checkpoint root returns empty list with warning."""
    discoveries = parse_adaptive_checkpoint(
        "/nonexistent/path", sample_config, set()
    )

    assert discoveries == []


def test_handle_missing_checkpoint_file(temp_checkpoint_root, sample_config):
    """Missing checkpoint file is skipped gracefully."""
    # Create iter dir but no state file
    iter_dir = temp_checkpoint_root / "iter_0001"
    iter_dir.mkdir(parents=True)

    discoveries = parse_adaptive_checkpoint(
        str(temp_checkpoint_root), sample_config, set()
    )

    assert discoveries == []


def test_handle_malformed_json(temp_checkpoint_root, sample_config):
    """Malformed JSON checkpoint is skipped with warning."""
    iter_dir = temp_checkpoint_root / "iter_0001"
    iter_dir.mkdir(parents=True)

    state_file = iter_dir / "adaptive_state.json"
    with open(state_file, "w") as f:
        f.write("{invalid json")

    old_time = time.time() - 10.0
    os.utime(state_file, (old_time, old_time))

    discoveries = parse_adaptive_checkpoint(
        str(temp_checkpoint_root), sample_config, set()
    )

    assert discoveries == []


def test_handle_non_dict_json(temp_checkpoint_root, sample_config):
    """Non-dict JSON checkpoint is skipped with warning."""
    iter_dir = temp_checkpoint_root / "iter_0001"
    iter_dir.mkdir(parents=True)

    state_file = iter_dir / "adaptive_state.json"
    with open(state_file, "w") as f:
        json.dump(["not", "a", "dict"], f)

    old_time = time.time() - 10.0
    os.utime(state_file, (old_time, old_time))

    discoveries = parse_adaptive_checkpoint(
        str(temp_checkpoint_root), sample_config, set()
    )

    assert discoveries == []


def test_handle_missing_proposals_key(temp_checkpoint_root, sample_config):
    """Checkpoint without proposals key is handled gracefully."""
    iter_dir = temp_checkpoint_root / "iter_0001"
    iter_dir.mkdir(parents=True)

    state_file = iter_dir / "adaptive_state.json"
    with open(state_file, "w") as f:
        json.dump({"metadata": {"iter_index": 1}}, f)

    old_time = time.time() - 10.0
    os.utime(state_file, (old_time, old_time))

    discoveries = parse_adaptive_checkpoint(
        str(temp_checkpoint_root), sample_config, set()
    )

    assert discoveries == []


def test_skip_null_run_dir(temp_checkpoint_root, sample_config):
    """Proposals with null or empty run_dir are skipped."""
    proposals = [
        {"run_dir": None, "bin_index": 5, "lam1_min": -10.0},
        {"run_dir": "", "bin_index": 8, "lam1_min": 5.0},
        {"run_dir": "/path/to/run3", "bin_index": 3, "lam1_min": 0.0},
    ]

    state_file = create_checkpoint(temp_checkpoint_root, 1, proposals)
    old_time = time.time() - 10.0
    os.utime(state_file, (old_time, old_time))

    discoveries = parse_adaptive_checkpoint(
        str(temp_checkpoint_root), sample_config, set()
    )

    assert len(discoveries) == 1
    assert discoveries[0].run_dir == "/path/to/run3"


def test_handle_missing_optional_fields(temp_checkpoint_root, sample_config):
    """Proposals with missing optional fields use fallback defaults."""
    proposals = [
        {"run_dir": "/path/to/run1"},
        {"run_dir": "/path/to/run2", "bin_index": 5, "lam1_min": -10.0},
    ]

    state_file = create_checkpoint(temp_checkpoint_root, 1, proposals)
    old_time = time.time() - 10.0
    os.utime(state_file, (old_time, old_time))

    discoveries = parse_adaptive_checkpoint(
        str(temp_checkpoint_root), sample_config, set()
    )

    assert len(discoveries) == 2
    assert discoveries[0].run_dir == "/path/to/run1"
    assert discoveries[0].tb == 0
    assert discoveries[0].lam1_raw == 0.0
    assert discoveries[1].run_dir == "/path/to/run2"


def test_skip_invalid_field_types(temp_checkpoint_root, sample_config):
    """Proposals with invalid field types are skipped."""
    proposals = [
        {"run_dir": "/path/to/run1", "bin_index": "not-a-number", "lam1_min": -10.0},
        {"run_dir": "/path/to/run2", "bin_index": 5, "lam1_min": -10.0},
    ]

    state_file = create_checkpoint(temp_checkpoint_root, 1, proposals)
    old_time = time.time() - 10.0
    os.utime(state_file, (old_time, old_time))

    discoveries = parse_adaptive_checkpoint(
        str(temp_checkpoint_root), sample_config, set()
    )

    assert len(discoveries) == 1
    assert discoveries[0].run_dir == "/path/to/run2"


def test_multiple_iterations_sorted(temp_checkpoint_root, sample_config):
    """Multiple iterations are parsed and sorted correctly."""
    proposals_iter1 = [
        {"run_dir": "/path/to/run1", "bin_index": 5, "lam1_min": -10.0},
    ]
    proposals_iter2 = [
        {"run_dir": "/path/to/run2", "bin_index": 3, "lam1_min": 0.0},
    ]

    file1 = create_checkpoint(temp_checkpoint_root, 1, proposals_iter1)
    file2 = create_checkpoint(temp_checkpoint_root, 2, proposals_iter2)

    old_time = time.time() - 10.0
    os.utime(file1, (old_time, old_time))
    os.utime(file2, (old_time, old_time))

    discoveries = parse_adaptive_checkpoint(
        str(temp_checkpoint_root), sample_config, set()
    )

    assert len(discoveries) == 2
    assert discoveries[0].iter_index == 1
    assert discoveries[1].iter_index == 2


def test_empty_checkpoint_returns_empty(temp_checkpoint_root, sample_config):
    """Checkpoint with empty proposals list returns empty discoveries."""
    state_file = create_checkpoint(temp_checkpoint_root, 1, [])
    old_time = time.time() - 10.0
    os.utime(state_file, (old_time, old_time))

    discoveries = parse_adaptive_checkpoint(
        str(temp_checkpoint_root), sample_config, set()
    )

    assert discoveries == []


def test_lam1_binning_integration(temp_checkpoint_root, sample_config):
    """Verify lam1 binning uses axis contract correctly."""
    proposals = [
        {"run_dir": "/path/to/run1", "bin_index": 0, "lam1_min": -20.0},  # bin 0
        {"run_dir": "/path/to/run2", "bin_index": 0, "lam1_min": -18.0},  # bin 0
        {"run_dir": "/path/to/run3", "bin_index": 0, "lam1_min": -10.0},  # bin 2
        {"run_dir": "/path/to/run4", "bin_index": 0, "lam1_min": 0.0},    # bin 5
        {"run_dir": "/path/to/run5", "bin_index": 0, "lam1_min": 19.9},   # bin 9
    ]

    state_file = create_checkpoint(temp_checkpoint_root, 1, proposals)
    old_time = time.time() - 10.0
    os.utime(state_file, (old_time, old_time))

    discoveries = parse_adaptive_checkpoint(
        str(temp_checkpoint_root), sample_config, set()
    )

    assert len(discoveries) == 5
    assert discoveries[0].lam1_bin == 0
    assert discoveries[1].lam1_bin == 0
    assert discoveries[2].lam1_bin == 2
    assert discoveries[3].lam1_bin == 5
    assert discoveries[4].lam1_bin == 9
