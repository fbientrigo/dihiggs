import json
import os
import tempfile
import time
from pathlib import Path
from unittest.mock import patch

import pytest

from autoresearch.harness.dihiggs_parsers import (
    BranchDiscovery,
    parse_branch_checkpoint,
)


@pytest.fixture
def config():
    return {
        "search": {
            "lam1": {
                "min": -20.0,
                "max": 20.0,
                "n_bins": 40,
            }
        }
    }


@pytest.fixture
def temp_checkpoint_dir():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


def make_file_stable(file_path: Path) -> None:
    old_time = time.time() - 10.0
    os.utime(file_path, (old_time, old_time))


def test_parse_valid_checkpoint_returns_discoveries(temp_checkpoint_dir, config):
    track_0001 = temp_checkpoint_dir / "track_0001"
    track_0001.mkdir()
    events_jsonl = track_0001 / "events.jsonl"

    events = [
        {
            "event_type": "ATTEMPT_COMPLETED",
            "utc": "2026-04-01T10:00:00Z",
            "payload": {
                "result": {"run_dir": "/runs/run_001"},
                "params": {"tb": 10000, "lam1": -12.5},
                "step_label": "seed",
            },
        },
        {
            "event_type": "ATTEMPT_COMPLETED",
            "utc": "2026-04-01T10:05:00Z",
            "payload": {
                "result": {"run_dir": "/runs/run_002"},
                "params": {"tb": 10000, "lam1": 5.0},
                "step_label": "expand",
            },
        },
    ]

    with open(events_jsonl, "w", encoding="utf-8") as f:
        for event in events:
            f.write(json.dumps(event) + "\n")

    # Ensure file is stable (modified >2s ago)
    old_time = time.time() - 10.0
    os.utime(events_jsonl, (old_time, old_time))

    discoveries = parse_branch_checkpoint(str(temp_checkpoint_dir), config, set())

    assert len(discoveries) == 2
    
    # Sort by (track_id, step_name), so "expand" comes before "seed" alphabetically
    assert discoveries[0].run_dir == "/runs/run_002"
    assert discoveries[0].tb == 10000
    assert discoveries[0].lam1_raw == 5.0
    assert discoveries[0].track_id == "track_0001"
    assert discoveries[0].step_name == "expand"
    assert discoveries[0].lam1_bin == 25

    assert discoveries[1].run_dir == "/runs/run_001"
    assert discoveries[1].tb == 10000
    assert discoveries[1].lam1_raw == -12.5
    assert discoveries[1].track_id == "track_0001"
    assert discoveries[1].step_name == "seed"
    assert discoveries[1].lam1_bin == 7


def test_parse_skips_known_run_dirs(temp_checkpoint_dir, config):
    track_0001 = temp_checkpoint_dir / "track_0001"
    track_0001.mkdir()
    events_jsonl = track_0001 / "events.jsonl"

    events = [
        {
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {"run_dir": "/runs/run_001"},
                "params": {"tb": 10000, "lam1": -12.5},
                "step_label": "seed",
            },
        },
        {
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {"run_dir": "/runs/run_002"},
                "params": {"tb": 10000, "lam1": 5.0},
                "step_label": "expand",
            },
        },
    ]

    with open(events_jsonl, "w", encoding="utf-8") as f:
        for event in events:
            f.write(json.dumps(event) + "\n")

    make_file_stable(events_jsonl)

    known = {"/runs/run_001"}
    discoveries = parse_branch_checkpoint(str(temp_checkpoint_dir), config, known)

    assert len(discoveries) == 1
    assert discoveries[0].run_dir == "/runs/run_002"


def test_parse_missing_checkpoint_file(temp_checkpoint_dir, config):
    track_0001 = temp_checkpoint_dir / "track_0001"
    track_0001.mkdir()

    discoveries = parse_branch_checkpoint(str(temp_checkpoint_dir), config, set())

    assert len(discoveries) == 0


def test_parse_missing_checkpoint_root(config):
    nonexistent = "/tmp/nonexistent_checkpoint_dir_12345"
    discoveries = parse_branch_checkpoint(nonexistent, config, set())

    assert len(discoveries) == 0


def test_parse_malformed_jsonl_lines(temp_checkpoint_dir, config):
    track_0001 = temp_checkpoint_dir / "track_0001"
    track_0001.mkdir()
    events_jsonl = track_0001 / "events.jsonl"

    with open(events_jsonl, "w", encoding="utf-8") as f:
        f.write('{"event_type": "ATTEMPT_COMPLETED", "payload": {"result": {"run_dir": "/runs/run_001"}, "params": {"tb": 10000, "lam1": -12.5}, "step_label": "seed"}}\n')
        f.write("THIS IS NOT JSON\n")
        f.write('{"event_type": "ATTEMPT_COMPLETED", "payload": {"result": {"run_dir": "/runs/run_002"}, "params": {"tb": 10000, "lam1": 5.0}, "step_label": "expand"}}\n')

    make_file_stable(events_jsonl)

    discoveries = parse_branch_checkpoint(str(temp_checkpoint_dir), config, set())

    assert len(discoveries) == 2
    assert discoveries[0].run_dir == "/runs/run_002"
    assert discoveries[1].run_dir == "/runs/run_001"


def test_parse_filters_by_event_type(temp_checkpoint_dir, config):
    track_0001 = temp_checkpoint_dir / "track_0001"
    track_0001.mkdir()
    events_jsonl = track_0001 / "events.jsonl"

    events = [
        {
            "event_type": "ATTEMPT_STARTED",
            "payload": {
                "result": {"run_dir": "/runs/run_001"},
                "params": {"tb": 10000, "lam1": -12.5},
                "step_label": "seed",
            },
        },
        {
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {"run_dir": "/runs/run_002"},
                "params": {"tb": 10000, "lam1": 5.0},
                "step_label": "expand",
            },
        },
        {
            "event_type": "ATTEMPT_FAILED",
            "payload": {
                "result": {"run_dir": "/runs/run_003"},
                "params": {"tb": 10000, "lam1": 10.0},
                "step_label": "refine",
            },
        },
    ]

    with open(events_jsonl, "w", encoding="utf-8") as f:
        for event in events:
            f.write(json.dumps(event) + "\n")

    make_file_stable(events_jsonl)

    discoveries = parse_branch_checkpoint(str(temp_checkpoint_dir), config, set())

    assert len(discoveries) == 1
    assert discoveries[0].run_dir == "/runs/run_002"


def test_parse_multiple_tracks_sorted(temp_checkpoint_dir, config):
    track_0002 = temp_checkpoint_dir / "track_0002"
    track_0002.mkdir()
    events_2 = track_0002 / "events.jsonl"
    with open(events_2, "w", encoding="utf-8") as f:
        f.write(json.dumps({
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {"run_dir": "/runs/run_002"},
                "params": {"tb": 20000, "lam1": 0.0},
                "step_label": "zzz_last",
            },
        }) + "\n")
    make_file_stable(events_2)

    track_0001 = temp_checkpoint_dir / "track_0001"
    track_0001.mkdir()
    events_1 = track_0001 / "events.jsonl"
    with open(events_1, "w", encoding="utf-8") as f:
        f.write(json.dumps({
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {"run_dir": "/runs/run_001"},
                "params": {"tb": 10000, "lam1": -5.0},
                "step_label": "aaa_first",
            },
        }) + "\n")
    make_file_stable(events_1)

    discoveries = parse_branch_checkpoint(str(temp_checkpoint_dir), config, set())

    assert len(discoveries) == 2
    assert discoveries[0].track_id == "track_0001"
    assert discoveries[1].track_id == "track_0002"


def test_parse_empty_checkpoint_returns_empty_list(temp_checkpoint_dir, config):
    track_0001 = temp_checkpoint_dir / "track_0001"
    track_0001.mkdir()
    events_jsonl = track_0001 / "events.jsonl"
    events_jsonl.touch()

    make_file_stable(events_jsonl)

    discoveries = parse_branch_checkpoint(str(temp_checkpoint_dir), config, set())

    assert len(discoveries) == 0


def test_parse_skips_null_run_dir(temp_checkpoint_dir, config):
    track_0001 = temp_checkpoint_dir / "track_0001"
    track_0001.mkdir()
    events_jsonl = track_0001 / "events.jsonl"

    events = [
        {
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {"run_dir": None},
                "params": {"tb": 10000, "lam1": -12.5},
                "step_label": "seed",
            },
        },
        {
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {"run_dir": ""},
                "params": {"tb": 10000, "lam1": 5.0},
                "step_label": "expand",
            },
        },
        {
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {"run_dir": "/runs/run_003"},
                "params": {"tb": 10000, "lam1": 10.0},
                "step_label": "refine",
            },
        },
    ]

    with open(events_jsonl, "w", encoding="utf-8") as f:
        for event in events:
            f.write(json.dumps(event) + "\n")

    make_file_stable(events_jsonl)

    discoveries = parse_branch_checkpoint(str(temp_checkpoint_dir), config, set())

    assert len(discoveries) == 1
    assert discoveries[0].run_dir == "/runs/run_003"


def test_parse_handles_missing_keys(temp_checkpoint_dir, config):
    track_0001 = temp_checkpoint_dir / "track_0001"
    track_0001.mkdir()
    events_jsonl = track_0001 / "events.jsonl"

    events = [
        {
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {},
                "params": {"tb": 10000, "lam1": -12.5},
                "step_label": "seed",
            },
        },
        {
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {"run_dir": "/runs/run_002"},
                "params": {},
                "step_label": "expand",
            },
        },
        {
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {"run_dir": "/runs/run_003"},
                "params": {"tb": 10000, "lam1": 10.0},
            },
        },
        {
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {"run_dir": "/runs/run_004"},
                "params": {"tb": 10000, "lam1": 15.0},
                "step_label": "good",
            },
        },
    ]

    with open(events_jsonl, "w", encoding="utf-8") as f:
        for event in events:
            f.write(json.dumps(event) + "\n")

    make_file_stable(events_jsonl)

    discoveries = parse_branch_checkpoint(str(temp_checkpoint_dir), config, set())

    assert len(discoveries) == 1
    assert discoveries[0].run_dir == "/runs/run_004"


def test_parse_handles_invalid_payload_types(temp_checkpoint_dir, config):
    track_0001 = temp_checkpoint_dir / "track_0001"
    track_0001.mkdir()
    events_jsonl = track_0001 / "events.jsonl"

    with open(events_jsonl, "w", encoding="utf-8") as f:
        f.write(json.dumps({
            "event_type": "ATTEMPT_COMPLETED",
            "payload": "not a dict",
        }) + "\n")
        f.write(json.dumps({
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": "not a dict",
                "params": {"tb": 10000, "lam1": -12.5},
                "step_label": "seed",
            },
        }) + "\n")
        f.write(json.dumps({
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {"run_dir": "/runs/run_003"},
                "params": "not a dict",
                "step_label": "expand",
            },
        }) + "\n")
        f.write(json.dumps({
            "event_type": "ATTEMPT_COMPLETED",
            "payload": {
                "result": {"run_dir": "/runs/run_004"},
                "params": {"tb": 10000, "lam1": 10.0},
                "step_label": "good",
            },
        }) + "\n")

    make_file_stable(events_jsonl)

    discoveries = parse_branch_checkpoint(str(temp_checkpoint_dir), config, set())

    assert len(discoveries) == 1
    assert discoveries[0].run_dir == "/runs/run_004"
