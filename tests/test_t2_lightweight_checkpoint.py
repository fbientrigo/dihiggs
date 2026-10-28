import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = (Path(__file__).resolve().parents[1] /
          "campaigns/t2_t4_20261008T043858Z/scripts/run_t2_lightweight.py")
SPEC = importlib.util.spec_from_file_location("run_t2_lightweight", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_checkpoint_replay_and_conflict_detection(tmp_path):
    checkpoint = tmp_path / "points.jsonl"
    assert module.generate_points(checkpoint) == 300
    original = checkpoint.read_text()
    assert module.generate_points(checkpoint) == 300
    assert checkpoint.read_text() == original

    records = [json.loads(line) for line in original.splitlines()]
    assert all("source_sha256" in row and "grid_sha256" in row for row in records)
    assert all({f"lambda{i}" for i in range(1, 8)} <= row.keys() for row in records)
    assert all({f"Z{i}" for i in range(1, 8)} <= row.keys() for row in records)

    records[0]["lambda7"] = "corrupted"
    checkpoint.write_text("\n".join(json.dumps(row) for row in records) + "\n")
    with pytest.raises(ValueError, match="conflicts"):
        module.generate_points(checkpoint)


def test_checkpoint_rejects_invalid_json(tmp_path):
    checkpoint = tmp_path / "points.jsonl"
    checkpoint.write_text("{truncated\n")
    with pytest.raises(ValueError, match="invalid checkpoint"):
        module.generate_points(checkpoint)
