"""T6 (mission Sec. 26): the daemon implementation must never modify the frozen
physics substrate or the pre-existing discovery layer.

`tests/fixtures/frozen_substrate_hashes.json` was captured immediately before
any daemon file was written, over every pre-existing .py file in
search_substrate/ and search_discovery/ (excluding the four new daemon files:
policy_adapter.py, daemon_checkpoint.py, daemon.py, rebuild_archive.py). This
regression test stays in the repo permanently to catch any future accidental
edit to those frozen files.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = Path(__file__).resolve().parent / "fixtures" / "frozen_substrate_hashes.json"


def test_pre_existing_substrate_and_discovery_files_are_byte_identical():
    expected = json.loads(FIXTURE.read_text(encoding="utf-8"))
    mismatches = []
    missing = []
    for relative_path, expected_hash in expected.items():
        path = ROOT / relative_path
        if not path.exists():
            missing.append(relative_path)
            continue
        actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual_hash != expected_hash:
            mismatches.append(relative_path)
    assert not missing, f"frozen files disappeared: {missing}"
    assert not mismatches, f"frozen files were modified since fixture capture: {mismatches}"


def test_new_daemon_files_are_additive_only():
    """The four new daemon files must exist and must not appear in the frozen fixture."""
    expected = json.loads(FIXTURE.read_text(encoding="utf-8"))
    new_files = [
        "search_discovery/policy_adapter.py",
        "search_discovery/daemon_checkpoint.py",
        "search_discovery/daemon.py",
        "search_discovery/rebuild_archive.py",
    ]
    for relative_path in new_files:
        assert relative_path not in expected, f"{relative_path} should not be in the frozen fixture"
        assert (ROOT / relative_path).exists(), f"{relative_path} should exist"
