from __future__ import annotations

import json
from pathlib import Path

import pytest

from search_discovery import daemon_checkpoint as ckpt


def test_fresh_checkpoint_has_expected_shape():
    state = ckpt.fresh("digest123")
    assert state["schema_version"] == ckpt.SCHEMA_VERSION
    assert state["cycle"] == 0
    assert state["total_evaluations"] == 0
    assert state["config_digest"] == "digest123"
    assert state["stopped"] is False


def test_load_returns_none_when_absent(tmp_path):
    assert ckpt.load(tmp_path) is None


def test_save_then_load_round_trips(tmp_path):
    state = ckpt.fresh("digestABC")
    state["cycle"] = 3
    state["total_evaluations"] = 42
    saved = ckpt.save(tmp_path, state)
    assert saved["last_checkpoint_utc"] > 0
    loaded = ckpt.load(tmp_path)
    assert loaded["cycle"] == 3
    assert loaded["total_evaluations"] == 42
    assert loaded["config_digest"] == "digestABC"


def test_save_is_atomic_no_tmp_file_left_behind(tmp_path):
    state = ckpt.fresh("d")
    ckpt.save(tmp_path, state)
    assert ckpt.path_for(tmp_path).exists()
    assert not ckpt.path_for(tmp_path).with_suffix(".json.tmp").exists()


def test_save_overwrites_cleanly_and_prior_content_never_partially_visible(tmp_path):
    state = ckpt.fresh("d")
    for cycle in range(5):
        state["cycle"] = cycle
        ckpt.save(tmp_path, state)
        loaded = json.loads(ckpt.path_for(tmp_path).read_text(encoding="utf-8"))
        assert loaded["cycle"] == cycle  # never a half-written mix of old/new
