"""Unit tests for the torn-ledger-tail repair used on resume after a hard kill.
No physics evaluated here."""
from __future__ import annotations

import json

from search_discovery import daemon


def test_repair_is_noop_on_a_clean_ledger(tmp_path):
    path = tmp_path / "ledger.jsonl"
    path.write_text('{"a": 1}\n{"a": 2}\n', encoding="utf-8")
    changed = daemon.repair_torn_ledger_tail(path)
    assert changed is False
    assert path.read_text(encoding="utf-8") == '{"a": 1}\n{"a": 2}\n'


def test_repair_is_noop_on_a_missing_file(tmp_path):
    assert daemon.repair_torn_ledger_tail(tmp_path / "missing.jsonl") is False


def test_repair_truncates_a_torn_trailing_line(tmp_path):
    path = tmp_path / "ledger.jsonl"
    complete = '{"a": 1}\n{"a": 2}\n'
    torn = '{"a": 3, "big_field": "some data that got cut of'  # no trailing newline, invalid JSON
    path.write_text(complete + torn, encoding="utf-8")
    changed = daemon.repair_torn_ledger_tail(path)
    assert changed is True
    remaining = path.read_text(encoding="utf-8")
    assert remaining == complete
    for line in remaining.splitlines():
        json.loads(line)  # every remaining line is valid JSON


def test_repair_never_touches_a_valid_line_missing_only_the_final_newline(tmp_path):
    """A file ending in valid JSON with no trailing newline is left alone --
    it is not proof of a torn write, just an edge case the check should not
    misfire on."""
    path = tmp_path / "ledger.jsonl"
    path.write_text('{"a": 1}\n{"a": 2}', encoding="utf-8")
    changed = daemon.repair_torn_ledger_tail(path)
    assert changed is False
