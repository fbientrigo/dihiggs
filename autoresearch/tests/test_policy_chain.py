"""#17: selector -> abstain -> agy -> rules; with agy down or answering garbage the daemon goes on with rules."""
import sys

import pytest

from autoresearch.harness import dihiggs_runner
from autoresearch.harness.dihiggs_runner import DiHiggsRunner
from autoresearch.tests.test_journal_budget import launches, world

A, B = "survivor_v2|tb=100.0|l7_0p01_mh125p20", "survivor_v2|tb=10000000.0|l7_0p01_mh125p20"
FAKE_AGY = {  # fake local agy CLI (no network): the prompt arrives as the last argument
    "garbage": ["-c", "print('I would pick the tan beta 100 one')"],
    "down": ["-c", "import sys; sys.exit(3)"],
    "slow": ["-c", "import time; time.sleep(5)"],
    "valid": ["-c", f"print('`{B}`')"],
}


def _daemon(tmp_path, monkeypatch, agy_cmd, max_calls=10):
    monkeypatch.setattr(dihiggs_runner, "run_all_preflight_checks", lambda cfg: {"overall": "pass"})
    monkeypatch.setenv("FAKE_FRESH", "0")  # the lake already covers everything: no arm has support
    cfg = world(tmp_path)
    cfg["arms"] = [{"id": "a", "action": A}, {"id": "b", "action": B}, {"id": "mesh", "action": "adaptive_m2_mesh|tb=all|l7_0p01_mh125p20"}]
    cfg["policy"] = {"agy_cmd": agy_cmd, "agy_timeout_sec": 1, "agy_max_calls": max_calls}
    return DiHiggsRunner(cfg, str(tmp_path / "out"))


def _hops(runner, attempt):
    return [h["hop"] for h in runner.journal.get(f"decision:{attempt}")["hops"]]


@pytest.mark.parametrize("mode, hop", [("garbage", "agy:invalid_output"), ("down", "agy:error"), ("slow", "agy:timeout"),
                                       ("missing", "agy:unavailable")])
def test_agy_failing_falls_back_to_rules_and_every_hop_is_journaled(tmp_path, monkeypatch, mode, hop):
    runner = _daemon(tmp_path, monkeypatch, ["/nonexistent/agy"] if mode == "missing" else [sys.executable, *FAKE_AGY[mode]])
    for _ in range(4):
        runner.run_adaptation_round()
    # first two: UCB1 plays the unplayed arms; the mesh arm is masked as unavailable, never launched
    assert _hops(runner, "a0000") == _hops(runner, "a0001") == ["masked", "ucb1:selected"]
    assert runner.journal.get("decision:a0000")["law"] == "EXACT"
    assert _hops(runner, "a0002") == ["masked", "abstain:no_support", hop, "rules:selected"]
    assert runner.journal.get("decision:a0002")["law"] == "EXACT" and launches(tmp_path) == 4  # the daemon kept going
    assert all(c["action"] != "adaptive_m2_mesh|tb=all|l7_0p01_mh125p20" for _, _, c in runner.journal.events("closed"))


def test_valid_agy_answer_is_used_with_unknown_law_and_its_budget_is_kept(tmp_path, monkeypatch):
    runner = _daemon(tmp_path, monkeypatch, [sys.executable, *FAKE_AGY["valid"]], max_calls=1)
    for _ in range(4):
        runner.run_adaptation_round()
    third = runner.journal.get("decision:a0002")
    assert _hops(runner, "a0002")[-2:] == ["abstain:no_support", "agy:selected"] and third["law"] == "UNKNOWN" and third["action"] == B
    assert _hops(runner, "a0003")[-2:] == ["agy:budget_exhausted", "rules:selected"]


def test_every_arm_saturated_stops_without_agy_or_launch(tmp_path, monkeypatch):
    runner = _daemon(tmp_path, monkeypatch, [sys.executable, *FAKE_AGY["valid"]], max_calls=0)
    for _ in range(4):  # a, b, then one more 0-row attempt each: both arms' last 2 attempts are empty
        runner.run_adaptation_round()
    with pytest.raises(dihiggs_runner.policy.Saturated):
        runner.run_adaptation_round()
    hops = _hops(runner, "a0004")
    assert hops[-1] == "stop:saturated" and not any(h.startswith("agy:") for h in hops) and launches(tmp_path) == 4
    assert runner.journal.get("decision:a0004")["arm_id"] is None and not runner.journal.get("intent:a0004")
