"""#16: write-ahead journal, budget never exceeded, kill -9 mid-action reconciles without duplicate or loss."""
import json
import os
import signal
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

import pytest
import yaml

from autoresearch.harness import boundary_actions as ba
from autoresearch.harness import dihiggs_runner
from autoresearch.harness.dihiggs_runner import DiHiggsRunner
from autoresearch.harness.journal import BudgetExhausted, Journal

KEY = "survivor_v2|tb=10000000.0|l7_0p01_mh125p20"
# fake boundary: records argv, sleeps $FAKE_SLEEP, writes points + manifest at the end (as the real one does)
FAKE_DHB = """import json, os, sys, time, yaml
from pathlib import Path
a = sys.argv; cfg = yaml.safe_load(open(a[a.index("--config") + 1]))
out = Path(a[a.index("--output-dir") + 1]) / cfg["campaign_id"]; out.mkdir(parents=True, exist_ok=True)
Path(os.environ["FAKE_LOG"]).open("a").write(json.dumps(a) + "\\n")
time.sleep(float(os.environ.get("FAKE_SLEEP", "0")))
(out / "raw_x_points.csv").write_text("theory_ok_v1\\n1\\n0\\n")
(out / "manifest.json").write_text(json.dumps({"budget": {"fresh_rows": int(os.environ.get("FAKE_FRESH", "2"))}, "explored_seed": {"seeded_attempts": 5}}))
"""
# fake lake: explored writes a csv, commit appends to commits.tsv, ls prints it (cid \\t ... \\t message)
FAKE_LAKE = """import sys
from pathlib import Path
log = Path(__file__).with_name("commits.tsv"); a = sys.argv[1:]
if a[0] == "explored":
    Path(a[a.index("--export") + 1]).write_text("x\\n")
elif a[0] == "commit":
    import os, time; time.sleep(float(os.environ.get("FAKE_COMMIT_SLEEP", "0")))
    n = len(log.read_text().splitlines()) if log.exists() else 0
    log.open("a").write(f"cid{n}\\tcandidate\\t{a[a.index('-m') + 1]}\\n"); print(f"committed cid{n}: 1 files")
elif a[0] == "ls":
    print(log.read_text() if log.exists() else "", end="")
"""


def world(tmp_path, **runtime):
    root = tmp_path / "dihiggs_boundary"
    for spec in ba.STRATEGIES.values():
        path = root / spec["config"]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump({"campaign_id": path.stem, "dihiggs_root": "../../dihiggs", "tan_beta_slices": [100.0, 1e7],
                                        "physics": {"lambda7": 0.0, "mh_GeV": 125.13}}))
    (root / "dhb").mkdir()
    (root / "dhb" / "__init__.py").write_text("")
    (root / "dhb" / "survivor_adaptive_campaign.py").write_text(FAKE_DHB)
    (tmp_path / "lake").mkdir()
    (tmp_path / "lake" / "lake.py").write_text(FAKE_LAKE)
    os.environ["FAKE_LOG"] = str(tmp_path / "launches.jsonl")
    return {"campaign_id": "t", "paths": {"boundary_root": str(root), "lake": str(tmp_path / "lake" / "lake.py")},
            "runtime": {"threads": 1, "timeout_sec": 30, **runtime}, "arms": [{"id": "a", "action": KEY}]}


def launches(tmp_path):
    log = tmp_path / "launches.jsonl"
    return len(log.read_text().splitlines()) if log.exists() else 0


def commits(tmp_path):
    log = tmp_path / "lake" / "commits.tsv"
    return log.read_text().splitlines() if log.exists() else []


def test_journal_is_append_only_and_idempotent(tmp_path):
    j = Journal(tmp_path / "j.db")
    assert j.append("x:1", "x", {"v": 1}) == {"v": 1} and j.append("x:1", "x", {"v": 2}) == {"v": 1}
    with pytest.raises(sqlite3.IntegrityError, match="append-only"):
        j.db.execute("UPDATE events SET kind='y'")
    with pytest.raises(sqlite3.IntegrityError, match="append-only"):
        j.db.execute("DELETE FROM events")


def test_budget_is_never_exceeded(tmp_path, monkeypatch):
    monkeypatch.setattr(dihiggs_runner, "run_all_preflight_checks", lambda cfg: {"overall": "pass"})
    cfg = {**world(tmp_path, disk_mb_per_action=1), "budget": {"disk_gb": 3 / 1024, "cpu_hours": 1, "wall_hours": 1}}
    runner = DiHiggsRunner(cfg, str(tmp_path / "out"))
    with pytest.raises(BudgetExhausted, match="disk_mb: used 3 \\+ reserve 1 > limit 3"):
        for _ in range(10):
            runner.run_adaptation_round()
            assert all(runner.journal.ledger()[r] <= runner.limits[r] for r in runner.limits)
    assert launches(tmp_path) == 3 and len(runner.journal.events("closed")) == 3 and not runner.journal.open_attempts()

    from autoresearch.harness.campaign_supervisor import CampaignSupervisor
    sup = CampaignSupervisor(cfg, outdir=tmp_path / "out", preflight_runner=lambda c: {"overall": "pass"})
    assert sup.run() == "COMPLETED" and "budget_exhausted" in json.loads(sup.state_path.read_text())["stop_reason"]
    assert launches(tmp_path) == 3


DAEMON = """import sys
sys.path.insert(0, {repo!r})
from autoresearch.harness import dihiggs_runner
dihiggs_runner.run_all_preflight_checks = lambda cfg: {{"overall": "pass"}}
dihiggs_runner.DiHiggsRunner({cfg!r}, {out!r}).run_adaptation_round()
"""


@pytest.mark.parametrize("kill_at, slow", [("started:a0000", "FAKE_SLEEP"), ("committing:a0000", "FAKE_COMMIT_SLEEP")])
def test_kill_9_mid_action_reconciles_without_duplicate_or_loss(tmp_path, monkeypatch, kill_at, slow):
    """A real daemon is SIGKILLed while boundary (or `lake commit`) runs; the restart waits out the orphan, commits once."""
    monkeypatch.setenv(slow, "2")
    monkeypatch.setattr(dihiggs_runner, "run_all_preflight_checks", lambda cfg: {"overall": "pass"})
    cfg, out = world(tmp_path), str(tmp_path / "out")
    repo = str(Path(__file__).resolve().parents[2])
    daemon = subprocess.Popen([sys.executable, "-c", DAEMON.format(repo=repo, cfg=cfg, out=out)])
    for _ in range(200):
        if Path(out, "journal.db").exists() and Journal(Path(out) / "journal.db").get(kill_at):
            break
        time.sleep(0.05)
    os.kill(daemon.pid, signal.SIGKILL)
    daemon.wait()
    j = Journal(Path(out) / "journal.db")
    assert j.open_attempts() == ["a0000"] and j.get(kill_at) and not j.get("committed:a0000")

    runner = DiHiggsRunner(cfg, out)
    result = runner.run_adaptation_round()
    assert result["reconciled"] and result["attempt"] == "a0000" and result["subprocess_status"] == "success"
    assert result["fresh_rows"] == 2 and result["lake_commit"] == "cid0" and launches(tmp_path) == 1 and len(commits(tmp_path)) == 1
    monkeypatch.setenv(slow, "0")
    assert runner.run_adaptation_round()["attempt"] == "a0001" and launches(tmp_path) == 2 and len(commits(tmp_path)) == 2
    # the round result of the reconciled attempt is in events.jsonl exactly once
    events = [json.loads(line) for line in Path(out, "events.jsonl").read_text().splitlines()]
    assert len(events) == 2


def test_crash_before_launch_or_after_commit_never_duplicates(tmp_path, monkeypatch):
    cfg, out = world(tmp_path), str(tmp_path / "out")
    runner = DiHiggsRunner(cfg, out)
    # crash after the intent, before launch: closed not_executed, reservation released, nothing ran
    runner.journal.reserve("a0000", {"cpu_core_s": 30, "wall_s": 30, "disk_mb": 1}, {}, {"arm_id": "a", "action": KEY,
                           "work": str(Path(out, "actions", "a0000")), "threads": 1, "timeout_sec": 30})
    assert runner._settle("a0000", False, reconcile=True)["subprocess_status"] == "not_executed"
    assert runner.journal.ledger() == {"cpu_core_s": 0, "wall_s": 0, "disk_mb": 0} and launches(tmp_path) == 0
    # crash after `lake commit`, before its journal row: the restart finds the commit by tag instead of committing again
    original = ba.lake_commit
    def commit_then_die(*args):
        original(*args)
        raise KeyboardInterrupt
    monkeypatch.setattr(ba, "lake_commit", commit_then_die)
    with pytest.raises(KeyboardInterrupt):
        runner._run_action("a", cfg["arms"][0], False)
    monkeypatch.setattr(ba, "lake_commit", original)
    result = DiHiggsRunner(cfg, out)._settle("a0001", False, reconcile=True)
    assert result["lake_commit"] == "cid0" and len(commits(tmp_path)) == 1 and launches(tmp_path) == 1
