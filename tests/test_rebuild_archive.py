"""T4 (mission Sec. 26): the discovery archive must be reconstructible from
the ledger alone. Builds a small, real (non-live) run-dir via a couple of
bounded daemon cycles, then rebuilds the archive into a scratch location and
compares it against the live one.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
DAEMON = ROOT / "benchmark_search_daemon.py"


def _run_bounded_campaign(run_dir: Path, cycles: int = 3) -> None:
    config = {
        "campaign": {"name": "rebuild_archive_test"},
        "run_dir": str(run_dir),
        "evaluator": str(ROOT / "dihiggs/app/DihiggsPointV2Evaluator"),
        "envelope_id": "E1_mixed_exploit",
        "runtime": {"workers": 2, "max_evaluations_per_cycle": 16,
                    "checkpoint_every": 1, "summary_every": 1, "shutdown_grace_s": 10},
        "budgets": {"max_total_evaluations": None, "max_cycle_walltime_s": 300, "max_total_walltime_s": None},
        "policy": {"seed": 20260825, "family_validate_per_cycle": 1,
                   "allocation": {"explore": 0.4, "refine": 0.3, "optimize_mixed": 0.3}},
        "stopping": {"patience_cycles": 100, "max_consecutive_evaluator_failures": 5},
        "safety": {"downstream_snapshots_read_only": True},
    }
    run_dir.parent.mkdir(parents=True, exist_ok=True)
    config_path = run_dir.parent / f"{run_dir.name}_config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(DAEMON), "--config", str(config_path), "--max-cycles", str(cycles)],
        cwd=ROOT, capture_output=True, text=True, timeout=90,
    )
    assert result.returncode == 0, result.stderr


def test_rebuilt_archive_matches_the_live_one(tmp_path):
    run_dir = tmp_path / "runs" / "rebuild_test"
    _run_bounded_campaign(run_dir, cycles=3)

    live_archive_path = run_dir / "discovery_archive.json"
    assert live_archive_path.exists(), "campaign should have produced at least one family validation"
    live_archive = json.loads(live_archive_path.read_text(encoding="utf-8"))
    assert sum(len(v) for v in live_archive.values()) > 0, "expected at least one archived family"

    from search_discovery.rebuild_archive import rebuild
    target_dir = tmp_path / "runs" / "rebuild_scratch"
    rebuilt = rebuild(run_dir, target_dir=target_dir)

    assert rebuilt == live_archive


def test_rebuild_does_not_modify_the_live_archive_file(tmp_path):
    run_dir = tmp_path / "runs" / "rebuild_no_mutate_test"
    _run_bounded_campaign(run_dir, cycles=2)
    live_archive_path = run_dir / "discovery_archive.json"
    before = live_archive_path.read_bytes() if live_archive_path.exists() else None

    from search_discovery.rebuild_archive import rebuild
    rebuild(run_dir, target_dir=tmp_path / "runs" / "scratch2")

    after = live_archive_path.read_bytes() if live_archive_path.exists() else None
    assert before == after


# ---------------------------------------------------------------- bounded fix
def test_rebuild_refuses_a_protected_run_dir_without_an_explicit_target(tmp_path):
    """BLOCKER fix: rebuild() with the default target_dir=None must not be
    usable to wipe a review-blocked live campaign's archive+history -- it
    must refuse via the same guard daemon.py uses, not reimplement it."""
    from search_discovery.daemon import DaemonConfigError
    from search_discovery.rebuild_archive import rebuild

    protected = tmp_path / "runs" / "discovery_v1"
    protected.mkdir(parents=True)
    (protected / "ledger.jsonl").write_text("", encoding="utf-8")
    with pytest.raises(DaemonConfigError, match="protected_run_dir"):
        rebuild(protected)


def test_rebuild_never_deletes_an_existing_history_file(tmp_path):
    """BLOCKER fix: discovery_archive_history.jsonl is append-only eviction
    evidence and must never be unlinked by rebuild, even for an in-place
    (non-protected) rebuild where target_dir == run_dir."""
    run_dir = tmp_path / "runs" / "rebuild_history_preserved_test"
    _run_bounded_campaign(run_dir, cycles=3)
    history_path = run_dir / "discovery_archive_history.jsonl"
    assert history_path.exists(), "campaign should have produced at least one archive decision"
    before_lines = history_path.read_text(encoding="utf-8").splitlines()
    assert before_lines, "expected at least one history line to test preservation against"

    from search_discovery.rebuild_archive import rebuild
    rebuild(run_dir)  # in-place: target_dir defaults to run_dir itself

    assert history_path.exists(), "rebuild must never unlink discovery_archive_history.jsonl"
    after_lines = history_path.read_text(encoding="utf-8").splitlines()
    # replay may append further ARCHIVE_CONSIDER rows on top, but every
    # pre-existing line must still be present, in order, byte-for-byte
    assert after_lines[: len(before_lines)] == before_lines
