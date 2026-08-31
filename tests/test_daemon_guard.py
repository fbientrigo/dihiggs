"""T12 (mission Sec. 26): unit-level checks that the daemon's path-allowlist
guard rejects writes outside the configured run-dir and refuses the live,
review-blocked campaign directories / deliverables/ by name. No physics is
evaluated in this file.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from search_discovery import daemon


def test_guard_run_dir_rejects_discovery_v1_by_name(tmp_path):
    protected = tmp_path / "runs" / "discovery_v1"
    with pytest.raises(daemon.DaemonConfigError, match="protected_run_dir"):
        daemon.guard_run_dir(protected)


def test_guard_run_dir_rejects_discovery_photonic_v1_by_name(tmp_path):
    protected = tmp_path / "runs" / "discovery_photonic_v1"
    with pytest.raises(daemon.DaemonConfigError, match="protected_run_dir"):
        daemon.guard_run_dir(protected)


def test_guard_run_dir_rejects_anything_under_deliverables(tmp_path):
    protected = tmp_path / "deliverables" / "some_new_subdir"
    with pytest.raises(daemon.DaemonConfigError, match="protected_path"):
        daemon.guard_run_dir(protected)


def test_guard_run_dir_allows_a_fresh_harvest_dir(tmp_path):
    ok = tmp_path / "runs" / "harvest_v1"
    daemon.guard_run_dir(ok)  # must not raise


# ---------------------------------------------------------------- bounded fix
def test_guard_run_dir_rejects_a_subdirectory_of_discovery_v1(tmp_path):
    """MAJOR fix: the guard used to match only run_dir's own basename, so
    `run_dir: runs/discovery_v1/anything` sailed through path containment."""
    protected_subdir = tmp_path / "runs" / "discovery_v1" / "some_subdir"
    with pytest.raises(daemon.DaemonConfigError, match="protected_run_dir"):
        daemon.guard_run_dir(protected_subdir)


def test_guard_run_dir_rejects_a_deeply_nested_subdirectory_of_discovery_photonic_v1(tmp_path):
    protected_subdir = tmp_path / "runs" / "discovery_photonic_v1" / "a" / "b" / "c"
    with pytest.raises(daemon.DaemonConfigError, match="protected_run_dir"):
        daemon.guard_run_dir(protected_subdir)


def test_guard_write_path_rejects_escaping_the_run_dir(tmp_path):
    run_dir = tmp_path / "runs" / "harvest_v1"
    run_dir.mkdir(parents=True)
    escaping = tmp_path / "runs" / "discovery_v1" / "ledger.jsonl"
    with pytest.raises(daemon.DaemonConfigError, match="daemon_write_outside_run_dir"):
        daemon.guard_write_path(escaping, run_dir)


def test_guard_write_path_allows_paths_under_the_run_dir(tmp_path):
    run_dir = tmp_path / "runs" / "harvest_v1"
    run_dir.mkdir(parents=True)
    daemon.guard_write_path(run_dir / "ledger.jsonl", run_dir)  # must not raise
    daemon.guard_write_path(run_dir / "attempts" / "x" / "evaluation.json", run_dir)  # must not raise


def test_run_refuses_to_start_against_protected_run_dir(tmp_path):
    config = daemon.DEFAULTS.copy()
    with pytest.raises(daemon.DaemonConfigError):
        daemon.run(config, tmp_path / "runs" / "discovery_v1", resume=False,
                  workers_override=1, dry_run=True)
