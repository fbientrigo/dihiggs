"""
test_thread_headroom.py — the orchestrator must leave CPU headroom by default.

Regression guard for the jun-2026 incident where a `--threads 12` run on a
12-core host pinned every core, starved sshd, and froze the box. The CLI now
defaults OMP_NUM_THREADS to (logical CPUs - 2) unless the user explicitly opts
out, so a bare invocation can never grab every core.
"""

from __future__ import annotations

import dihiggs.app.orchestrator.cli as cli
from dihiggs.app.orchestrator.cli import (
    _CPU_HEADROOM,
    _safe_default_threads,
    build_parser,
    resolve_omp_threads,
)


class TestSafeDefaultThreads:
    def test_reserves_headroom(self, monkeypatch):
        monkeypatch.setattr(cli.os, "cpu_count", lambda: 12)
        assert _safe_default_threads() == 12 - _CPU_HEADROOM

    def test_never_below_one(self, monkeypatch):
        """Tiny / unknown core counts still yield a runnable >=1 thread count."""
        monkeypatch.setattr(cli.os, "cpu_count", lambda: 1)
        assert _safe_default_threads() == 1
        monkeypatch.setattr(cli.os, "cpu_count", lambda: None)
        assert _safe_default_threads() == 1


class TestResolveOmpThreads:
    def test_default_is_capped(self, monkeypatch):
        monkeypatch.setattr(cli.os, "cpu_count", lambda: 12)
        # No flags -> auto headroom, NOT all cores.
        assert resolve_omp_threads(all_cores=False, threads=None) == 10

    def test_explicit_threads_passthrough(self):
        assert resolve_omp_threads(all_cores=False, threads=4) == 4

    def test_all_cores_opt_out(self):
        # Explicit opt-out -> None means "inherit env / use every core".
        assert resolve_omp_threads(all_cores=True, threads=None) is None

    def test_all_cores_beats_threads(self):
        assert resolve_omp_threads(all_cores=True, threads=8) is None


class TestParserFlags:
    def test_threads_defaults_none_at_parse_time(self):
        """Parser keeps None; the safe default is applied in resolution, so the
        'explicit value' vs 'unset' distinction is preserved."""
        args = build_parser().parse_args(["--engine", "lambda1"])
        assert args.threads is None
        assert args.all_cores is False

    def test_all_cores_flag_present(self):
        args = build_parser().parse_args(["--engine", "lambda1", "--all-cores"])
        assert args.all_cores is True
