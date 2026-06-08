"""
conftest.py — Shared test fixtures for test_orchestrator.

Provides:
- tmp_outdir        : temporary output directory (Path)
- fake_success_exe  : path to a bash script that simulates a successful C++ scan
- fake_failure_exe  : path to a bash script that simulates a failing C++ scan
- make_grid         : factory for ScanGrid instances
- make_fixed_lam1   : factory for FixedParams (lambda1 engine)
- make_fixed_m2     : factory for FixedParams (m2 engine, includes lambda1)
"""

from __future__ import annotations

import os
import stat
import textwrap
from pathlib import Path

import pytest

from dihiggs.app.orchestrator.grid import ScanGrid
from dihiggs.app.orchestrator.models import FixedParams


# ---------------------------------------------------------------------------
# Fake executables
# ---------------------------------------------------------------------------

@pytest.fixture()
def fake_success_exe(tmp_path: Path) -> Path:
    """
    A bash script that mimics a successful C++ scan:
    - Writes a minimal CSV to the last positional argument (output path).
    - Prints ``Total Attempts: 100`` and ``TRIPLE_OK_POINTS 42`` to stdout.
    - Exits 0.
    """
    script = tmp_path / "fake_success.sh"
    script.write_text(
        textwrap.dedent("""\
            #!/usr/bin/env bash
            # Last argument is the output CSV path
            OUTPUT="${@: -1}"
            echo "mphi,lam1_or_m2,ctau_m" > "$OUTPUT"
            echo "200.0,1.0,5.39e-4"      >> "$OUTPUT"
            echo "Total Attempts: 100"
            echo "TRIPLE_OK_POINTS 42"
        """),
        encoding="utf-8",
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return script


@pytest.fixture()
def fake_failure_exe(tmp_path: Path) -> Path:
    """
    A bash script that mimics a failing C++ scan:
    - Writes nothing.
    - Prints an error to stderr.
    - Exits 1.
    """
    script = tmp_path / "fake_failure.sh"
    script.write_text(
        textwrap.dedent("""\
            #!/usr/bin/env bash
            echo "Fatal: parameter out of range" >&2
            echo "Attempting cleanup..." >&1
            exit 1
        """),
        encoding="utf-8",
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return script


# ---------------------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------------------

@pytest.fixture()
def tmp_outdir(tmp_path: Path) -> Path:
    """Temporary root output directory for scan results."""
    d = tmp_path / "outdir"
    d.mkdir()
    return d


# ---------------------------------------------------------------------------
# Grid factories
# ---------------------------------------------------------------------------

def make_grid(
    *,
    mphi_min: float = 125.0,
    mphi_max: float = 500.0,
    n_mphi: int = 5,
    axis_min: float = 0.0,
    axis_max: float = 12.0,
    n_axis: int = 10,
) -> ScanGrid:
    """Factory for ScanGrid with sensible defaults."""
    return ScanGrid(
        mphi_min=mphi_min,
        mphi_max=mphi_max,
        n_mphi=n_mphi,
        axis_min=axis_min,
        axis_max=axis_max,
        n_axis=n_axis,
    )


def make_fixed_lam1(
    *,
    mA: float = 500.0,
    sin_ba: float = 1.0,
    tan_beta: float = 50.0,
    lambda6: float = 0.001,
    lambda7: float = 0.0,
) -> FixedParams:
    """Factory for FixedParams for the lambda1 engine (no lambda1 in fixed)."""
    return FixedParams(
        mA=mA,
        sin_ba=sin_ba,
        tan_beta=tan_beta,
        lambda6=lambda6,
        lambda7=lambda7,
        lambda1=None,  # lambda1 is the scan axis
    )


def make_fixed_m2(
    *,
    mA: float = 500.0,
    sin_ba: float = 1.0,
    tan_beta: float = 50.0,
    lambda6: float = 0.001,
    lambda7: float = 0.0,
    lambda1: float = 1.0,
) -> FixedParams:
    """Factory for FixedParams for the M2 engine (lambda1 is fixed constant)."""
    return FixedParams(
        mA=mA,
        sin_ba=sin_ba,
        tan_beta=tan_beta,
        lambda6=lambda6,
        lambda7=lambda7,
        lambda1=lambda1,  # lambda1 is a fixed constant for M2 scans
    )
