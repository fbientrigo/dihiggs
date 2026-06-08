"""
test_command_m2.py — Verify M2Engine command construction.

Tests that:
- The command has exactly 13 tokens.
- Positions 4–6 are M2_min, M2_max, N_M2 (GeV^2 values, NOT lambda1).
- fixed.lambda1 is required (ValueError if None).
- axis_metadata() unambiguously identifies the axis as M2 in GeV^2.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dihiggs.app.orchestrator.engines.m2 import M2Engine
from dihiggs.app.orchestrator.models import FixedParams
from tests.test_orchestrator.conftest import make_fixed_m2, make_grid


class TestM2CommandConstruction:
    """Unit tests for M2Engine.build_command()."""

    def setup_method(self):
        self.engine = M2Engine()
        # M2 axis: 0 to 500000 GeV^2
        self.grid = make_grid(
            mphi_min=125.0, mphi_max=500.0, n_mphi=50,
            axis_min=0.0, axis_max=500_000.0, n_axis=50,
        )
        self.fixed = make_fixed_m2(
            mA=500.0, sin_ba=1.0, tan_beta=50.0,
            lambda6=0.001, lambda7=0.0, lambda1=1.0,
        )
        self.output_csv = Path("/tmp/m2_scan_out.csv")
        self.exe = Path("/path/to/Phys_M2BoundaryScan")

    def _cmd(self) -> list[str]:
        return self.engine.build_command(
            executable=self.exe,
            grid=self.grid,
            fixed=self.fixed,
            output_csv=self.output_csv,
        )

    def test_token_count(self):
        """Command must have exactly 13 tokens."""
        cmd = self._cmd()
        assert len(cmd) == 13, f"Expected 13 tokens, got {len(cmd)}: {cmd}"

    def test_executable_is_first(self):
        cmd = self._cmd()
        assert cmd[0] == str(self.exe)

    def test_mphi_range_positions(self):
        """Tokens 1–3: mphi_min, mphi_max, n_mphi."""
        cmd = self._cmd()
        assert float(cmd[1]) == pytest.approx(125.0)
        assert float(cmd[2]) == pytest.approx(500.0)
        assert int(cmd[3]) == 50

    def test_m2_range_positions(self):
        """Tokens 4–6 must carry M^2 values (GeV^2), NOT lambda1 values."""
        cmd = self._cmd()
        # M2_min = 0.0 GeV^2
        assert float(cmd[4]) == pytest.approx(0.0), (
            "Position 4 should be M2_min (0.0 GeV^2)"
        )
        # M2_max = 500000.0 GeV^2 — this is clearly NOT a lambda1 value
        assert float(cmd[5]) == pytest.approx(500_000.0), (
            "Position 5 should be M2_max (500000 GeV^2), not a lambda1 value"
        )
        assert int(cmd[6]) == 50

    def test_fixed_params_positions(self):
        """Tokens 7–11: mA, sin_ba, tan_beta, lambda6, lambda7."""
        cmd = self._cmd()
        assert float(cmd[7]) == pytest.approx(500.0)   # mA
        assert float(cmd[8]) == pytest.approx(1.0)     # sin_ba
        assert float(cmd[9]) == pytest.approx(50.0)    # tan_beta
        assert float(cmd[10]) == pytest.approx(0.001)  # lambda6
        assert float(cmd[11]) == pytest.approx(0.0)    # lambda7

    def test_output_csv_is_last(self):
        """Token 12 must be the output CSV path."""
        cmd = self._cmd()
        assert cmd[12] == str(self.output_csv)

    def test_raises_if_lambda1_missing(self):
        """M2Engine must raise ValueError when fixed.lambda1 is None."""
        bad_fixed = FixedParams(
            mA=500.0, sin_ba=1.0, tan_beta=50.0,
            lambda6=0.001, lambda7=0.0,
            lambda1=None,  # Error: M2 scans require a fixed lambda1
        )
        with pytest.raises(ValueError, match="lambda1"):
            self.engine.build_command(
                executable=self.exe,
                grid=self.grid,
                fixed=bad_fixed,
                output_csv=self.output_csv,
            )

    def test_engine_name(self):
        assert self.engine.engine_name == "m2"

    def test_executable_basename(self):
        assert self.engine.executable_basename == "Phys_M2BoundaryScan"

    def test_axis_metadata_is_gev2(self):
        """axis_metadata must declare axis as M2 in GeV^2, not lambda1."""
        meta = self.engine.axis_metadata()
        assert meta["scan_axis"] == "M2"
        assert "GeV^2" in meta["axis_units"], (
            "M2 axis must be labelled as GeV^2"
        )
        assert "m12_sq" in meta["axis_label"] or "M2" in meta["axis_label"]

    def test_axis_metadata_documents_csv_column_ambiguity(self):
        """axis_metadata csv_column_note must warn about m12 vs m12_sq."""
        meta = self.engine.axis_metadata()
        note = meta["csv_column_note"].lower()
        # Must mention that 'm12' column stores GeV^2 (m12_sq)
        assert "m12" in note
        assert "gev" in note or "gev^2" in note or "m12_sq" in note, (
            f"csv_column_note should mention m12_sq or GeV^2; got: {meta['csv_column_note']}"
        )

    def test_m2_differs_from_lambda1_in_axis(self):
        """The M2 engine scan axis must differ from the lambda1 engine axis."""
        from dihiggs.app.orchestrator.engines.lambda1 import Lambda1Engine
        lam1_engine = Lambda1Engine()
        assert self.engine.scan_axis != lam1_engine.scan_axis
        assert self.engine.engine_name != lam1_engine.engine_name
