"""
test_command_lambda1.py — Verify Lambda1Engine command construction.

Tests that:
- The command has exactly 13 tokens.
- Positions 1–3 are mphi_min, mphi_max, n_mphi.
- Positions 4–6 are lam1_min, lam1_max, N_lam1 (NOT M2).
- Positions 7–11 are mA, sin_ba, tan_beta, lambda6, lambda7.
- Position 12 is the output CSV path.
- ValueError is raised if fixed.lambda1 is set (axis confusion guard).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dihiggs.app.orchestrator.engines.lambda1 import Lambda1Engine
from tests.test_orchestrator.conftest import make_fixed_lam1, make_grid


class TestLambda1CommandConstruction:
    """Unit tests for Lambda1Engine.build_command()."""

    def setup_method(self):
        self.engine = Lambda1Engine()
        self.grid = make_grid(
            mphi_min=130.0, mphi_max=290.0, n_mphi=15,
            axis_min=0.0, axis_max=12.0, n_axis=666,
        )
        self.fixed = make_fixed_lam1(
            mA=300.0, sin_ba=0.995, tan_beta=10000.0,
            lambda6=0.1, lambda7=0.0,
        )
        self.output_csv = Path("/tmp/scan_out.csv")
        self.exe = Path("/path/to/PhysScanWithFixings")

    def _cmd(self) -> list[str]:
        return self.engine.build_command(
            executable=self.exe,
            grid=self.grid,
            fixed=self.fixed,
            output_csv=self.output_csv,
        )

    def test_token_count(self):
        """Command must have exactly 13 tokens (exe + 12 positional args)."""
        cmd = self._cmd()
        assert len(cmd) == 13, f"Expected 13 tokens, got {len(cmd)}: {cmd}"

    def test_executable_is_first(self):
        cmd = self._cmd()
        assert cmd[0] == str(self.exe)

    def test_mphi_range_positions(self):
        """Tokens 1–3 must be mphi_min, mphi_max, n_mphi."""
        cmd = self._cmd()
        assert float(cmd[1]) == pytest.approx(130.0)
        assert float(cmd[2]) == pytest.approx(290.0)
        assert int(cmd[3]) == 15

    def test_lambda1_range_positions(self):
        """Tokens 4–6 must be lam1_min, lam1_max, N_lam1 (NOT M2 values)."""
        cmd = self._cmd()
        assert float(cmd[4]) == pytest.approx(0.0), (
            "Position 4 should be lam1_min; got something unexpected"
        )
        assert float(cmd[5]) == pytest.approx(12.0), (
            "Position 5 should be lam1_max; got something unexpected"
        )
        assert int(cmd[6]) == 666

    def test_fixed_params_positions(self):
        """Tokens 7–11 must be mA, sin_ba, tan_beta, lambda6, lambda7."""
        cmd = self._cmd()
        assert float(cmd[7]) == pytest.approx(300.0)   # mA
        assert float(cmd[8]) == pytest.approx(0.995)   # sin_ba
        assert float(cmd[9]) == pytest.approx(10000.0) # tan_beta
        assert float(cmd[10]) == pytest.approx(0.1)    # lambda6
        assert float(cmd[11]) == pytest.approx(0.0)    # lambda7

    def test_output_csv_is_last(self):
        """Token 12 (last) must be the output CSV path."""
        cmd = self._cmd()
        assert cmd[12] == str(self.output_csv)

    def test_raises_if_lambda1_in_fixed(self):
        """If fixed.lambda1 is set, Lambda1Engine must raise ValueError."""
        from dihiggs.app.orchestrator.models import FixedParams
        bad_fixed = FixedParams(
            mA=300.0, sin_ba=0.995, tan_beta=10000.0,
            lambda6=0.1, lambda7=0.0,
            lambda1=1.5,  # Error: lambda1 should not be in fixed for this engine
        )
        with pytest.raises(ValueError, match="lambda1"):
            self.engine.build_command(
                executable=self.exe,
                grid=self.grid,
                fixed=bad_fixed,
                output_csv=self.output_csv,
            )

    def test_engine_name(self):
        assert self.engine.engine_name == "lambda1"

    def test_executable_basename(self):
        assert self.engine.executable_basename == "PhysScanWithFixings"

    def test_axis_metadata_no_m2(self):
        """axis_metadata must declare the axis as lambda1, not M2."""
        meta = self.engine.axis_metadata()
        assert meta["scan_axis"] == "lambda1"
        assert "GeV^2" not in meta["axis_units"]
        assert "M2" not in meta["axis_label"]
        # Must warn explicitly about the distinction
        assert "lambda1" in meta["csv_column_note"].lower() or \
               "conflate" in meta["csv_column_note"].lower()
