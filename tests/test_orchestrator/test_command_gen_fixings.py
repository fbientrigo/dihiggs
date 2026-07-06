"""
test_command_gen_fixings.py — Verify GenFixingsEngine command construction.

Tests that:
- The command uses named flags (--bronze-csv, --output-csv, --calibration-*).
- fixed.bronze_shard_csv is required (ValueError if missing).
- Optional calibration flags are emitted only when set.
- axis_metadata() declares the GEN_FIXINGS (no generated grid) axis.
- expected_csv_columns() advertises the variation cloud column (variation_idx)
  and the chris cross-check / ctau columns the comparison report relies on.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dihiggs.app.orchestrator.engines.gen_fixings import GenFixingsEngine
from dihiggs.app.orchestrator.engines.lambda1 import Lambda1Engine
from dihiggs.app.orchestrator.models import FixedParams
from tests.test_orchestrator.conftest import make_grid


def _fixed(**overrides) -> FixedParams:
    base = dict(
        mA=300.0, sin_ba=1.0, tan_beta=1000.0, lambda6=1e-4, lambda7=0.0,
        bronze_shard_csv="/tmp/bronze_all.csv",
        calibration_n=10, calibration_frac=0.10, rng_seed=0,
    )
    base.update(overrides)
    return FixedParams(**base)


class TestGenFixingsCommandConstruction:
    def setup_method(self):
        self.engine = GenFixingsEngine()
        self.grid = make_grid()  # ignored by this engine, but required by the protocol
        self.output_csv = Path("/tmp/silver_out.csv")
        self.exe = Path("/path/to/GenScanWithFixings")

    def _cmd(self, fixed: FixedParams) -> list[str]:
        return self.engine.build_command(
            executable=self.exe, grid=self.grid, fixed=fixed, output_csv=self.output_csv,
        )

    def test_executable_is_first(self):
        cmd = self._cmd(_fixed())
        assert cmd[0] == str(self.exe)

    def test_named_flags_present(self):
        cmd = self._cmd(_fixed())
        assert "--bronze-csv=/tmp/bronze_all.csv" in cmd
        assert f"--output-csv={self.output_csv}" in cmd
        assert "--calibration-n=10" in cmd
        assert "--calibration-frac=0.1" in cmd
        assert "--rng-seed=0" in cmd

    def test_bronze_csv_required(self):
        with pytest.raises(ValueError):
            self._cmd(_fixed(bronze_shard_csv=None))

    def test_optional_flags_omitted_when_unset(self):
        cmd = self._cmd(_fixed(calibration_n=None, calibration_frac=None, rng_seed=None))
        assert not any(t.startswith("--calibration-n=") for t in cmd)
        assert not any(t.startswith("--calibration-frac=") for t in cmd)
        assert not any(t.startswith("--rng-seed=") for t in cmd)
        # bronze/output are still present
        assert any(t.startswith("--bronze-csv=") for t in cmd)
        assert any(t.startswith("--output-csv=") for t in cmd)

    def test_engine_identity(self):
        assert self.engine.engine_name == "gen_fixings"
        assert self.engine.executable_basename == "GenScanWithFixings"

    def test_axis_is_gen_fixings_not_lambda1(self):
        lam1 = Lambda1Engine()
        assert self.engine.scan_axis != lam1.scan_axis
        meta = self.engine.axis_metadata()
        assert meta["scan_axis"] == "gen_fixings"

    def test_expected_columns_include_variation_and_crosscheck(self):
        cols = self.engine.expected_csv_columns()
        # the variation cloud index (one row per +/-10% candidate)
        assert "variation_idx" in cols
        # chris cross-check + precomputed comparison columns the report uses
        for c in ("chris_width_gaga", "chris_ctau_mm", "ratio_ctau_mm",
                  "ratio_width_gaga", "br_gaga", "calibration_score", "stability_ok"):
            assert c in cols, f"missing expected column: {c}"
