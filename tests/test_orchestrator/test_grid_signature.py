"""
test_grid_signature.py — Verify that grid signatures differ across engines.

Key invariants:
- Same numeric ranges + same fixed params → DIFFERENT signatures for
  lambda1 vs M2 engines (engine name is included in the hash).
- Same engine + same ranges + same fixed → same signature (deterministic).
- Different numeric ranges → different signatures (even for same engine).
- fixed.lambda1 participates in the M2 signature (it's a fixed param).
"""

from __future__ import annotations

import pytest

from dihiggs.app.orchestrator.grid import ScanGrid, grid_signature
from dihiggs.app.orchestrator.models import FixedParams
from tests.test_orchestrator.conftest import make_fixed_lam1, make_fixed_m2, make_grid


class TestGridSignature:

    # ---- shared fixtures ----------------------------------------------------

    def _grid(self) -> ScanGrid:
        """A grid with the same numbers for both engine types."""
        return make_grid(
            mphi_min=125.0, mphi_max=500.0, n_mphi=10,
            axis_min=0.0, axis_max=12.0, n_axis=100,
        )

    def _fixed_lam1(self) -> FixedParams:
        return make_fixed_lam1(
            mA=500.0, sin_ba=1.0, tan_beta=50.0, lambda6=0.001, lambda7=0.0
        )

    def _fixed_m2(self) -> FixedParams:
        return make_fixed_m2(
            mA=500.0, sin_ba=1.0, tan_beta=50.0, lambda6=0.001, lambda7=0.0,
            lambda1=1.0,
        )

    # ---- tests --------------------------------------------------------------

    def test_different_engines_produce_different_signatures(self):
        """
        Critical: lambda1 and M2 signatures must differ even for identical
        numeric parameters.  This prevents resume collisions.
        """
        grid = self._grid()
        # Use the same numeric fixed params (excluding lambda1 field)
        fixed_lam1 = make_fixed_lam1(mA=500.0, sin_ba=1.0, tan_beta=50.0,
                                      lambda6=0.001, lambda7=0.0)
        fixed_m2 = make_fixed_m2(mA=500.0, sin_ba=1.0, tan_beta=50.0,
                                   lambda6=0.001, lambda7=0.0, lambda1=1.0)

        sig_lam1 = grid_signature("lambda1", grid, fixed_lam1)
        sig_m2 = grid_signature("m2", grid, fixed_m2)

        assert sig_lam1 != sig_m2, (
            "Lambda1 and M2 signatures must differ to prevent resume collisions. "
            f"Both produced: {sig_lam1}"
        )

    def test_same_engine_same_params_deterministic(self):
        """Same inputs must always produce the same signature."""
        grid = self._grid()
        fixed = self._fixed_lam1()

        sig1 = grid_signature("lambda1", grid, fixed)
        sig2 = grid_signature("lambda1", grid, fixed)

        assert sig1 == sig2

    def test_different_axis_ranges_produce_different_signatures(self):
        """Changing axis_max must change the signature."""
        fixed = self._fixed_lam1()
        grid_a = make_grid(axis_min=0.0, axis_max=12.0, n_axis=100)
        grid_b = make_grid(axis_min=0.0, axis_max=6.0, n_axis=100)

        sig_a = grid_signature("lambda1", grid_a, fixed)
        sig_b = grid_signature("lambda1", grid_b, fixed)

        assert sig_a != sig_b

    def test_different_n_axis_produces_different_signatures(self):
        """Changing n_axis must change the signature."""
        fixed = self._fixed_lam1()
        grid_a = make_grid(n_axis=100)
        grid_b = make_grid(n_axis=200)

        assert grid_signature("lambda1", grid_a, fixed) != \
               grid_signature("lambda1", grid_b, fixed)

    def test_different_mphi_range_produces_different_signatures(self):
        """Changing mphi_max must change the signature."""
        fixed = self._fixed_lam1()
        grid_a = make_grid(mphi_max=290.0)
        grid_b = make_grid(mphi_max=500.0)

        assert grid_signature("lambda1", grid_a, fixed) != \
               grid_signature("lambda1", grid_b, fixed)

    def test_different_tan_beta_produces_different_signatures(self):
        """Per-task signatures differ by tan_beta."""
        grid = self._grid()
        fixed_a = make_fixed_lam1(tan_beta=10000.0)
        fixed_b = make_fixed_lam1(tan_beta=20000.0)

        assert grid_signature("lambda1", grid, fixed_a) != \
               grid_signature("lambda1", grid, fixed_b)

    def test_lambda1_fixed_param_included_in_m2_signature(self):
        """Changing fixed.lambda1 must change the M2 signature."""
        grid = self._grid()
        fixed_a = make_fixed_m2(lambda1=1.0)
        fixed_b = make_fixed_m2(lambda1=2.0)

        sig_a = grid_signature("m2", grid, fixed_a)
        sig_b = grid_signature("m2", grid, fixed_b)

        assert sig_a != sig_b, (
            "fixed.lambda1 must participate in M2 grid signature"
        )

    def test_signature_length(self):
        """Signature is always exactly 16 hex characters."""
        sig = grid_signature("lambda1", self._grid(), self._fixed_lam1())
        assert len(sig) == 16
        assert all(c in "0123456789abcdef" for c in sig)

    def test_engine_name_in_signature_is_case_sensitive(self):
        """'lambda1' and 'Lambda1' produce different signatures."""
        grid = self._grid()
        fixed = self._fixed_lam1()
        assert grid_signature("lambda1", grid, fixed) != \
               grid_signature("Lambda1", grid, fixed)
