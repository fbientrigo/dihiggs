from __future__ import annotations

from pathlib import Path

from dihiggs.app.orchestrate_scans import FixedParams, build_run_dir


def test_fixed_dir_lambda6_precision_avoids_alias() -> None:
    base = Path("/tmp/lake_test")
    f1 = FixedParams(mA=500.0, sin_ba=1.0, lambda6=0.0019683, lambda7=0.0)
    f2 = FixedParams(mA=500.0, sin_ba=1.0, lambda6=0.0020273, lambda7=0.0)
    r1 = build_run_dir(outdir=base, lake_name="lake", campaign="c", fixed=f1, run_name="r")
    r2 = build_run_dir(outdir=base, lake_name="lake", campaign="c", fixed=f2, run_name="r")
    assert r1.parent != r2.parent
    assert "_l6=" in r1.parent.name and "_l6=" in r2.parent.name


def test_fixed_dir_lambda7_precision_avoids_alias() -> None:
    base = Path("/tmp/lake_test")
    f1 = FixedParams(mA=500.0, sin_ba=1.0, lambda6=0.001, lambda7=0.0019683)
    f2 = FixedParams(mA=500.0, sin_ba=1.0, lambda6=0.001, lambda7=0.0020273)
    r1 = build_run_dir(outdir=base, lake_name="lake", campaign="c", fixed=f1, run_name="r")
    r2 = build_run_dir(outdir=base, lake_name="lake", campaign="c", fixed=f2, run_name="r")
    assert r1.parent != r2.parent
    assert "_l7=" in r1.parent.name and "_l7=" in r2.parent.name


def test_build_run_dir_is_deterministic() -> None:
    base = Path("/tmp/lake_test")
    fixed = FixedParams(mA=500.0, sin_ba=1.0, lambda6=0.0019683, lambda7=0.0)
    a = build_run_dir(outdir=base, lake_name="lake", campaign="c", fixed=fixed, run_name="r")
    b = build_run_dir(outdir=base, lake_name="lake", campaign="c", fixed=fixed, run_name="r")
    assert a == b
