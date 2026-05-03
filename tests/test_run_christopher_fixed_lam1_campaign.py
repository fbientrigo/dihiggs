from pathlib import Path

from scripts.run_christopher_fixed_lam1_campaign import (
    WrapperConfig,
    build_orchestrator_command,
    normalize_campaign_name,
    parse_float_csv,
)


def test_parse_float_csv():
    assert parse_float_csv("0.4,1.0,12.566370614359172") == [0.4, 1.0, 12.566370614359172]


def test_campaign_name_not_double_prefixed():
    assert normalize_campaign_name("campaign=abc") == "abc"
    assert normalize_campaign_name("abc") == "abc"


def test_lam1_fixed_mode_command_has_equal_bounds_and_one_point(tmp_path):
    cfg = WrapperConfig(
        campaign="christopher_fixed_lam1_2026apr",
        outdir="/tmp",
        lake_name="dihiggs_lake",
        exec_path="/tmp/PhysScanWithFixings",
        threads=8,
        dry_run=False,
        force=True,
        mphi_min=125.0,
        mphi_max=300.0,
        n_mphi=400,
        tanbeta=10000.0,
        mA=300.0,
        sin_ba=1.0,
        lambda7=0.0,
        lambda1_values=[0.4],
        lambda6_values=[0.01],
    )
    cmd = build_orchestrator_command(
        cfg=cfg,
        lam1=0.4,
        lam6=0.01,
        run_name="run_name",
        orchestrator_path=Path("/tmp/orchestrate_scans.py"),
    )
    joined = " ".join(cmd)
    assert " --lam1-min 0.4 " in f" {joined} "
    assert " --lam1-max 0.4 " in f" {joined} "
    assert " --n-lam1 1 " in f" {joined} "
    assert " --campaign christopher_fixed_lam1_2026apr " in f" {joined} "


def test_command_contains_requested_flags():
    cfg = WrapperConfig(
        campaign="camp",
        outdir="/tmp",
        lake_name="dihiggs_lake",
        exec_path="/tmp/PhysScanWithFixings",
        threads=4,
        dry_run=True,
        force=True,
        mphi_min=125.0,
        mphi_max=300.0,
        n_mphi=400,
        tanbeta=10000.0,
        mA=300.0,
        sin_ba=1.0,
        lambda7=0.0,
        lambda1_values=[1.0],
        lambda6_values=[0.1],
    )
    cmd = build_orchestrator_command(cfg, 1.0, 0.1, "rn", Path("/tmp/orchestrate_scans.py"))
    assert "--dry-run" in cmd
    assert "--force" in cmd
    assert "--threads" in cmd
