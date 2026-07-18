from __future__ import annotations

from pathlib import Path

from dihiggs.app.orchestrator.cli import build_parser, main


def test_parser_defaults_to_canonical_lambda1_v2_and_exposes_tracker() -> None:
    parser = build_parser()
    args = parser.parse_args([])
    assert args.engine == "lambda1_v2"
    assert set((parser._option_string_actions["--engine"].choices)) >= {
        "lambda1_v2", "lambda1_legacy", "m2", "m2_tracker",
    }


def test_readme_uses_current_orchestrator_output_flag() -> None:
    readme = Path(__file__).resolve().parents[2] / "README.md"
    text = readme.read_text()
    assert "--outdir" in text
    assert "--output-dir" not in text


def test_m2_rejects_fixed_lambda1(capsys) -> None:
    assert main(["--engine", "m2", "--lambda1", "1", "--dry-run"]) == 2
    error = capsys.readouterr().err
    assert "rejected" in error
    assert "reconstructed output" in error


def test_lambda1_v2_dry_run_writes_manifest_and_exact_input(tmp_path: Path) -> None:
    assert main([
        "--engine", "lambda1_v2", "--dry-run", "--outdir", str(tmp_path),
        "--lake-name", "lake", "--campaign", "pilot", "--run-name", "run",
        "--mH-min", "130", "--mH-max", "130", "--n-mH", "1",
        "--axis-min", "1", "--axis-max", "2", "--n-axis", "2",
        "--tanbeta", "50,100",
    ]) == 0
    run = tmp_path / "lake/campaign=pilot/run"
    assert run.joinpath("input_v2.csv").read_text().splitlines()[0] == (
        "point_id,mh_gev,mH_gev,mA_gev,mHp_gev,sin_beta_minus_alpha,"
        "tan_beta,lambda1_target,lambda6_input,lambda7_input"
    )
