from pathlib import Path

from dihiggs.app.orchestrator.engines.m2 import M2Engine
from tests.test_orchestrator.conftest import make_fixed_m2, make_grid


def pairs(command: list[str]) -> dict[str, str]:
    return dict(zip(command[1::2], command[2::2]))


def test_m2_command_matches_real_named_cli():
    engine = M2Engine("campaign-17", "run-4")
    command = engine.build_command(
        Path("/repo/dihiggs/app/DihiggsPointV2Evaluator"),
        make_grid(mphi_min=125.13, mphi_max=500.0, n_mphi=7,
                  axis_min=0.0, axis_max=500_000.0, n_axis=11),
        make_fixed_m2(mA=510.0, sin_ba=0.999, tan_beta=50.0,
                      lambda6=0.001, lambda7=0.0),
        Path("/tmp/points.csv"),
    )
    assert command[0].endswith("DihiggsPointV2Evaluator")
    assert pairs(command) == {
        "--campaign-id": "campaign-17", "--run-id": "run-4", "--mh": "125.13",
        "--mH-min": "125.13", "--mH-max": "500", "--n-mH": "7",
        "--mA": "510", "--mHp": "510", "--yukawa-type": "1",
        "--sin-ba": "0.999", "--tan-beta": "50", "--M2-min": "0",
        "--M2-max": "500000", "--n-M2": "11", "--lambda6": "0.001",
        "--lambda7": "0", "--output": "/tmp/points.csv",
    }


def test_m2_engine_metadata_and_columns_are_canonical():
    engine = M2Engine()
    metadata = engine.axis_metadata()
    assert engine.executable_basename == "DihiggsPointV2Evaluator"
    assert metadata["schema_version"] == "dihiggs.point.v2"
    assert metadata["mass_convention"]["mh_GeV"] == 125.13
    assert "m12_sq /" in metadata["axis_description"]
    assert metadata["acceptance"]["theory_ok_v1"] == "triple_ok_legacy"
    assert {"point_id", "M2_input_GeV2", "experimental_ok"} <= set(engine.expected_csv_columns())


def test_m2_target_is_buildable():
    root = Path(__file__).resolve().parents[2]
    assert (root / "dihiggs/src/DihiggsPointV2Evaluator.cpp").is_file()
    assert "DihiggsPointV2Evaluator.cpp" in (root / "dihiggs/Makefile").read_text()
