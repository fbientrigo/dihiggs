import csv
import math
import os
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"
HEADER = (
    "schema_version,producer,producer_commit,producer_dirty,evaluator_api,campaign_id,run_id,point_id,"
    "yukawa_type,yukawa_type_installed,mh_input_GeV,mH_input_GeV,mA_input_GeV,mHp_input_GeV,sin_beta_minus_alpha_input,"
    "tan_beta_input,beta_input_rad,M2_input_GeV2,m12_sq_input_GeV2,lambda6_input,lambda7_input,"
    "lambda1_reconstructed,lambda2_reconstructed,lambda3_reconstructed,lambda4_reconstructed,"
    "lambda5_reconstructed,lambda6_reconstructed,lambda7_reconstructed,tan_beta_reconstructed,"
    "m12_sq_reconstructed_GeV2,M2_reconstructed_GeV2,construction_ok,numerical_ok,rejection_stage,rejection_reason,"
    "positivity_reported_ok,unitarity_ok,perturbativity_ok,stability_reported_ok,stability_dependency_alias,"
    "triple_ok_legacy,theory_ok_v1,experimental_evaluated,experimental_ok,g_hH2H2_GeV,width_bb_GeV,width_cc_GeV,"
    "width_tt_GeV,width_tautau_GeV,width_WW_GeV,width_ZZ_GeV,width_gammagamma_GeV,width_Zgamma_GeV,width_gg_GeV,"
    "width_hh_GeV,total_width_GeV,width_unaccounted_GeV,br_bb,br_cc,br_tt,br_tautau,br_WW,br_ZZ,br_gammagamma,br_Zgamma,br_gg,"
    "br_hh,width_ok,ctau_mm"
)


@pytest.fixture(scope="module", autouse=True)
def require_binary():
    if not BINARY.is_file():
        pytest.skip("build dihiggs/app/DihiggsPointV2Evaluator first")


def run_point(tmp_path: Path, name: str, *, mh=125.20, mH=130.0, mH_max=None, n_mH=1,
              mA=300.0, mHp=300.0, sba=0.999, tb=50.0, M2=16721.68154468371,
              M2_max=None, n_M2=1, l6=0.1, l7=0.0, campaign="pilot", run="run"):
    output = tmp_path / f"{name}.csv"
    command = [
        str(BINARY), "--campaign-id", campaign, "--run-id", run, "--mh", repr(mh),
        "--mH-min", repr(mH), "--mH-max", repr(mH if mH_max is None else mH_max),
        "--n-mH", str(n_mH), "--mA", repr(mA), "--mHp", repr(mHp),
        "--yukawa-type", "1", "--sin-ba", repr(sba), "--tan-beta", repr(tb),
        "--M2-min", repr(M2), "--M2-max", repr(M2 if M2_max is None else M2_max),
        "--n-M2", str(n_M2), "--lambda6", repr(l6), "--lambda7", repr(l7),
        "--output", str(output),
    ]
    env = {**os.environ, "DIHIGGS_GIT_COMMIT": "test-commit", "DIHIGGS_GIT_DIRTY": "no"}
    completed = subprocess.run(command, check=True, text=True, capture_output=True, env=env)
    with output.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return output, rows, completed


def test_exact_header_cardinality_and_provenance(tmp_path):
    output, rows, completed = run_point(
        tmp_path, "grid", mH=130.0, mH_max=131.0, n_mH=2,
        M2=16000.0, M2_max=17000.0, n_M2=3,
    )
    assert output.read_text().splitlines()[0] == HEADER
    assert len(rows) == 6
    assert completed.stdout.startswith("ROWS_WRITTEN 6\nTotal Attempts: 6\nTRIPLE_OK_POINTS ")
    assert {(row["schema_version"], row["campaign_id"], row["run_id"]) for row in rows} == {
        ("dihiggs.point.v2", "pilot", "run")
    }
    assert all("THDM::get_coupling_hhh" in row["evaluator_api"] for row in rows)


def test_adjacent_doubles_round_trip_and_get_distinct_ids(tmp_path):
    low = 16721.68154468371
    high = math.nextafter(low, math.inf)
    _, rows, _ = run_point(tmp_path, "adjacent", M2=low, M2_max=high, n_M2=2)
    assert [float(row["M2_input_GeV2"]) for row in rows] == [low, high]
    assert rows[0]["point_id"] != rows[1]["point_id"]
    assert all("e" in row["M2_input_GeV2"] for row in rows)


def test_ids_ignore_campaign_run_and_repeat_is_byte_deterministic(tmp_path):
    first, rows_a, _ = run_point(tmp_path, "first", campaign="c1", run="r1")
    second, rows_b, _ = run_point(tmp_path, "second", campaign="c2", run="r2")
    assert rows_a[0]["point_id"] == rows_b[0]["point_id"]
    third, _, _ = run_point(tmp_path, "third", campaign="c1", run="r1")
    assert first.read_bytes() == third.read_bytes()


def test_m2_reconstruction_width_br_lifetime_and_unevaluated_experiments(tmp_path):
    _, (row,), _ = run_point(tmp_path, "l01")
    beta = math.atan(float(row["tan_beta_input"]))
    expected_m12 = float(row["M2_input_GeV2"]) * math.sin(beta) * math.cos(beta)
    assert float(row["m12_sq_input_GeV2"]) == pytest.approx(expected_m12, rel=2e-16)
    assert float(row["M2_reconstructed_GeV2"]) == pytest.approx(float(row["M2_input_GeV2"]), rel=1e-12)
    assert math.isfinite(float(row["g_hH2H2_GeV"]))
    assert float(row["g_hH2H2_GeV"]) >= 0.0
    total = float(row["total_width_GeV"])
    assert row["width_ok"] == "1.00000000000000000e+00"
    assert float(row["br_gammagamma"]) == pytest.approx(float(row["width_gammagamma_GeV"]) / total)
    assert float(row["ctau_mm"]) == pytest.approx(1.973269804e-13 / total)
    assert row["yukawa_type_installed"] == "1.00000000000000000e+00"
    selected = sum(float(row[field]) for field in (
        "width_bb_GeV", "width_cc_GeV", "width_tt_GeV", "width_tautau_GeV", "width_WW_GeV",
        "width_ZZ_GeV", "width_gammagamma_GeV", "width_Zgamma_GeV", "width_gg_GeV", "width_hh_GeV",
    ))
    assert float(row["width_unaccounted_GeV"]) + selected == pytest.approx(total, rel=1e-14)
    assert row["experimental_evaluated"] == "0.00000000000000000e+00"
    assert row["experimental_ok"] == "nan"


def test_validated_h2_benchmark_exports_frozen_coupling(tmp_path):
    _, (row,), _ = run_point(
        tmp_path,
        "h2-benchmark",
        mh=125.20,
        mH=150.0,
        mA=450.0,
        mHp=450.0,
        sba=1.0,
        tb=300000.0,
        M2=22499.999999500335,
        l6=1e-10,
        l7=0.0,
    )
    assert row["construction_ok"] == "1"
    assert float(row["g_hH2H2_GeV"]) == pytest.approx(63.6625935034957138, rel=0.0, abs=1e-10)
    assert float(row["total_width_GeV"]) == pytest.approx(4.56119462052655178e-14, rel=0.0, abs=1e-24)
    assert float(row["br_bb"]) == pytest.approx(0.756737070372915044, rel=0.0, abs=1e-12)
    assert float(row["ctau_mm"]) == pytest.approx(4.32621268805276848, rel=0.0, abs=1e-10)


def test_top_pair_width_is_explicit_and_zero_below_threshold(tmp_path):
    _, (row,), _ = run_point(
        tmp_path, "below-tt", mH=250.0, mA=800.0, mHp=800.0, sba=1.0, tb=50.0,
        M2=16721.68154468371, l6=0.1, l7=0.0,
    )
    assert row["construction_ok"] == "1"
    assert float(row["width_tt_GeV"]) == 0.0
    assert float(row["br_tt"]) == 0.0
    assert float(row["total_width_GeV"]) == pytest.approx(4.15329937678322846e-06, rel=0.0, abs=1e-18)
    assert float(row["width_unaccounted_GeV"]) == pytest.approx(1.41522478536371598e-09, rel=0.0, abs=1e-20)


def test_top_pair_width_dominates_above_threshold(tmp_path):
    _, (row,), _ = run_point(
        tmp_path, "above-tt", mH=400.0, mA=800.0, mHp=800.0, sba=1.0, tb=50.0,
        M2=16721.68154468371, l6=0.1, l7=0.0,
    )
    assert row["construction_ok"] == "1"
    assert float(row["width_tt_GeV"]) == pytest.approx(1.67384971497812948e-03, rel=0.0, abs=1e-15)
    assert float(row["br_tt"]) == pytest.approx(0.966756638952538605, rel=0.0, abs=1e-12)
    assert float(row["width_tt_GeV"]) > 0.0
    # width_tt_GeV must never be silently folded into width_unaccounted_GeV.
    assert float(row["width_unaccounted_GeV"]) < 0.01 * float(row["width_tt_GeV"])


def test_ordering_boundary_emits_success_and_failure_with_nan_masks(tmp_path):
    below = math.nextafter(125.20, -math.inf)
    _, rows, _ = run_point(tmp_path, "ordering", mH=below, mH_max=125.20, n_mH=2)
    assert [row["construction_ok"] for row in rows] == ["0", "1"]
    assert rows[0]["rejection_reason"] == "mh_gt_mH"
    for field in ("numerical_ok", "lambda1_reconstructed", "theory_ok_v1", "g_hH2H2_GeV", "total_width_GeV", "width_ok"):
        assert rows[0][field] == "nan"
    assert rows[0]["yukawa_type_installed"] == "nan"


def test_construction_failure_preserves_row_without_decay_evaluation(tmp_path):
    _, (row,), _ = run_point(tmp_path, "construction-failure", mA=-1.0)
    assert row["construction_ok"] == "0"
    assert row["yukawa_type_installed"] == "nan"
    for field in ("g_hH2H2_GeV", "width_bb_GeV", "width_tautau_GeV", "total_width_GeV", "width_unaccounted_GeV", "ctau_mm"):
        assert row[field] == "nan"


@pytest.mark.parametrize(
    "name,kwargs,theory_ok",
    [
        ("L01", {}, 1.0),
        ("L05", {"M2": 16239.109978356435}, 0.0),
        ("L06", {"mH": 200.0, "mA": 500.0, "mHp": 500.0, "sba": 1.0,
                 "tb": 10000.0, "M2": 39999.9995713761, "l6": 1e-10}, 1.0),
    ],
)
def test_pilot_anchors_at_125_20(tmp_path, name, kwargs, theory_ok):
    _, (row,), _ = run_point(tmp_path, name, **kwargs)
    assert float(row["mh_input_GeV"]) == 125.20
    assert float(row["theory_ok_v1"]) == theory_ok
    if name == "L06":
        assert 5e-11 < float(row["total_width_GeV"]) < 7e-11
        assert 0.0004 < float(row["br_gammagamma"]) < 0.0006
        assert float(row["ctau_mm"]) < 0.01


def test_yukawa_type_i_is_active_and_type_ii_changes_fermionic_widths(tmp_path):
    _, (type_i,), _ = run_point(tmp_path, "type-i")
    output = tmp_path / "type-ii.csv"
    command = [
        str(BINARY), "--campaign-id", "pilot", "--run-id", "run", "--mh", "125.20",
        "--mH-min", "130", "--mH-max", "130", "--n-mH", "1", "--mA", "300", "--mHp", "300",
        "--yukawa-type", "2", "--sin-ba", "0.999", "--tan-beta", "50", "--M2-min", "16721.68154468371",
        "--M2-max", "16721.68154468371", "--n-M2", "1", "--lambda6", "0.1", "--lambda7", "0", "--output", str(output),
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    with output.open(newline="") as handle:
        type_ii = next(csv.DictReader(handle))
    assert type_i["yukawa_type_installed"] == "1.00000000000000000e+00"
    assert type_ii["yukawa_type_installed"] == "2.00000000000000000e+00"
    for field in ("construction_ok", "positivity_reported_ok", "unitarity_ok", "perturbativity_ok", "theory_ok_v1"):
        assert type_i[field] == type_ii[field]
    assert float(type_i["width_bb_GeV"]) != pytest.approx(float(type_ii["width_bb_GeV"]))


def test_unsupported_yukawa_type_is_rejected(tmp_path):
    output = tmp_path / "bad.csv"
    command = [str(BINARY), "--campaign-id", "pilot", "--run-id", "run", "--mh", "125.20",
               "--mH-min", "130", "--mH-max", "130", "--n-mH", "1", "--mA", "300", "--mHp", "300",
               "--yukawa-type", "0", "--sin-ba", "0.999", "--tan-beta", "50", "--M2-min", "1",
               "--M2-max", "1", "--n-M2", "1", "--lambda6", "0", "--lambda7", "0", "--output", str(output)]
    completed = subprocess.run(command, capture_output=True, text=True)
    assert completed.returncode == 2
    assert "unsupported_yukawa_type" in completed.stderr


def test_cli_and_output_failures_are_fatal(tmp_path):
    bad = subprocess.run([str(BINARY), "--campaign-id", "only"], text=True, capture_output=True)
    assert bad.returncode == 2
    output_directory = tmp_path / "is-a-directory.csv"
    output_directory.mkdir()
    with pytest.raises(subprocess.CalledProcessError):
        run_point(tmp_path, "is-a-directory")