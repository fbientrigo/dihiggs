import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "benchmarks/H2scan_mH150_tb300000_production_coupling.json"


def test_h2_production_coupling_replay_contract() -> None:
    data = json.loads(ARTIFACT.read_text())
    assert data["point_id"] == "H2scan_mH150_tb300000"
    assert data["construction"]["roundtrip_coordinate_rejected"] is True
    assert data["scalar_state_mapping"]["native_call"] == "THDM::get_coupling_hhh(1, 2, 2, c)"

    coupling = data["coupling"]
    assert math.isclose(coupling["two_hdmc_returned_imag_GeV"], -coupling["g_hH2H2_GeV"], rel_tol=0, abs_tol=1e-12)
    assert math.isclose(coupling["converted_GHphiphi_GeV"], coupling["two_hdmc_returned_imag_GeV"], rel_tol=0, abs_tol=1e-12)

    replay = data["replay"]
    assert replay["same_direct_set_param_phys_construction"] is True
    assert replay["width_and_branching_replay_status"] == "PASS"
    assert replay["total_width_GeV"] == 4.56118529862185007e-14
    assert replay["br_bb"] == 0.756737485808578692

    width = data["ufo_width_cross_check"]
    assert abs(width["reproduced_width_GeV"] - replay["total_width_GeV"]) <= width["absolute_tolerance_GeV"]
