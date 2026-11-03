import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "benchmarks/H2scan_mH150_tb300000_production_coupling.json"
COUPLINGS_F = ROOT / "2hdmc/MGME/2HDMC/couplings.f"
INTERACTIONS = ROOT / "2hdmc/MGME/2HDMC/interactions.dat"


def test_h2_production_coupling_replay_contract() -> None:
    data = json.loads(ARTIFACT.read_text())
    assert data["point_id"] == "H2scan_mH150_tb300000"
    assert data["construction"]["roundtrip_coordinate_rejected"] is True
    assert data["scalar_state_mapping"]["native_call"] == "THDM::get_coupling_hhh(1, 2, 2, c)"

    coupling = data["coupling"]
    assert coupling["two_hdmc_hphiphi_real_GeV"] == 0.0
    assert math.isclose(coupling["g_physical_hphiphi_GeV"], -coupling["two_hdmc_hphiphi_imag_GeV"], rel_tol=0, abs_tol=1e-12)
    assert coupling["g_physical_hphiphi_abs_GeV"] == abs(coupling["g_physical_hphiphi_GeV"])
    assert coupling["g_hH2H2_GeV"] == coupling["g_physical_hphiphi_abs_GeV"]
    opposite = data["opposite_sign_valid_point"]
    assert opposite["theory_ok_v1"] is True
    assert opposite["g_physical_hphiphi_GeV"] < 0.0 < coupling["g_physical_hphiphi_GeV"]
    assert opposite["g_physical_hphiphi_abs_GeV"] == abs(opposite["g_physical_hphiphi_GeV"])

    replay = data["replay"]
    assert replay["same_direct_set_param_phys_construction"] is True
    assert replay["width_and_branching_replay_status"] == "PASS"
    assert replay["total_width_GeV"] == 4.56118529862185007e-14
    assert replay["br_bb"] == 0.756737485808578692

    width = data["ufo_width_cross_check"]
    assert abs(width["reproduced_width_GeV"] - replay["total_width_GeV"]) <= width["absolute_tolerance_GeV"]

    assert data["evaluator_source_sha256"] == {
        "benchmarks/check_H2scan_mH150_tb300000.cpp":
        "3ce43838b366702e7ebeabd2a2a122ad1311df543100bfbad985c87faae8fae9"
    }


def test_vendored_ufo_mapping_is_verified_separately() -> None:
    data = json.loads(ARTIFACT.read_text())
    coupling = data["coupling"]
    assert "GH1H2H2=dcmplx(IMGH1H2H2,-REGH1H2H2)" in COUPLINGS_F.read_text()
    assert "h1  h2  h2 GH1H2H2 QED" in INTERACTIONS.read_text()

    returned = complex(
        coupling["two_hdmc_hphiphi_real_GeV"],
        coupling["two_hdmc_hphiphi_imag_GeV"],
    )
    ufo_parameter = complex(returned.imag, -returned.real)  # -i*c
    assert 1j * ufo_parameter == returned
