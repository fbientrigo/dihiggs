import pytest
from pathlib import Path

from dihiggs.app.orchestrator.aligned_lambda7 import build_aligned_lambda7_family
from dihiggs.app.orchestrator.analytic_bridge import evaluate_t4
from dihiggs.app.orchestrator.aligned_constraints import unitarity

POINT = build_aligned_lambda7_family(tan_beta="1000000", lambda1_ref="1.2",
    mh_GeV="125.13", mphi_GeV="150", mA_GeV="200", mHp_GeV="200",
    v_GeV="246", X_lambda6_tanbeta="1")


def test_bridge_rejects_unknown_mode_before_reading_point():
    with pytest.raises(ValueError, match="mode must be one of"):
        evaluate_t4({}, mode="unknown")

def test_bridge_validates_point_shape_before_dispatch():
    with pytest.raises(ValueError, match="input and value mappings"):
        evaluate_t4({}, mode="analytic")

def test_analytic_reordered_cp_even_state_is_inconclusive_without_downstream_calls():
    point = build_aligned_lambda7_family(tan_beta="100",lambda1_ref="1.2",
        mh_GeV="125.13",mphi_GeV="100",mA_GeV="200",mHp_GeV="200",
        v_GeV="246",X_lambda6_tanbeta="1")
    bridge=Path(__file__).resolve().parents[1]/"2hdmc/analytic-bridge"
    result = evaluate_t4(point, mode="analytic", executable=bridge)
    assert result["constructor_status"] == "INCONCLUSIVE"
    assert result["backend_constructor_status"] == "PASS"
    assert result["family_constructor_agreement_status"] == "NOT_TESTED"
    assert result["requested_h_index_backend"] == 2
    assert result["requested_phi_index_backend"] == 1
    degenerate=build_aligned_lambda7_family(tan_beta="1000000",lambda1_ref="1.2",
        mh_GeV="125.13",mphi_GeV="125.13",mA_GeV="200",mHp_GeV="200",
        v_GeV="246",X_lambda6_tanbeta="1")
    assert evaluate_t4(degenerate, mode="analytic")["constructor_status"] == "INCONCLUSIVE"
    reordered = build_aligned_lambda7_family(tan_beta="100",lambda1_ref="1.2",
        mh_GeV="125.13",mphi_GeV="100",mA_GeV="200",mHp_GeV="200",
        v_GeV="246",X_lambda6_tanbeta="1")
    with pytest.raises(ValueError, match="mphi_GeV >= mh_GeV"):
        evaluate_t4(reordered, mode="explicit", executable="/missing/bridge")


def test_bridge_reports_missing_executable_with_build_command(tmp_path):
    with pytest.raises(FileNotFoundError, match="make -C .* analytic-bridge"):
        evaluate_t4(POINT, mode="analytic", executable=tmp_path/"missing")

def test_backend_spectrum_is_exported_from_calc_unitarity():
    campaign=Path(__file__).resolve().parents[1]/"2hdmc/analytic-bridge"
    result=evaluate_t4(POINT,mode="analytic",executable=campaign)
    lambdas={f"lambda{i}":str(result[f"lambda{i}"]) for i in range(1,8)}
    independent=unitarity(lambdas,precision=80)
    for block in ("Y2_sigma1","Y2_sigma0","Y0_sigma1","Y0_sigma0"):
        backend=result[f"backend_eigenvalues_{block}"]
        reference=[float(x) for x in independent["eigenvalues"][block]]
        assert len(backend)==len(reference)
        assert max(abs(a-b) for a,b in zip(backend,reference))<5e-11


def test_analytic_reordered_state_mapping_preserves_t1_signed_conventions():
    point=build_aligned_lambda7_family(tan_beta="100",lambda1_ref="1.2",
        mh_GeV="125.13",mphi_GeV="100",mA_GeV="200",mHp_GeV="200",
        v_GeV="246",X_lambda6_tanbeta="1")
    bridge=Path(__file__).resolve().parents[1]/"2hdmc/analytic-bridge"
    result=evaluate_t4(point,mode="analytic",executable=bridge)
    assert result["physical_state_ordering"]=="MAPPED"
    assert result["requested_h_index_backend"]==2
    assert result["requested_phi_index_backend"]==1
    assert result["requested_phi_phase_factor"]==-result["requested_h_phase_factor"]
    assert abs(result["kappa_phi_top_backend_signed"]+0.01)<1e-12
    assert abs(result["g_hphi_phi_backend_GeV"]-2*100**2/246)<1e-8
    assert abs(result["g_phi_HpHm_backend_GeV"]-float(point.values["CphiHpHm_potential_GeV"]))<1e-8


def test_exact_cp_even_degeneracy_is_inconclusive_with_common_mode_metadata():
    point=build_aligned_lambda7_family(tan_beta="100",lambda1_ref="1.2",
        mh_GeV="125.13",mphi_GeV="125.13",mA_GeV="200",mHp_GeV="200",
        v_GeV="246",X_lambda6_tanbeta="1")
    result=evaluate_t4(point,mode="analytic")
    assert result["mode"]=="analytic"
    assert result["requested_mode"]=="analytic"
    assert result["constructor_status"]=="INCONCLUSIVE"
    assert result["backend_constructor_status"]=="NOT_TESTED"
    assert result["physical_state_ordering"]=="INCONCLUSIVE"
    assert result["widths"]=="NOT_TESTED"


def test_serialized_family_point_rejects_derived_field_mutation():
    bad=POINT.as_dict()
    bad["values"]["Z7"]="0"
    with pytest.raises(ValueError,match="derived field Z7"):
        evaluate_t4(bad,mode="analytic")
