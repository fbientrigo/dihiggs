import pytest
import sys
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
    bridge=Path(__file__).resolve().parents[1]/"campaigns/t2_t4_20261008T043858Z/scripts/analytic_bridge"
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
    campaign=Path(__file__).resolve().parents[1]/"campaigns/t2_t4_20261008T043858Z/scripts/analytic_bridge"
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
    bridge=Path(__file__).resolve().parents[1]/"campaigns/t2_t4_20261008T043858Z/scripts/analytic_bridge"
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


def test_native_diphoton_diagnostics_reconcile_qcd_components_and_width():
    campaign=Path(__file__).resolve().parents[1]/"campaigns/t2_t4_20261008T043858Z"
    sys.path.insert(0,str(campaign/"scripts"))
    from run_t4_bridge import _independent_diphoton, run_point
    rows=run_point(point_builder=lambda digits:POINT,precisions=(50,))
    analytic=next(row for row in rows if row["mode"]=="analytic")
    audit=analytic["independent_diphoton_components"]
    assert audit["status"]=="PASS"
    assert float(audit["absolute_width_residual_GeV"]) < 1e-25
    native=analytic["diphoton_native_diagnostics"]
    top=next(x for x in native["fermions"] if x["family"]=="up" and x["generation"]==3)
    before=complex(top["scalar_amplitude_before_qcd"]["real"],top["scalar_amplitude_before_qcd"]["imag"])
    after=complex(top["scalar_amplitude_after_qcd"]["real"],top["scalar_amplitude_after_qcd"]["imag"])
    assert abs(after-before)>0
    factor=complex(audit["top_qcd_component_factor"])
    assert 0 < factor.real < 1
    assert abs(after-before*factor)<1e-24
    qcd_off=evaluate_t4(POINT,mode="analytic",executable=campaign/"scripts/analytic_bridge",qcd_on=False)
    off=next(x for x in qcd_off["diphoton_native_diagnostics"]["fermions"] if x["family"]=="up" and x["generation"]==3)
    assert qcd_off["diphoton_native_diagnostics"]["qcd_enabled"] is False
    assert qcd_off["diphoton_native_diagnostics"]["loop_width_GeV"] > 0
    assert off["scalar_amplitude_after_qcd"]==off["scalar_amplitude_before_qcd"]
    assert qcd_off["width_diphoton_GeV"] == 0
    assert qcd_off["diphoton_width_channel_status"] == "NOT_TESTED"
    off_audit=_independent_diphoton(qcd_off)
    assert off_audit["status"]=="PASS"
    assert off_audit["production_width_channel_status"]=="NOT_TESTED"
    assert qcd_off["diphoton_native_diagnostics"]["scalar_amplitude_raw"] != {"real": 0, "imag": 0}


def test_independent_diphoton_zero_amplitude_uses_absolute_residuals():
    campaign=Path(__file__).resolve().parents[1]/"campaigns/t2_t4_20261008T043858Z"
    sys.path.insert(0,str(campaign/"scripts"))
    from run_t4_bridge import _independent_diphoton
    zero={"real":0,"imag":0}
    native={"state_phase":1,"mass_GeV":"125","vev_GeV":"246","alpha0":"0.0073",
        "mW_GeV":"80","mHp_GeV":"200","fermions":[],"heavy_damping_enabled":False,
        "qcd_enabled":True,"W_coupling_raw":zero,"charged_higgs_coupling_raw":zero,
        "scalar_amplitude_raw":zero,"pseudoscalar_amplitude_raw":zero,
        "alpha_s_at_mass":"0.1","W_amplitude_raw":zero,"charged_higgs_amplitude_raw":zero,
        "loop_width_GeV":"0"}
    result=_independent_diphoton({"diphoton_native_diagnostics":native,"width_diphoton_GeV":"0"})
    assert result["status"]=="PASS"
    assert result["scalar_amplitude_normalized_residual"] is None
    assert result["pseudoscalar_amplitude_normalized_residual"] is None
    assert result["absolute_width_residual_GeV"]=="0.0"
