#!/usr/bin/env python3
"""Reproducibly compare the two physical and one analytic 2HDMC routes."""
from __future__ import annotations

import json
import sys
from pathlib import Path
import time
import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from dihiggs.app.orchestrator.aligned_lambda7 import build_aligned_lambda7_family
from dihiggs.app.orchestrator.aligned_constraints import unitarity
from dihiggs.app.orchestrator.analytic_bridge import evaluate_t4

HERE = Path(__file__).resolve().parent
CAMPAIGN = HERE.parent
MODES=("legacy","explicit","analytic")
PRECISIONS=(50,80,120)
DECLARED_ACCEPTANCE_TOLERANCES={"dimensionless_lambda_Z_abs":"5e-12",
    "dimensionless_eigenvalue_abs":"5e-11","mass_GeV_abs":"1e-8",
    "coupling_GeV_abs":"1e-8","v_GeV_abs":"1e-10",
    "tan_beta_abs":"1e-8","signed_kappa_abs":"1e-12",
    "alignment_abs":"1e-12","m12_GeV2_abs":"1e-12",
    "diphoton_amplitude_rel":"1e-11","diphoton_component_abs_scale":"1e-11",
    "diphoton_width_rel":"1e-10","diphoton_width_abs_GeV":"1e-30",
    "running_mass_identity_abs_scale":"1e-12"}


def run_point(*, tan_beta="1000000", lambda1_ref="1.2", mphi="150",
              point_builder=None, precisions=(50,80,120)):
    records = []
    for precision in precisions:
        if point_builder is None:
            point = build_aligned_lambda7_family(
                tan_beta=tan_beta, lambda1_ref=lambda1_ref, mh_GeV="125.13",
                mphi_GeV=mphi, mA_GeV=str(int(mphi)+50), mHp_GeV=str(int(mphi)+50),
                v_GeV="246", X_lambda6_tanbeta="1", precision=precision).as_dict()
        else:
            constructed=point_builder(precision)
            point=constructed.as_dict() if hasattr(constructed,"as_dict") else constructed
        v, i = point["values"], point["inputs"]
        expected = {"lambda1": float(v["lambda1"]), **{
            f"lambda{k}": float(v[f"lambda{k}"]) for k in range(2,8)},
            **{f"Z{k}": float(v[f"Z{k}"]) for k in range(1,8)},
            "m12_sq_GeV2": float(v["m12_sq_GeV2"]),
            "m_h_backend_GeV": float(i["mh_GeV"]), "m_H_backend_GeV": float(i["mphi_GeV"]),
            "m_A_backend_GeV": float(i["mA_GeV"]), "m_Hp_backend_GeV": float(i["mHp_GeV"]),
            "tan_beta_backend": float(i["tan_beta"]), "backend_v_GeV": float(i["v_GeV"])}
        for mode in ("legacy", "explicit", "analytic"):
            started = time.perf_counter()
            result = evaluate_t4(point, mode=mode, executable=HERE/"analytic_bridge")
            result["declared_acceptance_tolerances"] = DECLARED_ACCEPTANCE_TOLERANCES
            result["calculator_elapsed_seconds"] = time.perf_counter()-started
            result["precision_digits_constructed"] = precision
            result["inputs_decimal"] = i | {"X_lambda6_tanbeta": v["X_lambda6_tanbeta"]}
            result["target_values_decimal"] = {key: v[key] for key in
                [*(f"lambda{k}" for k in range(1,8)), *(f"Z{k}" for k in range(1,8)),
                 "m12_sq_GeV2", "g_hphiphi_GeV", "CphiHpHm_potential_GeV"]}
            result["phase_convention_mapping"] = {
                "status": result.get("physical_state_ordering", "MAPPED"),
                "h_index": result.get("requested_h_index_backend"),
                "phi_index": result.get("requested_phi_index_backend"),
                "h_phase_factor": result.get("requested_h_phase_factor"),
                "phi_phase_factor": result.get("requested_phi_phase_factor"),
                "raw_h_vev_projection": result.get("raw_h_vev_projection"),
                "raw_phi_vev_projection": result.get("raw_phi_vev_projection"),
                "mapped_h_vev_projection": result.get("mapped_h_vev_projection"),
                "mapped_phi_vev_projection": result.get("mapped_phi_vev_projection"),
                "raw_sin_beta_minus_alpha": result.get("sin_beta_minus_alpha_raw"),
                "raw_cos_beta_minus_alpha": result.get("cos_beta_minus_alpha_raw"),
            }
            if result.get("constructor_ok"):
                diffs = {key: result[key]-target for key,target in expected.items()}
                result["backend_abs_residuals"] = {key: abs(value) for key,value in diffs.items()}
                result["backend_max_abs_lambda_residual"] = max(abs(diffs[f"lambda{k}"]) for k in range(1,8))
                result["backend_max_abs_Z_residual"] = max(abs(diffs[f"Z{k}"]) for k in range(1,8))
                result["backend_abs_residuals"]["g_hphi_phi_GeV"] = abs(result["g_hphi_phi_backend_GeV"]-float(v["g_hphiphi_GeV"]))
                result["backend_abs_residuals"]["g_phi_HpHm_GeV"] = abs(result["g_phi_HpHm_backend_GeV"]-float(v["CphiHpHm_potential_GeV"]))
                result["backend_abs_residuals"]["kappa_phi_top_signed"] = abs(result["kappa_phi_top_backend_signed"]+1/float(i["tan_beta"]))
                result["backend_abs_residuals"]["sin_beta_minus_alpha"] = abs(result["sin_beta_minus_alpha_backend"]-1.0)
                result["backend_abs_residuals"]["cos_beta_minus_alpha"] = abs(result["cos_beta_minus_alpha_backend"])
                result["signed_couplings"] = {
                    "expected_kappa_phi_top_Type_I": str(-1/float(i["tan_beta"])),
                    "backend_kappa_phi_top_signed": result["kappa_phi_top_backend_signed"],
                    "expected_CV_phiHpHm_GeV": v["CphiHpHm_potential_GeV"],
                    "backend_CV_phiHpHm_GeV": result["g_phi_HpHm_backend_GeV"],
                    "alignment_backend_sin": result["sin_beta_minus_alpha_backend"],
                    "alignment_backend_cos": result["cos_beta_minus_alpha_backend"]}
                if result.get("widths_computed",True):
                    channels = {key: result[f"width_{key}_GeV"] for key in
                                ("bb","cc","tautau","mumu","gg","diphoton","Zgamma","WW","ZZ","hh")}
                    result["width_channels_GeV"] = channels
                    total = result["width_total_phi_GeV"]
                    result["branching_fractions_reported_channels"] = {
                        key: (value/total if total else None) for key,value in channels.items()}
                    result["physical_lifetime"] = {"ctau_mm": result["ctau_mm"],
                        "status": "INCONCLUSIVE", "reason": "non-diphoton decay channels remain calculator-reported, without an independent all-channel width audit"}
                    result["decay_status"] = {"widths": "INCONCLUSIVE", "branching_fractions": "INCONCLUSIVE",
                        "displaced_photon_selection": "NOT_TESTED", "experimental_exclusions": "NOT_TESTED"}
                    result["appendix_diphoton_amplitude"] = _appendix_amplitude(v, i, result)
                    result["independent_diphoton_components"] = _independent_diphoton(result)
                    native_audit=result["independent_diphoton_components"]
                    appendix=result["appendix_diphoton_amplitude"]
                    if native_audit.get("status") in ("PASS","FAIL") and appendix.get("appendix_LO_width_GeV"):
                        appendix_width=mp.mpf(appendix["appendix_LO_width_GeV"])
                        native_width=mp.mpf(native_audit["independent_width_GeV"])
                        native_audit["appendix_LO_vs_running_mass_QCD_width"]={
                            "appendix_LO_pole_scheme_width_GeV":mp.nstr(appendix_width,50),
                            "independent_running_mass_QCD_width_GeV":mp.nstr(native_width,50),
                            "absolute_difference_GeV":mp.nstr(abs(appendix_width-native_width),50),
                            "relative_difference":mp.nstr(abs(appendix_width-native_width)/max(abs(appendix_width),abs(native_width),mp.mpf("1e-100")),50),
                            "interpretation":"scheme correction diagnostic; independent native-scheme reconstruction is compared separately against DecayTable components"}
                else:
                    result["constructor_status"]="INCONCLUSIVE"
                    result["exact_failure_reason"]=result.get("exact_failure_reason") or "CP-even state ordering is unresolved"
                    result["widths"]="NOT_TESTED"
                    result["branching_fractions"]="NOT_TESTED"
                    result["lifetime"]="NOT_TESTED"
                backend_lambdas = {f"lambda{k}": str(result[f"lambda{k}"]) for k in range(1,8)}
                independent_spectrum = unitarity(backend_lambdas, precision=80)
                spectrum_differences = {}
                for block, eigenvalues in independent_spectrum["eigenvalues"].items():
                    backend_values = result[f"backend_eigenvalues_{block}"]
                    spectrum_differences[block] = max(abs(float(x)-y) for x,y in zip(eigenvalues,backend_values))
                result["independent_full_unitarity"] = independent_spectrum
                result["backend_spectrum_max_abs_eigenvalue_residual"] = max(spectrum_differences.values())
                result["backend_spectrum_block_residuals"] = spectrum_differences
                result["double_backend_inputs"] = {key: float(value) for key,value in v.items()
                    if key.startswith("lambda") or key.startswith("Z") or key == "m12_sq_GeV2"}
                from decimal import Decimal
                result["input_serialization_abs_loss"] = {
                    key: str(abs(Decimal.from_float(float(value))-Decimal(str(value))))
                    for key,value in {**v,**i}.items()
                    if isinstance(value,(int,float,str)) and str(value).replace(".","",1).replace("-","",1).replace("+","",1).replace("E","",1).replace("e","",1).isdigit()}
                tolerances=DECLARED_ACCEPTANCE_TOLERANCES
                checks={}
                for key,residual in result["backend_abs_residuals"].items():
                    if key.startswith("lambda") or key.startswith("Z"):
                        tol=tolerances["dimensionless_lambda_Z_abs"]
                    elif "eigenvalue" in key:
                        tol=tolerances["dimensionless_eigenvalue_abs"]
                    elif key.startswith("m_") and key.endswith("_GeV"):
                        tol=tolerances["mass_GeV_abs"]
                    elif key.startswith("g_"):
                        tol=tolerances["coupling_GeV_abs"]
                    elif key in ("backend_v_GeV",):
                        tol=tolerances["v_GeV_abs"]
                    elif key=="tan_beta_backend":
                        tol=tolerances["tan_beta_abs"]
                    elif key=="m12_sq_GeV2":
                        tol=tolerances["m12_GeV2_abs"]
                    elif key=="kappa_phi_top_signed":
                        tol=tolerances["signed_kappa_abs"]
                    else:
                        tol=tolerances["alignment_abs"]
                    status="PASS" if residual<float(tol) else "FAIL"
                    reason=""
                    checks[key]={"status":status,"absolute_residual":residual,"declared_tolerance":tol,
                                 "reason":reason}
                checks["complete_unitarity_spectrum"]={"status":"PASS" if result["backend_spectrum_max_abs_eigenvalue_residual"] < float(tolerances["dimensionless_eigenvalue_abs"]) else "FAIL",
                    "absolute_residual":result["backend_spectrum_max_abs_eigenvalue_residual"],"declared_tolerance":tolerances["dimensionless_eigenvalue_abs"]}
                result["constructor_checks"]=checks
                statuses={item["status"] for item in checks.values()}
                result["constructor_status"]=("INCONCLUSIVE" if result.get("physical_state_ordering")=="INCONCLUSIVE" else "FAIL" if "FAIL" in statuses else
                    "INCONCLUSIVE" if "INCONCLUSIVE" in statuses else "PASS")
                result["backend_constructor_status"]="PASS"
                result["family_constructor_agreement_status"]=result["constructor_status"]
                result["exact_failure_reason"]=(result.get("exact_failure_reason","") if result.get("physical_state_ordering")=="INCONCLUSIVE" else "; ".join(
                    key+(f": {item['reason']}" if item.get("reason") else "")
                    for key,item in checks.items() if item["status"] in ("FAIL","INCONCLUSIVE")))
            else:
                if result.get("status")=="INCONCLUSIVE":
                    result["constructor_status"]="INCONCLUSIVE"
                    result.setdefault("exact_failure_reason","required downstream calculator observables unavailable")
                else:
                    result["constructor_status"] = "FAIL"
                    result["backend_constructor_status"] = "FAIL"
                    result["family_constructor_agreement_status"] = "FAIL"
                    result["exact_failure_reason"] = "2HDMC constructor returned false"
            records.append(result)
    return records


def _loop_f(tau):
    tau = mp.mpf(tau)
    if tau <= 1:
        return mp.asin(mp.sqrt(tau))**2
    root=mp.sqrt(1-1/tau)
    return -mp.mpf("0.25")*(mp.log((1+root)/(1-root))-1j*mp.pi)**2


def _f_sf(tau):
    return (1/tau)*(1+(1-1/tau)*_loop_f(tau))


def _f_pf(tau):
    return (1/tau)*_loop_f(tau)


def _f0(tau):
    return (1/tau)*(-1+(1/tau)*_loop_f(tau))


def _f1(tau):
    return 2+3/tau+3/tau*(2-1/tau)*_loop_f(tau)


def _complex_json(value):
    return mp.mpc(str(value["real"]),str(value["imag"]))


def _independent_diphoton(backend):
    """Rebuild 2HDMC's photon amplitude from exported masses and vertices."""
    native=backend.get("diphoton_native_diagnostics")
    if not native:
        return {"status":"NOT_TESTED","reason":"native DecayTable amplitude diagnostics unavailable"}
    with mp.workdps(50):
        phase=mp.mpf(native["state_phase"])
        mass=mp.mpf(str(native["mass_GeV"])); v=mp.mpf(str(native["vev_GeV"]))
        alpha=mp.mpf(str(native["alpha0"])); mW=mp.mpf(str(native["mW_GeV"]))
        mHp=mp.mpf(str(native["mHp_GeV"])); scalar=mp.mpc(0); pseudoscalar=mp.mpc(0)
        amplitude_scale=max(mp.mpf(1),abs(phase*_complex_json(native["scalar_amplitude_raw"])),
            abs(phase*_complex_json(native["pseudoscalar_amplitude_raw"])))
        entries=[]; component_checks=[]; max_mass_identity=mp.mpf(0); scalar_before_qcd=mp.mpf(0)
        def component_check(name,independent,reference):
            difference=abs(independent-reference)
            normalized=None if reference==0 else difference/abs(reference)
            # The exported backend is binary64: compare every component on the
            # declared dimensionless amplitude scale, while also retaining its
            # individual normalized residual for diagnosis. Exact-zero
            # references are evaluated by absolute residual only.
            ok=difference<=mp.mpf(DECLARED_ACCEPTANCE_TOLERANCES["diphoton_component_abs_scale"])*amplitude_scale
            component_checks.append({"component":name,"status":"PASS" if ok else "FAIL",
                "absolute_residual":mp.nstr(difference,50),
                "normalized_residual":None if normalized is None else mp.nstr(normalized,50),
                "comparison_scale":mp.nstr(amplitude_scale,50),
                "declared_absolute_tolerance":mp.nstr(mp.mpf(DECLARED_ACCEPTANCE_TOLERANCES["diphoton_component_abs_scale"])*amplitude_scale,50),
                "zero_reference_uses_absolute_tolerance":reference==0})
            return ok
        components_ok=True
        for fermion in native["fermions"]:
            if not fermion["included"]: continue
            fam=fermion["family"]
            charge={"down":-mp.mpf(1)/3,"up":mp.mpf(2)/3,"lepton":-1}[fam]
            color=3 if fam in ("down","up") else 1
            mhd=mp.mpf(str(fermion["heavy_damping_mass_GeV"]))
            mrun=mp.mpf(str(fermion["vertex_running_mass_GeV"]))
            tau=mass**2/(4*mhd**2)
            gs=phase*_complex_json(fermion["scalar_coupling"])
            gp=phase*_complex_json(fermion["pseudoscalar_coupling"])
            sf=2*color*charge**2*gs*v/mrun*_f_sf(tau)
            pf=2*color*charge**2*gp*v/mrun*_f_pf(tau)
            sf_after=sf
            scalar_before_qcd+=sf
            if fam=="up" and fermion["generation"]==3 and not native["heavy_damping_enabled"] and native["qcd_enabled"]:
                asval=mp.mpf(str(fermion["alpha_s"]))
                sf_after*=1-asval/mp.pi
            scalar+=sf_after; pseudoscalar+=pf
            native_sf=phase*_complex_json(fermion["scalar_amplitude_after_qcd"])
            native_pf=phase*_complex_json(fermion["pseudoscalar_amplitude_after_qcd"])
            components_ok &= component_check(f"{fam}{fermion['generation']}_scalar",sf_after,native_sf)
            components_ok &= component_check(f"{fam}{fermion['generation']}_pseudoscalar",pf,native_pf)
            if fam != "lepton":
                reconstructed=mhd  # populated below from the exported running masses
                mhalf=mp.mpf(str(fermion["running_mass_at_half_scale_GeV"]))
                mpole=mp.mpf(str(fermion["pole_mass_GeV"]))
                mpole_run=mp.mpf(str(fermion["running_mass_at_pole_scale_GeV"]))
                reconstructed=mhalf*mpole/mpole_run
                max_mass_identity=max(max_mass_identity,abs(reconstructed-mhd))
            entries.append({"family":fam,"generation":fermion["generation"],
                "scalar_amplitude_independent":mp.nstr(sf_after,50),
                "pseudoscalar_amplitude_independent":mp.nstr(pf,50),
                "scalar_amplitude_native_mapped":mp.nstr(phase*_complex_json(fermion["scalar_amplitude_after_qcd"]),50),
                "pseudoscalar_amplitude_native_mapped":mp.nstr(phase*_complex_json(fermion["pseudoscalar_amplitude_after_qcd"]),50)})
        w_coupling=phase*_complex_json(native["W_coupling_raw"])
        hp_coupling=phase*_complex_json(native["charged_higgs_coupling_raw"])
        w_amp=w_coupling*v/(2*mW**2)*_f1(mass**2/(4*mW**2))
        hp_amp=hp_coupling/v*v**2/(2*mHp**2)*_f0(mass**2/(4*mHp**2))
        components_ok &= component_check("W_scalar",w_amp,phase*_complex_json(native["W_amplitude_raw"]))
        components_ok &= component_check("charged_Higgs_scalar",hp_amp,
            phase*_complex_json(native["charged_higgs_amplitude_raw"]))
        scalar_no_qcd=scalar_before_qcd+w_amp+hp_amp
        scalar+=w_amp+hp_amp
        native_scalar=phase*_complex_json(native["scalar_amplitude_raw"])
        native_pseudo=phase*_complex_json(native["pseudoscalar_amplitude_raw"])
        width=mass**3*alpha**2/(256*mp.pi**3*v**2)*(abs(scalar)**2+abs(pseudoscalar)**2)
        width_no_qcd=mass**3*alpha**2/(256*mp.pi**3*v**2)*(abs(scalar_no_qcd)**2+abs(pseudoscalar)**2)
        native_width=mp.mpf(str(native["loop_width_GeV"]))
        scalar_diff=abs(scalar-native_scalar); pseudo_diff=abs(pseudoscalar-native_pseudo)
        comp=max(scalar_diff,pseudo_diff)
        scalar_norm=scalar_diff/abs(native_scalar) if native_scalar else None
        pseudo_norm=pseudo_diff/abs(native_pseudo) if native_pseudo else None
        width_diff=abs(width-native_width)
        production_width=mp.mpf(str(backend["width_diphoton_GeV"]))
        production_width_diff=abs(production_width-native_width)
        production_width_ok=(production_width_diff<=mp.mpf(DECLARED_ACCEPTANCE_TOLERANCES["diphoton_width_abs_GeV"]) if native_width==0 else
            production_width_diff/abs(native_width)<=mp.mpf(DECLARED_ACCEPTANCE_TOLERANCES["diphoton_width_rel"]))
        production_width_status=("NOT_TESTED" if not native["qcd_enabled"] else
            "PASS" if production_width_ok else "FAIL")
        amplitude_ok=all((diff<=mp.mpf(DECLARED_ACCEPTANCE_TOLERANCES["diphoton_component_abs_scale"])*amplitude_scale if reference==0 else
            diff/abs(reference)<=mp.mpf(DECLARED_ACCEPTANCE_TOLERANCES["diphoton_amplitude_rel"]))
            for diff,reference in ((scalar_diff,native_scalar),(pseudo_diff,native_pseudo)))
        width_ok=(width_diff<=mp.mpf(DECLARED_ACCEPTANCE_TOLERANCES["diphoton_width_abs_GeV"]) if native_width==0 else
            width_diff/abs(native_width)<=mp.mpf(DECLARED_ACCEPTANCE_TOLERANCES["diphoton_width_rel"]))
        mass_ok=max_mass_identity<=mp.mpf(DECLARED_ACCEPTANCE_TOLERANCES["running_mass_identity_abs_scale"])*max(1,mass)
        components_ok &= production_width_status in ("PASS","NOT_TESTED")
        return {"status":"PASS" if amplitude_ok and width_ok and components_ok and mass_ok else "FAIL",
            "method":"independent mpmath implementation of DecayTable Fsf/Fpf/F0/F1 and exported running-mass/QCD scheme",
            "precision_digits":50,"declared_diagnostic_tolerance":"amplitude <= 1e-11*max(1,|S|,|P|); width <= max(1e-30 GeV,1e-10*Gamma); running-mass identity <= 1e-12*max(1,M)",
            "independent_scalar_amplitude":mp.nstr(scalar,50),"independent_pseudoscalar_amplitude":mp.nstr(pseudoscalar,50),
            "independent_scalar_amplitude_before_top_qcd":mp.nstr(scalar_no_qcd,50),
            "independent_width_before_top_qcd_GeV":mp.nstr(width_no_qcd,50),
            "top_qcd_width_shift_GeV":mp.nstr(width-width_no_qcd,50),
            "native_scalar_amplitude_mapped":mp.nstr(native_scalar,50),"native_pseudoscalar_amplitude_mapped":mp.nstr(native_pseudo,50),
            "max_amplitude_component_abs_residual":mp.nstr(comp,50),
            "amplitude_component_checks":component_checks,
            "production_width_GeV":mp.nstr(production_width,50),
            "production_width_absolute_residual_GeV":mp.nstr(production_width_diff,50),
            "production_width_normalized_residual":(None if native_width==0 else mp.nstr(production_width_diff/abs(native_width),50)),
            "production_width_channel_status":production_width_status,
            "production_width_channel_note":("QCD-off preserves 2HDMC's legacy disabled width channel; compare the independently reconstructed loop width against loop_width_GeV, not the channel's zero return"
                if not native["qcd_enabled"] else "production channel is compared with the native diagnostic loop width"),
            "scalar_amplitude_normalized_residual":None if scalar_norm is None else mp.nstr(scalar_norm,50),
            "pseudoscalar_amplitude_normalized_residual":None if pseudo_norm is None else mp.nstr(pseudo_norm,50),
            "independent_width_GeV":mp.nstr(width,50),
            "native_width_GeV":mp.nstr(native_width,50),"absolute_width_residual_GeV":mp.nstr(width_diff,50),
            "max_running_mass_reconstruction_residual_GeV":mp.nstr(max_mass_identity,50),
            "fermion_components":entries,"W_amplitude_independent":mp.nstr(w_amp,50),
            "charged_higgs_amplitude_independent":mp.nstr(hp_amp,50),
            "qcd_on":native["qcd_enabled"],"heavy_damping_on":native["heavy_damping_enabled"],
            "top_qcd_component_factor":(mp.nstr(1-mp.mpf(str(native["alpha_s_at_mass"]))/mp.pi,50)
                if native["qcd_enabled"] and not native["heavy_damping_enabled"] else "1"),
            "native_W_amplitude_mapped":mp.nstr(phase*_complex_json(native["W_amplitude_raw"]),50),
            "native_charged_higgs_amplitude_mapped":mp.nstr(phase*_complex_json(native["charged_higgs_amplitude_raw"]),50)}


def _appendix_amplitude(values, inputs, backend):
    """Appendix Eq. phi_gammagamma_alignment, leading-order pole-mass check."""
    with mp.workdps(50):
        mass=mp.mpf(inputs["mphi_GeV"]); v=mp.mpf(inputs["v_GeV"])
        tb=mp.mpf(inputs["tan_beta"]); mhp=mp.mpf(inputs["mHp_GeV"])
        z7=mp.mpf(values["Z7"])
        fermions=((172.5,3,mp.mpf(2)/3),(4.75,3,mp.mpf(-1)/3),
                  (1.42,3,mp.mpf(2)/3),(1.77686,1,-1),(0.105658,1,-1))
        fermion_sum=mp.mpc(0)
        for mf,nc,charge in fermions:
            tau=mass**2/(4*mp.mpf(str(mf))**2)
            fermion_sum += nc*charge**2*2*(tau+(tau-1)*_loop_f(tau))/(tau**2)
        tau_h=mass**2/(4*mhp**2)
        A0=-(tau_h-_loop_f(tau_h))/(tau_h**2)
        kappa_f=-1/tb
        scalar_coefficient=-v**2*z7/(2*mhp**2)
        amplitude=kappa_f*fermion_sum+scalar_coefficient*A0
        alpha0=mp.mpf(1)/mp.mpf("137.0359997")
        width=alpha0**2*mass**3/(256*mp.pi**3*v**2)*abs(amplitude)**2
        calc=mp.mpf(str(backend["width_diphoton_GeV"]))
        return {"status":"INCONCLUSIVE", "method":"independent leading-order appendix amplitude, Type-I exact-alignment modifiers, pole-mass loop functions",
            "observable":"phi -> gamma gamma reduced amplitude and partial width", "tolerance":"no frozen QCD/running-mass correction budget; diagnostic comparison only",
            "precision":"mpmath 50 decimal digits; decimal pole masses and appendix coupling coefficients",
            "kappa_fermions_signed":mp.nstr(kappa_f,50), "kappa_V_signed":"0",
            "C_V_phiHpHm_GeV":mp.nstr(-v*z7,50), "fermion_loop_sum":mp.nstr(fermion_sum,50),
            "charged_scalar_A0":mp.nstr(A0,50), "charged_scalar_reduced_amplitude":mp.nstr(scalar_coefficient*A0,50),
            "reduced_amplitude":mp.nstr(amplitude,50), "appendix_LO_width_GeV":mp.nstr(width,50),
            "calculator_width_GeV":mp.nstr(calc,50), "absolute_width_difference_GeV":mp.nstr(abs(width-calc),50),
            "relative_width_difference":mp.nstr(abs(width-calc)/max(abs(width),abs(calc),mp.mpf("1e-100")),50),
            "evidence":"2HDMC uses running-mass and QCD conventions; scheme corrections have not been decomposed"}


if __name__ == "__main__":
    rows = run_point()
    output = {"status": "PARTIAL", "scope": "three actual calculator constructor routes, high-tanbeta cancellation anchor",
              "precision_constructions": [50,80,120], "route_count": len(rows),
              "unitarity_backend_note": "bridge exports backend GSL spectra for all four GI blocks and calc_unitarity maximum; mpmath spectra are compared block by block",
              "production_dispatch": {"status": "PARTIAL", "reason": "new isolated orchestration API implements legacy-default/explicit/analytic routes; existing production CLI and legacy schemas are not modified"},
              "diphoton_amplitude_crosscheck": {"status": "INCONCLUSIVE", "reason": "leading-order appendix amplitude is computed; scheme-level QCD/running-mass comparison remains unresolved"},
              "displaced_photon_selection": {"status": "NOT_TESTED", "reason": "no applicable validated selection configured"},
              "experimental_exclusions": {"status": "NOT_TESTED", "reason": "no validated method or production normalization in scope"},
              "routes": rows}
    from collections import Counter
    output["constructor_status_counts"]={mode:dict(Counter(x["constructor_status"] for x in rows if x["mode"]==mode))
        for mode in ("legacy","explicit","analytic")}
    # Freeze measured anchor-level budgets before any future boundary refinement.
    # These are measured numerical error components, separate from construction
    # residuals and acceptance tolerances.
    measured={}
    for mode in ("legacy","explicit","analytic"):
        group=[x for x in rows if x["mode"]==mode and x.get("constructor_ok")]
        if not group:
            measured[mode]={"status":"NOT_MEASURED"}
            continue
        serialization_keys=set.intersection(*(set(x["input_serialization_abs_loss"]) for x in group))
        measured[mode]={"status":"FROZEN_SINGLE_ANCHOR","precision_digits":[50,80,120],
            "max_input_decimal_to_binary64_loss_by_observable":{
                key:max(float(x["input_serialization_abs_loss"][key]) for x in group) for key in sorted(serialization_keys)},
            "max_observed_precision_spread_by_backend_output":{
                key:max(x[key] for x in group)-min(x[key] for x in group)
                for key in ("backend_v_GeV","lambda1","lambda2","lambda3","lambda4","lambda5","lambda6","lambda7",
                            "Z1","Z2","Z3","Z4","Z5","Z6","Z7","m12_sq_GeV2","m_h_backend_GeV","m_H_backend_GeV",
                            "m_A_backend_GeV","m_Hp_backend_GeV","g_hphi_phi_backend_GeV","g_phi_HpHm_backend_GeV")
                if all(key in x for x in group)},
            "max_eigensolver_error_by_scattering_block":{
                block:max(float(x["independent_full_unitarity"]["measured_eigensolver_error_components"][
                    "max_eigenpair_residual"]) for x in group)
                for block in group[0]["independent_full_unitarity"]["eigenvalues"]}}
    output["measured_numerical_error_budgets"]={"scope":"one high-t anchor only; frozen before any refinement; construction residuals are recorded separately",
        "components":"decimal-to-binary64 input serialization loss, measured backend output spread over 50/80/120 digit Decimal constructions, and independent scattering eigensolver residuals",
        "by_constructor_route":measured}
    (CAMPAIGN / "T4_integration_report.json").write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"output": str(CAMPAIGN / "T4_integration_report.json"),
                      "constructors_ok": sum(bool(row.get("constructor_ok")) for row in rows),
                      "routes": len(rows)}, indent=2))
