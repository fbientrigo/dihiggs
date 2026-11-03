"""Isolated 2HDMC backend API for legacy, explicit, and analytic construction.

``legacy`` remains the default. The API does not change the historical CLI or
its output schema; it calls a small campaign-owned C++ bridge into the same
THDM/SM/DecayTable library.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any
from decimal import Decimal, InvalidOperation


_MODES = {"legacy": "THDM::set_param_phys",
          "explicit": "THDM::set_param_phys_lam1",
          "analytic": "THDM::set_param_gen_exact (delegates to set_param_gen)"}


def _finite(name: str, raw: Any, *, positive: bool=False) -> str:
    try:
        value=Decimal(str(raw))
    except (InvalidOperation,ValueError,TypeError) as exc:
        raise ValueError(f"{name} must be a finite decimal number") from exc
    if not value.is_finite() or (positive and value<=0):
        rule="finite and > 0" if positive else "finite"
        raise ValueError(f"{name} must be {rule}")
    return str(value)


def validate_family_point(point: Any) -> dict[str, Any]:
    """Rebuild a serialized stable-family point and reject altered derived fields."""
    from dihiggs.app.orchestrator.aligned_lambda7 import build_aligned_lambda7_family

    data=point.as_dict() if hasattr(point,"as_dict") else point
    if not isinstance(data,dict) or not isinstance(data.get("inputs"),dict) or not isinstance(data.get("values"),dict):
        raise ValueError("point must contain input and value mappings from the aligned-family constructor")
    inputs=data["inputs"]
    provenance=data.get("provenance")
    if not isinstance(provenance,dict) or provenance.get("decimal_inputs_preserved") is not True:
        raise ValueError("point provenance must certify preserved decimal-string inputs")
    for key,value in inputs.items():
        if not isinstance(value,str):
            raise ValueError(f"point input {key} must remain a decimal string")
    for key,value in data["values"].items():
        if not isinstance(value,str):
            raise ValueError(f"point derived value {key} must remain a decimal string")
    common={key:inputs[key] for key in ("tan_beta","lambda1_ref","mh_GeV","mphi_GeV","mA_GeV","mHp_GeV","v_GeV")}
    precision=int(provenance.get("precision_digits",80))
    if precision<20 or precision>500:
        raise ValueError("point provenance precision_digits must be in [20, 500]")
    lambda6=inputs.get("lambda6_input")
    X=inputs.get("X_lambda6_tanbeta_input")
    trajectory=str(data.get("trajectory",""))
    fixed_l6_label=trajectory=="fixed_lambda6" or trajectory.startswith("fixed_lambda6_")
    fixed_X_label=trajectory=="fixed_X" or trajectory.startswith("fixed_X_") or trajectory=="Xstar"
    if not (fixed_l6_label or fixed_X_label):
        raise ValueError(f"unsupported trajectory label {trajectory!r}")
    if lambda6 not in (None,"") and X not in (None,""):
        raise ValueError("point cannot declare both fixed lambda6 and fixed X")
    if lambda6 not in (None,""):
        if not fixed_l6_label:
            raise ValueError("trajectory label does not match fixed lambda6 input")
        rebuilt=build_aligned_lambda7_family(**common,lambda6=lambda6,precision=precision)
    elif X not in (None,""):
        if not fixed_X_label:
            raise ValueError("trajectory label does not match fixed X input")
        rebuilt=build_aligned_lambda7_family(**common,X_lambda6_tanbeta=X,precision=precision)
        if trajectory=="Xstar":
            from dihiggs.app.orchestrator.aligned_lambda7 import xstar_lambda6_tanbeta
            expected_X=xstar_lambda6_tanbeta(**{key:common[key] for key in
                ("tan_beta","lambda1_ref","mh_GeV","mphi_GeV","v_GeV")},precision=precision)
            if Decimal(X)!=Decimal(expected_X):
                raise ValueError("Xstar trajectory X does not match its evaluated Xstar formula")
    else:
        raise ValueError("point must preserve either lambda6_input or X_lambda6_tanbeta_input")
    expected=rebuilt.as_dict()
    supplied=data["values"]
    for key in (*(f"lambda{i}" for i in range(1,8)),*(f"Z{i}" for i in range(1,8)),
                "m12_sq_GeV2","g_hphiphi_GeV","CphiHpHm_potential_GeV","X_lambda6_tanbeta"):
        if key not in supplied or Decimal(str(supplied[key])) != Decimal(str(expected["values"][key])):
            raise ValueError(f"serialized derived field {key} does not match constructor inputs")
    expected["trajectory"] = trajectory
    return expected


def evaluate_t4(point: Any, *, mode: str = "legacy", executable: str | Path | None = None,
                qcd_on: bool = True) -> dict[str, Any]:
    """Evaluate one complete stable-family point through an actual 2HDMC route.

    Inputs to the double backend are serialized decimal strings from the point
    and are retained in the result. ``explicit`` supplies lambda1/lambda6/lambda7
    to the calculator's lambda1-aware physical constructor and matches v;
    ``analytic`` sends
    all stable quartics directly through ``set_param_gen`` and matches the
    requested v with ``SM::set_GF``. ``legacy`` is the unchanged default
    physical constructor with its default SM parameters.
    """
    if mode not in _MODES:
        raise ValueError(f"mode must be one of {', '.join(_MODES)}")
    data = validate_family_point(point)
    inputs, values = data["inputs"], data["values"]
    clean_inputs={key:_finite(key,inputs[key],positive=True) for key in
                  ("tan_beta","mh_GeV","mphi_GeV","mA_GeV","mHp_GeV","v_GeV")}
    if mode != "analytic" and Decimal(clean_inputs["mphi_GeV"]) < Decimal(clean_inputs["mh_GeV"]):
        raise ValueError("physical constructor requires mphi_GeV >= mh_GeV")
    clean_lambdas={f"lambda{i}":_finite(f"lambda{i}",values[f"lambda{i}"]) for i in range(1,8)}
    clean_m12=_finite("m12_sq_GeV2",values["m12_sq_GeV2"])
    if Decimal(clean_inputs["mh_GeV"]) == Decimal(clean_inputs["mphi_GeV"]):
        reason="CP-even target masses are exactly degenerate, so the requested h/phi eigenstates and signed couplings are not identifiable"
        return {"mode":mode,"requested_mode":mode,"constructor_route":_MODES[mode],
            "constructor_ok":False,"backend_constructor_status":"NOT_TESTED",
            "family_constructor_agreement_status":"NOT_TESTED",
            "constructor_status":"INCONCLUSIVE","status":"INCONCLUSIVE",
            "physical_state_ordering":"INCONCLUSIVE","exact_failure_reason":reason,
            "widths":"NOT_TESTED","branching_fractions":"NOT_TESTED","lifetime":"NOT_TESTED",
            "displaced_photon_selection":"NOT_TESTED","experimental_exclusions":"NOT_TESTED",
            "diphoton_width_channel_status":"NOT_TESTED",
            "diphoton_width_channel_reason":"CP-even eigenstate is degenerate; decay observables are not assigned",
            "inputs_decimal":dict(inputs),"values_decimal":dict(values),
            "double_backend_decimal_inputs":{"inputs":dict(inputs),
                "lambdas":clean_lambdas,
                "Z":{f"Z{i}":values[f"Z{i}"] for i in range(1,8)},"m12_sq_GeV2":clean_m12},
            "requested_qcd_on":bool(qcd_on)}
    if executable is None:
        repo = Path(__file__).resolve().parents[3]
        executable = repo / "2hdmc/analytic-bridge"
    executable = Path(executable)
    if not executable.is_file():
        raise FileNotFoundError(f"2HDMC bridge executable not found: {executable}; build it with `make -C {executable.parent} analytic-bridge`")
    route_lambdas=dict(clean_lambdas)
    if mode=="legacy":
        route_lambdas["lambda7"]="0"
    argv = [str(executable), mode, clean_inputs["mh_GeV"], clean_inputs["mphi_GeV"], clean_inputs["mA_GeV"],
            clean_inputs["mHp_GeV"], clean_inputs["tan_beta"], *(route_lambdas[f"lambda{i}"] for i in range(1,8)),
            clean_m12, clean_inputs["v_GeV"], "1" if qcd_on else "0"]
    done = subprocess.run(argv, check=True, text=True, capture_output=True, timeout=60)
    result = json.loads(done.stdout)
    result["mode"] = mode
    result["requested_mode"] = mode
    result["constructor_route"] = _MODES[mode]
    result["inputs_decimal"] = dict(inputs)
    result["values_decimal"] = dict(values)
    result["requested_qcd_on"] = bool(qcd_on)
    if not result.get("constructor_ok"):
        result["backend_constructor_status"] = "FAIL"
        result["family_constructor_agreement_status"] = "FAIL"
        result["constructor_status"] = "FAIL"
        result["exact_failure_reason"] = f"{_MODES[mode]} returned false for the validated finite inputs"
    elif result.get("physical_state_ordering") == "INCONCLUSIVE":
        result["backend_constructor_status"] = "PASS"
        result["family_constructor_agreement_status"] = "NOT_TESTED"
        result["constructor_status"] = "INCONCLUSIVE"
        result["widths"] = "NOT_TESTED"
        result["branching_fractions"] = "NOT_TESTED"
        result["lifetime"] = "NOT_TESTED"
    else:
        result["backend_constructor_status"] = "PASS"
        result["family_constructor_agreement_status"] = "NOT_TESTED"
        result["constructor_status"] = "INCONCLUSIVE"
        result["exact_failure_reason"] = (
            "backend constructor succeeded and state mapping is resolved; per-observable family agreement "
            "has not been evaluated by the single-point API"
        )
        channels={key:result[f"width_{key}_GeV"] for key in
            ("bb","cc","tautau","mumu","gg","diphoton","Zgamma","WW","ZZ","hh")}
        total=result["width_total_phi_GeV"]
        result["width_channels_GeV"]=channels
        result["branching_fractions_reported_channels"]={key:(value/total if total else None)
            for key,value in channels.items()}
        result["physical_lifetime"]={"ctau_mm":result["ctau_mm"],
            "status":"INCONCLUSIVE","reason":"non-diphoton decay channels remain calculator-reported without an independent all-channel audit"}
    result["diphoton_width_channel_status"] = (
        "PASS" if qcd_on and result.get("backend_constructor_status") == "PASS" and result.get("widths_computed") else "NOT_TESTED"
    )
    if not qcd_on:
        result["diphoton_width_channel_reason"] = (
            "2HDMC preserves legacy behavior and disables the production diphoton width channel when QCD is off; "
            "the separate native diagnostics still report the loop amplitude"
        )
    result["double_backend_decimal_inputs"] = {"inputs": dict(inputs),
        "lambdas": route_lambdas,
        "Z": {f"Z{i}": values[f"Z{i}"] for i in range(1,8)},
        "m12_sq_GeV2": values["m12_sq_GeV2"]}
    return result
