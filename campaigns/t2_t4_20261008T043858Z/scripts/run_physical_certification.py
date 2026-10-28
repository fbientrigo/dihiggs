#!/usr/bin/env python3
"""Resumable exact-target and rounded-calculator physical certification pilot."""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from contextlib import nullcontext
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import signal
import sys
import subprocess
import time
from typing import Any

# Set thread caps before importing SymPy, mpmath, or the calculator bridge.
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

HERE = Path(__file__).resolve().parent
CAMPAIGN = HERE.parent
REPO = CAMPAIGN.parents[1]
OUT = CAMPAIGN / "physical_certification"
RESULTS = OUT / "physical_certification.jsonl"
CHECKPOINT = OUT / "physical_certification_checkpoint.json"
LOCK = CAMPAIGN / "physical_certification.lock"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))

import sympy as sp  # noqa: E402

from dihiggs.app.orchestrator.aligned_lambda7 import build_aligned_lambda7_family  # noqa: E402
from dihiggs.app.orchestrator.aligned_constraints import _analytic_family_Z  # noqa: E402
from dihiggs.app.orchestrator.physical_certification import (  # noqa: E402
    RealPotential,
    _pi_interval,
    certify_bfb,
    certify_global_minimum,
    certify_marginal_full_stability,
    certify_tree_unitarity,
    enumerate_stationary_candidates,
    rounded_input_potential,
)
import run_t4_scan as frozen_t4  # noqa: E402
from scan_lock import campaign_scan_lock  # noqa: E402

SCHEMA = "dihiggs.physical_certification.v1"
METHOD_VERSION = "physical-pilot-exact-target-rounded-binary64-v1"
SM_GET_V2_TEXT = "60515.999999999993"
SM_GET_V2_BINARY64 = float(SM_GET_V2_TEXT)
MAX_POINTS_PER_BATCH = 50
CPU_SECONDS = 900
ADDRESS_SPACE_BYTES = 9 * 1024**3
PRECISION_DIGITS = 80


def canonical_sha(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                            ensure_ascii=True).encode()).hexdigest()


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def q(value: Any) -> sp.Rational:
    return sp.Rational(str(value))


def _source_rows() -> list[dict[str, Any]]:
    # Reuse the frozen T4 pilot's row normalizer so point IDs, convention tags,
    # and exact row fingerprints use exactly the calculator-scan semantics.
    rows = frozen_t4.scan_rows()
    if len(rows) != 390:
        raise ValueError(f"expected 390 source pilot rows (300 T1 + 90 active), got {len(rows)}")
    keys = [(str(row["point_id"]), str(row["mass_convention"])) for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError("pilot grid contains duplicate point/convention keys")
    conventions = Counter(row["mass_convention"] for row in rows)
    if conventions != Counter({"mh125.13_T1": 300, "mh125.20_active": 90}):
        raise ValueError(f"mass-convention slices are mixed or incomplete: {conventions}")
    return rows


def _frozen_t4(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    identity = frozen_t4.make_scan_identity(rows)
    checkpoint_path = CAMPAIGN / "calculator_scan_checkpoint.json"
    checkpoint = json.loads(checkpoint_path.read_text())
    if checkpoint.get("status") != "COMPLETE":
        raise ValueError("frozen T4 calculator scan is not COMPLETE")
    if checkpoint.get("scan_identity") != identity["fingerprint"] or checkpoint.get("identity") != identity:
        raise ValueError("frozen T4 calculator identity no longer matches its source/grid files")
    output_path = CAMPAIGN / "calculator_scan_agreement.jsonl"
    done = frozen_t4.load_existing(rows, identity, output_path, checkpoint_path)
    expected = len(rows) * len(frozen_t4.PRECISIONS) * len(frozen_t4.MODES)
    if len(done) != expected or checkpoint.get("completed_route_precision_records") != expected:
        raise ValueError(f"frozen T4 output is incomplete: {len(done)}/{expected} records")
    selected: dict[str, dict[str, Any]] = {}
    wanted = set()
    for row in rows:
        key = (str(row["point_id"]), str(row["mass_convention"]), 80, "analytic")
        wanted.add(key)
    for line_no, line in enumerate(output_path.read_text().splitlines(), 1):
        try:
            record = json.loads(line)
        except json.JSONDecodeError as exc:
            raise ValueError(f"malformed frozen T4 JSONL line {line_no}") from exc
        key = (str(record.get("point_id")), str(record.get("mass_convention")),
               record.get("precision_digits_constructed"), record.get("mode"))
        if key not in wanted:
            continue
        point_key = key[0] + "|" + key[1]
        if point_key in selected:
            raise ValueError(f"duplicate frozen T4 analytic/80 record: {point_key}")
        selected[point_key] = record
    if len(selected) != len(rows):
        raise ValueError(f"missing frozen T4 analytic/80 joins: {len(selected)}/{len(rows)}")
    identity = {"checkpoint_scan_identity": identity["fingerprint"],
                "checkpoint_identity": identity,
                "checkpoint_sha256": file_sha(checkpoint_path),
                "output_sha256": file_sha(output_path),
                "output_records": expected,
                "selected_mode": "analytic", "selected_precision_digits": 80}
    return identity, selected


def _run_sm_probe(rows:list[dict[str,Any]])->dict[str,str]:
    requested_vs={q(row["v_GeV"]) for row in rows}
    if requested_vs!={sp.Integer(246)}:
        raise ValueError(f"isolated SM.get_v2 probe only covers v=246 GeV, source grid has {requested_vs}")
    probe=CAMPAIGN/"verification/sm_v2_rounding_probe"
    completed=subprocess.run([str(probe)],check=True,capture_output=True,text=True,timeout=10)
    observed={}
    for line in completed.stdout.splitlines():
        if "=" in line:
            key,value=line.split("=",1)
            observed[key]=value
    if observed.get("get_v2")!=SM_GET_V2_TEXT or observed.get("requested_v_squared")!="60516":
        raise ValueError(f"SM rounding probe output changed: {completed.stdout!r}")
    return {"stdout":completed.stdout.strip(),"GF":observed.get("GF",""),
            "get_v2":observed["get_v2"],"requested_v_squared":observed["requested_v_squared"]}


def exact_target_model(row: dict[str, Any]) -> tuple[RealPotential, dict[str, str]]:
    """Build exact rational T1 family coefficients from decimal input lexemes."""
    t, L = q(row["tan_beta"]), q(row["lambda1_ref"])
    mh, mp, ma, mhp, v = (q(row[key]) for key in ("mh_GeV", "mphi_GeV", "mA_GeV", "mHp_GeV", "v_GeV"))
    r = 1 / t
    a = mh**2 / v**2
    d = (mp**2 - mh**2/2) / v**2
    k = L-a
    b = 2*mp**2/v**2
    if row["trajectory"] == "fixed_lambda6_1e-10":
        l6 = q(row["lambda6"])
        X = l6*t
    elif row["trajectory"] == "Xstar":
        # Xstar is an exact trajectory definition.  The displayed/evaluated_X
        # is a rounded diagnostic and must not seed the exact target model.
        X = (b-L-(L-a)*r**2)/(1+r**2)
        l6 = X/t
    else:
        X = q(row["evaluated_X"])
        l6 = X/t
    l7 = -2*d*r+3*l6*r**2+2*k*r**3
    l2 = a+4*d*r**2-4*l6*r**3-3*k*r**4
    common = 2*l6*r+k*r**2
    l3 = 2*mhp**2/v**2-common
    l4 = (mp**2+ma**2-2*mhp**2)/v**2-common
    l5 = (mp**2-ma**2)/v**2-common
    z = (a,
         L-4*l6*r+2*(b-L)*r**2-4*l6*r**3-3*k*r**4,
         2*mhp**2/v**2,
         (mp**2+ma**2-2*mhp**2)/v**2,
         (mp**2-ma**2)/v**2,
         sp.S.Zero,
         -l6*(1+r**2)+(b-L)*r-k*r**3)
    lambdas = (L,l2,l3,l4,l5,l6,l7)
    values = {**{f"lambda{i}":str(value) for i,value in enumerate(lambdas,1)},
              **{f"Z{i}":str(value) for i,value in enumerate(z,1)},
              "Y1_GeV2":str(-mh**2/2), "Y2_GeV2":"0", "Y3_GeV2":"0",
              "v2_GeV2":str(v**2), "X_lambda6_tanbeta":str(X)}
    if row["trajectory"] == "Xstar":
        values["Xstar_exact"] = str(X)
        values["evaluated_X_decimal80"] = str(row["evaluated_X"])
    return RealPotential.from_mapping(values), values


def _constraint(status: str, method: str, observable: str, tolerance: str,
                precision: str, evidence: Any, failure_reason: str = "") -> dict[str, Any]:
    if status not in ("PASS", "FAIL", "INCONCLUSIVE", "NOT_TESTED"):
        raise ValueError(f"invalid constraint status {status!r}")
    if status in ("FAIL", "INCONCLUSIVE", "NOT_TESTED") and not failure_reason:
        if isinstance(evidence, str):
            failure_reason = evidence
        elif isinstance(evidence, dict):
            failure_reason = str(evidence.get("reason") or evidence.get("evidence") or
                                 evidence.get("exact_failure_reason") or
                                 json.dumps(evidence, sort_keys=True, default=str))
        else:
            failure_reason = str(evidence)
    return {"status":status,"method":method,"observable":observable,
            "tolerance":tolerance,"precision":precision,"evidence":evidence,
            "failure_reason":failure_reason}


def _copy_result(result: dict[str, Any]) -> dict[str, Any]:
    return _constraint(result["status"],result.get("method","exact independent certification"),
        result.get("observable","see evidence"),str(result.get("tolerance","exact comparisons")),
        str(result.get("precision","exact rational/algebraic arithmetic")),result.get("evidence"),
        result.get("failure_reason", ""))


def quartic_screen(lambdas: dict[str, Any], zvals: dict[str, Any]) -> dict[str, Any]:
    values={**lambdas,**zvals}
    maximum=max(abs(q(v)) for v in values.values())
    pi_lo,pi_hi=_pi_interval(100)
    lower,upper=4*pi_lo,4*pi_hi
    status="PASS" if maximum<lower else "FAIL" if maximum>=upper else "INCONCLUSIVE"
    return _constraint(status,"exact rational magnitude screen with rigorous Machin pi bounds",
        "all lambda1..7 and Z1..7","strict max(abs(q)) < 4*pi",
        "exact rational inputs and 100-term alternating arctangent bounds",
        {"maximum_absolute_quartic":str(maximum),"4pi_interval":[str(lower),str(upper)],
         "all_lambdas":lambdas,"all_Z":zvals})


def _exact_target_constraints(p: RealPotential, coefficients: dict[str, str]) -> dict[str, Any]:
    y1,y2,y3=p.Y;z1,z2,z3,z4,z5,z6,z7=p.Z
    v2=p.v2
    assert v2 is not None
    stationarity_res=(y1+z1*v2/2,y2,y3+z6*v2/2)
    stationary=all(item==0 for item in stationarity_res)
    hessian=sp.Matrix([[z1*v2,z6*v2],
                       [z6*v2,y2+(z3+z4+z5)*v2/2]])
    det=sp.factor(hessian.det())
    h11=z1*v2
    a2=y2+(z3+z4-z5)*v2/2
    hp2=y2+z3*v2/2
    local=(h11>0 and det>0 and a2>0 and hp2>0)
    mh2=q(coefficients["lambda1"])*0 + q(_row_value(coefficients,"mh_GeV"))**2 if "mh_GeV" in coefficients else h11
    closure=(h11==mh2 and det==mh2*q(coefficients.get("mphi_GeV",0))**2
             and a2==q(coefficients.get("mA_GeV",0))**2
             and hp2==q(coefficients.get("mHp_GeV",0))**2) if all(k in coefficients for k in ("mh_GeV","mphi_GeV","mA_GeV","mHp_GeV")) else None
    # Positive CP-even principal minors, pseudoscalar and charged eigenvalues
    # are exact field-Hessian tests. The target closure residuals retain the
    # mass identities independently.
    return {
        "stationarity":_constraint("PASS" if stationary else "FAIL",
            "exact Higgs-basis potential gradients", "Y1+Z1*v2/2, Y2, Y3+Z6*v2/2",
            "all three residuals exactly zero", "exact rational arithmetic",
            {"residuals_GeV2":[str(v) for v in stationarity_res]}),
        "alignment":_constraint("PASS" if z6==0 else "FAIL","exact target Z6 coefficient",
            "Z6 alignment","Z6=0 exactly","exact rational arithmetic",{"Z6":str(z6)}),
        "local_minimum_tree":_constraint("PASS" if local else "FAIL",
            "exact Hessian principal minors and physical odd/charged masses",
            "physical scalar Hessian eigenvalues","strictly positive; Goldstones excluded",
            "exact rational arithmetic",{"even_hessian_11_GeV2":str(h11),
                "even_hessian_determinant_GeV4":str(det),"A_squared_GeV2":str(a2),
                "Hplus_squared_GeV2":str(hp2)}),
        "spectrum_reconstruction":_constraint("PASS" if closure else ("INCONCLUSIVE" if closure is None else "FAIL"),
            "exact target spectrum identities","h, phi, A, H+ masses",
            "exact equality to each declared input mass squared","exact rational arithmetic",
            {"closure_exact":closure,"inputs":{k:coefficients.get(k) for k in ("mh_GeV","mphi_GeV","mA_GeV","mHp_GeV")}}),
    }


def _row_value(coefficients: dict[str,str], name:str):
    return coefficients[name]


def _certify_target(p: RealPotential, coefficients: dict[str,str], row:dict[str,Any]) -> dict[str,Any]:
    lam={f"lambda{i}":coefficients[f"lambda{i}"] for i in range(1,8)}
    z={f"Z{i}":coefficients[f"Z{i}"] for i in range(1,8)}
    constraints=_exact_target_constraints(p,coefficients)
    constraints["quartic_screen_4pi"]=quartic_screen(lam,z)
    constraints["full_tree_unitarity"]=_copy_result(certify_tree_unitarity(lam))
    constraints["BFB_general"]=_copy_result(certify_bfb(p))
    global_result=certify_global_minimum(p,(p.v2/2,0,0,0))
    enum=global_result.get("evidence",{}).get("stationary_enumeration")
    if global_result["status"]=="PASS" and (not isinstance(enum,dict) or enum.get("status")=="NOT_TESTED"):
        enum=enumerate_stationary_candidates(p)
        if not enum["complete"]:
            global_result={**global_result,"status":"INCONCLUSIVE",
                "evidence":{**global_result.get("evidence",{}),"stationary_enumeration":enum,
                    "reason":"The global decomposition passes, but requested competing-extrema enumeration is incomplete."}}
    constraints["global_EW_minimum"]=_copy_result(global_result)
    if isinstance(enum,dict) and enum.get("status")!="NOT_TESTED":
        constraints["stationary_candidate_enumeration"]=_constraint(
            enum["status"],enum["method"],enum["observable"],enum["tolerance"],enum["precision"],enum)
    else:
        constraints["stationary_candidate_enumeration"]=_constraint("NOT_TESTED",
            "not needed to establish exact negative global-minimum witness",
            "complete competing stationary point enumeration","required only for global PASS",
            "not run",enum or "global proof did not expose an enumeration")
    return {"constraints":constraints,"target_coefficients_exact":coefficients}


def _certify_rounded(calc:dict[str,Any],row:dict[str,Any]) -> dict[str,Any]:
    doubles=calc.get("double_backend_inputs")
    if not isinstance(doubles,dict) or not all(f"lambda{i}" in doubles for i in range(1,8)) or "m12_sq_GeV2" not in doubles:
        reason="frozen analytic/80 calculator record lacks complete rounded lambda and m12 inputs"
        nottested=_constraint("NOT_TESTED","calculator input extraction","rounded backend real potential",
            "all seven lambdas, m12^2, tan(beta), and requested v required","not available",reason)
        return {"constraints":{key:nottested for key in ("stationarity","alignment","local_minimum_tree",
            "quartic_screen_4pi","full_tree_unitarity","BFB_general","global_EW_minimum",
            "stationary_candidate_enumeration")},"rounded_backend_inputs":doubles}
    inputs={**{f"lambda{i}":doubles[f"lambda{i}"] for i in range(1,8)},
            "m12_sq_GeV2":doubles["m12_sq_GeV2"],
            "tan_beta":calc["double_backend_inputs"].get("tan_beta",calc["tan_beta_backend"]),
            "v_GeV":calc.get("requested_v_GeV",row["v_GeV"])}
    # The explicit/analytic bridge sets GF from requested v; SM.get_v2 then
    # returns this independently measured binary64 value.
    rounded=rounded_input_potential(inputs,v2_sm=SM_GET_V2_BINARY64)
    p=RealPotential.from_mapping(rounded["potential"])
    # Use the exact dyadic rationals that define the reconstructed rounded
    # model, rather than reparsing shortened float display strings.
    lam={f"lambda{i}":rounded["inputs_binary64_rational"][f"lambda{i}"] for i in range(1,8)}
    z={f"Z{i}":str(p.Z[i-1]) for i in range(1,8)}
    constraints={}
    constraints["stationarity"]=_constraint("PASS" if rounded["stationary"] else "FAIL",
        "exact reconstruction of rounded backend inputs","rounded-input tadpole identities",
        "Y1+Z1*v2/2=0 and Y3+Z6*v2/2=0 exactly","exact rationals from binary64 inputs",
        rounded["stationarity_residuals"])
    z6=sp.Rational(z["Z6"])
    align_tol=sp.Rational(str(frozen_t4.DECLARED_ACCEPTANCE_TOLERANCES["alignment_abs"]))
    constraints["alignment"]=_constraint("PASS" if abs(z6)<=align_tol else "FAIL",
        "rounded-input Higgs-basis transform","rounded-model alignment",
        f"abs(Z6) <= {align_tol}","exact rational coefficient; declared calculator tolerance",
        {"Z6":str(z6),"absolute_tolerance":str(align_tol)})
    constraints["local_minimum_tree"]=_constraint("PASS" if rounded["physical_positive"] else "FAIL",
        "exact rounded-input physical Hessian positivity","CP-even, CP-odd and charged masses",
        "strict positive principal minors/eigenvalues","exact rationals from binary64 inputs",
        rounded["physical_hessian_and_masses"])
    constraints["quartic_screen_4pi"]=quartic_screen(lam,z)
    constraints["full_tree_unitarity"]=_copy_result(certify_tree_unitarity(lam))
    constraints["BFB_general"]=_copy_result(certify_bfb(p))
    global_result=certify_global_minimum(p,(p.v2/2,0,0,0))
    enum=global_result.get("evidence",{}).get("stationary_enumeration")
    if global_result["status"]=="PASS" and (not isinstance(enum,dict) or enum.get("status")=="NOT_TESTED"):
        enum=enumerate_stationary_candidates(p)
        if not enum["complete"]:
            global_result={**global_result,"status":"INCONCLUSIVE",
                "evidence":{**global_result.get("evidence",{}),"stationary_enumeration":enum,
                    "reason":"The global decomposition passes, but requested competing-extrema enumeration is incomplete."}}
    constraints["global_EW_minimum"]=_copy_result(global_result)
    if isinstance(enum,dict) and enum.get("status")!="NOT_TESTED":
        constraints["stationary_candidate_enumeration"]=_constraint(enum["status"],enum["method"],
            enum["observable"],enum["tolerance"],enum["precision"],enum)
    else:
        constraints["stationary_candidate_enumeration"]=_constraint("NOT_TESTED",
            "no complete enumeration required for exact negative/unboundedness rejection",
            "competing stationary candidate enumeration","required for global PASS","not run",
            enum or "Global certificate did not expose candidate enumeration")
    return {"constraints":constraints,"potential":rounded["potential"],
            "rounded_reconstruction":rounded,"rounded_backend_inputs":inputs,
            "SM_get_v2_binary64":repr(SM_GET_V2_BINARY64)}


def _aggregate(constraints:dict[str,dict[str,Any]]) -> str:
    statuses={item["status"] for item in constraints.values()}
    if "FAIL" in statuses:return "FAIL"
    if statuses & {"INCONCLUSIVE","NOT_TESTED"}:return "INCONCLUSIVE"
    return "PASS"


def certify_point(row:dict[str,Any],calc:dict[str,Any],identity:dict[str,Any]) -> dict[str,Any]:
    target,target_coefficients=exact_target_model(row)
    target_coefficients.update({k:str(row[k]) for k in ("mh_GeV","mphi_GeV","mA_GeV","mHp_GeV","v_GeV")})
    exact=_certify_target(target,target_coefficients,row)
    rounded=_certify_rounded(calc,row)
    calculator_status=calc.get("constructor_status","NOT_TESTED")
    calc_reason=calc.get("exact_failure_reason","")
    calculator=_constraint(calculator_status,"frozen 2HDMC analytic constructor at 80 digits",
        "stable-family versus actual analytic calculator constructor agreement",
        json.dumps(calc.get("declared_acceptance_tolerances",{}),sort_keys=True),
        "frozen calculator record; binary64 backend",
        {"record_payload_sha256":calc["record_payload_sha256"],
         "constructor_checks":calc.get("constructor_checks"),
         "physical_state_ordering":calc.get("physical_state_ordering")},calc_reason)
    if calc.get("widths_computed"):
        phenomenology={
            "width_total_phi_GeV":calc.get("width_total_phi_GeV"),
            "width_channels_GeV":calc.get("width_channels_GeV"),
            "branching_fractions_reported_channels":calc.get("branching_fractions_reported_channels"),
            "physical_lifetime":calc.get("physical_lifetime"),
            "decay_status":calc.get("decay_status"),
        }
    else:
        phenomenology={"status":"NOT_TESTED","failure_reason":
            "The frozen analytic/80 calculator record has no completed decay calculation."}
    phenomenology["displaced_photon_selection"]={"status":"NOT_TESTED","failure_reason":
        "No applicable validated displaced-photon selection is implemented in this tree-level certification runner."}
    phenomenology["experimental_exclusions"]={"status":"NOT_TESTED","failure_reason":
        "No applicable validated experimental exclusion method is included in this campaign."}
    all_constraints={**{f"target.{k}":v for k,v in exact["constraints"].items()},
                     **{f"rounded.{k}":v for k,v in rounded["constraints"].items()},
                     "calculator.analytic_constructor":calculator}
    overall=_aggregate(all_constraints)
    if overall=="PASS":
        point_class="CERTIFIED_TREE_LEVEL"
    elif overall=="FAIL":
        point_class="REJECTED"
    else:
        point_class="CANDIDATE_VALIDITY_INCONCLUSIVE"
    failure_reasons={key:value["failure_reason"] for key,value in all_constraints.items()
                     if value["status"] in ("FAIL","INCONCLUSIVE","NOT_TESTED")}
    record={"schema":SCHEMA,"method_version":METHOD_VERSION,
        "physical_scan_identity":identity["fingerprint"],
        "point_id":row["point_id"],"mass_convention":row["mass_convention"],
        "slice":row["slice"],"trajectory":row["trajectory"],
        "point_input_sha256":frozen_t4.row_fingerprint(row),
        "frozen_t4_record_sha256":calc["record_payload_sha256"],
        "frozen_t4_scan_identity":calc["scan_identity"],
        "pilot_input_lexemes":{k:str(row[k]) for k in ("tan_beta","lambda1_ref","mh_GeV","mphi_GeV","mA_GeV","mHp_GeV","v_GeV","evaluated_X","lambda6")},
        "source_decimal80_lambda_Z":{**{f"lambda{i}":str(row[f"lambda{i}"]) for i in range(1,8)},
                                     **{f"Z{i}":str(row[f"Z{i}"]) for i in range(1,8)}},
        "target_exact":exact,"rounded_backend":rounded,
        "calculator_constructor":calculator,
        "phenomenology":phenomenology,
        "constraints":all_constraints,"overall_status":overall,"point_class":point_class,
        "exact_failure_reasons":failure_reasons,
        "provenance":{"point_input_sha256":frozen_t4.row_fingerprint(row),
            "frozen_t4_analytic80_payload_sha256":calc["record_payload_sha256"],
            "frozen_t4_identity":calc["scan_identity"],
            "SM_v2_probe_value":SM_GET_V2_TEXT}}
    record["point_payload_sha256"]=canonical_sha(record)
    return record


def _source_spec(rows:list[dict[str,Any]],t4_identity:dict[str,Any],sm_probe:dict[str,str]) -> dict[str,Any]:
    sources=(HERE/"run_physical_certification.py",REPO/"dihiggs/app/orchestrator/physical_certification.py",
        REPO/"dihiggs/app/orchestrator/aligned_lambda7.py",REPO/"dihiggs/app/orchestrator/aligned_constraints.py",
        HERE/"run_t4_scan.py",HERE/"run_t4_bridge.py",HERE/"scan_lock.py",HERE/"analytic_bridge",
        CAMPAIGN/"verification/sm_v2_rounding_probe.cpp",CAMPAIGN/"verification/sm_v2_rounding_probe",
        CAMPAIGN/"build_analytic/2hdmc/lib/lib2HDMC.a",REPO/"2hdmc/src/SM.cpp",REPO/"2hdmc/src/SM.h")
    source_hashes={str(path.relative_to(REPO)):file_sha(path) for path in sources}
    grid_paths=(CAMPAIGN/"t2_lightweight/trajectory_points.jsonl",CAMPAIGN/"pilot_scan_tables_mh125_20.csv")
    grid_hashes={str(path.relative_to(CAMPAIGN)):file_sha(path) for path in grid_paths}
    methods={"exact_target":"exact rational stable T1 coefficients from original decimal input lexemes",
        "rounded_model":"Fraction.from_float binary64 constructor inputs plus isolated SM.get_v2 measurement",
        "bfb":"exact rational witness or S-lemma principal-minor certificate",
        "global":"exact full-domain aligned decomposition or complete stationary enumeration and ordering",
        "unitarity":"exact GI principal minors with rigorous Machin rational bounds for pi",
        "quartic_screen":"strict max(abs(lambda,Z)) < 4*pi with rigorous pi bounds",
        "mass_slices":{"mh125.13_T1":"T1 reproduction","mh125.20_active":"active convention; kept as a separate point key"}}
    spec={"method_version":METHOD_VERSION,"schema":SCHEMA,"source_sha256":source_hashes,
        "grid_sha256":grid_hashes,"frozen_t4":t4_identity,
        "SM_get_v2":{"requested_v_GeV":"246","measured_binary64":SM_GET_V2_TEXT,
                      "probe_stdout":sm_probe["stdout"],"probe_GF_binary64":sm_probe["GF"],
                      "probe_source_sha256":source_hashes[str((CAMPAIGN/"verification/sm_v2_rounding_probe.cpp").relative_to(REPO))],
                      "probe_executable_sha256":source_hashes[str((CAMPAIGN/"verification/sm_v2_rounding_probe").relative_to(REPO))],
                      "SM_library_sha256":source_hashes[str((CAMPAIGN/"build_analytic/2hdmc/lib/lib2HDMC.a").relative_to(REPO))]},
        "pilot_rows":len(rows),"point_input_sha256":[frozen_t4.row_fingerprint(row) for row in rows],
        "precision_digits":PRECISION_DIGITS,"tolerances":{"quartic_screen":"strict < 4*pi",
            "unitarity":"strict |Lambda| < 8*pi; a0=Lambda/(16*pi)",
            "alignment_abs":frozen_t4.DECLARED_ACCEPTANCE_TOLERANCES["alignment_abs"],
            "algebraic_signs":"exact rational/root-isolation signs; unresolved signs are INCONCLUSIVE"},
        "methods":methods,"resource_policy":{"workers":1,"address_space_bytes":ADDRESS_SPACE_BYTES,
            "cpu_seconds_per_batch":CPU_SECONDS,"max_points_per_batch":MAX_POINTS_PER_BATCH}}
    return {**spec,"fingerprint":canonical_sha(spec)}


def _write_checkpoint(identity:dict[str,Any],done:dict[str,dict[str,Any]],rows:list[dict[str,Any]],
                      status:str,start_index:int,limit:int,pending_point:dict[str,Any]|None=None)->None:
    payload={"schema":SCHEMA,"status":status,"scan_identity":identity["fingerprint"],"identity":identity,
        "completed_points":len(done),"expected_points":len(rows),
        "completed_point_keys":sorted(done),"last_completed_point_key":next(reversed(done),None),
        "start_index":start_index,"limit_points":limit,"worker_count":1,
        "results_jsonl":RESULTS.name}
    if pending_point is not None:
        payload["pending_point"]=pending_point
    OUT.mkdir(parents=True,exist_ok=True)
    tmp=CHECKPOINT.with_suffix(".tmp")
    with tmp.open("w",encoding="utf-8") as stream:
        stream.write(json.dumps(payload,indent=2,sort_keys=True)+"\n")
        stream.flush();os.fsync(stream.fileno())
    tmp.replace(CHECKPOINT)


_TARGET_CONSTRAINTS={"stationarity","alignment","local_minimum_tree","spectrum_reconstruction",
    "quartic_screen_4pi","full_tree_unitarity","BFB_general","global_EW_minimum",
    "stationary_candidate_enumeration"}
_ROUNDED_CONSTRAINTS={"stationarity","alignment","local_minimum_tree","quartic_screen_4pi",
    "full_tree_unitarity","BFB_general","global_EW_minimum","stationary_candidate_enumeration"}
_CONSTRAINT_METADATA={"status","method","observable","tolerance","precision","evidence","failure_reason"}


def _validate_complete_record(record:dict[str,Any],row:dict[str,Any],calc:dict[str,Any],line_no:int)->None:
    prefix=f"physical record line {line_no}: "
    if not isinstance(record.get("target_exact"),dict) or not isinstance(record.get("rounded_backend"),dict):
        raise ValueError(prefix+"missing exact-target or rounded-model reconstruction")
    expected_source={**{f"lambda{i}":str(row[f"lambda{i}"]) for i in range(1,8)},
                     **{f"Z{i}":str(row[f"Z{i}"]) for i in range(1,8)}}
    if record.get("source_decimal80_lambda_Z")!=expected_source:
        raise ValueError(prefix+"source decimal80 lambda/Z values are incomplete or mismatched")
    for field in ("target_exact","rounded_backend"):
        if not isinstance(record.get(field),dict):
            raise ValueError(prefix+f"missing {field}")
    target_values=record["target_exact"].get("target_coefficients_exact")
    expected_target,expected_target_values=exact_target_model(row)
    del expected_target
    expected_target_values.update({k:str(row[k]) for k in ("mh_GeV","mphi_GeV","mA_GeV","mHp_GeV","v_GeV")})
    if target_values!=expected_target_values:
        raise ValueError(prefix+"exact target Y/Z/lambda coefficients or trajectory values are incomplete/mismatched")
    if record["target_exact"].get("constraints") is None:
        raise ValueError(prefix+"exact-target constraint evidence is missing")

    doubles=calc.get("double_backend_inputs")
    if not isinstance(doubles,dict):
        raise ValueError(prefix+"frozen calculator backend inputs are missing")
    expected_backend={**{f"lambda{i}":doubles[f"lambda{i}"] for i in range(1,8)},
        "m12_sq_GeV2":doubles["m12_sq_GeV2"],
        "tan_beta":doubles.get("tan_beta",calc["tan_beta_backend"]),
        "v_GeV":calc.get("requested_v_GeV",row["v_GeV"])}
    rounded=record["rounded_backend"]
    rounded_inputs=rounded.get("rounded_backend_inputs")
    if not isinstance(rounded_inputs,dict) or set(rounded_inputs)!=set(expected_backend):
        raise ValueError(prefix+"rounded backend inputs are incomplete")
    if any(q(rounded_inputs[k])!=q(expected_backend[k]) for k in expected_backend):
        raise ValueError(prefix+"rounded backend inputs mismatch the frozen T4 record")
    expected_rounded=rounded_input_potential(expected_backend,v2_sm=SM_GET_V2_BINARY64)
    if rounded.get("potential")!=expected_rounded["potential"] or rounded.get("rounded_reconstruction")!=expected_rounded:
        raise ValueError(prefix+"rounded exact-rational Y/Z potential reconstruction is incomplete or mismatched")
    target_lambdas={f"lambda{i}":expected_target_values[f"lambda{i}"] for i in range(1,8)}
    target_z={f"Z{i}":expected_target_values[f"Z{i}"] for i in range(1,8)}
    rounded_lambdas={f"lambda{i}":expected_rounded["inputs_binary64_rational"][f"lambda{i}"] for i in range(1,8)}
    rounded_potential=RealPotential.from_mapping(expected_rounded["potential"])
    rounded_z={f"Z{i}":str(v) for i,v in enumerate(rounded_potential.Z,1)}
    expected_screens={"target":quartic_screen(target_lambdas,target_z),
                      "rounded":quartic_screen(rounded_lambdas,rounded_z)}
    expected_unitarity={"target":_copy_result(certify_tree_unitarity(target_lambdas)),
                        "rounded":_copy_result(certify_tree_unitarity(rounded_lambdas))}

    expected_constraint_keys={*(f"target.{name}" for name in _TARGET_CONSTRAINTS),
        *(f"rounded.{name}" for name in _ROUNDED_CONSTRAINTS),"calculator.analytic_constructor"}
    constraints=record.get("constraints")
    if not isinstance(constraints,dict) or set(constraints)!=expected_constraint_keys:
        raise ValueError(prefix+"required physics constraint set is incomplete or contains unknown entries")
    for name,item in constraints.items():
        if not isinstance(item,dict) or not _CONSTRAINT_METADATA.issubset(item):
            raise ValueError(prefix+f"constraint {name} lacks required status/method/observable/tolerance/precision/evidence metadata")
        if item["status"] not in ("PASS","FAIL","INCONCLUSIVE","NOT_TESTED"):
            raise ValueError(prefix+f"invalid constraint status {name}")
        if not all(isinstance(item[k],str) and item[k] for k in ("method","observable","tolerance","precision")):
            raise ValueError(prefix+f"constraint {name} has empty method/observable/tolerance/precision")
        if item["evidence"] is None:
            raise ValueError(prefix+f"constraint {name} has no evidence")
        if item["status"]!="PASS" and not item["failure_reason"]:
            raise ValueError(prefix+f"constraint {name} has no exact failure/uncertainty reason")
    for name,item in record["target_exact"]["constraints"].items():
        if name not in _TARGET_CONSTRAINTS or constraints.get("target."+name)!=item:
            raise ValueError(prefix+"target nested constraints do not match complete campaign constraint evidence")
    for name,item in rounded.get("constraints",{}).items():
        if name not in _ROUNDED_CONSTRAINTS or constraints.get("rounded."+name)!=item:
            raise ValueError(prefix+"rounded nested constraints do not match complete campaign constraint evidence")
    if set(record["target_exact"]["constraints"])!=_TARGET_CONSTRAINTS or set(rounded.get("constraints",{}))!=_ROUNDED_CONSTRAINTS:
        raise ValueError(prefix+"nested target/rounded constraint set is incomplete")
    for side in ("target","rounded"):
        if record["constraints"][f"{side}.quartic_screen_4pi"]!=expected_screens[side]:
            raise ValueError(prefix+f"{side} quartic-screen evidence is inconsistent with its exact lambda/Z coefficients")
        if record["constraints"][f"{side}.full_tree_unitarity"]!=expected_unitarity[side]:
            raise ValueError(prefix+f"{side} unitarity evidence is inconsistent with its exact lambda coefficients")
    calculator=record.get("calculator_constructor")
    if not isinstance(calculator,dict) or constraints["calculator.analytic_constructor"].get("status")!=calc.get("constructor_status"):
        raise ValueError(prefix+"calculator constructor status evidence is missing or inconsistent")
    if record.get("exact_failure_reasons")!={k:v["failure_reason"] for k,v in constraints.items()
            if v["status"] in ("FAIL","INCONCLUSIVE","NOT_TESTED")}:
        raise ValueError(prefix+"exact failure/uncertainty reason index is inconsistent")
    overall=_aggregate(constraints)
    point_class=("CERTIFIED_TREE_LEVEL" if overall=="PASS" else "REJECTED" if overall=="FAIL"
                 else "CANDIDATE_VALIDITY_INCONCLUSIVE")
    if record.get("overall_status")!=overall or record.get("point_class")!=point_class:
        raise ValueError(prefix+"overall status/classification does not follow required constraint outcomes")


def _load_existing(rows:list[dict[str,Any]],identity:dict[str,Any],
                   t4_records:dict[str,dict[str,Any]])->dict[str,dict[str,Any]]:
    if not RESULTS.exists() and not CHECKPOINT.exists():return {}
    if RESULTS.exists() != CHECKPOINT.exists():
        raise ValueError("physical output/checkpoint pair is incomplete; preserve and inspect before replay")
    checkpoint=json.loads(CHECKPOINT.read_text())
    if checkpoint.get("schema")!=SCHEMA:
        raise ValueError("physical checkpoint schema mismatch")
    if checkpoint.get("scan_identity")!=identity["fingerprint"] or checkpoint.get("identity")!=identity:
        raise ValueError("physical checkpoint identity mismatch: source, grid, calculator, method, or tolerances changed")
    expected={row["point_id"]+"|"+row["mass_convention"]:row for row in rows}
    done={}
    for line_no,line in enumerate(RESULTS.read_text().splitlines(),1):
        try:record=json.loads(line)
        except json.JSONDecodeError as exc:raise ValueError(f"malformed physical record line {line_no}") from exc
        if record.get("schema")!=SCHEMA:raise ValueError(f"physical record schema mismatch at line {line_no}")
        key=str(record.get("point_id"))+"|"+str(record.get("mass_convention"))
        if key not in expected:raise ValueError(f"unknown point/convention key at line {line_no}: {key}")
        if key in done:raise ValueError(f"duplicate physical point record at line {line_no}: {key}")
        if record.get("physical_scan_identity")!=identity["fingerprint"]:
            raise ValueError(f"physical record identity mismatch at line {line_no}")
        row=expected[key]
        if record.get("point_input_sha256")!=frozen_t4.row_fingerprint(row):
            raise ValueError(f"physical record point-input hash mismatch at line {line_no}")
        calc=t4_records[key]
        if record.get("frozen_t4_record_sha256")!=calc["record_payload_sha256"] or record.get("frozen_t4_scan_identity")!=calc["scan_identity"]:
            raise ValueError(f"physical record frozen calculator join mismatch at line {line_no}")
        payload={k:v for k,v in record.items() if k!="point_payload_sha256"}
        if record.get("point_payload_sha256")!=canonical_sha(payload):
            raise ValueError(f"physical record payload checksum mismatch at line {line_no}")
        _validate_complete_record(record,row,calc,line_no)
        done[key]=record
    expected_done=checkpoint.get("completed_point_keys")
    if sorted(done)!=expected_done or checkpoint.get("completed_points")!=len(done):
        raise ValueError("checkpoint completion index does not match validated physical JSONL records")
    pending=checkpoint.get("pending_point")
    if pending is not None:
        pending_key=pending.get("point_key")
        if checkpoint.get("status")!="RESOURCE_LIMITED" or pending_key not in expected or pending_key in done:
            raise ValueError("checkpoint pending-point evidence is inconsistent with completion state")
    if checkpoint.get("status")=="CAMPAIGN_COMPLETE" and (len(done)!=len(rows) or pending is not None):
        raise ValueError("campaign checkpoint claims completion with missing or pending points")
    if checkpoint.get("status")=="BATCH_COMPLETE":
        start=checkpoint.get("start_index");limit=checkpoint.get("limit_points")
        if not isinstance(start,int) or not isinstance(limit,int) or start<0 or limit<1:
            raise ValueError("batch-complete checkpoint has an invalid selected range")
        selected=rows[start:start+limit]
        if pending is not None or any(r["point_id"]+"|"+r["mass_convention"] not in done for r in selected):
            raise ValueError("batch checkpoint claims completion with unfinished selected points")
    return done


class BatchCpuLimit(Exception):
    pass


def _enforce_resources()->None:
    old_as=resource.getrlimit(resource.RLIMIT_AS)
    hard_as=ADDRESS_SPACE_BYTES if old_as[1]==resource.RLIM_INFINITY else min(old_as[1],ADDRESS_SPACE_BYTES)
    soft_as=hard_as if old_as[0]==resource.RLIM_INFINITY else min(old_as[0],hard_as)
    resource.setrlimit(resource.RLIMIT_AS,(soft_as,hard_as))
    soft_cpu,hard_cpu=resource.getrlimit(resource.RLIMIT_CPU)
    hard=CPU_SECONDS if hard_cpu==resource.RLIM_INFINITY else min(hard_cpu,CPU_SECONDS)
    soft=min(hard,CPU_SECONDS-10)
    resource.setrlimit(resource.RLIMIT_CPU,(soft,hard))
    def cpu_limit(_signum,_frame):raise BatchCpuLimit("batch reached its enforced CPU budget")
    signal.signal(signal.SIGXCPU,cpu_limit)


def run(args:argparse.Namespace)->dict[str,Any]:
    rows=_source_rows()
    t4_identity,t4_records=_frozen_t4(rows)
    sm_probe=_run_sm_probe(rows)
    identity=_source_spec(rows,t4_identity,sm_probe)
    with campaign_scan_lock(LOCK,identity["fingerprint"]):
        OUT.mkdir(parents=True,exist_ok=True)
        done=_load_existing(rows,identity,t4_records)
        selected=rows[args.start_index:args.start_index+args.limit_points]
        if CHECKPOINT.exists():
            previous_checkpoint=json.loads(CHECKPOINT.read_text())
            pending=previous_checkpoint.get("pending_point")
            if pending is not None and pending.get("point_key") not in {
                    row["point_id"]+"|"+row["mass_convention"] for row in selected}:
                raise ValueError("resume range must include the previously interrupted pending point")
        profile=[]
        pending_point=None
        pending_profile=None
        _enforce_resources()
        for row in selected:
            key=row["point_id"]+"|"+row["mass_convention"]
            if key in done:continue
            start=time.perf_counter();before=resource.getrusage(resource.RUSAGE_SELF)
            try:
                record=certify_point(row,t4_records[key],identity)
            except BatchCpuLimit as exc:
                after=resource.getrusage(resource.RUSAGE_SELF)
                elapsed=time.perf_counter()-start
                pending_point={"point_key":key,"reason":str(exc),
                    "point_input_sha256":frozen_t4.row_fingerprint(row),
                    "frozen_t4_record_sha256":t4_records[key]["record_payload_sha256"],
                    "source_decimal80_lambda_Z":{**{f"lambda{i}":str(row.get(f"lambda{i}","")) for i in range(1,8)},
                        **{f"Z{i}":str(row.get(f"Z{i}","")) for i in range(1,8)}},
                    "calculator_decimal_inputs":t4_records[key].get("inputs_decimal"),
                    "calculator_backend_values_decimal":t4_records[key].get("values_decimal")}
                pending_profile={"point_id":row["point_id"],"mass_convention":row["mass_convention"],
                    "status":"INTERRUPTED_RETRY_REQUIRED","wall_seconds":elapsed,
                    "user_cpu_seconds":after.ru_utime-before.ru_utime,
                    "system_cpu_seconds":after.ru_stime-before.ru_stime,
                    "process_peak_rss_bytes":after.ru_maxrss*1024}
                break
            after=resource.getrusage(resource.RUSAGE_SELF)
            elapsed=time.perf_counter()-start
            record["execution_profile"]={"worker_count":1,"wall_seconds":elapsed,
                "user_cpu_seconds":after.ru_utime-before.ru_utime,
                "system_cpu_seconds":after.ru_stime-before.ru_stime,
                "process_peak_rss_bytes":after.ru_maxrss*1024,
                "address_space_limit_bytes":ADDRESS_SPACE_BYTES,"cpu_batch_limit_seconds":CPU_SECONDS}
            record.pop("point_payload_sha256",None)
            record["point_payload_sha256"]=canonical_sha(record)
            with RESULTS.open("a",encoding="utf-8") as stream:
                stream.write(json.dumps(record,sort_keys=True,separators=(",",":"))+"\n")
                stream.flush();os.fsync(stream.fileno())
            done[key]=record
            _write_checkpoint(identity,done,rows,"RUNNING",args.start_index,args.limit_points)
            profile.append({"point_id":row["point_id"],"mass_convention":row["mass_convention"],
                "overall_status":record["overall_status"],"point_class":record["point_class"],
                "wall_seconds":elapsed,"user_cpu_seconds":after.ru_utime-before.ru_utime,
                "process_peak_rss_bytes":after.ru_maxrss*1024})
            print(f"physical point {len(done)}/{len(rows)} {key}: {record['overall_status']} ({elapsed:.3f}s)",flush=True)
        batch_complete=all(row["point_id"]+"|"+row["mass_convention"] in done for row in selected)
        campaign_complete=len(done)==len(rows)
        status="RESOURCE_LIMITED" if pending_point is not None else "CAMPAIGN_COMPLETE" if campaign_complete else "BATCH_COMPLETE" if batch_complete else "RESOURCE_LIMITED"
        if pending_point is not None:
            # Keep the resumable output/checkpoint pair intact without turning
            # the interrupted point into a result row or completed key.
            RESULTS.touch(exist_ok=True)
        _write_checkpoint(identity,done,rows,status,args.start_index,args.limit_points,pending_point)
        batch_result={"status":status,"physical_scan_identity":identity,"completed_this_run":len(profile),
            "completed_total":len(done),"expected_total":len(rows),"worker_count":1,"point_profiles":profile,
            "pending_point":pending_point,"pending_point_profile":pending_profile,
            "output_jsonl":str(RESULTS.relative_to(CAMPAIGN)),"checkpoint_json":str(CHECKPOINT.relative_to(CAMPAIGN))}
        summary_path=OUT/"physical_certification_batch_profile.json"
        summary_path.write_text(json.dumps(batch_result,indent=2,sort_keys=True)+"\n")
        return batch_result


def main()->None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-index",type=int,default=0)
    parser.add_argument("--limit-points",type=int,default=50)
    args=parser.parse_args()
    if args.limit_points<1 or args.limit_points>MAX_POINTS_PER_BATCH:
        raise SystemExit(f"--limit-points must be in 1..{MAX_POINTS_PER_BATCH}")
    if args.start_index<0 or args.start_index>=390:
        raise SystemExit("--start-index must be in 0..389")
    if args.start_index+args.limit_points>390:
        args.limit_points=390-args.start_index
    result=run(args)
    print(json.dumps({k:v for k,v in result.items() if k!="physical_scan_identity"},sort_keys=True))


if __name__=="__main__":
    main()
