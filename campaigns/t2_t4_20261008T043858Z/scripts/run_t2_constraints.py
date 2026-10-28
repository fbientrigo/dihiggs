#!/usr/bin/env python3
"""Run serial exact-witness/sufficient-certificate screens on the saved grid."""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path
import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from dihiggs.app.orchestrator.aligned_lambda7 import build_aligned_lambda7_family, xstar_lambda6_tanbeta
from dihiggs.app.orchestrator.aligned_constraints import assess_family_point, unitarity

CAMPAIGN=Path(__file__).resolve().parents[1]
LIGHT=CAMPAIGN/"t2_lightweight"


def build_point(row):
    args={key:row[key] for key in ("tan_beta","lambda1_ref","mh_GeV","mphi_GeV","mA_GeV","mHp_GeV","v_GeV")}
    if row["trajectory"]=="fixed_lambda6_1e-10":
        return build_aligned_lambda7_family(**args,lambda6=row["lambda6"],precision=80).as_dict()
    return build_aligned_lambda7_family(**args,X_lambda6_tanbeta=row["evaluated_X"],precision=80).as_dict()


def main():
    source=LIGHT/"trajectory_points.jsonl"
    rows=[json.loads(line) for line in source.read_text().splitlines() if line.strip()]
    assessments=[]
    for index,row in enumerate(rows,1):
        point=build_point(row)
        assessment=assess_family_point(point)
        c=assessment["constraints"]
        core=("spectrum_reconstruction","stationarity","alignment","local_minimum_tree",
              "quartic_screen_4pi","BFB_general","global_EW_minimum","full_tree_unitarity")
        failures=[name for name in core if c[name]["status"]=="FAIL"]
        undecided=[name for name in core if c[name]["status"] not in ("PASS","FAIL")]
        candidate="CANDIDATE_VALID" if not failures and not undecided else "SCREENED_ONLY"
        reasons=[]
        if failures:
            for name in failures:
                evidence=c[name].get("evidence",{})
                why=evidence.get("quartic_value") or c[name].get("method") or name
                reasons.append(f"{name}: {why}")
        if undecided:
            reasons.append("INCONCLUSIVE/NOT_TESTED: "+", ".join(undecided))
        record={"point_id":row["point_id"],"slice":row["slice"],"trajectory":row["trajectory"],
                "inputs":point["inputs"],"values":point["values"],"constraints":c,
                "candidate_status":candidate,"benchmark_certification":"NOT_TESTED",
                "exact_failure_reason":"; ".join(reasons) if reasons else "",
                "calculator_route_agreement":"NOT_TESTED for this scan point",
                "phenomenology":"NOT_TESTED for this scan point"}
        assessments.append(record)
        if index%50==0:
            checkpoint=CAMPAIGN/"physical_constraints_checkpoint.jsonl"
            checkpoint.write_text("".join(json.dumps(x,sort_keys=True)+"\n" for x in assessments))
    out_jsonl=CAMPAIGN/"physical_constraints_report.jsonl"
    out_jsonl.write_text("".join(json.dumps(x,sort_keys=True)+"\n" for x in assessments))
    counts=Counter(x["candidate_status"] for x in assessments)
    status_counts={}
    for name in ("quartic_screen_4pi","BFB_general","global_EW_minimum","full_tree_unitarity",
                 "local_minimum_tree","alignment","spectrum_reconstruction","stationarity"):
        status_counts[name]=dict(Counter(x["constraints"][name]["status"] for x in assessments))
    precision_row=next(row for row in rows if row["slice"]=="illustrative_90"
        and row["trajectory"]=="fixed_X_1" and row["tan_beta"]=="1000000"
        and row["lambda1_ref"]=="1.2")
    precision_point=build_point(precision_row)
    precision_lambdas={f"lambda{i}":precision_point["values"][f"lambda{i}"] for i in range(1,8)}
    spectra={digits:unitarity(precision_lambdas,precision=digits) for digits in (50,80,120)}
    convergence={}
    for low,high in ((50,80),(80,120)):
        with mp.workdps(150):
            convergence[f"{low}_to_{high}"]={block:mp.nstr(max(abs(mp.mpf(spectra[high]["eigenvalues"][block][i])-mp.mpf(spectra[low]["eigenvalues"][block][i]))
                for i in range(len(spectra[high]["eigenvalues"][block]))),120) for block in spectra[low]["eigenvalues"]}
    lightweight=json.loads((LIGHT/"t2_lightweight_report.json").read_text())
    t4file=CAMPAIGN/"T4_integration_report.json"
    precision_report={"status":"MEASURED","scope":"one declared high-t cancellation point; complete GI matrices, 50/80/120 digits; double backend routes separately below",
        "point_inputs":precision_point["inputs"],"decimal_lambdas":precision_lambdas,
        "forward_constructor_inverse_arithmetic":lightweight["precision_audit"],
        "unitarity_spectra_by_precision":spectra,"spectrum_convergence_absolute_max_by_block":convergence,
        "eigensolver_error_budget_method":"maximum of eigensystem residual, eigenvector orthogonality error times matrix norm, and measured spectrum change at +20 decimal digits",
        "calculator_routes":json.loads(t4file.read_text())["routes"] if t4file.exists() else "NOT_RUN"}
    columns=["point_id","slice","trajectory","tan_beta","lambda1_ref","mh_GeV","mphi_GeV","mA_GeV","mHp_GeV","v_GeV","evaluated_X","lambda6",*(f"lambda{i}" for i in range(1,8)),*(f"Z{i}" for i in range(1,8)),"quartic_screen_4pi","BFB_general","global_EW_minimum","full_tree_unitarity","local_minimum_tree","alignment","candidate_status","benchmark_certification","exact_failure_reason"]
    # Evaluate the active mh=125.20 convention independently. Keep these rows
    # and their reports separate from the T1 reproduction table at mh=125.13.
    active_rows=[]
    for row in rows:
        if row["slice"]!="illustrative_90":
            continue
        active=dict(row,mh_GeV="125.20")
        if active["trajectory"]=="Xstar":
            active["evaluated_X"]=str(xstar_lambda6_tanbeta(**{k:active[k] for k in
                ("tan_beta","lambda1_ref","mh_GeV","mphi_GeV","v_GeV")}))
        active_rows.append(active)
    active_records=[]
    for row in active_rows:
        point=build_point(row)
        c=assess_family_point(point)["constraints"]
        core=("spectrum_reconstruction","stationarity","alignment","local_minimum_tree",
              "quartic_screen_4pi","BFB_general","global_EW_minimum","full_tree_unitarity")
        failures=[name for name in core if c[name]["status"]=="FAIL"]
        undecided=[name for name in core if c[name]["status"] not in ("PASS","FAIL")]
        reasons=[f"{name}: {c[name].get('method',name)}" for name in failures]
        if undecided: reasons.append("INCONCLUSIVE/NOT_TESTED: "+", ".join(undecided))
        active_records.append({"point_id":row["point_id"],"inputs":point["inputs"],
            "values":point["values"],"constraints":c,
            "candidate_status":"CANDIDATE_VALID" if not failures and not undecided else "SCREENED_ONLY",
            "benchmark_certification":"NOT_TESTED","exact_failure_reason":"; ".join(reasons)})
    active_jsonl=CAMPAIGN/"active_mh125_20_constraints.jsonl"
    active_jsonl.write_text("".join(json.dumps(x,sort_keys=True)+"\n" for x in active_records))
    active_csv=CAMPAIGN/"pilot_scan_tables_mh125_20.csv"
    with active_csv.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=columns);writer.writeheader()
        for row,record in zip(active_rows,active_records):
            output={k:row.get(k,"") for k in columns}
            for name in ("lambda1","lambda2","lambda3","lambda4","lambda5","lambda6","lambda7",
                         "Z1","Z2","Z3","Z4","Z5","Z6","Z7"):
                output[name]=record["values"][name]
            output["evaluated_X"]=record["values"]["X_lambda6_tanbeta"]
            output["lambda6"]=record["values"]["lambda6"]
            for name in ("quartic_screen_4pi","BFB_general","global_EW_minimum","full_tree_unitarity","local_minimum_tree","alignment"):
                output[name]=record["constraints"][name]["status"]
            output["candidate_status"]=record["candidate_status"]
            output["benchmark_certification"]=record["benchmark_certification"]
            output["exact_failure_reason"]=record["exact_failure_reason"]
            writer.writerow(output)
    precision_report["active_convention_mh125_20"]={"status":"SEPARATE_EVALUATION",
        "point_count":len(active_records),"pilot_csv":active_csv.name,"constraints_jsonl":active_jsonl.name,
        "point_classes":dict(Counter(x["candidate_status"] for x in active_records)),
        "mixed_with_T1_mh125_13":False}
    (CAMPAIGN/"numerical_precision_audit.json").write_text(json.dumps(precision_report,indent=2,sort_keys=True)+"\n")
    summary={"status":"PARTIAL","point_count":len(assessments),"point_classes":dict(counts),
        "constraint_status_counts":status_counts,"precision_digits":80,
        "certificate_scope":"exact decimal-rational bilinear witness or sufficient positive certificate for the generated aligned Z6=0 family; no optimizer PASS",
        "global_minimum_scope":"exact input-rational aligned Higgs-basis family certificate; serialized and backend closure residuals are separately recorded; no arbitrary-potential claim",
        "benchmark_validations":0,"benchmark_status":"no benchmark physically validated; actual calculator agreement and signed decay-width checks are incomplete",
        "point_report_jsonl":str(out_jsonl.relative_to(CAMPAIGN)),
        "unitarity_precision_study":"numerical_precision_audit.json",
        "active_mh125_20_table":"pilot_scan_tables_mh125_20.csv",
        "not_tested":"experimental exclusions and applicable displaced-photon selection"}
    (CAMPAIGN/"physical_constraints_report.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    with (CAMPAIGN/"pilot_scan_tables.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=columns);writer.writeheader()
        for row,record in zip(rows,assessments):
            output={k:row.get(k,"") for k in columns}
            for name in ("quartic_screen_4pi","BFB_general","global_EW_minimum","full_tree_unitarity","local_minimum_tree","alignment"):
                output[name]=record["constraints"][name]["status"]
            output["candidate_status"]=record["candidate_status"]
            output["benchmark_certification"]=record["benchmark_certification"]
            output["exact_failure_reason"]=record["exact_failure_reason"]
            writer.writerow(output)
    boundary=[]
    for dim,transform in (("tan_beta",lambda x:Decimal(x)),("lambda1_ref",lambda x:Decimal(x)),
                          ("mphi_GeV",lambda x:Decimal(x))):
        groups=defaultdict(list)
        for row,record in zip(rows,assessments):
            key=(row["slice"],row["trajectory"],*(row[k] for k in ("tan_beta","lambda1_ref","mphi_GeV") if k!=dim))
            groups[key].append((transform(row[dim]),record))
        for key,items in groups.items():
            items.sort(key=lambda pair:pair[0])
            for (x,a),(y,b) in zip(items,items[1:]):
                sa,sb=a["candidate_status"],b["candidate_status"]
                if sa!=sb:
                    boundary.append({"axis":dim,"slice_trajectory_and_fixed_coordinates":key,
                        "sample_a":{"coordinate":str(x),"point_id":a["point_id"],"status":sa,"failure_reason":a["exact_failure_reason"]},
                        "sample_b":{"coordinate":str(y),"point_id":b["point_id"],"status":sb,"failure_reason":b["exact_failure_reason"]},
                        "refinement_status":"NOT_PERFORMED","reason":"candidate-valid class requires simultaneous certified physical conditions; no validated calculator/decay error budget exists for interpolated points"})
    analysis={"status":"PARTIAL","sampled_adjacent_status_changes":boundary,
        "sampled_status_change_count":len(boundary),"boundary_brackets_are_certified":False,
        "unsampled_regions":"not classified","physical_regions":"only exact sufficient certificates and exact negative witnesses receive PASS/FAIL; remaining points are INCONCLUSIVE",
        "refinement":"no interpolated refinement performed while backend/decay agreement error budgets remain unresolved"}
    (CAMPAIGN/"physical_boundary_analysis.json").write_text(json.dumps(analysis,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"points":len(assessments),"classes":dict(counts),"boundary_status_changes":len(boundary)},indent=2))


if __name__=="__main__":main()
