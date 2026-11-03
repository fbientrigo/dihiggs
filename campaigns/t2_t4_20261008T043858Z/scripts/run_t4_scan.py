#!/usr/bin/env python3
"""Serial, resumable 50/80/120-digit calculator comparison over the pilots."""
from __future__ import annotations

import argparse
import csv
from collections import Counter,defaultdict
import hashlib
import json
import os
import resource
import sys
import time
from pathlib import Path
import tempfile

HERE=Path(__file__).resolve().parent
CAMPAIGN=HERE.parent
REPO=CAMPAIGN.parents[1]
sys.path.insert(0,str(REPO))
sys.path.insert(0,str(HERE))
from dihiggs.app.orchestrator.aligned_lambda7 import build_aligned_lambda7_family
from run_t4_bridge import DECLARED_ACCEPTANCE_TOLERANCES, MODES, PRECISIONS, run_point
from scan_lock import campaign_scan_lock

OUTPUT=CAMPAIGN/"calculator_scan_agreement.jsonl"
CHECKPOINT=CAMPAIGN/"calculator_scan_checkpoint.json"
SCAN_METHOD_VERSION="calculator-scan-agreement-v4-mapped-native-decay-audit-locked"
ROW_KEYS=("point_id","mass_convention","slice","trajectory")
INPUT_KEYS=("tan_beta","lambda1_ref","mh_GeV","mphi_GeV","mA_GeV","mHp_GeV","v_GeV","evaluated_X","lambda6")


def canonical_sha(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=True).encode()).hexdigest()


def row_inputs(row):
    return {key:str(row.get(key,"")) for key in INPUT_KEYS}


def row_fingerprint(row):
    return canonical_sha({"identity":{key:str(row[key]) for key in ROW_KEYS},"inputs":row_inputs(row)})


def make_scan_identity(rows):
    source_files=[HERE/"run_t4_scan.py",HERE/"run_t4_bridge.py",HERE/"scan_lock.py",HERE/"analytic_bridge",
        CAMPAIGN/"build_analytic/2hdmc/Makefile",CAMPAIGN/"build_analytic/2hdmc/lib/lib2HDMC.a",
        REPO/"dihiggs/app/orchestrator/aligned_lambda7.py",
        REPO/"dihiggs/app/orchestrator/aligned_constraints.py",
        REPO/"dihiggs/app/orchestrator/analytic_bridge.py",
        REPO/"2hdmc/src/THDM.cpp",REPO/"2hdmc/src/THDM.h",
        REPO/"2hdmc/src/DecayTable.cpp",REPO/"2hdmc/src/DecayTable.h",
        REPO/"2hdmc/src/SM.cpp",REPO/"2hdmc/src/SM.h",
        REPO/"2hdmc/src/AnalyticBridge.cpp"]
    source_hashes={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
    grid_files=[CAMPAIGN/"t2_lightweight/trajectory_points.jsonl",CAMPAIGN/"pilot_scan_tables_mh125_20.csv"]
    grid_hashes={str(p.relative_to(CAMPAIGN)):hashlib.sha256(p.read_bytes()).hexdigest() for p in grid_files}
    spec={"method_version":SCAN_METHOD_VERSION,"source_hashes":source_hashes,
        "calculator_executable_sha256":source_hashes[str((HERE/"analytic_bridge").relative_to(REPO))],
        "grid_hashes":grid_hashes,"row_fingerprints":[row_fingerprint(row) for row in rows],
        "precision_digits":list(PRECISIONS),"constructor_modes":list(MODES),
        "acceptance_tolerances":DECLARED_ACCEPTANCE_TOLERANCES}
    return {**spec,"fingerprint":canonical_sha(spec)}


def scan_rows():
    light=CAMPAIGN/"t2_lightweight/trajectory_points.jsonl"
    rows=[json.loads(x) for x in light.read_text().splitlines() if x.strip()]
    for row in rows:
        row["mass_convention"]="mh125.13_T1"
    active_path=CAMPAIGN/"pilot_scan_tables_mh125_20.csv"
    with active_path.open(newline="") as f:
        active=list(csv.DictReader(f))
    for row in active:
        row["mass_convention"]="mh125.20_active"
    return rows+active


def build_point(row,precision):
    args={key:row[key] for key in ("tan_beta","lambda1_ref","mh_GeV","mphi_GeV","mA_GeV","mHp_GeV","v_GeV")}
    if row["trajectory"]=="fixed_lambda6_1e-10":
        return build_aligned_lambda7_family(**args,lambda6=row["lambda6"],precision=precision)
    return build_aligned_lambda7_family(**args,X_lambda6_tanbeta=row["evaluated_X"],precision=precision)


def run_one(row,identity,precisions=PRECISIONS):
    def builder(digits):
        return build_point(row,digits)
    records=run_point(point_builder=builder,precisions=precisions)
    for record in records:
        record["point_id"]=row["point_id"]
        record["mass_convention"]=row["mass_convention"]
        record["pilot_slice"]=row["slice"]
        record["trajectory"]=row["trajectory"]
        record["pilot_point_metadata"]={key:row.get(key) for key in
            ("tan_beta","lambda1_ref","mh_GeV","mphi_GeV","mA_GeV","mHp_GeV","v_GeV","evaluated_X","lambda6")}
        record["scan_identity"]=identity["fingerprint"]
        record["point_input_sha256"]=row_fingerprint(row)
        record["record_payload_sha256"]=canonical_sha(record)
    return records


def profile(row,identity):
    start=time.perf_counter()
    before=resource.getrusage(resource.RUSAGE_SELF)
    records=run_one(row,identity)
    after=resource.getrusage(resource.RUSAGE_SELF)
    elapsed=time.perf_counter()-start
    profile_record={"status":"PROFILED_SERIAL_ONE_POINT","point_id":row["point_id"],
        "mass_convention":row["mass_convention"],"route_precision_invocations":len(records),
        "scan_identity":identity,
        "worker_count":1,"wall_seconds":elapsed,
        "user_cpu_seconds":after.ru_utime-before.ru_utime,
        "system_cpu_seconds":after.ru_stime-before.ru_stime,
        "process_peak_rss_bytes":after.ru_maxrss*1024,
        "measured_seconds_per_route_precision":elapsed/max(1,len(records)),
        "estimated_full_scan_wall_seconds_serial":elapsed*390,
        "resource_policy":{"serial_only":True,"safe_mem_ceiling":"11G","safe_high_ceiling":"9G","safe_cpu_ceiling":"900%"},
        "records_summary":[{"precision":r["precision_digits_constructed"],"mode":r["mode"],
            "constructor_status":r.get("constructor_status"),"exact_failure_reason":r.get("exact_failure_reason","")}
            for r in records]}
    (CAMPAIGN/"calculator_scan_profile.json").write_text(json.dumps(profile_record,indent=2,sort_keys=True)+"\n")
    return profile_record


def complete_point_count(done,rows):
    return sum(all((r["point_id"],r["mass_convention"],p,m) in done
        for p in PRECISIONS for m in MODES) for r in rows)


def load_existing(rows,identity,output_path,checkpoint_path):
    expected_rows={(str(r["point_id"]),str(r["mass_convention"])):r for r in rows}
    if len(expected_rows)!=len(rows):
        raise ValueError("pilot grid contains duplicate point/convention keys")
    checkpoint=None
    if checkpoint_path.exists():
        checkpoint=json.loads(checkpoint_path.read_text())
        if checkpoint.get("scan_identity")!=identity["fingerprint"] or checkpoint.get("identity")!=identity:
            raise ValueError("checkpoint scan identity mismatch: sources, executable, grid, precision, mode or tolerance changed")
    done=set()
    if not output_path.exists():
        return done
    for line_number,line in enumerate(output_path.read_text().splitlines(),1):
        try:
            record=json.loads(line)
        except json.JSONDecodeError as exc:
            raise ValueError(f"malformed scan record at line {line_number}") from exc
        required=("point_id","mass_convention","precision_digits_constructed","mode","scan_identity","point_input_sha256","pilot_point_metadata","record_payload_sha256")
        for key in required:
            if key not in record:
                raise ValueError(f"scan record line {line_number} lacks {key}")
        pointkey=(str(record["point_id"]),str(record["mass_convention"]))
        row=expected_rows.get(pointkey)
        if row is None:
            raise ValueError(f"unknown point/convention key at line {line_number}: {pointkey}")
        if record["scan_identity"]!=identity["fingerprint"]:
            raise ValueError(f"scan record line {line_number} source/executable/grid identity mismatch")
        if record["point_input_sha256"]!=row_fingerprint(row):
            raise ValueError(f"scan record line {line_number} input-grid fingerprint mismatch")
        if record["pilot_point_metadata"]!=row_inputs(row):
            raise ValueError(f"scan record line {line_number} inputs do not match the current grid")
        checksum=record["record_payload_sha256"]
        payload={key:value for key,value in record.items() if key!="record_payload_sha256"}
        if checksum!=canonical_sha(payload):
            raise ValueError(f"scan record line {line_number} scientific payload checksum mismatch")
        precision=record["precision_digits_constructed"]
        mode=record["mode"]
        if precision not in PRECISIONS or mode not in MODES:
            raise ValueError(f"unknown precision or constructor mode at line {line_number}")
        common_payload=("constructor_status","constructor_ok","constructor_route","exact_failure_reason",
            "inputs_decimal","target_values_decimal","requested_mode","declared_acceptance_tolerances")
        missing=[name for name in common_payload if name not in record]
        if missing:
            raise ValueError(f"scan record line {line_number} has incomplete result payload: {missing}")
        if (not isinstance(record["inputs_decimal"],dict) or not record["inputs_decimal"] or
                not isinstance(record["target_values_decimal"],dict) or not record["target_values_decimal"] or
                not isinstance(record["declared_acceptance_tolerances"],dict) or not record["declared_acceptance_tolerances"] or
                not isinstance(record["exact_failure_reason"],str)):
            raise ValueError(f"scan record line {line_number} has malformed result payload values")
        if record["constructor_status"] not in ("PASS","FAIL","INCONCLUSIVE","NOT_TESTED"):
            raise ValueError(f"scan record line {line_number} has invalid constructor status")
        if record["constructor_ok"]:
            required_success=("constructor_checks","backend_abs_residuals","independent_full_unitarity",
                "backend_spectrum_block_residuals",*(f"lambda{i}" for i in range(1,8)),
                *(f"Z{i}" for i in range(1,8)),"m12_sq_GeV2","backend_eigenvalues_Y2_sigma1",
                "backend_eigenvalues_Y2_sigma0","backend_eigenvalues_Y0_sigma1","backend_eigenvalues_Y0_sigma0",
                "backend_v_GeV","m_h_backend_GeV","m_H_backend_GeV","m_A_backend_GeV","m_Hp_backend_GeV",
                "g_hphi_phi_backend_GeV","g_phi_HpHm_backend_GeV","kappa_phi_top_backend_signed")
            if record.get("physical_state_ordering")=="MAPPED" and record.get("widths_computed"):
                required_success += ("width_total_phi_GeV","width_channels_GeV",
                    "branching_fractions_reported_channels","physical_lifetime")
            else:
                required_success += ("widths","branching_fractions","lifetime")
            missing=[name for name in required_success if name not in record]
            if missing:
                raise ValueError(f"scan record line {line_number} has incomplete successful scientific outputs: {missing}")
            if (not isinstance(record["constructor_checks"],dict) or not record["constructor_checks"] or
                    not isinstance(record["backend_abs_residuals"],dict) or
                    not isinstance(record["independent_full_unitarity"],dict) or
                    not isinstance(record["backend_spectrum_block_residuals"],dict)):
                raise ValueError(f"scan record line {line_number} has malformed successful scientific outputs")
        elif not record["exact_failure_reason"]:
            raise ValueError(f"scan record line {line_number} has no exact failure reason")
        key=(pointkey[0],pointkey[1],precision,mode)
        if key in done:
            raise ValueError(f"duplicate scan record at line {line_number}: {key}")
        done.add(key)
    if checkpoint and checkpoint.get("completed_route_precision_records",0)>len(done):
        raise ValueError("checkpoint claims more completed records than the JSONL contains")
    return done


def write_checkpoint(identity,done,rows,batch_size,last_row_index,status="RUNNING"):
    payload={"status":status,"scan_identity":identity["fingerprint"],"identity":identity,
        "completed_route_precision_records":len(done),
        "expected_route_precision_records":len(rows)*len(PRECISIONS)*len(MODES),
        "completed_pilot_points":complete_point_count(done,rows),"pilot_point_count":len(rows),
        "worker_count":1,"last_row_index":last_row_index,"batch_size":batch_size,
        "output_jsonl":OUTPUT.name}
    tmp=CHECKPOINT.with_suffix(".tmp")
    tmp.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    tmp.replace(CHECKPOINT)


def write_summary(rows,identity,done):
    records=[json.loads(line) for line in OUTPUT.read_text().splitlines() if line.strip()]
    expected=len(rows)*len(PRECISIONS)*len(MODES)
    if len(records)!=expected or len(done)!=expected:
        raise ValueError(f"cannot summarize incomplete scan: {len(records)}/{expected} records")
    by_status={mode:dict(Counter(r["constructor_status"] for r in records if r["mode"]==mode)) for mode in MODES}
    by_precision={mode:{str(digits):dict(Counter(r["constructor_status"] for r in records
        if r["mode"]==mode and r["precision_digits_constructed"]==digits)) for digits in PRECISIONS} for mode in MODES}
    reasons=Counter(r.get("exact_failure_reason","") for r in records if r.get("constructor_status") in ("FAIL","INCONCLUSIVE"))
    max_residuals={}
    spectrum_residuals={}
    spreads={}
    grouped=defaultdict(list)
    spread_keys=(*(f"lambda{i}" for i in range(1,8)),*(f"Z{i}" for i in range(1,8)),
        "m12_sq_GeV2","m_h_backend_GeV","m_H_backend_GeV","m_A_backend_GeV","m_Hp_backend_GeV",
        "backend_v_GeV","width_total_phi_GeV","width_diphoton_GeV","ctau_mm",
        "g_hphi_phi_backend_GeV","g_phi_HpHm_backend_GeV","kappa_phi_top_backend_signed")
    for record in records:
        route=(record["mode"],record["precision_digits_constructed"])
        max_residuals.setdefault(str(route),{})
        for key,value in record.get("backend_abs_residuals",{}).items():
            max_residuals[str(route)][key]=max(max_residuals[str(route)].get(key,0.0),float(value))
        spectrum_residuals.setdefault(str(route),{})
        for block,value in record.get("backend_spectrum_block_residuals",{}).items():
            spectrum_residuals[str(route)][block]=max(spectrum_residuals[str(route)].get(block,0.0),float(value))
        grouped[(record["point_id"],record["mass_convention"],record["mode"])].append(record)
    for (point_id,convention,mode),items in grouped.items():
        if len(items)!=len(PRECISIONS):
            raise ValueError(f"precision convergence rows missing for {point_id}/{convention}/{mode}")
        if len({item["precision_digits_constructed"] for item in items})!=len(PRECISIONS):
            raise ValueError(f"duplicate precision convergence rows for {point_id}/{convention}/{mode}")
        spreads.setdefault(mode,{})
        for key in spread_keys:
            values=[float(item[key]) for item in items if isinstance(item.get(key),(float,int))]
            if len(values)==len(items):
                spreads[mode][key]=max(spreads[mode].get(key,0.0),max(values)-min(values))
    payload={"status":"COMPLETE","scan_identity":identity,
        "pilot_point_count":len(rows),"route_precision_record_count":len(records),
        "point_convention_counts":dict(Counter(r["mass_convention"] for r in rows)),
        "constructor_status_counts_by_route":by_status,
        "constructor_status_counts_by_route_and_precision":by_precision,
        "exact_failure_reason_counts":dict(reasons),
        "maximum_backend_absolute_residual_by_route_precision":max_residuals,
        "maximum_actual_backend_spectrum_residual_by_route_precision":spectrum_residuals,
        "maximum_backend_output_spread_over_50_80_120_digits_by_route":spreads,
        "acceptance_tolerances_frozen_in_scan_identity":DECLARED_ACCEPTANCE_TOLERANCES,
        "record_schema_note":"Every route/precision row contains decimal constructor inputs, backend double inputs, complete actual-backend and independent spectra, per-observable residuals/tolerances, widths/BR/lifetime where state ordering is supported, status and exact failure reason.",
        "output_jsonl":OUTPUT.name,"output_sha256":hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
        "checkpoint_json":CHECKPOINT.name,"worker_count":1}
    (CAMPAIGN/"calculator_scan_summary.json").write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")


def execute_scan(args,rows,identity):
    if args.profile_only:
        pivot=next(r for r in rows if r["mass_convention"]=="mh125.13_T1" and
            r["trajectory"]=="fixed_X_1" and r["tan_beta"]=="1000000" and r["lambda1_ref"]=="1.2")
        print(json.dumps(profile(pivot,identity),indent=2,sort_keys=True))
        return
    done=load_existing(rows,identity,OUTPUT,CHECKPOINT)
    expected={(r["point_id"],r["mass_convention"],p,m) for r in rows for p in PRECISIONS for m in MODES}
    completed=len(done & expected)
    with OUTPUT.open("a",buffering=1) as output:
        pending=0
        for row_index,row in enumerate(rows):
            missing=[p for p in PRECISIONS if any((row["point_id"],row["mass_convention"],p,m) not in done for m in MODES)]
            if not missing:
                continue
            records=run_one(row,identity,missing)
            for record in records:
                key=(record["point_id"],record["mass_convention"],record["precision_digits_constructed"],record["mode"])
                if key in done:
                    continue
                output.write(json.dumps(record,sort_keys=True)+"\n")
                done.add(key);completed+=1
            output.flush();os.fsync(output.fileno())
            pending+=1
            if pending>=args.batch_size:
                write_checkpoint(identity,done,rows,args.batch_size,row_index)
                pending=0
                print(f"checkpoint: {completed}/{len(expected)} route/precision records",flush=True)
        write_checkpoint(identity,done,rows,args.batch_size,len(rows)-1,
            status="COMPLETE" if expected<=done else "PARTIAL")
    print(json.dumps({"completed":completed,"expected":len(expected),"output":str(OUTPUT),
        "status":"COMPLETE" if expected<=done else "PARTIAL"},indent=2))
    if expected<=done:
        write_summary(rows,identity,done)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--profile-only",action="store_true")
    parser.add_argument("--batch-size",type=int,default=10)
    args=parser.parse_args()
    if args.batch_size<1:
        raise SystemExit("--batch-size must be positive")
    rows=scan_rows()
    if len(rows)!=390:
        raise SystemExit(f"expected 390 pilot rows, got {len(rows)}")
    identity=make_scan_identity(rows)
    lock_path=Path(tempfile.gettempdir())/f"{CAMPAIGN.name}.calculator-scan.lock"
    with campaign_scan_lock(lock_path,identity["fingerprint"]):
        execute_scan(args,rows,identity)


if __name__=="__main__":
    main()
