#!/usr/bin/env python3
"""Campaign-only audit of complete native 2HDMC decay-width accounting."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import resource
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

HERE=Path(__file__).resolve().parent
CAMPAIGN=HERE.parent
REPO=CAMPAIGN.parents[1]
OUT=CAMPAIGN/"decay_accounting"
EXECUTABLE=OUT/"DecayAccountingAudit"
RESULTS=OUT/"results.jsonl"
CHECKPOINT=OUT/"checkpoint.json"
SUMMARY=OUT/"summary.json"
PROFILE=OUT/"profile.json"
METHOD_VERSION="native-decay-accounting-v1-78-terms"
HBAR_C_GEV_MM=1.973269804e-13
EPS=2.220446049250313e-16
EXPECTED_FAMILIES={"hdd":9,"huu":9,"hdu":9,"hll":9,"hln":9,
    "hvv":3,"hZga":1,"hgg":1,"hvh":12,"hhh":16}

sys.path.insert(0,str(HERE))
import run_t4_scan as frozen_scan
from scan_lock import campaign_scan_lock


def sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_sha(value) -> str:
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),
        ensure_ascii=True,allow_nan=False).encode()).hexdigest()


def expected_term_keys():
    keys=[]
    for i in range(1,4):
        for j in range(1,4):
            for family in ("hdd","huu","hdu","hll","hln"):
                keys.append((family,i,j))
    keys.extend(("hvv",i,0) for i in range(1,4))
    keys.extend((("hZga",0,0),("hgg",0,0)))
    keys.extend(("hvh",i,j) for i in range(1,4) for j in range(1,5))
    keys.extend(("hhh",i,j) for i in range(1,5) for j in range(1,5))
    return keys


def _within_ulps(a: float,b: float,ulps: int=64) -> tuple[bool,float,float]:
    if not math.isfinite(a) or not math.isfinite(b):
        return False,math.inf,0.0
    residual=abs(a-b)
    scale=max(abs(a),abs(b))
    if scale==0.0:
        return residual==0.0,residual,0.0
    tolerance=ulps*EPS*scale
    return residual<=tolerance,residual,tolerance


def audit_native_payload(native: dict,frozen: dict,*,audit_identity: str="fixture") -> dict:
    """Check channel completeness, native sum/BR closure and conversion identities."""
    result={"audit_identity":audit_identity,"scan_identity":frozen.get("scan_identity"),
        "point_id":frozen.get("point_id"),"mass_convention":frozen.get("mass_convention"),
        "precision_digits":frozen.get("precision_digits_constructed"),"mode":"analytic",
        "source_record_payload_sha256":frozen.get("record_payload_sha256"),
        "validation_scope":"native DecayTable accounting and algebraic closure only; no independent non-diphoton formula audit",
        "benchmark_physical_validation_status":"INCONCLUSIVE","checks":{},"exact_failure_reason":""}
    if native.get("constructor_ok") is not True:
        result.update(status="INCONCLUSIVE",exact_failure_reason=native.get("reason","native constructor failed"))
        return result
    if native.get("state_identity_ok") is not True:
        result.update(status="INCONCLUSIVE",exact_failure_reason=native.get("reason","state identity unresolved"))
        return result
    if frozen.get("constructor_status")!="PASS" or frozen.get("physical_state_ordering")!="MAPPED":
        result.update(status="INCONCLUSIVE",exact_failure_reason="frozen calculator constructor/state mapping did not PASS")
        return result
    gamma_audit=frozen.get("independent_diphoton_components",{})
    if gamma_audit.get("status")!="PASS":
        result.update(status="INCONCLUSIVE",exact_failure_reason="frozen native diphoton component audit did not PASS")
        return result
    for name,native_key,frozen_key,tolerance in (
        ("requested_h_index","h_index","requested_h_index_backend",0.0),
        ("requested_phi_index","phi_index","requested_phi_index_backend",0.0),
        ("mapped_h_mass","m_h_index_GeV","m_h_backend_GeV",1e-8),
        ("mapped_phi_mass","m_phi_index_GeV","m_H_backend_GeV",1e-8),
        ("backend_vev","backend_v_GeV","backend_v_GeV",1e-10),
        ("requested_vev","requested_v_GeV","requested_v_GeV",0.0)):
        actual=native.get(native_key); expected_value=frozen.get(frozen_key)
        if isinstance(actual,bool) or not isinstance(actual,(int,float)) or not math.isfinite(actual) or expected_value is None:
            result.update(status="INCONCLUSIVE",exact_failure_reason=f"missing/nonfinite frozen state identity field {name}")
            return result
        if tolerance==0.0:
            agrees=int(actual)==int(expected_value) if "index" in name else float(actual)==float(expected_value)
            residual=abs(float(actual)-float(expected_value))
        else:
            agrees,residual,_=_within_ulps(float(actual),float(expected_value),ulps=0)
            agrees=agrees or residual<=tolerance
        if not agrees:
            result.update(status="INCONCLUSIVE",exact_failure_reason=f"reconstructed {name} differs from frozen analytic record")
            result["state_identity_residual"]={"observable":name,"absolute_residual":residual,"declared_tolerance":tolerance}
            return result
    result["checks"]["frozen_state_identity"]={"status":"PASS",
        "method":"requested mapped indices, masses, and requested/backend electroweak scales match the frozen analytic 80-digit record"}

    terms=native.get("terms")
    expected=expected_term_keys()
    keys=[]
    bad=[]
    if not isinstance(terms,list):
        result.update(status="FAIL",exact_failure_reason="native term table is missing or malformed")
        return result
    for item in terms:
        if not isinstance(item,dict):
            bad.append("malformed_term_object")
            continue
        keys.append((item.get("family"),item.get("i"),item.get("j")))
        value=item.get("width_GeV")
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:
            bad.append(f"invalid_width:{item.get('family')}:{item.get('i')}:{item.get('j')}")
    if keys!=expected:
        result["checks"]["complete_indexed_channel_schema"]={"status":"FAIL",
            "expected_term_count":len(expected),"actual_term_count":len(keys),
            "missing_or_reordered_terms":True}
        result["exact_failure_reason"]="native indexed channel table must contain all 78 terms in DecayTable summation order"
        result["status"]="FAIL"
        return result
    if bad:
        result.update(status="FAIL",exact_failure_reason="; ".join(bad))
        return result

    family_terms={family:[] for family in EXPECTED_FAMILIES}
    for item in terms:
        family_terms[item["family"]].append(float(item["width_GeV"]))
    expected_counts=EXPECTED_FAMILIES
    actual_counts={key:len(value) for key,value in family_terms.items()}
    result["checks"]["complete_indexed_channel_schema"]={"status":"PASS",
        "expected_term_count":78,"actual_term_count":len(terms),
        "family_term_counts":actual_counts}

    reported_families=native.get("family_totals_GeV")
    if not isinstance(reported_families,dict) or set(reported_families)!=set(EXPECTED_FAMILIES):
        result.update(status="FAIL",exact_failure_reason="native family total schema is incomplete")
        return result
    family_residuals={}
    all_checks=["complete_indexed_channel_schema"]
    for family,values in family_terms.items():
        if actual_counts[family]!=expected_counts[family]:
            result.update(status="FAIL",exact_failure_reason=f"native family {family} has incomplete indexed terms")
            return result
        subtotal=0.0
        for value in values:
            subtotal+=value
        ok,residual,tolerance=_within_ulps(subtotal,float(reported_families[family]),ulps=64)
        family_residuals[family]={"absolute_residual_GeV":residual,"tolerance_GeV":tolerance}
        if not ok:
            all_checks.append(f"family_sum_{family}")
            result["checks"][f"family_sum_{family}"]={"status":"FAIL",**family_residuals[family]}
        else:
            result["checks"][f"family_sum_{family}"]={"status":"PASS",**family_residuals[family]}
        all_checks.append(f"family_sum_{family}")

    ordered_sum=0.0
    for item in terms:
        ordered_sum+=float(item["width_GeV"])
    total=native.get("total_width_GeV")
    native_ordered=native.get("ordered_term_sum_GeV")
    frozen_total=frozen.get("width_total_phi_GeV")
    diphoton=native.get("diphoton_width_GeV")
    frozen_diphoton=frozen.get("width_diphoton_GeV")
    ctau=native.get("ctau_mm")
    if any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or x<0
           for x in (total,native_ordered,frozen_total,diphoton,frozen_diphoton)):
        result.update(status="FAIL",exact_failure_reason="total or diphoton width is nonfinite/negative")
        return result
    if isinstance(ctau,bool) or not isinstance(ctau,(int,float)) or not math.isfinite(ctau) or ctau<0:
        result.update(status="FAIL",exact_failure_reason="lifetime is nonfinite/negative")
        return result
    result["total_width_GeV"]=total
    result["ordered_channel_sum_GeV"]=ordered_sum
    result["diphoton_width_GeV"]=diphoton
    result["ctau_mm"]=ctau
    result["complete_branching_fractions"]={
        f"{term['family']}:{term['i']}:{term['j']}":float(term["width_GeV"])/total if total>0 else None
        for term in terms}
    result["family_branching_fractions"]={family:float(reported_families[family])/total if total>0 else None
        for family in EXPECTED_FAMILIES}
    if total<=0:
        result.update(status="FAIL",exact_failure_reason="total decay width must be positive for BR/lifetime validation")
        return result

    for name,a,b in (("python_term_sum_vs_native_ordered_sum",ordered_sum,float(native_ordered)),
                     ("native_total_vs_all_terms",ordered_sum,total),
                     ("native_total_vs_frozen_t4",total,float(frozen_total)),
                     ("native_diphoton_vs_frozen_t4",diphoton,float(frozen_diphoton))):
        ok,residual,tolerance=_within_ulps(float(a),float(b),ulps=64)
        result["checks"][name]={"status":"PASS" if ok else "FAIL",
            "absolute_residual_GeV":residual,"tolerance_GeV":tolerance}
        all_checks.append(name)
    br_sum=sum(result["complete_branching_fractions"].values())
    br_tolerance=8*len(terms)*EPS
    result["checks"]["complete_branching_fraction_normalization"]={
        "status":"PASS" if abs(br_sum-1.0)<=br_tolerance else "FAIL",
        "sum":br_sum,"absolute_residual":abs(br_sum-1.0),"tolerance":br_tolerance}
    all_checks.append("complete_branching_fraction_normalization")
    expected_ctau=HBAR_C_GEV_MM/float(total)
    ok,residual,tolerance=_within_ulps(float(ctau),expected_ctau,ulps=8)
    result["checks"]["lifetime_conversion"]={"status":"PASS" if ok else "FAIL",
        "expected_hbar_c_GeV_mm":HBAR_C_GEV_MM,"formula_ctau_mm":expected_ctau,
        "absolute_residual_mm":residual,"tolerance_mm":tolerance}
    all_checks.append("lifetime_conversion")
    result["checks"]["native_diphoton_components"]={"status":"PASS",
        "source_record_status":gamma_audit.get("status"),
        "scheme_note":"native DecayTable running-mass/QCD component audit; appendix pole-mass LO comparison remains a separate scheme diagnostic"}
    all_checks.append("native_diphoton_components")
    failed=[name for name,check in result["checks"].items() if check["status"]!="PASS"]
    result["status"]="FAIL" if failed else "PASS"
    result["exact_failure_reason"]="; ".join(failed)
    result["decay_status"]={"total_width":"PASS" if not failed else "FAIL",
        "branching_fractions":"PASS" if not failed else "FAIL",
        "lifetime":"PASS" if not failed else "FAIL",
        "scope":"CALCULATOR_ACCOUNTING_ONLY",
        "independent_non_diphoton_formulas":"NOT_TESTED",
        "benchmark_physical_validation":"INCONCLUSIVE"}
    return result


def load_frozen_scan():
    rows=frozen_scan.scan_rows()
    if len(rows)!=390:
        raise ValueError(f"frozen T4 pilot grid has {len(rows)} rows, expected 390")
    checkpoint=json.loads(frozen_scan.CHECKPOINT.read_text())
    identity=checkpoint.get("identity")
    if not isinstance(identity,dict) or identity.get("fingerprint")!=checkpoint.get("scan_identity"):
        raise ValueError("frozen T4 checkpoint identity is malformed")
    spec={key:value for key,value in identity.items() if key!="fingerprint"}
    if canonical_sha(spec)!=identity["fingerprint"]:
        raise ValueError("frozen T4 identity fingerprint does not match its recorded specification")
    for relative,expected in identity.get("source_hashes",{}).items():
        source=REPO/relative
        if not source.is_file() or sha(source)!=expected:
            raise ValueError(f"frozen T4 source/executable changed since scan: {relative}")
    done=frozen_scan.load_existing(rows,identity,frozen_scan.OUTPUT,frozen_scan.CHECKPOINT)
    expected={(str(row["point_id"]),str(row["mass_convention"]),precision,mode)
        for row in rows for precision in frozen_scan.PRECISIONS for mode in frozen_scan.MODES}
    if len(done)!=3510 or done!=expected or checkpoint.get("status")!="COMPLETE":
        raise ValueError("frozen T4 scan is not a complete 3510-record replay")
    records=[json.loads(line) for line in frozen_scan.OUTPUT.read_text().splitlines() if line.strip()]
    selected=[record for record in records if record.get("mode")=="analytic" and
        record.get("precision_digits_constructed")==80]
    if len(selected)!=390:
        raise ValueError(f"expected 390 analytic/80-digit records, found {len(selected)}")
    keys=set()
    for record in selected:
        key=(str(record["point_id"]),str(record["mass_convention"]))
        if key in keys:
            raise ValueError(f"duplicate analytic/80 record for {key}")
        keys.add(key)
        if record.get("constructor_status")!="PASS" or record.get("physical_state_ordering")!="MAPPED":
            raise ValueError(f"frozen constructor/mapping is not PASS for {key}")
        if record.get("independent_diphoton_components",{}).get("status")!="PASS":
            raise ValueError(f"frozen native diphoton audit is not PASS for {key}")
    return rows,identity,selected


def build_executable():
    OUT.mkdir(parents=True,exist_ok=True)
    build=CAMPAIGN/"build_analytic/2hdmc"
    command=["g++","-std=c++11","-Wall","-Wextra","-O0","-fopenmp",
        "-I",str(build/"src"),str(HERE/"DecayAccountingAudit.cpp"),
        "-L",str(build/"lib"),"-l2HDMC","-lgsl","-lgslcblas","-lm","-o",str(EXECUTABLE)]
    env=os.environ.copy(); env.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1")
    subprocess.run(command,cwd=REPO,env=env,check=True,timeout=60,preexec_fn=_preexec_limiter)
    return command


def executable_identity():
    build=CAMPAIGN/"build_analytic/2hdmc"
    compiler=subprocess.run(["g++","--version"],capture_output=True,text=True,check=True,timeout=10).stdout.splitlines()[0]
    files={
        "audit_cpp_sha256":sha(HERE/"DecayAccountingAudit.cpp"),
        "runner_sha256":sha(HERE/"run_decay_accounting_audit.py"),
        "isolated_library_sha256":sha(build/"lib/lib2HDMC.a"),
        "native_decay_source_sha256":sha(REPO/"2hdmc/src/DecayTable.cpp"),
        "native_decay_header_sha256":sha(REPO/"2hdmc/src/DecayTable.h"),
        "generic_model_source_sha256":sha(REPO/"2hdmc/src/THDM.cpp"),
        "generic_model_header_sha256":sha(REPO/"2hdmc/src/THDM.h"),
        "analytic_source_sha256":sha(REPO/"2hdmc/src/AnalyticBridge.cpp"),
    }
    if EXECUTABLE.is_file():
        files["audit_executable_sha256"]=sha(EXECUTABLE)
    return {"compiler":compiler,"files":files,"method_version":METHOD_VERSION,
        "build_command_template":["g++","-std=c++11","-Wall","-Wextra","-O0","-fopenmp",
            "-I <campaign>/build_analytic/2hdmc/src","<campaign>/scripts/DecayAccountingAudit.cpp",
            "-L <campaign>/build_analytic/2hdmc/lib","-l2HDMC","-lgsl","-lgslcblas","-lm"],
        "term_counts":EXPECTED_FAMILIES,"total_terms":78,
        "tolerances":{"dimensionless_branching_normalization":"8*N_terms*binary64_epsilon",
            "sum_and_replay_width_comparisons":"64*binary64_epsilon*relative_scale",
            "lifetime_conversion":"8*binary64_epsilon*relative_scale"},
        "resource_policy":{"workers":1,"threads":1,"RLIMIT_AS_bytes":9*1024**3,"RLIMIT_CPU_seconds":900}}


def scan_identity(scan_spec,selected,base):
    return {**base,"frozen_scan_identity":scan_spec["fingerprint"],
        "frozen_jsonl_sha256":sha(frozen_scan.OUTPUT),
        "frozen_checkpoint_sha256":sha(frozen_scan.CHECKPOINT),
        "selected_record_payload_sha256":[x["record_payload_sha256"] for x in selected],
        "selected_point_count":len(selected),
        "fingerprint":canonical_sha({**base,"frozen_scan_identity":scan_spec["fingerprint"],
            "frozen_jsonl_sha256":sha(frozen_scan.OUTPUT),
            "frozen_checkpoint_sha256":sha(frozen_scan.CHECKPOINT),
            "selected_record_payload_sha256":[x["record_payload_sha256"] for x in selected],
            "selected_point_count":len(selected)})}


def _preexec_limiter():
    resource.setrlimit(resource.RLIMIT_AS,(9*1024**3,9*1024**3))
    resource.setrlimit(resource.RLIMIT_CPU,(900,900))


def run_native(record):
    args=[str(EXECUTABLE)]
    args.extend(str(record[f"lambda{i}"]) for i in range(1,8))
    args.extend((str(record["m12_sq_GeV2"]),str(record["tan_beta_backend"]),
        str(record["requested_v_GeV"]),str(record["requested_h_index_backend"]),
        str(record["requested_phi_index_backend"]),str(record["m_h_backend_GeV"]),
        str(record["m_H_backend_GeV"]),str(record["backend_v_GeV"])))
    env=os.environ.copy(); env.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1")
    start=time.perf_counter()
    proc=subprocess.run(args,cwd=REPO,capture_output=True,text=True,check=True,timeout=20,
        env=env,preexec_fn=_preexec_limiter)
    elapsed=time.perf_counter()-start
    try:
        native=json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        raise ValueError(f"native audit emitted malformed JSON: {proc.stdout[:300]}") from exc
    return native,elapsed


def run_record(frozen,audit_id):
    native,elapsed=run_native(frozen)
    result=audit_native_payload(native,frozen,audit_identity=audit_id)
    result["calculator_elapsed_seconds"]=elapsed
    result["native_payload"]=native
    result["record_payload_sha256"]=canonical_sha(result)
    return result


def load_existing(selected,identity):
    if CHECKPOINT.exists():
        checkpoint=json.loads(CHECKPOINT.read_text())
        if checkpoint.get("audit_identity")!=identity:
            raise ValueError("decay-accounting checkpoint identity mismatch")
        if checkpoint.get("status") not in ("RUNNING","PARTIAL","COMPLETE"):
            raise ValueError("decay-accounting checkpoint has invalid status")
    else:
        checkpoint=None
    expected={(str(x["point_id"]),str(x["mass_convention"])):x for x in selected}
    if len(expected)!=len(selected):
        raise ValueError("frozen audit input has duplicate point/convention keys")
    done=set()
    if not RESULTS.exists():
        if checkpoint and (checkpoint.get("completed_records",0)>0 or checkpoint.get("status")=="COMPLETE"):
            raise ValueError("checkpoint claims completed records but the output file is missing")
        return done
    if checkpoint and checkpoint.get("output_sha256") is not None:
        data=RESULTS.read_bytes()
        prefix_bytes=checkpoint.get("output_bytes")
        if not isinstance(prefix_bytes,int) or prefix_bytes<0 or prefix_bytes>len(data):
            raise ValueError("checkpoint output byte boundary is invalid")
        prefix=data[:prefix_bytes]
        if prefix and not prefix.endswith(b"\n"):
            raise ValueError("checkpoint output boundary does not end at a complete JSONL record")
        if hashlib.sha256(prefix).hexdigest()!=checkpoint["output_sha256"]:
            raise ValueError("checkpoint output prefix checksum mismatch")
    for line_number,line in enumerate(RESULTS.read_text().splitlines(),1):
        try: row=json.loads(line)
        except json.JSONDecodeError as exc: raise ValueError(f"malformed audit record on line {line_number}") from exc
        key=(str(row.get("point_id")),str(row.get("mass_convention")))
        source=expected.get(key)
        if source is None: raise ValueError(f"unknown audit key at line {line_number}: {key}")
        if row.get("audit_identity")!=identity["fingerprint"]: raise ValueError(f"audit identity mismatch on line {line_number}")
        if row.get("source_record_payload_sha256")!=source.get("record_payload_sha256"):
            raise ValueError(f"frozen source-record mismatch on line {line_number}")
        checksum=row.get("record_payload_sha256")
        payload={k:v for k,v in row.items() if k!="record_payload_sha256"}
        if not isinstance(checksum,str) or checksum!=canonical_sha(payload):
            raise ValueError(f"audit payload checksum mismatch on line {line_number}")
        if row.get("status") not in ("PASS","FAIL","INCONCLUSIVE","NOT_TESTED"):
            raise ValueError(f"invalid audit status on line {line_number}")
        if not isinstance(row.get("checks"),dict) or not isinstance(row.get("native_payload"),dict):
            raise ValueError(f"incomplete accounting payload on line {line_number}")
        expected_result=audit_native_payload(row["native_payload"],source,
            audit_identity=identity["fingerprint"])
        for name,value in expected_result.items():
            if row.get(name)!=value:
                raise ValueError(f"audit-derived result field {name} disagrees with native payload on line {line_number}")
        if key in done: raise ValueError(f"duplicate audit record at line {line_number}: {key}")
        done.add(key)
    if checkpoint and checkpoint.get("completed_records",0)>len(done):
        raise ValueError("checkpoint claims records missing from JSONL")
    if checkpoint and checkpoint.get("status")=="COMPLETE" and len(done)!=len(selected):
        raise ValueError("checkpoint is COMPLETE but the output is missing expected records")
    if checkpoint and checkpoint.get("status")=="COMPLETE" and checkpoint.get("completed_records")!=len(done):
        raise ValueError("COMPLETE checkpoint record count differs from JSONL")
    return done


def write_checkpoint(identity,done,selected,batch_size,last_index,status):
    output_bytes=RESULTS.stat().st_size if RESULTS.exists() else 0
    payload={"status":status,"audit_identity":identity,"completed_records":len(done),
        "expected_records":len(selected),"worker_count":1,"batch_size":batch_size,
        "last_row_index":last_index,"output":RESULTS.name,
        "output_sha256":sha(RESULTS) if RESULTS.exists() else hashlib.sha256(b"").hexdigest(),
        "output_bytes":output_bytes}
    tmp=CHECKPOINT.with_suffix(".tmp")
    tmp.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    tmp.replace(CHECKPOINT)


def write_summary(selected,identity,done):
    rows=[json.loads(line) for line in RESULTS.read_text().splitlines() if line.strip()]
    if len(rows)!=len(selected) or len(done)!=len(selected):
        raise ValueError("cannot summarize an incomplete decay-accounting audit")
    statuses=Counter(row["status"] for row in rows)
    check_counts={}
    max_residuals={}
    for row in rows:
        for key,check in row["checks"].items():
            check_counts.setdefault(key,Counter())[check["status"]]+=1
            for field,value in check.items():
                if field.startswith("absolute_residual") and isinstance(value,(int,float)):
                    max_residuals[key]=max(max_residuals.get(key,0.0),float(value))
    summary={"status":"COMPLETE","method_scope":"native 2HDMC channel accounting and identity checks; not independent formulas for all non-diphoton channels",
        "audit_identity":identity,"point_count":len(selected),"record_count":len(rows),
        "status_counts":dict(statuses),"check_status_counts":{k:dict(v) for k,v in check_counts.items()},
        "maximum_absolute_residuals":max_residuals,"benchmark_physical_validation":"INCONCLUSIVE",
        "results_jsonl":RESULTS.name,"results_sha256":sha(RESULTS),"checkpoint_json":CHECKPOINT.name}
    tmp=SUMMARY.with_suffix(".tmp")
    tmp.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    tmp.replace(SUMMARY)


def _reversed_profile_record(executable):
    sys.path.insert(0,str(REPO))
    from dihiggs.app.orchestrator.aligned_lambda7 import build_aligned_lambda7_family
    from run_t4_bridge import run_point
    point=build_aligned_lambda7_family(tan_beta="10000",lambda1_ref="1.2",mh_GeV="125.13",
        mphi_GeV="100",mA_GeV="160",mHp_GeV="160",v_GeV="246",
        X_lambda6_tanbeta="1",precision=80).as_dict()
    routes=run_point(point_builder=lambda digits:point,precisions=(80,))
    match=next((row for row in routes if row.get("mode")=="analytic"),None)
    if not match or match.get("constructor_status")!="PASS" or match.get("physical_state_ordering")!="MAPPED":
        raise ValueError("lighter-phi profile anchor did not construct/map through frozen analytic bridge")
    return match


def profile_two_anchors(selected,identity):
    before=resource.getrusage(resource.RUSAGE_CHILDREN)
    start=time.perf_counter()
    pilot=selected[0]
    reversed_point=_reversed_profile_record(EXECUTABLE)
    snapshots=[]
    for label,point in (("pilot_mphi150",pilot),("reversed_phi_mphi100",reversed_point)):
        native,elapsed=run_native(point)
        audited=audit_native_payload(native,point,audit_identity=identity["fingerprint"])
        snapshots.append({"label":label,"source_point_id":point.get("point_id"),
            "mass_GeV":point.get("m_H_backend_GeV"),"phi_index":point.get("requested_phi_index_backend"),
            "status":audited["status"],"exact_failure_reason":audited["exact_failure_reason"],
            "checks":audited["checks"],"route_wall_seconds":elapsed})
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    payload={"status":"PROFILED_TWO_ANCHORS","audit_identity":identity,
        "anchors":snapshots,"worker_count":1,"wall_seconds":time.perf_counter()-start,
        "child_user_cpu_seconds":after.ru_utime-before.ru_utime,
        "child_system_cpu_seconds":after.ru_stime-before.ru_stime,
        "child_peak_rss_bytes":after.ru_maxrss*1024,
        "resource_policy":{"RLIMIT_AS":"9 GiB per child","RLIMIT_CPU":"900 s per child",
            "OMP_NUM_THREADS":1,"OPENBLAS_NUM_THREADS":1}}
    PROFILE.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    return payload


def execute_full(selected,identity,batch_size):
    done=load_existing(selected,identity)
    expected={(str(row["point_id"]),str(row["mass_convention"])) for row in selected}
    completed=len(done)
    with RESULTS.open("a",buffering=1) as output:
        pending=0
        for index,point in enumerate(selected):
            key=(str(point["point_id"]),str(point["mass_convention"]))
            if key in done:
                continue
            row=run_record(point,identity["fingerprint"])
            output.write(json.dumps(row,sort_keys=True,allow_nan=False)+"\n")
            output.flush();os.fsync(output.fileno())
            done.add(key);completed+=1;pending+=1
            if pending>=batch_size:
                write_checkpoint(identity,done,selected,batch_size,index,"RUNNING")
                pending=0
                print(f"checkpoint: {completed}/{len(expected)} decay-accounting records",flush=True)
        write_checkpoint(identity,done,selected,batch_size,len(selected)-1,
            "COMPLETE" if expected<=done else "PARTIAL")
    if expected<=done:
        write_summary(selected,identity,done)
    return {"status":"COMPLETE" if expected<=done else "PARTIAL",
        "completed":completed,"expected":len(expected),"audit_identity":identity["fingerprint"]}


def main():
    parser=argparse.ArgumentParser()
    action=parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--profile-only",action="store_true")
    action.add_argument("--run",action="store_true",help="execute/resume the full analytic/80 pilot; coordinator must authorize this phase")
    parser.add_argument("--build",action="store_true",help="build the campaign-only executable against the existing isolated library")
    parser.add_argument("--batch-size",type=int,default=50)
    args=parser.parse_args()
    if args.batch_size<1: raise SystemExit("--batch-size must be positive")
    OUT.mkdir(parents=True,exist_ok=True)
    frozen_source_spec={"frozen_scan_jsonl_sha256":sha(frozen_scan.OUTPUT),
        "frozen_checkpoint_sha256":sha(frozen_scan.CHECKPOINT),
        "audit_cpp_sha256":sha(HERE/"DecayAccountingAudit.cpp"),
        "runner_sha256":sha(HERE/"run_decay_accounting_audit.py"),
        "isolated_library_sha256":sha(CAMPAIGN/"build_analytic/2hdmc/lib/lib2HDMC.a"),
        "method_version":METHOD_VERSION}
    lock_id=canonical_sha(frozen_source_spec)
    lock_path=OUT/"audit.lock"
    with campaign_scan_lock(lock_path,lock_id):
        if args.build:
            compile_command=build_executable()
        else:
            compile_command=None
        if not EXECUTABLE.is_file(): raise SystemExit("audit executable missing; rerun with --build")
        rows,source_identity,selected=load_frozen_scan()
        identity=scan_identity(source_identity,selected,executable_identity())
        identity_path=OUT/"identity.json"
        if identity_path.exists() and json.loads(identity_path.read_text())!=identity:
            raise ValueError("decay-accounting identity.json does not match current source/grid/executable identity")
        identity_path.write_text(json.dumps(identity,indent=2,sort_keys=True)+"\n")
        if args.profile_only:
            profile=profile_two_anchors(selected,identity)
            if compile_command is not None: profile["compile_command"]=compile_command
            print(json.dumps(profile,indent=2,sort_keys=True))
            return
        print(json.dumps(execute_full(selected,identity,args.batch_size),indent=2,sort_keys=True))


if __name__=="__main__": main()
