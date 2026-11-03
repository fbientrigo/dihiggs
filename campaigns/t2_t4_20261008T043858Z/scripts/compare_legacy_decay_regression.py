#!/usr/bin/env python3
"""Compare historical evaluator output with the same source linked to the T4 build."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import subprocess
import tempfile
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
CAMPAIGN=HERE.parent
REPO=CAMPAIGN.parents[1]
sys.path.insert(0,str(REPO))
OLD=REPO/"dihiggs/app/DihiggsPointV2Evaluator"
HIST=CAMPAIGN/"verification/historical"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_point(base_command: list[str], name: str, new_executable: Path, temp: Path):
    old_cmd=list(base_command)
    new_cmd=list(base_command)
    output_index=old_cmd.index("--output")+1
    old_cmd[output_index]=str(temp/f"{name}.old.csv")
    new_cmd[0]=str(new_executable)
    new_cmd[output_index]=str(temp/f"{name}.new.csv")
    env=os.environ|{"OMP_NUM_THREADS":"1","OPENBLAS_NUM_THREADS":"1"}
    subprocess.run(old_cmd,check=True,capture_output=True,text=True,env=env,timeout=60)
    subprocess.run(new_cmd,check=True,capture_output=True,text=True,env=env,timeout=60)
    with Path(old_cmd[output_index]).open(newline="") as f: old=list(csv.DictReader(f))
    with Path(new_cmd[output_index]).open(newline="") as f: new=list(csv.DictReader(f))
    if len(old)!=1 or len(new)!=1 or old[0].keys()!=new[0].keys():
        raise RuntimeError(f"{name}: output row count or CSV schema changed")
    old,new=old[0],new[0]
    changes={}; mismatches=[]
    for key in old:
        a,b=old[key],new[key]
        try:
            af,bf=float(a),float(b)
        except ValueError:
            if a!=b: mismatches.append({"field":key,"old":a,"new":b})
            continue
        if math.isnan(af) and math.isnan(bf): continue
        difference=abs(af-bf)
        scale=max(abs(af),abs(bf))
        tolerance=max(1e-30,1e-12*scale)
        changes[key]={"absolute_difference":difference,"tolerance":tolerance}
        if not math.isfinite(difference) or difference>tolerance:
            mismatches.append({"field":key,"old":a,"new":b,"absolute_difference":difference,"tolerance":tolerance})
    return {"status":"PASS" if not mismatches else "FAIL","rows":1,
        "old_command":old_cmd,"relinked_command":new_cmd,
        "max_absolute_difference":max((x["absolute_difference"] for x in changes.values()),default=0.),
        "mismatch_count":len(mismatches),"mismatches":mismatches,
        "csv_schema_sha256":hashlib.sha256(",".join(old.keys()).encode()).hexdigest()}


def active_anchor(new_executable: Path, temp: Path):
    from dihiggs.app.orchestrator.aligned_lambda7 import build_aligned_lambda7_family
    point=build_aligned_lambda7_family(tan_beta="100",lambda1_ref="1",mh_GeV="125.20",
        mphi_GeV="150",mA_GeV="200",mHp_GeV="200",v_GeV="246",
        X_lambda6_tanbeta="0",precision=80)
    tb=100.
    m12=float(point.values["m12_sq_GeV2"])
    M2=m12/(tb/(1.+tb*tb))
    cmd=[str(OLD),"--campaign-id","t2_t4_phase2_regression","--run-id","active12520",
        "--mh","125.20","--mH-min","150","--mH-max","150","--n-mH","1",
        "--mA","200","--mHp","200","--yukawa-type","1","--sin-ba","1",
        "--tan-beta","100","--M2-min",format(M2,".17g"),"--M2-max",format(M2,".17g"),
        "--n-M2","1","--lambda6","0","--lambda7","0","--output","unused"]
    result=run_point(cmd,"active_mh12520_fixed_X0_t100",new_executable,temp)
    result["scientific_inputs"]={"mh_GeV":"125.20","mH_GeV":"150","mA_GeV":"200",
        "mHp_GeV":"200","tan_beta":"100","lambda1_ref":"1","trajectory":"fixed_X_0",
        "evaluated_X":"0","M2_GeV2":format(M2,".17g"),"lambda6":"0","lambda7":"0"}
    return result


def main():
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument("--relinked-executable",type=Path,required=True)
    p.add_argument("--output",type=Path,default=CAMPAIGN/"legacy_decay_regression.json")
    args=p.parse_args()
    temp=Path(tempfile.mkdtemp(prefix="t4-default-decay-"))
    anchors={}
    for name in ("historical450","historical500"):
        historical=json.loads((HIST/f"{name}.json").read_text())
        anchors[name]=run_point(historical["command"],name,args.relinked_executable,temp)
    anchors["active_mh12520_t100"] = active_anchor(args.relinked_executable,temp)
    result={"status":"PASS" if all(x["status"]=="PASS" for x in anchors.values()) else "FAIL",
        "scope":"default DihiggsPointV2Evaluator source linked once to the pre-phase2 executable and once to the rebuilt campaign 2HDMC library",
        "old_executable_sha256":sha(OLD),"relinked_executable_sha256":sha(args.relinked_executable),
        "source_sha256":{str(path.relative_to(REPO)):sha(path) for path in (
            REPO/"dihiggs/src/DihiggsPointV2Evaluator.cpp",REPO/"dihiggs/src/ReplaySafeOutput.cpp",
            REPO/"2hdmc/src/DecayTable.cpp",REPO/"2hdmc/src/DecayTable.h",
            REPO/"2hdmc/src/THDM.cpp",REPO/"2hdmc/src/THDM.h")},
        "comparison_tolerance":"absolute difference <= max(1e-30, 1e-12 * max(|old|, |new|)) for finite numeric CSV cells; exact CSV schema and nonnumeric cells",
        "anchors":anchors,"stdout_capture":"commands used --output paths in an isolated /tmp directory"}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":result["status"],"anchors":{k:v["status"] for k,v in anchors.items()},
        "output":str(args.output)},indent=2))
    if result["status"]!="PASS": raise SystemExit(1)

if __name__=="__main__": main()
