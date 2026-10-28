#!/usr/bin/env python3
"""Write a checksummed duplicate audit to a separate campaign snapshot only."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

HERE=Path(__file__).resolve().parent
CAMPAIGN=HERE.parent
DEFAULT_INPUT=CAMPAIGN/"snapshots/concurrent_writer_collision_20261008T/calculator_scan_agreement_5292_lines.jsonl"
CANONICAL_OUTPUT=CAMPAIGN/"calculator_scan_agreement.jsonl"
DEFAULT_OUTPUT=CAMPAIGN/"snapshots/concurrent_writer_collision_20261008T/calculator_scan_agreement_deduplicated.jsonl"


def canonical_sha(value):
    raw=json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=True).encode()
    return hashlib.sha256(raw).hexdigest()


def canonicalize(input_path: Path, output_path: Path, *, expected_unique: int = 3510,
                 canonical_output: Path = CANONICAL_OUTPUT):
    input_path=Path(input_path).resolve()
    output_path=Path(output_path).resolve()
    canonical_output=Path(canonical_output).resolve()
    snapshots=(CAMPAIGN/"snapshots").resolve()
    if output_path==canonical_output:
        raise ValueError("duplicate audit output cannot be the canonical scan path")
    if output_path==input_path:
        raise ValueError("duplicate audit must preserve its input and write a separate snapshot")
    try:
        output_path.relative_to(snapshots)
    except ValueError as exc:
        raise ValueError("duplicate audit output must be inside the campaign snapshots directory") from exc

    groups=defaultdict(list)
    identities=set()
    total=0
    for line_number,line in enumerate(input_path.read_text().splitlines(),1):
        if not line.strip():
            continue
        try:
            record=json.loads(line)
        except json.JSONDecodeError as exc:
            raise ValueError(f"line {line_number} is malformed JSON") from exc
        if not isinstance(record,dict):
            raise ValueError(f"line {line_number} is not a scan record object")
        checksum=record.get("record_payload_sha256")
        payload={key:value for key,value in record.items() if key!="record_payload_sha256"}
        if not isinstance(checksum,str) or checksum!=canonical_sha(payload):
            raise ValueError(f"line {line_number} has an invalid scientific payload checksum")
        identity=record.get("scan_identity")
        if not isinstance(identity,str) or not identity:
            raise ValueError(f"line {line_number} has no scan identity")
        identities.add(identity)
        key=(record.get("point_id"),record.get("mass_convention"),
             record.get("precision_digits_constructed"),record.get("mode"))
        if any(value is None for value in key):
            raise ValueError(f"line {line_number} has an incomplete scan key")
        groups[key].append(record)
        total+=1
    if len(identities)!=1:
        raise ValueError(f"input mixes scan identities ({len(identities)} identities)")

    kept=[]
    removed=0
    for key,records in groups.items():
        reference=dict(records[0])
        reference.pop("calculator_elapsed_seconds",None)
        reference.pop("record_payload_sha256",None)
        for record in records[1:]:
            comparison=dict(record)
            comparison.pop("calculator_elapsed_seconds",None)
            comparison.pop("record_payload_sha256",None)
            if comparison!=reference:
                raise ValueError(f"duplicate key {key} has differing scientific payload")
        selected=dict(records[0])
        selected["record_payload_sha256"]=canonical_sha(
            {key:value for key,value in selected.items() if key!="record_payload_sha256"})
        kept.append(selected)
        removed+=len(records)-1
    if len(kept)!=expected_unique:
        raise ValueError(f"expected {expected_unique} unique route/precision keys; found {len(kept)}")
    output_path.parent.mkdir(parents=True,exist_ok=True)
    temp=output_path.with_suffix(output_path.suffix+".tmp")
    try:
        with temp.open("w") as stream:
            for record in kept:
                stream.write(json.dumps(record,sort_keys=True)+"\n")
        temp.replace(output_path)
    finally:
        if temp.exists():
            temp.unlink()
    return {"status":"CANONICALIZED_IDENTICAL_DUPLICATES","input_records":total,
        "unique_records":len(kept),"identical_duplicate_rows_removed":removed,
        "identity":next(iter(identities)),"output":str(output_path)}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,default=DEFAULT_INPUT)
    parser.add_argument("--output",type=Path,default=DEFAULT_OUTPUT)
    args=parser.parse_args()
    print(json.dumps(canonicalize(args.input,args.output),indent=2))


if __name__=="__main__":
    main()
