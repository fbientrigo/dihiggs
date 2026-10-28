import hashlib
import json
import sys
from pathlib import Path

import pytest

CAMPAIGN=Path(__file__).resolve().parents[1]/"campaigns/t2_t4_20261008T043858Z"
sys.path.insert(0,str(CAMPAIGN/"scripts"))
import run_decay_accounting_audit as audit


def _frozen():
    return {"point_id":"fixture-point","mass_convention":"fixture-slice",
        "scan_identity":"frozen-identity","record_payload_sha256":"source-payload",
        "precision_digits_constructed":80,"constructor_status":"PASS",
        "physical_state_ordering":"MAPPED","requested_h_index_backend":1,
        "requested_phi_index_backend":2,"m_h_backend_GeV":125.13,
        "m_H_backend_GeV":150.0,"backend_v_GeV":246.0,"requested_v_GeV":246.0,
        "width_total_phi_GeV":1.0,
        "width_diphoton_GeV":0.25,
        "independent_diphoton_components":{"status":"PASS"}}


def _native():
    terms=[]
    weights={audit.expected_term_keys()[0]:0.6,
        ("hvv",1,0):0.25,audit.expected_term_keys()[-1]:0.15}
    families={family:0.0 for family in audit.EXPECTED_FAMILIES}
    ordered=0.0
    for family,i,j in audit.expected_term_keys():
        width=weights.get((family,i,j),0.0)
        terms.append({"family":family,"i":i,"j":j,"width_GeV":width})
        families[family]+=width
        ordered+=width
    return {"constructor_ok":True,"state_identity_ok":True,
        "h_index":1,"phi_index":2,"m_h_index_GeV":125.13,"m_phi_index_GeV":150.0,
        "backend_v_GeV":246.0,"requested_v_GeV":246.0,
        "total_width_GeV":ordered,"ordered_term_sum_GeV":ordered,
        "diphoton_width_GeV":0.25,"ctau_mm":audit.HBAR_C_GEV_MM/ordered,
        "family_totals_GeV":families,"terms":terms}


def _complete_audit_row(identity="audit-identity"):
    source=_frozen()
    native=_native()
    row=audit.audit_native_payload(native,source,audit_identity=identity)
    row["calculator_elapsed_seconds"]=0.01
    row["native_payload"]=native
    row["record_payload_sha256"]=audit.canonical_sha(row)
    return source,row


def test_native_accounting_accepts_all_78_terms_and_normalizes_br_lifetime():
    result=audit.audit_native_payload(_native(),_frozen())
    assert len(result["complete_branching_fractions"])==78
    assert result["status"]=="PASS"
    assert result["decay_status"]["scope"]=="CALCULATOR_ACCOUNTING_ONLY"
    assert result["decay_status"]["benchmark_physical_validation"]=="INCONCLUSIVE"
    assert result["checks"]["complete_branching_fraction_normalization"]["status"]=="PASS"
    assert result["checks"]["lifetime_conversion"]["status"]=="PASS"


def test_native_accounting_rejects_omitted_term_and_negative_or_nonfinite_widths():
    native=_native(); native["terms"].pop()
    assert audit.audit_native_payload(native,_frozen())["status"]=="FAIL"
    for invalid in (-0.1,float("inf")):
        native=_native(); native["terms"][0]["width_GeV"]=invalid
        assert audit.audit_native_payload(native,_frozen())["status"]=="FAIL"


def test_native_accounting_rejects_total_and_lifetime_mismatches():
    native=_native(); native["total_width_GeV"]*=1.01
    frozen=_frozen(); frozen["width_total_phi_GeV"]=native["total_width_GeV"]
    result=audit.audit_native_payload(native,frozen)
    assert result["status"]=="FAIL"
    assert result["checks"]["native_total_vs_all_terms"]["status"]=="FAIL"
    native=_native(); native["ctau_mm"]*=2
    result=audit.audit_native_payload(native,_frozen())
    assert result["checks"]["lifetime_conversion"]["status"]=="FAIL"


def test_native_accounting_keeps_missing_constructor_or_diphoton_evidence_inconclusive():
    frozen=_frozen(); frozen["independent_diphoton_components"]={"status":"INCONCLUSIVE"}
    result=audit.audit_native_payload(_native(),frozen)
    assert result["status"]=="INCONCLUSIVE"
    assert "diphoton" in result["exact_failure_reason"]


def test_native_accounting_accepts_explicit_reversed_cp_even_index_mapping():
    native=_native(); frozen=_frozen()
    native.update(h_index=2,phi_index=1,m_phi_index_GeV=100.0)
    frozen.update(requested_h_index_backend=2,requested_phi_index_backend=1,m_H_backend_GeV=100.0)
    result=audit.audit_native_payload(native,frozen)
    assert result["status"]=="PASS"
    assert result["checks"]["frozen_state_identity"]["status"]=="PASS"


def test_strict_audit_replay_accepts_and_rejects_corruption(tmp_path,monkeypatch):
    source,row=_complete_audit_row()
    output=tmp_path/"results.jsonl"; checkpoint=tmp_path/"checkpoint.json"
    monkeypatch.setattr(audit,"RESULTS",output); monkeypatch.setattr(audit,"CHECKPOINT",checkpoint)
    output.write_text(json.dumps(row,sort_keys=True)+"\n")
    identity={"fingerprint":"audit-identity"}
    checkpoint.write_text(json.dumps({"status":"COMPLETE","audit_identity":identity,
        "completed_records":1,"expected_records":1,"output_sha256":hashlib.sha256(output.read_bytes()).hexdigest(),
        "output_bytes":output.stat().st_size}))
    assert audit.load_existing([source],identity)=={("fixture-point","fixture-slice")}
    altered=dict(row); altered["native_payload"]={**row["native_payload"],"total_width_GeV":2.0}
    altered["record_payload_sha256"]=audit.canonical_sha({k:v for k,v in altered.items() if k!="record_payload_sha256"})
    output.write_text(json.dumps(altered,sort_keys=True)+"\n")
    checkpoint.write_text(json.dumps({"status":"COMPLETE","audit_identity":identity,
        "completed_records":1,"expected_records":1,"output_sha256":hashlib.sha256(output.read_bytes()).hexdigest(),
        "output_bytes":output.stat().st_size}))
    with pytest.raises(ValueError,match="derived result field"):
        audit.load_existing([source],identity)


def test_strict_audit_replay_rejects_stale_prefix_checksum(tmp_path,monkeypatch):
    source,row=_complete_audit_row()
    output=tmp_path/"results.jsonl"; checkpoint=tmp_path/"checkpoint.json"
    monkeypatch.setattr(audit,"RESULTS",output); monkeypatch.setattr(audit,"CHECKPOINT",checkpoint)
    output.write_text(json.dumps(row,sort_keys=True)+"\n")
    checkpoint.write_text(json.dumps({"status":"COMPLETE","audit_identity":{"fingerprint":"audit-identity"},
        "completed_records":1,"expected_records":1,"output_sha256":"bad-checksum",
        "output_bytes":output.stat().st_size}))
    with pytest.raises(ValueError,match="prefix checksum"):
        audit.load_existing([source],{"fingerprint":"audit-identity"})


def test_strict_audit_replay_rejects_duplicate_records(tmp_path,monkeypatch):
    source,row=_complete_audit_row()
    output=tmp_path/"results.jsonl"; checkpoint=tmp_path/"checkpoint.json"
    monkeypatch.setattr(audit,"RESULTS",output); monkeypatch.setattr(audit,"CHECKPOINT",checkpoint)
    output.write_text(json.dumps(row,sort_keys=True)+"\n"+json.dumps(row,sort_keys=True)+"\n")
    identity={"fingerprint":"audit-identity"}
    checkpoint.write_text(json.dumps({"status":"RUNNING","audit_identity":identity,
        "completed_records":1,"expected_records":1,"output_sha256":hashlib.sha256(b"").hexdigest(),
        "output_bytes":0}))
    with pytest.raises(ValueError,match="duplicate audit record"):
        audit.load_existing([source],identity)


def test_complete_replay_does_not_call_native_executable(tmp_path,monkeypatch):
    source,row=_complete_audit_row()
    output=tmp_path/"results.jsonl"; checkpoint=tmp_path/"checkpoint.json"; summary=tmp_path/"summary.json"
    monkeypatch.setattr(audit,"RESULTS",output); monkeypatch.setattr(audit,"CHECKPOINT",checkpoint)
    monkeypatch.setattr(audit,"SUMMARY",summary)
    output.write_text(json.dumps(row,sort_keys=True)+"\n")
    identity={"fingerprint":"audit-identity"}
    checkpoint.write_text(json.dumps({"status":"COMPLETE","audit_identity":identity,
        "completed_records":1,"expected_records":1,"output_sha256":hashlib.sha256(output.read_bytes()).hexdigest(),
        "output_bytes":output.stat().st_size}))
    monkeypatch.setattr(audit,"run_record",lambda *a,**k: (_ for _ in ()).throw(AssertionError("native call on replay")))
    result=audit.execute_full([source],identity,50)
    assert result["status"]=="COMPLETE" and result["completed"]==1
    assert json.loads(summary.read_text())["record_count"]==1


def test_replay_rejects_missing_output_claimed_complete(tmp_path,monkeypatch):
    source,row=_complete_audit_row()
    output=tmp_path/"missing.jsonl"; checkpoint=tmp_path/"checkpoint.json"
    monkeypatch.setattr(audit,"RESULTS",output); monkeypatch.setattr(audit,"CHECKPOINT",checkpoint)
    checkpoint.write_text(json.dumps({"status":"COMPLETE","audit_identity":{"fingerprint":"audit-identity"},
        "completed_records":1,"expected_records":1}))
    with pytest.raises(ValueError,match="output file is missing"):
        audit.load_existing([source],{"fingerprint":"audit-identity"})
