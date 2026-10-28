import json
import subprocess
import sys
from pathlib import Path

import pytest

CAMPAIGN = Path(__file__).resolve().parents[1] / "campaigns/t2_t4_20261008T043858Z"
sys.path.insert(0, str(CAMPAIGN / "scripts"))
import run_t4_scan as scan
import canonicalize_identical_scan_duplicates as duplicate_audit


def _record(row, identity):
    record = {
        "point_id": row["point_id"],
        "mass_convention": row["mass_convention"],
        "precision_digits_constructed": scan.PRECISIONS[0],
        "mode": scan.MODES[0],
        "scan_identity": identity["fingerprint"],
        "point_input_sha256": scan.row_fingerprint(row),
        "pilot_point_metadata": scan.row_inputs(row),
        "constructor_status": "FAIL",
        "constructor_ok": False,
        "constructor_route": "legacy",
        "exact_failure_reason": "synthetic replay fixture failure",
        "inputs_decimal": {"tan_beta": row["tan_beta"]},
        "target_values_decimal": {"lambda1": "1.0"},
        "requested_mode": "legacy",
        "declared_acceptance_tolerances": scan.DECLARED_ACCEPTANCE_TOLERANCES,
    }
    record["record_payload_sha256"] = scan.canonical_sha(record)
    return record


def _load(tmp_path, record, identity, duplicate=False):
    output = tmp_path / "results.jsonl"
    rows = [{
        "point_id": "scan-fixture", "mass_convention": "fixture",
        "slice": "fixture", "trajectory": "fixture",
        "tan_beta": "2", "lambda1_ref": "1", "mh_GeV": "125.13",
        "mphi_GeV": "150", "mA_GeV": "200", "mHp_GeV": "200",
        "v_GeV": "246", "evaluated_X": "1", "lambda6": "",
    }]
    text = json.dumps(record) + "\n"
    output.write_text(text + text if duplicate else text)
    return scan.load_existing(rows, identity, output, tmp_path / "checkpoint.json")


def test_checkpoint_replay_accepts_complete_failure_payload(tmp_path):
    row = {"point_id": "scan-fixture", "mass_convention": "fixture", "slice": "fixture",
        "trajectory": "fixture", "tan_beta": "2", "lambda1_ref": "1", "mh_GeV": "125.13",
        "mphi_GeV": "150", "mA_GeV": "200", "mHp_GeV": "200", "v_GeV": "246",
        "evaluated_X": "1", "lambda6": ""}
    identity = {"fingerprint": "fixture-identity"}
    assert len(_load(tmp_path, _record(row, identity), identity)) == 1


def test_checkpoint_replay_rejects_mutated_scientific_payload(tmp_path):
    row = {"point_id": "scan-fixture", "mass_convention": "fixture", "slice": "fixture",
        "trajectory": "fixture", "tan_beta": "2", "lambda1_ref": "1", "mh_GeV": "125.13",
        "mphi_GeV": "150", "mA_GeV": "200", "mHp_GeV": "200", "v_GeV": "246",
        "evaluated_X": "1", "lambda6": ""}
    identity = {"fingerprint": "fixture-identity"}
    record = _record(row, identity)
    record["exact_failure_reason"] = "mutated after checkpoint"
    with pytest.raises(ValueError, match="payload checksum"):
        _load(tmp_path, record, identity)


def test_checkpoint_replay_rejects_metadata_only_and_duplicates(tmp_path):
    row = {"point_id": "scan-fixture", "mass_convention": "fixture", "slice": "fixture",
        "trajectory": "fixture", "tan_beta": "2", "lambda1_ref": "1", "mh_GeV": "125.13",
        "mphi_GeV": "150", "mA_GeV": "200", "mHp_GeV": "200", "v_GeV": "246",
        "evaluated_X": "1", "lambda6": ""}
    identity = {"fingerprint": "fixture-identity"}
    record = _record(row, identity)
    record["inputs_decimal"] = None
    record["record_payload_sha256"] = scan.canonical_sha(
        {key: value for key, value in record.items() if key != "record_payload_sha256"})
    with pytest.raises(ValueError, match="malformed result payload"):
        _load(tmp_path, record, identity)
    record = _record(row, identity)
    with pytest.raises(ValueError, match="duplicate scan record"):
        _load(tmp_path, record, identity, duplicate=True)


def test_scan_lock_blocks_a_concurrent_writer_and_records_owner(tmp_path):
    lock_path=tmp_path/"scan.lock"
    command=("import sys\nfrom pathlib import Path\n"
        "sys.path.insert(0,sys.argv[1])\nfrom scan_lock import campaign_scan_lock\n"
        "try:\n with campaign_scan_lock(Path(sys.argv[2]),'second-identity'):\n  print('UNEXPECTED_ACQUIRE')\n"
        "except RuntimeError as exc:\n print(exc)\n sys.exit(3)\n")
    with scan.campaign_scan_lock(lock_path,"first-identity"):
        result=subprocess.run([sys.executable,"-c",command,str(CAMPAIGN/"scripts"),str(lock_path)],
            capture_output=True,text=True,timeout=10)
        assert result.returncode==3
        assert "pid" in result.stdout and "first-identity" in result.stdout


def _audit_record(point_id, identity="fixture-identity"):
    record={"point_id":point_id,"mass_convention":"fixture","precision_digits_constructed":50,
        "mode":"analytic","scan_identity":identity,"constructor_status":"PASS",
        "calculator_elapsed_seconds":0.1}
    record["record_payload_sha256"]=duplicate_audit.canonical_sha(record)
    return record


def _write_audit_input(path, records):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text("".join(json.dumps(record)+"\n" for record in records))


def test_duplicate_audit_writes_only_snapshot_and_verifies_payloads(tmp_path):
    records=[_audit_record("p1"),_audit_record("p1"),_audit_record("p2")]
    source=tmp_path/"input.jsonl"
    output=CAMPAIGN/"snapshots"/f"test-dedup-{tmp_path.name}.jsonl"
    _write_audit_input(source,records)
    try:
        result=duplicate_audit.canonicalize(source,output,expected_unique=2)
        assert result["input_records"]==3 and result["unique_records"]==2
        saved=[json.loads(line) for line in output.read_text().splitlines()]
        assert len(saved)==2
        assert all(record["record_payload_sha256"]==duplicate_audit.canonical_sha(
            {key:value for key,value in record.items() if key!="record_payload_sha256"}) for record in saved)
    finally:
        output.unlink(missing_ok=True)


def test_duplicate_audit_rejects_canonical_output(tmp_path):
    with pytest.raises(ValueError,match="cannot be the canonical scan path"):
        duplicate_audit.canonicalize(tmp_path/"missing.jsonl",duplicate_audit.CANONICAL_OUTPUT,
            expected_unique=0)


def test_duplicate_audit_rejects_invalid_input_payload_checksum(tmp_path):
    record=_audit_record("p1")
    record["constructor_status"]="FAIL"
    source=tmp_path/"corrupt.jsonl"
    _write_audit_input(source,[record])
    output=CAMPAIGN/"snapshots"/f"test-corrupt-{tmp_path.name}.jsonl"
    with pytest.raises(ValueError,match="invalid scientific payload checksum"):
        duplicate_audit.canonicalize(source,output,expected_unique=1)
    assert not output.exists()


def test_duplicate_audit_rejects_mixed_identities_across_keys(tmp_path):
    source=tmp_path/"mixed.jsonl"
    _write_audit_input(source,[_audit_record("p1","identity-a"),_audit_record("p2","identity-b")])
    output=CAMPAIGN/"snapshots"/f"test-mixed-{tmp_path.name}.jsonl"
    with pytest.raises(ValueError,match="mixes scan identities"):
        duplicate_audit.canonicalize(source,output,expected_unique=2)
    assert not output.exists()
