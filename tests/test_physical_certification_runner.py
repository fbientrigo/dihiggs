import json
import subprocess
import sys
from pathlib import Path

import pytest

CAMPAIGN = Path(__file__).resolve().parents[1] / "campaigns/t2_t4_20261008T043858Z"
sys.path.insert(0, str(CAMPAIGN / "scripts"))
import run_physical_certification as runner


@pytest.fixture(scope="module")
def complete_point():
    rows = runner._source_rows()
    t4_identity, t4_records = runner._frozen_t4(rows)
    sm_probe = runner._run_sm_probe(rows)
    identity = runner._source_spec(rows, t4_identity, sm_probe)
    row = rows[0]
    key = row["point_id"] + "|" + row["mass_convention"]
    record = runner.certify_point(row,t4_records[key],identity)
    return rows, t4_records, row, key, identity, record


def _install_replay_files(tmp_path, monkeypatch, complete_point, record=None):
    rows, t4_records, row, key, identity, fixture_record = complete_point
    output, checkpoint = tmp_path / "physical.jsonl", tmp_path / "checkpoint.json"
    monkeypatch.setattr(runner, "RESULTS", output)
    monkeypatch.setattr(runner, "CHECKPOINT", checkpoint)
    current = json.loads(json.dumps(record if record is not None else fixture_record))
    current["point_payload_sha256"] = runner.canonical_sha(
        {k: v for k, v in current.items() if k != "point_payload_sha256"})
    output.write_text(json.dumps(current) + "\n")
    checkpoint.write_text(json.dumps({"schema": runner.SCHEMA, "status": "BATCH_COMPLETE",
        "scan_identity": identity["fingerprint"], "identity": identity,
        "completed_point_keys": [key], "completed_points": 1,
        "expected_points": len(rows), "start_index": 0, "limit_points": 1}) + "\n")
    return rows, t4_records, row, key, identity


def _rewrite_record(path, mutate):
    record = json.loads(path.read_text())
    mutate(record)
    record["point_payload_sha256"] = runner.canonical_sha(
        {k: v for k, v in record.items() if k != "point_payload_sha256"})
    path.write_text(json.dumps(record) + "\n")


def _mutate_rounded_screen(record):
    for path in (record["constraints"]["rounded.quartic_screen_4pi"],
                 record["rounded_backend"]["constraints"]["quartic_screen_4pi"]):
        path["evidence"]["all_lambdas"]["lambda2"] = "1"


def test_exact_xstar_uses_trajectory_definition_not_rounded_evaluated_x():
    row = next(r for r in runner._source_rows() if r["trajectory"] == "Xstar")
    _potential, coefficients = runner.exact_target_model(row)
    assert coefficients["Z7"] == "0"
    assert coefficients["Xstar_exact"] == coefficients["X_lambda6_tanbeta"]
    assert coefficients["evaluated_X_decimal80"] == row["evaluated_X"]
    assert coefficients["X_lambda6_tanbeta"] != str(runner.q(row["evaluated_X"]))


def test_source_grid_join_and_sm_probe_are_frozen_separately():
    rows = runner._source_rows()
    assert len(rows) == 390
    assert {r["mass_convention"] for r in rows} == {"mh125.13_T1", "mh125.20_active"}
    identity, records = runner._frozen_t4(rows)
    assert len(records) == 390
    assert identity["selected_mode"] == "analytic"
    assert identity["selected_precision_digits"] == 80
    measured = runner._run_sm_probe(rows)
    assert measured["get_v2"] == runner.SM_GET_V2_TEXT
    assert "requested_v_squared=60516" in measured["stdout"]
    rows[0]["v_GeV"] = "250"
    with pytest.raises(ValueError, match="only covers v=246"):
        runner._run_sm_probe(rows)


def test_strict_replay_accepts_complete_real_point(complete_point, tmp_path, monkeypatch):
    rows, t4_records, _row, _key, identity, _record = complete_point
    _install_replay_files(tmp_path, monkeypatch, complete_point)
    done = runner._load_existing(rows, identity, t4_records)
    assert len(done) == 1
    record = next(iter(done.values()))
    assert record["point_class"] == "CERTIFIED_TREE_LEVEL"
    assert set(record["source_decimal80_lambda_Z"]) == {
        *(f"lambda{i}" for i in range(1, 8)), *(f"Z{i}" for i in range(1, 8))}


@pytest.mark.parametrize("mutation,reason", [
    (lambda r: r["constraints"].pop("target.BFB_general"), "required physics constraint set"),
    (lambda r: r["source_decimal80_lambda_Z"].pop("lambda7"), "source decimal80 lambda/Z"),
    (lambda r: r["target_exact"]["target_coefficients_exact"].pop("Z7"), "exact target Y/Z/lambda"),
    (lambda r: r["rounded_backend"]["potential"].update({"Z7": "0"}), "rounded exact-rational"),
    (_mutate_rounded_screen, "rounded quartic-screen evidence"),
    (lambda r: r.update({"overall_status": "FAIL", "point_class": "REJECTED"}), "overall status/classification"),
    (lambda r: r["constraints"]["target.BFB_general"].pop("method"), "required status/method"),
])
def test_strict_replay_rejects_rechecksummed_scientific_mutations(
        complete_point, tmp_path, monkeypatch, mutation, reason):
    rows, t4_records, _row, _key, identity, record = complete_point
    _install_replay_files(tmp_path, monkeypatch, complete_point, record)
    _rewrite_record(runner.RESULTS, mutation)
    with pytest.raises(ValueError, match=reason):
        runner._load_existing(rows, identity, t4_records)


def test_strict_replay_rejects_duplicate_and_join_mismatch(complete_point, tmp_path, monkeypatch):
    rows, t4_records, _row, _key, identity, _record = complete_point
    _install_replay_files(tmp_path, monkeypatch, complete_point)
    data = runner.RESULTS.read_text()
    runner.RESULTS.write_text(data + data)
    with pytest.raises(ValueError, match="duplicate"):
        runner._load_existing(rows, identity, t4_records)
    runner.RESULTS.write_text(data)
    bad_t4 = dict(t4_records)
    key = next(iter(runner._load_existing(rows, identity, t4_records)))
    bad_t4[key] = {**bad_t4[key], "record_payload_sha256": "wrong"}
    with pytest.raises(ValueError, match="calculator join"):
        runner._load_existing(rows, identity, bad_t4)


def test_checkpoint_identity_schema_and_pending_state_are_strict(complete_point, tmp_path, monkeypatch):
    rows, t4_records, _row, _key, identity, _record = complete_point
    _install_replay_files(tmp_path, monkeypatch, complete_point)
    checkpoint = json.loads(runner.CHECKPOINT.read_text())
    checkpoint["schema"] = "wrong"
    runner.CHECKPOINT.write_text(json.dumps(checkpoint))
    with pytest.raises(ValueError, match="checkpoint schema"):
        runner._load_existing(rows, identity, t4_records)


def test_campaign_lock_rejects_concurrent_identity(tmp_path):
    lock_path = tmp_path / "physical.lock"
    command = ("import sys\nfrom pathlib import Path\n"
        "sys.path.insert(0,sys.argv[1])\nfrom scan_lock import campaign_scan_lock\n"
        "try:\n with campaign_scan_lock(Path(sys.argv[2]),'second'):\n  print('UNEXPECTED')\n"
        "except RuntimeError as exc:\n print(exc)\n sys.exit(3)\n")
    with runner.campaign_scan_lock(lock_path, "first"):
        result = subprocess.run([sys.executable, "-c", command,
            str(CAMPAIGN / "scripts"), str(lock_path)], capture_output=True, text=True, timeout=10)
        assert result.returncode == 3
        assert "first" in result.stdout


def test_final_selected_cpu_interruption_is_pending_and_never_complete(tmp_path, monkeypatch):
    row = {"point_id": "last", "mass_convention": "fixture", "slice": "fixture",
        "trajectory": "fixture", "v_GeV": "246", "lambda1": "1"}
    key = "last|fixture"
    identity = {"fingerprint": "cpu-fixture"}
    calc = {"record_payload_sha256": "t4-sha", "scan_identity": "t4-identity",
        "inputs_decimal": {"lambda1": "1"}, "values_decimal": {"lambda1": "1"}}
    out = tmp_path / "out"
    monkeypatch.setattr(runner, "CAMPAIGN", tmp_path)
    monkeypatch.setattr(runner, "OUT", out)
    monkeypatch.setattr(runner, "RESULTS", out / "results.jsonl")
    monkeypatch.setattr(runner, "CHECKPOINT", out / "checkpoint.json")
    monkeypatch.setattr(runner, "LOCK", tmp_path / "lock")
    monkeypatch.setattr(runner, "_source_rows", lambda: [row])
    monkeypatch.setattr(runner, "_frozen_t4", lambda _rows: ({}, {key: calc}))
    monkeypatch.setattr(runner, "_run_sm_probe", lambda _rows: {})
    monkeypatch.setattr(runner, "_source_spec", lambda *_args: identity)
    monkeypatch.setattr(runner, "_enforce_resources", lambda: None)
    calls=[]
    def interrupt(*_args):
        calls.append("attempt")
        raise runner.BatchCpuLimit("simulated final-row CPU interruption")
    monkeypatch.setattr(runner, "certify_point", interrupt)
    result = runner.run(type("Args", (), {"start_index": 0, "limit_points": 1})())
    checkpoint = json.loads((out / "checkpoint.json").read_text())
    assert result["status"] == "RESOURCE_LIMITED"
    assert result["completed_total"] == 0
    assert (out / "results.jsonl").read_text()==""
    assert checkpoint["status"] == "RESOURCE_LIMITED"
    assert checkpoint["pending_point"]["point_key"] == key
    assert checkpoint["completed_points"] == 0
    assert runner._load_existing([row],identity,{key:calc})=={}
    retry = runner.run(type("Args", (), {"start_index": 0, "limit_points": 1})())
    assert retry["status"]=="RESOURCE_LIMITED"
    assert len(calls)==2


def test_constraint_aggregation_fails_closed():
    assert runner._aggregate({"a": {"status": "PASS"}}) == "PASS"
    assert runner._aggregate({"a": {"status": "NOT_TESTED"}}) == "INCONCLUSIVE"
    assert runner._aggregate({"a": {"status": "INCONCLUSIVE"},
                              "b": {"status": "FAIL"}}) == "FAIL"
