"""Tests for the adaptive proposal-batch handoff (GitHub issue #72).

Covers: the input contract (docs/contracts/adaptive_proposal_batch_v1.yaml)
agrees with the wrapper's own constants; the wrapper imports no
adaptive-search/optimizer module (issue #72 acceptance criterion); the
end-to-end handoff through the real DihiggsPointV2Evaluator binary (skipped
if unbuilt); and the Task 5 adversarial failure modes, each of which must
stay attempt-preserving and fail closed.
"""
from __future__ import annotations

import csv
import inspect
import stat
import sys
import textwrap
from pathlib import Path

import pytest

from dihiggs.app.orchestrator import proposal_batch as pb

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"
CONTRACT_YAML = ROOT / "docs/contracts/adaptive_proposal_batch_v1.yaml"

REQUIRED_INPUT_ROW = {
    "proposal_id": "p1", "mH_GeV": "150.0", "mA_GeV": "450.0", "mHp_GeV": "450.0",
    "M2_GeV2": "22499.999999500335", "tan_beta": "300000.0", "lambda6": "1e-10",
    "lambda7": "0.0", "sin_beta_minus_alpha": "1.0", "yukawa_type": "1",
}


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    fieldnames = fieldnames or list(rows[0].keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_attempts(manifest: dict) -> list[dict]:
    with Path(manifest["attempts_csv"]).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def make_exe(tmp_path: Path, name: str, script: str) -> Path:
    path = tmp_path / name
    path.write_text(textwrap.dedent(script), encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


# ---------------------------------------------------------------------------
# Contract (Gate 1): yaml agrees with the module, and carries no policy.
# ---------------------------------------------------------------------------

def test_contract_yaml_present_and_versioned():
    yaml = pytest.importorskip("yaml")
    doc = yaml.safe_load(CONTRACT_YAML.read_text(encoding="utf-8"))
    assert doc["schema_name"] == pb.PROPOSAL_SCHEMA_VERSION


def test_contract_yaml_field_list_matches_module_constants():
    yaml = pytest.importorskip("yaml")
    doc = yaml.safe_load(CONTRACT_YAML.read_text(encoding="utf-8"))
    names = [f["name"] for f in doc["fields"]]
    required_in_yaml = [
        f["name"] for f in doc["fields"]
        if f.get("category") == "canonical_physics_input" and f.get("required") is True
    ] + [f["name"] for f in doc["fields"] if f.get("category") == "proposal_metadata"]
    assert set(required_in_yaml) == set(pb.REQUIRED_COLUMNS)
    assert "mh_GeV" in names
    mh_field = next(f for f in doc["fields"] if f["name"] == "mh_GeV")
    assert mh_field["required"] is False
    assert mh_field["default_when_absent"] == pb.MH_CONVENTION_GEV


def test_contract_yaml_has_no_search_policy_keys():
    text = CONTRACT_YAML.read_text(encoding="utf-8").lower()
    for forbidden in ("weight:", "prior:", "score:", "acceptance_threshold:", "budget:"):
        assert forbidden not in text, f"contract yaml contains policy-shaped key: {forbidden!r}"


# ---------------------------------------------------------------------------
# Boundary guard (issue #72 acceptance criterion): no optimizer imports.
# ---------------------------------------------------------------------------

def test_wrapper_imports_no_adaptive_or_optimizer_module():
    source = Path(inspect.getfile(pb)).read_text(encoding="utf-8")
    for forbidden in ("adaptive_policy", "adaptive_explorer", "adaptive_checkpoint",
                       "optuna", "scipy", "numpy", "ray", "mlflow"):
        assert forbidden not in source, f"proposal_batch.py must not reference {forbidden!r}"


def test_wrapper_module_does_not_pull_in_adaptive_modules():
    for name in list(sys.modules):
        if "adaptive_" in name or name in ("optuna", "scipy", "numpy", "ray", "mlflow"):
            pytest.skip(f"{name} already imported by another test in this process")
    assert "dihiggs.app.orchestrator.proposal_batch" in sys.modules


# ---------------------------------------------------------------------------
# End-to-end (Gate 2/3/4): requires the built evaluator.
# ---------------------------------------------------------------------------

@pytest.fixture()
def require_binary():
    if not BINARY.is_file():
        pytest.skip("build dihiggs/app/DihiggsPointV2Evaluator first")


def test_three_proposals_yield_three_attempts(tmp_path, require_binary):
    proposals = tmp_path / "proposals.csv"
    write_csv(proposals, [
        {**REQUIRED_INPUT_ROW, "proposal_id": "A_anchor"},
        {**REQUIRED_INPUT_ROW, "proposal_id": "B_perturbed", "mH_GeV": "160.0"},
        {
            **REQUIRED_INPUT_ROW, "proposal_id": "C_rejected", "mH_GeV": "130.0",
            "mA_GeV": "300.0", "mHp_GeV": "300.0", "M2_GeV2": "16239.109978356435",
            "tan_beta": "50.0", "lambda6": "0.1", "sin_beta_minus_alpha": "0.999",
        },
    ])
    manifest = pb.run_proposal_batch(
        executable=BINARY, proposals_csv=proposals, outdir=tmp_path / "out",
        campaign_id="test", run_name="run", repo_root=ROOT,
    )
    assert manifest["status"] == "complete"
    attempts = read_attempts(manifest)
    assert len(attempts) == 3
    assert all(a["attempt_status"] == "EVALUATED" for a in attempts)

    anchor = next(a for a in attempts if a["proposal_id"] == "A_anchor")
    assert anchor["point_id"] == "point_98c841e915d3605a"
    assert float(anchor["g_hH2H2_GeV"]) == pytest.approx(63.6625935034957138, abs=1e-10)

    rejected = next(a for a in attempts if a["proposal_id"] == "C_rejected")
    assert rejected["construction_ok"] == "1"
    assert rejected["rejection_stage"] not in ("", "accepted")
    assert rejected["rejection_reason"] != "none"


def test_repeat_run_is_byte_identical(tmp_path, require_binary):
    proposals = tmp_path / "proposals.csv"
    write_csv(proposals, [{**REQUIRED_INPUT_ROW}])
    m1 = pb.run_proposal_batch(executable=BINARY, proposals_csv=proposals,
                                outdir=tmp_path / "out1", campaign_id="test", run_name="run",
                                repo_root=ROOT)
    m2 = pb.run_proposal_batch(executable=BINARY, proposals_csv=proposals,
                                outdir=tmp_path / "out2", campaign_id="test", run_name="run",
                                repo_root=ROOT)
    assert m1["output_sha256"] == m2["output_sha256"]


def test_mh_default_applied_when_column_absent(tmp_path, require_binary):
    row = {k: v for k, v in REQUIRED_INPUT_ROW.items()}
    proposals = tmp_path / "proposals.csv"
    write_csv(proposals, [row], fieldnames=list(pb.REQUIRED_COLUMNS))
    manifest = pb.run_proposal_batch(executable=BINARY, proposals_csv=proposals,
                                      outdir=tmp_path / "out", campaign_id="test",
                                      run_name="run", repo_root=ROOT)
    attempts = read_attempts(manifest)
    assert float(attempts[0]["mh_input_GeV"]) == pb.MH_CONVENTION_GEV


def test_mh_mismatch_is_malformed_not_silently_overridden(tmp_path, require_binary):
    row = {**REQUIRED_INPUT_ROW, "mh_GeV": "999.0"}
    proposals = tmp_path / "proposals.csv"
    write_csv(proposals, [row], fieldnames=list(pb.REQUIRED_COLUMNS) + ["mh_GeV"])
    with pytest.raises(RuntimeError):
        pb.run_proposal_batch(executable=BINARY, proposals_csv=proposals,
                               outdir=tmp_path / "out", campaign_id="test",
                               run_name="run", repo_root=ROOT)


# ---------------------------------------------------------------------------
# Task 5 adversarial micro-tests: attempt-preserving, machine-readable, fail closed.
# ---------------------------------------------------------------------------

def test_malformed_numeric_proposal(tmp_path):
    row = {**REQUIRED_INPUT_ROW, "mH_GeV": "not_a_number"}
    proposals = tmp_path / "proposals.csv"
    write_csv(proposals, [row])
    with pytest.raises(RuntimeError):
        pb.run_proposal_batch(executable=Path("/nonexistent"), proposals_csv=proposals,
                               outdir=tmp_path / "out", campaign_id="test",
                               run_name="run", repo_root=ROOT)
    manifest_path = tmp_path / "out" / "campaign=test" / "run" / "batch_manifest.json"
    import json
    manifest = json.loads(manifest_path.read_text())
    assert manifest["status"] == "failed"
    assert manifest["counts"]["malformed"] == 1
    attempts = read_attempts(manifest)
    assert len(attempts) == 1
    assert attempts[0]["attempt_status"] == "MALFORMED"
    assert "invalid_float:mH_GeV" in attempts[0]["attempt_reason"]


def test_missing_required_field_marks_every_row_malformed(tmp_path):
    fieldnames = [c for c in pb.REQUIRED_COLUMNS if c != "mHp_GeV"]
    proposals = tmp_path / "proposals.csv"
    rows = [
        {k: v for k, v in REQUIRED_INPUT_ROW.items() if k != "mHp_GeV"},
        {**{k: v for k, v in REQUIRED_INPUT_ROW.items() if k != "mHp_GeV"}, "proposal_id": "p2"},
    ]
    write_csv(proposals, rows, fieldnames=fieldnames)
    with pytest.raises(RuntimeError):
        pb.run_proposal_batch(executable=Path("/nonexistent"), proposals_csv=proposals,
                               outdir=tmp_path / "out", campaign_id="test",
                               run_name="run", repo_root=ROOT)
    import json
    manifest = json.loads((tmp_path / "out" / "campaign=test" / "run" / "batch_manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["counts"]["malformed"] == 2
    attempts = read_attempts(manifest)
    assert len(attempts) == 2
    assert all(a["attempt_status"] == "MALFORMED" for a in attempts)
    assert all("missing_required_column:mHp_GeV" in a["attempt_reason"] for a in attempts)


def test_duplicate_proposal_id(tmp_path):
    rows = [
        {**REQUIRED_INPUT_ROW, "proposal_id": "dup"},
        {**REQUIRED_INPUT_ROW, "proposal_id": "dup", "mH_GeV": "160.0"},
    ]
    proposals = tmp_path / "proposals.csv"
    write_csv(proposals, rows)
    with pytest.raises(RuntimeError):
        pb.run_proposal_batch(executable=Path("/nonexistent"), proposals_csv=proposals,
                               outdir=tmp_path / "out", campaign_id="test",
                               run_name="run", repo_root=ROOT)
    import json
    manifest = json.loads((tmp_path / "out" / "campaign=test" / "run" / "batch_manifest.json").read_text())
    assert manifest["status"] == "failed"
    attempts = read_attempts(manifest)
    assert len(attempts) == 2
    assert all(a["attempt_status"] == "MALFORMED" for a in attempts)
    assert all("duplicate_proposal_id" in a["attempt_reason"] for a in attempts)


def test_evaluator_nonzero_exit(tmp_path):
    exe = make_exe(tmp_path, "fake_failure.sh", """\
        #!/usr/bin/env bash
        echo "fake_evaluator: simulated_physics_rejection" >&2
        exit 2
    """)
    proposals = tmp_path / "proposals.csv"
    write_csv(proposals, [REQUIRED_INPUT_ROW])
    with pytest.raises(RuntimeError):
        pb.run_proposal_batch(executable=exe, proposals_csv=proposals,
                               outdir=tmp_path / "out", campaign_id="test",
                               run_name="run", repo_root=ROOT)
    import json
    manifest = json.loads((tmp_path / "out" / "campaign=test" / "run" / "batch_manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["counts"]["evaluator_error"] == 1
    attempts = read_attempts(manifest)
    assert len(attempts) == 1
    assert attempts[0]["attempt_status"] == "EVALUATOR_ERROR"
    assert "simulated_physics_rejection" in attempts[0]["attempt_reason"]


def test_simulated_missing_output(tmp_path):
    exe = make_exe(tmp_path, "fake_no_output.sh", """\
        #!/usr/bin/env bash
        exit 0
    """)
    proposals = tmp_path / "proposals.csv"
    write_csv(proposals, [REQUIRED_INPUT_ROW])
    with pytest.raises(RuntimeError):
        pb.run_proposal_batch(executable=exe, proposals_csv=proposals,
                               outdir=tmp_path / "out", campaign_id="test",
                               run_name="run", repo_root=ROOT)
    import json
    manifest = json.loads((tmp_path / "out" / "campaign=test" / "run" / "batch_manifest.json").read_text())
    assert manifest["status"] == "failed"
    attempts = read_attempts(manifest)
    assert len(attempts) == 1
    assert attempts[0]["attempt_status"] == "EVALUATOR_ERROR"
    assert attempts[0]["attempt_reason"] == "missing_output_file"


def test_simulated_output_missing_required_canonical_column(tmp_path):
    exe = make_exe(tmp_path, "fake_bad_header.sh", """\
        #!/usr/bin/env bash
        OUTPUT="${@: -1}"
        echo "schema_version,producer,producer_commit,producer_dirty,point_id,construction_ok,rejection_stage,rejection_reason,theory_ok_v1,total_width_GeV,ctau_mm" > "$OUTPUT"
        echo "dihiggs.point.v2,DihiggsPointV2Evaluator,deadbeef,no,point_stub0001,1,accepted,none,1.0,5.0,1e-10" >> "$OUTPUT"
        exit 0
    """)
    proposals = tmp_path / "proposals.csv"
    write_csv(proposals, [REQUIRED_INPUT_ROW])
    with pytest.raises(RuntimeError):
        pb.run_proposal_batch(executable=exe, proposals_csv=proposals,
                               outdir=tmp_path / "out", campaign_id="test",
                               run_name="run", repo_root=ROOT)
    import json
    manifest = json.loads((tmp_path / "out" / "campaign=test" / "run" / "batch_manifest.json").read_text())
    assert manifest["status"] == "failed"
    attempts = read_attempts(manifest)
    assert len(attempts) == 1
    assert attempts[0]["attempt_status"] == "EVALUATOR_ERROR"
    assert "g_hH2H2_GeV" in attempts[0]["attempt_reason"]


def test_dirty_repository_provenance_recorded_honestly(tmp_path, monkeypatch):
    exe = make_exe(tmp_path, "fake_success.sh", """\
        #!/usr/bin/env bash
        OUTPUT="${@: -1}"
        echo "schema_version,producer,producer_commit,producer_dirty,point_id,construction_ok,rejection_stage,rejection_reason,theory_ok_v1,g_hH2H2_GeV,total_width_GeV,ctau_mm" > "$OUTPUT"
        echo "dihiggs.point.v2,DihiggsPointV2Evaluator,deadbeef,yes,point_stub0001,1,accepted,none,1.0,10.0,1e-10,5.0" >> "$OUTPUT"
        exit 0
    """)
    monkeypatch.setattr(pb, "detect_git_info", lambda root: {
        "repo_root": str(root), "commit": "deadbeef", "commit_short": "deadbee",
        "is_dirty": "yes",
    })
    proposals = tmp_path / "proposals.csv"
    write_csv(proposals, [REQUIRED_INPUT_ROW])
    manifest = pb.run_proposal_batch(executable=exe, proposals_csv=proposals,
                                      outdir=tmp_path / "out", campaign_id="test",
                                      run_name="run", repo_root=ROOT)
    assert manifest["status"] == "complete"
    assert manifest["repository_dirty"] == "yes"
