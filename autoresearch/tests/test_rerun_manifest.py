from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from autoresearch.harness.rerun_manifest import (
    build_rerun_manifest,
    load_campaign_state,
    write_rerun_manifest_json,
    write_rerun_manifest_markdown,
)


def _campaign_state_fixture() -> dict[str, object]:
    return {
        "schema_version": "1.0",
        "campaign": "phase5_0d_transition",
        "proposals": [
            {
                "proposal_id": "p1",
                "operator": "autoresearch",
                "bounds": {"lambda6": [0.0006, 0.0006], "lambda1": [4.0, 5.0], "mH": [250.0, 260.0], "tan_beta": [10000.0, 40000.0]},
                "scan_values": {"tan_beta": [30000.0]},
                "fixed": {"sin_ba": 1.0, "lambda7": 0.0, "mA": 300.0, "mHp": 300.0, "yukawa_type": 1},
                "resolution": {"mH": 2, "lambda1": 2},
            }
        ],
        "runs": [
            {
                "run_dir": "/tmp/run1",
                "proposal_id": "p1",
                "health_status": "failed",
                "dry_run": False,
                "recommended_action": "needs_rerun_or_inspection",
                "csv_count": 1,
            },
            {
                "run_dir": "/tmp/run2",
                "proposal_id": "p2",
                "health_status": "ok",
                "dry_run": False,
                "recommended_action": "ready_for_review",
                "csv_count": 1,
            },
            {
                "run_dir": "/tmp/run3",
                "proposal_id": "p3",
                "health_status": "failed",
                "dry_run": True,
                "recommended_action": "needs_rerun_or_inspection",
                "csv_count": 1,
            },
        ],
        "slices": [
            {
                "source_csv": "/tmp/run1/tb_1000/scan_tb_1000.csv",
                "run_dir": "/tmp/run1",
                "proposal_id": "p1",
                "scoreable": False,
                "is_empty": False,
                "missing_columns": ["m_phi"],
            },
            {
                "source_csv": "/tmp/run2/tb_1000/scan_tb_1000.csv",
                "run_dir": "/tmp/run2",
                "proposal_id": "p2",
                "scoreable": True,
                "is_empty": False,
                "missing_columns": [],
            },
            {
                "source_csv": "/tmp/run3/tb_1000/scan_tb_1000.csv",
                "run_dir": "/tmp/run3",
                "proposal_id": "p3",
                "scoreable": False,
                "is_empty": False,
                "missing_columns": ["m_phi"],
                "dry_run": True,
            },
        ],
    }


def test_build_rerun_manifest_selection_and_reason_codes() -> None:
    state = _campaign_state_fixture()
    out = build_rerun_manifest(state, source_campaign_state="/tmp/campaign_state.json")

    assert out["schema_version"] == "1.0"
    assert out["selection_policy"]["include_actions"] == ["needs_rerun_or_inspection"]

    selected_runs = [r for r in out["run_actions"] if r["selected"]]
    assert len(selected_runs) == 1
    assert selected_runs[0]["run_dir"] == "/tmp/run1"
    assert "missing_required_columns" in selected_runs[0]["reason_codes"]
    assert selected_runs[0]["affected_slices"]

    run3 = [r for r in out["run_actions"] if r["run_dir"] == "/tmp/run3"][0]
    assert run3["selected"] is False
    assert "dry_run_rejected" in run3["reason_codes"]

    selected_slices = [s for s in out["slice_actions"] if s["selected"]]
    assert len(selected_slices) == 1
    assert selected_slices[0]["source_csv"] == "/tmp/run1/tb_1000/scan_tb_1000.csv"
    assert "missing_required_columns" in selected_slices[0]["reason_codes"]
    assert selected_slices[0]["evidence_paths"]


def test_write_json_markdown_and_required_sentence(tmp_path: Path) -> None:
    manifest = build_rerun_manifest(_campaign_state_fixture(), source_campaign_state="/tmp/campaign_state.json")
    out_json = tmp_path / "rerun_manifest.json"
    out_md = tmp_path / "rerun_manifest.md"

    write_rerun_manifest_json(manifest, out_json)
    write_rerun_manifest_markdown(manifest, out_md)

    loaded = json.loads(out_json.read_text(encoding="utf-8"))
    assert loaded["summary"]["n_runs_selected"] == 1

    md = out_md.read_text(encoding="utf-8")
    assert "Coverage and novelty are not implemented yet" in md
    assert "## Selected runs" in md
    assert "## Selected slices" in md


def test_load_campaign_state_errors_on_invalid_json(tmp_path: Path) -> None:
    p = tmp_path / "bad.json"
    p.write_text("{bad", encoding="utf-8")
    try:
        load_campaign_state(p)
    except ValueError as exc:
        assert "failed to load campaign_state" in str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_cli_generates_outputs(tmp_path: Path) -> None:
    state_path = tmp_path / "campaign_state.json"
    state_path.write_text(json.dumps(_campaign_state_fixture()), encoding="utf-8")
    out_json = tmp_path / "rerun_manifest.json"
    out_md = tmp_path / "rerun_manifest.md"

    p = subprocess.run(
        [
            sys.executable,
            "-m",
            "autoresearch.harness.rerun_manifest",
            "--campaign-state",
            str(state_path),
            "--out-json",
            str(out_json),
            "--out-md",
            str(out_md),
        ],
        capture_output=True,
        text=True,
    )
    assert p.returncode == 0
    assert out_json.exists()
    assert out_md.exists()


def test_command_contract_missing_payload_is_unknown() -> None:
    state = _campaign_state_fixture()
    state["runs"][0]["proposal_id"] = "p_missing"
    out = build_rerun_manifest(state, source_campaign_state="/tmp/campaign_state.json")
    run1 = [r for r in out["run_actions"] if r["run_dir"] == "/tmp/run1"][0]
    assert run1["command_contract_ready"] == "unknown"
    assert "proposal_payload_missing" in run1["command_contract_reason_codes"]


def test_command_contract_fields_reported_for_known_payload() -> None:
    out = build_rerun_manifest(_campaign_state_fixture(), source_campaign_state="/tmp/campaign_state.json")
    run1 = [r for r in out["run_actions"] if r["run_dir"] == "/tmp/run1"][0]
    assert run1["lambda6_singleton"] is True
    assert run1["has_scan_values_tan_beta"] is True
    assert "command_contract_ready" in run1
    assert run1["command_contract_hint"]["effective_cpp_omp_threads"] == 1
    assert "not a broad safety" in run1["command_contract_hint"]["thread_policy_note"]


def test_command_contract_reports_explicit_cpp_omp_override() -> None:
    state = _campaign_state_fixture()
    state["proposals"][0]["runtime_config"] = {"cpp_omp_threads": 2}
    out = build_rerun_manifest(state, source_campaign_state="/tmp/campaign_state.json")
    run1 = [r for r in out["run_actions"] if r["run_dir"] == "/tmp/run1"][0]
    assert run1["command_contract_hint"]["effective_cpp_omp_threads"] == 2
