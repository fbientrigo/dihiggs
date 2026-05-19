from __future__ import annotations

import csv
import json
from pathlib import Path

from autoresearch.harness.run_health import required_physics_columns
from autoresearch.harness.safe_automation_layer import apply_gates, build_commands, run_safe_campaign, triple_ok_only_analysis, validate_contract


def _mkrow(*, triple_ok: bool, br_gaga: float, width_bb: float, total_width: float) -> list[str]:
    ok = "1" if triple_ok else "0"
    return [
        "260.0", "0.001", "0.0", "30000", "4.2", ok, ok, ok,
        f"{width_bb}", "1e-14", "1e-15", "1e-13", "1e-13", "1e-13", "1e-13", f"{total_width}", f"{br_gaga}",
    ]


def _write_csv(path: Path, rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as h:
        w = csv.writer(h)
        w.writerow(required_physics_columns())
        w.writerows(rows)


def _proposal(pid: str) -> dict[str, object]:
    return {
        "proposal_id": pid,
        "operator": "autoresearch",
        "bounds": {"lambda6": [0.001, 0.001], "tan_beta": [30000, 30000], "lambda1": [4.0, 4.5], "mH": [260.0, 262.0]},
        "scan_values": {"tan_beta": [30000]},
        "fixed": {"sin_ba": 1.0, "lambda7": 0.0, "mA": 300.0, "mHp": 300.0, "yukawa_type": 1},
        "resolution": {"mH": 2, "lambda1": 2},
    }


def _contract(tmp_path: Path) -> dict[str, object]:
    outdir = tmp_path / "out"
    return {
        "enable_proposal_generation": False,
        "allow_broad_scans": False,
        "modify_physics_semantics": False,
        "modify_scoring_semantics": False,
        "modify_triple_ok_definition": False,
        "use_all_row_br_metrics_as_physics_claims": False,
        "investigate_gsl_root_cause": False,
        "require_cpp_omp_threads_1": True,
        "analysis": {"triple_ok_only": True},
        "gates": {"min_triple_ok_rate": 0.2},
        "evidence_mode": "fixture",
        "runtime": {"threads": 1, "exec_path": "/bin/true", "outdir": str(outdir), "lake_name": "lake", "dry_run": True, "force": False},
        "campaign": {"name": "safe_v1_test", "max_runs": 1, "expansion_max_runs": 0, "proposals": [_proposal("p1")]},
    }


def test_contract_validation_strict_policy() -> None:
    c = validate_contract(_contract(Path("/tmp")))
    assert c["require_cpp_omp_threads_1"] is True
    bad = _contract(Path("/tmp"))
    bad["enable_proposal_generation"] = True
    try:
        validate_contract(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("expected contract rejection")


def test_threads_1_enforced_in_commands(tmp_path: Path) -> None:
    cmds = build_commands(_contract(tmp_path))
    assert len(cmds) == 1
    cmd = cmds[0]
    assert "--threads" in cmd
    assert cmd[cmd.index("--threads") + 1] == "1"


def test_triple_ok_only_metric_computation(tmp_path: Path) -> None:
    run = tmp_path / "run1"
    _write_csv(run / "tb_30000" / "scan_tb_30000.csv", [
        _mkrow(triple_ok=True, br_gaga=1e-4, width_bb=2e-3, total_width=1e-2),
        _mkrow(triple_ok=False, br_gaga=9e-4, width_bb=9e-3, total_width=1e-2),
        _mkrow(triple_ok=True, br_gaga=2e-4, width_bb=1e-3, total_width=2e-2),
    ])
    out = triple_ok_only_analysis([run])
    assert out["n_rows"] == 3
    assert out["n_triple_ok"] == 2
    assert out["ranking_basis"] == "triple_ok_rows_only"


def test_apply_gates_missing_manifest_fails_in_production() -> None:
    campaign_state = {
        "summary": {"n_health_failed": 0},
        "runs": [{"run_manifest_path": "/tmp/nope/run_manifest.json", "orchestrator_log_path": "/tmp/nope/orchestrator.log", "task_summary_path": "/tmp/nope/task_summary.jsonl"}],
        "slices": [],
    }
    out = apply_gates(
        commands=[["python", "x.py", "--threads", "1"]],
        run_results=[],
        campaign_state=campaign_state,
        rerun_manifest={"summary": {"n_runs_selected": 0, "n_slices_selected": 0}},
        analysis={"triple_ok_rate": 1.0, "ranking_basis": "triple_ok_rows_only"},
        min_triple_ok_rate=0.2,
        production_validation=True,
        evidence_mode="production",
    )
    assert out["passed"] is False
    assert any(f["code"] == "gate_failed_due_to_missing_evidence" for f in out["failures"])


def test_apply_gates_header_only_csv_fails(tmp_path: Path) -> None:
    csv_path = tmp_path / "a.csv"
    _write_csv(csv_path, [])
    out = apply_gates(
        commands=[["python", "x.py", "--threads", "1"]],
        run_results=[],
        campaign_state={"summary": {"n_health_failed": 0}, "runs": [], "slices": [{"source_csv": str(csv_path), "row_count": 0}]},
        rerun_manifest={"summary": {"n_runs_selected": 0, "n_slices_selected": 0}},
        analysis={"triple_ok_rate": 1.0, "ranking_basis": "triple_ok_rows_only"},
        min_triple_ok_rate=0.2,
        production_validation=False,
        evidence_mode="fixture",
    )
    assert any(f["code"] == "gate_failed_due_to_header_only_csv" for f in out["failures"])


def test_stop_report_no_execute_is_non_production(tmp_path: Path) -> None:
    c = _contract(tmp_path)
    campaign_dir = Path(c["runtime"]["outdir"]) / "lake" / "safe_v1_test" / "r1" / "run"
    _write_csv(campaign_dir / "tb_30000" / "scan_tb_30000.csv", [_mkrow(triple_ok=True, br_gaga=1e-4, width_bb=1e-3, total_width=1e-2)])
    (campaign_dir / "task_summary.jsonl").write_text('{"event":"done"}\n', encoding="utf-8")
    (campaign_dir / "run_manifest.json").write_text(json.dumps({"runtime": {"omp_num_threads": 1}}), encoding="utf-8")
    (campaign_dir / "orchestrator.log").write_text("[CONF] OMP_NUM_THREADS = 1\n", encoding="utf-8")

    out = run_safe_campaign(c, workdir=tmp_path, execute=False)
    stop = out["stop_report"]
    assert stop["production_validation"] is False
    assert stop["evidence_mode"] == "fixture"
    assert stop["expansion_allowed"] is False


def test_min_triple_ok_rate_gate() -> None:
    out = apply_gates(
        commands=[["python", "x.py", "--threads", "1"]],
        run_results=[],
        campaign_state={"summary": {"n_health_failed": 0}, "runs": [], "slices": []},
        rerun_manifest={"summary": {"n_runs_selected": 0, "n_slices_selected": 0}},
        analysis={"triple_ok_rate": 0.1, "ranking_basis": "triple_ok_rows_only"},
        min_triple_ok_rate=0.2,
        production_validation=False,
        evidence_mode="fixture",
    )
    assert any(f["code"] == "gate_failed_due_to_physics_viability" for f in out["failures"])


def test_existing_artifacts_mode_can_validate_as_production(tmp_path: Path) -> None:
    c = _contract(tmp_path)
    c["evidence_mode"] = "existing_artifacts"
    campaign_dir = Path(c["runtime"]["outdir"]) / "lake" / "safe_v1_test" / "r1" / "run"
    _write_csv(campaign_dir / "tb_30000" / "scan_tb_30000.csv", [_mkrow(triple_ok=True, br_gaga=1e-4, width_bb=1e-3, total_width=1e-2)])
    (campaign_dir / "task_summary.jsonl").write_text('{"event":"done"}\n', encoding="utf-8")
    (campaign_dir / "run_manifest.json").write_text(json.dumps({"runtime": {"omp_num_threads": 1}}), encoding="utf-8")
    (campaign_dir / "orchestrator.log").write_text("[CONF] OMP_NUM_THREADS = 1\n", encoding="utf-8")

    out = run_safe_campaign(c, workdir=tmp_path, execute=False)
    stop = out["stop_report"]
    assert stop["evidence_mode"] == "existing_artifacts"
    assert stop["gates"]["passed"] is True
    assert stop["production_validation"] is True


def test_run_name_scoping_excludes_unrelated_artifacts(tmp_path: Path) -> None:
    c = _contract(tmp_path)
    c["evidence_mode"] = "existing_artifacts"
    c["run_filter"] = {"run_name": "r_good"}

    base = Path(c["runtime"]["outdir"]) / "lake" / "safe_v1_test"
    good = base / "r_good" / "run"
    bad = base / "r_bad" / "run"

    _write_csv(good / "tb_30000" / "scan_tb_30000.csv", [_mkrow(triple_ok=True, br_gaga=1e-4, width_bb=1e-3, total_width=1e-2)])
    (good / "task_summary.jsonl").write_text('{"event":"done"}\n', encoding="utf-8")
    (good / "run_manifest.json").write_text(json.dumps({"runtime": {"omp_num_threads": 1}}), encoding="utf-8")
    (good / "orchestrator.log").write_text("[CONF] OMP_NUM_THREADS = 1\n", encoding="utf-8")

    _write_csv(bad / "tb_30000" / "scan_tb_30000.csv", [])
    (bad / "task_summary.jsonl").write_text('{"event":"error"}\n', encoding="utf-8")
    (bad / "run_manifest.json").write_text(json.dumps({"runtime": {"omp_num_threads": 8}}), encoding="utf-8")
    (bad / "orchestrator.log").write_text("bad\n", encoding="utf-8")

    out = run_safe_campaign(c, workdir=tmp_path, execute=False)
    assert out["stop_report"]["gates"]["passed"] is True
    roots = out["stop_report"]["artifact_roots_used_for_validation"]
    assert roots["run_filter"]["run_name"] == "r_good"
