from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path

from autoresearch.harness.campaign_state import (
    build_campaign_state,
    build_run_record,
    build_slice_records,
    detect_report_staleness,
    discover_run_dirs,
    file_sha256,
    summarize_score,
    write_campaign_state_csv,
    write_campaign_state_json,
    write_campaign_state_markdown,
)
from autoresearch.harness.run_health import required_physics_columns


def _write_csv(path: Path, rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as h:
        w = csv.writer(h)
        w.writerow(required_physics_columns())
        w.writerows(rows)


def _minimal_row() -> list[str]:
    return [
        "250.0", "0.001", "0.0", "30000", "5.0", "1", "1", "1", "1e-12", "1e-14",
        "1e-15", "1e-13", "1e-13", "1e-13", "1e-13", "2e-12", "0.01",
    ]


def _make_run(run_dir: Path, *, ok: bool = True, dry: bool = False, empty: bool = False) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    if dry:
        (run_dir / "run_manifest.json").write_text(json.dumps({"runtime": {"dry_run": True}}), encoding="utf-8")
    (run_dir / "task_summary.jsonl").write_text('{"event":"done"}\n', encoding="utf-8")
    if ok or empty:
        rows = [] if empty else [_minimal_row()]
        _write_csv(run_dir / "tb_30000" / "scan_tb_30000.csv", rows)
    return run_dir


def _write_registry(path: Path, proposals: list[dict[str, object]]) -> None:
    lines = []
    for p in proposals:
        lines.append(json.dumps({"event_type": "create_proposal", "proposal": p}, sort_keys=True))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _proposal(pid: str, run_dir: str | None) -> dict[str, object]:
    return {
        "schema_version": "1.0",
        "proposal_id": pid,
        "parent_id": None,
        "generation": 0,
        "operator": "autoresearch",
        "bounds": {"m_phi": [200.0, 290.0], "tan_beta": [10000.0, 40000.0], "lambda6": [0.0006, 0.0006], "lambda7": [0.0, 0.0], "lam1": [4.0, 12.0]},
        "fixed": {"sin_ba": 1.0, "mA": 300.0, "mHp": 300.0, "yukawa_type": 2, "lambda7": 0.0},
        "resolution": {"m_phi": 4, "tan_beta": 4, "lambda6": 2, "lambda7": 1, "lam1": 4},
        "budget": {"requested_points": 100, "max_walltime_sec": 60.0},
        "status": "planned",
        "created_utc": "2026-01-01T00:00:00Z",
        "updated_utc": "2026-01-01T00:00:00Z",
        "run_dir": run_dir,
        "metrics": {},
        "error": None,
    }


def test_file_sha256_behaviors(tmp_path: Path) -> None:
    p = tmp_path / "x.txt"
    assert file_sha256(p) is None
    p.write_text("abc", encoding="utf-8")
    h1 = file_sha256(p)
    h2 = file_sha256(p)
    assert h1 is not None and h1 == h2
    p.write_text("abcd", encoding="utf-8")
    assert file_sha256(p) != h1


def test_discover_run_dirs_explicit_sorted_dedup_and_campaign_bounded(tmp_path: Path) -> None:
    run_a = _make_run(tmp_path / "camp" / "a")
    run_b = _make_run(tmp_path / "camp" / "nested" / "b")
    unrelated = tmp_path / "camp" / "zzz"
    unrelated.mkdir(parents=True)
    (unrelated / "deep" / "x").mkdir(parents=True)
    dirs = discover_run_dirs(explicit_run_dirs=[run_b, run_a, run_b], campaign_dir=tmp_path / "camp")
    assert dirs == sorted(set(dirs), key=lambda p: str(p))
    assert run_a.resolve() in dirs
    assert run_b.resolve() in dirs
    assert all("deep/x" not in str(d) for d in dirs)


def test_build_slice_records_one_csv_hash_rows_and_failures(tmp_path: Path) -> None:
    run = _make_run(tmp_path / "run1")
    slices = build_slice_records(run, proposal_id="p1")
    assert len(slices) == 1
    s = slices[0]
    assert s["csv_hash"]
    assert s["row_count"] == 1

    bad = tmp_path / "run_bad"
    (bad / "tb_30000").mkdir(parents=True)
    with (bad / "tb_30000" / "scan_tb_30000.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.writer(h)
        header = required_physics_columns()
        header.remove("br_gaga")
        w.writerow(header)
        w.writerow(_minimal_row()[:-1])
    s_bad = build_slice_records(bad)[0]
    assert s_bad["scoreable"] is False
    assert s_bad["reason_codes"]
    assert s_bad["suggested_next_check"]

    empty = _make_run(tmp_path / "run_empty", empty=True)
    s_empty = build_slice_records(empty)[0]
    assert s_empty["is_empty"] is True
    assert s_empty["health_status"] == "failed"


def test_summarize_score_shapes() -> None:
    flat = summarize_score({"rankable": True, "csv_rows": 3, "br_gaga_q95": 0.1, "warnings": []})
    assert flat["score_status"] == "scored"
    nested = summarize_score({"rankable": False, "aggregate_metrics": {"csv_rows": 2, "br_bb_median": 0.2}})
    assert nested["csv_rows"] == 2
    none = summarize_score(None)
    assert none["score_status"] == "missing"


def test_build_run_record_recommended_actions_and_limitations(tmp_path: Path) -> None:
    dry = _make_run(tmp_path / "dry", dry=True)
    rec = build_run_record(dry)
    assert rec["dry_run"] is True
    assert rec["recommended_action"] == "ignore_for_scoring"

    failed = _make_run(tmp_path / "failed", empty=True)
    rec_f = build_run_record(failed)
    assert rec_f["recommended_action"] == "needs_rerun_or_inspection"
    assert rec_f["reason_codes"]
    assert rec_f["evidence_paths"]
    assert "first_fail_event" in rec_f
    assert "orchestrator_log_signature_count" in rec_f

    ok = _make_run(tmp_path / "ok")
    score = {"rankable": True, "aggregate_metrics": {"csv_rows": 1}}
    rec_ok = build_run_record(ok, score=score)
    assert rec_ok["recommended_action"] == "ready_for_review"

    rec_missing_score = build_run_record(ok, score=None)
    assert rec_missing_score["recommended_action"] in {"score", "review_warning_then_score_or_accept"}
    assert "coverage_not_computed" in rec_missing_score["limitations"]
    assert "novelty_not_computed" in rec_missing_score["limitations"]


def test_detect_report_staleness_cases(tmp_path: Path) -> None:
    src = tmp_path / "src.txt"
    src.write_text("x", encoding="utf-8")
    missing = detect_report_staleness(tmp_path / "missing.md", source_paths=[src])
    assert missing["exists"] is False

    rep = tmp_path / "report.md"
    time.sleep(0.01)
    rep.write_text("r", encoding="utf-8")
    fresh = detect_report_staleness(rep, source_paths=[src])
    assert fresh["is_stale"] is False

    time.sleep(0.01)
    src.write_text("new", encoding="utf-8")
    stale = detect_report_staleness(rep, source_paths=[src])
    assert stale["is_stale"] is True


def test_build_campaign_state_mixed_and_association(tmp_path: Path) -> None:
    campaign = tmp_path / "campaign"
    run_dry = _make_run(campaign / "run_p_dry", dry=True)
    run_ok = _make_run(campaign / "run_p_good")
    run_fail = _make_run(campaign / "run_orphan_fail", empty=True)

    reg = tmp_path / "proposal_registry.jsonl"
    p1 = _proposal("p_dry", str(run_dry))
    p2 = _proposal("p_good", str(run_ok))
    p3 = _proposal("p_planned_only", None)
    _write_registry(reg, [p1, p2, p3])

    report = tmp_path / "old_report.md"
    report.write_text("old", encoding="utf-8")
    time.sleep(0.01)
    reg.write_text(reg.read_text(encoding="utf-8") + "", encoding="utf-8")

    state = build_campaign_state(
        campaign="synthetic",
        registry_path=reg,
        campaign_dir=campaign,
        report_paths=[report],
        score_runs=True,
    )

    assert "proposals" in state and "runs" in state and "slices" in state and "reports" in state
    summary = state["summary"]
    assert summary["n_proposals"] == 3
    assert summary["n_runs_found"] == 3
    assert summary["n_health_failed"] >= 2
    assert any("coverage_not_implemented" in w for w in state["warnings"])
    assert any("novelty_not_implemented" in w for w in state["warnings"])
    assert any("proposal_generation_disabled" in w for w in state["warnings"])
    assert any("orphan_run_dir" in w for w in state["warnings"])


def test_writers_and_markdown_guard(tmp_path: Path) -> None:
    run = _make_run(tmp_path / "r")
    state = build_campaign_state(campaign="x", explicit_run_dirs=[run], score_runs=False)
    p_json = tmp_path / "state.json"
    p_csv = tmp_path / "state.csv"
    p_md = tmp_path / "state.md"
    write_campaign_state_json(state, p_json)
    write_campaign_state_csv(state, p_csv)
    write_campaign_state_markdown(state, p_md)

    loaded = json.loads(p_json.read_text(encoding="utf-8"))
    assert loaded["campaign"] == "x"

    with p_csv.open("r", encoding="utf-8", newline="") as h:
        hdr = next(csv.reader(h))
    for col in [
        "proposal_id", "operator", "run_dir", "dry_run", "health_status", "health_scoreable", "health_recommended_action",
        "score_status", "report_status", "staleness_status", "csv_rows", "triple_ok_rate", "br_gaga_q95", "br_bb_median", "total_width_median", "recommended_action", "limitations",
    ]:
        assert col in hdr

    md = p_md.read_text(encoding="utf-8")
    assert "Coverage and novelty are not implemented yet" in md
    assert "| Metric | Value |" in md


def test_cli_outputs_and_exit_codes(tmp_path: Path) -> None:
    bad = _make_run(tmp_path / "bad", empty=True)
    outj = tmp_path / "o.json"
    outm = tmp_path / "o.md"
    outc = tmp_path / "o.csv"
    p = subprocess.run(
        [
            sys.executable,
            "-m",
            "autoresearch.harness.campaign_state",
            "--campaign",
            "cli_bad",
            "--run-dir",
            str(bad),
            "--out-json",
            str(outj),
            "--out-md",
            str(outm),
            "--out-csv",
            str(outc),
            "--no-score",
        ],
        capture_output=True,
        text=True,
    )
    assert p.returncode == 2
    assert outj.exists() and outm.exists() and outc.exists()

    good = _make_run(tmp_path / "good")
    p2 = subprocess.run(
        [sys.executable, "-m", "autoresearch.harness.campaign_state", "--campaign", "cli_ok", "--run-dir", str(good), "--no-score"],
        capture_output=True,
        text=True,
    )
    assert p2.returncode == 0


def test_power_tool_smoke_questions_answerable(tmp_path: Path) -> None:
    run_ok = _make_run(tmp_path / "run_ok")
    run_fail = _make_run(tmp_path / "run_fail", empty=True)
    state = build_campaign_state(campaign="q", explicit_run_dirs=[run_ok, run_fail], score_runs=True)
    summary = state["summary"]
    assert summary["n_runs_found"] == 2  # what exists / what ran
    assert summary["n_health_failed"] >= 1  # what failed
    assert summary["n_scoreable"] >= 1  # what is scoreable
    assert summary["n_scored"] >= 1  # what is scored
    assert summary["n_needs_rerun"] >= 1  # what needs rerun
    assert isinstance(summary["n_ready_for_review"], int)  # ready for review signal
    assert state["capabilities"]["proposal_generation_enabled"] is False  # why proposal generation disabled
