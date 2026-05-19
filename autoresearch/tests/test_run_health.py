from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

from autoresearch.harness.run_health import (
    detect_dry_run,
    inspect_csv,
    inspect_run_dir,
    inspect_task_summary,
    required_physics_columns,
    summarize_health,
)


def _write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def _minimal_header() -> list[str]:
    return required_physics_columns()


def _minimal_row() -> list[str]:
    return [
        "250.0", "0.001", "0.0", "30000", "5.0", "1", "1", "1", "1e-12", "1e-14",
        "1e-15", "1e-13", "1e-13", "1e-13", "1e-13", "2e-12", "0.01",
    ]


def test_inspect_csv_missing_file_failed(tmp_path: Path) -> None:
    rep = inspect_csv(tmp_path / "missing.csv")
    assert rep["status"] == "failed"
    assert rep["is_scoreable"] is False
    assert rep["exists"] is False


def test_inspect_csv_empty_file_failed(tmp_path: Path) -> None:
    p = tmp_path / "x.csv"
    p.write_text("", encoding="utf-8")
    rep = inspect_csv(p)
    assert rep["status"] == "failed"
    assert rep["is_empty"] is True


def test_inspect_csv_header_only_failed(tmp_path: Path) -> None:
    p = tmp_path / "h.csv"
    _write_csv(p, _minimal_header(), [])
    rep = inspect_csv(p)
    assert rep["status"] == "failed"
    assert rep["row_count"] == 0
    assert "empty_csv" in rep["reason_codes"]
    assert rep["suggested_next_check"]


def test_inspect_csv_valid_minimal_ok_scoreable(tmp_path: Path) -> None:
    p = tmp_path / "ok.csv"
    _write_csv(p, _minimal_header(), [_minimal_row()])
    rep = inspect_csv(p)
    assert rep["status"] == "ok"
    assert rep["is_schema_valid"] is True
    assert rep["is_scoreable"] is True


def test_inspect_csv_missing_required_column_failed(tmp_path: Path) -> None:
    p = tmp_path / "bad.csv"
    header = _minimal_header()
    header.remove("br_gaga")
    row = _minimal_row()[:-1]
    _write_csv(p, header, [row])
    rep = inspect_csv(p)
    assert rep["status"] == "failed"
    assert "br_gaga" in rep["missing_columns"]
    assert "missing_required_columns" in rep["reason_codes"]
    assert "br_gaga" not in str(rep["evidence_snippet"])


def test_inspect_csv_extra_columns_allowed(tmp_path: Path) -> None:
    p = tmp_path / "extra.csv"
    header = _minimal_header() + ["extra_col"]
    row = _minimal_row() + ["x"]
    _write_csv(p, header, [row])
    rep = inspect_csv(p)
    assert rep["status"] == "ok"
    assert "extra_col" in rep["extra_columns"]


def test_inspect_csv_min_rows_enforced(tmp_path: Path) -> None:
    p = tmp_path / "minrows.csv"
    _write_csv(p, _minimal_header(), [_minimal_row()])
    rep = inspect_csv(p, min_rows=2)
    assert rep["status"] == "failed"
    assert rep["is_scoreable"] is False


def test_inspect_task_summary_missing_warning(tmp_path: Path) -> None:
    rep = inspect_task_summary(tmp_path / "task_summary.jsonl")
    assert rep["status"] == "warning"
    assert rep["exists"] is False


def test_inspect_task_summary_valid_done_ok(tmp_path: Path) -> None:
    p = tmp_path / "task_summary.jsonl"
    p.write_text('{"event":"done"}\n', encoding="utf-8")
    rep = inspect_task_summary(p)
    assert rep["status"] == "ok"
    assert rep["valid_jsonl"] is True
    assert rep["has_done"] is True


def test_inspect_task_summary_invalid_json_failed(tmp_path: Path) -> None:
    p = tmp_path / "task_summary.jsonl"
    p.write_text('{"event":"done"}\n{bad}\n', encoding="utf-8")
    rep = inspect_task_summary(p)
    assert rep["status"] == "failed"
    assert rep["valid_jsonl"] is False


def test_inspect_task_summary_fail_or_crash_flagged(tmp_path: Path) -> None:
    p = tmp_path / "task_summary.jsonl"
    p.write_text('{"event":"crash"}\n', encoding="utf-8")
    rep = inspect_task_summary(p)
    assert rep["has_fail_or_crash"] is True
    assert rep["status"] == "warning"


def test_first_fail_event_from_fail_record_and_bounded_snippet(tmp_path: Path) -> None:
    p = tmp_path / "task_summary.jsonl"
    long_reason = "x" * 1200
    p.write_text(
        '{"event":"start"}\n'
        + json.dumps({"event": "fail", "reason": long_reason, "stdout_path": "o.log", "stderr_path": "e.log", "tan_beta": 30000})
        + "\n",
        encoding="utf-8",
    )
    rep = inspect_task_summary(p)
    assert rep["first_fail_event"] == "fail"
    assert rep["first_fail_event_index"] == 1
    assert rep["first_fail_stdout_path"] == "o.log"
    assert rep["first_fail_stderr_path"] == "e.log"
    assert isinstance(rep["first_fail_payload_snippet"], str)
    assert len(str(rep["first_fail_payload_snippet"])) <= 500


def test_first_fail_event_from_crash_exception_record(tmp_path: Path) -> None:
    p = tmp_path / "task_summary.jsonl"
    p.write_text(json.dumps({"event": "crash", "exception": "segfault"}) + "\n", encoding="utf-8")
    rep = inspect_task_summary(p)
    assert rep["first_fail_event"] == "crash"
    assert rep["first_fail_reason"] == "segfault"


def test_orchestrator_log_signatures_are_bounded(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    tb = run_dir / "tb_30000"
    tb.mkdir(parents=True)
    (run_dir / "task_summary.jsonl").write_text('{"event":"done"}\n', encoding="utf-8")
    _write_csv(tb / "scan_tb_30000.csv", _minimal_header(), [_minimal_row()])
    (run_dir / "orchestrator.log").write_text(
        "ok\nERROR bad thing\nmore\nTraceback: fail\n",
        encoding="utf-8",
    )
    rep = inspect_run_dir(run_dir)
    assert rep["orchestrator_log_signature_count"] == 2
    assert len(rep["orchestrator_log_tail_relevant"]) <= 8
    assert len(rep["orchestrator_log_signatures"]) <= 10


def test_detect_dry_run_manifest_flag_true(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "run_manifest.json").write_text(json.dumps({"runtime": {"dry_run": True}}), encoding="utf-8")
    assert detect_dry_run(run_dir) is True


def test_detect_dry_run_task_summary_event_true(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "task_summary.jsonl").write_text('{"event":"dry_run"}\n', encoding="utf-8")
    assert detect_dry_run(run_dir) is True


def test_detect_dry_run_normal_done_with_rows_false(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    csv_dir = run_dir / "tb_30000"
    csv_dir.mkdir(parents=True)
    (run_dir / "task_summary.jsonl").write_text('{"event":"done"}\n', encoding="utf-8")
    _write_csv(csv_dir / "scan_tb_30000.csv", _minimal_header(), [_minimal_row()])
    assert detect_dry_run(run_dir) is False


def test_inspect_run_dir_missing_failed(tmp_path: Path) -> None:
    rep = inspect_run_dir(tmp_path / "missing")
    assert rep["health_status"] == "failed"
    assert rep["is_scoreable"] is False
    assert rep["recommended_action"] == "missing_run_dir"


def test_inspect_run_dir_dry_run_not_scoreable(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "run_manifest.json").write_text(json.dumps({"runtime": {"dry_run": True}}), encoding="utf-8")
    rep = inspect_run_dir(run_dir)
    assert rep["dry_run"] is True
    assert rep["is_scoreable"] is False
    assert rep["recommended_action"] == "reject_dry_run"
    assert "dry_run_rejected" in rep["reason_codes"]


def test_inspect_run_dir_no_csvs_failed(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "task_summary.jsonl").write_text('{"event":"done"}\n', encoding="utf-8")
    rep = inspect_run_dir(run_dir)
    assert rep["health_status"] == "failed"
    assert rep["csv_count"] == 0


def test_inspect_run_dir_one_empty_csv_among_valid_failed(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    tb1 = run_dir / "tb_30000"
    tb2 = run_dir / "tb_40000"
    tb1.mkdir(parents=True)
    tb2.mkdir(parents=True)
    (run_dir / "task_summary.jsonl").write_text('{"event":"done"}\n', encoding="utf-8")
    _write_csv(tb1 / "scan_tb_30000.csv", _minimal_header(), [_minimal_row()])
    _write_csv(tb2 / "scan_tb_40000.csv", _minimal_header(), [])
    rep = inspect_run_dir(run_dir)
    assert rep["health_status"] == "failed"
    assert rep["recommended_action"] == "rerun_empty_or_invalid_csvs"
    assert len(rep["empty_csvs"]) == 1
    assert rep["reason_codes"]
    assert rep["evidence_paths"]


def test_inspect_run_dir_all_valid_missing_task_summary_warning_but_scoreable(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    tb = run_dir / "tb_30000"
    tb.mkdir(parents=True)
    _write_csv(tb / "scan_tb_30000.csv", _minimal_header(), [_minimal_row()])
    rep = inspect_run_dir(run_dir)
    assert rep["health_status"] == "warning"
    assert rep["is_scoreable"] is True
    assert rep["recommended_action"] == "score"


def test_inspect_run_dir_all_valid_with_done_ok(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    tb = run_dir / "tb_30000"
    tb.mkdir(parents=True)
    (run_dir / "task_summary.jsonl").write_text('{"event":"done"}\n', encoding="utf-8")
    _write_csv(tb / "scan_tb_30000.csv", _minimal_header(), [_minimal_row()])
    rep = inspect_run_dir(run_dir)
    assert rep["health_status"] == "ok"
    assert rep["is_scoreable"] is True
    assert rep["recommended_action"] == "score"
    assert "healthy_run" in rep["reason_codes"]


def test_task_summary_fail_marker_has_evidence_fields(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    tb = run_dir / "tb_30000"
    tb.mkdir(parents=True)
    (run_dir / "task_summary.jsonl").write_text('{"event":"timeout"}\n', encoding="utf-8")
    _write_csv(tb / "scan_tb_30000.csv", _minimal_header(), [_minimal_row()])
    rep = inspect_run_dir(run_dir)
    assert rep["health_status"] == "warning"
    assert rep["task_summary_events_seen"]
    assert rep["fail_or_crash_markers"]


def test_summarize_health_aggregate_counts() -> None:
    reports = [
        {"run_dir": "a", "health_status": "ok", "is_scoreable": True, "dry_run": False, "empty_csvs": [], "recommended_action": "score"},
        {"run_dir": "b", "health_status": "warning", "is_scoreable": True, "dry_run": False, "empty_csvs": [], "recommended_action": "inspect_task_summary"},
        {"run_dir": "c", "health_status": "failed", "is_scoreable": False, "dry_run": True, "empty_csvs": ["x"], "recommended_action": "reject_dry_run"},
    ]
    s = summarize_health(reports)
    assert s["n_runs"] == 3
    assert s["n_ok"] == 1
    assert s["n_warning"] == 1
    assert s["n_failed"] == 1
    assert s["n_scoreable"] == 2
    assert s["n_dry_run"] == 1
    assert s["n_empty_csvs"] == 1


def test_cli_scoreable_exit0_and_json_output(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    tb = run_dir / "tb_30000"
    tb.mkdir(parents=True)
    (run_dir / "task_summary.jsonl").write_text('{"event":"done"}\n', encoding="utf-8")
    _write_csv(tb / "scan_tb_30000.csv", _minimal_header(), [_minimal_row()])

    out_json = tmp_path / "out.json"
    cmd = [sys.executable, "-m", "autoresearch.harness.run_health", "--run-dir", str(run_dir), "--json", str(out_json)]
    p = subprocess.run(cmd, text=True, capture_output=True)
    assert p.returncode == 0
    payload = json.loads(out_json.read_text(encoding="utf-8"))
    assert payload["is_scoreable"] is True


def test_cli_non_scoreable_exit2(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    cmd = [sys.executable, "-m", "autoresearch.harness.run_health", "--run-dir", str(run_dir)]
    p = subprocess.run(cmd, text=True, capture_output=True)
    assert p.returncode == 2
