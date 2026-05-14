from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any


_FAIL_EVENT_TOKENS = {"fail", "failed", "crash", "error", "exception", "timeout"}
_LOG_SIGNATURE_RE = re.compile(r"error|exception|failed|fail|timeout|traceback|crash", re.IGNORECASE)


def _bounded_text(value: object, *, limit: int = 500) -> str | None:
    if value is None:
        return None
    text = str(value)
    if len(text) <= limit:
        return text
    return text[: limit - 1] + "…"


def _event_kind(event: Mapping[str, object]) -> str:
    return str(event.get("event") or event.get("status") or event.get("outcome") or "").strip().lower()


def _looks_like_failure_event(kind: str, event: Mapping[str, object]) -> bool:
    if any(tok in kind for tok in _FAIL_EVENT_TOKENS):
        return True
    returncode = event.get("returncode")
    try:
        if returncode is not None and int(returncode) != 0:
            return True
    except Exception:
        pass
    for key in ("reason", "exception", "stderr_path", "stdout_path"):
        val = event.get(key)
        if isinstance(val, str) and any(tok in val.lower() for tok in _FAIL_EVENT_TOKENS):
            return True
    return False


def _extract_first_failure_event(events: list[Mapping[str, object]]) -> dict[str, object]:
    out: dict[str, object] = {
        "first_fail_event": None,
        "first_fail_event_index": None,
        "first_fail_reason": None,
        "first_fail_payload_snippet": None,
        "first_fail_stdout_path": None,
        "first_fail_stderr_path": None,
        "first_fail_output_csv": None,
        "first_fail_tanbeta": None,
        "first_fail_tb_tag": None,
    }
    for idx, ev in enumerate(events):
        kind = _event_kind(ev)
        if not _looks_like_failure_event(kind, ev):
            continue
        out["first_fail_event"] = kind or "unknown"
        out["first_fail_event_index"] = idx
        out["first_fail_reason"] = _bounded_text(ev.get("reason") or ev.get("exception") or ev.get("error"))
        out["first_fail_payload_snippet"] = _bounded_text(json.dumps(ev, sort_keys=True, default=str))
        out["first_fail_stdout_path"] = ev.get("stdout_path")
        out["first_fail_stderr_path"] = ev.get("stderr_path")
        out["first_fail_output_csv"] = ev.get("output_csv") or ev.get("csv_path")
        out["first_fail_tanbeta"] = ev.get("tan_beta")
        out["first_fail_tb_tag"] = ev.get("tb_tag")
        break
    return out


def inspect_orchestrator_log(path: str | Path, *, max_signatures: int = 10, tail_lines: int = 8) -> dict[str, object]:
    log_path = Path(path)
    out: dict[str, object] = {
        "orchestrator_log_path": str(log_path),
        "orchestrator_log_signature_count": 0,
        "orchestrator_log_signatures": [],
        "orchestrator_log_tail_relevant": [],
    }
    if not log_path.exists() or not log_path.is_file():
        return out
    hits: list[dict[str, object]] = []
    try:
        with log_path.open("r", encoding="utf-8", errors="replace") as handle:
            for line_no, line in enumerate(handle, start=1):
                line_stripped = line.rstrip("\n")
                if _LOG_SIGNATURE_RE.search(line_stripped):
                    hits.append({"line": line_no, "text": _bounded_text(line_stripped, limit=300)})
    except Exception:
        return out
    tail = hits[-tail_lines:]
    out["orchestrator_log_signature_count"] = len(hits)
    out["orchestrator_log_signatures"] = hits[-max_signatures:]
    out["orchestrator_log_tail_relevant"] = [f"L{item['line']}: {item['text']}" for item in tail]
    return out


def required_physics_columns() -> list[str]:
    return [
        "m_phi",
        "lambda6",
        "lambda7",
        "tan_beta",
        "lam1",
        "positivity_ok",
        "unitarity_ok",
        "perturbativity_ok",
        "width_bb",
        "width_gaga",
        "width_Zga",
        "width_gg",
        "width_WW",
        "width_ZZ",
        "width_hh",
        "total_width",
        "br_gaga",
    ]


def _new_csv_report(path: Path) -> dict[str, object]:
    return {
        "source_csv": str(path),
        "path": str(path),
        "exists": path.exists(),
        "is_file": path.is_file(),
        "has_header": False,
        "header_seen": [],
        "expected_columns": required_physics_columns(),
        "row_count": 0,
        "missing_columns": [],
        "extra_columns": [],
        "is_empty": False,
        "is_schema_valid": False,
        "is_scoreable": False,
        "status": "failed",
        "reason_codes": [],
        "first_error": None,
        "evidence_snippet": None,
        "suggested_next_check": "inspect_csv_schema_and_rows",
        "suggested_owner_module": "run_health",
        "warnings": [],
        "errors": [],
    }


def inspect_csv(path: str | Path, *, min_rows: int = 1) -> dict[str, object]:
    csv_path = Path(path)
    out = _new_csv_report(csv_path)

    warnings = out["warnings"]
    errors = out["errors"]
    assert isinstance(warnings, list)
    assert isinstance(errors, list)

    reason_codes: list[str] = []

    def _finish() -> dict[str, object]:
        out["reason_codes"] = sorted(set(reason_codes))
        out["first_error"] = errors[0] if errors else None
        if errors:
            out["status"] = "failed"
            out["is_scoreable"] = False
        return out

    if not csv_path.exists():
        errors.append(f"missing CSV file: {csv_path}")
        reason_codes.append("missing_csv")
        out["is_empty"] = True
        return _finish()
    if not csv_path.is_file():
        errors.append(f"CSV path is not a file: {csv_path}")
        reason_codes.append("csv_not_file")
        out["is_empty"] = True
        return _finish()

    try:
        with csv_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames
            if fieldnames is None:
                out["is_empty"] = True
                errors.append("CSV is empty or missing header")
                reason_codes.append("missing_header")
                return _finish()
            out["has_header"] = True
            out["header_seen"] = list(fieldnames)
            out["evidence_snippet"] = ",".join(fieldnames)

            required = required_physics_columns()
            missing = [c for c in required if c not in fieldnames]
            extra = [c for c in fieldnames if c not in required]
            out["missing_columns"] = missing
            out["extra_columns"] = extra

            row_count = sum(1 for _ in reader)
            out["row_count"] = row_count
            out["is_empty"] = row_count == 0

            if missing:
                errors.append(f"missing required columns: {', '.join(missing)}")
                reason_codes.append("missing_required_columns")
            if row_count < min_rows:
                errors.append(f"row_count={row_count} is below min_rows={min_rows}")
                reason_codes.append("row_count_below_min")
            if row_count == 0:
                reason_codes.append("empty_csv")
            if extra:
                warnings.append(f"extra columns present: {', '.join(extra)}")

            schema_valid = not missing and out["has_header"]
            out["is_schema_valid"] = schema_valid
            out["is_scoreable"] = schema_valid and row_count >= min_rows
            out["status"] = "ok" if out["is_scoreable"] else "failed"
            out["suggested_next_check"] = "score_csv" if out["status"] == "ok" else "fix_csv_schema_or_rerun_slice"
            if out["status"] == "ok":
                reason_codes.append("schema_valid")
            return _finish()
    except Exception as exc:
        out["is_empty"] = True
        errors.append(f"failed to parse CSV: {exc}")
        reason_codes.append("csv_parse_error")
        return _finish()


def inspect_task_summary(path: str | Path) -> dict[str, object]:
    summary_path = Path(path)
    out: dict[str, object] = {
        "path": str(summary_path),
        "exists": summary_path.exists(),
        "valid_jsonl": False,
        "line_count": 0,
        "events_seen": [],
        "has_done": False,
        "has_fail_or_crash": False,
        "status": "warning",
        "warnings": [],
        "errors": [],
        "first_fail_event": None,
        "first_fail_event_index": None,
        "first_fail_reason": None,
        "first_fail_payload_snippet": None,
        "first_fail_stdout_path": None,
        "first_fail_stderr_path": None,
        "first_fail_output_csv": None,
        "first_fail_tanbeta": None,
        "first_fail_tb_tag": None,
    }
    warnings = out["warnings"]
    errors = out["errors"]
    events_seen = out["events_seen"]
    assert isinstance(warnings, list)
    assert isinstance(errors, list)
    assert isinstance(events_seen, list)

    if not summary_path.exists():
        warnings.append(f"missing task_summary.jsonl: {summary_path}")
        out["status"] = "warning"
        return out
    if not summary_path.is_file():
        errors.append(f"task_summary path is not a file: {summary_path}")
        out["status"] = "failed"
        return out

    line_count = 0
    has_done = False
    has_fail_or_crash = False
    parsed_events: list[Mapping[str, object]] = []

    try:
        with summary_path.open("r", encoding="utf-8") as handle:
            for idx, line in enumerate(handle, start=1):
                text = line.strip()
                if text == "":
                    continue
                line_count += 1
                try:
                    event = json.loads(text)
                except json.JSONDecodeError as exc:
                    errors.append(f"invalid JSONL at line {idx}: {exc.msg}")
                    out["line_count"] = line_count
                    out["status"] = "failed"
                    out["valid_jsonl"] = False
                    return out

                if isinstance(event, Mapping):
                    kind = _event_kind(event)
                    parsed_events.append(event)
                else:
                    kind = ""
                if kind:
                    events_seen.append(kind)
                if kind == "done":
                    has_done = True
                if kind in {"failed", "fail", "crash", "error", "timeout"}:
                    has_fail_or_crash = True

        out["line_count"] = line_count
        out["valid_jsonl"] = True
        out["has_done"] = has_done
        out["has_fail_or_crash"] = has_fail_or_crash
        out.update(_extract_first_failure_event(parsed_events))

        if has_fail_or_crash:
            warnings.append("task_summary includes fail/crash/timeout events")
            out["status"] = "warning"
        else:
            out["status"] = "ok"
        return out
    except Exception as exc:
        errors.append(f"failed to read task_summary.jsonl: {exc}")
        out["status"] = "failed"
        out["valid_jsonl"] = False
        return out


def detect_dry_run(run_dir: str | Path) -> bool:
    run_path = Path(run_dir)

    manifest_path = run_path / "run_manifest.json"
    if manifest_path.exists() and manifest_path.is_file():
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            runtime = payload.get("runtime") if isinstance(payload, Mapping) else None
            if isinstance(runtime, Mapping) and bool(runtime.get("dry_run")):
                return True
        except Exception:
            pass

    task_path = run_path / "task_summary.jsonl"
    has_dry_event = False
    if task_path.exists() and task_path.is_file():
        try:
            with task_path.open("r", encoding="utf-8") as handle:
                for line in handle:
                    text = line.strip()
                    if text == "":
                        continue
                    try:
                        event = json.loads(text)
                    except json.JSONDecodeError:
                        continue
                    if not isinstance(event, Mapping):
                        continue
                    kind = str(event.get("event") or event.get("status") or event.get("outcome") or "").strip().lower()
                    if kind == "dry_run":
                        has_dry_event = True
                        break
        except Exception:
            pass
    if has_dry_event:
        return True

    csv_paths = sorted(run_path.glob("tb_*/scan_tb_*.csv"))
    total_rows = 0
    for csv_path in csv_paths:
        rep = inspect_csv(csv_path, min_rows=0)
        total_rows += int(rep.get("row_count", 0) or 0)
    if total_rows == 0 and has_dry_event:
        return True

    return False


def inspect_run_dir(run_dir: str | Path, *, min_rows_per_csv: int = 1) -> dict[str, object]:
    run_path = Path(run_dir)
    out: dict[str, object] = {
        "run_dir": str(run_path),
        "run_manifest_path": str(run_path / "run_manifest.json"),
        "task_summary_path": str(run_path / "task_summary.jsonl"),
        "exists": run_path.exists(),
        "has_manifest": (run_path / "run_manifest.json").exists(),
        "dry_run": False,
        "csv_paths": [],
        "csv_count": 0,
        "csv_ok_count": 0,
        "csv_failed_count": 0,
        "empty_csvs": [],
        "missing_column_files": [],
        "task_summary": inspect_task_summary(run_path / "task_summary.jsonl"),
        "csv_reports": [],
        "health_status": "failed",
        "is_scoreable": False,
        "reason_codes": [],
        "task_summary_events_seen": [],
        "fail_or_crash_markers": [],
        "evidence_paths": [],
        "source_code_hints": [
            "run_health.inspect_run_dir",
            "run_health.inspect_csv",
            "run_health.inspect_task_summary",
        ],
        "recommended_action": "missing_run_dir",
        "orchestrator_log_path": str(run_path / "orchestrator.log"),
        "orchestrator_log_signature_count": 0,
        "orchestrator_log_signatures": [],
        "orchestrator_log_tail_relevant": [],
        "first_fail_event": None,
        "first_fail_event_index": None,
        "first_fail_reason": None,
        "first_fail_payload_snippet": None,
        "first_fail_stdout_path": None,
        "first_fail_stderr_path": None,
        "first_fail_output_csv": None,
        "first_fail_tanbeta": None,
        "first_fail_tb_tag": None,
    }

    reason_codes: list[str] = []

    def _finish() -> dict[str, object]:
        task_summary = out.get("task_summary")
        if isinstance(task_summary, Mapping):
            events = task_summary.get("events_seen")
            if isinstance(events, list):
                out["task_summary_events_seen"] = list(events)
            if bool(task_summary.get("has_fail_or_crash")):
                out["fail_or_crash_markers"] = ["task_summary_fail_or_crash"]
            for key in (
                "first_fail_event",
                "first_fail_event_index",
                "first_fail_reason",
                "first_fail_payload_snippet",
                "first_fail_stdout_path",
                "first_fail_stderr_path",
                "first_fail_output_csv",
                "first_fail_tanbeta",
                "first_fail_tb_tag",
            ):
                out[key] = task_summary.get(key)

        out.update(inspect_orchestrator_log(run_path / "orchestrator.log"))

        ev_paths: list[str] = []
        for p in [out.get("run_manifest_path"), out.get("task_summary_path"), out.get("orchestrator_log_path")]:
            if isinstance(p, str) and Path(p).exists():
                ev_paths.append(p)
        for p in out.get("csv_paths", []):
            if isinstance(p, str) and Path(p).exists():
                ev_paths.append(p)
        out["evidence_paths"] = sorted(set(ev_paths))
        out["reason_codes"] = sorted(set(reason_codes))
        return out

    if not run_path.exists() or not run_path.is_dir():
        reason_codes.append("missing_run_dir")
        out["health_status"] = "failed"
        out["is_scoreable"] = False
        out["recommended_action"] = "missing_run_dir"
        return _finish()

    dry = detect_dry_run(run_path)
    out["dry_run"] = dry
    if dry:
        reason_codes.append("dry_run_rejected")
        out["health_status"] = "failed"
        out["is_scoreable"] = False
        out["recommended_action"] = "reject_dry_run"
        return _finish()

    csv_paths = sorted(run_path.glob("tb_*/scan_tb_*.csv"))
    out["csv_paths"] = [str(p) for p in csv_paths]
    out["csv_count"] = len(csv_paths)
    if len(csv_paths) == 0:
        reason_codes.append("missing_csv")
        out["health_status"] = "failed"
        out["is_scoreable"] = False
        out["recommended_action"] = "rerun_empty_or_invalid_csvs"
        return _finish()

    csv_reports: list[dict[str, object]] = [inspect_csv(path, min_rows=min_rows_per_csv) for path in csv_paths]
    out["csv_reports"] = csv_reports

    failed_reports = [rep for rep in csv_reports if rep.get("status") != "ok"]
    ok_reports = [rep for rep in csv_reports if rep.get("status") == "ok"]

    out["csv_ok_count"] = len(ok_reports)
    out["csv_failed_count"] = len(failed_reports)
    out["empty_csvs"] = [str(rep.get("path")) for rep in csv_reports if bool(rep.get("is_empty"))]
    out["missing_column_files"] = [
        str(rep.get("path"))
        for rep in csv_reports
        if isinstance(rep.get("missing_columns"), list) and len(rep.get("missing_columns")) > 0
    ]

    if failed_reports:
        for rep in failed_reports:
            if isinstance(rep.get("reason_codes"), list):
                reason_codes.extend([str(code) for code in rep["reason_codes"]])
        out["health_status"] = "failed"
        out["is_scoreable"] = False
        out["recommended_action"] = "rerun_empty_or_invalid_csvs"
        return _finish()

    task_summary = out["task_summary"]
    if isinstance(task_summary, Mapping) and task_summary.get("status") == "failed":
        reason_codes.append("task_summary_parse_failed")
        out["health_status"] = "failed"
        out["is_scoreable"] = False
        out["recommended_action"] = "inspect_task_summary"
        return _finish()

    if isinstance(task_summary, Mapping) and task_summary.get("exists") is False:
        reason_codes.append("task_summary_missing")
        out["health_status"] = "warning"
        out["is_scoreable"] = True
        out["recommended_action"] = "score"
        return _finish()

    if isinstance(task_summary, Mapping) and task_summary.get("status") == "warning":
        reason_codes.append("task_summary_warning")
        out["health_status"] = "warning"
        out["is_scoreable"] = True
        out["recommended_action"] = "inspect_task_summary"
        return _finish()

    out["health_status"] = "ok"
    out["is_scoreable"] = True
    out["recommended_action"] = "score"
    reason_codes.append("healthy_run")
    return _finish()


def summarize_health(reports: Iterable[Mapping[str, object]]) -> dict[str, object]:
    reports_list = list(reports)
    n_ok = sum(1 for r in reports_list if str(r.get("health_status")) == "ok")
    n_warning = sum(1 for r in reports_list if str(r.get("health_status")) == "warning")
    n_failed = sum(1 for r in reports_list if str(r.get("health_status")) == "failed")
    n_scoreable = sum(1 for r in reports_list if bool(r.get("is_scoreable")))
    n_dry = sum(1 for r in reports_list if bool(r.get("dry_run")))

    n_empty_csvs = 0
    failed_run_dirs: list[str] = []
    actions = Counter()
    for rep in reports_list:
        empty_csvs = rep.get("empty_csvs")
        if isinstance(empty_csvs, list):
            n_empty_csvs += len(empty_csvs)
        if str(rep.get("health_status")) == "failed":
            failed_run_dirs.append(str(rep.get("run_dir")))
        actions[str(rep.get("recommended_action"))] += 1

    return {
        "n_runs": len(reports_list),
        "n_ok": n_ok,
        "n_warning": n_warning,
        "n_failed": n_failed,
        "n_scoreable": n_scoreable,
        "n_dry_run": n_dry,
        "n_empty_csvs": n_empty_csvs,
        "failed_run_dirs": failed_run_dirs,
        "recommended_actions_count": dict(actions),
    }


def _cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Inspect DiHiggs autoresearch run health")
    parser.add_argument("--run-dir", required=True, help="Path to run_dir")
    parser.add_argument("--json", dest="json_out", default=None, help="Optional JSON output path")
    parser.add_argument("--min-rows-per-csv", type=int, default=1)
    args = parser.parse_args(argv)

    report = inspect_run_dir(args.run_dir, min_rows_per_csv=args.min_rows_per_csv)

    concise = {
        "run_dir": report["run_dir"],
        "health_status": report["health_status"],
        "is_scoreable": report["is_scoreable"],
        "recommended_action": report["recommended_action"],
        "csv_count": report["csv_count"],
        "csv_failed_count": report["csv_failed_count"],
        "dry_run": report["dry_run"],
    }
    print(json.dumps(concise, indent=2))

    if args.json_out:
        out_path = Path(args.json_out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    if bool(report.get("is_scoreable")):
        if str(report.get("health_status")) == "ok":
            return 0
        return 1
    return 2


if __name__ == "__main__":
    raise SystemExit(_cli())
