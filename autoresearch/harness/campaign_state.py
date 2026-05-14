from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from autoresearch.harness.dihiggs_physics_score import score_run_dir_checked
from autoresearch.harness.proposal_registry import load_events, replay_registry
from autoresearch.harness.run_health import inspect_csv, inspect_run_dir


DISABLED_CAPABILITY_WARNING = (
    "Coverage and novelty are not implemented yet. Existing results must not be interpreted as "
    "representative coverage of the full parameter space. Autonomous proposal generation is disabled."
)


_REQUIRED_CSV_COLUMNS = [
    "proposal_id",
    "operator",
    "run_dir",
    "dry_run",
    "health_status",
    "health_scoreable",
    "health_recommended_action",
    "score_status",
    "report_status",
    "staleness_status",
    "csv_rows",
    "triple_ok_rate",
    "br_gaga_q95",
    "br_bb_median",
    "total_width_median",
    "recommended_action",
    "first_fail_event",
    "first_fail_reason",
    "first_fail_stdout_path",
    "first_fail_stderr_path",
    "orchestrator_log_path",
    "orchestrator_log_signature_count",
    "limitations",
]


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def file_sha256(path: str | Path) -> str | None:
    p = Path(path)
    if not p.exists() or not p.is_file():
        return None
    try:
        h = hashlib.sha256()
        with p.open("rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return None


def load_registry_state(registry_path: str | Path | None, contract: Mapping[str, object] | None = None) -> dict[str, dict]:
    if registry_path is None:
        return {}
    p = Path(registry_path)
    if not p.exists():
        return {}
    try:
        return replay_registry(load_events(p), contract=contract)
    except Exception as exc:
        raise ValueError(f"failed to load/replay registry {p}: {exc}") from exc


def _is_plausible_run_dir(path: Path) -> bool:
    if not path.exists() or not path.is_dir():
        return False
    if (path / "run_manifest.json").exists() or (path / "task_summary.jsonl").exists():
        return True
    return any(path.glob("tb_*/scan_tb_*.csv"))


def discover_run_dirs(
    *,
    explicit_run_dirs: Iterable[str | Path] | None = None,
    campaign_dir: str | Path | None = None,
) -> list[Path]:
    found: set[Path] = set()

    for item in explicit_run_dirs or ():
        p = Path(item)
        if p.exists() and p.is_dir():
            found.add(p.resolve())

    if campaign_dir is not None:
        base = Path(campaign_dir)
        if base.exists() and base.is_dir():
            candidates: list[Path] = [base]
            candidates.extend([x for x in base.iterdir() if x.is_dir()])
            for child in [x for x in base.iterdir() if x.is_dir()]:
                for gc in child.iterdir():
                    if gc.is_dir():
                        candidates.append(gc)
            for c in candidates:
                if _is_plausible_run_dir(c):
                    found.add(c.resolve())

    return sorted(found, key=lambda p: str(p))


def _derive_tan_beta(csv_path: Path) -> float | None:
    tb_dir = csv_path.parent.name
    m = re.search(r"tb_([0-9eE+\-.]+)", tb_dir)
    if m:
        try:
            return float(m.group(1))
        except Exception:
            pass
    m = re.search(r"scan_tb_([0-9eE+\-.]+)\.csv$", csv_path.name)
    if m:
        try:
            return float(m.group(1))
        except Exception:
            pass
    try:
        with csv_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            first = next(reader, None)
            if first and first.get("tan_beta") not in (None, ""):
                return float(str(first["tan_beta"]))
    except Exception:
        return None
    return None


def build_slice_records(run_dir: str | Path, proposal_id: str | None = None) -> list[dict[str, object]]:
    run_path = Path(run_dir)
    rows: list[dict[str, object]] = []
    for csv_path in sorted(run_path.glob("tb_*/scan_tb_*.csv")):
        rep = inspect_csv(csv_path, min_rows=1)
        scoreable = bool(rep.get("is_scoreable"))
        is_empty = bool(rep.get("is_empty"))
        health_status = str(rep.get("status") or "failed")
        if is_empty or not scoreable:
            action = "needs_rerun_or_inspection"
        else:
            action = "score_or_include"
        rows.append(
            {
                "proposal_id": proposal_id,
                "run_dir": str(run_path),
                "source_csv": str(csv_path),
                "tb_dir": csv_path.parent.name,
                "tan_beta": _derive_tan_beta(csv_path),
                "row_count": int(rep.get("row_count", 0) or 0),
                "header_seen": list(rep.get("header_seen", [])) if isinstance(rep.get("header_seen"), list) else [],
                "expected_columns": list(rep.get("expected_columns", [])) if isinstance(rep.get("expected_columns"), list) else [],
                "reason_codes": list(rep.get("reason_codes", [])) if isinstance(rep.get("reason_codes"), list) else [],
                "first_error": rep.get("first_error"),
                "evidence_snippet": rep.get("evidence_snippet"),
                "suggested_next_check": rep.get("suggested_next_check"),
                "suggested_owner_module": rep.get("suggested_owner_module"),
                "csv_hash": file_sha256(csv_path),
                "health_status": health_status,
                "scoreable": scoreable,
                "missing_columns": list(rep.get("missing_columns", [])) if isinstance(rep.get("missing_columns"), list) else [],
                "extra_columns": list(rep.get("extra_columns", [])) if isinstance(rep.get("extra_columns"), list) else [],
                "is_empty": is_empty,
                "recommended_action": action,
                "evidence_paths": [str(csv_path)],
            }
        )
    return rows


def summarize_score(score: Mapping[str, object] | None) -> dict[str, object]:
    if score is None:
        return {
            "score_status": "missing",
            "rankable": False,
            "csv_rows": 0,
            "triple_ok_rate": None,
            "br_gaga_q95": None,
            "br_gaga_median": None,
            "br_bb_median": None,
            "total_width_median": None,
            "objective_score": None,
            "warnings": ["score_missing"],
        }
    src = score.get("aggregate_metrics") if isinstance(score.get("aggregate_metrics"), Mapping) else score
    assert isinstance(src, Mapping)
    warnings = src.get("warnings")
    return {
        "score_status": "scored",
        "rankable": bool(score.get("rankable", src.get("rankable"))),
        "csv_rows": int(src.get("csv_rows", 0) or 0),
        "triple_ok_rate": src.get("triple_ok_rate_over_csv_rows", src.get("triple_ok_rate")),
        "br_gaga_q95": src.get("br_gaga_q95"),
        "br_gaga_median": src.get("br_gaga_median"),
        "br_bb_median": src.get("br_bb_median"),
        "total_width_median": src.get("total_width_median"),
        "objective_score": src.get("objective_score"),
        "warnings": list(warnings) if isinstance(warnings, list) else [],
    }


def build_run_record(
    run_dir: str | Path,
    *,
    proposal_id: str | None = None,
    proposal: Mapping[str, object] | None = None,
    score: Mapping[str, object] | None = None,
    min_rows_per_csv: int = 1,
    allow_warnings: bool = True,
) -> dict[str, object]:
    run_path = Path(run_dir)
    health = inspect_run_dir(run_path, min_rows_per_csv=min_rows_per_csv)
    score_s = summarize_score(score)

    dry_run = bool(health.get("dry_run"))
    health_status = str(health.get("health_status") or "failed")
    health_scoreable = bool(health.get("is_scoreable"))
    report_status = "missing"
    staleness_status = "unknown"

    if dry_run:
        recommended_action = "ignore_for_scoring"
    elif health_status == "failed":
        recommended_action = "needs_rerun_or_inspection"
    elif health_status == "warning" and health_scoreable:
        recommended_action = "review_warning_then_score_or_accept"
    elif health_scoreable and score is None:
        recommended_action = "score"
    elif score is not None:
        recommended_action = "ready_for_review"
    else:
        recommended_action = "inspect"

    limitations = ["coverage_not_computed", "novelty_not_computed", "proposal_generation_disabled"]

    cmd_fingerprint = None
    manifest_path = run_path / "run_manifest.json"
    if manifest_path.exists() and manifest_path.is_file():
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            cmd = None
            if isinstance(payload, Mapping):
                cmd = payload.get("command")
                if cmd is None and isinstance(payload.get("runtime"), Mapping):
                    cmd = payload["runtime"].get("command")  # type: ignore[index]
            if cmd is not None:
                cmd_fingerprint = hashlib.sha256(json.dumps(cmd, sort_keys=True, default=str).encode("utf-8")).hexdigest()
        except Exception:
            cmd_fingerprint = None

    out: dict[str, object] = {
        "proposal_id": proposal_id,
        "operator": proposal.get("operator") if isinstance(proposal, Mapping) else None,
        "run_dir": str(run_path),
        "run_name": run_path.name,
        "run_manifest_path": health.get("run_manifest_path"),
        "task_summary_path": health.get("task_summary_path"),
        "csv_paths": list(health.get("csv_paths", [])) if isinstance(health.get("csv_paths"), list) else [],
        "dry_run": dry_run,
        "health_status": health_status,
        "health_scoreable": health_scoreable,
        "health_recommended_action": health.get("recommended_action"),
        "reason_codes": list(health.get("reason_codes", [])) if isinstance(health.get("reason_codes"), list) else [],
        "task_summary_events_seen": list(health.get("task_summary_events_seen", [])) if isinstance(health.get("task_summary_events_seen"), list) else [],
        "fail_or_crash_markers": list(health.get("fail_or_crash_markers", [])) if isinstance(health.get("fail_or_crash_markers"), list) else [],
        "first_fail_event": health.get("first_fail_event"),
        "first_fail_event_index": health.get("first_fail_event_index"),
        "first_fail_reason": health.get("first_fail_reason"),
        "first_fail_payload_snippet": health.get("first_fail_payload_snippet"),
        "first_fail_stdout_path": health.get("first_fail_stdout_path"),
        "first_fail_stderr_path": health.get("first_fail_stderr_path"),
        "first_fail_output_csv": health.get("first_fail_output_csv"),
        "first_fail_tanbeta": health.get("first_fail_tanbeta"),
        "first_fail_tb_tag": health.get("first_fail_tb_tag"),
        "orchestrator_log_path": health.get("orchestrator_log_path"),
        "orchestrator_log_signature_count": int(health.get("orchestrator_log_signature_count", 0) or 0),
        "orchestrator_log_tail_relevant": list(health.get("orchestrator_log_tail_relevant", [])) if isinstance(health.get("orchestrator_log_tail_relevant"), list) else [],
        "source_code_hints": list(health.get("source_code_hints", [])) if isinstance(health.get("source_code_hints"), list) else [],
        "evidence_paths": list(health.get("evidence_paths", [])) if isinstance(health.get("evidence_paths"), list) else [],
        "score_status": "missing" if score is None else score_s["score_status"],
        "report_status": report_status,
        "staleness_status": staleness_status,
        "csv_count": int(health.get("csv_count", 0) or 0),
        "csv_ok_count": int(health.get("csv_ok_count", 0) or 0),
        "csv_failed_count": int(health.get("csv_failed_count", 0) or 0),
        "csv_rows": score_s.get("csv_rows", 0),
        "triple_ok_rate": score_s.get("triple_ok_rate"),
        "br_gaga_q95": score_s.get("br_gaga_q95"),
        "br_bb_median": score_s.get("br_bb_median"),
        "total_width_median": score_s.get("total_width_median"),
        "decision_label": None,
        "recommended_action": recommended_action,
        "missing_artifacts": [],
        "limitations": limitations,
        "command_fingerprint": cmd_fingerprint,
    }
    if out["csv_count"] == 0:
        out["missing_artifacts"] = ["scan_csvs"]
    return out


def build_proposal_record(proposal: Mapping[str, object]) -> dict[str, object]:
    metrics = proposal.get("metrics")
    status = str(proposal.get("status") or "")
    lifecycle_hint = "active"
    if status in {"planned", "validated", "submitted", "running"}:
        lifecycle_hint = "in_flight"
    elif status in {"done", "failed", "timeout", "skipped", "rejected"}:
        lifecycle_hint = "terminal"
    return {
        "proposal_id": proposal.get("proposal_id"),
        "parent_id": proposal.get("parent_id"),
        "generation": proposal.get("generation"),
        "operator": proposal.get("operator"),
        "status": proposal.get("status"),
        "bounds": proposal.get("bounds"),
        "fixed": proposal.get("fixed"),
        "scan_values": proposal.get("scan_values"),
        "resolution": proposal.get("resolution"),
        "budget": proposal.get("budget"),
        "run_dir": proposal.get("run_dir"),
        "metrics_keys": sorted(list(metrics.keys())) if isinstance(metrics, Mapping) else [],
        "lifecycle_hint": lifecycle_hint,
    }


def detect_report_staleness(
    report_path: str | Path | None,
    *,
    source_paths: Iterable[str | Path] = (),
) -> dict[str, object]:
    rp = Path(report_path) if report_path is not None else None
    if rp is None or not rp.exists() or not rp.is_file():
        return {
            "report_path": str(rp) if rp is not None else None,
            "exists": False,
            "is_stale": None,
            "stale_reason": "missing_report",
            "report_mtime": None,
            "newest_source_mtime": None,
            "hash_based_staleness_status": "not_implemented",
        }
    report_mtime = rp.stat().st_mtime
    src_mtimes = [Path(s).stat().st_mtime for s in source_paths if Path(s).exists() and Path(s).is_file()]
    if not src_mtimes:
        return {
            "report_path": str(rp),
            "exists": True,
            "is_stale": False,
            "stale_reason": "no_source_paths",
            "report_mtime": report_mtime,
            "newest_source_mtime": None,
            "hash_based_staleness_status": "not_implemented",
        }
    newest = max(src_mtimes)
    stale = report_mtime < newest
    return {
        "report_path": str(rp),
        "exists": True,
        "is_stale": stale,
        "stale_reason": "report_older_than_source" if stale else "fresh",
        "report_mtime": report_mtime,
        "newest_source_mtime": newest,
        "hash_based_staleness_status": "not_implemented",
    }


def _associate_proposals_to_runs(
    proposals: Mapping[str, Mapping[str, object]], run_dirs: list[Path]
) -> tuple[dict[str, str | None], list[str]]:
    warnings: list[str] = []
    mapped: dict[str, str | None] = {}
    run_dirs_s = [str(p) for p in run_dirs]
    assigned_run = set()

    for pid, prop in sorted(proposals.items()):
        rd = prop.get("run_dir") if isinstance(prop, Mapping) else None
        if isinstance(rd, str) and rd.strip() != "":
            normalized = str(Path(rd).resolve()) if Path(rd).exists() else str(Path(rd))
            for r in run_dirs:
                if str(r) == normalized or str(r) == str(Path(rd)):
                    mapped[str(r)] = pid
                    assigned_run.add(str(r))
                    break

    for r in run_dirs_s:
        if r in mapped:
            continue
        found_pid = None
        for pid in sorted(proposals):
            if pid in r:
                found_pid = pid
                break
        mapped[r] = found_pid
        if found_pid is None:
            warnings.append(f"orphan_run_dir:{r}")

    for pid, prop in sorted(proposals.items()):
        pr = prop.get("run_dir") if isinstance(prop, Mapping) else None
        if not pr:
            warnings.append(f"proposal_without_run_dir:{pid}")

    return mapped, warnings


def build_campaign_state(
    *,
    campaign: str,
    registry_path: str | Path | None = None,
    campaign_dir: str | Path | None = None,
    explicit_run_dirs: Iterable[str | Path] | None = None,
    report_paths: Iterable[str | Path] | None = None,
    score_runs: bool = True,
    min_rows_per_csv: int = 1,
    allow_health_warnings: bool = True,
) -> dict[str, object]:
    proposals_state = load_registry_state(registry_path)
    run_dirs = discover_run_dirs(explicit_run_dirs=explicit_run_dirs, campaign_dir=campaign_dir)
    association, assoc_warnings = _associate_proposals_to_runs(proposals_state, run_dirs)

    runs: list[dict[str, object]] = []
    slices: list[dict[str, object]] = []
    warnings = list(assoc_warnings)

    reports: list[dict[str, object]] = []
    for rp in report_paths or ():
        rep = detect_report_staleness(rp, source_paths=[registry_path] if registry_path is not None else [])
        reports.append(rep)
        if rep.get("is_stale") is True:
            warnings.append(f"report_stale:{rep.get('report_path')}")

    for run_dir in run_dirs:
        pid = association.get(str(run_dir))
        prop = proposals_state.get(pid) if pid else None
        health = inspect_run_dir(run_dir, min_rows_per_csv=min_rows_per_csv)
        score = None
        if score_runs and (bool(health.get("is_scoreable")) or (allow_health_warnings and health.get("health_status") == "warning")):
            score = score_run_dir_checked(
                run_dir,
                allow_warnings=allow_health_warnings,
                allow_failed=False,
                min_rows_per_csv=min_rows_per_csv,
            )
        run_record = build_run_record(
            run_dir,
            proposal_id=pid,
            proposal=prop,
            score=score,
            min_rows_per_csv=min_rows_per_csv,
            allow_warnings=allow_health_warnings,
        )
        run_report_candidates = [r for r in reports if isinstance(r.get("report_path"), str) and str(run_dir) in str(r.get("report_path"))]
        if run_report_candidates:
            rr = run_report_candidates[0]
            run_record["report_status"] = "present" if rr.get("exists") else "missing"
            run_record["staleness_status"] = "stale" if rr.get("is_stale") else "fresh"
            if rr.get("is_stale"):
                run_record["recommended_action"] = "regenerate_report"
        elif run_record["recommended_action"] == "ready_for_review":
            run_record["report_status"] = "missing"
            run_record["staleness_status"] = "unknown"

        if run_record["dry_run"]:
            warnings.append(f"dry_run_present:{run_record['run_dir']}")
        if run_record["health_status"] == "failed":
            warnings.append(f"health_failed:{run_record['run_dir']}")

        runs.append(run_record)
        slices.extend(build_slice_records(run_dir, proposal_id=pid))

    proposals = [build_proposal_record(proposals_state[k]) for k in sorted(proposals_state)]

    capabilities = {
        "coverage_status": "not_implemented",
        "novelty_status": "not_implemented",
        "proposal_generation_enabled": False,
        "rerun_planning_status": "not_implemented",
        "registry_runner_status": "not_implemented",
    }
    warnings.extend([
        "coverage_not_implemented",
        "novelty_not_implemented",
        "proposal_generation_disabled",
    ])

    state: dict[str, object] = {
        "schema_version": "1.0",
        "generated_utc": utc_now_iso(),
        "campaign": campaign,
        "inputs": {
            "registry_path": str(registry_path) if registry_path is not None else None,
            "campaign_dir": str(campaign_dir) if campaign_dir is not None else None,
            "explicit_run_dirs": [str(Path(x)) for x in (explicit_run_dirs or [])],
            "report_paths": [str(Path(x)) for x in (report_paths or [])],
            "score_runs": score_runs,
            "min_rows_per_csv": min_rows_per_csv,
            "allow_health_warnings": allow_health_warnings,
        },
        "summary": {},
        "capabilities": capabilities,
        "proposals": proposals,
        "runs": sorted(runs, key=lambda x: str(x.get("run_dir"))),
        "slices": sorted(slices, key=lambda x: str(x.get("source_csv"))),
        "reports": reports,
        "warnings": sorted(set(warnings)),
    }
    state["summary"] = build_campaign_summary(state)
    return state


def build_campaign_summary(state: Mapping[str, object]) -> dict[str, object]:
    proposals = state.get("proposals") if isinstance(state.get("proposals"), list) else []
    runs = state.get("runs") if isinstance(state.get("runs"), list) else []
    slices = state.get("slices") if isinstance(state.get("slices"), list) else []
    reports = state.get("reports") if isinstance(state.get("reports"), list) else []

    def _cnt(items: list[object], pred: Any) -> int:
        return sum(1 for x in items if pred(x))

    n_stale = _cnt(reports, lambda r: isinstance(r, Mapping) and r.get("is_stale") is True)
    return {
        "n_proposals": len(proposals),
        "n_runs_found": len(runs),
        "n_dry_run": _cnt(runs, lambda r: isinstance(r, Mapping) and r.get("dry_run") is True),
        "n_health_ok": _cnt(runs, lambda r: isinstance(r, Mapping) and r.get("health_status") == "ok"),
        "n_health_warning": _cnt(runs, lambda r: isinstance(r, Mapping) and r.get("health_status") == "warning"),
        "n_health_failed": _cnt(runs, lambda r: isinstance(r, Mapping) and r.get("health_status") == "failed"),
        "n_scoreable": _cnt(runs, lambda r: isinstance(r, Mapping) and r.get("health_scoreable") is True),
        "n_scored": _cnt(runs, lambda r: isinstance(r, Mapping) and r.get("score_status") == "scored"),
        "n_reported": _cnt(runs, lambda r: isinstance(r, Mapping) and r.get("report_status") == "present"),
        "n_stale_reports": n_stale,
        "n_missing_run_dir": _cnt(proposals, lambda p: isinstance(p, Mapping) and not p.get("run_dir")),
        "n_needs_rerun": _cnt(runs, lambda r: isinstance(r, Mapping) and r.get("recommended_action") == "needs_rerun_or_inspection"),
        "n_ready_for_review": _cnt(runs, lambda r: isinstance(r, Mapping) and r.get("recommended_action") == "ready_for_review"),
        "n_slices": len(slices),
        "n_slices_scoreable": _cnt(slices, lambda s: isinstance(s, Mapping) and s.get("scoreable") is True),
        "n_slices_failed": _cnt(slices, lambda s: isinstance(s, Mapping) and s.get("health_status") != "ok"),
        "coverage_status": "not_implemented",
        "proposal_generation_status": "disabled",
    }


def write_campaign_state_json(state: Mapping[str, object], path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(dict(state), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_campaign_state_csv(state: Mapping[str, object], path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    runs = state.get("runs") if isinstance(state.get("runs"), list) else []
    with p.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_REQUIRED_CSV_COLUMNS)
        writer.writeheader()
        for run in runs:
            if not isinstance(run, Mapping):
                continue
            row: dict[str, object] = {}
            for col in _REQUIRED_CSV_COLUMNS:
                val = run.get(col)
                if isinstance(val, (dict, list)):
                    row[col] = json.dumps(val, separators=(",", ":"), sort_keys=True)
                else:
                    row[col] = val
            writer.writerow(row)


def write_campaign_state_markdown(state: Mapping[str, object], path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    summary = state.get("summary") if isinstance(state.get("summary"), Mapping) else {}
    caps = state.get("capabilities") if isinstance(state.get("capabilities"), Mapping) else {}
    runs = state.get("runs") if isinstance(state.get("runs"), list) else []
    warnings = state.get("warnings") if isinstance(state.get("warnings"), list) else []

    lines: list[str] = []
    lines.append(f"# Campaign State: {state.get('campaign')}")
    lines.append("")
    lines.append("## Executive status")
    lines.append(DISABLED_CAPABILITY_WARNING)
    lines.append("")
    lines.append("## Summary")
    lines.append("| Metric | Value |")
    lines.append("| --- | --- |")
    for k, v in summary.items():
        lines.append(f"| {k} | {v} |")
    lines.append("")
    lines.append("## Capability status")
    lines.append("| Capability | Status |")
    lines.append("| --- | --- |")
    for k, v in caps.items():
        lines.append(f"| {k} | {v} |")
    lines.append("")
    lines.append("## Runs")
    lines.append("| proposal_id | run_dir | dry_run | health_status | first_fail_event | first_fail_reason | orchestrator_log_signature_count | score_status | report_status | staleness_status | recommended_action |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for run in runs:
        if not isinstance(run, Mapping):
            continue
        lines.append(
            "| "
            + " | ".join(
                str(
                    run.get(k)
                )
                for k in [
                    "proposal_id",
                    "run_dir",
                    "dry_run",
                    "health_status",
                    "first_fail_event",
                    "first_fail_reason",
                    "orchestrator_log_signature_count",
                    "score_status",
                    "report_status",
                    "staleness_status",
                    "recommended_action",
                ]
            )
            + " |"
        )
    lines.append("")
    lines.append("## Warnings")
    for w in warnings:
        lines.append(f"- {w}")
    lines.append("")
    lines.append("## Limitations / non-covered capabilities")
    lines.append("- coverage_not_computed")
    lines.append("- novelty_not_computed")
    lines.append("- proposal_generation_disabled")
    lines.append("- rerun_planning_not_implemented")
    lines.append("")
    lines.append("## Next recommended operational actions")
    actions = sorted({str(run.get("recommended_action")) for run in runs if isinstance(run, Mapping)})
    for action in actions:
        lines.append(f"- {action}")

    p.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _cli(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Build deterministic campaign state")
    ap.add_argument("--campaign", required=True)
    ap.add_argument("--registry", default=None)
    ap.add_argument("--campaign-dir", default=None)
    ap.add_argument("--run-dir", action="append", dest="run_dirs", default=[])
    ap.add_argument("--report", action="append", dest="reports", default=[])
    ap.add_argument("--out-json", default=None)
    ap.add_argument("--out-md", default=None)
    ap.add_argument("--out-csv", default=None)
    ap.add_argument("--no-score", action="store_true")
    args = ap.parse_args(argv)

    if not args.registry and not args.campaign_dir and not args.run_dirs:
        ap.error("At least one of --registry, --campaign-dir, or --run-dir is required")

    state = build_campaign_state(
        campaign=args.campaign,
        registry_path=args.registry,
        campaign_dir=args.campaign_dir,
        explicit_run_dirs=args.run_dirs,
        report_paths=args.reports,
        score_runs=not args.no_score,
    )

    if args.out_json:
        write_campaign_state_json(state, args.out_json)
    if args.out_md:
        write_campaign_state_markdown(state, args.out_md)
    if args.out_csv:
        write_campaign_state_csv(state, args.out_csv)

    summary = state.get("summary", {})
    print(
        json.dumps(
            {
                "campaign": state.get("campaign"),
                "n_runs_found": summary.get("n_runs_found"),
                "n_health_failed": summary.get("n_health_failed"),
                "n_stale_reports": summary.get("n_stale_reports"),
                "n_ready_for_review": summary.get("n_ready_for_review"),
            }
        )
    )

    n_failed = int(summary.get("n_health_failed", 0) or 0)
    n_stale = int(summary.get("n_stale_reports", 0) or 0)
    if n_failed > 0:
        return 2
    if n_stale > 0:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
