from __future__ import annotations

import argparse
import csv
import json
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from autoresearch.harness.dihiggs_physics_score import score_run_dir_checked
from autoresearch.harness.proposal_registry import load_events, replay_registry


DEFAULT_THRESHOLDS: dict[str, object] = {
    "min_csv_rows_for_rankable": 10,
    "min_triple_ok_rate": 0.01,
    "low_triple_ok_rate": 0.001,
    "high_br_gaga_q95": 0.01,
    "low_br_bb_median": 0.1,
    "suspicious_total_width_low": 1e-15,
    "suspicious_total_width_high": 1e3,
    "outlier_gap_warning": 0.1,
}

DECISION_PRIORITY: dict[str, int] = {
    "exploit": 0,
    "sanity_check": 1,
    "needs_more_points": 2,
    "explore": 3,
    "discard": 4,
    "not_rankable": 5,
}

NEXT_ACTIONS: dict[str, str] = {
    "exploit": "run local refinement around this region",
    "explore": "keep as exploratory candidate",
    "needs_more_points": "rerun or extend statistics before ranking",
    "sanity_check": "inspect widths and boundary behavior before exploiting",
    "discard": "do not allocate more budget unless used as control",
    "not_rankable": "obtain valid scored outputs first",
}


def _merged_thresholds(thresholds: Mapping[str, object] | None) -> dict[str, object]:
    out = dict(DEFAULT_THRESHOLDS)
    if thresholds is not None:
        out.update(dict(thresholds))
    return out


def _safe_float(value: object) -> float | None:
    if isinstance(value, bool) or value is None:
        return None
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _safe_int(value: object) -> int | None:
    parsed = _safe_float(value)
    if parsed is None:
        return None
    return int(parsed)


def _extract_metrics(score: Mapping[str, object] | None) -> dict[str, object]:
    if score is None:
        return {}
    if isinstance(score.get("aggregate_metrics"), Mapping):
        aggregate = dict(score["aggregate_metrics"])  # type: ignore[index]
    else:
        aggregate = dict(score)

    if isinstance(score.get("task_summary_metrics"), Mapping):
        summary = score["task_summary_metrics"]  # type: ignore[index]
        attempts = summary.get("attempts_total")
        triple_attempt = summary.get("triple_ok_rate_over_attempts")
        if attempts is not None:
            aggregate["attempts_total"] = attempts
        if triple_attempt is not None:
            aggregate["triple_ok_rate_over_attempts"] = triple_attempt

    if "rankable" in score:
        aggregate["rankable"] = bool(score.get("rankable"))
    if "run_dir" in score and score.get("run_dir") is not None:
        aggregate["run_dir"] = str(score.get("run_dir"))
    if isinstance(score.get("health"), Mapping):
        aggregate["health"] = dict(score["health"])  # type: ignore[index]
    return aggregate


def classify_decision(metrics: Mapping[str, object], *, thresholds: Mapping[str, object] | None = None) -> dict[str, object]:
    thr = _merged_thresholds(thresholds)

    min_rows = int(thr["min_csv_rows_for_rankable"])
    min_triple = float(thr["min_triple_ok_rate"])
    low_triple = float(thr["low_triple_ok_rate"])
    high_br_gaga = float(thr["high_br_gaga_q95"])
    low_br_bb = float(thr["low_br_bb_median"])
    suspicious_w_low = float(thr["suspicious_total_width_low"])
    suspicious_w_high = float(thr["suspicious_total_width_high"])
    outlier_gap = float(thr["outlier_gap_warning"])

    csv_rows = _safe_int(metrics.get("csv_rows")) or 0
    rankable = bool(metrics.get("rankable"))
    triple_ok_rate = _safe_float(metrics.get("triple_ok_rate"))
    br_gaga_q95 = _safe_float(metrics.get("br_gaga_q95"))
    br_gaga_median = _safe_float(metrics.get("br_gaga_median"))
    br_bb_median = _safe_float(metrics.get("br_bb_median"))
    total_width_median = _safe_float(metrics.get("total_width_median"))

    has_useful_br = br_gaga_q95 is not None or br_gaga_median is not None or br_bb_median is not None

    triggered_rules: list[str] = []
    warnings: list[str] = []

    if br_gaga_q95 is not None and br_gaga_median is not None:
        gap = br_gaga_q95 - br_gaga_median
        if gap > outlier_gap:
            warnings.append(f"large br_gaga tail gap: {gap:.6g}")

    suspicious_width = (
        total_width_median is not None
        and (total_width_median < suspicious_w_low or total_width_median > suspicious_w_high)
    )

    decision_label = "explore"
    decision_reason = "baseline exploratory candidate"

    if not rankable:
        decision_label = "not_rankable"
        decision_reason = "rankable flag is false"
        triggered_rules.append("rule:not_rankable_flag_false")
    elif csv_rows < min_rows:
        if has_useful_br:
            decision_label = "needs_more_points"
            decision_reason = "insufficient CSV rows for stable ranking"
            triggered_rules.append("rule:insufficient_rows_with_metrics")
        else:
            decision_label = "not_rankable"
            decision_reason = "insufficient rows and no useful BR metrics"
            triggered_rules.append("rule:insufficient_rows_no_metrics")
    elif br_gaga_q95 is not None and br_gaga_q95 >= high_br_gaga and suspicious_width:
        decision_label = "sanity_check"
        decision_reason = "high br_gaga with suspicious total_width scale"
        triggered_rules.append("rule:high_br_gaga_suspicious_width")
    elif (
        br_gaga_q95 is not None
        and br_gaga_q95 >= high_br_gaga
        and triple_ok_rate is not None
        and triple_ok_rate < low_triple
    ):
        decision_label = "sanity_check"
        decision_reason = "high br_gaga but extremely low triple_ok_rate"
        triggered_rules.append("rule:high_br_gaga_very_low_triple_ok")
    elif (
        br_gaga_q95 is not None
        and br_gaga_q95 >= high_br_gaga
        and br_bb_median is not None
        and br_bb_median < low_br_bb
        and triple_ok_rate is not None
        and triple_ok_rate >= min_triple
    ):
        decision_label = "exploit"
        decision_reason = "high br_gaga, low br_bb, and acceptable triple_ok_rate"
        triggered_rules.append("rule:high_br_gaga_low_br_bb_feasible")
    elif (
        br_gaga_q95 is not None
        and br_gaga_q95 < high_br_gaga
        and br_bb_median is not None
        and br_bb_median >= low_br_bb
        and triple_ok_rate is not None
        and triple_ok_rate < low_triple
    ):
        decision_label = "discard"
        decision_reason = "low br_gaga, high br_bb, and low triple_ok_rate"
        triggered_rules.append("rule:weak_signal_low_feasibility")
    else:
        decision_label = "explore"
        decision_reason = "metrics are mixed without strong exploit/discard signal"
        triggered_rules.append("rule:mixed_metrics_explore")

    if not triggered_rules:
        triggered_rules.append("rule:fallback")

    return {
        "decision_label": decision_label,
        "decision_reason": decision_reason,
        "triggered_rules": triggered_rules,
        "warnings": warnings,
    }


def compute_derived_decision_metrics(score: Mapping[str, object]) -> dict[str, object]:
    metrics = _extract_metrics(score)

    feasibility_score = _safe_float(metrics.get("triple_ok_rate_over_attempts"))
    if feasibility_score is None:
        feasibility_score = _safe_float(metrics.get("triple_ok_rate_over_csv_rows"))

    br_gaga_q95 = _safe_float(metrics.get("br_gaga_q95"))
    br_gaga_median = _safe_float(metrics.get("br_gaga_median"))
    br_bb_median = _safe_float(metrics.get("br_bb_median"))

    if br_gaga_q95 is not None and br_bb_median is not None:
        physics_score = br_gaga_q95 - 0.1 * br_bb_median
    else:
        physics_score = None

    if br_gaga_q95 is not None and br_gaga_median is not None:
        robustness_gap = br_gaga_q95 - br_gaga_median
    else:
        robustness_gap = None

    attempts_total = _safe_float(metrics.get("attempts_total"))
    csv_rows = _safe_float(metrics.get("csv_rows"))
    if attempts_total is not None and attempts_total > 0 and csv_rows is not None:
        efficiency_score = csv_rows / attempts_total
    else:
        efficiency_score = None

    return {
        "feasibility_score": feasibility_score,
        "physics_score": physics_score,
        "robustness_gap_br_gaga": robustness_gap,
        "efficiency_score": efficiency_score,
        "boundary_warning": None,
        "novelty_score": None,
        "coverage_warning": "not_computed_yet",
    }


def summarize_proposal_state(proposal: Mapping[str, object]) -> dict[str, object]:
    return {
        "proposal_id": proposal.get("proposal_id"),
        "parent_id": proposal.get("parent_id"),
        "generation": proposal.get("generation"),
        "operator": proposal.get("operator"),
        "status": proposal.get("status"),
        "bounds": proposal.get("bounds"),
        "fixed": proposal.get("fixed"),
        "resolution": proposal.get("resolution"),
        "budget": proposal.get("budget"),
        "run_dir": proposal.get("run_dir"),
        "metrics": proposal.get("metrics"),
    }


def build_decision_row(
    proposal: Mapping[str, object],
    score: Mapping[str, object] | None = None,
    thresholds: Mapping[str, object] | None = None,
) -> dict[str, object]:
    base = summarize_proposal_state(proposal)
    metric_source = _extract_metrics(score)
    derived = compute_derived_decision_metrics(metric_source)

    triple_ok_rate = _safe_float(metric_source.get("triple_ok_rate_over_attempts"))
    if triple_ok_rate is None:
        triple_ok_rate = _safe_float(metric_source.get("triple_ok_rate_over_csv_rows"))

    decision_input = {
        "rankable": bool(metric_source.get("rankable")) if metric_source else False,
        "csv_rows": metric_source.get("csv_rows"),
        "triple_ok_rate": triple_ok_rate,
        "br_gaga_q95": metric_source.get("br_gaga_q95"),
        "br_gaga_median": metric_source.get("br_gaga_median"),
        "br_bb_median": metric_source.get("br_bb_median"),
        "total_width_median": metric_source.get("total_width_median"),
    }
    decision = classify_decision(decision_input, thresholds=thresholds)

    health = metric_source.get("health") if isinstance(metric_source.get("health"), Mapping) else None
    health_status = health.get("health_status") if isinstance(health, Mapping) else None
    health_scoreable = health.get("is_scoreable") if isinstance(health, Mapping) else None
    health_recommended_action = health.get("recommended_action") if isinstance(health, Mapping) else None
    dry_run = health.get("dry_run") if isinstance(health, Mapping) else None

    if isinstance(health, Mapping) and health_scoreable is False:
        rules = list(decision.get("triggered_rules", []))
        rules.append("rule:health_not_scoreable")
        decision["triggered_rules"] = rules
        warnings = list(decision.get("warnings", []))
        warnings.append("run_dir_not_scoreable")
        decision["warnings"] = warnings
        decision["decision_label"] = "not_rankable"
        decision["decision_reason"] = "health gate marked run_dir as not scoreable"

    budget = proposal.get("budget") if isinstance(proposal.get("budget"), Mapping) else {}
    run_dir = metric_source.get("run_dir") if metric_source.get("run_dir") is not None else proposal.get("run_dir")

    return {
        "proposal_id": base.get("proposal_id"),
        "parent_id": base.get("parent_id"),
        "generation": base.get("generation"),
        "operator": base.get("operator"),
        "status": base.get("status"),
        "run_dir": run_dir,
        "bounds": base.get("bounds"),
        "requested_points": budget.get("requested_points") if isinstance(budget, Mapping) else None,
        "max_walltime_sec": budget.get("max_walltime_sec") if isinstance(budget, Mapping) else None,
        "rankable": bool(metric_source.get("rankable")) if metric_source else False,
        "csv_rows": metric_source.get("csv_rows"),
        "triple_ok_rate": triple_ok_rate,
        "br_gaga_q95": metric_source.get("br_gaga_q95"),
        "br_gaga_median": metric_source.get("br_gaga_median"),
        "br_bb_median": metric_source.get("br_bb_median"),
        "total_width_median": metric_source.get("total_width_median"),
        "objective_score": metric_source.get("objective_score"),
        "feasibility_score": derived["feasibility_score"],
        "physics_score": derived["physics_score"],
        "robustness_gap_br_gaga": derived["robustness_gap_br_gaga"],
        "efficiency_score": derived["efficiency_score"],
        "novelty_score": derived["novelty_score"],
        "coverage_warning": derived["coverage_warning"],
        "health_status": health_status,
        "health_scoreable": health_scoreable,
        "health_recommended_action": health_recommended_action,
        "dry_run": dry_run,
        "decision_label": decision["decision_label"],
        "decision_reason": decision["decision_reason"],
        "triggered_rules": decision["triggered_rules"],
        "warnings": decision["warnings"],
        "recommended_next_action": NEXT_ACTIONS[str(decision["decision_label"])],
    }


def load_registry_state(registry_path: str | Path, contract: Mapping[str, object] | None = None) -> dict[str, dict]:
    events = load_events(registry_path)
    return replay_registry(events, contract=contract)


def map_scores_by_run_dir(score_inputs: Iterable[Mapping[str, object]]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for item in score_inputs:
        run_dir = item.get("run_dir")
        if run_dir is None:
            continue
        normalized = str(Path(str(run_dir)))
        out[normalized] = dict(item)
    return out


def _priority_for_label(label: object) -> int:
    return DECISION_PRIORITY.get(str(label), 999)


def _score_for_sort(value: object) -> float:
    parsed = _safe_float(value)
    if parsed is None:
        return float("-inf")
    return parsed


def build_decision_report(
    proposals: Mapping[str, Mapping[str, object]],
    scores_by_run_dir: Mapping[str, Mapping[str, object]],
    thresholds: Mapping[str, object] | None = None,
) -> dict[str, object]:
    rows: list[dict[str, object]] = []

    for proposal_id in sorted(proposals):
        proposal = proposals[proposal_id]
        run_dir_raw = proposal.get("run_dir")
        score = None
        if run_dir_raw is not None:
            score = scores_by_run_dir.get(str(Path(str(run_dir_raw))))
        rows.append(build_decision_row(proposal, score=score, thresholds=thresholds))

    rows.sort(
        key=lambda row: (
            _priority_for_label(row.get("decision_label")),
            -_score_for_sort(row.get("physics_score")),
            str(row.get("proposal_id") or ""),
        )
    )

    n_by_label = {label: 0 for label in DECISION_PRIORITY}
    n_score_rankable = 0
    n_decision_rankable = 0
    n_health_ok = 0
    n_health_warning = 0
    n_health_failed = 0
    n_health_scoreable = 0
    n_health_dry_run = 0
    summary_warnings: list[str] = []
    best_physics_id = None
    best_physics = None
    best_objective_id = None
    best_objective = None

    for row in rows:
        label = str(row.get("decision_label"))
        if label in n_by_label:
            n_by_label[label] += 1
        if bool(row.get("rankable")):
            n_score_rankable += 1
        if label != "not_rankable":
            n_decision_rankable += 1

        hstatus = row.get("health_status")
        if hstatus == "ok":
            n_health_ok += 1
        elif hstatus == "warning":
            n_health_warning += 1
        elif hstatus == "failed":
            n_health_failed += 1
        if row.get("health_scoreable") is True:
            n_health_scoreable += 1
        if row.get("dry_run") is True:
            n_health_dry_run += 1

        physics = _safe_float(row.get("physics_score"))
        if physics is not None and (best_physics is None or physics > best_physics):
            best_physics = physics
            best_physics_id = row.get("proposal_id")

        objective = _safe_float(row.get("objective_score"))
        if objective is not None and (best_objective is None or objective > best_objective):
            best_objective = objective
            best_objective_id = row.get("proposal_id")

    if not rows:
        summary_warnings.append("empty proposal set")

    summary = {
        "n_proposals": len(rows),
        "n_rankable": n_decision_rankable,
        "n_score_rankable": n_score_rankable,
        "n_decision_rankable": n_decision_rankable,
        "n_exploit": n_by_label["exploit"],
        "n_explore": n_by_label["explore"],
        "n_needs_more_points": n_by_label["needs_more_points"],
        "n_sanity_check": n_by_label["sanity_check"],
        "n_discard": n_by_label["discard"],
        "n_not_rankable": n_by_label["not_rankable"],
        "n_health_ok": n_health_ok,
        "n_health_warning": n_health_warning,
        "n_health_failed": n_health_failed,
        "n_health_scoreable": n_health_scoreable,
        "n_health_dry_run": n_health_dry_run,
        "best_proposal_id_by_physics_score": best_physics_id,
        "best_proposal_id_by_objective_score": best_objective_id,
        "warnings": summary_warnings,
    }

    return {
        "schema_version": "phase4_decision_report_v1",
        "summary": summary,
        "rows": rows,
    }


def write_decision_report_json(report: Mapping[str, object], path: str | Path) -> None:
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(dict(report), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_decision_report_markdown(report: Mapping[str, object], path: str | Path) -> None:
    summary = report.get("summary") if isinstance(report.get("summary"), Mapping) else {}
    rows = report.get("rows") if isinstance(report.get("rows"), list) else []

    def _table(headers: list[str], body_rows: list[list[object]]) -> str:
        lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
        for row in body_rows:
            lines.append("| " + " | ".join(str(x) for x in row) + " |")
        return "\n".join(lines)

    top_candidates = [r for r in rows if isinstance(r, Mapping) and r.get("decision_label") == "exploit"][:10]
    sanity_rows = [r for r in rows if isinstance(r, Mapping) and r.get("decision_label") == "sanity_check"][:10]
    dropped = [r for r in rows if isinstance(r, Mapping) and r.get("decision_label") in {"discard", "not_rankable"}]

    summary_headers = ["Metric", "Value"]
    summary_rows = [[k, summary.get(k)] for k in [
        "n_proposals",
        "n_rankable",
        "n_exploit",
        "n_explore",
        "n_needs_more_points",
        "n_sanity_check",
        "n_discard",
        "n_not_rankable",
        "best_proposal_id_by_physics_score",
        "best_proposal_id_by_objective_score",
    ]]

    top_headers = ["proposal_id", "physics_score", "objective_score", "triple_ok_rate", "br_gaga_q95", "decision_reason"]
    top_body = [
        [
            r.get("proposal_id"),
            r.get("physics_score"),
            r.get("objective_score"),
            r.get("triple_ok_rate"),
            r.get("br_gaga_q95"),
            r.get("decision_reason"),
        ]
        for r in top_candidates
    ]

    sanity_headers = ["proposal_id", "total_width_median", "br_gaga_q95", "triple_ok_rate", "decision_reason"]
    sanity_body = [
        [
            r.get("proposal_id"),
            r.get("total_width_median"),
            r.get("br_gaga_q95"),
            r.get("triple_ok_rate"),
            r.get("decision_reason"),
        ]
        for r in sanity_rows
    ]

    dropped_headers = ["proposal_id", "decision_label", "decision_reason", "triggered_rules"]
    dropped_body = [
        [r.get("proposal_id"), r.get("decision_label"), r.get("decision_reason"), json.dumps(r.get("triggered_rules"), separators=(",", ":"))]
        for r in dropped
    ]

    sections = [
        "# DiHiggs Autoresearch Decision Report",
        "",
        "## Summary",
        _table(summary_headers, summary_rows),
        "",
        "## Top Candidates",
        _table(top_headers, top_body),
        "",
        "## Sanity Checks",
        _table(sanity_headers, sanity_body),
        "",
        "## Discarded and Not Rankable",
        _table(dropped_headers, dropped_body),
        "",
        "## Notes",
        "- novelty_score is currently a placeholder (None).",
        "- coverage_warning is currently set to not_computed_yet.",
        "- Decisions are deterministic rule outputs from reported metrics.",
    ]

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(sections) + "\n", encoding="utf-8")


def write_decision_rows_csv(report: Mapping[str, object], path: str | Path) -> None:
    rows = report.get("rows") if isinstance(report.get("rows"), list) else []
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        out_path.write_text("", encoding="utf-8")
        return

    fieldnames: list[str] = []
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(str(key))

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            if not isinstance(row, Mapping):
                continue
            serialized: dict[str, object] = {}
            for key in fieldnames:
                value = row.get(key)
                if isinstance(value, (dict, list)):
                    serialized[key] = json.dumps(value, sort_keys=True, separators=(",", ":"))
                else:
                    serialized[key] = value
            writer.writerow(serialized)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build deterministic DiHiggs decision reports")
    parser.add_argument("--registry", required=True, help="Path to proposal_registry.jsonl")
    parser.add_argument("--run-dir", action="append", default=[], help="Run directory to score (repeatable)")
    parser.add_argument("--out-json", default=None)
    parser.add_argument("--out-md", default=None)
    parser.add_argument("--out-csv", default=None)
    parser.add_argument(
        "--allow-health-warning",
        choices=["true", "false"],
        default="true",
        help="Allow run_dirs with health_status=warning to be scored/rankable",
    )
    parser.add_argument(
        "--allow-health-failed",
        choices=["true", "false"],
        default="false",
        help="Allow non-scoreable/failed run_dirs to proceed to raw scoring",
    )
    parser.add_argument(
        "--min-rows-per-csv",
        type=int,
        default=1,
        help="Minimum required rows per CSV for run_health inspect",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    proposals = load_registry_state(args.registry)
    allow_health_warning = str(args.allow_health_warning).lower() != "false"
    allow_health_failed = str(args.allow_health_failed).lower() == "true"
    score_inputs = [
        score_run_dir_checked(
            run_dir,
            allow_warnings=allow_health_warning,
            allow_failed=allow_health_failed,
            min_rows_per_csv=max(1, int(args.min_rows_per_csv)),
        )
        for run_dir in args.run_dir
    ]
    scores_by_run_dir = map_scores_by_run_dir(score_inputs)
    report = build_decision_report(proposals, scores_by_run_dir)

    if args.out_json:
        write_decision_report_json(report, args.out_json)
    if args.out_md:
        write_decision_report_markdown(report, args.out_md)
    if args.out_csv:
        write_decision_rows_csv(report, args.out_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
