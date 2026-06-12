from __future__ import annotations

import csv
import json
import math
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

from autoresearch.harness.dihiggs_validators import compute_triple_ok, derive_br_bb
from autoresearch.harness.run_health import inspect_run_dir


def read_physics_csv(path: str | Path) -> list[dict[str, object]]:
    csv_path = Path(path)
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return []
        return [dict(row) for row in reader]


def finite_float(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    if value is None:
        return None
    if isinstance(value, str) and value.strip() == "":
        return None
    if isinstance(value, (int, float, str)):
        try:
            out = float(value)
        except (TypeError, ValueError):
            return None
        if math.isfinite(out):
            return out
    return None


def quantile(values: Sequence[float], q: float) -> float | None:
    if not 0.0 <= q <= 1.0:
        raise ValueError("q must be within [0, 1]")
    if len(values) == 0:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    pos = (len(ordered) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return ordered[lo]
    frac = pos - lo
    return ordered[lo] + frac * (ordered[hi] - ordered[lo])


def _base_row_metrics() -> dict[str, object]:
    return {
        "csv_rows": 0,
        "triple_ok_count": 0,
        "triple_ok_rate_over_csv_rows": None,
        "br_gaga_count": 0,
        "br_gaga_max": None,
        "br_gaga_q95": None,
        "br_gaga_median": None,
        "br_bb_count": 0,
        "br_bb_min": None,
        "br_bb_q05": None,
        "br_bb_median": None,
        "total_width_count": 0,
        "total_width_median": None,
        "rankable": False,
        "objective_score": None,
        "warnings": [],
    }


def score_rows(rows: Iterable[Mapping[str, object]]) -> dict[str, object]:
    out = _base_row_metrics()
    rows_list = [dict(r) for r in rows]
    csv_rows = len(rows_list)
    out["csv_rows"] = csv_rows

    triple_ok_count = sum(1 for row in rows_list if compute_triple_ok(row))
    out["triple_ok_count"] = triple_ok_count
    out["triple_ok_rate_over_csv_rows"] = (triple_ok_count / csv_rows) if csv_rows > 0 else None

    br_gaga_values: list[float] = []
    br_bb_values: list[float] = []
    total_width_values: list[float] = []

    for row in rows_list:
        br_gaga = finite_float(row.get("br_gaga"))
        if br_gaga is not None:
            br_gaga_values.append(br_gaga)

        br_bb = derive_br_bb(row)
        if br_bb is not None:
            br_bb_values.append(br_bb)

        total_width = finite_float(row.get("total_width"))
        if total_width is not None and total_width > 0.0:
            total_width_values.append(total_width)

    out["br_gaga_count"] = len(br_gaga_values)
    out["br_gaga_max"] = max(br_gaga_values) if br_gaga_values else None
    out["br_gaga_q95"] = quantile(br_gaga_values, 0.95)
    out["br_gaga_median"] = quantile(br_gaga_values, 0.5)

    out["br_bb_count"] = len(br_bb_values)
    out["br_bb_min"] = min(br_bb_values) if br_bb_values else None
    out["br_bb_q05"] = quantile(br_bb_values, 0.05)
    out["br_bb_median"] = quantile(br_bb_values, 0.5)

    out["total_width_count"] = len(total_width_values)
    out["total_width_median"] = quantile(total_width_values, 0.5)

    rankable = csv_rows > 0 and len(br_gaga_values) > 0 and len(br_bb_values) > 0
    out["rankable"] = rankable

    br_gaga_q95 = out["br_gaga_q95"]
    br_bb_median = out["br_bb_median"]
    if br_gaga_q95 is not None and br_bb_median is not None:
        out["objective_score"] = br_gaga_q95 - 0.1 * br_bb_median

    warnings = out["warnings"]
    assert isinstance(warnings, list)
    if csv_rows == 0:
        warnings.append("no CSV rows available")
    if not br_gaga_values:
        warnings.append("no finite br_gaga values")
    if not br_bb_values:
        warnings.append("no derivable br_bb values")
    if not total_width_values:
        warnings.append("no finite positive total_width values")
    if out["objective_score"] is None:
        warnings.append("objective_score unavailable due to missing quantiles")

    return out


def score_csv(path: str | Path) -> dict[str, object]:
    csv_path = Path(path)
    rows = read_physics_csv(csv_path)
    out = score_rows(rows)
    out["csv_path"] = str(csv_path)
    return out


def _as_int_like(value: object) -> int:
    parsed = finite_float(value)
    if parsed is None:
        return 0
    return int(parsed)


# Status / event vocabularies accepted from real orchestrator task summaries.
# The orchestrator emits status="done"|"fail"|"crash"|"skip"|"dry_run", while
# older/aggregated summaries used "failed"/"skipped"/"timeout".  Accept both so
# no real failure or skip event is silently dropped.
_DONE_STATUSES = frozenset({"done", "success", "succeeded", "complete", "completed", "ok"})
_FAIL_STATUSES = frozenset({"fail", "failed", "crash", "crashed", "error"})
_SKIP_STATUSES = frozenset({"skip", "skipped"})
_TIMEOUT_STATUSES = frozenset({"timeout", "timed_out", "timedout"})


def _summary_attempts(event: Mapping[str, object]) -> int:
    """Read the attempt count using the documented fallback order.

    Order: ``total_attempts`` → ``attempts`` → ``metrics.total_attempts`` →
    ``metrics.attempts``.
    """
    for key in ("total_attempts", "attempts"):
        if event.get(key) is not None:
            return _as_int_like(event.get(key))
    metrics = event.get("metrics")
    if isinstance(metrics, Mapping):
        for key in ("total_attempts", "attempts"):
            if metrics.get(key) is not None:
                return _as_int_like(metrics.get(key))
    return 0


def _summary_triple_ok(event: Mapping[str, object]) -> int:
    triple_ok_points = event.get("triple_ok_points")
    if triple_ok_points is None and isinstance(event.get("metrics"), Mapping):
        metrics = event["metrics"]
        assert isinstance(metrics, Mapping)
        triple_ok_points = metrics.get("triple_ok_points")
    return _as_int_like(triple_ok_points)


def _summary_return_code(event: Mapping[str, object]) -> int | None:
    for key in ("returncode", "return_code", "rc", "exit_code"):
        if event.get(key) is not None:
            parsed = finite_float(event.get(key))
            if parsed is not None:
                return int(parsed)
    return None


def _new_task_summary_metrics() -> dict[str, object]:
    return {
        "tasks_total": 0,
        "tasks_done": 0,
        "tasks_failed": 0,
        "tasks_skipped": 0,
        "tasks_timeout": 0,
        "attempts_total": 0,
        "triple_ok_total": 0,
        "triple_ok_rate_over_attempts": None,
        "warnings": [],
    }


def score_task_summary(path: str | Path) -> dict[str, object]:
    summary_path = Path(path)
    out = _new_task_summary_metrics()
    warnings = out["warnings"]
    assert isinstance(warnings, list)

    if not summary_path.exists():
        warnings.append(f"missing task summary file: {summary_path}")
        return out

    with summary_path.open("r", encoding="utf-8") as handle:
        for idx, line in enumerate(handle, start=1):
            text = line.strip()
            if text == "":
                continue
            try:
                event = json.loads(text)
            except json.JSONDecodeError as exc:
                raise ValueError(f"invalid JSON on line {idx}") from exc

            if not isinstance(event, Mapping):
                out["tasks_total"] = int(out["tasks_total"]) + 1
                continue

            out["tasks_total"] = int(out["tasks_total"]) + 1

            status_raw = event.get(
                "status",
                event.get("task_status", event.get("event", event.get("outcome", event.get("result")))),
            )
            status = str(status_raw).strip().lower() if status_raw is not None else ""
            return_code = _summary_return_code(event)
            if status in _DONE_STATUSES:
                out["tasks_done"] = int(out["tasks_done"]) + 1
            elif status in _SKIP_STATUSES:
                out["tasks_skipped"] = int(out["tasks_skipped"]) + 1
            elif status in _TIMEOUT_STATUSES:
                out["tasks_timeout"] = int(out["tasks_timeout"]) + 1
            elif status in _FAIL_STATUSES:
                out["tasks_failed"] = int(out["tasks_failed"]) + 1
            elif return_code is not None and return_code != 0:
                # No recognized status but the runner reported a nonzero exit:
                # treat it as a failure rather than dropping it.
                out["tasks_failed"] = int(out["tasks_failed"]) + 1

            out["attempts_total"] = int(out["attempts_total"]) + _summary_attempts(event)
            out["triple_ok_total"] = int(out["triple_ok_total"]) + _summary_triple_ok(event)

    attempts_total = int(out["attempts_total"])
    triple_ok_total = int(out["triple_ok_total"])
    out["triple_ok_rate_over_attempts"] = (triple_ok_total / attempts_total) if attempts_total > 0 else None
    return out


def score_run_dir(run_dir: str | Path) -> dict[str, object]:
    run_path = Path(run_dir)
    summary_metrics = score_task_summary(run_path / "task_summary.jsonl")

    csv_paths = sorted(run_path.glob("tb_*/scan_tb_*.csv"))
    csv_metrics = [score_csv(path) for path in csv_paths]

    all_rows: list[dict[str, object]] = []
    for path in csv_paths:
        all_rows.extend(read_physics_csv(path))

    aggregate = score_rows(all_rows)

    return {
        "run_dir": str(run_path),
        "task_summary_metrics": summary_metrics,
        "csv_metrics": csv_metrics,
        "aggregate_metrics": aggregate,
        "rankable": bool(aggregate.get("rankable")),
    }


def score_run_dir_checked(
    run_dir: str | Path,
    *,
    allow_warnings: bool = True,
    allow_failed: bool = False,
    min_rows_per_csv: int = 1,
) -> dict[str, object]:
    run_path = Path(run_dir)
    health = inspect_run_dir(run_path, min_rows_per_csv=min_rows_per_csv)

    if not bool(health.get("is_scoreable")) and not allow_failed:
        aggregate = _base_row_metrics()
        warnings = list(aggregate.get("warnings", []))
        warnings.append("run_dir_not_scoreable")
        aggregate["warnings"] = warnings
        return {
            "run_dir": str(run_path),
            "task_summary_metrics": {},
            "csv_metrics": [],
            "aggregate_metrics": aggregate,
            "rankable": False,
            "health": health,
        }

    if str(health.get("health_status")) == "warning" and not allow_warnings:
        aggregate = _base_row_metrics()
        warnings = list(aggregate.get("warnings", []))
        warnings.append("run_dir_health_warning_blocked")
        aggregate["warnings"] = warnings
        return {
            "run_dir": str(run_path),
            "task_summary_metrics": {},
            "csv_metrics": [],
            "aggregate_metrics": aggregate,
            "rankable": False,
            "health": health,
        }

    scored = score_run_dir(run_path)
    scored["health"] = health
    return scored
