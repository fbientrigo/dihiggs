from __future__ import annotations

import argparse
import json
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path

from autoresearch.harness.orchestrator_adapter import (
    build_orchestrator_command,
    normalize_runtime_config,
    validate_executable_proposal,
)

DISABLED_CAPABILITY_WARNING = (
    "Coverage and novelty are not implemented yet. Existing results must not be interpreted as "
    "representative coverage of the full parameter space. Autonomous proposal generation is disabled."
)


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _as_list(value: object) -> list[object]:
    return value if isinstance(value, list) else []


def load_campaign_state(path: str | Path) -> dict[str, object]:
    p = Path(path)
    try:
        payload = json.loads(p.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ValueError(f"failed to load campaign_state from {p}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"campaign_state must be a JSON object: {p}")
    return payload


def _slice_reason_codes(slice_rec: Mapping[str, object]) -> list[str]:
    codes: list[str] = []
    if _as_list(slice_rec.get("missing_columns")):
        codes.append("missing_required_columns")
    if bool(slice_rec.get("is_empty")):
        codes.append("empty_csv")
    if bool(slice_rec.get("scoreable")) is False and not codes:
        codes.append("not_scoreable")
    return sorted(set(codes))


def _run_reason_codes(run_rec: Mapping[str, object], slices_by_run: Mapping[str, list[Mapping[str, object]]]) -> list[str]:
    codes: list[str] = []
    if isinstance(run_rec.get("reason_codes"), list):
        codes.extend([str(x) for x in run_rec.get("reason_codes", [])])
    if str(run_rec.get("health_status")) == "failed":
        codes.append("health_failed")
    if bool(run_rec.get("dry_run")):
        codes.append("dry_run_rejected")
    if int(run_rec.get("csv_count", 0) or 0) == 0:
        codes.append("missing_csv")
    run_dir = str(run_rec.get("run_dir") or "")
    for s in slices_by_run.get(run_dir, []):
        if _as_list(s.get("missing_columns")):
            codes.append("missing_required_columns")
            break
    if not codes:
        codes.append("needs_inspection")
    return sorted(set(codes))


def _proposal_index(state: Mapping[str, object]) -> dict[str, Mapping[str, object]]:
    proposals = [x for x in _as_list(state.get("proposals")) if isinstance(x, Mapping)]
    out: dict[str, Mapping[str, object]] = {}
    for p in proposals:
        pid = p.get("proposal_id")
        if isinstance(pid, str) and pid:
            out[pid] = p
    return out


def _is_lambda6_singleton(payload: Mapping[str, object]) -> bool | None:
    bounds = payload.get("bounds")
    if not isinstance(bounds, Mapping):
        return None
    l6 = bounds.get("lambda6")
    if isinstance(l6, (int, float, str)):
        return True
    if isinstance(l6, list) and len(l6) == 2:
        return str(l6[0]) == str(l6[1])
    return None


def _command_contract_eval(proposal_payload: Mapping[str, object] | None) -> dict[str, object]:
    if proposal_payload is None:
        return {
            "command_contract_ready": "unknown",
            "command_contract_reason_codes": ["proposal_payload_missing"],
            "command_contract_missing_fields": ["proposal_payload"],
            "command_contract_hint": {"deterministic_command_buildable": "unknown"},
            "lambda6_singleton": None,
            "has_scan_values_tan_beta": None,
            "has_runtime_config": None,
        }

    has_scan_values_tan_beta = isinstance(proposal_payload.get("scan_values"), Mapping) and (
        proposal_payload.get("scan_values") or {}
    ).get("tan_beta") is not None
    lambda6_singleton = _is_lambda6_singleton(proposal_payload)
    has_runtime_config = isinstance(proposal_payload.get("runtime_config"), Mapping)

    missing_fields: list[str] = []
    if not has_scan_values_tan_beta:
        missing_fields.append("scan_values.tan_beta")
    if lambda6_singleton is False:
        missing_fields.append("bounds.lambda6_singleton")

    reason_codes: list[str] = []
    try:
        validation_errors = validate_executable_proposal(proposal_payload)
    except Exception as exc:
        validation_errors = [f"adapter_validation_exception:{exc}"]

    if validation_errors:
        reason_codes.append("proposal_contract_invalid")
    if missing_fields:
        reason_codes.append("proposal_fields_missing_or_invalid")

    deterministic_buildable = False
    runtime_cfg = proposal_payload.get("runtime_config") if has_runtime_config else None
    effective_cpp_omp_threads: int | None = None
    try:
        cfg = normalize_runtime_config(runtime_cfg if isinstance(runtime_cfg, Mapping) else None)
        effective_cpp_omp_threads = int(cfg["cpp_omp_threads"])
        _ = build_orchestrator_command(proposal_payload, cfg, dry_run=True, force=False)
        deterministic_buildable = True
    except Exception:
        deterministic_buildable = False
        reason_codes.append("deterministic_command_not_buildable")

    ready = deterministic_buildable and not validation_errors and not missing_fields
    return {
        "command_contract_ready": ready,
        "command_contract_reason_codes": sorted(set(reason_codes)) if reason_codes else ["command_contract_ready"],
        "command_contract_missing_fields": sorted(set(missing_fields)),
        "command_contract_hint": {
            "deterministic_command_buildable": deterministic_buildable,
            "effective_cpp_omp_threads": effective_cpp_omp_threads,
            "thread_policy_note": "C++ subprocess OpenMP threads; not a broad safety or root-cause claim.",
            "validation_errors": validation_errors[:5],
        },
        "lambda6_singleton": lambda6_singleton,
        "has_scan_values_tan_beta": bool(has_scan_values_tan_beta),
        "has_runtime_config": bool(has_runtime_config),
    }


def build_rerun_manifest(
    state: Mapping[str, object],
    *,
    source_campaign_state: str | Path,
) -> dict[str, object]:
    runs = [x for x in _as_list(state.get("runs")) if isinstance(x, Mapping)]
    slices = [x for x in _as_list(state.get("slices")) if isinstance(x, Mapping)]

    slices_by_run: dict[str, list[Mapping[str, object]]] = {}
    for s in slices:
        run_dir = str(s.get("run_dir") or "")
        slices_by_run.setdefault(run_dir, []).append(s)

    run_actions: list[dict[str, object]] = []
    blocked_items: list[dict[str, object]] = []
    proposals_by_id = _proposal_index(state)

    for r in sorted(runs, key=lambda x: str(x.get("run_dir") or "")):
        run_dir = str(r.get("run_dir") or "")
        reason_codes = _run_reason_codes(r, slices_by_run)
        selected = str(r.get("recommended_action")) == "needs_rerun_or_inspection"
        if bool(r.get("dry_run")):
            selected = False
        affected_slices = [str(s.get("source_csv")) for s in slices_by_run.get(run_dir, []) if isinstance(s.get("source_csv"), str)]
        proposal_id = r.get("proposal_id")
        proposal_payload = proposals_by_id.get(str(proposal_id)) if isinstance(proposal_id, str) and proposal_id else None
        contract_eval = _command_contract_eval(proposal_payload)
        action = {
            "run_dir": run_dir,
            "proposal_id": proposal_id,
            "health_status": r.get("health_status"),
            "reason_codes": reason_codes,
            "selected": selected,
            "blocked_reason": None if selected else "selection_policy_excluded_or_dry_run",
            "priority": "high" if selected else "normal",
            "affected_slices": affected_slices,
            "preflight_required": [
                "csv_schema_contract",
                "task_summary_present",
                "run_manifest_present",
            ],
            "recommended_operator_action": "rerun_with_correct_schema" if selected else "inspect",
            "recommended_action": "rerun_with_correct_schema" if selected else "inspect",
            "evidence_paths": list(r.get("evidence_paths", [])) if isinstance(r.get("evidence_paths"), list) else [],
            "command_contract_ready": contract_eval["command_contract_ready"],
            "command_contract_reason_codes": contract_eval["command_contract_reason_codes"],
            "command_contract_missing_fields": contract_eval["command_contract_missing_fields"],
            "command_contract_hint": contract_eval["command_contract_hint"],
            "lambda6_singleton": contract_eval["lambda6_singleton"],
            "has_scan_values_tan_beta": contract_eval["has_scan_values_tan_beta"],
            "has_runtime_config": contract_eval["has_runtime_config"],
        }
        run_actions.append(action)
        if not selected:
            blocked_items.append(
                {
                    "item_type": "run",
                    "run_dir": run_dir,
                    "reason_codes": reason_codes,
                    "blocked_reason": "selection_policy_excluded_or_dry_run",
                    "evidence_paths": list(r.get("evidence_paths", [])) if isinstance(r.get("evidence_paths"), list) else [],
                }
            )

    slice_actions: list[dict[str, object]] = []
    for s in sorted(slices, key=lambda x: str(x.get("source_csv") or "")):
        selected = bool(s.get("scoreable")) is False
        if bool(s.get("dry_run")):
            selected = False
        reason_codes = _slice_reason_codes(s)
        if bool(s.get("dry_run")):
            reason_codes = sorted(set(reason_codes + ["dry_run_rejected"]))
        slice_actions.append(
            {
                "source_csv": s.get("source_csv"),
                "run_dir": s.get("run_dir"),
                "proposal_id": s.get("proposal_id"),
                "scoreable": s.get("scoreable"),
                "reason_codes": reason_codes,
                "selected": selected,
                "blocked_reason": None if selected else "slice_scoreable_or_dry_run",
                "affected_slices": [s.get("source_csv")],
                "evidence_paths": list(s.get("evidence_paths", [])) if isinstance(s.get("evidence_paths"), list) else [s.get("source_csv")],
                "recommended_operator_action": "rerun_slice" if selected else "inspect",
                "recommended_action": "rerun_slice" if selected else "inspect",
            }
        )

    selected_runs = [x for x in run_actions if bool(x.get("selected"))]
    selected_slices = [x for x in slice_actions if bool(x.get("selected"))]

    warnings = [
        "coverage_not_implemented",
        "novelty_not_implemented",
        "proposal_generation_disabled",
    ]

    return {
        "schema_version": "1.0",
        "generated_utc": utc_now_iso(),
        "campaign": state.get("campaign"),
        "source_campaign_state": str(Path(source_campaign_state)),
        "summary": {
            "n_runs_considered": len(runs),
            "n_runs_selected": len(selected_runs),
            "n_slices_selected": len(selected_slices),
            "n_blocked": len(blocked_items),
        },
        "selection_policy": {
            "include_actions": ["needs_rerun_or_inspection"],
            "exclude_dry_run": True,
            "exclude_health_ok": True,
        },
        "run_actions": run_actions,
        "slice_actions": slice_actions,
        "blocked_items": blocked_items,
        "limitations": [
            "coverage_not_computed",
            "novelty_not_computed",
            "proposal_generation_disabled",
            "rerun_execution_not_implemented",
        ],
        "warnings": warnings,
    }


def write_rerun_manifest_json(manifest: Mapping[str, object], path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(dict(manifest), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_rerun_manifest_markdown(manifest: Mapping[str, object], path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)

    summary = manifest.get("summary") if isinstance(manifest.get("summary"), Mapping) else {}
    run_actions = [x for x in _as_list(manifest.get("run_actions")) if isinstance(x, Mapping) and bool(x.get("selected"))]
    slice_actions = [x for x in _as_list(manifest.get("slice_actions")) if isinstance(x, Mapping) and bool(x.get("selected"))]

    lines: list[str] = []
    lines.append(f"# Rerun Manifest: {manifest.get('campaign')}")
    lines.append("")
    lines.append("## Executive status")
    lines.append("Deterministic rerun planning only. No rerun execution performed.")
    lines.append("")
    lines.append("## Summary")
    lines.append("| Metric | Value |")
    lines.append("| --- | --- |")
    for k in ["n_runs_considered", "n_runs_selected", "n_slices_selected", "n_blocked"]:
        lines.append(f"| {k} | {summary.get(k)} |")
    lines.append("")

    lines.append("## Selected runs")
    lines.append("| run_dir | proposal_id | health_status | reason_codes | preflight_required | recommended_operator_action |")
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for r in run_actions:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(r.get("run_dir")),
                    str(r.get("proposal_id")),
                    str(r.get("health_status")),
                    json.dumps(r.get("reason_codes"), separators=(",", ":")),
                    json.dumps(r.get("preflight_required"), separators=(",", ":")),
                    str(r.get("recommended_operator_action")),
                ]
            )
            + " |"
        )
    lines.append("")

    lines.append("## Selected slices")
    lines.append("| source_csv | run_dir | proposal_id | scoreable | reason_codes | recommended_operator_action |")
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for s in slice_actions:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(s.get("source_csv")),
                    str(s.get("run_dir")),
                    str(s.get("proposal_id")),
                    str(s.get("scoreable")),
                    json.dumps(s.get("reason_codes"), separators=(",", ":")),
                    str(s.get("recommended_operator_action")),
                ]
            )
            + " |"
        )
    lines.append("")
    lines.append("## Limitations")
    lines.append("- coverage_not_computed")
    lines.append("- novelty_not_computed")
    lines.append("- proposal_generation_disabled")
    lines.append("- rerun_execution_not_implemented")
    lines.append(f"- {DISABLED_CAPABILITY_WARNING}")
    lines.append("")
    lines.append("## Next operator actions (manual, non-autonomous)")
    lines.append("- Validate CSV schema contract against required physics columns before rerun")
    lines.append("- Regenerate run manifests/task summaries where missing")
    lines.append("- Queue selected runs and selected slices for supervised rerun")
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _cli(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Build deterministic rerun manifest from campaign_state")
    ap.add_argument("--campaign-state", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args(argv)

    state = load_campaign_state(args.campaign_state)
    manifest = build_rerun_manifest(state, source_campaign_state=args.campaign_state)
    write_rerun_manifest_json(manifest, args.out_json)
    write_rerun_manifest_markdown(manifest, args.out_md)

    print(
        json.dumps(
            {
                "campaign": manifest.get("campaign"),
                "n_runs_selected": manifest.get("summary", {}).get("n_runs_selected") if isinstance(manifest.get("summary"), Mapping) else None,
                "n_slices_selected": manifest.get("summary", {}).get("n_slices_selected") if isinstance(manifest.get("summary"), Mapping) else None,
            }
        )
    )

    n_selected = manifest.get("summary", {}).get("n_runs_selected") if isinstance(manifest.get("summary"), Mapping) else 0
    return 0 if int(n_selected or 0) >= 0 else 1


if __name__ == "__main__":
    raise SystemExit(_cli())
