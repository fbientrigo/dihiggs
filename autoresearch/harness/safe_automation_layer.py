from __future__ import annotations

import argparse
import csv
import json
import subprocess
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from autoresearch.harness.campaign_state import build_campaign_state, write_campaign_state_json
from autoresearch.harness.dihiggs_physics_score import finite_float, quantile
from autoresearch.harness.dihiggs_validators import compute_triple_ok, derive_br_bb
from autoresearch.harness.orchestrator_adapter import build_orchestrator_command
from autoresearch.harness.rerun_manifest import build_rerun_manifest, write_rerun_manifest_json

GSL_ABORT_SIGNATURES = ("gsl", "sigabrt", "returncode=-6", "abort")
OMP_AUDIT_LINE = "[CONF] OMP_NUM_THREADS = 1"


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _load_yaml(path: Path) -> dict[str, Any]:
    try:
        import yaml  # type: ignore
    except Exception as exc:
        raise ValueError("YAML contract requires PyYAML") from exc
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("contract must deserialize to a mapping")
    return payload


def load_contract(path: str | Path) -> dict[str, Any]:
    p = Path(path)
    if not p.exists():
        raise ValueError(f"contract not found: {p}")
    if p.suffix.lower() in {".yaml", ".yml"}:
        return _load_yaml(p)
    payload = json.loads(p.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("contract must deserialize to a mapping")
    return payload


def _require_bool(val: object, key: str) -> bool:
    if not isinstance(val, bool):
        raise ValueError(f"{key} must be boolean")
    return val


def _require_str(val: object, key: str) -> str:
    if not isinstance(val, str) or not val.strip():
        raise ValueError(f"{key} must be non-empty string")
    return val.strip()


def _require_int(val: object, key: str, *, min_value: int = 0) -> int:
    if not isinstance(val, int):
        raise ValueError(f"{key} must be int")
    if val < min_value:
        raise ValueError(f"{key} must be >= {min_value}")
    return val


def _require_mapping(val: object, key: str) -> dict[str, Any]:
    if not isinstance(val, Mapping):
        raise ValueError(f"{key} must be mapping")
    return dict(val)


def _require_sequence(val: object, key: str) -> list[Any]:
    if not isinstance(val, Sequence) or isinstance(val, (str, bytes)):
        raise ValueError(f"{key} must be sequence")
    return list(val)


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    c = dict(contract)
    required_false = [
        "enable_proposal_generation",
        "allow_broad_scans",
        "modify_physics_semantics",
        "modify_scoring_semantics",
        "modify_triple_ok_definition",
        "use_all_row_br_metrics_as_physics_claims",
        "investigate_gsl_root_cause",
    ]
    for k in required_false:
        if _require_bool(c.get(k), k):
            raise ValueError(f"{k} must be false")

    if not _require_bool(c.get("require_cpp_omp_threads_1"), "require_cpp_omp_threads_1"):
        raise ValueError("require_cpp_omp_threads_1 must be true")

    analysis = _require_mapping(c.get("analysis"), "analysis")
    if not _require_bool(analysis.get("triple_ok_only"), "analysis.triple_ok_only"):
        raise ValueError("analysis.triple_ok_only must be true")

    gates = _require_mapping(c.get("gates"), "gates")
    _ = float(gates.get("min_triple_ok_rate", 0.0))

    evidence_mode = str(c.get("evidence_mode", "production")).strip().lower()
    if evidence_mode not in {"production", "existing_artifacts", "fixture", "no_execute"}:
        raise ValueError("evidence_mode must be one of production|existing_artifacts|fixture|no_execute")

    runtime = _require_mapping(c.get("runtime"), "runtime")
    runtime_threads = runtime.get("threads", runtime.get("cpp_omp_threads"))
    if _require_int(runtime_threads, "runtime.threads", min_value=1) != 1:
        raise ValueError("runtime.threads must be 1")

    campaign = _require_mapping(c.get("campaign"), "campaign")
    _require_str(campaign.get("name"), "campaign.name")
    _require_int(campaign.get("max_runs"), "campaign.max_runs", min_value=1)
    if _require_int(campaign.get("expansion_max_runs"), "campaign.expansion_max_runs", min_value=0) != 0:
        raise ValueError("campaign.expansion_max_runs must be 0 in v1.1")

    proposals = _require_sequence(campaign.get("proposals"), "campaign.proposals")
    if not proposals:
        raise ValueError("campaign.proposals must be non-empty")

    out = dict(c)
    out["campaign"] = campaign
    out["runtime"] = runtime
    out["analysis"] = analysis
    out["gates"] = gates
    return out


def build_commands(contract: Mapping[str, Any]) -> list[list[str]]:
    campaign = _require_mapping(contract.get("campaign"), "campaign")
    runtime = _require_mapping(contract.get("runtime"), "runtime")
    proposals = _require_sequence(campaign.get("proposals"), "campaign.proposals")

    commands: list[list[str]] = []
    for item in proposals:
        if not isinstance(item, Mapping):
            raise ValueError("each campaign.proposals item must be mapping")
        proposal = dict(item)
        rcfg = dict(proposal.get("runtime_config", {})) if isinstance(proposal.get("runtime_config"), Mapping) else {}
        rcfg["cpp_omp_threads"] = 1
        rcfg["threads"] = 1
        rcfg["campaign"] = campaign["name"]
        rcfg["exec_path"] = runtime["exec_path"]
        rcfg["outdir"] = runtime["outdir"]
        rcfg["lake_name"] = runtime["lake_name"]
        cmd = build_orchestrator_command(proposal, rcfg, dry_run=bool(runtime.get("dry_run", False)), force=bool(runtime.get("force", False)))
        if "--threads" not in cmd or cmd[cmd.index("--threads") + 1] != "1":
            raise ValueError("constructed command violates --threads 1")
        commands.append(cmd)
    return commands


def _run_command(cmd: Sequence[str], workdir: Path) -> dict[str, Any]:
    start = datetime.now(timezone.utc)
    proc = subprocess.run(list(cmd), cwd=str(workdir), text=True, capture_output=True)
    end = datetime.now(timezone.utc)
    return {
        "command": list(cmd),
        "returncode": int(proc.returncode),
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "runtime_seconds": (end - start).total_seconds(),
    }


def _scan_csv_rows(run_dir: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for csv_path in sorted(run_dir.glob("tb_*/scan_tb_*.csv")):
        with csv_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for r in reader:
                rr = dict(r)
                rr["_source_csv"] = str(csv_path)
                rows.append(rr)
    return rows


def triple_ok_only_analysis(run_dirs: Sequence[str | Path], *, top_k: int = 5) -> dict[str, Any]:
    all_rows: list[dict[str, Any]] = []
    for rd in run_dirs:
        all_rows.extend(_scan_csv_rows(Path(rd)))

    triple_rows = [r for r in all_rows if compute_triple_ok(r)]
    n_rows = len(all_rows)
    n_triple_ok = len(triple_rows)
    br_gaga_vals = [v for v in (finite_float(r.get("br_gaga")) for r in triple_rows) if v is not None]
    br_bb_vals = [v for v in (derive_br_bb(r) for r in triple_rows) if v is not None]
    width_vals = [v for v in (finite_float(r.get("total_width")) for r in triple_rows) if v is not None]

    def _top(rows: list[dict[str, Any]], key: str, reverse: bool) -> list[dict[str, Any]]:
        scored: list[tuple[float, dict[str, Any]]] = []
        for r in rows:
            val = derive_br_bb(r) if key == "br_bb" else finite_float(r.get(key))
            if val is not None:
                scored.append((val, r))
        scored.sort(key=lambda x: x[0], reverse=reverse)
        return [{"value": val, "source_csv": row.get("_source_csv")} for val, row in scored[:top_k]]

    return {
        "n_rows": n_rows,
        "n_triple_ok": n_triple_ok,
        "triple_ok_rate": (n_triple_ok / n_rows) if n_rows > 0 else None,
        "br_gaga_q50": quantile(br_gaga_vals, 0.50),
        "br_gaga_q90": quantile(br_gaga_vals, 0.90),
        "br_gaga_q95": quantile(br_gaga_vals, 0.95),
        "br_bb_median": quantile(br_bb_vals, 0.50),
        "br_bb_q90": quantile(br_bb_vals, 0.90),
        "total_width_median": quantile(width_vals, 0.50),
        "top_br_gaga": _top(triple_rows, "br_gaga", True),
        "top_low_br_bb": _top(triple_rows, "br_bb", False),
        "ranking_basis": "triple_ok_rows_only",
    }


@dataclass
class GateFailure:
    code: str
    severity: str
    explanation: str
    owner: str
    next_check: str
    artifact_path: str | None = None


def _task_summary_has_fail(path: Path) -> tuple[bool, str | None]:
    if not path.exists() or not path.is_file():
        return True, "missing task_summary"
    try:
        with path.open("r", encoding="utf-8") as h:
            for i, line in enumerate(h, start=1):
                t = line.strip()
                if not t:
                    continue
                ev = json.loads(t)
                if isinstance(ev, Mapping):
                    kind = str(ev.get("event") or ev.get("status") or ev.get("outcome") or "").lower()
                    if any(tok in kind for tok in ("fail", "crash", "error", "timeout")):
                        return True, f"failure event at line {i}: {kind}"
        return False, None
    except Exception as exc:
        return True, f"invalid task_summary.jsonl: {exc}"


def apply_gates(*, commands: Sequence[Sequence[str]], run_results: Sequence[Mapping[str, Any]], campaign_state: Mapping[str, Any], rerun_manifest: Mapping[str, Any], analysis: Mapping[str, Any], min_triple_ok_rate: float, production_validation: bool, evidence_mode: str, operational_only: bool = False) -> dict[str, Any]:
    failures: list[GateFailure] = []
    warnings: list[GateFailure] = []

    def fail(code: str, explanation: str, owner: str, next_check: str, artifact: str | None = None, severity: str = "error") -> None:
        obj = GateFailure(code=code, severity=severity, explanation=explanation, owner=owner, next_check=next_check, artifact_path=artifact)
        (failures if severity == "error" else warnings).append(obj)

    all_threads_flag = all(("--threads" in cmd and cmd[cmd.index("--threads") + 1] == "1") for cmd in commands)
    if not all_threads_flag:
        fail("gate_failed_due_to_runtime_policy", "Not all commands include --threads 1", "safe_automation_layer", "inspect build_commands")

    if analysis.get("ranking_basis") != "triple_ok_rows_only":
        fail("gate_failed_due_to_non_triple_ok_ranking_basis", "ranking_basis must be triple_ok_rows_only", "safe_automation_layer", "inspect triple_ok_only_analysis")

    if production_validation:
        runs = campaign_state.get("runs", []) if isinstance(campaign_state.get("runs"), list) else []
        for run in runs:
            if not isinstance(run, Mapping):
                continue
            manifest_path = Path(str(run.get("run_manifest_path") or ""))
            if not manifest_path.exists():
                fail("gate_failed_due_to_missing_evidence", "missing run_manifest.json", "campaign_state/run_health", "check run artifact generation", str(manifest_path))
            else:
                try:
                    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
                    runtime = payload.get("runtime") if isinstance(payload, Mapping) else None
                    if not isinstance(runtime, Mapping):
                        fail("gate_failed_due_to_invalid_evidence", "manifest missing runtime object", "run_manifest", "inspect run_manifest schema", str(manifest_path))
                    elif "omp_num_threads" not in runtime:
                        fail("gate_failed_due_to_missing_evidence", "manifest missing runtime.omp_num_threads", "run_manifest", "record omp_num_threads in manifest", str(manifest_path))
                    elif int(runtime.get("omp_num_threads")) != 1:
                        fail("gate_failed_due_to_runtime_policy", "runtime.omp_num_threads must equal 1", "run_manifest", "enforce OMP=1 at runtime", str(manifest_path))
                except Exception as exc:
                    fail("gate_failed_due_to_invalid_evidence", f"invalid manifest JSON: {exc}", "run_manifest", "regenerate manifest", str(manifest_path))

            log_path = Path(str(run.get("orchestrator_log_path") or ""))
            if not log_path.exists():
                fail("gate_failed_due_to_missing_evidence", "missing orchestrator.log", "orchestrator", "ensure log is emitted", str(log_path))
            else:
                txt = log_path.read_text(encoding="utf-8", errors="replace")
                if OMP_AUDIT_LINE not in txt:
                    fail("gate_failed_due_to_runtime_policy", "orchestrator.log missing OMP audit line", "orchestrator", "emit [CONF] OMP_NUM_THREADS = 1", str(log_path))
                low = txt.lower()
                if any(sig in low for sig in GSL_ABORT_SIGNATURES):
                    fail("gate_failed_due_to_gsl_signature", "GSL/SIGABRT signature detected in orchestrator.log", "run_health", "inspect stderr + stack evidence", str(log_path))

            task_path = Path(str(run.get("task_summary_path") or ""))
            has_fail, detail = _task_summary_has_fail(task_path)
            if has_fail:
                code = "gate_failed_due_to_missing_evidence" if not task_path.exists() else "gate_failed_due_to_invalid_evidence"
                if detail and "failure event" in detail:
                    code = "gate_failed_due_to_health"
                fail(code, detail or "task_summary invalid", "run_health", "inspect task_summary.jsonl", str(task_path))

    for slice_rec in campaign_state.get("slices", []):
        if not isinstance(slice_rec, Mapping):
            continue
        csv_path = Path(str(slice_rec.get("source_csv") or ""))
        if not csv_path.exists():
            fail("gate_failed_due_to_missing_evidence", "missing expected CSV", "run_health", "rerun missing slice", str(csv_path))
            continue
        row_count = int(slice_rec.get("row_count", 0) or 0)
        if row_count <= 0:
            fail("gate_failed_due_to_header_only_csv", "CSV is empty or header-only", "run_health", "rerun slice and verify rows", str(csv_path))

    combined = "\n".join(str(r.get("stderr", "")) + "\n" + str(r.get("stdout", "")) for r in run_results).lower()
    if any(sig in combined for sig in GSL_ABORT_SIGNATURES):
        fail("gate_failed_due_to_gsl_signature", "GSL/SIGABRT signature detected in execution outputs", "safe_automation_layer", "inspect subprocess stderr")

    summary = campaign_state.get("summary") if isinstance(campaign_state.get("summary"), Mapping) else {}
    health_failed = int(summary.get("n_health_failed", 0) or 0)
    if health_failed != 0:
        fail("gate_failed_due_to_health", f"health_failed must be 0 (got {health_failed})", "campaign_state", "inspect failed runs")

    rsum = rerun_manifest.get("summary") if isinstance(rerun_manifest.get("summary"), Mapping) else {}
    n_runs_sel = int(rsum.get("n_runs_selected", 0) or 0)
    n_slices_sel = int(rsum.get("n_slices_selected", 0) or 0)
    if n_runs_sel != 0 or n_slices_sel != 0:
        fail("gate_failed_due_to_rerun_manifest", f"rerun selections must be 0 (runs={n_runs_sel}, slices={n_slices_sel})", "rerun_manifest", "resolve run/slice failures")

    tok_rate = analysis.get("triple_ok_rate")
    tok_ok = isinstance(tok_rate, (int, float)) and float(tok_rate) >= float(min_triple_ok_rate)
    if not operational_only and not tok_ok:
        fail("gate_failed_due_to_physics_viability", f"triple_ok_rate={tok_rate} below threshold={min_triple_ok_rate}", "safe_automation_layer", "improve triple_ok viability")

    return {
        "passed": len(failures) == 0,
        "failures": [f.__dict__ for f in failures],
        "warnings": [w.__dict__ for w in warnings],
        "production_validation": production_validation,
        "evidence_mode": evidence_mode,
        "all_commands_have_threads_1": all_threads_flag,
    }


def run_safe_campaign(contract: Mapping[str, Any], *, workdir: str | Path, execute: bool = True) -> dict[str, Any]:
    validated = validate_contract(contract)
    commands = build_commands(validated)
    campaign = validated["campaign"]
    outdir = Path(validated["runtime"]["outdir"])
    campaign_dir = outdir / str(validated["runtime"]["lake_name"]) / str(campaign["name"])

    run_results: list[dict[str, Any]] = []
    if execute:
        for cmd in commands[: int(campaign["max_runs"])]:
            run_results.append(_run_command(cmd, Path(workdir)))

    run_dirs = [p for p in campaign_dir.glob("*/*") if p.is_dir()] if campaign_dir.exists() else []
    run_filter = validated.get("run_filter") if isinstance(validated.get("run_filter"), Mapping) else {}
    run_name_filter = str(run_filter.get("run_name", "")).strip() if isinstance(run_filter, Mapping) else ""
    if run_name_filter:
        run_dirs = [p for p in run_dirs if p.name == run_name_filter or p.parent.name == run_name_filter]

    state = build_campaign_state(campaign=str(campaign["name"]), campaign_dir=campaign_dir, score_runs=False)
    if run_name_filter:
        state = dict(state)
        state["runs"] = [r for r in state.get("runs", []) if isinstance(r, Mapping) and str(r.get("run_name", "")) == run_name_filter]
        scoped_run_dirs = {str(r.get("run_dir")) for r in state["runs"] if isinstance(r, Mapping)}
        state["slices"] = [s for s in state.get("slices", []) if isinstance(s, Mapping) and str(s.get("run_dir", "")) in scoped_run_dirs]
        summary = dict(state.get("summary", {})) if isinstance(state.get("summary"), Mapping) else {}
        summary["n_runs"] = len(state["runs"])
        summary["n_slices"] = len(state["slices"])
        summary["n_health_failed"] = sum(1 for r in state["runs"] if isinstance(r, Mapping) and str(r.get("health_status", "")) != "ok")
        state["summary"] = summary
    rerun = build_rerun_manifest(state, source_campaign_state="inline")
    analysis = triple_ok_only_analysis(run_dirs) if run_dirs else {
        "n_rows": 0,
        "n_triple_ok": 0,
        "triple_ok_rate": None,
        "ranking_basis": "triple_ok_rows_only",
        "top_br_gaga": [],
        "top_low_br_bb": [],
    }

    contract_evidence_mode = str(validated.get("evidence_mode", "")).strip().lower()
    if execute:
        evidence_mode = "production"
    else:
        if contract_evidence_mode in {"existing_artifacts", "production", "fixture", "no_execute"}:
            evidence_mode = contract_evidence_mode
        else:
            evidence_mode = "fixture" if run_dirs else "no_execute"

    enforce_production_gates = execute or evidence_mode in {"production", "existing_artifacts"}

    gate_eval = apply_gates(
        commands=commands,
        run_results=run_results,
        campaign_state=state,
        rerun_manifest=rerun,
        analysis=analysis,
        min_triple_ok_rate=float(validated["gates"]["min_triple_ok_rate"]),
        production_validation=enforce_production_gates,
        evidence_mode=evidence_mode,
        operational_only=bool(validated.get("gates", {}).get("operational_only", False)),
    )

    production_validation = bool(gate_eval["passed"] and evidence_mode in {"production", "existing_artifacts"})

    stop_report = {
        "schema_version": "safe_automation_v1_1_stop_report",
        "generated_utc": utc_now_iso(),
        "campaign": campaign["name"],
        "bounded": {"max_runs": campaign["max_runs"], "expansion_max_runs": campaign["expansion_max_runs"]},
        "commands": commands,
        "command_count": len(commands),
        "execution": run_results,
        "campaign_state_summary": state.get("summary"),
        "artifact_roots_used_for_validation": {
            "campaign_dir": str(campaign_dir),
            "run_filter": {"run_name": run_name_filter} if run_name_filter else {},
            "scoped_run_dirs": [str(p) for p in run_dirs],
        },
        "rerun_manifest_summary": rerun.get("summary"),
        "triple_ok_analysis": analysis,
        "gates": gate_eval,
        "evidence_mode": evidence_mode,
        "production_validation": production_validation,
        "expansion_allowed": bool(gate_eval["passed"] and campaign["expansion_max_runs"] > 0),
        "ranking_basis": analysis.get("ranking_basis"),
        "all_commands_have_threads_1": gate_eval["all_commands_have_threads_1"],
        "all_manifests_have_omp_threads_1": all("runtime.omp_num_threads must equal 1" not in f["explanation"] and "missing runtime.omp_num_threads" not in f["explanation"] for f in gate_eval["failures"]),
        "all_logs_have_omp_audit": all("missing OMP audit line" not in f["explanation"] for f in gate_eval["failures"]),
        "all_csvs_have_rows": all("header-only" not in f["explanation"] for f in gate_eval["failures"]),
    }

    return {
        "commands": commands,
        "campaign_state": state,
        "rerun_manifest": rerun,
        "triple_ok_analysis": analysis,
        "gates": gate_eval,
        "stop_report": stop_report,
    }


def write_outputs(result: Mapping[str, Any], output_dir: str | Path) -> dict[str, str]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    p_state = out / "campaign_state.json"
    p_rerun = out / "rerun_manifest.json"
    p_analysis = out / "triple_ok_analysis.json"
    p_stop = out / "stop_report.json"
    write_campaign_state_json(result["campaign_state"], p_state)
    write_rerun_manifest_json(result["rerun_manifest"], p_rerun)
    p_analysis.write_text(json.dumps(result["triple_ok_analysis"], indent=2, sort_keys=True) + "\n", encoding="utf-8")
    p_stop.write_text(json.dumps(result["stop_report"], indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"campaign_state": str(p_state), "rerun_manifest": str(p_rerun), "triple_ok_analysis": str(p_analysis), "stop_report": str(p_stop)}


def _cli(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Safe Automation Layer v1.1 runner")
    ap.add_argument("--contract", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--no-execute", action="store_true")
    args = ap.parse_args(argv)

    contract = load_contract(args.contract)
    result = run_safe_campaign(contract, workdir=args.workdir, execute=not args.no_execute)
    paths = write_outputs(result, args.output_dir)
    print(json.dumps({"ok": True, "gates_passed": result["gates"]["passed"], "paths": paths}, sort_keys=True))
    return 0 if result["gates"]["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(_cli())
