from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class OrchestratorAdapter(Protocol):
    def dispatch(self, payload: Mapping[str, Any]) -> Mapping[str, Any]:
        ...


@dataclass(frozen=True)
class DispatchRecord:
    goal_id: str
    dispatched_at: str
    metadata: Mapping[str, Any] = field(default_factory=dict)
    status: str = "dry_run"


def _as_float(value: object) -> float | None:
    try:
        return float(value)
    except Exception:
        return None


def _as_int(value: object) -> int | None:
    try:
        iv = int(value)
    except Exception:
        return None
    return iv


def validate_proposer_payload(payload: Mapping[str, Any]) -> list[str]:
    errors = []
    if "status" not in payload:
        errors.append("Missing required key: 'status'")

    status = payload.get("status")
    if status not in ("success", "saturated"):
        errors.append(f"Invalid status: {status}")

    if "proposals" not in payload:
        errors.append("Missing required key: 'proposals'")
    else:
        proposals = payload.get("proposals")
        if not isinstance(proposals, Sequence) or isinstance(proposals, (str, bytes)):
            errors.append("'proposals' must be a sequence")
        else:
            for i, p in enumerate(proposals):
                if not isinstance(p, Mapping):
                    errors.append(f"Proposal at index {i} is not a mapping")
                    continue
                if "goal_id" not in p:
                    errors.append(f"Proposal at index {i} missing 'goal_id'")
                if "axes_binned" not in p:
                    errors.append(f"Proposal at index {i} missing 'axes_binned'")

    return errors


def normalize_runtime_config(raw: Mapping[str, Any] | None = None) -> dict[str, Any]:
    cfg = dict(raw or {})

    cpp_threads = _as_int(cfg.get("cpp_omp_threads"))
    legacy_threads = _as_int(cfg.get("threads"))
    effective_threads = cpp_threads if cpp_threads is not None else legacy_threads
    if effective_threads is None or effective_threads < 1:
        effective_threads = 1

    outdir = str(cfg.get("outdir", "/home/fabi/dihiggs/scripts/out"))
    lake_name = str(cfg.get("lake_name", "dihiggs_lake"))
    campaign = str(cfg.get("campaign", "autoresearch_campaign"))
    exec_path = str(cfg.get("exec_path", "/home/fabi/dihiggs/dihiggs/app/PhysScanWithFixings"))

    return {
        "exec_path": exec_path,
        "outdir": outdir,
        "lake_name": lake_name,
        "campaign": campaign,
        "cpp_omp_threads": effective_threads,
        "threads": effective_threads,
        "resume_scope": str(cfg.get("resume_scope", "run")),
        "materialize": bool(cfg.get("materialize", False)),
    }


def validate_executable_proposal(proposal: Mapping[str, Any], contract: Mapping[str, Any] | None = None) -> list[str]:
    errors: list[str] = []

    bounds = proposal.get("bounds")
    fixed = proposal.get("fixed")
    resolution = proposal.get("resolution")
    scan_values = proposal.get("scan_values")

    if not isinstance(bounds, Mapping):
        return ["bounds must be mapping"]
    if not isinstance(fixed, Mapping):
        return ["fixed must be mapping"]
    if not isinstance(resolution, Mapping):
        return ["resolution must be mapping"]
    if not isinstance(scan_values, Mapping):
        return ["scan_values must be mapping"]

    for k in ("sin_ba", "lambda7", "mA"):
        if _as_float(fixed.get(k)) is None:
            errors.append(f"fixed.{k} missing or non-numeric")

    tan_beta = scan_values.get("tan_beta")
    tanbeta_vals: list[float] = []
    if isinstance(tan_beta, Sequence) and not isinstance(tan_beta, (str, bytes)):
        for v in tan_beta:
            fv = _as_float(v)
            if fv is None:
                errors.append("scan_values.tan_beta contains non-numeric value")
                continue
            tanbeta_vals.append(fv)
    else:
        fv = _as_float(tan_beta)
        if fv is None:
            errors.append("scan_values.tan_beta missing or non-numeric")
        else:
            tanbeta_vals = [fv]
    if not tanbeta_vals:
        errors.append("scan_values.tan_beta must be non-empty")

    lambda6 = bounds.get("lambda6")
    if isinstance(lambda6, Sequence) and not isinstance(lambda6, (str, bytes)):
        vals = list(lambda6)
        if len(vals) != 2:
            errors.append("bounds.lambda6 must have exactly 2 entries")
        else:
            lo = _as_float(vals[0])
            hi = _as_float(vals[1])
            if lo is None or hi is None:
                errors.append("bounds.lambda6 values must be numeric")
            elif abs(lo - hi) > 0:
                errors.append("bounds.lambda6 must be singleton [x, x]")
    else:
        if _as_float(lambda6) is None:
            errors.append("bounds.lambda6 missing or non-numeric")

    for k in ("lambda1", "mH"):
        rng = bounds.get(k)
        if not (isinstance(rng, Sequence) and not isinstance(rng, (str, bytes)) and len(list(rng)) == 2):
            errors.append(f"bounds.{k} must be [min,max]")

    for k in ("mH", "lambda1"):
        iv = _as_int(resolution.get(k))
        if iv is None or iv < 1:
            errors.append(f"resolution.{k} must be integer >= 1")

    return errors


def build_orchestrator_command(
    proposal: Mapping[str, Any],
    runtime_config: Mapping[str, Any] | None = None,
    *,
    dry_run: bool = False,
    force: bool = False,
) -> list[str]:
    cfg = normalize_runtime_config(runtime_config)
    errors = validate_executable_proposal(proposal)
    if errors:
        raise ValueError("invalid executable proposal: " + "; ".join(errors))

    bounds = proposal["bounds"]
    fixed = proposal["fixed"]
    resolution = proposal["resolution"]
    scan_values = proposal["scan_values"]

    tan_beta = scan_values["tan_beta"]
    if isinstance(tan_beta, Sequence) and not isinstance(tan_beta, (str, bytes)):
        tanbeta_arg = ",".join(str(float(v)).rstrip("0").rstrip(".") if "." in str(float(v)) else str(int(float(v))) for v in tan_beta)
    else:
        tanbeta_arg = str(tan_beta)

    lambda6_bounds = bounds.get("lambda6")
    if isinstance(lambda6_bounds, Sequence) and not isinstance(lambda6_bounds, (str, bytes)):
        l6_val = float(list(lambda6_bounds)[0])
    else:
        l6_val = float(lambda6_bounds)

    lam1_min, lam1_max = list(bounds["lambda1"])
    mphi_min, mphi_max = list(bounds["mH"])

    run_name = str(proposal.get("proposal_id") or proposal.get("goal_id") or "autoresearch_run")

    cmd: list[str] = [
        "python",
        "dihiggs/app/orchestrate_scans.py",
        "--exec",
        str(cfg["exec_path"]),
        "--outdir",
        str(cfg["outdir"]),
        "--lake-name",
        str(cfg["lake_name"]),
        "--campaign",
        str(cfg["campaign"]),
        "--run-name",
        run_name,
        "--sin-ba",
        str(float(fixed["sin_ba"])),
        "--lambda6",
        str(l6_val),
        "--lambda7",
        str(float(fixed["lambda7"])),
        "--mA",
        str(float(fixed["mA"])),
        "--tanbeta",
        tanbeta_arg,
        "--mphi-min",
        str(float(mphi_min)),
        "--mphi-max",
        str(float(mphi_max)),
        "--n-mphi",
        str(int(resolution["mH"])),
        "--lam1-min",
        str(float(lam1_min)),
        "--lam1-max",
        str(float(lam1_max)),
        "--n-lam1",
        str(int(resolution["lambda1"])),
        "--threads",
        str(int(cfg["cpp_omp_threads"])),
        "--resume-scope",
        str(cfg["resume_scope"]),
    ]

    if bool(cfg.get("materialize", False)):
        cmd.append("--materialize")
    if dry_run:
        cmd.append("--dry-run")
    if force:
        cmd.append("--force")

    return cmd


class DryRunOrchestratorAdapter:
    def dispatch(self, payload: Mapping[str, Any]) -> Mapping[str, Any]:
        errors = validate_proposer_payload(payload)
        if errors:
            return {
                "status": "error",
                "errors": errors,
                "timestamp": datetime.now(timezone.utc).isoformat(),
            }

        proposals = payload.get("proposals", [])
        dispatch_records = []

        for proposal in proposals:
            record = DispatchRecord(
                goal_id=str(proposal.get("goal_id")),
                status="dry_run",
                dispatched_at=datetime.now(timezone.utc).isoformat(),
                metadata={
                    "axes_binned": proposal.get("axes_binned"),
                    "rationale": proposal.get("rationale"),
                },
            )
            dispatch_records.append(
                {
                    "goal_id": record.goal_id,
                    "status": record.status,
                    "dispatched_at": record.dispatched_at,
                    "metadata": record.metadata,
                }
            )

        return {
            "status": "success",
            "dispatch_count": len(dispatch_records),
            "records": dispatch_records,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
