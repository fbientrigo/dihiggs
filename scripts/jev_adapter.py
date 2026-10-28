#!/usr/bin/env python3
"""Canonical execution adapter for Jev-selected DiHiggs parameter-space templates.

Bridges the Jev policy layer to the canonical Lambda1EvaluatorV2 producer (schema
dihiggs.lambda1.v2, engine lambda1_v2).

Key invariants:
- Canonical lifetime unit is millimetres (ctau_mm). Explicit conversion to/from metres.
- Evaluator is dihiggs/app/Lambda1EvaluatorV2. Legacy PhysScanWithFixings is strictly rejected.
- Each template is decomposed into fixed (mA, lambda6) subruns using Lambda1Fixed and cartesian_rows.
- For ctau_mphi_mA_refine: 9 subruns of 120 points = 1080 total requested points.
- All outputs belong under RUNS_ROOT / "jev_harness" (never in the source tree).
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from dihiggs.app.orchestrator.io_utils import detect_git_info
from dihiggs.app.orchestrator.lambda1_v2 import (
    INPUT_HEADER,
    Lambda1Fixed,
    SCHEMA_VERSION as CANONICAL_SCHEMA,
    cartesian_rows,
    grid_values,
    run_lambda1_v2,
    validate_output,
)
from dihiggs.app.orchestrator.physics_authority import active_mass_gev, convention_provenance

WORKSPACE_ROOT = REPO_ROOT.parent
RUNS_ROOT = WORKSPACE_ROOT / "runs"
DEFAULT_OUTPUT_ROOT = RUNS_ROOT / "jev_harness"
CANONICAL_PRODUCER = "Lambda1EvaluatorV2"
CANONICAL_ENGINE = "lambda1_v2"
DEFAULT_EXECUTABLE = REPO_ROOT / "dihiggs" / "app" / CANONICAL_PRODUCER

DEFAULT_TAN_BETA_CENTER = 124416.0
DEFAULT_LAMBDA6_CENTER = 0.0019683

# Declared reference anchor from GEMINI.md section 4 / golden case L07
# ("mphi=200, mA=500, tan_beta=126904, lambda6~0.0019, lambda1=1.0"). This is
# NOT a recovered historical-best artifact (none could be found); m_h always
# comes from active_mass_gev(), never from this dict.
REFERENCE_ANCHOR_L07 = {
    "mA": 500.0,
    "mHp": 500.0,
    "sin_ba": 1.0,
    "lambda6": 0.0019,
    "lambda7": 0.0,
    "mH": 200.0,
    "lambda1": 1.0,
    "tan_beta": 126904.0,
}

TEMPLATES: dict[str, dict[str, Any]] = {
    "ctau_mphi_mA_refine": {
        "mA_values": [450.0, 475.0, 500.0],
        "lambda6_local_multipliers": [0.97, 1.0, 1.03],
        "tan_beta_local_multipliers": [0.98, 1.0, 1.02],
        "mphi_min": 180.0,
        "mphi_max": 220.0,
        "n_mphi": 40,
        "n_lam1": 1,
    },
    "mphi_lower_edge_probe": {
        "mA_values": [500.0],
        "lambda6_local_multipliers": [0.97, 1.0, 1.03],
        "tan_beta_local_multipliers": [0.96, 0.98, 1.0, 1.02, 1.04],
        "mphi_min": 160.0,
        "mphi_max": 220.0,
        "n_mphi": 80,
        "n_lam1": 1,
    },
    "mA_upper_edge_probe": {
        "mA_values": [450.0, 475.0, 500.0],
        "lambda6_local_multipliers": [1.0],
        "tan_beta_local_multipliers": [0.98, 1.0, 1.02],
        "mphi_min": 180.0,
        "mphi_max": 220.0,
        "n_mphi": 120,
        "n_lam1": 1,
    },
    "sensitivity_probe": {
        "mA_values": [500.0],
        "lambda6_local_multipliers": [1.0],
        "tan_beta_local_multipliers": [1.0],
        "mphi_min": 180.0,
        "mphi_max": 220.0,
        "n_mphi": 334,
        "n_lam1": 1,
        "allow_lambda1_sensitivity": False,
    },
    "control_box": {
        "mA_values": [500.0],
        "lambda6_local_multipliers": [1.0],
        "tan_beta_local_multipliers": [1.0],
        "mphi_min": 195.0,
        "mphi_max": 205.0,
        "n_mphi": 1000,
        "n_lam1": 1,
    },
}

ALLOWED_OVERRIDES = {
    "mphi_min",
    "mphi_max",
    "n_mphi",
    "mA_values",
    "tan_beta_center",
    "lambda6_center",
    "tan_beta_local_multipliers",
    "lambda6_local_multipliers",
    "variable",
    "relative_steps",
    "anchor",
    "allow_lambda1_sensitivity",
}

FORBIDDEN_OVERRIDES = {
    "lambda1",
    "n_lam1",
    "scoring",
    "triple_ok",
    "envelope",
    "command",
    "shell",
    "bash",
    "exec_path",
    "force",
}


def ctau_m_to_mm(m: float) -> float:
    """Explicit unit conversion: metres to millimetres."""
    return float(m) * 1000.0


def ctau_mm_to_m(mm: float) -> float:
    """Explicit unit conversion: millimetres to metres."""
    return float(mm) / 1000.0


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _inside_range(vals: list[float], lo: float, hi: float) -> bool:
    return all(lo <= float(v) <= hi for v in vals)


def resolve_template_grid(
    template_name: str,
    envelope: Mapping[str, Any],
    overrides: Mapping[str, Any] | None = None,
    min_points: int = 1000,
    max_points: int = 3000,
    tan_beta_center_default: float = DEFAULT_TAN_BETA_CENTER,
    lambda6_center_default: float = DEFAULT_LAMBDA6_CENTER,
) -> dict[str, Any]:
    """Resolve and validate the parameter grid for an allowlisted template."""
    if template_name not in TEMPLATES:
        raise ValueError(f"unknown template: {template_name}")
    spec = dict(TEMPLATES[template_name])
    ov = dict(overrides or {})

    for k in ov:
        if k in FORBIDDEN_OVERRIDES:
            raise ValueError(f"forbidden override: {k}")
        if k not in ALLOWED_OVERRIDES:
            raise ValueError(f"override not allowed: {k}")

    spec.update({k: ov[k] for k in ov if k in ALLOWED_OVERRIDES})

    if template_name == "sensitivity_probe":
        variable = str(ov.get("variable", "tan_beta")).strip()
        if variable not in {"tan_beta", "lambda6", "mA", "lambda1"}:
            raise ValueError("sensitivity_probe variable must be one of tan_beta|lambda6|mA|lambda1")
        rel_steps_raw = ov.get("relative_steps", [-0.10, 0.0, 0.10])
        if not isinstance(rel_steps_raw, list) or len(rel_steps_raw) == 0:
            raise ValueError("sensitivity_probe relative_steps must be non-empty list")
        rel_steps = [float(v) for v in rel_steps_raw]
        anchor = dict(ov.get("anchor", {})) if isinstance(ov.get("anchor"), dict) else {}
        allow_l1 = bool(ov.get("allow_lambda1_sensitivity", spec.get("allow_lambda1_sensitivity", False)))

        anchor_tan_beta = float(anchor.get("tan_beta", ov.get("tan_beta_center", tan_beta_center_default)))
        anchor_lambda6 = float(anchor.get("lambda6", ov.get("lambda6_center", lambda6_center_default)))
        anchor_mA = float(anchor.get("mA", spec.get("mA_values", [500.0])[0]))
        anchor_lambda1 = float(anchor.get("lambda1", 1.0))
        spec["mphi_min"] = float(anchor.get("mphi_min", spec.get("mphi_min", 180.0)))
        spec["mphi_max"] = float(anchor.get("mphi_max", spec.get("mphi_max", 220.0)))
        spec["n_lam1"] = 1

        def _rel_values(base: float) -> list[float]:
            return sorted({base * (1.0 + s) for s in rel_steps})

        if variable == "tan_beta":
            spec["tan_beta_values"] = _rel_values(anchor_tan_beta)
            spec["lambda6_values"] = [anchor_lambda6]
            spec["mA_values"] = [anchor_mA]
        elif variable == "lambda6":
            spec["tan_beta_values"] = [anchor_tan_beta]
            spec["lambda6_values"] = _rel_values(anchor_lambda6)
            spec["mA_values"] = [anchor_mA]
        elif variable == "mA":
            spec["tan_beta_values"] = [anchor_tan_beta]
            spec["lambda6_values"] = [anchor_lambda6]
            spec["mA_values"] = _rel_values(anchor_mA)
        else:
            if not allow_l1:
                raise ValueError("sensitivity_probe lambda1 variation disabled unless allow_lambda1_sensitivity=true")
            spec["tan_beta_values"] = [anchor_tan_beta]
            spec["lambda6_values"] = [anchor_lambda6]
            spec["mA_values"] = [anchor_mA]
            l1_values = _rel_values(anchor_lambda1)
            if len(l1_values) < 2:
                raise ValueError("sensitivity_probe lambda1 variation requires at least two distinct lambda1 values")
            spec["lambda1_values"] = l1_values
            spec["n_lam1"] = len(l1_values)

    mA_values = [float(v) for v in spec["mA_values"]]
    lambda6_local_multipliers = [float(v) for v in spec.get("lambda6_local_multipliers", [1.0])]
    tan_beta_local_multipliers = [float(v) for v in spec.get("tan_beta_local_multipliers", [1.0])]
    mphi_min = float(spec["mphi_min"])
    mphi_max = float(spec["mphi_max"])
    n_mphi = int(spec["n_mphi"])
    lambda1_values = [float(v) for v in spec.get("lambda1_values", [1.0])]

    if mphi_min > mphi_max:
        raise ValueError("mphi_min > mphi_max")
    if not _inside_range(mA_values, float(envelope["mA"][0]), float(envelope["mA"][1])):
        raise ValueError("mA_values outside envelope")
    if not (float(envelope["mphi"][0]) <= mphi_min <= float(envelope["mphi"][1])):
        raise ValueError("mphi_min outside envelope")
    if not (float(envelope["mphi"][0]) <= mphi_max <= float(envelope["mphi"][1])):
        raise ValueError("mphi_max outside envelope")

    tan_beta_values = [float(v) for v in spec.get("tan_beta_values", [])]
    lambda6_values = [float(v) for v in spec.get("lambda6_values", [])]

    if tan_beta_values:
        if not _inside_range(tan_beta_values, float(envelope["tan_beta"][0]), float(envelope["tan_beta"][1])):
            raise ValueError("tan_beta_values outside envelope")
        tan_beta_center = float(tan_beta_values[len(tan_beta_values) // 2])
    else:
        tan_beta_center = float(ov.get("tan_beta_center", tan_beta_center_default))
        if not (float(envelope["tan_beta"][0]) <= tan_beta_center <= float(envelope["tan_beta"][1])):
            raise ValueError("tan_beta_center outside envelope")

    if lambda6_values:
        if not _inside_range(lambda6_values, float(envelope["lambda6"][0]), float(envelope["lambda6"][1])):
            raise ValueError("lambda6_values outside envelope")
        lambda6_center = float(lambda6_values[len(lambda6_values) // 2])
    else:
        lambda6_center = float(ov.get("lambda6_center", lambda6_center_default))
        if not (float(envelope["lambda6"][0]) <= lambda6_center <= float(envelope["lambda6"][1])):
            raise ValueError("lambda6_center outside envelope")

    if not tan_beta_values:
        tan_beta_values = sorted({
            tan_beta_center * m
            for m in tan_beta_local_multipliers
            if float(envelope["tan_beta"][0]) <= tan_beta_center * m <= float(envelope["tan_beta"][1])
        })
    if not lambda6_values:
        lambda6_values = sorted({
            lambda6_center * m
            for m in lambda6_local_multipliers
            if float(envelope["lambda6"][0]) <= lambda6_center * m <= float(envelope["lambda6"][1])
        })

    if not tan_beta_values:
        raise ValueError("tan_beta materialization produced no in-envelope values")
    if not lambda6_values:
        raise ValueError("lambda6 materialization produced no in-envelope values")

    points_per_mphi = int(len(mA_values) * len(lambda6_values) * len(tan_beta_values) * len(lambda1_values))
    real_points = int(points_per_mphi * n_mphi)
    if real_points < int(min_points):
        needed_n_mphi = max(1, (int(min_points) + points_per_mphi - 1) // points_per_mphi)
        n_mphi = int(needed_n_mphi)
        real_points = int(points_per_mphi * n_mphi)

    if real_points < int(min_points):
        raise ValueError(f"template real_points ({real_points}) below min_points_per_subrun ({min_points})")
    if real_points > int(max_points):
        raise ValueError(f"template real_points ({real_points}) above max_points_per_subrun ({max_points})")

    return {
        "template_name": template_name,
        "mA_values": mA_values,
        "lambda6_values": lambda6_values,
        "tan_beta_values": tan_beta_values,
        "mphi_min": mphi_min,
        "mphi_max": mphi_max,
        "n_mphi": n_mphi,
        "lambda1_values": lambda1_values,
        "tan_beta_center": tan_beta_center,
        "lambda6_center": lambda6_center,
        "real_points": real_points,
    }


@dataclass(frozen=True)
class CanonicalSubrun:
    subrun_index: int
    subrun_name: str
    fixed: Lambda1Fixed
    mH_values: list[float]
    lambda1_values: list[float]
    tan_beta_values: list[float]
    row_count: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "subrun_index": self.subrun_index,
            "subrun_name": self.subrun_name,
            "row_count": self.row_count,
            "fixed": {
                "mh": self.fixed.mh,
                "mA": self.fixed.mA,
                "mHp": self.fixed.mHp,
                "sin_ba": self.fixed.sin_ba,
                "lambda6": self.fixed.lambda6,
                "lambda7": self.fixed.lambda7,
            },
            "tan_beta_values": self.tan_beta_values,
            "mH_values": self.mH_values,
            "lambda1_target_values": self.lambda1_values,
        }


@dataclass
class CanonicalPlan:
    template_name: str
    producer: str
    producer_schema: str
    engine: str
    total_requested_points: int
    total_row_count: int
    subrun_count: int
    subruns: list[CanonicalSubrun]
    fixed_controls: dict[str, Any]
    scanned_parameters: dict[str, Any]
    output_root: Path
    plan_hash: str
    plan_hash_short: str
    repository_provenance: dict[str, Any]
    mass_convention: dict[str, Any]
    iteration: int
    cycle_id: str
    iteration_id: str
    created_at: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": "jev_canonical_plan_v1",
            "template_name": self.template_name,
            "producer": self.producer,
            "producer_schema": self.producer_schema,
            "engine": self.engine,
            "total_requested_points": self.total_requested_points,
            "total_row_count": self.total_row_count,
            "subrun_count": self.subrun_count,
            "fixed_controls": self.fixed_controls,
            "scanned_parameters": self.scanned_parameters,
            "subruns": [s.to_dict() for s in self.subruns],
            "output_root": str(self.output_root),
            "plan_hash": self.plan_hash,
            "plan_hash_short": self.plan_hash_short,
            "repository_provenance": self.repository_provenance,
            "mass_convention": self.mass_convention,
            "iteration": self.iteration,
            "cycle_id": self.cycle_id,
            "iteration_id": self.iteration_id,
            "created_at": self.created_at,
        }


def plan_template_v2(
    template_name: str,
    cfg: Mapping[str, Any],
    *,
    iteration: int = 1,
    cycle_id: str = "cycle01",
    iteration_id: str = "iter01",
    overrides: Mapping[str, Any] | None = None,
    output_root: Path | None = None,
    repo_root: Path | None = None,
) -> CanonicalPlan:
    """Create the structured canonical V2 execution plan decomposing into fixed subruns."""
    resolved_repo = repo_root or REPO_ROOT
    resolved_output = output_root or DEFAULT_OUTPUT_ROOT
    envelope = cfg["envelope"]
    min_points = int(cfg.get("min_points_per_subrun", 1000))
    max_points = int(cfg.get("max_points_per_subrun", 3000))
    tb_center = float(cfg.get("tan_beta_center", DEFAULT_TAN_BETA_CENTER))
    l6_center = float(cfg.get("lambda6_center", DEFAULT_LAMBDA6_CENTER))

    active_mh = active_mass_gev()

    # Fail closed if an explicit mass override is provided and mismatches active authority
    explicit_mh = None
    if overrides:
        for k in ("mh", "mh_gev", "m_h_GeV", "m_h"):
            if k in overrides:
                explicit_mh = float(overrides[k])
                break
    if explicit_mh is None and cfg:
        for k in ("mh", "mh_gev", "m_h_GeV", "m_h"):
            if k in cfg:
                explicit_mh = float(cfg[k])
                break

    if explicit_mh is not None and abs(explicit_mh - active_mh) > 1e-4:
        raise ValueError(
            f"Mismatched explicit Higgs mass {explicit_mh} GeV does not match "
            f"active new-production authority {active_mh} GeV. "
            f"Fail closed: new production must use active_mass_gev()."
        )

    grid = resolve_template_grid(
        template_name=template_name,
        envelope=envelope,
        overrides=overrides,
        min_points=min_points,
        max_points=max_points,
        tan_beta_center_default=tb_center,
        lambda6_center_default=l6_center,
    )

    git = detect_git_info(resolved_repo)
    mass_conv = convention_provenance(
        source_commit=git.get("commit"),
        requested_mass_gev=active_mh,
    )

    mH_values = grid_values(grid["mphi_min"], grid["mphi_max"], grid["n_mphi"])
    lambda1_values = grid["lambda1_values"]
    tan_beta_values = grid["tan_beta_values"]

    subruns: list[CanonicalSubrun] = []
    idx = 0
    for mA in grid["mA_values"]:
        for l6 in grid["lambda6_values"]:
            fixed = Lambda1Fixed(
                mh=active_mh,
                mA=mA,
                mHp=mA,
                sin_ba=1.0,
                lambda6=l6,
                lambda7=0.0,
            )
            count = len(tan_beta_values) * len(mH_values) * len(lambda1_values)
            name = f"subrun_{idx:02d}_mA{int(mA)}_l6_{l6:.6f}".replace(".", "p")
            subruns.append(
                CanonicalSubrun(
                    subrun_index=idx,
                    subrun_name=name,
                    fixed=fixed,
                    mH_values=mH_values,
                    lambda1_values=lambda1_values,
                    tan_beta_values=tan_beta_values,
                    row_count=count,
                )
            )
            idx += 1

    total_rows = sum(s.row_count for s in subruns)
    if total_rows != grid["real_points"]:
        raise ValueError(
            f"Decomposition row count mismatch: subruns sum to {total_rows}, "
            f"grid computed {grid['real_points']}"
        )

    identity_payload = {
        "template_name": template_name,
        "producer": CANONICAL_PRODUCER,
        "schema": CANONICAL_SCHEMA,
        "engine": CANONICAL_ENGINE,
        "mA_values": grid["mA_values"],
        "lambda6_values": grid["lambda6_values"],
        "tan_beta_values": grid["tan_beta_values"],
        "mphi_min": grid["mphi_min"],
        "mphi_max": grid["mphi_max"],
        "n_mphi": grid["n_mphi"],
        "lambda1_values": grid["lambda1_values"],
        "total_rows": total_rows,
        "iteration": iteration,
        "cycle_id": cycle_id,
        "iteration_id": iteration_id,
    }
    plan_hash = hashlib.sha256(json.dumps(identity_payload, sort_keys=True).encode("utf-8")).hexdigest()

    return CanonicalPlan(
        template_name=template_name,
        producer=CANONICAL_PRODUCER,
        producer_schema=CANONICAL_SCHEMA,
        engine=CANONICAL_ENGINE,
        total_requested_points=total_rows,
        total_row_count=total_rows,
        subrun_count=len(subruns),
        subruns=subruns,
        fixed_controls={
            "mA": grid["mA_values"],
            "mHp": grid["mA_values"],
            "lambda6": grid["lambda6_values"],
            "lambda7": 0.0,
            "sin_ba": 1.0,
            "mh_gev": active_mh,
            "mass_convention_id": mass_conv["convention_id"],
        },
        scanned_parameters={
            "tan_beta": grid["tan_beta_values"],
            "mH_gev": [grid["mphi_min"], grid["mphi_max"], grid["n_mphi"]],
            "lambda1_target": grid["lambda1_values"],
        },
        output_root=resolved_output,
        plan_hash=plan_hash,
        plan_hash_short=plan_hash[:8],
        repository_provenance={
            "commit": git.get("commit", "unknown"),
            "is_dirty": git.get("is_dirty", False),
        },
        mass_convention=mass_conv,
        iteration=iteration,
        cycle_id=cycle_id,
        iteration_id=iteration_id,
        created_at=now_iso(),
    )


def materialize_canonical_plan(
    plan: CanonicalPlan,
    plan_dir: Path,
    *,
    dry_run: bool = True,
    executable: Path | None = None,
    repo_root: Path | None = None,
) -> dict[str, Any]:
    """Materialize the canonical plan to disk without running C++ when dry_run=True."""
    resolved_repo = repo_root or REPO_ROOT
    exec_path = executable or DEFAULT_EXECUTABLE
    plan_dir.mkdir(parents=True, exist_ok=True)
    plan_json_path = plan_dir / "canonical_plan.json"
    plan_json_path.write_text(json.dumps(plan.to_dict(), indent=2) + "\n", encoding="utf-8")

    subrun_summaries = []
    campaign_name = f"{plan.cycle_id}_{plan.iteration_id}"

    for subrun in plan.subruns:
        subrun_rows = cartesian_rows(
            fixed=subrun.fixed,
            mH_values=subrun.mH_values,
            lambda1_values=subrun.lambda1_values,
            tan_beta_values=subrun.tan_beta_values,
        )
        manifest = run_lambda1_v2(
            executable=exec_path,
            rows=subrun_rows,
            outdir=plan_dir,
            campaign=campaign_name,
            run_name=subrun.subrun_name,
            repo_root=resolved_repo,
            dry_run=True,
            force=True,
        )
        subrun_summaries.append({
            "subrun_index": subrun.subrun_index,
            "subrun_name": subrun.subrun_name,
            "rows": len(subrun_rows),
            "input_csv": manifest["input_csv"],
            "input_sha256": manifest["input_sha256"],
        })

    return {
        "template_name": plan.template_name,
        "plan_path": str(plan_json_path),
        "contract_path": str(plan_json_path),
        "plan_hash": plan.plan_hash,
        "plan_hash_short": plan.plan_hash_short,
        "contract_hash_short": plan.plan_hash_short,
        "planned_points": plan.total_row_count,
        "materialized_points": plan.total_row_count,
        "subrun_count": plan.subrun_count,
        "subruns": subrun_summaries,
        "plancheck_ok": True,
        "plancheck_rejection_reasons": [],
        "plancheck_stderr_tail": "",
        "plan_dir": str(plan_dir),
        "producer": plan.producer,
        "producer_schema": plan.producer_schema,
        "engine": plan.engine,
        "mass_convention": plan.mass_convention,
        "status": "planned",
    }


def validate_canonical_executable(executable: Path) -> Path:
    """Validate that the evaluator binary is the canonical Lambda1EvaluatorV2."""
    resolved = executable.expanduser().resolve()
    base = resolved.name
    if "PhysScanWithFixings" in base or "legacy" in base.lower():
        raise ValueError(
            f"Non-canonical executor rejected: {base} is legacy compatibility only. "
            f"New production must use {CANONICAL_PRODUCER}."
        )
    if not resolved.exists():
        raise FileNotFoundError(f"Canonical executable not found: {resolved}")
    if not os.access(resolved, os.X_OK):
        raise PermissionError(f"Canonical executable is not executable: {resolved}")
    return resolved


def execute_canonical_plan(
    plan: CanonicalPlan,
    plan_dir: Path,
    *,
    executable: Path | None = None,
    repo_root: Path | None = None,
    force: bool = False,
    timeout_s: float | None = None,
) -> dict[str, Any]:
    """Execute the canonical Lambda1EvaluatorV2 producer across all subruns in the plan."""
    resolved_repo = repo_root or REPO_ROOT
    exec_path = validate_canonical_executable(executable or DEFAULT_EXECUTABLE)

    manifests = []
    output_csvs: list[Path] = []
    campaign_name = f"{plan.cycle_id}_{plan.iteration_id}"

    for subrun in plan.subruns:
        subrun_rows = cartesian_rows(
            fixed=subrun.fixed,
            mH_values=subrun.mH_values,
            lambda1_values=subrun.lambda1_values,
            tan_beta_values=subrun.tan_beta_values,
        )
        manifest = run_lambda1_v2(
            executable=exec_path,
            rows=subrun_rows,
            outdir=plan_dir,
            campaign=campaign_name,
            run_name=subrun.subrun_name,
            repo_root=resolved_repo,
            dry_run=False,
            force=force,
            timeout_s=timeout_s,
        )
        manifests.append(manifest)
        output_csvs.append(Path(str(manifest["output_csv"])))

    outcome = extract_campaign_outcome(output_csvs)
    outcome["campaign"] = campaign_name
    outcome["template_name"] = plan.template_name
    outcome["plan_hash"] = plan.plan_hash
    outcome["status"] = "executed"
    return {
        "status": "executed",
        "subruns_executed": len(manifests),
        "manifests": manifests,
        "outcome": outcome,
    }


def extract_subrun_outcome(csv_path: Path) -> dict[str, Any]:
    """Extract theory acceptance, lifetime, and best point from a canonical output CSV."""
    if not csv_path.exists():
        raise FileNotFoundError(f"Output CSV not found: {csv_path}")

    total_evaluated_rows = 0
    construction_ok_rows = 0
    triple_ok_rows = 0
    width_ok_rows = 0
    valid_rows = 0
    failures_by_stage: dict[str, int] = {}
    failures_by_reason: dict[str, int] = {}

    best_ctau_mm: float | None = None
    best_point: dict[str, Any] | None = None

    with csv_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            total_evaluated_rows += 1
            schema = row.get("schema_version")
            if schema != CANONICAL_SCHEMA:
                raise ValueError(f"Schema mismatch: expected {CANONICAL_SCHEMA}, got {schema}")

            c_ok = int(row.get("construction_ok", 0) or 0) == 1
            t_ok = int(row.get("triple_ok", 0) or 0) == 1
            w_ok = int(row.get("width_ok", 0) or 0) == 1

            if c_ok:
                construction_ok_rows += 1
            if t_ok:
                triple_ok_rows += 1
            if w_ok:
                width_ok_rows += 1

            stage = str(row.get("rejection_stage", "unknown"))
            reason = str(row.get("rejection_reason", "none"))
            failures_by_stage[stage] = failures_by_stage.get(stage, 0) + 1
            if reason != "none":
                failures_by_reason[reason] = failures_by_reason.get(reason, 0) + 1

            ctau_raw = row.get("ctau_mm", "")
            ctau_val: float | None = None
            if ctau_raw and ctau_raw not in ("nan", "NaN", "null"):
                try:
                    parsed = float(ctau_raw)
                    if math.isfinite(parsed) and parsed > 0:
                        ctau_val = parsed
                except ValueError:
                    ctau_val = None

            is_valid = c_ok and t_ok and w_ok and (ctau_val is not None)
            if is_valid:
                valid_rows += 1
                assert ctau_val is not None
                if best_ctau_mm is None or ctau_val > best_ctau_mm:
                    best_ctau_mm = ctau_val
                    best_point = {
                        "point_id": row.get("point_id"),
                        "mh_gev": float(row.get("mh_input_gev", 0)),
                        "mH_gev": float(row.get("mH_input_gev", 0)),
                        "mA_gev": float(row.get("mA_input_gev", 0)),
                        "mHp_gev": float(row.get("mHp_input_gev", 0)),
                        "sin_beta_minus_alpha": float(row.get("sin_beta_minus_alpha_input", 0)),
                        "tan_beta": float(row.get("tan_beta_input", 0)),
                        "lambda1_target": float(row.get("lambda1_target", 0)),
                        "lambda6": float(row.get("lambda6_input", 0)),
                        "lambda7": float(row.get("lambda7_input", 0)),
                        "ctau_mm": ctau_val,
                        "ctau_m": ctau_mm_to_m(ctau_val),
                        "total_width_gev": float(row.get("total_width_gev", 0)),
                        "triple_ok": 1,
                        "width_ok": 1,
                        "rejection_stage": stage,
                    }

    triple_ok_fraction = (triple_ok_rows / total_evaluated_rows) if total_evaluated_rows > 0 else 0.0
    valid_fraction = (valid_rows / total_evaluated_rows) if total_evaluated_rows > 0 else 0.0

    return {
        "csv_path": str(csv_path),
        "total_evaluated_rows": total_evaluated_rows,
        "construction_ok_rows": construction_ok_rows,
        "triple_ok_rows": triple_ok_rows,
        "width_ok_rows": width_ok_rows,
        "valid_rows": valid_rows,
        "triple_ok_fraction": triple_ok_fraction,
        "valid_fraction": valid_fraction,
        "failures_by_stage": failures_by_stage,
        "failures_by_reason": failures_by_reason,
        "best_ctau_mm": best_ctau_mm,
        "best_ctau_m": ctau_mm_to_m(best_ctau_mm) if best_ctau_mm is not None else None,
        "best_point": best_point,
        "ctau_units": "mm",
        "unit_conversion_rule": "ctau_mm = ctau_m * 1000.0; ctau_m = ctau_mm / 1000.0",
    }


def extract_campaign_outcome(csv_paths: Sequence[Path]) -> dict[str, Any]:
    """Aggregate outcomes across all subruns in a campaign."""
    if not csv_paths:
        raise ValueError("No CSV paths provided for campaign outcome extraction")

    total_evaluated_rows = 0
    construction_ok_rows = 0
    triple_ok_rows = 0
    width_ok_rows = 0
    valid_rows = 0
    failures_by_stage: dict[str, int] = {}
    failures_by_reason: dict[str, int] = {}

    overall_best_ctau_mm: float | None = None
    overall_best_point: dict[str, Any] | None = None

    subrun_outcomes = []
    for p in csv_paths:
        sub = extract_subrun_outcome(p)
        subrun_outcomes.append(sub)
        total_evaluated_rows += sub["total_evaluated_rows"]
        construction_ok_rows += sub["construction_ok_rows"]
        triple_ok_rows += sub["triple_ok_rows"]
        width_ok_rows += sub["width_ok_rows"]
        valid_rows += sub["valid_rows"]
        for st, count in sub["failures_by_stage"].items():
            failures_by_stage[st] = failures_by_stage.get(st, 0) + count
        for re, count in sub["failures_by_reason"].items():
            failures_by_reason[re] = failures_by_reason.get(re, 0) + count

        sub_best_mm = sub["best_ctau_mm"]
        if sub_best_mm is not None:
            if overall_best_ctau_mm is None or sub_best_mm > overall_best_ctau_mm:
                overall_best_ctau_mm = sub_best_mm
                overall_best_point = sub["best_point"]

    triple_ok_fraction = (triple_ok_rows / total_evaluated_rows) if total_evaluated_rows > 0 else 0.0
    valid_fraction = (valid_rows / total_evaluated_rows) if total_evaluated_rows > 0 else 0.0

    return {
        "producer": CANONICAL_PRODUCER,
        "producer_schema": CANONICAL_SCHEMA,
        "total_evaluated_rows": total_evaluated_rows,
        "construction_ok_rows": construction_ok_rows,
        "triple_ok_rows": triple_ok_rows,
        "width_ok_rows": width_ok_rows,
        "valid_rows": valid_rows,
        "triple_ok_fraction": triple_ok_fraction,
        "valid_fraction": valid_fraction,
        "failures_by_stage": failures_by_stage,
        "failures_by_reason": failures_by_reason,
        "best_ctau_mm": overall_best_ctau_mm,
        "best_ctau_m": ctau_mm_to_m(overall_best_ctau_mm) if overall_best_ctau_mm is not None else None,
        "best_point": overall_best_point,
        "subrun_count": len(csv_paths),
        "ctau_units": "mm",
        "unit_conversion_rule": "ctau_mm = ctau_m * 1000.0; ctau_m = ctau_mm / 1000.0",
    }


def update_scientific_state(
    current_state: Mapping[str, Any],
    config: Mapping[str, Any],
    outcome: Mapping[str, Any],
    *,
    iteration: int,
    cycle_id: str,
    template_name: str,
) -> dict[str, Any]:
    """Produce the updated ScientificState after an executed outcome."""
    new_state = dict(current_state)
    total_eval = int(outcome.get("total_evaluated_rows", 0) or 0)
    valid_cnt = int(outcome.get("valid_rows", 0) or 0)
    rejected_cnt = total_eval - valid_cnt

    new_state["executed_points_total"] = int(new_state.get("executed_points_total", 0) or 0) + total_eval
    new_state["rejected_points_total"] = int(new_state.get("rejected_points_total", 0) or 0) + rejected_cnt
    new_state["iterations_completed"] = iteration
    new_state["cycle_id"] = cycle_id

    # Resolve previous best with explicit units
    if "best_ctau_mm" in new_state:
        prev_best_mm = float(new_state["best_ctau_mm"])
    elif "best_ctau_m" in new_state:
        prev_best_mm = ctau_m_to_mm(float(new_state["best_ctau_m"]))
    elif "baseline_best_ctau_mm" in config:
        prev_best_mm = float(config["baseline_best_ctau_mm"])
    else:
        prev_best_mm = ctau_m_to_mm(float(config["baseline_best_ctau_m"]))

    outcome_best_mm = outcome.get("best_ctau_mm")
    improved = False
    if outcome_best_mm is not None and float(outcome_best_mm) > prev_best_mm:
        improved = True
        new_best_mm = float(outcome_best_mm)
        new_state["best_ctau_mm"] = new_best_mm
        new_state["best_ctau_m"] = ctau_mm_to_m(new_best_mm)
        new_state["best_point"] = outcome.get("best_point")
        new_state["no_improvement_count"] = 0
    else:
        new_state["best_ctau_mm"] = prev_best_mm
        new_state["best_ctau_m"] = ctau_mm_to_m(prev_best_mm)
        new_state["no_improvement_count"] = int(new_state.get("no_improvement_count", 0) or 0) + 1

    summary = {
        "iteration": iteration,
        "cycle_id": cycle_id,
        "template_name": template_name,
        "executed_points": total_eval,
        "valid_points": valid_cnt,
        "triple_ok_fraction": outcome.get("triple_ok_fraction", 0.0),
        "best_ctau_mm": outcome_best_mm,
        "best_ctau_m": outcome.get("best_ctau_m"),
        "improved": improved,
    }
    new_state["last_iteration_summary_for_prompt"] = summary

    history = list(new_state.get("history", []))
    history.append(summary)
    new_state["history"] = history
    new_state["status"] = "in_progress"
    return new_state


def run_canonical_microrun(
    outdir: Path,
    executable: Path | None = None,
    repo_root: Path | None = None,
    point: Mapping[str, float] | None = None,
    campaign: str = "integration_smoke",
    run_name: str = "microrun",
) -> dict[str, Any]:
    """Run an intentionally tiny (1-point) canonical V2 smoke run to verify end-to-end integration.

    ``point`` overrides the default smoke coordinates (e.g. pass
    ``REFERENCE_ANCHOR_L07`` to replay the declared reference anchor). m_h is
    always taken from ``active_mass_gev()``, never from ``point``.
    """
    resolved_repo = repo_root or REPO_ROOT
    exec_path = validate_canonical_executable(executable or DEFAULT_EXECUTABLE)
    p = dict(point) if point else {}

    active_mh = active_mass_gev()
    fixed = Lambda1Fixed(
        mh=active_mh,
        mA=p.get("mA", 500.0),
        mHp=p.get("mHp", 500.0),
        sin_ba=p.get("sin_ba", 1.0),
        lambda6=p.get("lambda6", 0.0019),
        lambda7=p.get("lambda7", 0.0),
    )
    rows = cartesian_rows(
        fixed=fixed,
        mH_values=[p.get("mH", 200.0)],
        lambda1_values=[p.get("lambda1", 1.0)],
        tan_beta_values=[p.get("tan_beta", 124416.0)],
    )
    manifest = run_lambda1_v2(
        executable=exec_path,
        rows=rows,
        outdir=outdir,
        campaign=campaign,
        run_name=run_name,
        repo_root=resolved_repo,
        dry_run=False,
        force=True,
    )
    out_csv = Path(str(manifest["output_csv"]))
    outcome = extract_subrun_outcome(out_csv)
    return {
        "status": "complete",
        "manifest": manifest,
        "outcome": outcome,
    }
