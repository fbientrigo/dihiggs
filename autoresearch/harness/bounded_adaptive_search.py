from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from dataclasses import dataclass, replace
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

from autoresearch.harness.dihiggs_physics_score import finite_float, quantile
from autoresearch.harness.dihiggs_validators import compute_triple_ok, derive_br_bb
from autoresearch.harness.orchestrator_adapter import build_orchestrator_command
from autoresearch.harness.safe_automation_layer import run_safe_campaign


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def load_json(path: str | Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("contract must deserialize to object")
    return payload


def _pair(v: Any, key: str) -> tuple[float, float]:
    if not isinstance(v, list) or len(v) != 2:
        raise ValueError(f"{key} must be [lo, hi]")
    lo = float(v[0])
    hi = float(v[1])
    if lo > hi:
        raise ValueError(f"{key} must satisfy lo <= hi")
    return lo, hi


def validate_envelope(search_envelope: dict[str, Any]) -> dict[str, tuple[float, float]]:
    required = ["sin_ba", "mA", "lambda7", "lambda6", "tan_beta", "mphi", "lambda1"]
    out: dict[str, tuple[float, float]] = {}
    for key in required:
        out[key] = _pair(search_envelope[key], f"search_envelope.{key}")
    return out


def validate_budget(budget: dict[str, Any]) -> dict[str, Any]:
    out = {
        "max_iterations": int(budget["max_iterations"]),
        "min_raw_points_per_subrun": int(budget.get("min_raw_points_per_subrun", budget.get("min_points_per_subcampaign", 1000))),
        "max_raw_points_per_subrun": int(budget.get("max_raw_points_per_subrun", budget.get("max_points_per_subcampaign", 10000))),
        "max_total_points": int(budget["max_total_points"]),
        "max_walltime_minutes": budget.get("max_walltime_minutes"),
        "max_per_subrun_minutes": int(budget.get("max_per_subrun_minutes", 45)),
        "min_remaining_minutes_to_start": int(budget.get("min_remaining_minutes_to_start", 10)),
        "max_failed_subcampaigns": int(budget.get("max_failed_subcampaigns", 2)),
        "preferred_min_n_lam1": int(budget.get("preferred_min_n_lam1", 100)),
        "production_target_n_lam1": int(budget.get("production_target_n_lam1", 500)),
        "min_n_mphi_when_varied": int(budget.get("min_n_mphi_when_varied", 10)),
        "operational_smoke": bool(budget.get("operational_smoke", False)),
    }
    if out["min_raw_points_per_subrun"] < 1 or out["max_raw_points_per_subrun"] < out["min_raw_points_per_subrun"]:
        raise ValueError("invalid raw-point limits")
    return out


def _clamp_pair(pair: tuple[float, float], envelope: tuple[float, float]) -> tuple[float, float]:
    lo = max(pair[0], envelope[0])
    hi = min(pair[1], envelope[1])
    if lo > hi:
        raise ValueError("proposal escaped envelope")
    return lo, hi


def proposal_inside_envelope(ranges: dict[str, tuple[float, float]], envelope: dict[str, tuple[float, float]]) -> bool:
    return all(ranges[k][0] >= envelope[k][0] and ranges[k][1] <= envelope[k][1] for k in ranges)


def calculate_raw_points(tanbeta_values: list[float], n_mphi: int, n_lam1: int) -> int:
    return len(tanbeta_values) * int(n_mphi) * int(n_lam1)


def calculate_raw_points_multiplicative(
    tanbeta_values: list[float],
    lambda6_values: list[float],
    n_mphi: int,
    n_lam1: int,
) -> int:
    return len(tanbeta_values) * len(lambda6_values) * int(n_mphi) * int(n_lam1)


HBARC_GEV_M = 1.973269804e-16


def enforce_raw_point_rules(*, tanbeta_values: list[float], n_mphi: int, n_lam1: int, mphi_varies: bool, budget: dict[str, Any]) -> tuple[int, int, int, str | None]:
    min_pts = int(budget["min_raw_points_per_subrun"])
    max_pts = int(budget["max_raw_points_per_subrun"])
    preferred_min_n_lam1 = int(budget["preferred_min_n_lam1"])
    production_target_n_lam1 = int(budget["production_target_n_lam1"])
    min_n_mphi_when_varied = int(budget["min_n_mphi_when_varied"])

    n_lam1 = max(n_lam1, preferred_min_n_lam1)
    if mphi_varies:
        n_mphi = max(n_mphi, min_n_mphi_when_varied)

    pts = calculate_raw_points(tanbeta_values, n_mphi, n_lam1)
    if pts < min_pts:
        while pts < min_pts and n_lam1 < production_target_n_lam1:
            n_lam1 += 10
            pts = calculate_raw_points(tanbeta_values, n_mphi, n_lam1)
        while pts < min_pts and mphi_varies:
            n_mphi += 1
            pts = calculate_raw_points(tanbeta_values, n_mphi, n_lam1)
    if pts < min_pts and not bool(budget.get("operational_smoke", False)):
        return n_mphi, n_lam1, pts, "insufficient_points"

    while pts > max_pts and n_lam1 > preferred_min_n_lam1:
        n_lam1 -= 10
        pts = calculate_raw_points(tanbeta_values, n_mphi, n_lam1)
    while pts > max_pts and n_mphi > 1:
        n_mphi -= 1
        pts = calculate_raw_points(tanbeta_values, n_mphi, n_lam1)
    if pts > max_pts:
        return n_mphi, n_lam1, pts, "exceeds_max_points"
    return n_mphi, n_lam1, pts, None


@dataclass
class Subcampaign:
    subcampaign_id: str
    parent_iteration: int
    strategy: str
    variable_ranges: dict[str, tuple[float, float]]
    tanbeta_values: list[float]
    n_mphi: int
    n_lam1: int
    estimated_raw_points: int
    reason: str
    expected_output_campaign_name: str
    safety_checks: list[str]
    inside_envelope: bool
    lambda6_local_values: list[float] | None = None
    tan_beta_values: list[float] | None = None
    mA_values: list[float] | None = None
    expanded_subruns: list[dict[str, Any]] | None = None
    command_count: int | None = None
    template_name: str | None = None
    planned_points: int | None = None
    path_step: int | None = None
    tan_beta_center: float | None = None
    lambda6_center: float | None = None
    lambda1_fixed: float | None = None
    tan_beta_factor: float | None = None
    lambda6_factor: float | None = None
    objective_mode: str | None = None
    expected_direction: str | None = None
    lifecycle_state: str = "PENDING"

    def to_dict(self) -> dict[str, Any]:
        return {
            "subcampaign_id": self.subcampaign_id,
            "parent_iteration": self.parent_iteration,
            "strategy": self.strategy,
            "variable_ranges": {k: [v[0], v[1]] for k, v in self.variable_ranges.items()},
            "tanbeta_values": self.tanbeta_values,
            "tan_beta_values": self.tan_beta_values if self.tan_beta_values is not None else self.tanbeta_values,
            "lambda6_local_values": self.lambda6_local_values,
            "mA_values": self.mA_values,
            "expanded_subruns": self.expanded_subruns,
            "command_count": self.command_count,
            "template_name": self.template_name,
            "planned_points": self.planned_points if self.planned_points is not None else self.estimated_raw_points,
            "n_mphi": self.n_mphi,
            "n_lam1": self.n_lam1,
            "estimated_raw_points": self.estimated_raw_points,
            "reason": self.reason,
            "expected_output_campaign_name": self.expected_output_campaign_name,
            "safety_checks": self.safety_checks,
            "inside_envelope": self.inside_envelope,
            "path_step": self.path_step,
            "tan_beta_center": self.tan_beta_center,
            "lambda6_center": self.lambda6_center,
            "lambda1_fixed": self.lambda1_fixed,
            "tan_beta_factor": self.tan_beta_factor,
            "lambda6_factor": self.lambda6_factor,
            "objective_mode": self.objective_mode,
            "expected_direction": self.expected_direction,
            "lifecycle_state": self.lifecycle_state,
        }


def _local_refine_half_width(var_name: str, env: tuple[float, float], cfg: dict[str, Any]) -> float:
    span = max(0.0, float(env[1]) - float(env[0]))
    frac = float(cfg.get("refine_width_fraction", 0.1))
    abs_cfg = cfg.get("absolute_widths", {}) if isinstance(cfg.get("absolute_widths"), dict) else {}
    min_cfg = cfg.get("min_widths", {}) if isinstance(cfg.get("min_widths"), dict) else {}
    max_cfg = cfg.get("max_widths", {}) if isinstance(cfg.get("max_widths"), dict) else {}

    width = float(abs_cfg[var_name]) if var_name in abs_cfg else span * frac
    width = max(width, float(min_cfg.get(var_name, 0.0)))
    if var_name in max_cfg:
        width = min(width, float(max_cfg[var_name]))
    width = min(width, span)
    return max(0.0, width / 2.0)


def detect_ctau_metric_source(rows: list[dict[str, Any]]) -> str:
    if not rows:
        return "unavailable"
    keys = {k.lower() for r in rows for k in r.keys()}
    if any(k in keys for k in ("ctau", "c_tau", "decay_length", "lifetime")):
        return "explicit_column"
    if "total_width" in keys:
        return "derived_from_total_width"
    return "unavailable"


def _ctau_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    src = detect_ctau_metric_source(rows)
    tok = [r for r in rows if compute_triple_ok(r)]
    width_vals = [v for v in (finite_float(r.get("total_width")) for r in tok) if v is not None]
    out: dict[str, Any] = {
        "ctau_metric_source": src,
        "total_width_q50": quantile(width_vals, 0.50),
        "total_width_q10": quantile(width_vals, 0.10),
        "total_width_q05": quantile(width_vals, 0.05),
    }
    if src == "explicit_column":
        vals: list[float] = []
        for r in tok:
            for key in ("ctau", "c_tau", "decay_length", "lifetime"):
                v = finite_float(r.get(key))
                if v is not None:
                    vals.append(v)
                    break
        out.update({
            "ctau_m_q50": quantile(vals, 0.50),
            "ctau_m_q90": quantile(vals, 0.90),
            "ctau_m_q95": quantile(vals, 0.95),
            "ctau_m_max": max(vals) if vals else None,
        })
        return out
    if src == "derived_from_total_width":
        vals = [HBARC_GEV_M / v for v in width_vals if v > 0.0 and math.isfinite(v)]
        out.update({
            "ctau_m_q50": quantile(vals, 0.50),
            "ctau_m_q90": quantile(vals, 0.90),
            "ctau_m_q95": quantile(vals, 0.95),
            "ctau_m_max": max(vals) if vals else None,
        })
        return out
    return {"ctau_metric_source": "unavailable"}


def build_subcampaign(
    iteration: int,
    strategy: str,
    envelope: dict[str, tuple[float, float]],
    anchor: dict[str, float],
    budget: dict[str, Any],
    refine_cfg: dict[str, Any] | None = None,
    strategy_cfg: dict[str, Any] | None = None,
) -> Subcampaign:
    ranges = dict(envelope)
    cfg = dict(refine_cfg or {})
    s_cfg = dict(strategy_cfg or {})
    path_step = None
    tan_beta_center = None
    lambda6_center = None
    lambda1_fixed = None
    tan_beta_factor = None
    lambda6_factor = None
    tanbeta_locals: list[float] = []
    lambda6_locals: list[float] = []
    mA_values: list[float] = []
    expanded_subruns: list[dict[str, Any]] | None = None
    template_name: str | None = None
    planned_points: int | None = None
    objective_mode = str(s_cfg.get("objective_mode", "ctau_viability")) if s_cfg else None
    expected_direction = None

    if strategy == "coordinate_sweep":
        for k in ("tan_beta", "mphi", "lambda1"):
            c = anchor.get(k, (envelope[k][0] + envelope[k][1]) / 2)
            half = (envelope[k][1] - envelope[k][0]) * 0.2
            ranges[k] = _clamp_pair((c - half, c + half), envelope[k])
    elif strategy == "local_refine":
        for k in ("lambda6", "tan_beta", "mphi", "lambda1"):
            c = anchor.get(k, (envelope[k][0] + envelope[k][1]) / 2)
            half = _local_refine_half_width(k, envelope[k], cfg)
            ranges[k] = _clamp_pair((c - half, c + half), envelope[k])
    elif strategy == "box_partition":
        for k in ("lambda6", "tan_beta", "mphi", "lambda1"):
            lo, hi = envelope[k]
            ranges[k] = (lo, (lo + hi) / 2)
    elif strategy == "control_box":
        ranges["lambda6"] = _clamp_pair((0.0027, 0.0029), envelope["lambda6"])
        ranges["tan_beta"] = _clamp_pair((50000.0, 56000.0), envelope["tan_beta"])
        ranges["mphi"] = _clamp_pair((286.0, 292.0), envelope["mphi"])
        ranges["lambda1"] = _clamp_pair((4.1, 4.7), envelope["lambda1"])
    elif strategy == "multiplicative_ctau_path":
        path_step = int(s_cfg.get("path_step", max(0, iteration - 1)))
        tan_beta_factor = float(s_cfg.get("tan_beta_factor", 1.2))
        lambda6_factor = float(s_cfg.get("lambda6_factor", 0.9))
        anchor_cfg = dict(s_cfg.get("anchor", {})) if isinstance(s_cfg.get("anchor"), dict) else {}
        anchor_tan = float(anchor_cfg.get("tan_beta", anchor.get("tan_beta", (envelope["tan_beta"][0] + envelope["tan_beta"][1]) / 2)))
        anchor_l6 = float(anchor_cfg.get("lambda6", anchor.get("lambda6", (envelope["lambda6"][0] + envelope["lambda6"][1]) / 2)))
        tan_beta_center = anchor_tan * (tan_beta_factor ** path_step)
        lambda6_center = anchor_l6 * (lambda6_factor ** path_step)
        lambda1_fixed = float(s_cfg.get("lambda1_fixed", 1.0))
        if lambda1_fixed != 1.0:
            raise ValueError("multiplicative_ctau_path requires lambda1_fixed == 1.0")
        if not (envelope["tan_beta"][0] <= tan_beta_center <= envelope["tan_beta"][1]):
            raise ValueError("multiplicative_ctau_path tan_beta escaped envelope")
        if not (envelope["lambda6"][0] <= lambda6_center <= envelope["lambda6"][1]):
            raise ValueError("multiplicative_ctau_path lambda6 escaped envelope")
        l6_mults = [float(v) for v in s_cfg.get("lambda6_local_multiplier_values", [1.0])]
        tb_mults = [float(v) for v in s_cfg.get("tan_beta_local_multiplier_values", [1.0])]
        lambda6_locals = sorted({lambda6_center * m for m in l6_mults if envelope["lambda6"][0] <= lambda6_center * m <= envelope["lambda6"][1]})
        tanbeta_locals = sorted({tan_beta_center * m for m in tb_mults if envelope["tan_beta"][0] <= tan_beta_center * m <= envelope["tan_beta"][1]})
        if not lambda6_locals:
            raise ValueError("multiplicative_ctau_path has no in-envelope lambda6 local values")
        if not tanbeta_locals:
            raise ValueError("multiplicative_ctau_path has no in-envelope tan_beta local values")
        ranges["tan_beta"] = (min(tanbeta_locals), max(tanbeta_locals))
        ranges["lambda6"] = (min(lambda6_locals), max(lambda6_locals))
        ranges["lambda1"] = (lambda1_fixed, lambda1_fixed)
        if isinstance(anchor_cfg.get("mphi"), list) and len(anchor_cfg["mphi"]) == 2:
            ranges["mphi"] = _clamp_pair((float(anchor_cfg["mphi"][0]), float(anchor_cfg["mphi"][1])), envelope["mphi"])
        expected_direction = "larger_ctau_or_smaller_total_width"
    elif strategy == "template_grid_probe":
        template_name = str(s_cfg.get("template_name", "template_grid_probe"))
        lambda1_fixed = float(s_cfg.get("lambda1_fixed", 1.0))
        if lambda1_fixed != 1.0:
            raise ValueError("template_grid_probe requires lambda1_fixed == 1.0")
        n_lam1_req = int(s_cfg.get("n_lam1", 1))
        if n_lam1_req != 1:
            raise ValueError("template_grid_probe requires n_lam1 == 1")
        tan_beta_center = float(s_cfg.get("tan_beta_center", anchor.get("tan_beta", (envelope["tan_beta"][0] + envelope["tan_beta"][1]) / 2)))
        lambda6_center = float(s_cfg.get("lambda6_center", anchor.get("lambda6", (envelope["lambda6"][0] + envelope["lambda6"][1]) / 2)))
        mA_values = [float(v) for v in s_cfg.get("mA_values", [envelope["mA"][1]])]
        if not mA_values:
            raise ValueError("template_grid_probe requires non-empty mA_values")
        if any(v < envelope["mA"][0] or v > envelope["mA"][1] for v in mA_values):
            raise ValueError("template_grid_probe mA_values escaped envelope")
        if "lambda6_values" in s_cfg:
            lambda6_locals = sorted({float(v) for v in s_cfg.get("lambda6_values", []) if envelope["lambda6"][0] <= float(v) <= envelope["lambda6"][1]})
        else:
            l6_mults = [float(v) for v in s_cfg.get("lambda6_local_multiplier_values", [1.0])]
            lambda6_locals = sorted({lambda6_center * m for m in l6_mults if envelope["lambda6"][0] <= lambda6_center * m <= envelope["lambda6"][1]})
        if "tan_beta_values" in s_cfg:
            tanbeta_locals = sorted({float(v) for v in s_cfg.get("tan_beta_values", []) if envelope["tan_beta"][0] <= float(v) <= envelope["tan_beta"][1]})
        else:
            tb_mults = [float(v) for v in s_cfg.get("tan_beta_local_multiplier_values", [1.0])]
            tanbeta_locals = sorted({tan_beta_center * m for m in tb_mults if envelope["tan_beta"][0] <= tan_beta_center * m <= envelope["tan_beta"][1]})
        if not lambda6_locals:
            raise ValueError("template_grid_probe has no in-envelope lambda6 values")
        if not tanbeta_locals:
            raise ValueError("template_grid_probe has no in-envelope tan_beta values")
        mphi_min = float(s_cfg.get("mphi_min", envelope["mphi"][0]))
        mphi_max = float(s_cfg.get("mphi_max", envelope["mphi"][1]))
        n_mphi_req = int(s_cfg.get("n_mphi", 10))
        if not (envelope["mphi"][0] <= mphi_min <= envelope["mphi"][1] and envelope["mphi"][0] <= mphi_max <= envelope["mphi"][1] and mphi_min <= mphi_max):
            raise ValueError("template_grid_probe mphi range escaped envelope")
        ranges["mA"] = (min(mA_values), max(mA_values))
        ranges["lambda6"] = (min(lambda6_locals), max(lambda6_locals))
        ranges["tan_beta"] = (min(tanbeta_locals), max(tanbeta_locals))
        ranges["lambda1"] = (1.0, 1.0)
        ranges["mphi"] = (mphi_min, mphi_max)
        expanded_subruns = []
        for ma in mA_values:
            for l6 in lambda6_locals:
                expanded_subruns.append({"mA": float(ma), "lambda6": float(l6), "tanbeta_values": [float(v) for v in tanbeta_locals], "mphi_min": mphi_min, "mphi_max": mphi_max, "n_mphi": n_mphi_req, "n_lam1": 1})
        planned_points = int(len(mA_values) * len(lambda6_locals) * len(tanbeta_locals) * int(n_mphi_req) * 1)
        expected_direction = "larger_ctau_or_smaller_total_width"
    else:
        raise ValueError(f"unknown strategy {strategy}")

    tanbeta_values = [ranges["tan_beta"][0], ranges["tan_beta"][1]]
    lambda6_local_values: list[float] | None = None
    tan_beta_values: list[float] | None = None
    if strategy == "multiplicative_ctau_path":
        tanbeta_values = tanbeta_locals
        tan_beta_values = tanbeta_locals
        lambda6_local_values = lambda6_locals
        n_mphi = int(s_cfg.get("n_mphi", 10))
        n_lam1 = 1
        raw_points = calculate_raw_points_multiplicative(tanbeta_values, lambda6_local_values, n_mphi, n_lam1)
        reject = None
        if raw_points < int(budget["min_raw_points_per_subrun"]) and not bool(budget.get("operational_smoke", False)):
            reject = "insufficient_points"
        if raw_points > int(budget["max_raw_points_per_subrun"]):
            reject = "exceeds_max_points"
    elif strategy == "template_grid_probe":
        tanbeta_values = tanbeta_locals
        tan_beta_values = tanbeta_locals
        lambda6_local_values = lambda6_locals
        n_mphi = int(s_cfg.get("n_mphi", 10))
        n_lam1 = 1
        raw_points = int(planned_points or 0)
        reject = None
        if raw_points < int(budget["min_raw_points_per_subrun"]) and not bool(budget.get("operational_smoke", False)):
            reject = "insufficient_points"
        if raw_points > int(budget["max_raw_points_per_subrun"]):
            reject = "exceeds_max_points"
    else:
        n_mphi, n_lam1, raw_points, reject = enforce_raw_point_rules(
            tanbeta_values=tanbeta_values,
            n_mphi=10,
            n_lam1=100,
            mphi_varies=ranges["mphi"][0] != ranges["mphi"][1],
            budget=budget,
        )
    subcampaign_suffix = str(s_cfg.get("run_identity_suffix", "")).strip()
    subcampaign_id = f"iter{iteration:02d}_{strategy}"
    if subcampaign_suffix:
        subcampaign_id = f"{subcampaign_id}__{subcampaign_suffix}"

    return Subcampaign(
        subcampaign_id=subcampaign_id,
        parent_iteration=iteration,
        strategy=strategy,
        variable_ranges=ranges,
        tanbeta_values=tanbeta_values,
        tan_beta_values=tan_beta_values,
        lambda6_local_values=lambda6_local_values,
        mA_values=mA_values or None,
        expanded_subruns=expanded_subruns,
        command_count=len(expanded_subruns) if isinstance(expanded_subruns, list) else (len(lambda6_locals) if strategy == "multiplicative_ctau_path" else 1),
        template_name=template_name,
        planned_points=planned_points,
        n_mphi=n_mphi,
        n_lam1=n_lam1,
        estimated_raw_points=raw_points,
        reason=reject or f"deterministic_{strategy}",
        expected_output_campaign_name=f"bounded_v1_iter{iteration:02d}_{strategy}",
        safety_checks=["inside_envelope", "raw_points_bounds", "safe_layer_v1_1_gates", "triple_ok_rows_only"],
        inside_envelope=proposal_inside_envelope(ranges, envelope),
        path_step=path_step,
        tan_beta_center=tan_beta_center,
        lambda6_center=lambda6_center,
        lambda1_fixed=lambda1_fixed,
        tan_beta_factor=tan_beta_factor,
        lambda6_factor=lambda6_factor,
        objective_mode=objective_mode,
        expected_direction=expected_direction,
        lifecycle_state="PENDING",
    )


def next_thread_attempts(policy: dict[str, Any], fixed_threads: int | None = None) -> list[int]:
    if fixed_threads is not None:
        return [int(fixed_threads)]
    attempts = [int(policy.get("initial_threads", 8))]
    min_t = int(policy.get("min_threads", 1))
    backoff = int(policy.get("backoff_factor", 2))
    max_r = int(policy.get("max_retries_per_subcampaign", 4))
    while attempts[-1] > min_t and len(attempts) < max_r + 1:
        attempts.append(max(min_t, attempts[-1] // backoff))
    return attempts


def _proposal_for_subrun(sc: Subcampaign, *, run_name: str) -> dict[str, Any]:
    return {
        "proposal_id": run_name,
        "operator": "bounded_adaptive",
        "bounds": {
            "lambda6": [sc.variable_ranges["lambda6"][0], sc.variable_ranges["lambda6"][0]],
            "tan_beta": [min(sc.tanbeta_values), max(sc.tanbeta_values)],
            "lambda1": [sc.variable_ranges["lambda1"][0], sc.variable_ranges["lambda1"][1]],
            "mH": [sc.variable_ranges["mphi"][0], sc.variable_ranges["mphi"][1]],
        },
        "scan_values": {"tan_beta": sc.tanbeta_values},
        "fixed": {"sin_ba": sc.variable_ranges["sin_ba"][0], "lambda7": sc.variable_ranges["lambda7"][0], "mA": sc.variable_ranges["mA"][0], "mHp": sc.variable_ranges["mA"][0], "yukawa_type": 1},
        "resolution": {"mH": sc.n_mphi, "lambda1": sc.n_lam1},
    }


def build_orchestrator_command_for_subrun(sc: Subcampaign, runtime: dict[str, Any], threads: int, run_name: str) -> list[str]:
    proposal = _proposal_for_subrun(sc, run_name=run_name)
    rcfg = {
        "exec_path": runtime["exec_path"],
        "outdir": runtime["outdir"],
        "lake_name": runtime.get("lake_name", "dihiggs_lake"),
        "campaign": runtime["campaign"],
        "cpp_omp_threads": threads,
        "resume_scope": "run",
        "materialize": False,
    }
    cmd = build_orchestrator_command(proposal, rcfg, dry_run=False, force=False)
    cmd.insert(1, "-u")
    return cmd


def build_orchestrator_commands_for_subcampaign(sc: Subcampaign, runtime: dict[str, Any], threads: int, run_name_base: str) -> list[tuple[str, list[str], dict[str, Any]]]:
    if sc.strategy == "template_grid_probe" and sc.expanded_subruns:
        out: list[tuple[str, list[str], dict[str, Any]]] = []
        for i, sub in enumerate(sc.expanded_subruns, start=1):
            rn = f"{run_name_base}_mA{i:02d}_l6{i:02d}"
            sc_i = replace(
                sc,
                variable_ranges={
                    **sc.variable_ranges,
                    "mA": (float(sub["mA"]), float(sub["mA"])),
                    "lambda6": (float(sub["lambda6"]), float(sub["lambda6"])),
                    "mphi": (float(sub["mphi_min"]), float(sub["mphi_max"])),
                },
                tanbeta_values=[float(v) for v in sub["tanbeta_values"]],
                tan_beta_values=[float(v) for v in sub["tanbeta_values"]],
                n_mphi=int(sub["n_mphi"]),
                n_lam1=int(sub["n_lam1"]),
            )
            cmd = build_orchestrator_command_for_subrun(sc_i, runtime, threads, rn)
            out.append((rn, cmd, _proposal_for_subrun(sc_i, run_name=rn)))
        return out
    if sc.strategy != "multiplicative_ctau_path" or not sc.lambda6_local_values:
        rn = run_name_base
        cmd = build_orchestrator_command_for_subrun(sc, runtime, threads, rn)
        return [(rn, cmd, _proposal_for_subrun(sc, run_name=rn))]
    out: list[tuple[str, list[str], dict[str, Any]]] = []
    for i, l6 in enumerate(sc.lambda6_local_values, start=1):
        rn = f"{run_name_base}_l6{i:02d}"
        sc_i = replace(sc, variable_ranges={**sc.variable_ranges, "lambda6": (float(l6), float(l6))})
        cmd = build_orchestrator_command_for_subrun(sc_i, runtime, threads, rn)
        out.append((rn, cmd, _proposal_for_subrun(sc_i, run_name=rn)))
    return out


def _scan_rows_for_run(lake_dir: Path, campaign_name: str, run_name: str) -> tuple[list[dict[str, Any]], list[str]]:
    rows: list[dict[str, Any]] = []
    csvs: list[str] = []
    marker = f"campaign={campaign_name}"
    run_marker = f"run={run_name}"
    for csv_path in sorted(lake_dir.glob("**/tb_*/scan_tb_*.csv")):
        sp = str(csv_path)
        run_scoped = (run_marker in sp) or (f"/{run_name}__" in sp) or (f"/{run_name}/" in sp)
        if marker not in sp or not run_scoped:
            continue
        csvs.append(sp)
        with csv_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                rec = dict(row)
                rec["_source_csv"] = sp
                rows.append(rec)
    return rows, csvs


def _top_rows(rows: list[dict[str, Any]], *, key: str, reverse: bool, top_k: int = 10) -> list[dict[str, Any]]:
    scored: list[tuple[float, dict[str, Any]]] = []
    for row in rows:
        val = derive_br_bb(row) if key == "br_bb" else finite_float(row.get(key))
        if val is None:
            continue
        scored.append((val, row))
    scored.sort(key=lambda x: x[0], reverse=reverse)
    out: list[dict[str, Any]] = []
    for val, row in scored[:top_k]:
        canon = canonicalize_physical_row(row)
        out.append(
            {
                "value": val,
                "lambda6": canon["lambda6"],
                "tan_beta": canon["tan_beta"],
                "mphi": canon["mphi"],
                "lambda1": canon["lambda1"],
                "source_csv": row.get("_source_csv"),
            }
        )
    return out


def attempt_isolated_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    tok = [r for r in rows if compute_triple_ok(r)]
    br_gaga_vals = [v for v in (finite_float(r.get("br_gaga")) for r in tok) if v is not None]
    br_bb_vals = [v for v in (derive_br_bb(r) for r in tok) if v is not None]
    width_vals = [v for v in (finite_float(r.get("total_width")) for r in tok) if v is not None]
    out = {
        "isolated_n_rows": len(rows),
        "isolated_n_triple_ok": len(tok),
        "isolated_triple_ok_rate": (len(tok) / len(rows)) if rows else 0.0,
        "isolated_br_gaga_q50": quantile(br_gaga_vals, 0.50),
        "isolated_br_gaga_q90": quantile(br_gaga_vals, 0.90),
        "isolated_br_gaga_q95": quantile(br_gaga_vals, 0.95),
        "isolated_br_bb_q50": quantile(br_bb_vals, 0.50),
        "isolated_br_bb_q90": quantile(br_bb_vals, 0.90),
        "isolated_total_width_q50": quantile(width_vals, 0.50),
        "isolated_best_rows_top10_br_gaga": _top_rows(tok, key="br_gaga", reverse=True, top_k=10),
        "isolated_best_rows_top10_low_br_bb": _top_rows(tok, key="br_bb", reverse=False, top_k=10),
        "ranking_basis": "triple_ok_rows_only",
    }
    out.update(_ctau_summary(rows))
    return out


def _cumulative_summary(accepted_isolated: list[dict[str, Any]]) -> dict[str, Any]:
    if not accepted_isolated:
        return {
            "cumulative_n_rows": 0,
            "cumulative_n_triple_ok": 0,
            "cumulative_triple_ok_rate": 0.0,
            "cumulative_br_gaga_q95": None,
            "cumulative_br_bb_q50": None,
        }
    n_rows = sum(int(s.get("isolated_n_rows", 0) or 0) for s in accepted_isolated)
    n_tok = sum(int(s.get("isolated_n_triple_ok", 0) or 0) for s in accepted_isolated)
    return {
        "cumulative_n_rows": n_rows,
        "cumulative_n_triple_ok": n_tok,
        "cumulative_triple_ok_rate": (n_tok / n_rows) if n_rows else 0.0,
        "cumulative_br_gaga_q95": accepted_isolated[-1].get("isolated_br_gaga_q95"),
        "cumulative_br_bb_q50": accepted_isolated[-1].get("isolated_br_bb_q50"),
    }


def _delta_summary(isolated: dict[str, Any], cumulative: dict[str, Any]) -> dict[str, Any]:
    return {
        "delta_n_rows_vs_cumulative": int(isolated.get("isolated_n_rows", 0) or 0) - int(cumulative.get("cumulative_n_rows", 0) or 0),
        "delta_n_triple_ok_vs_cumulative": int(isolated.get("isolated_n_triple_ok", 0) or 0) - int(cumulative.get("cumulative_n_triple_ok", 0) or 0),
        "delta_triple_ok_rate_vs_cumulative": float(isolated.get("isolated_triple_ok_rate", 0.0) or 0.0) - float(cumulative.get("cumulative_triple_ok_rate", 0.0) or 0.0),
    }


def _anchor_score(row: dict[str, Any]) -> float:
    br_gaga = finite_float(row.get("br_gaga")) or 0.0
    br_bb = derive_br_bb(row)
    bb = br_bb if br_bb is not None else 1.0
    return br_gaga / math.sqrt(max(bb, 1e-12))


def _first_finite(row: dict[str, Any], keys: tuple[str, ...]) -> float | None:
    for key in keys:
        val = finite_float(row.get(key))
        if val is not None:
            return val
    return None


def canonicalize_physical_row(row: dict[str, Any]) -> dict[str, float | None]:
    total_width = _first_finite(row, ("total_width",))
    width_bb = _first_finite(row, ("width_bb",))
    br_bb = _first_finite(row, ("br_bb",))
    if br_bb is None and width_bb is not None and total_width is not None and abs(total_width) > 0.0:
        br_bb = width_bb / total_width
    ctau_m = _first_finite(row, ("ctau", "c_tau", "decay_length", "lifetime"))
    if ctau_m is None and total_width is not None and total_width > 0.0 and math.isfinite(total_width):
        ctau_m = HBARC_GEV_M / total_width
    return {
        "lambda6": _first_finite(row, ("lambda6",)),
        "tan_beta": _first_finite(row, ("tan_beta",)),
        "mphi": _first_finite(row, ("mphi", "m_phi", "mH")),
        "lambda1": _first_finite(row, ("lambda1", "lam1")),
        "mA": _first_finite(row, ("mA",)),
        "br_gaga": _first_finite(row, ("br_gaga",)),
        "br_bb": br_bb,
        "total_width": total_width,
        "ctau_m": ctau_m,
    }


def anchor_is_complete(anchor: dict[str, Any] | None) -> bool:
    if not isinstance(anchor, dict):
        return False
    required = ("lambda6", "tan_beta", "mphi", "lambda1", "br_gaga", "br_bb", "total_width")
    return all(finite_float(anchor.get(k)) is not None for k in required)


def select_best_anchor_row(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    tok = [r for r in rows if compute_triple_ok(r)]
    if not tok:
        return None
    best = max(tok, key=_anchor_score)
    canon = canonicalize_physical_row(best)
    return {
        "lambda6": canon["lambda6"],
        "tan_beta": canon["tan_beta"],
        "mphi": canon["mphi"],
        "lambda1": canon["lambda1"],
        "br_gaga": canon["br_gaga"],
        "br_bb": canon["br_bb"],
        "total_width": canon["total_width"],
        "anchor_score": _anchor_score(best),
        "source_csv": best.get("_source_csv"),
        "ranking_basis": "triple_ok_rows_only",
    }


def normalize(x: float | None, lo: float, hi: float) -> float:
    if x is None or not math.isfinite(x) or hi <= lo:
        return 0.0
    return max(0.0, min(1.0, (x - lo) / (hi - lo)))


def compute_subcampaign_score(summary: dict[str, Any], weights: dict[str, float] | None = None) -> float:
    w = {"w_gaga": 1.0, "w_bb": 1.0, "w_ok": 1.0, "w_health": 1.0}
    if weights:
        w.update({k: float(v) for k, v in weights.items()})
    return (
        w["w_gaga"] * normalize(finite_float(summary.get("isolated_br_gaga_q95")), 0.0, 1e-3)
        - w["w_bb"] * normalize(finite_float(summary.get("isolated_br_bb_q50")), 0.0, 1.0)
        + w["w_ok"] * normalize(finite_float(summary.get("isolated_triple_ok_rate")), 0.0, 1.0)
        - w["w_health"] * float(summary.get("health_penalty", 0.0) or 0.0)
    )


def should_retry_due_to_gate(stop_report: dict[str, Any], policy: dict[str, Any]) -> bool:
    if not bool(policy.get("retry_on_gsl_abort", True)):
        return False
    gates = stop_report.get("gates", {}) if isinstance(stop_report, dict) else {}
    failures = gates.get("failures", []) if isinstance(gates, dict) else []
    for f in failures:
        if not isinstance(f, dict):
            continue
        code = str(f.get("code", ""))
        exp = str(f.get("explanation", "")).lower()
        if code in {"gate_failed_due_to_gsl_signature", "gate_failed_due_to_header_only_csv", "gate_failed_due_to_health"}:
            return True
        if "sigabrt" in exp or "gsl" in exp or "header-only" in exp:
            return True
    return False


def select_strategy(iteration: int, accepted_summaries: list[dict[str, Any]]) -> str:
    if iteration == 1:
        return "box_partition"
    if not accepted_summaries:
        return "coordinate_sweep"
    best = max(accepted_summaries, key=lambda x: float(x.get("score", -1e9)))
    if int(best.get("isolated_n_triple_ok", best.get("n_triple_ok", 0)) or 0) <= 0:
        return "coordinate_sweep"
    if bool(best.get("local_refine_unavailable", False)):
        return "coordinate_sweep"
    return "local_refine"


def _append_jsonl(path: Path, obj: dict[str, Any]) -> None:
    with path.open("a", encoding="utf-8") as h:
        h.write(json.dumps(obj, sort_keys=True) + "\n")


def _write_progress(state: dict[str, Any], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "adaptive_state.json").write_text(json.dumps(state, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md = ["# Bounded Adaptive Search State", "", f"Generated: {state.get('generated_utc')}", f"Total points executed: {state.get('total_points_executed', 0)}"]
    for it in state.get("iterations", []):
        md.append(f"- Iteration {it['iteration']}: {it['status']} ({it['strategy']})")
        if isinstance(it.get("best_row_anchor"), dict):
            md.append(f"  - best_row_anchor: {json.dumps(it.get('best_row_anchor'), sort_keys=True)}")
        if isinstance(it.get("summary"), dict):
            md.append(f"  - isolated_summary: {json.dumps(it['summary'].get('isolated_summary', {}), sort_keys=True)}")
            md.append(f"  - cumulative_summary: {json.dumps(it['summary'].get('cumulative_summary', {}), sort_keys=True)}")
    (outdir / "adaptive_state.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    dlines: list[str] = ["# Iteration decisions"]
    for it in state.get("iterations", []):
        dlines.append(f"## Iteration {it['iteration']}")
        dlines.append(f"- Strategy: {it['strategy']}")
        dlines.append(f"- Status: {it['status']}")
        dlines.append(f"- Decision: {json.dumps(it.get('decision', {}), sort_keys=True)}")
        dlines.append(f"- Envelope clipping: {json.dumps(it.get('range_clipping', {}), sort_keys=True)}")
    (outdir / "decisions.md").write_text("\n".join(dlines) + "\n", encoding="utf-8")
    (outdir / "bounded_adaptive_report.md").write_text("# Bounded Adaptive Search v1.1 report\n\nContains isolated_summary, cumulative_summary, delta_summary, anchor selection, envelope clipping, budget remaining, and effective threads per iteration.\n", encoding="utf-8")


def _write_incremental_progress(outdir: Path, record: dict[str, Any]) -> None:
    progress_jsonl = outdir / "progress.jsonl"
    _append_jsonl(progress_jsonl, record)

    lines = ["# Incremental progress", ""]
    for raw in progress_jsonl.read_text(encoding="utf-8").splitlines()[-50:]:
        if not raw.strip():
            continue
        rec = json.loads(raw)
        lines.append(f"## {rec.get('timestamp', 'unknown')} | {rec.get('subcampaign_id', 'na')}")
        lines.append(f"- path_step: {rec.get('path_step')}")
        lines.append(f"- effective_threads: {rec.get('effective_threads')}")
        lines.append(f"- gates_passed: {rec.get('gates_passed')}")
        lines.append(f"- n_rows: {rec.get('n_rows')} | n_triple_ok: {rec.get('n_triple_ok')} | triple_ok_rate: {rec.get('triple_ok_rate')}")
        lines.append(f"- best_ctau_m: {rec.get('best_ctau_m')} | ctau_q95: {rec.get('ctau_q95')}")
        lines.append(f"- best_point: {json.dumps(rec.get('best_point', {}), sort_keys=True)}")
        lines.append(f"- recommendation: {rec.get('current_recommendation')}")
        lines.append(f"- budget_remaining_points: {rec.get('budget_remaining_points')}")
        lines.append("")
    (outdir / "progress.md").write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def _safe_layer_contract(runtime: dict[str, Any], campaign: str, run_name: str, proposal: dict[str, Any]) -> dict[str, Any]:
    return {
        "enable_proposal_generation": False,
        "allow_broad_scans": False,
        "modify_physics_semantics": False,
        "modify_scoring_semantics": False,
        "modify_triple_ok_definition": False,
        "use_all_row_br_metrics_as_physics_claims": False,
        "investigate_gsl_root_cause": False,
        "require_cpp_omp_threads_1": True,
        "analysis": {"triple_ok_only": True},
        "gates": {"min_triple_ok_rate": 0.0, "operational_only": True},
        "evidence_mode": "existing_artifacts",
        "runtime": {"threads": 1, "exec_path": runtime["exec_path"], "outdir": runtime["outdir"], "lake_name": runtime.get("lake_name", "dihiggs_lake"), "dry_run": True, "force": False},
        "campaign": {"name": campaign, "max_runs": 1, "expansion_max_runs": 0, "proposals": [proposal]},
        "run_filter": {"run_name": run_name},
    }


def run_bounded_search(contract: dict[str, Any], *, execute: bool, plan_only: bool, output_dir: str | Path, max_iterations_override: int | None = None, resume_state: dict[str, Any] | None = None) -> dict[str, Any]:
    envelope = validate_envelope(dict(contract["search_envelope"]))
    budget = validate_budget(dict(contract["budget"]))
    runtime = dict(contract.get("runtime", {}))
    runtime.setdefault("exec_path", "/home/fabi/dihiggs/dihiggs/app/PhysScanWithFixings")
    runtime.setdefault("outdir", "scripts/out")
    runtime.setdefault("lake_name", "dihiggs_lake")
    runtime.setdefault("campaign", f"bounded_adaptive_{datetime.now().strftime('%Y%m%d_%H%M%S')}")

    weights = dict(contract.get("score_weights", {"w_gaga": 1.0, "w_bb": 1.0, "w_ok": 1.0, "w_health": 1.0}))
    refine_cfg = dict(contract.get("local_refine", {}))
    thread_policy = dict(contract.get("runtime_thread_policy", {}))
    fixed_threads = contract.get("fixed_threads")
    allow_fixed_retry = bool(contract.get("allow_fixed_thread_retries", False))

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    for name in ["attempts.jsonl", "commands.jsonl", "subcampaigns.jsonl"]:
        (outdir / name).touch()

    started = datetime.now(timezone.utc)
    deadline = None
    if budget["max_walltime_minutes"] is not None:
        deadline = started + timedelta(minutes=int(budget["max_walltime_minutes"]))

    state = dict(resume_state or {"iterations": [], "total_points_executed": 0, "failed_subcampaigns": 0, "campaign": runtime["campaign"]})
    accepted_summaries: list[dict[str, Any]] = []
    anchor = {k: (v[0] + v[1]) / 2 for k, v in envelope.items()}
    strategy_cfg = dict(contract.get("strategy", {})) if isinstance(contract.get("strategy"), dict) else {}
    configured_strategy = str(strategy_cfg.get("name", "")).strip()

    max_iter = int(max_iterations_override or budget["max_iterations"])
    for iteration in range(len(state["iterations"]) + 1, max_iter + 1):
        remaining_minutes = None
        if deadline is not None:
            remaining_minutes = (deadline - datetime.now(timezone.utc)).total_seconds() / 60.0
            if remaining_minutes < float(budget["min_remaining_minutes_to_start"]):
                break

        strategy = configured_strategy or select_strategy(iteration, accepted_summaries)
        if strategy == "multiplicative_ctau_path":
            strategy_cfg["path_step"] = iteration - 1
        try:
            sc = build_subcampaign(iteration, strategy, envelope, anchor, budget, refine_cfg=refine_cfg, strategy_cfg=strategy_cfg)
        except ValueError as exc:
            state["iterations"].append({"iteration": iteration, "strategy": strategy, "status": "stopped", "reason": str(exc)})
            _write_progress(state, outdir)
            break
        _append_jsonl(outdir / "subcampaigns.jsonl", sc.to_dict())

        if not sc.inside_envelope or sc.reason in {"insufficient_points", "exceeds_max_points"}:
            state["iterations"].append({"iteration": iteration, "strategy": strategy, "status": "rejected", "reason": sc.reason, "subcampaign": sc.to_dict()})
            _write_progress(state, outdir)
            continue

        if state["total_points_executed"] + sc.estimated_raw_points > int(budget["max_total_points"]):
            break

        thread_plan = next_thread_attempts(thread_policy, int(fixed_threads) if fixed_threads is not None else None)
        if fixed_threads is not None and not allow_fixed_retry:
            thread_plan = thread_plan[:1]

        iter_rec: dict[str, Any] = {
            "iteration": iteration,
            "strategy": strategy,
            "status": "technical_failure",
            "subcampaign": sc.to_dict(),
            "attempts": [],
            "path_step": sc.path_step,
            "tan_beta_factor": sc.tan_beta_factor,
            "lambda6_factor": sc.lambda6_factor,
            "lambda1_fixed": sc.lambda1_fixed,
            "ctau_target_window_m": strategy_cfg.get("target_ctau_window_m"),
            "objective_mode": sc.objective_mode,
        }
        accepted_summary: dict[str, Any] | None = None

        for aidx, threads in enumerate(thread_plan, start=1):
            run_name = f"{sc.subcampaign_id}_t{threads}_a{aidx}"
            run_cmds = build_orchestrator_commands_for_subcampaign(sc, runtime, threads, run_name)
            rc, stdout, stderr = 0, "", ""
            all_rows: list[dict[str, Any]] = []
            all_csv_paths: list[str] = []
            stop_report = {"gates": {"passed": True, "failures": []}, "production_validation": False}
            for ridx, (run_name_i, cmd, proposal) in enumerate(run_cmds, start=1):
                _append_jsonl(outdir / "commands.jsonl", {"iteration": iteration, "subcampaign_id": sc.subcampaign_id, "attempt_index": aidx, "run_index": ridx, "command": cmd, "threads": threads})
                if execute and not plan_only:
                    proc = subprocess.run(cmd, cwd=".", text=True, capture_output=True, timeout=int(budget["max_per_subrun_minutes"]) * 60)
                    rc_i, so_i, se_i = int(proc.returncode), proc.stdout, proc.stderr
                    rc = max(rc, rc_i)
                    stdout += so_i
                    stderr += se_i
                safe_contract = _safe_layer_contract(runtime, state["campaign"], run_name_i, proposal)
                safe = run_safe_campaign(safe_contract, workdir=".", execute=False)
                stop_report = safe["stop_report"]
                rows_i, csv_paths_i = _scan_rows_for_run(Path(runtime["outdir"]) / runtime.get("lake_name", "dihiggs_lake"), state["campaign"], run_name_i)
                all_rows.extend(rows_i)
                all_csv_paths.extend(csv_paths_i)
            stop_path = outdir / f"stop_report_iter{iteration:02d}_attempt{aidx}.json"
            stop_path.write_text(json.dumps(stop_report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

            campaign_state_path = outdir / f"campaign_state_iter{iteration:02d}_attempt{aidx}.json"
            rerun_manifest_path = outdir / f"rerun_manifest_iter{iteration:02d}_attempt{aidx}.json"
            campaign_state_path.write_text(json.dumps(safe["campaign_state"], indent=2, sort_keys=True) + "\n", encoding="utf-8")
            rerun_manifest_path.write_text(json.dumps(safe["rerun_manifest"], indent=2, sort_keys=True) + "\n", encoding="utf-8")

            rows, csv_paths = all_rows, all_csv_paths
            isolated = attempt_isolated_summary(rows)
            accepted = bool(stop_report.get("gates", {}).get("passed", False)) and rc == 0
            isolated["health_penalty"] = 0.0 if accepted else 1.0
            isolated["score"] = compute_subcampaign_score(isolated, weights)

            attempt_rec = {
                "iteration": iteration,
                "subcampaign_id": sc.subcampaign_id,
                "attempt_index": aidx,
                "run_name": run_name,
                "command": cmd,
                "threads": threads,
                "return_code": rc,
                "run_dir": str(Path(runtime["outdir"]) / runtime.get("lake_name", "dihiggs_lake") / state["campaign"]),
                "artifact_roots_used_for_validation": {
                    "run_filter": {"run_name": run_name},
                    "attempt_csv_paths": csv_paths,
                    "campaign_state_path": str(campaign_state_path),
                    "rerun_manifest_path": str(rerun_manifest_path),
                    "safe_automation_stop_report_path": str(stop_path),
                },
                "gates_passed": bool(stop_report.get("gates", {}).get("passed", False)),
                "production_validation": bool(stop_report.get("production_validation", False)),
                "failure_reason_codes": [f.get("code") for f in stop_report.get("gates", {}).get("failures", []) if isinstance(f, dict)],
                "isolated_n_rows": isolated["isolated_n_rows"],
                "isolated_n_triple_ok": isolated["isolated_n_triple_ok"],
                "isolated_triple_ok_rate": isolated["isolated_triple_ok_rate"],
                "accepted_for_physics_analysis": accepted,
                "stdout_tail": stdout[-5000:],
                "stderr_tail": stderr[-5000:],
            }
            iter_rec["attempts"].append(attempt_rec)
            _append_jsonl(outdir / "attempts.jsonl", attempt_rec)

            progress_rec = {
                "timestamp": utc_now_iso(),
                "iteration": iteration,
                "subcampaign_id": sc.subcampaign_id,
                "path_step": sc.path_step,
                "effective_threads": threads,
                "gates_passed": bool(stop_report.get("gates", {}).get("passed", False)),
                "n_rows": isolated.get("isolated_n_rows"),
                "n_triple_ok": isolated.get("isolated_n_triple_ok"),
                "triple_ok_rate": isolated.get("isolated_triple_ok_rate"),
                "best_ctau_m": isolated.get("ctau_m_max"),
                "ctau_q95": isolated.get("ctau_m_q95"),
                "best_point": select_best_anchor_row(rows) or {},
                "current_recommendation": "continue" if accepted else "retry_or_switch_axis",
                "budget_remaining_points": int(budget["max_total_points"]) - int(state["total_points_executed"]),
            }
            _write_incremental_progress(outdir, progress_rec)

            if accepted:
                iter_rec["status"] = "accepted"
                iter_rec["accepted_attempt_index"] = aidx
                accepted_summary = isolated
                break

            if aidx < len(thread_plan) and should_retry_due_to_gate(stop_report, thread_policy):
                iter_rec["status"] = "retrying"
                _write_progress(state | {"iterations": state["iterations"] + [iter_rec], "generated_utc": utc_now_iso()}, outdir)
                continue

        if accepted_summary is None:
            state["failed_subcampaigns"] += 1
            iter_rec["decision"] = {"continue": state["failed_subcampaigns"] < int(budget["max_failed_subcampaigns"]), "reason": "technical_failure"}
        else:
            accepted_summary["excluded_failed_attempts"] = True
            accepted_summaries.append(accepted_summary)
            cumulative = _cumulative_summary(accepted_summaries)
            iter_rec["summary"] = {
                "isolated_summary": accepted_summary,
                "cumulative_summary": cumulative,
                "delta_summary": _delta_summary(accepted_summary, cumulative),
            }
            state["total_points_executed"] += sc.estimated_raw_points

            best_row_anchor = select_best_anchor_row(rows)
            iter_rec["best_row_anchor"] = best_row_anchor
            if best_row_anchor is None:
                accepted_summary["local_refine_unavailable"] = True
                iter_rec["decision"] = {"continue": True, "reason": "accepted_no_triple_ok_no_local_refine"}
            elif not anchor_is_complete(best_row_anchor):
                accepted_summary["local_refine_unavailable"] = True
                iter_rec["decision"] = {
                    "continue": True,
                    "reason": "local_refine_unavailable",
                    "reason_code": "missing_anchor_fields",
                    "fallback_strategy": "coordinate_sweep",
                }
            else:
                accepted_summary["local_refine_unavailable"] = False
                clip_info: dict[str, bool] = {}
                for k in ("lambda6", "tan_beta", "mphi", "lambda1"):
                    val = float(best_row_anchor[k])
                    clip_val = min(max(val, envelope[k][0]), envelope[k][1])
                    clip_info[k] = (clip_val != val)
                    anchor[k] = clip_val
                iter_rec["range_clipping"] = clip_info
                iter_rec["decision"] = {"continue": True, "reason": "accepted_best_row_anchor_selected"}

        state["iterations"].append(iter_rec)
        state["generated_utc"] = utc_now_iso()
        _write_progress(state, outdir)
        if state["failed_subcampaigns"] >= int(budget["max_failed_subcampaigns"]):
            break

    state["generated_utc"] = utc_now_iso()
    _write_progress(state, outdir)
    return state


def _cli(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Bounded Adaptive Search v1")
    ap.add_argument("--contract", required=True)
    ap.add_argument("--output-dir", default="scripts/out/bounded_adaptive_v1_demo")
    ap.add_argument("--execute", action="store_true")
    ap.add_argument("--plan-only", action="store_true")
    ap.add_argument("--max-iterations", type=int)
    ap.add_argument("--resume")
    args = ap.parse_args(argv)

    contract = load_json(args.contract)
    resume = load_json(args.resume) if args.resume else None
    state = run_bounded_search(contract, execute=args.execute, plan_only=args.plan_only or not args.execute, output_dir=args.output_dir, max_iterations_override=args.max_iterations, resume_state=resume)
    print(json.dumps({"ok": True, "iterations": len(state.get("iterations", [])), "output_dir": str(args.output_dir)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
