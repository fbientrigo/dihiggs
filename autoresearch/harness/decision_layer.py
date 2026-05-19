from __future__ import annotations

import math
from typing import Any


def _f(x: Any) -> float | None:
    try:
        v = float(x)
    except Exception:
        return None
    return v if math.isfinite(v) else None


def reachability_analysis(
    *,
    path_trend_summary: dict[str, Any] | None,
    sensitivity_summary: dict[str, Any] | None,
    envelope: dict[str, list[float]] | None = None,
    target_ctau_m: float = 1.0,
    short_term_factor: float = 2.0,
) -> dict[str, Any]:
    p = path_trend_summary or {}
    s = sensitivity_summary or {}

    best_vals = [
        _f(s.get("new_best_ctau_m")),
        _f(s.get("best_ctau_m")),
        _f(p.get("ctau_m_max")),
        _f(p.get("best_ctau_m")),
    ]
    best_ctau_m = max([v for v in best_vals if v is not None], default=None)
    if best_ctau_m is None or best_ctau_m <= 0:
        return {
            "best_ctau_m": None,
            "target_ctau_m": target_ctau_m,
            "improvement_factor_to_target": None,
            "short_term_target_ctau_m": None,
            "estimated_steps_to_target": None,
            "recommendation": "stop_no_progress",
            "rationale": "No valid best_ctau_m in summaries.",
        }

    slope = _f(p.get("log10_ctau_slope_per_step"))
    if slope is None:
        slope = _f(p.get("slope_dex_per_step"))

    steps_to_target = None
    if slope is not None and slope > 0:
        steps_to_target = max(0.0, math.log10(target_ctau_m / best_ctau_m) / slope)

    steps_available = None
    projected_boundary = None
    if envelope and p:
        tb = envelope.get("tan_beta") if isinstance(envelope.get("tan_beta"), list) else None
        l6 = envelope.get("lambda6") if isinstance(envelope.get("lambda6"), list) else None
        tb_now = _f(p.get("tan_beta_current") or p.get("tan_beta_center") or p.get("tan_beta"))
        l6_now = _f(p.get("lambda6_current") or p.get("lambda6_center") or p.get("lambda6"))
        tb_factor = _f(p.get("tan_beta_factor") or p.get("tan_beta_step_factor"))
        l6_factor = _f(p.get("lambda6_factor") or p.get("lambda6_step_factor"))

        candidates: list[float] = []
        if tb and tb_now and tb_factor and tb_factor > 1.0 and tb[1] > tb_now:
            candidates.append(math.log(tb[1] / tb_now, tb_factor))
        if l6 and l6_now and l6_factor and l6_factor > 1.0 and l6[1] > l6_now:
            candidates.append(math.log(l6[1] / l6_now, l6_factor))
        if candidates:
            steps_available = max(0.0, min(candidates))
            if slope and slope > 0:
                projected_boundary = best_ctau_m * (10 ** (slope * steps_available))

    near_target = best_ctau_m * short_term_factor
    improv_factor = target_ctau_m / best_ctau_m

    recommendation = "continue_path"
    rationale = "Path continuation is currently viable."
    if steps_to_target is not None and steps_available is not None and steps_to_target > steps_available:
        recommendation = "widen_envelope"
        rationale = "Projected steps to target exceed remaining path steps before envelope exit."
    if projected_boundary is not None and projected_boundary < near_target:
        recommendation = "refine_mphi"
        rationale = "Path likely cannot reach short-term 2x before boundary; prioritize new lever."

    return {
        "best_ctau_m": best_ctau_m,
        "target_ctau_m": target_ctau_m,
        "improvement_factor_to_target": improv_factor,
        "short_term_target_ctau_m": near_target,
        "log10_slope_per_step": slope,
        "estimated_steps_to_target": steps_to_target,
        "steps_available_before_envelope_exit": steps_available,
        "projected_ctau_at_envelope_boundary": projected_boundary,
        "recommendation": recommendation,
        "rationale": rationale,
    }


def next_axis_selector(
    *,
    per_mA_summary: list[dict[str, Any]] | None,
    global_top_by_ctau: list[dict[str, Any]] | None,
    path_trend_summary: dict[str, Any] | None,
    triple_ok_rate: float | None,
    mphi_bounds: tuple[float, float] | None,
    mA_bounds: tuple[float, float] | None,
    target_ctau_m: float,
    best_ctau_m: float,
) -> dict[str, Any]:
    tops = global_top_by_ctau or []
    lower_edge_share = 0.0
    upper_mA_share = 0.0
    if tops and mphi_bounds:
        lo = float(mphi_bounds[0])
        lower_edge_share = sum(1 for r in tops if _f(r.get("m_phi")) is not None and _f(r.get("m_phi")) <= lo + 5.0) / len(tops)
    if tops and mA_bounds:
        hi = float(mA_bounds[1])
        upper_mA_share = sum(1 for r in tops if _f(r.get("mA") or r.get("mA_path")) is not None and _f(r.get("mA") or r.get("mA_path")) >= hi - 1e-9) / len(tops)

    slope = _f((path_trend_summary or {}).get("log10_ctau_slope_per_step") or (path_trend_summary or {}).get("slope_dex_per_step"))
    cannot_reach = False
    if slope and slope > 0 and best_ctau_m > 0:
        steps = math.log10(target_ctau_m / best_ctau_m) / slope
        cannot_reach = steps > 20

    if triple_ok_rate is not None and triple_ok_rate < 0.80:
        return {
            "action": "backoff",
            "axis": "stability",
            "proposed_envelope_patch": {},
            "reason": "triple_ok_rate collapsed; prioritize stability before further refinement.",
            "confidence": 0.9,
            "expected_min_points": 300,
            "safety_notes": ["Do not claim physics from failed attempts.", "Keep bounded scope."],
        }

    if lower_edge_share >= 0.6:
        return {
            "action": "refine",
            "axis": "mphi",
            "proposed_envelope_patch": {"mphi": "shift lower within declared bounds"},
            "reason": "Top ctau rows pile up at lower mphi edge.",
            "confidence": 0.8,
            "expected_min_points": 600,
            "safety_notes": ["Remain inside declared mphi envelope."],
        }

    if upper_mA_share >= 0.6:
        return {
            "action": "refine",
            "axis": "mA",
            "proposed_envelope_patch": {"mA": "densify high-mA region within bounds"},
            "reason": "Top ctau rows pile up at upper mA edge.",
            "confidence": 0.78,
            "expected_min_points": 600,
            "safety_notes": ["Do not widen mA without explicit approval."],
        }

    if cannot_reach:
        return {
            "action": "switch",
            "axis": "mixed_refine",
            "proposed_envelope_patch": {"focus": ["mphi", "mA"]},
            "reason": "Positive path slope exists but target is unreachable in practical remaining steps.",
            "confidence": 0.7,
            "expected_min_points": 800,
            "safety_notes": ["Do not continue path-only by default."],
        }

    return {
        "action": "stop",
        "axis": "control_box",
        "proposed_envelope_patch": {},
        "reason": "No clear improvement axis signal.",
        "confidence": 0.55,
        "expected_min_points": 200,
        "safety_notes": ["Run small control box or stop."],
    }
