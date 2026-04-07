#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import adaptive_artifacts
import orchestrate_scans


logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] BranchContinuationV2: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Path defaults
# -----------------------------------------------------------------------------
def _default_exec_path() -> str:
    sibling = Path(__file__).with_name("PhysScanWithFixings")
    if sibling.exists():
        return str(sibling)
    return "./PhysScanWithFixings"


# -----------------------------------------------------------------------------
# Data models
# -----------------------------------------------------------------------------
@dataclass(frozen=True)
class FixedParams:
    mA: float
    sin_ba: float
    lambda7: float


@dataclass(frozen=True)
class SearchWindow:
    mphi_min: float
    mphi_max: float
    n_mphi: int
    lam1_min: float
    lam1_max: float
    n_lam1: int
    label: str
    expansion_index: int


@dataclass(frozen=True)
class ExternalPoint:
    tanbeta: float
    lambda6: float
    source: str
    approx_level_shape: float | None
    approx_level_ratio: float | None
    note: str | None = None


@dataclass(frozen=True)
class CandidatePoint:
    m_phi: float
    lam1: float
    lam2: float
    m12: float
    width_bb: float
    width_gaga: float
    br_gaga: float
    total_width: float
    positivity_ok: int
    unitarity_ok: int
    perturbativity_ok: int
    csv_path: str
    row_index: int
    raw: dict[str, Any]

    @property
    def triple_ok(self) -> bool:
        return self.positivity_ok == 1 and self.unitarity_ok == 1 and self.perturbativity_ok == 1

    @property
    def r_exact(self) -> float | None:
        if self.width_bb <= 0.0 or self.width_gaga < 0.0:
            return None
        return float(self.width_gaga / self.width_bb)


@dataclass(frozen=True)
class BranchState:
    lambda6: float
    tanbeta: float
    m_phi: float
    lam1: float
    lam2: float
    m12: float
    width_bb: float
    width_gaga: float
    br_gaga: float
    total_width: float
    r_exact: float | None
    approx_level_shape: float | None
    approx_level_ratio: float | None
    csv_path: str
    row_index: int
    selection_reason: str


@dataclass(frozen=True)
class AttemptResult:
    status: str
    command: list[str]
    run_dir: str | None
    csv_paths: list[str]
    n_candidates_total: int
    n_candidates_triple_ok: int
    selected: BranchState | None
    selection_score: float | None
    window: SearchWindow
    elapsed_sec: float
    notes: list[str]


@dataclass(frozen=True)
class DistanceWeights:
    mphi: float
    lam1: float
    lam2: float
    m12: float
    width_bb: float
    width_gaga: float
    br_gaga: float
    total_width: float
    r_exact: float


@dataclass(frozen=True)
class TrackDefinition:
    track_id: str
    label: str
    mode: str
    seed_point: ExternalPoint
    forward_points: list[ExternalPoint]
    backward_points: list[ExternalPoint]
    metadata: dict[str, Any]


# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------
def _float_or_raise(value: str, field: str) -> float:
    try:
        out = float(value)
    except Exception as exc:
        raise ValueError(f"Invalid float for {field}: {value!r}") from exc
    if not math.isfinite(out):
        raise ValueError(f"Non-finite float for {field}: {value!r}")
    return out


def _parse_csv_floats(text: str | None) -> list[float]:
    if text is None:
        return []
    values: list[float] = []
    for raw in text.split(","):
        s = raw.strip()
        if not s:
            continue
        values.append(_float_or_raise(s, "csv list"))
    return values


def _parse_seed_pairs(text: str | None) -> list[tuple[float, float]]:
    if text is None:
        return []
    pairs: list[tuple[float, float]] = []
    normalized = text.replace(";", ",")
    for raw in normalized.split(","):
        s = raw.strip()
        if not s:
            continue
        if ":" not in s:
            raise ValueError(f"Invalid seed pair {s!r}; expected tb:lambda6")
        tb_raw, l6_raw = s.split(":", 1)
        tb = _float_or_raise(tb_raw.strip(), "seed tb")
        l6 = _float_or_raise(l6_raw.strip(), "seed lambda6")
        pairs.append((tb, l6))
    return pairs


def _parse_manual_points(text: str | None) -> list[tuple[float, float]]:
    return _parse_seed_pairs(text)


def _make_range(min_value: float, max_value: float, count: int, scale: str) -> list[float]:
    if count <= 0:
        raise ValueError("count must be > 0")
    if count == 1:
        return [float(min_value)]
    if scale == "linear":
        step = (max_value - min_value) / float(count - 1)
        return [float(min_value + i * step) for i in range(count)]
    if scale == "log":
        if min_value <= 0.0 or max_value <= 0.0:
            raise ValueError("log scale requires positive min/max")
        log_min = math.log10(min_value)
        log_max = math.log10(max_value)
        step = (log_max - log_min) / float(count - 1)
        return [float(10 ** (log_min + i * step)) for i in range(count)]
    raise ValueError(f"Unknown scale: {scale}")


def _resolve_value_grid(
    *,
    explicit_values: str | None,
    min_value: float | None,
    max_value: float | None,
    count: int | None,
    scale: str,
    field_name: str,
) -> list[float]:
    values = _parse_csv_floats(explicit_values)
    if values:
        return values
    if min_value is None or max_value is None or count is None:
        raise ValueError(
            f"{field_name}: provide either --{field_name}-values or the trio --{field_name}-min/--{field_name}-max/--n-{field_name}"
        )
    return _make_range(min_value=float(min_value), max_value=float(max_value), count=int(count), scale=scale)


def _jsonable(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, tuple):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    return obj


def _json_dump(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(_jsonable(payload), f, indent=2, sort_keys=True)


def _jsonl_append(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(_jsonable(payload), sort_keys=True))
        f.write("\n")


def _centered_window(
    *,
    center: float,
    half_span: float,
    global_min: float,
    global_max: float,
) -> tuple[float, float]:
    lo = center - half_span
    hi = center + half_span
    if lo < global_min:
        shift = global_min - lo
        lo += shift
        hi += shift
    if hi > global_max:
        shift = hi - global_max
        lo -= shift
        hi -= shift
    lo = max(lo, global_min)
    hi = min(hi, global_max)
    if hi < lo:
        hi = lo
    return float(lo), float(hi)


def _make_window_around_state(
    *,
    prev: BranchState,
    mphi_half_span: float,
    lam1_half_span: float,
    mphi_global_min: float,
    mphi_global_max: float,
    lam1_global_min: float,
    lam1_global_max: float,
    n_mphi: int,
    n_lam1: int,
    label: str,
    expansion_index: int,
) -> SearchWindow:
    mphi_min, mphi_max = _centered_window(
        center=prev.m_phi,
        half_span=mphi_half_span,
        global_min=mphi_global_min,
        global_max=mphi_global_max,
    )
    lam1_min, lam1_max = _centered_window(
        center=prev.lam1,
        half_span=lam1_half_span,
        global_min=lam1_global_min,
        global_max=lam1_global_max,
    )
    return SearchWindow(
        mphi_min=mphi_min,
        mphi_max=mphi_max,
        n_mphi=int(n_mphi),
        lam1_min=lam1_min,
        lam1_max=lam1_max,
        n_lam1=int(n_lam1),
        label=label,
        expansion_index=int(expansion_index),
    )


def _global_window(
    *,
    mphi_min: float,
    mphi_max: float,
    n_mphi: int,
    lam1_min: float,
    lam1_max: float,
    n_lam1: int,
    label: str,
) -> SearchWindow:
    return SearchWindow(
        mphi_min=float(mphi_min),
        mphi_max=float(mphi_max),
        n_mphi=int(n_mphi),
        lam1_min=float(lam1_min),
        lam1_max=float(lam1_max),
        n_lam1=int(n_lam1),
        label=label,
        expansion_index=0,
    )


def _dedupe_preserve(values: list[float], *, rel_tol: float = 1e-15, abs_tol: float = 0.0) -> list[float]:
    out: list[float] = []
    for v in values:
        if not any(math.isclose(v, prev, rel_tol=rel_tol, abs_tol=abs_tol) for prev in out):
            out.append(v)
    return out


# -----------------------------------------------------------------------------
# Approximate level-curve model
# -----------------------------------------------------------------------------
def level_shape_tb(tb: float) -> float:
    if tb <= 0.0:
        raise ValueError("tanbeta must be > 0 for level-shape model")
    return float(((1.0 + tb * tb) ** 2) / (tb**5))


def lambda6_from_level_curve(tb: float, tb0: float, lambda6_0: float) -> float:
    return float(lambda6_0 * level_shape_tb(tb) / level_shape_tb(tb0))


def approx_level_invariant(tb: float, lambda6: float) -> float:
    # Constant along the approximate level curve if lambda6 scales with level_shape_tb.
    shape = level_shape_tb(tb)
    if shape == 0.0:
        raise ValueError("level_shape_tb returned zero")
    return float(lambda6 / shape)


# -----------------------------------------------------------------------------
# CSV parsing + candidate selection
# -----------------------------------------------------------------------------
def _flag_to_int(value: str | None) -> int:
    if value is None:
        return 0
    try:
        return 1 if float(value) == 1.0 else 0
    except Exception:
        return 0


def _safe_float_from_row(row: dict[str, str], key: str) -> float | None:
    value = row.get(key)
    if value is None:
        return None
    try:
        out = float(value)
    except Exception:
        return None
    if not math.isfinite(out):
        return None
    return out


def _load_candidates_from_csv(csv_path: Path) -> list[CandidatePoint]:
    candidates: list[CandidatePoint] = []
    if not csv_path.exists():
        return candidates

    with csv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader):
            m_phi = _safe_float_from_row(row, "m_phi")
            lam1 = _safe_float_from_row(row, "lam1")
            lam2 = _safe_float_from_row(row, "lam2")
            m12 = _safe_float_from_row(row, "m12")
            width_bb = _safe_float_from_row(row, "width_bb")
            width_gaga = _safe_float_from_row(row, "width_gaga")
            br_gaga = _safe_float_from_row(row, "br_gaga")
            total_width = _safe_float_from_row(row, "total_width")
            if None in (m_phi, lam1, lam2, m12, width_bb, width_gaga, br_gaga, total_width):
                continue
            candidates.append(
                CandidatePoint(
                    m_phi=float(m_phi),
                    lam1=float(lam1),
                    lam2=float(lam2),
                    m12=float(m12),
                    width_bb=float(width_bb),
                    width_gaga=float(width_gaga),
                    br_gaga=float(br_gaga),
                    total_width=float(total_width),
                    positivity_ok=_flag_to_int(row.get("positivity_ok")),
                    unitarity_ok=_flag_to_int(row.get("unitarity_ok")),
                    perturbativity_ok=_flag_to_int(row.get("perturbativity_ok")),
                    csv_path=str(csv_path),
                    row_index=int(idx),
                    raw={str(k): v for k, v in row.items()},
                )
            )
    return candidates


def _collect_csv_paths_from_run_dir(run_dir: Path) -> list[Path]:
    records = adaptive_artifacts.parse_task_summary_jsonl(run_dir / "task_summary.jsonl")
    paths: list[Path] = []
    for rec in records:
        raw_csv = rec.get("output_csv") or rec.get("previous_csv")
        if isinstance(raw_csv, str) and raw_csv.strip():
            p = Path(raw_csv)
            if p not in paths:
                paths.append(p)
    return paths


def _safe_log_distance(current: float, previous: float) -> float:
    if current <= 0.0 or previous <= 0.0:
        return float("inf")
    return abs(math.log(current) - math.log(previous))


def _normalized_delta(current: float, previous: float) -> float:
    scale = max(abs(previous), 1.0)
    return abs(current - previous) / scale


def _branch_distance(candidate: CandidatePoint, prev: BranchState, weights: DistanceWeights) -> float:
    score = 0.0
    score += weights.mphi * (_normalized_delta(candidate.m_phi, prev.m_phi) ** 2)
    score += weights.lam1 * (_normalized_delta(candidate.lam1, prev.lam1) ** 2)
    score += weights.lam2 * (_normalized_delta(candidate.lam2, prev.lam2) ** 2)
    score += weights.m12 * (_normalized_delta(candidate.m12, prev.m12) ** 2)

    log_bb = _safe_log_distance(candidate.width_bb, prev.width_bb)
    log_gaga = _safe_log_distance(candidate.width_gaga, prev.width_gaga)
    log_br = _safe_log_distance(candidate.br_gaga, prev.br_gaga)
    log_tot = _safe_log_distance(candidate.total_width, prev.total_width)

    score += weights.width_bb * ((log_bb**2) if math.isfinite(log_bb) else 1.0e6)
    score += weights.width_gaga * ((log_gaga**2) if math.isfinite(log_gaga) else 1.0e6)
    score += weights.br_gaga * ((log_br**2) if math.isfinite(log_br) else 1.0e6)
    score += weights.total_width * ((log_tot**2) if math.isfinite(log_tot) else 1.0e6)

    cand_r = candidate.r_exact
    prev_r = prev.r_exact
    if weights.r_exact > 0.0:
        if cand_r is None or prev_r is None:
            score += weights.r_exact * 1.0e6
        else:
            log_r = _safe_log_distance(cand_r, prev_r)
            score += weights.r_exact * ((log_r**2) if math.isfinite(log_r) else 1.0e6)

    return float(score)


def _seed_center_score(candidate: CandidatePoint, window: SearchWindow) -> float:
    center_mphi = 0.5 * (window.mphi_min + window.mphi_max)
    center_lam1 = 0.5 * (window.lam1_min + window.lam1_max)
    span_mphi = max(window.mphi_max - window.mphi_min, 1.0)
    span_lam1 = max(window.lam1_max - window.lam1_min, 1.0)
    dmphi = (candidate.m_phi - center_mphi) / span_mphi
    dlam1 = (candidate.lam1 - center_lam1) / span_lam1
    return float(dmphi * dmphi + dlam1 * dlam1)


def _select_seed_candidate(
    candidates: list[CandidatePoint],
    *,
    seed_policy: str,
    window: SearchWindow,
) -> tuple[CandidatePoint | None, float | None, str]:
    if not candidates:
        return None, None, "no-candidates"

    triple_ok = [c for c in candidates if c.triple_ok]
    pool = triple_ok if triple_ok else candidates
    if not pool:
        return None, None, "no-pool"

    if seed_policy == "center":
        ranked = sorted(pool, key=lambda c: (_seed_center_score(c, window), -c.br_gaga, -c.width_gaga, c.row_index))
        chosen = ranked[0]
        return chosen, _seed_center_score(chosen, window), "seed-center"

    if seed_policy == "max-br-gaga":
        ranked = sorted(pool, key=lambda c: (-c.br_gaga, -c.width_gaga, _seed_center_score(c, window), c.row_index))
        chosen = ranked[0]
        return chosen, float(-chosen.br_gaga), "seed-max-br-gaga"

    if seed_policy == "max-width-gaga":
        ranked = sorted(pool, key=lambda c: (-c.width_gaga, -c.br_gaga, _seed_center_score(c, window), c.row_index))
        chosen = ranked[0]
        return chosen, float(-chosen.width_gaga), "seed-max-width-gaga"

    raise ValueError(f"Unknown seed policy: {seed_policy}")


def _select_branch_candidate(
    candidates: list[CandidatePoint],
    *,
    prev: BranchState | None,
    weights: DistanceWeights,
    seed_policy: str,
    window: SearchWindow,
) -> tuple[CandidatePoint | None, float | None, str]:
    if prev is None:
        return _select_seed_candidate(candidates, seed_policy=seed_policy, window=window)

    triple_ok = [c for c in candidates if c.triple_ok]
    pool = triple_ok if triple_ok else candidates
    if not pool:
        return None, None, "no-pool"

    ranked = sorted(pool, key=lambda c: (_branch_distance(c, prev, weights), -c.br_gaga, -c.width_gaga, c.row_index))
    chosen = ranked[0]
    score = _branch_distance(chosen, prev, weights)
    return chosen, score, "branch-distance"


def _candidate_to_state(
    candidate: CandidatePoint,
    *,
    point: ExternalPoint,
    selection_reason: str,
) -> BranchState:
    return BranchState(
        lambda6=float(point.lambda6),
        tanbeta=float(point.tanbeta),
        m_phi=float(candidate.m_phi),
        lam1=float(candidate.lam1),
        lam2=float(candidate.lam2),
        m12=float(candidate.m12),
        width_bb=float(candidate.width_bb),
        width_gaga=float(candidate.width_gaga),
        br_gaga=float(candidate.br_gaga),
        total_width=float(candidate.total_width),
        r_exact=candidate.r_exact,
        approx_level_shape=point.approx_level_shape,
        approx_level_ratio=point.approx_level_ratio,
        csv_path=str(candidate.csv_path),
        row_index=int(candidate.row_index),
        selection_reason=selection_reason,
    )


# -----------------------------------------------------------------------------
# Orchestrator bridge
# -----------------------------------------------------------------------------
def _orchestrator_path() -> Path:
    return Path(__file__).with_name("orchestrate_scans.py")


def _build_run_name(base_prefix: str, point: ExternalPoint, step_label: str, attempt_index: int, window_label: str) -> str:
    l6_tag = orchestrate_scans.sanitize_for_path(orchestrate_scans.format_float_tag(point.lambda6, 8))
    tb_tag = orchestrate_scans.sanitize_for_path(orchestrate_scans.format_float_tag(point.tanbeta, 8))
    return f"{base_prefix}_l6_{l6_tag}_tb_{tb_tag}_{step_label}_try_{attempt_index:02d}_{window_label}"


def _build_orchestrator_command(
    *,
    exec_path: str,
    outdir: str,
    lake_name: str,
    campaign: str,
    run_name: str,
    threads: int | None,
    force: bool,
    dry_run: bool,
    fixed: FixedParams,
    point: ExternalPoint,
    window: SearchWindow,
    resume_scope: str,
    materialize: bool,
) -> list[str]:
    cmd: list[str] = [
        sys.executable,
        str(_orchestrator_path()),
        "--exec",
        exec_path,
        "--outdir",
        outdir,
        "--lake-name",
        lake_name,
        "--campaign",
        campaign,
        "--run-name",
        run_name,
        "--mA",
        str(fixed.mA),
        "--sin-ba",
        str(fixed.sin_ba),
        "--lambda6",
        str(point.lambda6),
        "--lambda7",
        str(fixed.lambda7),
        "--mphi-min",
        str(window.mphi_min),
        "--mphi-max",
        str(window.mphi_max),
        "--n-mphi",
        str(window.n_mphi),
        "--lam1-min",
        str(window.lam1_min),
        "--lam1-max",
        str(window.lam1_max),
        "--n-lam1",
        str(window.n_lam1),
        "--tanbeta",
        str(point.tanbeta),
        "--resume-scope",
        resume_scope,
    ]
    if threads is not None:
        cmd.extend(["--threads", str(int(threads))])
    if force:
        cmd.append("--force")
    if dry_run:
        cmd.append("--dry-run")
    if materialize:
        cmd.append("--materialize")
    return cmd


def _run_orchestrator(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


# -----------------------------------------------------------------------------
# Track builders
# -----------------------------------------------------------------------------
def _check_lambda6_policy(value: float, *, lo: float | None, hi: float | None, policy: str, where: str) -> tuple[bool, str | None]:
    if lo is not None and value < lo:
        msg = f"lambda6={value:g} is below lower bound {lo:g} at {where}"
        if policy == "error":
            raise ValueError(msg)
        if policy == "skip":
            return False, msg
        return True, msg
    if hi is not None and value > hi:
        msg = f"lambda6={value:g} is above upper bound {hi:g} at {where}"
        if policy == "error":
            raise ValueError(msg)
        if policy == "skip":
            return False, msg
        return True, msg
    return True, None


def _build_level_curve_points_from_seed(
    *,
    tb0: float,
    lambda6_0: float,
    tb_values: list[float] | None,
    tb_factor: float,
    n_up: int,
    n_down: int,
    lambda6_bound_min: float | None,
    lambda6_bound_max: float | None,
    lambda6_range_policy: str,
) -> tuple[ExternalPoint, list[ExternalPoint], list[ExternalPoint], dict[str, Any]]:
    if tb0 <= 0.0:
        raise ValueError("tb0 must be > 0")
    if tb_factor <= 0.0:
        raise ValueError("tb-factor must be > 0")
    if math.isclose(tb_factor, 1.0, rel_tol=1e-15, abs_tol=0.0) and (n_up > 0 or n_down > 0) and not tb_values:
        raise ValueError("tb-factor must differ from 1.0 when generating geometric steps")

    seed_point = ExternalPoint(
        tanbeta=float(tb0),
        lambda6=float(lambda6_0),
        source="level-curve-seed",
        approx_level_shape=level_shape_tb(tb0),
        approx_level_ratio=1.0,
        note="seed",
    )

    raw_up: list[float] = []
    raw_down: list[float] = []
    if tb_values:
        values = _dedupe_preserve([float(v) for v in tb_values])
        above = sorted([v for v in values if v > tb0])
        below = sorted([v for v in values if v < tb0], reverse=True)
        raw_up = above
        raw_down = below
    else:
        raw_up = [float(tb0 * (tb_factor**i)) for i in range(1, n_up + 1)]
        raw_down = [float(tb0 / (tb_factor**i)) for i in range(1, n_down + 1)]

    def make_points(seq: list[float], branch_label: str) -> list[ExternalPoint]:
        out: list[ExternalPoint] = []
        for idx, tb in enumerate(seq, start=1):
            l6 = lambda6_from_level_curve(tb, tb0, lambda6_0)
            keep, note = _check_lambda6_policy(
                l6,
                lo=lambda6_bound_min,
                hi=lambda6_bound_max,
                policy=lambda6_range_policy,
                where=f"{branch_label}[{idx}] tb={tb:g}",
            )
            if not keep:
                logger.warning("Skipping point: %s", note)
                continue
            out.append(
                ExternalPoint(
                    tanbeta=float(tb),
                    lambda6=float(l6),
                    source=f"level-curve-{branch_label}",
                    approx_level_shape=level_shape_tb(tb),
                    approx_level_ratio=(l6 / lambda6_0 if lambda6_0 != 0.0 else None),
                    note=note,
                )
            )
        return out

    metadata = {
        "tb0": tb0,
        "lambda6_0": lambda6_0,
        "tb_factor": tb_factor,
        "n_up": n_up,
        "n_down": n_down,
        "tb_values_input": tb_values,
        "approx_level_invariant": approx_level_invariant(tb0, lambda6_0),
        "lambda6_bound_min": lambda6_bound_min,
        "lambda6_bound_max": lambda6_bound_max,
        "lambda6_range_policy": lambda6_range_policy,
    }
    return seed_point, make_points(raw_up, "up"), make_points(raw_down, "down"), metadata


def _build_level_curve_tracks(args: argparse.Namespace) -> list[TrackDefinition]:
    pairs = _parse_seed_pairs(args.seed_pairs)
    if not pairs:
        if args.tb0 is None or args.lambda6_0 is None:
            raise ValueError("level-curve mode requires --tb0 and --lambda6-0, or --seed-pairs tb:l6,...")
        pairs = [(float(args.tb0), float(args.lambda6_0))]

    explicit_tb_values = _parse_csv_floats(args.tb_values) if args.tb_values else None
    tracks: list[TrackDefinition] = []
    for idx, (tb0, lambda6_0) in enumerate(pairs, start=1):
        seed_point, up_points, down_points, metadata = _build_level_curve_points_from_seed(
            tb0=float(tb0),
            lambda6_0=float(lambda6_0),
            tb_values=explicit_tb_values,
            tb_factor=float(args.tb_factor),
            n_up=int(args.n_up),
            n_down=int(args.n_down),
            lambda6_bound_min=args.lambda6_bound_min,
            lambda6_bound_max=args.lambda6_bound_max,
            lambda6_range_policy=str(args.lambda6_range_policy),
        )
        tracks.append(
            TrackDefinition(
                track_id=f"level_curve_{idx:03d}",
                label=f"level_curve_seed_tb_{orchestrate_scans.sanitize_for_path(orchestrate_scans.format_float_tag(tb0, 8))}_l6_{orchestrate_scans.sanitize_for_path(orchestrate_scans.format_float_tag(lambda6_0, 8))}",
                mode="level-curve",
                seed_point=seed_point,
                forward_points=up_points,
                backward_points=down_points,
                metadata=metadata,
            )
        )
    return tracks


def _build_manual_grid_tracks(args: argparse.Namespace) -> list[TrackDefinition]:
    lambda6_values = _resolve_value_grid(
        explicit_values=args.lambda6_values,
        min_value=args.lambda6_min,
        max_value=args.lambda6_max,
        count=args.n_lambda6,
        scale=args.lambda6_scale,
        field_name="lambda6",
    )
    tanbeta_values = _resolve_value_grid(
        explicit_values=args.tanbeta_values,
        min_value=args.tanbeta_min,
        max_value=args.tanbeta_max,
        count=args.n_tanbeta,
        scale=args.tanbeta_scale,
        field_name="tanbeta",
    )
    if not tanbeta_values:
        raise ValueError("manual-grid mode needs at least one tanbeta value")
    ordered_tb = sorted(tanbeta_values)
    tracks: list[TrackDefinition] = []
    for idx, l6 in enumerate(lambda6_values, start=1):
        seed_point = ExternalPoint(
            tanbeta=float(ordered_tb[0]),
            lambda6=float(l6),
            source="manual-grid-seed",
            approx_level_shape=(level_shape_tb(ordered_tb[0]) if ordered_tb[0] > 0 else None),
            approx_level_ratio=None,
            note=None,
        )
        forward = [
            ExternalPoint(
                tanbeta=float(tb),
                lambda6=float(l6),
                source="manual-grid-forward",
                approx_level_shape=(level_shape_tb(tb) if tb > 0 else None),
                approx_level_ratio=None,
                note=None,
            )
            for tb in ordered_tb[1:]
        ]
        tracks.append(
            TrackDefinition(
                track_id=f"manual_grid_{idx:03d}",
                label=f"manual_grid_l6_{orchestrate_scans.sanitize_for_path(orchestrate_scans.format_float_tag(l6, 8))}",
                mode="manual-grid",
                seed_point=seed_point,
                forward_points=forward,
                backward_points=[],
                metadata={"lambda6": l6, "tanbeta_values": ordered_tb},
            )
        )
    return tracks


def _build_manual_track_tracks(args: argparse.Namespace) -> list[TrackDefinition]:
    points_raw = _parse_manual_points(args.points)
    if not points_raw:
        raise ValueError("manual-track mode requires --points tb:l6,tb:l6,...")
    points = [
        ExternalPoint(
            tanbeta=float(tb),
            lambda6=float(l6),
            source="manual-track",
            approx_level_shape=(level_shape_tb(tb) if tb > 0 else None),
            approx_level_ratio=None,
            note=None,
        )
        for tb, l6 in points_raw
    ]
    return [
        TrackDefinition(
            track_id="manual_track_001",
            label="manual_track",
            mode="manual-track",
            seed_point=points[0],
            forward_points=points[1:],
            backward_points=[],
            metadata={"points": [asdict(p) for p in points]},
        )
    ]


def _build_track_definitions(args: argparse.Namespace) -> list[TrackDefinition]:
    if args.mode == "level-curve":
        return _build_level_curve_tracks(args)
    if args.mode == "manual-grid":
        return _build_manual_grid_tracks(args)
    if args.mode == "manual-track":
        return _build_manual_track_tracks(args)
    raise ValueError(f"Unknown mode: {args.mode}")


# -----------------------------------------------------------------------------
# Main continuation loop
# -----------------------------------------------------------------------------
def _make_windows_for_step(
    *,
    prev: BranchState | None,
    global_window: SearchWindow,
    local_mphi_half_span: float,
    local_lam1_half_span: float,
    local_n_mphi: int,
    local_n_lam1: int,
    max_expansions: int,
    expansion_factor: float,
    do_global_fallback: bool,
) -> list[SearchWindow]:
    windows: list[SearchWindow] = []
    if prev is None:
        windows.append(global_window)
        return windows

    for idx in range(max_expansions + 1):
        factor = expansion_factor**idx
        windows.append(
            _make_window_around_state(
                prev=prev,
                mphi_half_span=local_mphi_half_span * factor,
                lam1_half_span=local_lam1_half_span * factor,
                mphi_global_min=global_window.mphi_min,
                mphi_global_max=global_window.mphi_max,
                lam1_global_min=global_window.lam1_min,
                lam1_global_max=global_window.lam1_max,
                n_mphi=local_n_mphi,
                n_lam1=local_n_lam1,
                label=f"local_x{factor:.3g}",
                expansion_index=idx,
            )
        )

    if do_global_fallback:
        windows.append(
            SearchWindow(
                mphi_min=global_window.mphi_min,
                mphi_max=global_window.mphi_max,
                n_mphi=global_window.n_mphi,
                lam1_min=global_window.lam1_min,
                lam1_max=global_window.lam1_max,
                n_lam1=global_window.n_lam1,
                label="global_fallback",
                expansion_index=max_expansions + 1,
            )
        )
    return windows


def _evaluate_one_attempt(
    *,
    exec_path: str,
    outdir: str,
    lake_name: str,
    campaign: str,
    run_name_prefix: str,
    threads: int | None,
    force: bool,
    dry_run: bool,
    fixed: FixedParams,
    point: ExternalPoint,
    step_label: str,
    attempt_index: int,
    window: SearchWindow,
    prev: BranchState | None,
    weights: DistanceWeights,
    seed_policy: str,
    resume_scope: str,
    materialize: bool,
) -> AttemptResult:
    run_name = _build_run_name(
        base_prefix=run_name_prefix,
        point=point,
        step_label=step_label,
        attempt_index=attempt_index,
        window_label=window.label,
    )
    cmd = _build_orchestrator_command(
        exec_path=exec_path,
        outdir=outdir,
        lake_name=lake_name,
        campaign=campaign,
        run_name=run_name,
        threads=threads,
        force=force,
        dry_run=dry_run,
        fixed=fixed,
        point=point,
        window=window,
        resume_scope=resume_scope,
        materialize=materialize,
    )

    t0 = time.time()
    proc = _run_orchestrator(cmd)
    elapsed = time.time() - t0

    notes: list[str] = []
    stdout = proc.stdout or ""
    stderr = proc.stderr or ""
    if proc.returncode != 0:
        notes.append(f"orchestrator_returncode={proc.returncode}")
        if stderr.strip():
            notes.append(f"stderr_tail={stderr.strip()[-300:]}")
        if stdout.strip():
            notes.append(f"stdout_tail={stdout.strip()[-300:]}")
        return AttemptResult(
            status="orchestrator_failed",
            command=cmd,
            run_dir=str(adaptive_artifacts.parse_run_dir_from_orchestrator_output(stdout)) if stdout else None,
            csv_paths=[],
            n_candidates_total=0,
            n_candidates_triple_ok=0,
            selected=None,
            selection_score=None,
            window=window,
            elapsed_sec=float(elapsed),
            notes=notes,
        )

    run_dir = adaptive_artifacts.parse_run_dir_from_orchestrator_output(stdout)
    if run_dir is None:
        notes.append("run_dir_not_found")
        if stdout.strip():
            notes.append(f"stdout_tail={stdout.strip()[-300:]}")
        return AttemptResult(
            status="missing_run_dir",
            command=cmd,
            run_dir=None,
            csv_paths=[],
            n_candidates_total=0,
            n_candidates_triple_ok=0,
            selected=None,
            selection_score=None,
            window=window,
            elapsed_sec=float(elapsed),
            notes=notes,
        )

    if dry_run:
        return AttemptResult(
            status="dry_run",
            command=cmd,
            run_dir=str(run_dir),
            csv_paths=[],
            n_candidates_total=0,
            n_candidates_triple_ok=0,
            selected=None,
            selection_score=None,
            window=window,
            elapsed_sec=float(elapsed),
            notes=notes,
        )

    csv_paths = _collect_csv_paths_from_run_dir(run_dir)
    candidates: list[CandidatePoint] = []
    for csv_path in csv_paths:
        candidates.extend(_load_candidates_from_csv(csv_path))
    n_candidates_total = len(candidates)
    n_candidates_triple_ok = sum(1 for c in candidates if c.triple_ok)

    selected_candidate, score, selection_reason = _select_branch_candidate(
        candidates,
        prev=prev,
        weights=weights,
        seed_policy=seed_policy,
        window=window,
    )

    selected_state = None
    if selected_candidate is not None:
        selected_state = _candidate_to_state(
            selected_candidate,
            point=point,
            selection_reason=selection_reason,
        )

    status = "selected" if selected_state is not None else "no_candidate_selected"
    if n_candidates_total == 0:
        status = "no_candidates"
    elif n_candidates_triple_ok == 0:
        notes.append("no_triple_ok_candidates")

    return AttemptResult(
        status=status,
        command=cmd,
        run_dir=str(run_dir),
        csv_paths=[str(p) for p in csv_paths],
        n_candidates_total=n_candidates_total,
        n_candidates_triple_ok=n_candidates_triple_ok,
        selected=selected_state,
        selection_score=score,
        window=window,
        elapsed_sec=float(elapsed),
        notes=notes,
    )


def _run_point_with_windows(
    *,
    point: ExternalPoint,
    step_label: str,
    track_dir: Path,
    prev_state: BranchState | None,
    exec_path: str,
    outdir: str,
    lake_name: str,
    campaign: str,
    run_name_prefix: str,
    threads: int | None,
    force: bool,
    dry_run: bool,
    fixed: FixedParams,
    global_window: SearchWindow,
    local_mphi_half_span: float,
    local_lam1_half_span: float,
    local_n_mphi: int,
    local_n_lam1: int,
    max_expansions: int,
    expansion_factor: float,
    do_global_fallback: bool,
    weights: DistanceWeights,
    seed_policy: str,
    resume_scope: str,
    materialize: bool,
) -> tuple[AttemptResult | None, BranchState | None, list[dict[str, Any]]]:
    events_jsonl = track_dir / "events.jsonl"
    windows = _make_windows_for_step(
        prev=prev_state,
        global_window=global_window,
        local_mphi_half_span=local_mphi_half_span,
        local_lam1_half_span=local_lam1_half_span,
        local_n_mphi=local_n_mphi,
        local_n_lam1=local_n_lam1,
        max_expansions=max_expansions,
        expansion_factor=expansion_factor,
        do_global_fallback=do_global_fallback,
    )

    attempts_payload: list[dict[str, Any]] = []
    chosen_attempt: AttemptResult | None = None
    new_prev_state = prev_state

    for attempt_index, window in enumerate(windows, start=1):
        attempt_result = _evaluate_one_attempt(
            exec_path=exec_path,
            outdir=outdir,
            lake_name=lake_name,
            campaign=campaign,
            run_name_prefix=run_name_prefix,
            threads=threads,
            force=force,
            dry_run=dry_run,
            fixed=fixed,
            point=point,
            step_label=step_label,
            attempt_index=attempt_index,
            window=window,
            prev=prev_state,
            weights=weights,
            seed_policy=seed_policy,
            resume_scope=resume_scope,
            materialize=materialize,
        )
        attempt_dict = {
            "step_label": step_label,
            "point": asdict(point),
            "prev_state": (asdict(prev_state) if prev_state is not None else None),
            "attempt_index": attempt_index,
            "result": asdict(attempt_result),
        }
        attempts_payload.append(attempt_dict)
        _jsonl_append(events_jsonl, attempt_dict)

        if attempt_result.status == "selected":
            chosen_attempt = attempt_result
            new_prev_state = attempt_result.selected
            break
        if attempt_result.status == "dry_run":
            chosen_attempt = attempt_result
            new_prev_state = None
            break

    return chosen_attempt, new_prev_state, attempts_payload


def _run_branch_sequence(
    *,
    branch_name: str,
    points: list[ExternalPoint],
    initial_prev_state: BranchState | None,
    track_dir: Path,
    exec_path: str,
    outdir: str,
    lake_name: str,
    campaign: str,
    run_name_prefix: str,
    threads: int | None,
    force: bool,
    dry_run: bool,
    fixed: FixedParams,
    global_window: SearchWindow,
    local_mphi_half_span: float,
    local_lam1_half_span: float,
    local_n_mphi: int,
    local_n_lam1: int,
    max_expansions: int,
    expansion_factor: float,
    do_global_fallback: bool,
    weights: DistanceWeights,
    seed_policy: str,
    resume_scope: str,
    materialize: bool,
) -> dict[str, Any]:
    prev_state = initial_prev_state
    steps_payload: list[dict[str, Any]] = []
    status = "completed"

    for idx, point in enumerate(points, start=1):
        step_label = f"{branch_name}_{idx:04d}"
        logger.info("Track branch=%s | step=%s | tb=%s | l6=%s", branch_name, idx, point.tanbeta, point.lambda6)
        chosen_attempt, prev_state, attempts_payload = _run_point_with_windows(
            point=point,
            step_label=step_label,
            track_dir=track_dir,
            prev_state=prev_state,
            exec_path=exec_path,
            outdir=outdir,
            lake_name=lake_name,
            campaign=campaign,
            run_name_prefix=run_name_prefix,
            threads=threads,
            force=force,
            dry_run=dry_run,
            fixed=fixed,
            global_window=global_window,
            local_mphi_half_span=local_mphi_half_span,
            local_lam1_half_span=local_lam1_half_span,
            local_n_mphi=local_n_mphi,
            local_n_lam1=local_n_lam1,
            max_expansions=max_expansions,
            expansion_factor=expansion_factor,
            do_global_fallback=do_global_fallback,
            weights=weights,
            seed_policy=seed_policy,
            resume_scope=resume_scope,
            materialize=materialize,
        )
        step_payload = {
            "step_label": step_label,
            "point": asdict(point),
            "attempts": attempts_payload,
            "selected": (asdict(chosen_attempt.selected) if (chosen_attempt and chosen_attempt.selected) else None),
            "status": (chosen_attempt.status if chosen_attempt is not None else "no_attempts"),
        }
        steps_payload.append(step_payload)

        if chosen_attempt is None:
            status = "no_attempts"
            break
        if chosen_attempt.status not in {"selected", "dry_run"}:
            status = chosen_attempt.status
            break
        if dry_run:
            prev_state = None

    return {
        "branch_name": branch_name,
        "status": status,
        "n_requested_steps": len(points),
        "n_recorded_steps": len(steps_payload),
        "steps": steps_payload,
    }


def _run_track_definition(
    *,
    track: TrackDefinition,
    exec_path: str,
    outdir: str,
    lake_name: str,
    campaign: str,
    checkpoint_root: Path,
    run_name_prefix: str,
    threads: int | None,
    force: bool,
    dry_run: bool,
    fixed: FixedParams,
    global_window: SearchWindow,
    local_mphi_half_span: float,
    local_lam1_half_span: float,
    local_n_mphi: int,
    local_n_lam1: int,
    max_expansions: int,
    expansion_factor: float,
    do_global_fallback: bool,
    weights: DistanceWeights,
    seed_policy: str,
    resume_scope: str,
    materialize: bool,
) -> dict[str, Any]:
    track_dir = checkpoint_root / track.track_id
    track_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Starting track %s (%s)", track.track_id, track.label)

    seed_attempt, seed_state, seed_attempts = _run_point_with_windows(
        point=track.seed_point,
        step_label="seed",
        track_dir=track_dir,
        prev_state=None,
        exec_path=exec_path,
        outdir=outdir,
        lake_name=lake_name,
        campaign=campaign,
        run_name_prefix=run_name_prefix,
        threads=threads,
        force=force,
        dry_run=dry_run,
        fixed=fixed,
        global_window=global_window,
        local_mphi_half_span=local_mphi_half_span,
        local_lam1_half_span=local_lam1_half_span,
        local_n_mphi=local_n_mphi,
        local_n_lam1=local_n_lam1,
        max_expansions=max_expansions,
        expansion_factor=expansion_factor,
        do_global_fallback=do_global_fallback,
        weights=weights,
        seed_policy=seed_policy,
        resume_scope=resume_scope,
        materialize=materialize,
    )

    seed_payload = {
        "point": asdict(track.seed_point),
        "attempts": seed_attempts,
        "selected": (asdict(seed_attempt.selected) if (seed_attempt and seed_attempt.selected) else None),
        "status": (seed_attempt.status if seed_attempt is not None else "no_attempts"),
    }

    branches: list[dict[str, Any]] = []
    if seed_attempt is not None and seed_attempt.status in {"selected", "dry_run"}:
        if track.forward_points:
            branches.append(
                _run_branch_sequence(
                    branch_name="forward",
                    points=track.forward_points,
                    initial_prev_state=seed_state,
                    track_dir=track_dir,
                    exec_path=exec_path,
                    outdir=outdir,
                    lake_name=lake_name,
                    campaign=campaign,
                    run_name_prefix=run_name_prefix,
                    threads=threads,
                    force=force,
                    dry_run=dry_run,
                    fixed=fixed,
                    global_window=global_window,
                    local_mphi_half_span=local_mphi_half_span,
                    local_lam1_half_span=local_lam1_half_span,
                    local_n_mphi=local_n_mphi,
                    local_n_lam1=local_n_lam1,
                    max_expansions=max_expansions,
                    expansion_factor=expansion_factor,
                    do_global_fallback=do_global_fallback,
                    weights=weights,
                    seed_policy=seed_policy,
                    resume_scope=resume_scope,
                    materialize=materialize,
                )
            )
        if track.backward_points:
            branches.append(
                _run_branch_sequence(
                    branch_name="backward",
                    points=track.backward_points,
                    initial_prev_state=seed_state,
                    track_dir=track_dir,
                    exec_path=exec_path,
                    outdir=outdir,
                    lake_name=lake_name,
                    campaign=campaign,
                    run_name_prefix=run_name_prefix,
                    threads=threads,
                    force=force,
                    dry_run=dry_run,
                    fixed=fixed,
                    global_window=global_window,
                    local_mphi_half_span=local_mphi_half_span,
                    local_lam1_half_span=local_lam1_half_span,
                    local_n_mphi=local_n_mphi,
                    local_n_lam1=local_n_lam1,
                    max_expansions=max_expansions,
                    expansion_factor=expansion_factor,
                    do_global_fallback=do_global_fallback,
                    weights=weights,
                    seed_policy=seed_policy,
                    resume_scope=resume_scope,
                    materialize=materialize,
                )
            )

    track_summary = {
        "track_id": track.track_id,
        "label": track.label,
        "mode": track.mode,
        "metadata": track.metadata,
        "seed": seed_payload,
        "branches": branches,
    }
    _json_dump(track_dir / "track_summary.json", track_summary)
    return track_summary


# -----------------------------------------------------------------------------
# Pre-mortem
# -----------------------------------------------------------------------------
def _premortem(args: argparse.Namespace, tracks: list[TrackDefinition]) -> dict[str, Any]:
    risks: list[dict[str, Any]] = []

    if args.mode == "level-curve" and args.w_r_exact <= 0.0:
        risks.append(
            {
                "risk": "target_drift_not_penalized",
                "severity": "medium",
                "detail": "The search preserves branch continuity, but exact R drift is not part of the continuity score because w-r-exact <= 0.",
                "mitigation": "Use --w-r-exact 0.05 or 0.10 if exact ratio stability should influence branch selection.",
            }
        )

    if args.mode == "level-curve" and args.tb_factor > 3.0:
        risks.append(
            {
                "risk": "steps_too_aggressive",
                "severity": "medium",
                "detail": "Large geometric jumps in tanbeta can push the branch outside the local continuation window and force global fallback.",
                "mitigation": "Prefer tb-factor around 1.5-2.0 for first tests.",
            }
        )

    if args.max_expansions == 0 and args.no_global_fallback:
        risks.append(
            {
                "risk": "fragile_local_search",
                "severity": "high",
                "detail": "Without expansions and without global fallback, the branch can be lost after a single local miss.",
                "mitigation": "Keep at least one expansion or allow global fallback.",
            }
        )

    if args.resume_scope == "fixed" and not args.force:
        risks.append(
            {
                "risk": "cross_run_reuse",
                "severity": "low",
                "detail": "The orchestrator may reuse previous CSVs under the same fixed hyperparameters if the effective grid matches.",
                "mitigation": "Use --force when you want a fresh recomputation for every step.",
            }
        )

    for track in tracks:
        all_points = [track.seed_point] + list(track.forward_points) + list(track.backward_points)
        if len(all_points) <= 1:
            risks.append(
                {
                    "risk": "degenerate_track",
                    "severity": "medium",
                    "detail": f"Track {track.track_id} contains only the seed point, so it cannot test a level relation.",
                    "mitigation": "Increase n-up / n-down or provide more tb-values.",
                }
            )
        if track.mode == "level-curve":
            vals = [p.lambda6 for p in all_points]
            if any(v == 0.0 for v in vals):
                risks.append(
                    {
                        "risk": "zero_lambda6_on_track",
                        "severity": "medium",
                        "detail": f"Track {track.track_id} contains lambda6=0. This can kill the approximate guide and some comparisons in log space.",
                        "mitigation": "Use non-zero seed lambda6 or filter out pathological steps.",
                    }
                )

    return {
        "mode": args.mode,
        "n_tracks": len(tracks),
        "risks": risks,
        "risk_count": len(risks),
    }


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "2HDM branch-continuation explorer focused on (tanbeta, lambda6). "
            "Default mode is level-curve: generate external points from the approximate R relation, "
            "then use the current C++ motor via orchestrate_scans.py to preserve the same internal branch."
        )
    )

    # High-level mode
    p.add_argument("--mode", type=str, default="level-curve", choices=["level-curve", "manual-grid", "manual-track"], help="External-point generation mode.")

    # Runtime / paths
    p.add_argument("--exec", type=str, default=_default_exec_path(), help="Path to compiled C++ executable.")
    p.add_argument("--outdir", type=str, default=orchestrate_scans.DEFAULT_CERNBOX_OUTDIR, help="Root output directory.")
    p.add_argument("--lake-name", type=str, default=orchestrate_scans.DEFAULT_LAKE_DIRNAME, help="Lake directory name.")
    p.add_argument("--campaign", type=str, default="level_curve_branch", help="Campaign label.")
    p.add_argument("--checkpoint-root", type=str, default=None, help="Directory for explorer artifacts. Default: <outdir>/checkpoints/<campaign>")
    p.add_argument("--run-name-prefix", type=str, default="branchcurve", help="Prefix for orchestrator run_name.")
    p.add_argument("--threads", type=int, default=None, help="OMP_NUM_THREADS to forward to orchestrate_scans.py.")
    p.add_argument("--force", action="store_true", help="Force orchestrator to overwrite existing CSVs.")
    p.add_argument("--dry-run", action="store_true", help="Plan and invoke orchestrator in dry-run mode only.")
    p.add_argument("--resume-scope", type=str, default="fixed", choices=["run", "fixed"], help="Forwarded to orchestrator.")
    p.add_argument("--materialize", action="store_true", help="Forwarded to orchestrator.")

    # Fixed physical params
    p.add_argument("--mA", type=float, default=300.0, help="Fixed mA. Default tuned for common use.")
    p.add_argument("--sin-ba", type=float, default=1.0, help="Fixed sin(b-a). Default: 1.0")
    p.add_argument("--lambda7", type=float, default=0.0, help="Fixed lambda7. Default: 0.0")

    # Level-curve mode inputs
    p.add_argument("--tb0", type=float, default=None, help="Seed tanbeta for level-curve mode.")
    p.add_argument("--lambda6-0", type=float, default=None, help="Seed lambda6 at tb0 for level-curve mode.")
    p.add_argument("--seed-pairs", type=str, default=None, help="Multiple level-curve seeds as tb:l6,tb:l6,...")
    p.add_argument("--tb-factor", type=float, default=2.0, help="Geometric step factor in level-curve mode.")
    p.add_argument("--n-up", type=int, default=1, help="Number of upward tanbeta steps from the seed.")
    p.add_argument("--n-down", type=int, default=1, help="Number of downward tanbeta steps from the seed.")
    p.add_argument("--tb-values", type=str, default=None, help="Explicit tanbeta values for level-curve mode. Seed tb0 is still the anchor; values are split above/below tb0.")
    p.add_argument("--lambda6-bound-min", type=float, default=None, help="Optional lower bound for derived lambda6 in level-curve mode.")
    p.add_argument("--lambda6-bound-max", type=float, default=None, help="Optional upper bound for derived lambda6 in level-curve mode.")
    p.add_argument("--lambda6-range-policy", type=str, default="warn", choices=["warn", "skip", "error"], help="How to handle derived lambda6 values outside the allowed range.")

    # Manual-grid mode inputs
    p.add_argument("--lambda6-values", type=str, default=None, help="Explicit lambda6 list, comma-separated.")
    p.add_argument("--lambda6-min", type=float, default=None, help="lambda6 min when not using --lambda6-values.")
    p.add_argument("--lambda6-max", type=float, default=None, help="lambda6 max when not using --lambda6-values.")
    p.add_argument("--n-lambda6", type=int, default=None, help="Number of lambda6 points.")
    p.add_argument("--lambda6-scale", type=str, default="linear", choices=["linear", "log"], help="Grid scale for lambda6 range mode.")

    p.add_argument("--tanbeta-values", type=str, default=None, help="Explicit tanbeta list, comma-separated.")
    p.add_argument("--tanbeta-min", type=float, default=None, help="tanbeta min when not using --tanbeta-values.")
    p.add_argument("--tanbeta-max", type=float, default=None, help="tanbeta max when not using --tanbeta-values.")
    p.add_argument("--n-tanbeta", type=int, default=None, help="Number of tanbeta points.")
    p.add_argument("--tanbeta-scale", type=str, default="log", choices=["linear", "log"], help="Grid scale for tanbeta range mode.")

    # Manual-track mode
    p.add_argument("--points", type=str, default=None, help="Manual sequence as tb:l6,tb:l6,... used in manual-track mode.")

    # Global search window
    p.add_argument("--mphi-min", type=float, default=130.0, help="Global m_phi min.")
    p.add_argument("--mphi-max", type=float, default=290.0, help="Global m_phi max.")
    p.add_argument("--n-mphi-global", type=int, default=15, help="Global N_mphi for first search / fallback.")
    p.add_argument("--lam1-min", type=float, default=0.0, help="Global lambda1 min.")
    p.add_argument("--lam1-max", type=float, default=12.3, help="Global lambda1 max.")
    p.add_argument("--n-lam1-global", type=int, default=80, help="Global N_lam1 for first search / fallback.")

    # Local branch-continuation window
    p.add_argument("--local-mphi-half-span", type=float, default=8.0, help="Half-span around previous m_phi.")
    p.add_argument("--local-lam1-half-span", type=float, default=0.5, help="Half-span around previous lam1.")
    p.add_argument("--n-mphi-local", type=int, default=7, help="N_mphi for local windows.")
    p.add_argument("--n-lam1-local", type=int, default=25, help="N_lam1 for local windows.")
    p.add_argument("--max-expansions", type=int, default=2, help="Number of local window expansions before fallback.")
    p.add_argument("--expansion-factor", type=float, default=2.0, help="Half-span multiplier between expansions.")
    p.add_argument("--no-global-fallback", action="store_true", help="Disable final global fallback window.")

    # Selection controls
    p.add_argument(
        "--seed-policy",
        type=str,
        default="center",
        choices=["center", "max-br-gaga", "max-width-gaga"],
        help="How to choose the first branch seed when there is no previous state.",
    )
    p.add_argument("--w-mphi", type=float, default=1.0, help="Distance weight for m_phi continuity.")
    p.add_argument("--w-lam1", type=float, default=1.0, help="Distance weight for lam1 continuity.")
    p.add_argument("--w-lam2", type=float, default=0.75, help="Distance weight for lam2 continuity.")
    p.add_argument("--w-m12", type=float, default=0.25, help="Distance weight for m12 continuity.")
    p.add_argument("--w-width-bb", type=float, default=0.35, help="Distance weight for width_bb continuity.")
    p.add_argument("--w-width-gaga", type=float, default=0.35, help="Distance weight for width_gaga continuity.")
    p.add_argument("--w-br-gaga", type=float, default=0.20, help="Distance weight for BR_gaga continuity.")
    p.add_argument("--w-total-width", type=float, default=0.10, help="Distance weight for total_width continuity.")
    p.add_argument("--w-r-exact", type=float, default=0.0, help="Distance weight for exact ratio R = width_gaga / width_bb continuity.")

    return p


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    exec_path = str(Path(args.exec).expanduser())
    if not Path(exec_path).exists():
        raise SystemExit(f"Executable not found: {exec_path}")

    if args.expansion_factor < 1.0:
        raise SystemExit("--expansion-factor must be >= 1.0")
    if args.local_mphi_half_span <= 0.0:
        raise SystemExit("--local-mphi-half-span must be > 0")
    if args.local_lam1_half_span <= 0.0:
        raise SystemExit("--local-lam1-half-span must be > 0")
    if args.n_mphi_local < 1 or args.n_lam1_local < 1:
        raise SystemExit("local grid sizes must be >= 1")
    if args.n_mphi_global < 1 or args.n_lam1_global < 1:
        raise SystemExit("global grid sizes must be >= 1")
    if args.mphi_max <= args.mphi_min:
        raise SystemExit("mphi range must satisfy max > min")
    if args.lam1_max <= args.lam1_min:
        raise SystemExit("lam1 range must satisfy max > min")

    if args.checkpoint_root is None:
        checkpoint_root = Path(args.outdir).expanduser().resolve() / "checkpoints" / orchestrate_scans.sanitize_for_path(args.campaign)
    else:
        checkpoint_root = Path(args.checkpoint_root).expanduser().resolve()
    checkpoint_root.mkdir(parents=True, exist_ok=True)

    fixed = FixedParams(mA=float(args.mA), sin_ba=float(args.sin_ba), lambda7=float(args.lambda7))
    global_window = _global_window(
        mphi_min=float(args.mphi_min),
        mphi_max=float(args.mphi_max),
        n_mphi=int(args.n_mphi_global),
        lam1_min=float(args.lam1_min),
        lam1_max=float(args.lam1_max),
        n_lam1=int(args.n_lam1_global),
        label="global_seed",
    )
    weights = DistanceWeights(
        mphi=float(args.w_mphi),
        lam1=float(args.w_lam1),
        lam2=float(args.w_lam2),
        m12=float(args.w_m12),
        width_bb=float(args.w_width_bb),
        width_gaga=float(args.w_width_gaga),
        br_gaga=float(args.w_br_gaga),
        total_width=float(args.w_total_width),
        r_exact=float(args.w_r_exact),
    )

    tracks = _build_track_definitions(args)
    premortem = _premortem(args, tracks)

    manifest = {
        "created_unix": time.time(),
        "script": str(Path(__file__).resolve()),
        "mode": args.mode,
        "exec": exec_path,
        "outdir": args.outdir,
        "lake_name": args.lake_name,
        "campaign": args.campaign,
        "checkpoint_root": str(checkpoint_root),
        "dry_run": args.dry_run,
        "force": args.force,
        "resume_scope": args.resume_scope,
        "materialize": args.materialize,
        "fixed": asdict(fixed),
        "global_window": asdict(global_window),
        "local_window": {
            "local_mphi_half_span": args.local_mphi_half_span,
            "local_lam1_half_span": args.local_lam1_half_span,
            "n_mphi_local": args.n_mphi_local,
            "n_lam1_local": args.n_lam1_local,
            "max_expansions": args.max_expansions,
            "expansion_factor": args.expansion_factor,
            "global_fallback": not args.no_global_fallback,
        },
        "weights": asdict(weights),
        "seed_policy": args.seed_policy,
        "tracks": [
            {
                "track_id": t.track_id,
                "label": t.label,
                "mode": t.mode,
                "seed_point": asdict(t.seed_point),
                "forward_points": [asdict(p) for p in t.forward_points],
                "backward_points": [asdict(p) for p in t.backward_points],
                "metadata": t.metadata,
            }
            for t in tracks
        ],
        "premortem": premortem,
    }
    _json_dump(checkpoint_root / "manifest.json", manifest)
    _json_dump(checkpoint_root / "premortem.json", premortem)

    logger.info("mode = %s", args.mode)
    logger.info("checkpoint_root = %s", checkpoint_root)
    logger.info("n_tracks = %s", len(tracks))
    if premortem["risk_count"]:
        logger.warning("Pre-mortem flagged %d risk(s). See %s", premortem["risk_count"], checkpoint_root / "premortem.json")

    all_tracks: list[dict[str, Any]] = []
    for track in tracks:
        track_summary = _run_track_definition(
            track=track,
            exec_path=exec_path,
            outdir=args.outdir,
            lake_name=args.lake_name,
            campaign=args.campaign,
            checkpoint_root=checkpoint_root,
            run_name_prefix=args.run_name_prefix,
            threads=args.threads,
            force=args.force,
            dry_run=args.dry_run,
            fixed=fixed,
            global_window=global_window,
            local_mphi_half_span=float(args.local_mphi_half_span),
            local_lam1_half_span=float(args.local_lam1_half_span),
            local_n_mphi=int(args.n_mphi_local),
            local_n_lam1=int(args.n_lam1_local),
            max_expansions=int(args.max_expansions),
            expansion_factor=float(args.expansion_factor),
            do_global_fallback=not bool(args.no_global_fallback),
            weights=weights,
            seed_policy=str(args.seed_policy),
            resume_scope=str(args.resume_scope),
            materialize=bool(args.materialize),
        )
        all_tracks.append(track_summary)

    summary = {
        "campaign": args.campaign,
        "mode": args.mode,
        "n_tracks": len(all_tracks),
        "tracks": all_tracks,
        "premortem": premortem,
    }
    _json_dump(checkpoint_root / "summary.json", summary)

    logger.info("Done. Summary written to %s", checkpoint_root / "summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
