#!/usr/bin/env python3
"""Adaptive M2 continuation in mH2 for fixed heavy-state slices mA=mHp.

Traces a continuous, understandable sequence of theory-valid 2HDM points as
mH2 rises from just above the SM-like Higgs mass toward mA, on two heavy
slices (mA=mHp=300, 350 GeV). Calls the existing canonical
DihiggsPointV2Evaluator (schema dihiggs.point.v2) as a subprocess for every
physics evaluation. Does not modify, wrap, or reimplement 2HDMC physics.

Method (adaptive continuation, not a scan):
  0. Seed calibration (once per slice, at the lowest target mass): the
     starting sin(beta-alpha)=1.0, lambda6=1e-10 combination is held, but
     tan_beta must be searched downward from the requested 300000 anchor
     because that value is theory-invalid everywhere on these two slices
     (confirmed by direct wide scan -- see REPORT.md). The search finds the
     largest tan_beta admitting a "workable" band (>= MIN_BAND_POINTS valid
     M2 grid points) via a coarse geometric sweep then a log-space bisection
     refinement, and that tan_beta is then held fixed going forward exactly
     like lambda6/lambda7/sin_ba/Yukawa type.
  1. For each new target mH2: try the previous accepted point's exact
     (M2, tan_beta, lambda6) unchanged.
  2. If that fails: local M2 window search around a linear/constant
     prediction from the last 1-2 accepted points, expanding up to 3 times.
  3. If a window has valid points: refine once around the nearest-to-
     prediction candidate for a cleaner M2 (not just a coarse grid node).
  4. If M2 search fails entirely: try tan_beta x{0.5, 2.0} (relative to the
     CURRENT, possibly-already-relaxed tan_beta), each with one local M2
     window.
  5. If that fails: try lambda6 in {1e-9, 1e-11} (relative to the current
     lambda6), each with one local M2 window.
  6. If everything fails, record the mass as FAILED (with the dominant
     rejection_stage) and continue to the next mass using the last accepted
     point as the seed.

Every use of rung 4/5 is stamped in the output row (relaxation_used,
moved_coordinate) since the goal asks specifically when this extra freedom
becomes necessary.
"""
import argparse
import csv
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
BINARY_A = REPO_ROOT / "dihiggs" / "app" / "DihiggsPointV2Evaluator"
BINARY_B = REPO_ROOT / "dihiggs" / "app" / "Lambda1EvaluatorV2"

# --- fixed physics inputs (per the approved plan) --------------------------
MH = 125.20
SIN_BA = 1.0
LAMBDA6_SEED = 1e-10
LAMBDA7 = 0.0
YUKAWA_TYPE = 1
TAN_BETA_REQUESTED = 300000.0
TAN_BETA_FLOOR = 1.01
MIN_BAND_POINTS = 5          # "workable band", not a knife-edge single point
MZ = 91.15349
MW = 80.36951
HBAR_C_GEV_MM = 1.973269804e-13

# --- continuation ladder tuning ---------------------------------------------
LOCAL_N = 25
REFINE_N = 41
MAX_EXPANSIONS = 3
EXPANSION_FACTOR = 2.0
TAN_BETA_RELAX_FACTORS = (0.5, 2.0)
LAMBDA6_RELAX_VALUES = (1e-9, 1e-11)

VALID_TOKEN = "1.00000000000000000e+00"


def git_commit():
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True,
        capture_output=True, check=True,
    ).stdout.strip()


def git_dirty():
    out = subprocess.run(
        ["git", "status", "--porcelain"], cwd=REPO_ROOT, text=True,
        capture_output=True, check=True,
    ).stdout
    return "yes" if out.strip() else "no"


_ENV_CACHE = None


def eval_env():
    global _ENV_CACHE
    if _ENV_CACHE is None:
        env = dict(os.environ)
        env["DIHIGGS_GIT_COMMIT"] = git_commit()
        env["DIHIGGS_GIT_DIRTY"] = git_dirty()
        _ENV_CACHE = env
    return _ENV_CACHE


def run_grid(mH, M2_lo, M2_hi, n_M2, tan_beta, lambda6, mA, mHp,
             campaign_id, run_id, out_path, mh=MH, sin_ba=SIN_BA,
             lambda7=LAMBDA7, yukawa_type=YUKAWA_TYPE):
    if n_M2 == 1 and M2_lo != M2_hi:
        M2_hi = M2_lo
    cmd = [
        str(BINARY_A),
        "--campaign-id", campaign_id, "--run-id", run_id,
        "--mh", repr(mh),
        "--mH-min", repr(mH), "--mH-max", repr(mH), "--n-mH", "1",
        "--mA", repr(mA), "--mHp", repr(mHp),
        "--yukawa-type", str(yukawa_type),
        "--sin-ba", repr(sin_ba), "--tan-beta", repr(tan_beta),
        "--M2-min", repr(M2_lo), "--M2-max", repr(M2_hi), "--n-M2", str(n_M2),
        "--lambda6", repr(lambda6), "--lambda7", repr(lambda7),
        "--output", str(out_path),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True, env=eval_env())
    with open(out_path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def is_valid(row):
    return row["theory_ok_v1"] == VALID_TOKEN


def dominant_rejection(rows):
    stages = [r["rejection_stage"] for r in rows if r["rejection_stage"] != "none"]
    if not stages:
        return "none"
    counts = {}
    for s in stages:
        counts[s] = counts.get(s, 0) + 1
    return max(counts.items(), key=lambda kv: kv[1])[0]


def pick_nearest(valid_rows, target):
    if target is None:
        target = 0.0
    return min(valid_rows, key=lambda r: abs(float(r["M2_input_GeV2"]) - target))


# --- seed calibration --------------------------------------------------------

def bootstrap_tan_beta(mA, mHp, mH2_seed, campaign_id, workdir, logf):
    """Find the largest tan_beta (<=TAN_BETA_REQUESTED) admitting a workable
    M2 band at mH2_seed, holding sin_ba/lambda6/lambda7 at their seed values.
    Coarse geometric sweep (/3 per step) then log-space bisection refine."""
    scratch = workdir / "_scratch_bootstrap.csv"
    M2_half = 4.0 * mA * mA

    def probe(tb):
        rows = run_grid(mH2_seed, -M2_half, M2_half, 401, tb, LAMBDA6_SEED,
                         mA, mHp, campaign_id, "bootstrap_tb", scratch)
        valid = [r for r in rows if is_valid(r)]
        logf.write(f"[seed-cal mA={mA:g}] tan_beta={tb:.6g} valid_M2_count={len(valid)}\n")
        logf.flush()
        return valid

    candidates = []
    tb = TAN_BETA_REQUESTED
    while tb > TAN_BETA_FLOOR:
        candidates.append(tb)
        tb /= 3.0
    candidates.append(TAN_BETA_FLOOR)

    history = []
    passing_idx = None
    for i, tb in enumerate(candidates):
        valid = probe(tb)
        history.append((tb, len(valid)))
        if len(valid) >= MIN_BAND_POINTS:
            passing_idx = i
            break

    if passing_idx is None:
        # nothing passed even at the floor -- fall back to whatever had the
        # most valid points, however small.
        history.sort(key=lambda t: -t[1])
        best_tb, best_n = history[0]
        valid = probe(best_tb)
        logf.write(f"[seed-cal mA={mA:g}] WARNING: no candidate reached "
                    f"MIN_BAND_POINTS; using best-effort tan_beta={best_tb:.6g} "
                    f"(n_valid={best_n})\n")
        return best_tb, valid

    tb_pass = candidates[passing_idx]
    if passing_idx == 0:
        return tb_pass, probe(tb_pass)
    tb_fail = candidates[passing_idx - 1]

    for _ in range(6):
        mid = math.sqrt(tb_fail * tb_pass)
        valid_mid = probe(mid)
        if len(valid_mid) >= MIN_BAND_POINTS:
            tb_pass = mid
        else:
            tb_fail = mid

    final_valid = probe(tb_pass)
    logf.write(f"[seed-cal mA={mA:g}] SELECTED tan_beta={tb_pass!r} "
               f"(bracket_fail={tb_fail!r}, n_valid={len(final_valid)})\n")
    logf.flush()
    return tb_pass, final_valid


# --- continuation state -------------------------------------------------------

class TrajectoryPoint:
    __slots__ = ("mH2", "M2", "tan_beta", "lambda6", "row",
                 "n_attempt_rounds", "n_M2_probes", "moved_coordinate",
                 "relaxation_used", "M2_prediction_error", "failed")

    def __init__(self, mH2, M2, tan_beta, lambda6, row, n_attempt_rounds,
                 n_M2_probes, moved_coordinate, relaxation_used,
                 M2_prediction_error, failed=False):
        self.mH2 = mH2
        self.M2 = M2
        self.tan_beta = tan_beta
        self.lambda6 = lambda6
        self.row = row
        self.n_attempt_rounds = n_attempt_rounds
        self.n_M2_probes = n_M2_probes
        self.moved_coordinate = moved_coordinate
        self.relaxation_used = relaxation_used
        self.M2_prediction_error = M2_prediction_error
        self.failed = failed


def predict_M2(history):
    if not history:
        return None
    if len(history) == 1:
        return history[-1].M2
    a, b = history[-2], history[-1]
    if b.mH2 == a.mH2:
        return b.M2
    slope = (b.M2 - a.M2) / (b.mH2 - a.mH2)
    return b.M2  # caller adds slope*(target-b.mH2); kept for clarity below


def predicted_M2_for(history, target_mH2):
    if not history:
        return None
    if len(history) == 1:
        return history[-1].M2
    a, b = history[-2], history[-1]
    if b.mH2 == a.mH2:
        return b.M2
    slope = (b.M2 - a.M2) / (b.mH2 - a.mH2)
    return b.M2 + slope * (target_mH2 - b.mH2)


def try_mass(target_mH2, mA, mHp, history, baseline_tan_beta, baseline_lambda6,
             campaign_id, workdir, logf, attempt_log_prefix):
    scratch = workdir / "_scratch_ladder.csv"
    n_probes_total = 0
    n_attempt_rounds = 0
    last_rows_for_reason = []

    cur_tan_beta = history[-1].tan_beta if history else baseline_tan_beta
    cur_lambda6 = history[-1].lambda6 if history else baseline_lambda6
    predicted = predicted_M2_for(history, target_mH2)

    def log_attempt(tag, tan_beta, lambda6, best_row, valid_found):
        if best_row is not None:
            status = ("theory=OK" if valid_found else
                       f"theory=FAIL:{best_row['rejection_stage']}")
            constr = "OK" if best_row["construction_ok"] not in ("0", "") else "FAIL"
            logf.write(
                f"mH2={target_mH2:g} attempt={n_attempt_rounds} tag={tag} "
                f"M2={float(best_row['M2_input_GeV2']):.6g} "
                f"tanbeta={tan_beta:.6g} lambda6={lambda6:.3g} "
                f"construction={constr} {status}\n"
            )
        else:
            logf.write(
                f"mH2={target_mH2:g} attempt={n_attempt_rounds} tag={tag} "
                f"tanbeta={tan_beta:.6g} lambda6={lambda6:.3g} NO_ROWS\n"
            )
        logf.flush()

    def do_call(M2_lo, M2_hi, n_M2, tan_beta, lambda6, tag):
        nonlocal n_probes_total, n_attempt_rounds, last_rows_for_reason
        n_attempt_rounds += 1
        rows = run_grid(target_mH2, M2_lo, M2_hi, n_M2, tan_beta, lambda6,
                         mA, mHp, campaign_id, f"{attempt_log_prefix}_{tag}",
                         scratch)
        n_probes_total += len(rows)
        last_rows_for_reason = rows
        valid = [r for r in rows if is_valid(r)]
        center = M2_lo if M2_lo == M2_hi else 0.5 * (M2_lo + M2_hi)
        best = pick_nearest(valid, predicted if predicted is not None else center) if valid \
            else (min(rows, key=lambda r: abs(float(r["M2_input_GeV2"]) - center)) if rows else None)
        log_attempt(tag, tan_beta, lambda6, best, bool(valid))
        return rows, valid

    def accept(row, moved_coordinate, relaxation_used, tan_beta, lambda6,
               pred_err):
        M2 = float(row["M2_input_GeV2"])
        logf.write(
            f"mH2={target_mH2:g} attempt={n_attempt_rounds} VALID "
            f"M2={M2:.6g} ctau={row['ctau_mm']} BRbb={row['br_bb']}\n"
        )
        logf.flush()
        return TrajectoryPoint(target_mH2, M2, tan_beta, lambda6, row,
                                n_attempt_rounds, n_probes_total,
                                moved_coordinate, relaxation_used, pred_err)

    # --- attempt 1: previous accepted parameters, unchanged ---
    if history:
        prev = history[-1]
        rows, valid = do_call(prev.M2, prev.M2, 1, cur_tan_beta, cur_lambda6, "same")
        if valid:
            return accept(valid[0], "none", "none", cur_tan_beta, cur_lambda6, 0.0)

    # --- attempts 2..: local M2 window, expanding ---
    center = predicted if predicted is not None else target_mH2 ** 2
    half_width = max(0.05 * abs(predicted), 0.25 * target_mH2 ** 2) \
        if predicted is not None else 0.5 * target_mH2 ** 2

    window_valid = None
    for expansion in range(MAX_EXPANSIONS + 1):
        rows, valid = do_call(center - half_width, center + half_width,
                               LOCAL_N, cur_tan_beta, cur_lambda6,
                               f"win{expansion}")
        if valid:
            window_valid = valid
            break
        half_width *= EXPANSION_FACTOR

    if window_valid:
        best = pick_nearest(window_valid, predicted if predicted is not None else center)
        best_M2 = float(best["M2_input_GeV2"])
        refine_half = max(half_width / LOCAL_N * 3.0, 1.0)
        rows, valid = do_call(best_M2 - refine_half, best_M2 + refine_half,
                               REFINE_N, cur_tan_beta, cur_lambda6, "refine")
        if valid:
            best = pick_nearest(valid, predicted if predicted is not None else best_M2)
        pred_err = (abs(float(best["M2_input_GeV2"]) - predicted) / abs(predicted)
                    if predicted not in (None, 0.0) else 0.0)
        return accept(best, "M2", "none", cur_tan_beta, cur_lambda6, pred_err)

    # --- rung 4: tan_beta relaxation ---
    for factor in TAN_BETA_RELAX_FACTORS:
        tb_try = cur_tan_beta * factor
        rows, valid = do_call(center - half_width / (EXPANSION_FACTOR ** MAX_EXPANSIONS),
                               center + half_width / (EXPANSION_FACTOR ** MAX_EXPANSIONS),
                               LOCAL_N, tb_try, cur_lambda6, f"tb_x{factor:g}")
        if valid:
            best = pick_nearest(valid, predicted if predicted is not None else center)
            return accept(best, "tan_beta,M2", "tan_beta", tb_try, cur_lambda6, None)

    # --- rung 5: lambda6 relaxation ---
    for l6_try in LAMBDA6_RELAX_VALUES:
        rows, valid = do_call(center - half_width / (EXPANSION_FACTOR ** MAX_EXPANSIONS),
                               center + half_width / (EXPANSION_FACTOR ** MAX_EXPANSIONS),
                               LOCAL_N, cur_tan_beta, l6_try, f"l6_{l6_try:g}")
        if valid:
            best = pick_nearest(valid, predicted if predicted is not None else center)
            return accept(best, "lambda6,M2", "lambda6", cur_tan_beta, l6_try, None)

    # --- total failure for this mass ---
    reason = dominant_rejection(last_rows_for_reason)
    logf.write(f"mH2={target_mH2:g} attempt={n_attempt_rounds} FAILED "
               f"dominant_rejection={reason}\n")
    logf.flush()
    return TrajectoryPoint(target_mH2, None, cur_tan_beta, cur_lambda6, None,
                            n_attempt_rounds, n_probes_total, "none", "none",
                            None, failed=True)


# --- output row construction --------------------------------------------------

TRAJECTORY_HEADER = [
    "schema_version", "producer", "campaign_id",
    "m_H2_GeV", "m_A_GeV", "m_Hp_GeV", "Delta_heavy_GeV",
    "mh_GeV", "sin_beta_minus_alpha", "tan_beta", "lambda6", "lambda7", "lambda1",
    "M2_GeV2", "m12_sq_GeV2",
    "construction_ok", "numerical_ok", "positivity_ok", "unitarity_ok",
    "perturbativity_ok", "stability_ok", "theory_ok",
    "rejection_stage", "rejection_reason",
    "g_hH2H2_GeV", "total_width_GeV", "ctau_physical_mm",
    "BR_bb", "BR_tautau", "BR_WW", "BR_ZZ", "BR_gg", "BR_gammagamma", "BR_hh", "BR_tt",
    "H2_to_AZ_open", "H2_to_HpW_open", "H2_to_AA_open", "H2_to_HpHm_open",
    "H2_to_hh_open", "H2_to_tt_open",
    "n_attempt_rounds", "n_M2_probes", "moved_coordinate", "relaxation_used",
    "M2_prediction_error",
]


def build_row(tp, mA, mHp, campaign_id):
    mH2 = tp.mH2
    common = {
        "schema_version": "h2_m2_continuation.v1",
        "producer": "DihiggsPointV2Evaluator",
        "campaign_id": campaign_id,
        "m_H2_GeV": f"{mH2:.6f}",
        "m_A_GeV": f"{mA:.6f}",
        "m_Hp_GeV": f"{mHp:.6f}",
        "Delta_heavy_GeV": f"{mA - mH2:.6f}",
        "mh_GeV": f"{MH:.6f}",
        "sin_beta_minus_alpha": f"{SIN_BA:.6f}",
        "tan_beta": f"{tp.tan_beta!r}",
        "lambda6": f"{tp.lambda6!r}",
        "lambda7": f"{LAMBDA7:.6f}",
        "n_attempt_rounds": tp.n_attempt_rounds,
        "n_M2_probes": tp.n_M2_probes,
        "moved_coordinate": tp.moved_coordinate,
        "relaxation_used": tp.relaxation_used,
        "M2_prediction_error": (f"{tp.M2_prediction_error:.6g}"
                                 if tp.M2_prediction_error is not None else ""),
        "H2_to_AZ_open": mH2 > mA + MZ,
        "H2_to_HpW_open": mH2 > mHp + MW,
        "H2_to_AA_open": mH2 > 2 * mA,
        "H2_to_HpHm_open": mH2 > 2 * mHp,
        "H2_to_hh_open": mH2 > 2 * MH,
        "H2_to_tt_open": mH2 > 2 * 172.5,
    }
    if tp.failed or tp.row is None:
        common.update({
            "M2_GeV2": "", "m12_sq_GeV2": "", "lambda1": "",
            "construction_ok": "0", "numerical_ok": "", "positivity_ok": "",
            "unitarity_ok": "", "perturbativity_ok": "", "stability_ok": "",
            "theory_ok": "0", "rejection_stage": "FAILED_ALL_RUNGS",
            "rejection_reason": "no_valid_point_found_by_continuation",
            "g_hH2H2_GeV": "", "total_width_GeV": "", "ctau_physical_mm": "",
            "BR_bb": "", "BR_tautau": "", "BR_WW": "", "BR_ZZ": "", "BR_gg": "",
            "BR_gammagamma": "", "BR_hh": "", "BR_tt": "",
        })
    else:
        r = tp.row
        common.update({
            "M2_GeV2": r["M2_input_GeV2"],
            "m12_sq_GeV2": r["m12_sq_reconstructed_GeV2"],
            "lambda1": r["lambda1_reconstructed"],
            "construction_ok": r["construction_ok"],
            "numerical_ok": r["numerical_ok"],
            "positivity_ok": r["positivity_reported_ok"],
            "unitarity_ok": r["unitarity_ok"],
            "perturbativity_ok": r["perturbativity_ok"],
            "stability_ok": r["stability_reported_ok"],
            "theory_ok": r["theory_ok_v1"],
            "rejection_stage": r["rejection_stage"],
            "rejection_reason": r["rejection_reason"],
            "g_hH2H2_GeV": r["g_hH2H2_GeV"],
            "total_width_GeV": r["total_width_GeV"],
            "ctau_physical_mm": r["ctau_mm"],
            "BR_bb": r["br_bb"], "BR_tautau": r["br_tautau"],
            "BR_WW": r["br_WW"], "BR_ZZ": r["br_ZZ"], "BR_gg": r["br_gg"],
            "BR_gammagamma": r["br_gammagamma"], "BR_hh": r["br_hh"],
            "BR_tt": r["br_tt"],
        })
    return common


def fmt_elapsed(seconds):
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def run_slice(mA, mHp, mH2_start, mH2_stop, campaign_id, workdir, results_dir, logf):
    workdir.mkdir(parents=True, exist_ok=True)
    logf.write(f"\n===== SLICE mA=mHp={mA:g} GeV =====\n")
    logf.flush()

    baseline_tan_beta, boot_valid = bootstrap_tan_beta(mA, mHp, mH2_start, campaign_id, workdir, logf)
    baseline_lambda6 = LAMBDA6_SEED
    seed_best = pick_nearest(boot_valid, mH2_start ** 2)

    logf.write(f"[seed-cal mA={mA:g}] baseline_tan_beta={baseline_tan_beta!r} "
               f"baseline_lambda6={baseline_lambda6!r} "
               f"seed_M2={seed_best['M2_input_GeV2']}\n")
    logf.flush()

    seed_tp = TrajectoryPoint(mH2_start, float(seed_best["M2_input_GeV2"]),
                               baseline_tan_beta, baseline_lambda6, seed_best,
                               n_attempt_rounds=len(boot_valid) and 1 or 1,
                               n_M2_probes=401, moved_coordinate="tan_beta"
                               if baseline_tan_beta != TAN_BETA_REQUESTED else "none",
                               relaxation_used="tan_beta(seed_calibration)"
                               if baseline_tan_beta != TAN_BETA_REQUESTED else "none",
                               M2_prediction_error=0.0)

    history = [seed_tp]
    rows_out = [build_row(seed_tp, mA, mHp, campaign_id)]
    logf.write(f"mH2={mH2_start:g} attempt=1 VALID(seed) M2={seed_tp.M2:.6g} "
               f"tanbeta={baseline_tan_beta:.6g} ctau={seed_best['ctau_mm']} "
               f"BRbb={seed_best['br_bb']}\n")
    logf.flush()

    t_start = time.time()
    last_valid_mH2 = mH2_start
    recent_outcomes = []

    masses = list(range(int(mH2_start) + 1, int(mH2_stop) + 1))
    for i, mH2 in enumerate(masses):
        tp = try_mass(float(mH2), mA, mHp, history, baseline_tan_beta,
                       baseline_lambda6, campaign_id, workdir, logf,
                       attempt_log_prefix=f"m{mH2}")
        rows_out.append(build_row(tp, mA, mHp, campaign_id))
        if not tp.failed:
            history.append(tp)
            last_valid_mH2 = mH2
            recent_outcomes.append(1)
        else:
            recent_outcomes.append(0)
        recent_outcomes = recent_outcomes[-20:]
        frac = sum(recent_outcomes) / len(recent_outcomes)
        elapsed = fmt_elapsed(time.time() - t_start)
        logf.write(
            f"[mA={mA:g}] target={mH2} last_valid={last_valid_mH2} "
            f"attempts={tp.n_attempt_rounds} recent_success={frac:.2f} "
            f"elapsed={elapsed}\n"
        )
        logf.flush()

    out_csv = results_dir / f"trajectory_mA{int(mA)}.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=TRAJECTORY_HEADER)
        w.writeheader()
        for row in rows_out:
            w.writerow(row)
    logf.write(f"[mA={mA:g}] DONE. wrote {len(rows_out)} rows to {out_csv}\n")
    logf.flush()
    return out_csv, history


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", nargs="+", type=float, default=[300.0, 350.0])
    ap.add_argument("--mH2-start", type=float, default=126.0)
    ap.add_argument("--results-dir", type=Path,
                     default=REPO_ROOT / "results" / "h2_m2_continuation")
    ap.add_argument("--workdir", type=Path,
                     default=REPO_ROOT / "results" / "h2_m2_continuation" / "_scratch")
    args = ap.parse_args()

    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.workdir.mkdir(parents=True, exist_ok=True)
    log_path = args.results_dir / "attempts.log"

    histories = {}
    with open(log_path, "a", encoding="utf-8") as logf:
        logf.write(f"\n########## RUN {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} "
                    f"commit={git_commit()} dirty={git_dirty()} ##########\n")
        for mA in args.slices:
            mH2_stop = mA - 1.0
            _, history = run_slice(mA, mA, args.mH2_start, mH2_stop,
                                    campaign_id="h2_m2_continuation_v1_mh12520",
                                    workdir=args.workdir,
                                    results_dir=args.results_dir, logf=logf)
            histories[mA] = history

    print("Done. Trajectories written to", args.results_dir)


if __name__ == "__main__":
    main()
