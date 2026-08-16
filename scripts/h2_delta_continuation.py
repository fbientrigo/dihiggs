#!/usr/bin/env python3
"""Fixed-Delta_heavy M2 continuation in mH2, from the canonical LLP anchor.

Starts at the frozen, verified anchor point `H2scan_mH150_tb300000`
(mH2=150, mA=mHp=450, tan_beta=300000, lambda6=1e-10, lambda7=0,
sin(beta-alpha)=1, Type I, mh=125.13 GeV -- see
benchmarks/H2scan_mH150_tb300000_production_coupling.json) and increases
mH2 while holding Delta_heavy = mA - mH2 = mHp - mH2 fixed at 300 GeV.

Primary continuation invariants (never changed automatically):
    tan_beta, lambda6, lambda7, sin(beta-alpha), Yukawa type, mh.
Only M2 is searched at each new target mass. A mass that cannot be reached
by moving M2 alone is recorded FAILED -- it is physics information, not a
software retry target. Relaxing tan_beta/lambda6, or varying Delta_heavy,
only happens in the explicitly separate Phase 2 / secondary diagnostics
below, never silently merged into the primary trajectory.

Calls the existing canonical DihiggsPointV2Evaluator (schema
dihiggs.point.v2) as a subprocess for every physics evaluation. Does not
modify, wrap, or reimplement 2HDMC physics. Does not touch
h2_m2_continuation.py (a separate, already-completed investigation of
fixed *absolute* mA=300/350 slices with a different mh -- not the
question this script answers).
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
BINARY = REPO_ROOT / "dihiggs" / "app" / "DihiggsPointV2Evaluator"

# --- fixed physics invariants (the canonical anchor family) ----------------
MH = 125.13
SIN_BA = 1.0
LAMBDA6 = 1e-10
LAMBDA7 = 0.0
YUKAWA_TYPE = 1
TAN_BETA = 300000.0

# --- canonical anchor (H2scan_mH150_tb300000) -------------------------------
ANCHOR_MH2 = 150.0
ANCHOR_DELTA = 300.0
ANCHOR_M2 = 22499.9999995003345
ANCHOR_EXPECTED = {
    "g_hH2H2_GeV": 63.5914252007596588,
    "total_width_GeV": 4.56118529862185007e-14,
    "ctau_mm": 4.32622152973311191,
    "br_bb": 0.756737485808578692,
}

VALID_TOKEN = "1.00000000000000000e+00"

# --- continuation ladder tuning ---------------------------------------------
COARSE_STEP = 5.0
FINE_STEP = 1.0
LOCAL_N = 21
REFINE_N = 31
MAX_EXPANSIONS = 4
EXPANSION_FACTOR = 2.0
MAX_CONSECUTIVE_FINE_FAILURES = 5


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


def evaluate_grid(mH2, M2_lo, M2_hi, n_M2, mA, mHp, campaign_id, run_id,
                   out_path, tan_beta=TAN_BETA, lambda6=LAMBDA6, mh=MH,
                   sin_ba=SIN_BA, lambda7=LAMBDA7, yukawa_type=YUKAWA_TYPE):
    """Call the canonical evaluator once and return its rows (list of dict)."""
    if n_M2 == 1 and M2_lo != M2_hi:
        M2_hi = M2_lo
    cmd = [
        str(BINARY),
        "--campaign-id", campaign_id, "--run-id", run_id,
        "--mh", repr(mh),
        "--mH-min", repr(mH2), "--mH-max", repr(mH2), "--n-mH", "1",
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


def evaluate_point(mH2, M2, mA, mHp, workdir, campaign_id="scratch",
                    run_id="single", tan_beta=TAN_BETA, lambda6=LAMBDA6, **kw):
    """Evaluate exactly one (mH2, M2) point; returns the single row dict."""
    out_path = workdir / "_scratch_single.csv"
    rows = evaluate_grid(mH2, M2, M2, 1, mA, mHp, campaign_id, run_id,
                          out_path, tan_beta=tan_beta, lambda6=lambda6, **kw)
    return rows[0]


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
    """Among several valid candidate rows, choose the one nearest `target`."""
    if target is None:
        target = 0.0
    return min(valid_rows, key=lambda r: abs(float(r["M2_input_GeV2"]) - target))


GF = 1.16637e-5  # Fermi constant, GeV^-2 (matches benchmarks/check_H2scan_mH150_tb300000.py)


def analytic_M2_lambda1_one(mH2, mh=MH, sin_ba=SIN_BA, lambda6=LAMBDA6,
                             lambda7=LAMBDA7, tan_beta=TAN_BETA):
    """Closed-form M2 that reconstructs lambda1_target=1.0 exactly, given the
    other physics invariants held fixed. This is how the canonical anchor
    itself was constructed (see benchmarks/H2scan_mH150_tb300000_production_
    coupling.json: "lambda1_target": 1.0) and empirically it is also the
    *only* place in M2-space where theory validity holds at tan_beta=3e5:
    a direct scan found the valid M2 band around the anchor is only
    ~1e-6 GeV^2 wide (relative width ~1e-10) -- far too thin for any grid
    search to hit by chance. Positivity/unitarity at this near-decoupling,
    huge-tan_beta corner effectively pins lambda1 to ~1, so tracking this
    analytic curve *is* the adaptive continuation for this branch."""
    beta = math.atan(tan_beta)
    cb, tb = math.cos(beta), math.tan(beta)
    alpha = -math.asin(sin_ba) + beta
    sa, ca = math.sin(alpha), math.cos(alpha)
    v2 = 1.0 / (math.sqrt(2.0) * GF)
    bracket = 1.0 + 1.5 * lambda6 * tb - 0.5 * lambda7 * tb ** 3
    m12_sq = (mH2 ** 2 * ca ** 2 + mh ** 2 * sa ** 2 - v2 * cb ** 2 * bracket) / tb
    return m12_sq / (math.sin(beta) * math.cos(beta))


def predict_M2(history, target_mH2, tan_beta=TAN_BETA, lambda6=LAMBDA6):
    """Primary predictor: the analytic lambda1=1 curve (see
    analytic_M2_lambda1_one), evaluated at the *caller's* tan_beta/lambda6 --
    important for the secondary diagnostic, which perturbs those away from
    their primary-trajectory values. Linear extrapolation from history is
    not used: it is not accurate enough at this tan_beta to land inside the
    theory-valid band (see run.log PHASE0 exploration)."""
    return analytic_M2_lambda1_one(target_mH2, tan_beta=tan_beta, lambda6=lambda6)


def predict_M2_linear(history, target_mH2):
    if not history:
        return None
    if len(history) == 1:
        return history[-1]["M2"]
    a, b = history[-2], history[-1]
    if b["mH2"] == a["mH2"]:
        return b["M2"]
    slope = (b["M2"] - a["M2"]) / (b["mH2"] - a["mH2"])
    return b["M2"] + slope * (target_mH2 - b["mH2"])


def search_local_M2(target_mH2, mA, mHp, history, workdir, logf,
                     attempt_tag, tan_beta=TAN_BETA, lambda6=LAMBDA6,
                     campaign_id="continuation"):
    """M2-only search for a theory-valid point at target_mH2. Never touches
    tan_beta/lambda6. Returns (accepted_row_or_None, n_probes,
    n_attempt_rounds, moved_coordinate, dominant_rejection_reason).

    Ladder: (1) the analytic lambda1=1 prediction, evaluated exactly --
    this is the real continuation for this branch, see predict_M2; (2) the
    previous accepted M2 unchanged, in case the branch is locally flat;
    (3) a widening grid search around the analytic prediction, starting at
    an extremely small relative half-width because the theory-valid band at
    this tan_beta was measured (see run.log PHASE0/exploration) to be only
    ~1e-10 relative wide -- a coarse grid would silently miss it."""
    scratch = workdir / "_scratch_ladder.csv"
    n_probes = 0
    n_rounds = 0
    last_rows = []
    predicted = predict_M2(history, target_mH2, tan_beta=tan_beta, lambda6=lambda6)

    def call(M2_lo, M2_hi, n_M2, tag):
        nonlocal n_probes, n_rounds, last_rows
        n_rounds += 1
        rows = evaluate_grid(target_mH2, M2_lo, M2_hi, n_M2, mA, mHp,
                              campaign_id, f"{attempt_tag}_{tag}", scratch,
                              tan_beta=tan_beta, lambda6=lambda6)
        n_probes += len(rows)
        last_rows = rows
        valid = [r for r in rows if is_valid(r)]
        center = M2_lo if M2_lo == M2_hi else 0.5 * (M2_lo + M2_hi)
        best = pick_nearest(valid, predicted if predicted is not None else center) \
            if valid else None
        status = "OK" if best is not None else (
            f"FAIL:{rows[-1]['rejection_stage']}" if rows else "NO_ROWS")
        logf.write(f"mH2={target_mH2:g} attempt={n_rounds} tag={tag} "
                    f"M2={float((best or rows[0])['M2_input_GeV2']):.6g} "
                    f"theory={status}\n" if rows else
                    f"mH2={target_mH2:g} attempt={n_rounds} tag={tag} NO_ROWS\n")
        logf.flush()
        return rows, valid

    # 1: analytic lambda1=1 prediction, evaluated exactly
    rows, valid = call(predicted, predicted, 1, "analytic")
    if valid:
        return valid[0], n_probes, n_rounds, "M2(analytic)", "none"

    # 2: previous accepted M2, unchanged
    if history:
        prev_M2 = history[-1]["M2"]
        rows, valid = call(prev_M2, prev_M2, 1, "same")
        if valid:
            return valid[0], n_probes, n_rounds, "none", "none"

    # 3: widening grid search around the analytic prediction, starting at a
    # tiny relative half-width (the valid band is razor-thin at this
    # tan_beta) and expanding geometrically toward an ordinary-sized window.
    center = predicted
    half_width = max(1e-8 * abs(center), 1e-6)
    window_valid = None
    for expansion in range(MAX_EXPANSIONS + 1):
        rows, valid = call(center - half_width, center + half_width, LOCAL_N,
                            f"win{expansion}")
        if valid:
            window_valid = valid
            break
        half_width *= 100.0

    if window_valid:
        best = pick_nearest(window_valid, predicted if predicted is not None else center)
        best_M2 = float(best["M2_input_GeV2"])
        refine_half = max(half_width / LOCAL_N * 3.0, 0.5)
        rows, valid = call(best_M2 - refine_half, best_M2 + refine_half, REFINE_N,
                            "refine")
        if valid:
            best = pick_nearest(valid, predicted if predicted is not None else best_M2)
        return best, n_probes, n_rounds, "M2", "none"

    reason = dominant_rejection(last_rows)
    logf.write(f"mH2={target_mH2:g} attempt={n_rounds} FAILED "
                f"dominant_rejection={reason}\n")
    logf.flush()
    return None, n_probes, n_rounds, "none", reason


TRAJECTORY_HEADER = [
    "schema_version", "producer", "campaign_id",
    "m_H2_GeV", "m_A_GeV", "m_Hp_GeV", "Delta_heavy_GeV",
    "mh_GeV", "sin_beta_minus_alpha", "tan_beta", "lambda6", "lambda7",
    "lambda1", "M2_GeV2", "m12_sq_GeV2",
    "construction_ok", "theory_ok", "rejection_stage", "rejection_reason",
    "g_hH2H2_GeV", "total_width_GeV", "ctau_physical_mm",
    "BR_bb", "BR_tautau", "BR_WW", "BR_ZZ", "BR_gg", "BR_gammagamma",
    "BR_hh", "BR_tt",
    "H2_to_AZ_open", "H2_to_HpW_open", "H2_to_AA_open", "H2_to_HpHm_open",
    "H2_to_hh_open", "H2_to_tt_open",
    "grid_mode", "n_attempt_rounds", "n_M2_probes", "moved_coordinate",
]

MZ = 91.15349
MW = 80.36951
MTOP = 172.5


def build_row(mH2, mA, mHp, row, campaign_id, grid_mode, n_rounds, n_probes,
              moved_coordinate, rejection_reason, tan_beta=TAN_BETA,
              lambda6=LAMBDA6):
    common = {
        "schema_version": "h2_delta_continuation.v1",
        "producer": "DihiggsPointV2Evaluator",
        "campaign_id": campaign_id,
        "m_H2_GeV": f"{mH2:.6f}",
        "m_A_GeV": f"{mA:.6f}",
        "m_Hp_GeV": f"{mHp:.6f}",
        "Delta_heavy_GeV": f"{mA - mH2:.6f}",
        "mh_GeV": f"{MH:.6f}",
        "sin_beta_minus_alpha": f"{SIN_BA:.6f}",
        "tan_beta": f"{tan_beta!r}",
        "lambda6": f"{lambda6!r}",
        "lambda7": f"{LAMBDA7:.6f}",
        "grid_mode": grid_mode,
        "n_attempt_rounds": n_rounds,
        "n_M2_probes": n_probes,
        "moved_coordinate": moved_coordinate,
        "H2_to_AZ_open": mH2 > mA + MZ,
        "H2_to_HpW_open": mH2 > mHp + MW,
        "H2_to_AA_open": mH2 > 2 * mA,
        "H2_to_HpHm_open": mH2 > 2 * mHp,
        "H2_to_hh_open": mH2 > 2 * MH,
        "H2_to_tt_open": mH2 > 2 * MTOP,
    }
    if row is None:
        common.update({
            "M2_GeV2": "", "m12_sq_GeV2": "", "lambda1": "",
            "construction_ok": "0", "theory_ok": "0",
            "rejection_stage": "FAILED_ALL_RUNGS",
            "rejection_reason": rejection_reason,
            "g_hH2H2_GeV": "", "total_width_GeV": "", "ctau_physical_mm": "",
            "BR_bb": "", "BR_tautau": "", "BR_WW": "", "BR_ZZ": "", "BR_gg": "",
            "BR_gammagamma": "", "BR_hh": "", "BR_tt": "",
        })
    else:
        common.update({
            "M2_GeV2": row["M2_input_GeV2"],
            "m12_sq_GeV2": row["m12_sq_reconstructed_GeV2"],
            "lambda1": row["lambda1_reconstructed"],
            "construction_ok": row["construction_ok"],
            "theory_ok": row["theory_ok_v1"],
            "rejection_stage": row["rejection_stage"],
            "rejection_reason": row["rejection_reason"],
            "g_hH2H2_GeV": row["g_hH2H2_GeV"],
            "total_width_GeV": row["total_width_GeV"],
            "ctau_physical_mm": row["ctau_mm"],
            "BR_bb": row["br_bb"], "BR_tautau": row["br_tautau"],
            "BR_WW": row["br_WW"], "BR_ZZ": row["br_ZZ"], "BR_gg": row["br_gg"],
            "BR_gammagamma": row["br_gammagamma"], "BR_hh": row["br_hh"],
            "BR_tt": row["br_tt"],
        })
    return common


def verify_anchor(workdir, logf, tol_rel=1e-9):
    row = evaluate_point(ANCHOR_MH2, ANCHOR_M2, ANCHOR_MH2 + ANCHOR_DELTA,
                          ANCHOR_MH2 + ANCHOR_DELTA, workdir,
                          campaign_id="anchor_verify", run_id="phase0")
    ok = row["construction_ok"] == "1" and row["theory_ok_v1"] == VALID_TOKEN \
        and row["width_ok"] == VALID_TOKEN
    diffs = {}
    for key, expected in ANCHOR_EXPECTED.items():
        got = float(row[key])
        diffs[key] = abs(got - expected) / max(abs(expected), 1e-300)
        ok = ok and diffs[key] < tol_rel
    logf.write(f"[PHASE0] anchor construction_ok={row['construction_ok']} "
                f"theory_ok_v1={row['theory_ok_v1']} width_ok={row['width_ok']} "
                f"diffs={diffs} PASS={ok}\n")
    logf.flush()
    return ok, row, diffs


def continue_trajectory(delta_heavy, mH2_start, M2_start, mH2_hard_cap,
                         campaign_id, workdir, results_dir, logf,
                         attempts_jsonl):
    seed_mA = mH2_start + delta_heavy
    seed_row = evaluate_point(mH2_start, M2_start, seed_mA, seed_mA, workdir,
                               campaign_id=campaign_id, run_id="seed")
    assert is_valid(seed_row), "seed point must be theory-valid before continuation starts"

    history = [{"mH2": mH2_start, "M2": M2_start}]
    rows_out = [build_row(mH2_start, seed_mA, seed_mA, seed_row, campaign_id,
                           "seed", 1, 1, "none", "")]
    attempts_jsonl.write(json.dumps({
        "mH2": mH2_start, "Delta_heavy": delta_heavy, "grid_mode": "seed",
        "M2": M2_start, "theory": "PASS", "ctau_mm": seed_row["ctau_mm"],
    }) + "\n")
    attempts_jsonl.flush()

    logf.write(f"mH2={mH2_start:g} Delta={delta_heavy:g} attempt=1 SEED "
                f"M2={M2_start:.6g} theory=PASS ctau={seed_row['ctau_mm']}\n")
    logf.flush()

    mode = "coarse"
    step = COARSE_STEP
    target = mH2_start + step
    last_valid = mH2_start
    fine_return_target = None
    consecutive_fine_failures = 0
    t_start = time.time()
    boundary_mH2 = None
    boundary_reason = None
    first_failure_mH2 = None

    while target <= mH2_hard_cap:
        mA = target + delta_heavy
        mHp = mA
        best, n_probes, n_rounds, moved, reason = search_local_M2(
            target, mA, mHp, history, workdir, logf,
            attempt_tag=f"m{target:g}", campaign_id=campaign_id)

        attempts_jsonl.write(json.dumps({
            "mH2": target, "Delta_heavy": delta_heavy, "grid_mode": mode,
            "M2": (float(best["M2_input_GeV2"]) if best else None),
            "theory": "PASS" if best else f"FAIL:{reason}",
            "n_attempt_rounds": n_rounds, "n_M2_probes": n_probes,
        }) + "\n")
        attempts_jsonl.flush()

        rows_out.append(build_row(target, mA, mHp, best, campaign_id, mode,
                                   n_rounds, n_probes, moved, reason))

        if best is not None:
            M2_val = float(best["M2_input_GeV2"])
            history.append({"mH2": target, "M2": M2_val})
            last_valid = target
            consecutive_fine_failures = 0
            logf.write(f"mH2={target:g} Delta={delta_heavy:g} attempt={n_rounds} "
                        f"M2={M2_val:.6g} theory=PASS ctau={best['ctau_mm']} "
                        f"mode={mode} elapsed={time.time()-t_start:.1f}s\n")
            logf.flush()
            if mode == "fine" and fine_return_target is not None \
                    and target >= fine_return_target:
                mode = "coarse"
                step = COARSE_STEP
                target = last_valid + step
                fine_return_target = None
                continue
            target += step
        else:
            if first_failure_mH2 is None or target < first_failure_mH2:
                first_failure_mH2 = target
            logf.write(f"mH2={target:g} Delta={delta_heavy:g} attempt={n_rounds} "
                        f"theory=FAIL:{reason} mode={mode} "
                        f"elapsed={time.time()-t_start:.1f}s\n")
            logf.flush()
            if mode == "coarse":
                mode = "fine"
                step = FINE_STEP
                fine_return_target = target
                target = last_valid + step
                consecutive_fine_failures = 0
            else:
                consecutive_fine_failures += 1
                target += step
                if consecutive_fine_failures >= MAX_CONSECUTIVE_FINE_FAILURES:
                    boundary_mH2 = last_valid
                    boundary_reason = reason
                    logf.write(f"[BOUNDARY] Delta={delta_heavy:g} "
                                f"last_valid_mH2={last_valid:g} "
                                f"dominant_rejection={boundary_reason}\n")
                    logf.flush()
                    break

    out_csv = results_dir / f"trajectory_delta{int(delta_heavy)}.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=TRAJECTORY_HEADER)
        w.writeheader()
        for row in rows_out:
            w.writerow(row)
    logf.write(f"[Delta={delta_heavy:g}] DONE. wrote {len(rows_out)} rows to {out_csv}. "
                f"last_valid_mH2={last_valid:g} first_failure_mH2={first_failure_mH2}\n")
    logf.flush()
    return out_csv, history, last_valid, boundary_mH2, boundary_reason, first_failure_mH2


def delta_heavy_diagnostic(target_mH2, deltas, seed_history, workdir,
                            logf, campaign_id, results_dir):
    """Phase 2: local diagnostic at the first mH2 where fixed Delta=300
    genuinely failed. For each candidate Delta_heavy, search M2 locally
    (no further mH2 continuation) to see whether the failure is about
    mH2 itself or about the chosen heavy-state separation. `seed_history`
    (typically the last accepted trajectory point) only supplies the M2
    prediction anchor for the local search -- it is not re-evaluated."""
    rows_out = []
    for delta in deltas:
        mA = target_mH2 + delta
        pseudo_history = [seed_history] if seed_history else []
        best, n_probes, n_rounds, moved, reason = search_local_M2(
            target_mH2, mA, mA, pseudo_history, workdir, logf,
            attempt_tag=f"phase2_delta{delta:g}", campaign_id=campaign_id)
        status = "PASS" if best else f"FAIL:{reason}"
        logf.write(f"[PHASE2] mH2={target_mH2:g} Delta={delta:g} mA={mA:g} "
                    f"{status}\n")
        logf.flush()
        rows_out.append(build_row(target_mH2, mA, mA, best, campaign_id,
                                   "phase2", n_rounds, n_probes, moved, reason))
    out_csv = results_dir / "phase2_delta_diagnostic.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=TRAJECTORY_HEADER)
        w.writeheader()
        for row in rows_out:
            w.writerow(row)
    return out_csv, rows_out


def secondary_diagnostic_tan_beta_lambda6(target_mH2, delta_heavy, workdir,
                                           logf, campaign_id, results_dir):
    """SECONDARY DIAGNOSTIC -- explicitly separate from the primary
    fixed-parameter trajectory (Phase 1) and from the Phase 2 Delta_heavy
    scan. Only invoked once a boundary has already been documented at fixed
    tan_beta=300000, lambda6=1e-10. Tests small relative variations in
    tan_beta and lambda6 (one knob at a time) at the fixed Delta=300 mA/mHp
    to answer: can the LLP branch be rescued nearby at the boundary mass?
    This never overwrites or merges into trajectory_delta300.csv."""
    mA = target_mH2 + delta_heavy
    rows_out = []

    def probe(tag, tan_beta, lambda6):
        best, n_probes, n_rounds, moved, reason = search_local_M2(
            target_mH2, mA, mA, [], workdir, logf,
            attempt_tag=f"secondary_{tag}", tan_beta=tan_beta, lambda6=lambda6,
            campaign_id=campaign_id)
        status = "PASS" if best else f"FAIL:{reason}"
        logf.write(f"[SECONDARY] mH2={target_mH2:g} Delta={delta_heavy:g} "
                    f"tag={tag} tan_beta={tan_beta:.6g} lambda6={lambda6:.3g} "
                    f"{status}\n")
        logf.flush()
        out_row = build_row(target_mH2, mA, mA, best, campaign_id, "secondary",
                             n_rounds, n_probes, moved, reason,
                             tan_beta=tan_beta, lambda6=lambda6)
        rows_out.append(out_row)
        return best is not None

    tan_beta_factors = [0.999, 0.99, 0.9, 0.5, 1.001, 1.01, 1.1, 2.0]
    for factor in tan_beta_factors:
        probe(f"tb_x{factor:g}", TAN_BETA * factor, LAMBDA6)

    lambda6_values = [5e-11, 2e-10, 1e-9, 1e-11, -1e-10]
    for l6 in lambda6_values:
        probe(f"l6_{l6:g}", TAN_BETA, l6)

    out_csv = results_dir / "secondary_diagnostic_tan_beta_lambda6.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=TRAJECTORY_HEADER)
        w.writeheader()
        for row in rows_out:
            w.writerow(row)
    return out_csv, rows_out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mH2-hard-cap", type=float, default=900.0)
    ap.add_argument("--results-dir", type=Path,
                     default=REPO_ROOT / "results" / "h2_delta_continuation")
    ap.add_argument("--workdir", type=Path,
                     default=REPO_ROOT / "results" / "h2_delta_continuation" / "_scratch")
    ap.add_argument("--skip-phase2", action="store_true")
    ap.add_argument("--run-secondary-diagnostic", action="store_true",
                     help="Explicitly separate tan_beta/lambda6 rescue probe "
                          "at the boundary mass. Only meaningful after a "
                          "boundary has been documented; never runs by default.")
    args = ap.parse_args()

    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.workdir.mkdir(parents=True, exist_ok=True)
    log_path = args.results_dir / "run.log"
    attempts_path = args.results_dir / "attempts.jsonl"

    with open(log_path, "a", encoding="utf-8") as logf, \
         open(attempts_path, "a", encoding="utf-8") as attempts_jsonl:
        logf.write(f"\n########## RUN {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} "
                    f"commit={git_commit()} dirty={git_dirty()} ##########\n")
        logf.flush()

        ok, row, diffs = verify_anchor(args.workdir, logf)
        if not ok:
            logf.write("[PHASE0] ANCHOR VERIFICATION FAILED -- STOPPING.\n")
            logf.flush()
            print("ANCHOR VERIFICATION FAILED. See run.log. Stopping per mission contract.")
            sys.exit(1)
        print("Phase 0: anchor reproduced. Starting Phase 1 continuation.")

        out_csv, history, last_valid, boundary_mH2, boundary_reason, first_failure_mH2 = \
            continue_trajectory(
                delta_heavy=ANCHOR_DELTA, mH2_start=ANCHOR_MH2, M2_start=ANCHOR_M2,
                mH2_hard_cap=args.mH2_hard_cap,
                campaign_id="h2_delta_continuation_v1_mh12513",
                workdir=args.workdir, results_dir=args.results_dir, logf=logf,
                attempts_jsonl=attempts_jsonl)

        print(f"Phase 1 done. last_valid_mH2={last_valid}, boundary_mH2={boundary_mH2}, "
              f"first_failure_mH2={first_failure_mH2}, reason={boundary_reason}")

        if first_failure_mH2 is not None and not args.skip_phase2:
            seed_hist = history[-1] if history else None
            deltas = [250.0, 300.0, 350.0, 400.0]
            phase2_csv, phase2_rows = delta_heavy_diagnostic(
                first_failure_mH2, deltas, seed_hist, args.workdir, logf,
                "h2_delta_continuation_v1_mh12513_phase2", args.results_dir)
            print(f"Phase 2 diagnostic written to {phase2_csv}")

        if first_failure_mH2 is not None and args.run_secondary_diagnostic:
            sec_csv, sec_rows = secondary_diagnostic_tan_beta_lambda6(
                first_failure_mH2, ANCHOR_DELTA, args.workdir, logf,
                "h2_delta_continuation_v1_mh12513_secondary", args.results_dir)
            print(f"Secondary (tan_beta/lambda6) diagnostic written to {sec_csv}")

    manifest = {
        "schema": "h2_delta_continuation.manifest.v1",
        "producer": "DihiggsPointV2Evaluator",
        "git_commit": git_commit(),
        "git_dirty": git_dirty(),
        "anchor": {
            "mH2": ANCHOR_MH2, "Delta_heavy": ANCHOR_DELTA, "M2_GeV2": ANCHOR_M2,
            "expected": ANCHOR_EXPECTED,
        },
        "fixed_invariants": {
            "mh_GeV": MH, "sin_beta_minus_alpha": SIN_BA, "lambda6": LAMBDA6,
            "lambda7": LAMBDA7, "yukawa_type": YUKAWA_TYPE, "tan_beta": TAN_BETA,
        },
        "phase1_last_valid_mH2": last_valid,
        "phase1_boundary_mH2": boundary_mH2,
        "phase1_first_failure_mH2": first_failure_mH2,
        "phase1_boundary_reason": boundary_reason,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    with open(args.results_dir / "manifest.json", "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)

    print("Done. Results written to", args.results_dir)


if __name__ == "__main__":
    main()
