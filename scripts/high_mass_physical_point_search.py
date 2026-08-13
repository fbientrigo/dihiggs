#!/usr/bin/env python3
"""PHYSICAL_POINT_SCAN search for theory-valid high-mass H2 points.

Per docs/HIGH_MASS_H2_CONTRACT.md and docs/contracts/high_mass_campaign_template.yaml:
this is the PHYSICAL_POINT_SCAN experiment class (tan_beta, M2, lambda6 allowed
to vary to locate theory-valid points at each (m_H2, Delta_heavy) node), run
progressively region by region as required by the mission. Uses the canonical
DihiggsPointV2Evaluator binary exclusively; does not reimplement any 2HDMC
constraint check. Every attempted point (accepted or rejected) is retained.

Output: writes raw per-call CSVs to a scratch dir, then aggregates into
results/high_mass_valid_points/{high_mass_attempted_points.csv,
high_mass_valid_points.csv} plus rejection_summary.csv.

Not launched from pytest; this is a one-off research script, run manually.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"
RESULTS_DIR = ROOT / "results/high_mass_valid_points"
SCRATCH_DIR = Path("/tmp/claude-1000/-home-fabian-atlas-dihiggs/8d8e0b8a-cfad-43e3-a3a9-72d556a6c6ce/scratchpad/hm_scan_raw")

M_H = 125.13
M_W = 80.36951
M_Z = 91.15349
M_T_POLE = 172.5
HBAR_C_GEV_MM = 1.973269804e-13
YUKAWA_TYPE = 1  # Type I, documented default control point (Gate B)
LAMBDA7 = 0.0    # frozen control coordinate, exact numerical zero (Gate B allows 0.0 or 1e-10)
MAX_WORKERS = 6  # "a handful is fine" per mission constraints; machine has 8 cores

# --- Mass-region ladder (mission section: progressive targets ~150/200/250-300/
# 350-400/500/600/800 GeV). Delta_heavy chosen so mA=mHp<=2000 GeV always. ---
REGIONS = [
    ("R150_anchor", 150.0, 300.0),
    ("R200_below_hh", 200.0, 300.0),
    ("R250_near_hh_threshold", 250.0, 300.0),
    ("R300_above_hh", 300.0, 500.0),
    ("R350_near_tt_threshold", 350.0, 450.0),
    ("R400_above_tt", 400.0, 400.0),
    ("R500", 500.0, 500.0),
    ("R600", 600.0, 600.0),
    ("R800_near2000", 800.0, 1200.0),
]

# --- PHYSICAL_POINT_SCAN search grid (Section 6 of the contract: tan_beta,
# M2, lambda6 vary to locate theory-valid points; sin_ba recorded explicitly
# per point, not silently defaulted). Grid informed by the pilot's own
# finding that the only theory_ok=1 anchor (P0, 150 GeV) used near-exact
# alignment (sin_ba=1.0), very high tan_beta, and tiny lambda6; 2HDMC's
# physical-basis relation v^2*lambda5 = M2 - mA^2 (see THDM manual) makes
# M2 near mA^2 the natural decoupling coordinate, so M2 is swept as a
# fraction of mA^2 using the evaluator's own internal grid capability.
TAN_BETA_GRID = [2.0, 5.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 10000.0, 100000.0, 300000.0]
SIN_BA_GRID = [1.0, 0.999, 0.995]
LAMBDA6_GRID = [0.0, 1e-10, 0.1, -0.1]
M2_FRAC_MIN, M2_FRAC_MAX, M2_N = -1.0, 3.0, 51  # M2 = frac * mA^2, linear grid

WIDTH_FIELDS = (
    "width_bb_GeV", "width_cc_GeV", "width_tt_GeV", "width_tautau_GeV",
    "width_WW_GeV", "width_ZZ_GeV", "width_gammagamma_GeV", "width_Zgamma_GeV",
    "width_gg_GeV", "width_hh_GeV",
)


def git_head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, text=True, capture_output=True
    ).stdout.strip()


def git_dirty() -> str:
    out = subprocess.run(
        ["git", "status", "--porcelain"], cwd=ROOT, check=True, text=True, capture_output=True
    ).stdout
    return "yes" if out.strip() else "no"


def cascade_flags(mH2: float, mA: float, mHp: float) -> dict:
    return {
        "H2_to_AZ_open": mH2 > mA + M_Z,
        "H2_to_HpW_open": mH2 > mHp + M_W,
        "H2_to_AA_open": mH2 > 2 * mA,
        "H2_to_HpHm_open": mH2 > 2 * mHp,
        "H2_to_hh_open": mH2 > 2 * M_H,
        "H2_to_tt_open": mH2 > 2 * M_T_POLE,
    }


FORBIDDEN = ("H2_to_AZ_open", "H2_to_HpW_open", "H2_to_AA_open", "H2_to_HpHm_open")


def build_calls(region_label: str, mH2: float, delta: float, commit: str, dirty: str):
    """Yield (call_id, argv, meta) for every (sin_ba, tan_beta, lambda6) combo
    in this region; each call internally grids M2 via the evaluator's own
    --M2-min/--M2-max/--n-M2 (fraction of mA^2)."""
    mA = mH2 + delta
    calls = []
    idx = 0
    for sba, tb, l6 in itertools.product(SIN_BA_GRID, TAN_BETA_GRID, LAMBDA6_GRID):
        idx += 1
        run_id = f"{region_label}_c{idx:04d}"
        m2_min = M2_FRAC_MIN * mA * mA
        m2_max = M2_FRAC_MAX * mA * mA
        output = SCRATCH_DIR / f"{run_id}.csv"
        argv = [
            str(BINARY), "--campaign-id", "high_mass_h2_physical_point_scan_v1",
            "--run-id", run_id,
            "--mh", repr(M_H), "--mH-min", repr(mH2), "--mH-max", repr(mH2), "--n-mH", "1",
            "--mA", repr(mA), "--mHp", repr(mA), "--yukawa-type", str(YUKAWA_TYPE),
            "--sin-ba", repr(sba), "--tan-beta", repr(tb),
            "--M2-min", repr(m2_min), "--M2-max", repr(m2_max), "--n-M2", str(M2_N),
            "--lambda6", repr(l6), "--lambda7", repr(LAMBDA7), "--output", str(output),
        ]
        calls.append((run_id, argv, output, {
            "region_label": region_label, "mH2_target": mH2, "Delta_heavy_target": delta,
            "mA_target": mA, "sin_ba": sba, "tan_beta": tb, "lambda6": l6,
        }))
    return calls


def run_one(run_id, argv, output, meta, commit, dirty):
    env = {**os.environ, "DIHIGGS_GIT_COMMIT": commit, "DIHIGGS_GIT_DIRTY": dirty}
    proc = subprocess.run(argv, text=True, capture_output=True, env=env)
    rows = []
    if proc.returncode == 0 and output.is_file():
        with output.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                row.update({
                    "mass_region": meta["region_label"],
                    "Delta_heavy_target_GeV": meta["Delta_heavy_target"],
                })
                rows.append(row)
        output.unlink()
    else:
        rows.append({
            "mass_region": meta["region_label"], "Delta_heavy_target_GeV": meta["Delta_heavy_target"],
            "construction_ok": "0", "rejection_stage": "cli_invocation",
            "rejection_reason": f"nonzero_exit_or_missing_output:{proc.returncode}:{proc.stderr.strip()[:200]}",
            "mH_input_GeV": meta["mH2_target"], "mA_input_GeV": meta["mA_target"],
            "mHp_input_GeV": meta["mA_target"], "sin_beta_minus_alpha_input": meta["sin_ba"],
            "tan_beta_input": meta["tan_beta"], "lambda6_input": meta["lambda6"],
            "theory_ok_v1": "nan",
        })
    return rows


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    commit = git_head()
    dirty = git_dirty()
    print(f"producer_commit={commit} producer_dirty={dirty}", file=sys.stderr)

    all_calls = []
    for region_label, mH2, delta in REGIONS:
        all_calls.extend(build_calls(region_label, mH2, delta, commit, dirty))
    print(f"Total CLI invocations: {len(all_calls)}", file=sys.stderr)

    attempted_rows = []
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {
            pool.submit(run_one, run_id, argv, output, meta, commit, dirty): run_id
            for run_id, argv, output, meta in all_calls
        }
        done = 0
        for fut in as_completed(futures):
            rows = fut.result()
            attempted_rows.extend(rows)
            done += 1
            if done % 20 == 0 or done == len(futures):
                elapsed = time.time() - t0
                print(f"[{done}/{len(futures)}] calls done, {len(attempted_rows)} rows, {elapsed:.1f}s elapsed",
                      file=sys.stderr)

    print(f"Total rows collected: {len(attempted_rows)}", file=sys.stderr)

    # --- augment every row with cascade flags + classification scaffolding ---
    for row in attempted_rows:
        try:
            mH2 = float(row["mH_input_GeV"])
            mA = float(row["mA_input_GeV"])
            mHp = float(row["mHp_input_GeV"])
        except (KeyError, ValueError):
            mH2 = mA = mHp = float("nan")
        if math.isfinite(mH2) and math.isfinite(mA) and math.isfinite(mHp):
            flags = cascade_flags(mH2, mA, mHp)
        else:
            flags = {k: None for k in ("H2_to_AZ_open", "H2_to_HpW_open", "H2_to_AA_open",
                                        "H2_to_HpHm_open", "H2_to_hh_open", "H2_to_tt_open")}
        row.update({f"cascade_{k}": v for k, v in flags.items()})
        row["experiment_class"] = "PHYSICAL_POINT_SCAN"
        row["model_variant"] = "PHYSICAL_DECAYS_NO_HEAVY_CASCADES"

        theory_ok = row.get("theory_ok_v1", "nan")
        try:
            theory_ok_bool = float(theory_ok) == 1.0
        except ValueError:
            theory_ok_bool = False
        forbidden_open = any(flags.get(k) for k in FORBIDDEN) if flags.get(FORBIDDEN[0]) is not None else True
        accepted = theory_ok_bool and not forbidden_open
        row["accepted"] = "1" if accepted else "0"

    # --- write attempted_points.csv (native evaluator columns + augmentation) ---
    attempted_csv = RESULTS_DIR / "high_mass_attempted_points.csv"
    fieldnames = []
    seen = set()
    for row in attempted_rows:
        for key in row.keys():
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with attempted_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, restval="")
        writer.writeheader()
        writer.writerows(attempted_rows)
    print(f"Wrote {attempted_csv} ({len(attempted_rows)} rows)", file=sys.stderr)

    # --- build valid_points.csv (accepted==1) ---
    accepted_rows = [r for r in attempted_rows if r["accepted"] == "1"]
    valid_records = []
    for row in accepted_rows:
        mH2 = float(row["mH_input_GeV"])
        mA = float(row["mA_input_GeV"])
        mHp = float(row["mHp_input_GeV"])
        total_width = float(row["total_width_GeV"])
        widths_sum = sum(float(row[f]) for f in WIDTH_FIELDS)
        unaccounted = float(row["width_unaccounted_GeV"])
        closure_ok = math.isclose(widths_sum + unaccounted, total_width, rel_tol=1e-9, abs_tol=1e-30) \
            if total_width != 0 else True
        record = {
            "schema_version": "dihiggs.high_mass_point.v1",
            "producer": row["producer"], "producer_commit": row["producer_commit"],
            "producer_dirty": row["producer_dirty"], "twohdmc_version": "1.8.0",
            "campaign_id": row["campaign_id"], "run_id": row["run_id"], "point_id": row["point_id"],
            "mass_region": row["mass_region"],
            "model_variant": row["model_variant"], "experiment_class": row["experiment_class"],
            "m_h_GeV": float(row["mh_input_GeV"]), "m_H2_GeV": mH2, "m_A_GeV": mA, "m_Hp_GeV": mHp,
            "Delta_heavy_GeV": mA - mH2, "Delta_heavy_target_GeV": row["Delta_heavy_target_GeV"],
            "sin_beta_minus_alpha": float(row["sin_beta_minus_alpha_input"]),
            "tan_beta": float(row["tan_beta_input"]),
            "M2_GeV2": float(row["M2_input_GeV2"]), "m12_sq_GeV2": float(row["m12_sq_input_GeV2"]),
            "lambda1": float(row["lambda1_reconstructed"]), "lambda2": float(row["lambda2_reconstructed"]),
            "lambda3": float(row["lambda3_reconstructed"]), "lambda4": float(row["lambda4_reconstructed"]),
            "lambda5": float(row["lambda5_reconstructed"]), "lambda6": float(row["lambda6_reconstructed"]),
            "lambda7": float(row["lambda7_reconstructed"]),
            "yukawa_type": int(float(row["yukawa_type"])),
            "yukawa_type_installed": row["yukawa_type_installed"],
            "g_hH2H2_GeV": float(row["g_hH2H2_GeV"]),
            "width_bb_GeV": float(row["width_bb_GeV"]), "width_cc_GeV": float(row["width_cc_GeV"]),
            "width_tt_GeV": float(row["width_tt_GeV"]), "width_tautau_GeV": float(row["width_tautau_GeV"]),
            "width_WW_GeV": float(row["width_WW_GeV"]), "width_ZZ_GeV": float(row["width_ZZ_GeV"]),
            "width_gg_GeV": float(row["width_gg_GeV"]), "width_gammagamma_GeV": float(row["width_gammagamma_GeV"]),
            "width_Zgamma_GeV": float(row["width_Zgamma_GeV"]), "width_hh_GeV": float(row["width_hh_GeV"]),
            "total_width_GeV": total_width, "width_unaccounted_GeV": unaccounted,
            "width_closure_ok": closure_ok,
            "BR_bb": float(row["br_bb"]), "BR_cc": float(row["br_cc"]), "BR_tt": float(row["br_tt"]),
            "BR_tautau": float(row["br_tautau"]), "BR_WW": float(row["br_WW"]), "BR_ZZ": float(row["br_ZZ"]),
            "BR_gg": float(row["br_gg"]), "BR_gammagamma": float(row["br_gammagamma"]),
            "BR_Zgamma": float(row["br_Zgamma"]), "BR_hh": float(row["br_hh"]),
            "ctau_physical_mm": float(row["ctau_mm"]),
            "construction_ok": row["construction_ok"], "numerical_ok": row["numerical_ok"],
            "positivity_ok": row["positivity_reported_ok"], "unitarity_ok": row["unitarity_ok"],
            "perturbativity_ok": row["perturbativity_ok"], "stability_ok": row["stability_reported_ok"],
            "theory_ok": row["theory_ok_v1"],
            "rejection_stage": row["rejection_stage"], "rejection_reason": row["rejection_reason"],
            "H2_to_AZ_open": row["cascade_H2_to_AZ_open"], "H2_to_HpW_open": row["cascade_H2_to_HpW_open"],
            "H2_to_AA_open": row["cascade_H2_to_AA_open"], "H2_to_HpHm_open": row["cascade_H2_to_HpHm_open"],
            "H2_to_hh_open": row["cascade_H2_to_hh_open"], "H2_to_tt_open": row["cascade_H2_to_tt_open"],
            "classification": "THEORY_VALID",  # updated later after HB/HS attempt
        }
        valid_records.append(record)

    valid_csv = RESULTS_DIR / "high_mass_valid_points.csv"
    if valid_records:
        vfieldnames = list(valid_records[0].keys())
        with valid_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=vfieldnames)
            writer.writeheader()
            writer.writerows(valid_records)
    else:
        with valid_csv.open("w", newline="", encoding="utf-8") as handle:
            handle.write("# no theory-valid points found\n")
    print(f"Wrote {valid_csv} ({len(valid_records)} accepted points)", file=sys.stderr)

    # --- rejection_summary.csv ---
    from collections import Counter
    counts = Counter()
    for row in attempted_rows:
        key = (row.get("mass_region", "?"), row.get("rejection_stage", "?"), row.get("rejection_reason", "?"))
        counts[key] += 1
    summary_csv = RESULTS_DIR / "rejection_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mass_region", "rejection_stage", "rejection_reason", "count"])
        for (region, stage, reason), n in sorted(counts.items()):
            writer.writerow([region, stage, reason, n])
    print(f"Wrote {summary_csv}", file=sys.stderr)

    # --- search-parameter manifest fragment (consumed by build_manifest.py) ---
    search_manifest = {
        "search_script": "scripts/high_mass_physical_point_search.py",
        "producer_commit": commit, "producer_dirty": dirty,
        "evaluator_binary": str(BINARY.relative_to(ROOT)),
        "regions": [{"label": lbl, "m_H2_GeV": mh2, "Delta_heavy_GeV": d, "m_A_GeV": mh2 + d}
                    for lbl, mh2, d in REGIONS],
        "search_grid": {
            "tan_beta": TAN_BETA_GRID, "sin_beta_minus_alpha": SIN_BA_GRID,
            "lambda6": LAMBDA6_GRID, "lambda7_fixed": LAMBDA7,
            "M2_fraction_of_mA_squared_range": [M2_FRAC_MIN, M2_FRAC_MAX], "M2_grid_points": M2_N,
            "yukawa_type_fixed": YUKAWA_TYPE,
        },
        "total_cli_invocations": len(all_calls),
        "total_attempted_points": len(attempted_rows),
        "total_accepted_points": len(valid_records),
        "wall_time_seconds": round(time.time() - t0, 1),
        "config_hash_sha256": hashlib.sha256(json.dumps({
            "regions": REGIONS, "tan_beta": TAN_BETA_GRID, "sin_ba": SIN_BA_GRID,
            "lambda6": LAMBDA6_GRID, "lambda7": LAMBDA7, "m2_frac": [M2_FRAC_MIN, M2_FRAC_MAX, M2_N],
            "yukawa_type": YUKAWA_TYPE,
        }, sort_keys=True).encode()).hexdigest(),
    }
    with (RESULTS_DIR / "search_manifest_fragment.json").open("w") as handle:
        json.dump(search_manifest, handle, indent=2)
    print(f"Wrote {RESULTS_DIR / 'search_manifest_fragment.json'}", file=sys.stderr)
    print(json.dumps({"accepted": len(valid_records), "attempted": len(attempted_rows)}))


if __name__ == "__main__":
    main()
