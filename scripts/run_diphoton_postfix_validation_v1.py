#!/usr/bin/env python3
"""Run the bounded post-Yukawa-fix diphoton validation campaign.

This is intentionally a small bookkeeping/orchestration script.  It does not
change the canonical evaluator or construct replacement physics points.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO.parent
EVAL = REPO / "dihiggs" / "app" / "DihiggsPointV2Evaluator"
OUT = ROOT / "runs" / "diphoton_postfix_validation_v1"
CAMPAIGN = "diphoton_postfix_validation_v1"

MH = 125.13
SBA = 1.0
L1 = 1.0
L7 = 0.0
YUKAWA_TYPE = 1
GF = 1.16637e-5
L6_HIST = 1.0e-2
TB_HIST = 1.0e5
L6_LLP = 1.0e-10
TB_LLP = 3.0e5
MASSES = [150.0, 200.0, 250.0]

EXPERIMENTS = [
    ("exp01_sqrt10_local", "local_half_decade_sqrt10", [10.0 ** -0.5, 1.0, 10.0 ** 0.5], [10.0 ** -0.5, 1.0, 10.0 ** 0.5]),
    ("exp02_log_decade", "broad_log_decade", [0.1, 1.0, 10.0], [0.1, 1.0, 10.0]),
    ("exp03_linear_local", "local_linear", [0.5, 1.0, 1.5], [0.5, 1.0, 1.5]),
    ("exp04_historical_to_llp_bridge", "historical_to_llp_bridge", [1.0, math.sqrt(L6_LLP / L6_HIST), L6_LLP / L6_HIST], [1.0, math.sqrt(TB_LLP / TB_HIST), TB_LLP / TB_HIST]),
]

LAMBDA1_CLOSURE_TOLERANCE = 2.0e-4

HISTORICAL_SOURCE = "dihiggs/benchmarks/run_first_h2_bounded_scan.py; dihiggs/benchmarks/first_h2_bounded_scan_manifest.json; git history commit 6bfad7662fd87750d838bf2fe0bd7ac00ee2326a ancestry"


def sci(value: float) -> str:
    # C++ scientific + setprecision(max_digits10=17): 17 significant digits.
    return format(float(value), ".16e")


def faithful_m12_sq(mh: float, mH: float, tan_beta: float, lambda6: float, lambda7: float, lambda1: float) -> float:
    """The historical set_param_phys_lam1 inversion, copied as bookkeeping."""
    beta = math.atan(tan_beta)
    cb = math.cos(beta)
    tb = math.tan(beta)
    alpha = -math.asin(SBA) + beta
    sa = math.sin(alpha)
    ca = math.cos(alpha)
    v2 = 1.0 / (math.sqrt(2.0) * GF)
    bracket = lambda1 + 1.5 * lambda6 * tb - 0.5 * lambda7 * tb ** 3
    return (mH * mH * ca * ca + mh * mh * sa * sa - v2 * cb * cb * bracket) / tb


def point_inputs(lambda6: float, tan_beta: float, mH: float) -> dict[str, float]:
    mA = mH + 300.0
    m12 = faithful_m12_sq(MH, mH, tan_beta, lambda6, L7, L1)
    beta = math.atan(tan_beta)
    M2 = m12 / (math.sin(beta) * math.cos(beta))
    return {
        "m_h": MH,
        "m_H2": mH,
        "m_A": mA,
        "m_Hplus": mA,
        "sin_beta_minus_alpha": SBA,
        "tan_beta": tan_beta,
        "lambda6": lambda6,
        "lambda7": L7,
        "lambda1": L1,
        "M2_GeV2": M2,
        "m12_sq_GeV2": m12,
        "yukawa_type": YUKAWA_TYPE,
    }


def physical_key(p: dict[str, float]) -> tuple[str, ...]:
    return tuple(sci(p[k]) for k in ("m_h", "m_H2", "m_A", "m_Hplus", "sin_beta_minus_alpha", "tan_beta", "M2_GeV2", "lambda6", "lambda7", "yukawa_type"))


def fnv1a64(text: str) -> str:
    h = 14695981039346656037
    for b in text.encode():
        h ^= b
        h = (h * 1099511628211) & ((1 << 64) - 1)
    return f"point_{h:016x}"


def expected_point_id(p: dict[str, float]) -> str:
    return fnv1a64(",".join(sci(p[k]) for k in ("m_h", "m_H2", "m_A", "m_Hplus", "sin_beta_minus_alpha", "tan_beta", "M2_GeV2", "lambda6", "lambda7", "yukawa_type")))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def experiment_grid(exp_id: str, l6_scales: list[float], tb_scales: list[float]) -> tuple[list[float], list[float]]:
    lambda6_values = [L6_HIST * scale for scale in l6_scales]
    tan_beta_values = [TB_HIST * scale for scale in tb_scales]
    if exp_id == "exp04_historical_to_llp_bridge":
        expected_l6 = [1.0e-2, 1.0e-6, 1.0e-10]
        if not all(math.isclose(a, b, rel_tol=1.0e-14, abs_tol=0.0) for a, b in zip(lambda6_values, expected_l6)):
            raise AssertionError(f"Experiment 4 lambda6 bridge mismatch: {lambda6_values}")
        # Normalize only binary-rounding residue so the requested physical
        # inputs, and therefore point identities, are the exact decimal grid.
        lambda6_values = expected_l6
        expected_tb = [1.0e5, math.sqrt(3.0) * 1.0e5, 3.0e5]
        if tan_beta_values != expected_tb:
            raise AssertionError(f"Experiment 4 tan_beta bridge mismatch: {tan_beta_values}")
    return lambda6_values, tan_beta_values


def print_grids() -> None:
    total = 0
    for exp_id, label, l6_scales, tb_scales in EXPERIMENTS:
        l6_values, tb_values = experiment_grid(exp_id, l6_scales, tb_scales)
        print(json.dumps({"experiment_id": exp_id, "label": label, "lambda6": l6_values, "tan_beta": tb_values, "m_H2_GeV": MASSES, "slots": 27}, separators=(",", ":")))
        total += 27
    if total != 108:
        raise AssertionError(f"requested grid slots={total}, expected 108")
    if [L6_HIST * s for s in EXPERIMENTS[3][2]][:1] != [1.0e-2]:
        raise AssertionError("historical anchor changed")


def read_one(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise RuntimeError(f"expected exactly one evaluator row in {path}, got {len(rows)}")
    return rows[0]


def f(row: dict[str, str], *names: str) -> float:
    for name in names:
        if name in row:
            try:
                return float(row[name])
            except ValueError:
                return math.nan
    return math.nan


def b(row: dict[str, str], *names: str) -> int | str:
    for name in names:
        if name in row:
            value = row[name]
            if value in {"1", "1.0", "1.00000000000000000e+00"}:
                return 1
            if value in {"0", "0.0", "0.00000000000000000e+00"}:
                return 0
            return value
    return "MISSING"


def ratio(value: float, anchor: float) -> float:
    if not (math.isfinite(value) and math.isfinite(anchor)) or anchor == 0.0:
        return math.nan
    return value / anchor


def build_designs() -> tuple[dict[str, Any], list[dict[str, Any]], dict[tuple[str, ...], dict[str, Any]]]:
    slots: list[dict[str, Any]] = []
    unique: dict[tuple[str, ...], dict[str, Any]] = {}
    for exp_id, label, l6_scales, tb_scales in EXPERIMENTS:
        exp_dir = OUT / exp_id
        exp_dir.mkdir(parents=True, exist_ok=True)
        lambda6_values, tan_beta_values = experiment_grid(exp_id, l6_scales, tb_scales)
        design = {
            "experiment_id": exp_id,
            "label": label,
            "purpose": {
                "exp01_sqrt10_local": "moderately broad multiplicative neighborhood",
                "exp02_log_decade": "one full decade around the historical point",
                "exp03_linear_local": "local linear sensitivity",
                "exp04_historical_to_llp_bridge": "historical diphoton point to LLP trajectory",
            }[exp_id],
            "grid": {
                "lambda6": lambda6_values,
                "tan_beta": tan_beta_values,
                "m_H2_GeV": MASSES,
            },
            "grid_scales": {"lambda6": l6_scales, "tan_beta": tb_scales},
            "fixed_physics_contract": {
                "m_h_GeV": MH,
                "m_A_prescription": "m_A = m_Hplus = m_H2 + 300 GeV",
                "m_Hplus_prescription": "m_Hplus = m_A",
                "sin_beta_minus_alpha": SBA,
                "lambda1_prescription": "lambda1 = 1; solve set_param_phys_lam1 inversion for m12_sq, then M2 = m12_sq/(sin(beta)cos(beta))",
                "lambda1": L1,
                "lambda7": L7,
                "yukawa_type": YUKAWA_TYPE,
                "GF_GeV_minus_2": GF,
            },
            "planned_slots": 27,
        }
        write_json(exp_dir / "design.json", design)
        for li, l6 in enumerate(design["grid"]["lambda6"]):
            for ti, tb in enumerate(design["grid"]["tan_beta"]):
                for mi, mass in enumerate(MASSES):
                    p = point_inputs(l6, tb, mass)
                    key = physical_key(p)
                    if key not in unique:
                        unique[key] = {"point_id": expected_point_id(p), "inputs": p, "key": key, "experiments": [], "slots": []}
                    slot = {
                        "experiment_id": exp_id,
                        "lambda6_index": li,
                        "tan_beta_index": ti,
                        "mass_index": mi,
                        "point_id": unique[key]["point_id"],
                        **p,
                    }
                    slots.append(slot)
                    unique[key]["experiments"].append(exp_id)
                    unique[key]["slots"].append(slot)
        write_csv(exp_dir / "points_27.csv", [s for s in slots if s["experiment_id"] == exp_id])
    manifest = {
        "campaign_id": CAMPAIGN,
        "scope": "108 requested slots, deduplicated by complete physical input tuple; no MadGraph/Pythia/recast/detector simulation",
        "historical_anchor": {"lambda6": L6_HIST, "tan_beta": TB_HIST},
        "experiments": [e[0] for e in EXPERIMENTS],
        "requested_grid_slots": len(slots),
        "unique_physical_points": len(unique),
        "reused_points": len(slots) - len(unique),
        "masses_GeV": MASSES,
        "evaluator": str(EVAL),
        "evaluator_sha256": sha256(EVAL),
        "evaluator_source": "dihiggs/dihiggs/src/DihiggsPointV2Evaluator.cpp",
        "evaluator_source_sha256": sha256(REPO / "dihiggs" / "src" / "DihiggsPointV2Evaluator.cpp"),
        "yukawa_fix_commit": "6bfad7662fd87750d838bf2fe0bd7ac00ee2326a",
        "yukawa_fix_ancestry_check": "git merge-base --is-ancestor 6bfad7662fd87750d838bf2fe0bd7ac00ee2326a HEAD",
    }
    write_json(OUT / "manifest.json", manifest)
    if len(slots) != 108 or len(unique) != 99:
        raise AssertionError(f"campaign design expected 108 slots and 99 unique points, got {len(slots)} and {len(unique)}")
    return manifest, slots, unique


def run_unique(unique: dict[tuple[str, ...], dict[str, Any]], manifest: dict[str, Any]) -> list[dict[str, Any]]:
    raw_dir = OUT / "raw_points"
    raw_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    env = os.environ.copy()
    env["DIHIGGS_GIT_COMMIT"] = subprocess.check_output(["git", "-C", str(ROOT / "dihiggs"), "rev-parse", "HEAD"], text=True).strip()
    env["DIHIGGS_GIT_DIRTY"] = "no"
    for index, item in enumerate(unique.values(), start=1):
        p = item["inputs"]
        expected_id = item["point_id"]
        temporary_output = raw_dir / f".tmp_unique_{index:03d}.csv"
        if not temporary_output.exists():
            cmd = [
                str(EVAL), "--campaign-id", CAMPAIGN, "--run-id", f"unique_{index:03d}",
                "--mh", sci(p["m_h"]), "--mH-min", sci(p["m_H2"]), "--mH-max", sci(p["m_H2"]), "--n-mH", "1",
                "--mA", sci(p["m_A"]), "--mHp", sci(p["m_Hplus"]), "--yukawa-type", str(YUKAWA_TYPE),
                "--sin-ba", sci(p["sin_beta_minus_alpha"]), "--tan-beta", sci(p["tan_beta"]),
                "--M2-min", sci(p["M2_GeV2"]), "--M2-max", sci(p["M2_GeV2"]), "--n-M2", "1",
                "--lambda6", sci(p["lambda6"]), "--lambda7", sci(p["lambda7"]), "--output", str(temporary_output),
            ]
            proc = subprocess.run(cmd, cwd=ROOT / "dihiggs", env=env, text=True, capture_output=True)
            if proc.returncode != 0:
                raise RuntimeError(f"evaluator failed for {expected_id}: {proc.stdout}\n{proc.stderr}")
        row = read_one(temporary_output)
        point_id = row.get("point_id", "")
        if not point_id:
            raise RuntimeError(f"evaluator did not emit point_id for {expected_id}")
        output = raw_dir / f"{point_id}.csv"
        if output.exists():
            output.unlink()
        temporary_output.rename(output)
        item["point_id"] = point_id
        for slot in item["slots"]:
            slot["point_id"] = point_id
        derived = {
            "point_id": point_id,
            **p,
            "raw_file": str(output.relative_to(OUT)),
            "requested_experiment_ids": ";".join(sorted(set(item["experiments"]))),
            "requested_slot_count": len(item["slots"]),
            "construction_ok": b(row, "construction_ok"),
            "positivity_ok": b(row, "positivity_reported_ok", "positivity_ok"),
            "unitarity_ok": b(row, "unitarity_ok"),
            "perturbativity_ok": b(row, "perturbativity_ok"),
            "theory_ok": b(row, "theory_ok_v1", "theory_ok"),
            "width_ok": b(row, "width_ok"),
            "rejection_stage": row.get("rejection_stage", "MISSING"),
            "rejection_reason": row.get("rejection_reason", "MISSING"),
            "HB_HS_classification": "MISSING",
            "lambda1_target": L1,
            "lambda1_reconstructed": f(row, "lambda1_reconstructed"),
            "lambda1_closure_error": abs(f(row, "lambda1_reconstructed") - L1),
            "Gamma_total": f(row, "total_width_GeV", "total_width_gev"),
            "Gamma_bb": f(row, "width_bb_GeV", "width_bb_gev"),
            "Gamma_cc": f(row, "width_cc_GeV", "width_cc_gev"),
            "Gamma_tt": f(row, "width_tt_GeV", "width_tt_gev"),
            "Gamma_tautau": f(row, "width_tautau_GeV", "width_tautau_gev"),
            "Gamma_gg": f(row, "width_gg_GeV", "width_gg_gev"),
            "Gamma_gamma_gamma": f(row, "width_gammagamma_GeV", "width_gammagamma_gev"),
            "Gamma_Z_gamma": f(row, "width_Zgamma_GeV", "width_Zgamma_gev"),
            "Gamma_WW": f(row, "width_WW_GeV", "width_WW_gev"),
            "Gamma_ZZ": f(row, "width_ZZ_GeV", "width_ZZ_gev"),
            "Gamma_hh": f(row, "width_hh_GeV", "width_hh_gev"),
            "BR_bb": f(row, "br_bb"),
            "BR_cc": f(row, "br_cc"),
            "BR_tt": f(row, "br_tt"),
            "BR_tautau": f(row, "br_tautau"),
            "BR_gg": f(row, "br_gg"),
            "BR_gamma_gamma": f(row, "br_gammagamma"),
            "BR_Z_gamma": f(row, "br_Zgamma"),
            "BR_WW": f(row, "br_WW"),
            "BR_ZZ": f(row, "br_ZZ"),
            "BR_hh": f(row, "br_hh"),
            "BR_tt": f(row, "br_tt"),
            "ctau_mm": f(row, "ctau_mm"),
            "g_hH2H2": f(row, "g_hH2H2_GeV"),
            "g_H2HpHm": "MISSING",
            "status": "VALID" if b(row, "construction_ok") == 1 and b(row, "theory_ok_v1", "theory_ok") == 1 and b(row, "width_ok") == 1 else "INVALID",
        }
        rows.append(derived)
    return rows


def rewrite_authoritative_slot_files(slots: list[dict[str, Any]]) -> None:
    """Rewrite slot manifests after evaluator point IDs become authoritative."""
    for exp_id, _, _, _ in EXPERIMENTS:
        write_csv(OUT / exp_id / "points_27.csv", [s for s in slots if s["experiment_id"] == exp_id])
    write_csv(OUT / "comparison" / "all_requested_slots.csv", slots)


def enrich(rows: list[dict[str, Any]]) -> None:
    anchors = {(r["m_H2"], r["tan_beta"], r["lambda6"]): r for r in rows if r["tan_beta"] == TB_HIST and r["lambda6"] == L6_HIST}
    for r in rows:
        a = anchors.get((r["m_H2"], TB_HIST, L6_HIST))
        for name, source in [("BR_gamma_gamma", "BR_gamma_gamma"), ("Gamma_gamma_gamma", "Gamma_gamma_gamma"), ("Gamma_total", "Gamma_total"), ("ctau_mm", "ctau_mm")]:
            r[f"R_{name}_vs_postfix_anchor"] = ratio(r[name], a[source]) if a else math.nan
        r["theory_ok"] = r["theory_ok"]


def make_experiment_outputs(slots: list[dict[str, Any]], rows_by_id: dict[str, dict[str, Any]]) -> None:
    table_fields = ["m_H2", "m_A", "M2_GeV2", "m12_sq_GeV2", "Gamma_total", "ctau_mm", "Gamma_gamma_gamma", "BR_gamma_gamma", "BR_Z_gamma", "BR_bb", "BR_gg", "g_hH2H2", "g_H2HpHm", "theory_ok", "width_ok", "status"]
    matrix_fields = ["BR_gamma_gamma", "Gamma_gamma_gamma", "BR_Z_gamma", "BR_bb", "Gamma_total", "ctau_mm", "g_hH2H2", "g_H2HpHm", "theory_ok"]
    for exp_id, _, _, _ in EXPERIMENTS:
        exp_slots = [s for s in slots if s["experiment_id"] == exp_id]
        exp_dir = OUT / exp_id
        tables = exp_dir / "tables"
        matrices = exp_dir / "matrices"
        tables.mkdir(exist_ok=True)
        matrices.mkdir(exist_ok=True)
        by_cell: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
        for s in exp_slots:
            by_cell[(s["lambda6_index"], s["tan_beta_index"])].append(rows_by_id[s["point_id"]])
        for (li, ti), cell in sorted(by_cell.items()):
            cell = sorted(cell, key=lambda r: r["m_H2"])
            csv_path = tables / f"cell_l6_{li}_tb_{ti}.csv"
            write_csv(csv_path, [{k: r.get(k, "MISSING") for k in table_fields} for r in cell], table_fields)
            lines = [f"# {exp_id} cell lambda6_index={li}, tan_beta_index={ti}", "", "| " + " | ".join(table_fields) + " |", "|" + "|".join(["---"] * len(table_fields)) + "|"]
            for r in cell:
                lines.append("| " + " | ".join(str(r.get(k, "MISSING")) for k in table_fields) + " |")
            (tables / f"cell_l6_{li}_tb_{ti}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
        for mass in MASSES:
            for field in matrix_fields:
                path = matrices / f"mH2_{int(mass)}_{field}.csv"
                exp_design = next(e for e in EXPERIMENTS if e[0] == exp_id)
                l6_values = [L6_HIST * x for x in exp_design[2]]
                tb_values = [TB_HIST * x for x in exp_design[3]]
                mat_rows = []
                for li, l6 in enumerate(l6_values):
                    row = {"lambda6": l6}
                    for ti, tb in enumerate(tb_values):
                        slot = next(s for s in exp_slots if s["mass_index"] == MASSES.index(mass) and s["lambda6_index"] == li and s["tan_beta_index"] == ti)
                        row[f"tan_beta_{ti}"] = rows_by_id[slot["point_id"]].get(field, "MISSING")
                    mat_rows.append(row)
                write_csv(path, mat_rows)
        valid = [rows_by_id[s["point_id"]] for s in exp_slots if rows_by_id[s["point_id"]]["status"] == "VALID"]
        pool = [rows_by_id[s["point_id"]] for s in exp_slots]
        summary = {
            "experiment_id": exp_id,
            "requested_slots": len(exp_slots),
            "valid": len(valid),
            "invalid": len(exp_slots) - len(valid),
            "BR_gamma_gamma_range": [min((r["BR_gamma_gamma"] for r in valid), default=math.nan), max((r["BR_gamma_gamma"] for r in valid), default=math.nan)],
            "Gamma_gamma_gamma_range": [min((r["Gamma_gamma_gamma"] for r in valid), default=math.nan), max((r["Gamma_gamma_gamma"] for r in valid), default=math.nan)],
            "ctau_mm_range": [min((r["ctau_mm"] for r in valid), default=math.nan), max((r["ctau_mm"] for r in valid), default=math.nan)],
            "status_counts": {"VALID": len(valid), "INVALID": len(exp_slots) - len(valid)},
        }
        write_json(exp_dir / "status_summary.md.json", summary)
        (exp_dir / "status_summary.md").write_text("# " + exp_id + "\n\n" + json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")


def make_comparisons(manifest: dict[str, Any], slots: list[dict[str, Any]], unique: dict[tuple[str, ...], dict[str, Any]], rows: list[dict[str, Any]]) -> None:
    by_id = {r["point_id"]: r for r in rows}
    write_csv(OUT / "point_cache.csv", rows)
    write_csv(OUT / "comparison" / "all_unique_points.csv", rows)
    write_json(OUT / "comparison" / "duplicate_reuse_report.json", {
        "requested_grid_slots": len(slots), "unique_physical_points": len(rows), "reused_points": len(slots) - len(rows),
        "actual_evaluator_calls": len(rows),
        "duplicates": [{"point_id": item["point_id"], "slot_count": len(item["slots"]), "experiments": sorted(set(item["experiments"]))} for item in unique.values() if len(item["slots"]) > 1],
    })
    valid = [r for r in rows if r["status"] == "VALID"]
    best = sorted(valid, key=lambda r: r["BR_gamma_gamma"], reverse=True)
    candidate_fields = ["point_id", "m_H2", "lambda6", "tan_beta", "M2_GeV2", "m12_sq_GeV2", "lambda1_target", "lambda1_reconstructed", "lambda1_closure_error", "BR_gamma_gamma", "Gamma_gamma_gamma", "Gamma_total", "ctau_mm", "BR_bb", "BR_gg", "BR_Z_gamma", "g_hH2H2", "g_H2HpHm", "theory_ok", "width_ok", "status"]
    write_csv(OUT / "comparison" / "best_gamma_candidates.csv", [{k: r.get(k, "MISSING") for k in candidate_fields} for r in best[:10]], candidate_fields)
    llp = [r for r in valid if r["ctau_mm"] >= 0.1]
    llp = sorted(llp, key=lambda r: (r["BR_gamma_gamma"], r["Gamma_gamma_gamma"]), reverse=True)
    llp_fields = candidate_fields + ["ctau_ge_1mm", "ctau_ge_0p1mm"]
    write_csv(OUT / "comparison" / "llp_gamma_candidates.csv", [{**{k: r.get(k, "MISSING") for k in candidate_fields}, "ctau_ge_1mm": r["ctau_mm"] >= 1.0, "ctau_ge_0p1mm": r["ctau_mm"] >= 0.1} for r in llp[:10]], llp_fields)
    failures = [r for r in rows if r["status"] == "INVALID"]
    write_csv(OUT / "comparison" / "failures.csv", failures)
    exp_rows = []
    for exp_id, _, _, _ in EXPERIMENTS:
        exp = [by_id[s["point_id"]] for s in slots if s["experiment_id"] == exp_id]
        exp_valid = [r for r in exp if r["status"] == "VALID"]
        exp_rows.append({"experiment_id": exp_id, "requested_slots": 27, "valid": len(exp_valid), "invalid": 27 - len(exp_valid), "BR_gamma_gamma_min": min((r["BR_gamma_gamma"] for r in exp_valid), default=math.nan), "BR_gamma_gamma_max": max((r["BR_gamma_gamma"] for r in exp_valid), default=math.nan), "Gamma_gamma_gamma_min": min((r["Gamma_gamma_gamma"] for r in exp_valid), default=math.nan), "Gamma_gamma_gamma_max": max((r["Gamma_gamma_gamma"] for r in exp_valid), default=math.nan), "ctau_mm_min": min((r["ctau_mm"] for r in exp_valid), default=math.nan), "ctau_mm_max": max((r["ctau_mm"] for r in exp_valid), default=math.nan)})
    write_csv(OUT / "comparison" / "experiment_comparison.csv", exp_rows)


def write_anchor_artifacts(rows: list[dict[str, Any]]) -> None:
    anchor = OUT / "00_historical_anchor"
    anchor.mkdir(parents=True, exist_ok=True)
    write_json(anchor / "recovered_anchor.json", {
        "historical_lambda6": L6_HIST, "historical_tan_beta": TB_HIST, "m_h_GeV": MH,
        "m_A_prescription": "m_A = m_Hplus = m_H2 + 300 GeV", "sin_beta_minus_alpha": SBA,
        "lambda7": L7, "lambda1_prescription": "lambda1=1 via set_param_phys_lam1 inversion", "lambda1": L1,
        "M2_prescription": "M2 = m12_sq / (sin(beta)cos(beta))", "yukawa_type": "Type I (1)",
        "mass_grid_GeV": MASSES, "source_files": HISTORICAL_SOURCE,
        "ambiguity": "The repository preserves the l6=0.01, tan_beta=1e5 convention in historical validation scripts, but no single pre-fix diphoton result table with this exact full slice was found. The pair is therefore the best-supported historical study anchor, not a uniquely archived post-fix benchmark row.",
    })
    evidence = """# Historical evidence\n\n- `dihiggs/benchmarks/run_first_h2_bounded_scan.py` freezes `m_h=125.13`, `sin(beta-alpha)=1`, `lambda1=1`, `lambda7=0`, Type I, and `m_A=m_H+=m_H2+300 GeV`.\n- `dihiggs/benchmarks/check_H2scan_mH150_tb300000.py` records the exact `set_param_phys_lam1` inversion and `M2` conversion.\n- The older lifetime validation script explicitly uses `lambda6=0.01` and documents `tan_beta=10^5` as a study case.\n- The canonical evaluator fix is commit `6bfad7662fd87750d838bf2fe0bd7ac00ee2326a`; current campaign rows are post-fix only.\n- No exact pre-fix diphoton table for the complete requested 150/200/250 slice was found in tracked artifacts; missing values remain labelled reference-only.\n"""
    (anchor / "historical_evidence.md").write_text(evidence, encoding="utf-8")
    old_fields = ["provenance", "m_H2", "Gamma_bb", "Gamma_gg", "Gamma_gamma_gamma", "Gamma_Z_gamma", "Gamma_total", "BR_bb", "BR_gg", "BR_gamma_gamma", "BR_Z_gamma", "ctau_mm"]
    old_rows = [{"provenance": "HISTORICAL_PRE_FIX_REFERENCE_ONLY:MISSING", "m_H2": mass, **{k: "MISSING" for k in old_fields[2:]}} for mass in MASSES]
    write_csv(anchor / "old_values_reference.csv", old_rows, old_fields)
    postfix = [r for r in rows if r["lambda6"] == L6_HIST and r["tan_beta"] == TB_HIST]
    replay_fields = ["m_H2", "M2_GeV2", "m12_sq_GeV2", "Gamma_bb", "Gamma_gg", "Gamma_gamma_gamma", "Gamma_Z_gamma", "Gamma_total", "BR_bb", "BR_gg", "BR_gamma_gamma", "BR_Z_gamma", "ctau_mm", "status", "producer_provenance"]
    write_csv(anchor / "postfix_anchor_replay.csv", [{k: r.get(k, "MISSING") for k in replay_fields} | {"producer_provenance": "POSTFIX_CORRECTED_DihiggsPointV2Evaluator"} for r in sorted(postfix, key=lambda x: x["m_H2"])], replay_fields)
    write_csv(anchor / "old_vs_postfix.csv", [{"m_H2": r["m_H2"], "old_provenance": "HISTORICAL_PRE_FIX_REFERENCE_ONLY:MISSING", "postfix_provenance": "POSTFIX_CORRECTED", **{f"old_{k}": "MISSING" for k in ["Gamma_bb", "Gamma_gg", "Gamma_gamma_gamma", "Gamma_Z_gamma", "Gamma_total", "BR_bb", "BR_gg", "BR_gamma_gamma", "BR_Z_gamma", "ctau_mm"]}, **{f"postfix_{k}": r[k] for k in ["Gamma_bb", "Gamma_gg", "Gamma_gamma_gamma", "Gamma_Z_gamma", "Gamma_total", "BR_bb", "BR_gg", "BR_gamma_gamma", "BR_Z_gamma", "ctau_mm"]}} for r in sorted(postfix, key=lambda x: x["m_H2"])])


def run_prefix_reference() -> None:
    """Replay only the three historical central points with the pre-fix binary."""
    prefix_eval = Path(os.environ.get("DIHIGGS_PREFIX_EVAL", "/tmp/atlas_dihiggs_prefix_ref/dihiggs/app/DihiggsPointV2Evaluator"))
    if not prefix_eval.exists():
        raise RuntimeError(f"missing pre-fix reference evaluator: {prefix_eval}")
    cache_path = OUT / "point_cache.csv"
    with cache_path.open(newline="", encoding="utf-8") as stream:
        current = list(csv.DictReader(stream))
    central = sorted((r for r in current if float(r["lambda6"]) == L6_HIST and float(r["tan_beta"]) == TB_HIST), key=lambda r: float(r["m_H2"]))
    old_rows: list[dict[str, Any]] = []
    for index, r in enumerate(central, start=1):
        output = Path(f"/tmp/diphoton_prefix_reference_{index}.csv")
        p = [
            str(prefix_eval), "--campaign-id", "diphoton_postfix_validation_v1_prefix_reference", "--run-id", f"mH2_{int(float(r['m_H2']))}",
            "--mh", sci(MH), "--mH-min", sci(float(r["m_H2"])), "--mH-max", sci(float(r["m_H2"])), "--n-mH", "1",
            "--mA", sci(float(r["m_A"])), "--mHp", sci(float(r["m_Hplus"])), "--yukawa-type", "1", "--sin-ba", sci(SBA),
            "--tan-beta", sci(TB_HIST), "--M2-min", sci(float(r["M2_GeV2"])), "--M2-max", sci(float(r["M2_GeV2"])), "--n-M2", "1",
            "--lambda6", sci(L6_HIST), "--lambda7", sci(L7), "--output", str(output),
        ]
        proc = subprocess.run(p, cwd=ROOT / "dihiggs", env={**os.environ, "DIHIGGS_GIT_COMMIT": "cbb5079"}, text=True, capture_output=True)
        if proc.returncode != 0:
            raise RuntimeError(f"pre-fix reference failed: {proc.stdout}\n{proc.stderr}")
        old = read_one(output)
        fields = ["Gamma_bb", "Gamma_gg", "Gamma_gamma_gamma", "Gamma_Z_gamma", "Gamma_total", "BR_bb", "BR_gg", "BR_gamma_gamma", "BR_Z_gamma", "ctau_mm"]
        values = {
            "Gamma_bb": f(old, "width_bb_GeV"), "Gamma_gg": f(old, "width_gg_GeV"), "Gamma_gamma_gamma": f(old, "width_gammagamma_GeV"),
            "Gamma_Z_gamma": f(old, "width_Zgamma_GeV"), "Gamma_total": f(old, "total_width_GeV"), "BR_bb": f(old, "br_bb"),
            "BR_gg": f(old, "br_gg"), "BR_gamma_gamma": f(old, "br_gammagamma"), "BR_Z_gamma": f(old, "br_Zgamma"), "ctau_mm": f(old, "ctau_mm"),
        }
        old_rows.append({"provenance": "HISTORICAL_PRE_FIX_REFERENCE_ONLY", "source_commit": "cbb5079", "m_H2": float(r["m_H2"]), **values})
    anchor = OUT / "00_historical_anchor"
    old_fields = ["provenance", "source_commit", "m_H2", "Gamma_bb", "Gamma_gg", "Gamma_gamma_gamma", "Gamma_Z_gamma", "Gamma_total", "BR_bb", "BR_gg", "BR_gamma_gamma", "BR_Z_gamma", "ctau_mm"]
    write_csv(anchor / "old_values_reference.csv", old_rows, old_fields)
    old_by_mass = {r["m_H2"]: r for r in old_rows}
    comparison: list[dict[str, Any]] = []
    for post in central:
        old = old_by_mass[float(post["m_H2"])]
        out = {"m_H2": float(post["m_H2"]), "old_provenance": old["provenance"], "old_source_commit": old["source_commit"], "postfix_provenance": "POSTFIX_CORRECTED", "postfix_source_commit": post.get("producer_commit", "ea2681069d2edbe9b059e7b2e17f3169d5a4653d")}
        for name in ["Gamma_bb", "Gamma_gg", "Gamma_gamma_gamma", "Gamma_Z_gamma", "Gamma_total", "BR_bb", "BR_gg", "BR_gamma_gamma", "BR_Z_gamma", "ctau_mm"]:
            out[f"old_{name}"] = old[name]
            out[f"postfix_{name}"] = float(post[name])
            out[f"old_over_postfix_{name}"] = ratio(old[name], float(post[name]))
        comparison.append(out)
    write_csv(anchor / "old_vs_postfix.csv", comparison)


def validate(manifest: dict[str, Any], slots: list[dict[str, Any]], rows: list[dict[str, Any]]) -> None:
    issues: list[str] = []
    if len(slots) != 108: issues.append(f"requested slots={len(slots)}")
    if len(rows) != manifest["unique_physical_points"]: issues.append("unique count mismatch")
    if sorted({s["m_H2"] for s in slots}) != MASSES: issues.append("mass grid mismatch")
    for exp_id, _, _, _ in EXPERIMENTS:
        if len([s for s in slots if s["experiment_id"] == exp_id]) != 27: issues.append(f"{exp_id} slot count")
    referenced_ids = {s["point_id"] for s in slots}
    raw_paths = sorted((OUT / "raw_points").glob("point_*.csv"))
    raw_ids: list[str] = []
    for raw_path in raw_paths:
        raw_row = read_one(raw_path)
        raw_id = raw_row.get("point_id", "")
        if raw_id != raw_path.stem:
            issues.append(f"raw point filename/id mismatch {raw_path.name}:{raw_id}")
        raw_ids.append(raw_id)
    if len(raw_ids) != len(set(raw_ids)):
        issues.append("duplicate raw point IDs")
    if referenced_ids != set(raw_ids):
        issues.append(f"slot/raw point-id map mismatch referenced={len(referenced_ids)} raw={len(raw_ids)}")
    for r in rows:
        if r["width_ok"] == 1:
            expected = 1.973269804e-13 / r["Gamma_total"]
            if not math.isfinite(r["ctau_mm"]) or abs(r["ctau_mm"] - expected) / expected > 1e-10: issues.append(f"ctau mismatch {r['point_id']}")
        for name in ["BR_bb", "BR_cc", "BR_tt", "BR_tautau", "BR_gg", "BR_gamma_gamma", "BR_Z_gamma", "BR_WW", "BR_ZZ", "BR_hh"]:
            value = r[name]
            if math.isfinite(value) and not 0.0 <= value <= 1.0 + 1e-12: issues.append(f"BR range {r['point_id']}:{name}")
        if r["width_ok"] not in {0, 1, "MISSING"}: issues.append(f"width_ok not preserved {r['point_id']}")
        if not math.isfinite(r["lambda1_reconstructed"]):
            issues.append(f"lambda1 reconstruction missing {r['point_id']}")
        elif r["lambda1_closure_error"] >= LAMBDA1_CLOSURE_TOLERANCE:
            issues.append(f"lambda1 closure {r['point_id']}={r['lambda1_closure_error']}")
    lambda1_reconstructed = [r["lambda1_reconstructed"] for r in rows if math.isfinite(r["lambda1_reconstructed"])]
    lambda1_errors = [r["lambda1_closure_error"] for r in rows if math.isfinite(r["lambda1_closure_error"])]
    report = {"passed": not issues, "issues": issues, "requested_grid_slots": len(slots), "unique_physical_points": len(rows), "actual_evaluator_calls": len(rows), "reused_points": len(slots) - len(rows), "masses_GeV": MASSES, "lambda1_target": L1, "lambda1_closure_tolerance_abs": LAMBDA1_CLOSURE_TOLERANCE, "lambda1_reconstructed": {"min": min(lambda1_reconstructed, default=math.nan), "max": max(lambda1_reconstructed, default=math.nan)}, "lambda1_closure_error": {"max": max(lambda1_errors, default=math.nan), "all_below_tolerance": all(error < LAMBDA1_CLOSURE_TOLERANCE for error in lambda1_errors)}, "lambda1_closure_max_observed": max(lambda1_errors, default=math.nan), "yukawa_fix_commit": manifest["yukawa_fix_commit"], "madgraph_mixed": False, "pythia_mixed": False, "recast_mixed": False, "detector_quantity_mixed": False, "g_H2HpHm_status": "MISSING; no existing canonical evaluator field/API emitted; adding it would require an out-of-scope canonical mechanism change"}
    write_json(OUT / "comparison" / "validation_report.json", report)
    if issues:
        raise RuntimeError("validation failed: " + "; ".join(issues))


def main() -> int:
    if "--dry-run" in sys.argv[1:]:
        print_grids()
        return 0
    if "--reference-only" in sys.argv[1:]:
        run_prefix_reference()
        return 0
    if not EVAL.exists():
        raise SystemExit(f"missing evaluator binary: {EVAL}")
    OUT.mkdir(parents=True, exist_ok=True)
    manifest, slots, unique = build_designs()
    rows = run_unique(unique, manifest)
    rewrite_authoritative_slot_files(slots)
    enrich(rows)
    rows_by_id = {r["point_id"]: r for r in rows}
    write_anchor_artifacts(rows)
    make_experiment_outputs(slots, rows_by_id)
    make_comparisons(manifest, slots, unique, rows)
    validate(manifest, slots, rows)
    manifest["actual_evaluator_calls"] = len(rows)
    manifest["valid_points"] = sum(r["status"] == "VALID" for r in rows)
    manifest["invalid_points"] = sum(r["status"] == "INVALID" for r in rows)
    write_json(OUT / "manifest.json", manifest)
    print(json.dumps({"requested_slots": len(slots), "unique_points": len(rows), "actual_calls": len(rows), "reused_points": len(slots) - len(rows), "valid": manifest["valid_points"], "invalid": manifest["invalid_points"]}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
