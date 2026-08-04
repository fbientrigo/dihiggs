#!/usr/bin/env python3
"""Verify ctau_mm = hbar_c_GeV_mm / total_width for every width_ok row in the
post-Yukawa-fix engineering pilots and the bounded H2 benchmark scan, against
the invariant declared in hep_cross/contracts/model_point_to_llp_recast_v1.yaml.

Usage: python3 benchmarks/verify_pilot_ctau_invariant.py
Exits non-zero if any width_ok=1 row violates the invariant beyond rel_tol.
"""
import csv
import math
import sys
from pathlib import Path

HBAR_C_GEV_MM = 1.973269804e-13
REL_TOL = 1.0e-4

PILOT_FILES = [
    ("docs/pilots/dihiggs_point_v2_v1/L01_accepted_anchor.csv", "total_width_GeV", "ctau_mm", "width_ok", "point_id"),
    ("docs/pilots/dihiggs_point_v2_v1/L05_theory_rejected_anchor.csv", "total_width_GeV", "ctau_mm", "width_ok", "point_id"),
    ("docs/pilots/dihiggs_point_v2_v1/L06_llp_anchor.csv", "total_width_GeV", "ctau_mm", "width_ok", "point_id"),
    ("docs/pilots/dihiggs_point_v2_v1/ordering_boundary.csv", "total_width_GeV", "ctau_mm", "width_ok", "point_id"),
    ("docs/pilots/dihiggs_point_v2_v1/construction_failure.csv", "total_width_GeV", "ctau_mm", "width_ok", "point_id"),
    ("docs/pilots/lambda1_v2_yukawa_fix_v1/lambda1_v2_pilot.csv", "total_width_gev", "ctau_mm", "width_ok", "point_id"),
    ("benchmarks/first_h2_bounded_scan.csv", "total_width_gev", "ctau_mm", "width_ok", "point_id"),
]


def check_file(repo_root: Path, rel_path: str, width_col: str, ctau_col: str, ok_col: str, id_col: str) -> list[str]:
    failures = []
    checked = 0
    with open(repo_root / rel_path, newline="") as f:
        for row in csv.DictReader(f):
            try:
                width_ok = float(row[ok_col])
            except (TypeError, ValueError):
                continue
            if not math.isfinite(width_ok) or width_ok != 1.0:
                continue
            try:
                width = float(row[width_col])
            except (TypeError, ValueError):
                failures.append(f"{rel_path}:{row[id_col]} total width is non-finite: {row[width_col]!r}")
                continue
            if not math.isfinite(width):
                failures.append(f"{rel_path}:{row[id_col]} total width is non-finite: {width!r}")
                continue
            if width <= 0.0:
                failures.append(f"{rel_path}:{row[id_col]} total width <= 0: {width!r}")
                continue
            try:
                ctau = float(row[ctau_col])
            except (TypeError, ValueError):
                failures.append(f"{rel_path}:{row[id_col]} ctau is non-finite: {row[ctau_col]!r}")
                continue
            if not math.isfinite(ctau):
                failures.append(f"{rel_path}:{row[id_col]} ctau is non-finite: {ctau!r}")
                continue
            if ctau <= 0.0:
                failures.append(f"{rel_path}:{row[id_col]} ctau <= 0: {ctau!r}")
                continue
            predicted = HBAR_C_GEV_MM / width
            if not math.isfinite(predicted):
                failures.append(f"{rel_path}:{row[id_col]} predicted ctau is non-finite: {predicted!r}")
                continue
            rel_err = abs(predicted - ctau) / ctau
            if not math.isfinite(rel_err):
                failures.append(f"{rel_path}:{row[id_col]} relative error is non-finite: {rel_err!r}")
                continue
            checked += 1
            if rel_err > REL_TOL:
                failures.append(
                    f"{rel_path}:{row[id_col]} ctau_mm={ctau!r} predicted={predicted!r} rel_err={rel_err!r}"
                )
    print(f"{rel_path}: {checked} width_ok=1 row(s) checked")
    return failures


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    all_failures = []
    for rel_path, width_col, ctau_col, ok_col, id_col in PILOT_FILES:
        all_failures.extend(check_file(repo_root, rel_path, width_col, ctau_col, ok_col, id_col))
    if all_failures:
        print("FAIL: ctau invariant violated:")
        for line in all_failures:
            print(f"  {line}")
        return 1
    print("PASS: ctau_mm == hbar_c_GeV_mm / total_width for all width_ok=1 pilot rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
