#!/usr/bin/env python3
"""Verify ctau_mm = hbar_c_GeV_mm / total_width for every width_ok row in the
post-Yukawa-fix engineering pilots, against the invariant declared in
hep_cross/contracts/model_point_to_llp_recast_v1.yaml.

Usage: python3 benchmarks/verify_pilot_ctau_invariant.py
Exits non-zero if any width_ok=1 row violates the invariant beyond rel_tol.
"""
import csv
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
]


def check_file(repo_root: Path, rel_path: str, width_col: str, ctau_col: str, ok_col: str, id_col: str) -> list[str]:
    failures = []
    checked = 0
    with open(repo_root / rel_path, newline="") as f:
        for row in csv.DictReader(f):
            if row[ok_col] in ("", "nan") or float(row[ok_col]) != 1.0:
                continue
            width = float(row[width_col])
            ctau = float(row[ctau_col])
            predicted = HBAR_C_GEV_MM / width
            rel_err = abs(predicted - ctau) / ctau
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
