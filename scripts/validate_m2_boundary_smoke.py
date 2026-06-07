#!/usr/bin/env python3
"""
validate_m2_boundary_smoke.py

Validates a CSV produced by Phys_M2BoundaryScan against the phys_m2_boundary_v1
contract:

  - header contains M2_input and m12_sq_input
  - header does NOT contain computed_lam1
  - m12_sq_input matches M2_input * sin_beta_cos_beta within double-precision tolerance
  - m12_sq_out is finite for every data row
  - lam1_out and lam2_out are finite for every data row

Usage:
    python3 scripts/validate_m2_boundary_smoke.py <csv_file>

Exit codes:
    0  all checks passed
    1  one or more checks failed
"""

import csv
import math
import sys


# ---------------------------------------------------------------------------
# Tolerance for the M2_input * sin_beta_cos_beta == m12_sq_input check.
# We use a relative tolerance generous enough to allow for the long double ->
# double cast in the engine but strict enough to catch any gross derivation error.
# ---------------------------------------------------------------------------
REL_TOL = 1e-12
ABS_TOL = 1e-200   # fallback when the reference value is near zero


def validate(csv_path: str) -> bool:
    ok = True

    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []

        # --- Header checks ---
        if "M2_input" not in fieldnames:
            print("[FAIL] Header missing column: M2_input")
            ok = False

        if "m12_sq_input" not in fieldnames:
            print("[FAIL] Header missing column: m12_sq_input")
            ok = False

        if "computed_lam1" in fieldnames:
            print("[FAIL] Header contains forbidden column: computed_lam1")
            ok = False

        required_cols = {
            "M2_input", "m12_sq_input", "sin_beta_cos_beta",
            "m12_sq_out", "lam1_out", "lam2_out",
        }
        missing = required_cols - set(fieldnames)
        if missing:
            print(f"[FAIL] Header missing required columns: {sorted(missing)}")
            ok = False

        # Bail early if we cannot read data columns.
        if not ok:
            return False

        # --- Row checks ---
        n_rows = 0
        n_m12_mismatch = 0
        n_m12_out_bad = 0
        n_lam1_bad = 0
        n_lam2_bad = 0

        for row in reader:
            n_rows += 1

            try:
                M2_input          = float(row["M2_input"])
                m12_sq_input      = float(row["m12_sq_input"])
                sin_beta_cos_beta = float(row["sin_beta_cos_beta"])
                m12_sq_out        = float(row["m12_sq_out"])
                lam1_out          = float(row["lam1_out"])
                lam2_out          = float(row["lam2_out"])
            except (KeyError, ValueError) as exc:
                print(f"[FAIL] Could not parse row {n_rows}: {exc}")
                ok = False
                continue

            # m12_sq_input == M2_input * sin_beta_cos_beta
            expected = M2_input * sin_beta_cos_beta
            rel_err = (
                abs(m12_sq_input - expected) / max(abs(expected), ABS_TOL)
                if expected != 0.0 else abs(m12_sq_input - expected)
            )
            if rel_err > REL_TOL:
                n_m12_mismatch += 1
                if n_m12_mismatch <= 3:
                    print(
                        f"[FAIL] Row {n_rows}: m12_sq_input mismatch. "
                        f"got={m12_sq_input:.17e}, expected={expected:.17e}, "
                        f"rel_err={rel_err:.3e}"
                    )

            # Finiteness
            if not math.isfinite(m12_sq_out):
                n_m12_out_bad += 1
            if not math.isfinite(lam1_out):
                n_lam1_bad += 1
            if not math.isfinite(lam2_out):
                n_lam2_bad += 1

        # Summarise row-level failures
        if n_m12_mismatch > 0:
            print(f"[FAIL] m12_sq_input mismatch in {n_m12_mismatch}/{n_rows} rows.")
            ok = False
        if n_m12_out_bad > 0:
            print(f"[FAIL] m12_sq_out not finite in {n_m12_out_bad}/{n_rows} rows.")
            ok = False
        if n_lam1_bad > 0:
            print(f"[FAIL] lam1_out not finite in {n_lam1_bad}/{n_rows} rows.")
            ok = False
        if n_lam2_bad > 0:
            print(f"[FAIL] lam2_out not finite in {n_lam2_bad}/{n_rows} rows.")
            ok = False

        if n_rows == 0:
            print("[INFO] CSV has no data rows (empty output - physics may have filtered all).")
        else:
            print(f"[INFO] Checked {n_rows} data rows.")

    if ok:
        print("[PASS] All contract checks passed.")
    return ok


def main() -> None:
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <csv_file>")
        sys.exit(1)

    csv_path = sys.argv[1]
    passed = validate(csv_path)
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
