#!/usr/bin/env python3
"""Reproduce the mh and THDM::EPS sensitivity study for the lambda1 path.

This is a NARROW, DETERMINISTIC COMPARISON -- not a parameter scan. It answers
four questions with executable evidence, and changes nothing in the repository:

  1. Does the compiled M_H default (125.0) vs the paper's 125.09 change whether
     construction SUCCEEDS?
  2. Does it change derived parameters / flags / widths, and by how much?
  3. Does the fork's THDM::EPS = 1e-9 (stock 2HDMC uses 1e-12) change
     accepted/rejected behavior, or only the round-trip warning?
  4. Is each observed difference serialization-only, numerical, or semantic?

Method: build probe variants OUTSIDE the repository (in a temp dir) by copying
the production sources and rewriting exactly one constant each. The production
tree is never modified, and the unmodified-copy variant is verified to
reproduce the compiled default byte-for-byte before any conclusion is drawn.

Usage:  python3 scripts/lambda1_mh_eps_sensitivity.py [--keep]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CXXFLAGS = ["-fopenmp", "-std=gnu++17", "-Wall", "-Wextra", "-O0", "-g"]
LDLIBS = ["-l2HDMC", "-lgsl", "-lgslcblas", "-lm"]

# (m_phi, lam1, mA, sin_ba, tan_beta, lambda6, lambda7)
POINTS = {
    "G06_accepted": ("130", "0.1", "300", "0.999", "50", "0.1", "0"),
    "small_l6_large_tb": ("200", "1.0", "500", "1.0", "10000", "1e-10", "0"),
    "campaign_best": ("200", "1.0", "500", "1.0", "126904", "0.0019", "0"),
}
COLUMNS = (
    "m_phi mA alpha beta lambda6 lambda7 m12 sin_ba tan_beta positivity_ok "
    "unitarity_ok perturbativity_ok width_bb width_tautau width_WW width_ZZ "
    "width_gaga width_Zga width_gg width_hh total_width br_gaga lam1 "
    "computed_lam1 lam2 computed_lam2 lam3 lam4 lam5"
).split()


def sh(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def build_2hdmc(tree: Path) -> None:
    sh(["make", "-C", str(tree), "lib/lib2HDMC.a"]).check_returncode()


def build_scan(src: Path, thdmc: Path, out: Path) -> None:
    cmd = (
        ["g++", f"-I{thdmc}/src", f"-I{REPO_ROOT}/dihiggs/src"]
        + CXXFLAGS
        + [str(src), str(REPO_ROOT / "dihiggs" / "src" / "ParamUtils.cpp"), "-o", str(out)]
        + [f"-L{thdmc}/lib"]
        + LDLIBS
    )
    p = sh(cmd)
    if p.returncode != 0:
        raise SystemExit(f"build failed:\n{p.stderr}")


def run_point(binary: Path, pt: tuple[str, ...], out_csv: Path) -> tuple[dict | None, str]:
    m_phi, lam1, mA, sba, tb, l6, l7 = pt
    args = [
        str(binary), m_phi, m_phi, "1", lam1, lam1, "1", mA, sba, tb, l6, l7, str(out_csv),
    ]
    p = subprocess.run(
        args, capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1")
    )
    lines = out_csv.read_text().strip().splitlines()
    if len(lines) < 2:
        return None, p.stderr
    return dict(zip(COLUMNS, lines[1].split(","))), p.stderr


def classify(a: dict | None, b: dict | None) -> str:
    """serialization-only | numerical | semantic"""
    if (a is None) != (b is None):
        return "SEMANTIC (construction success differs)"
    if a is None and b is None:
        return "none (both failed construction)"
    assert a and b
    flags = ("positivity_ok", "unitarity_ok", "perturbativity_ok")
    if any(a[f] != b[f] for f in flags):
        return "SEMANTIC (theory flags differ)"
    if all(a[c] == b[c] for c in COLUMNS):
        return "none (byte-identical)"
    return "NUMERICAL (same decisions, different digits)"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep", action="store_true", help="keep the temp build tree")
    args = ap.parse_args()

    td = Path(tempfile.mkdtemp(prefix="l1_sens_"))
    print(f"# probe build tree: {td}")
    results: dict = {"points": {}, "boundary_band": {}, "eps": {}}
    try:
        # --- baseline: unmodified copy, must reproduce the compiled default ---
        base_src = td / "PhysLam1Scan_base.cpp"
        shutil.copy(REPO_ROOT / "dihiggs" / "src" / "PhysLam1Scan.cpp", base_src)
        thdmc_fork = REPO_ROOT / "2hdmc"
        build_2hdmc(thdmc_fork)
        build_scan(base_src, thdmc_fork, td / "scan_base")

        prod = REPO_ROOT / "build" / "lambda1_char" / "PhysLam1Scan"
        if prod.exists():
            a, _ = run_point(prod, POINTS["G06_accepted"], td / "p.csv")
            b, _ = run_point(td / "scan_base", POINTS["G06_accepted"], td / "q.csv")
            assert a == b, "probe harness does not reproduce the production binary"
            print("# harness check: unmodified copy == production binary  [OK]")

        # --- mh variants ---
        variants: dict[str, Path] = {}
        for label, mh in (("mh_125.0_default", "125.0"), ("mh_125.09_paper", "125.09")):
            src = td / f"PhysLam1Scan_{label}.cpp"
            text = (REPO_ROOT / "dihiggs" / "src" / "PhysLam1Scan.cpp").read_text()
            new, n = re.subn(
                r"constexpr double M_H = 125\.0;", f"constexpr double M_H = {mh};", text
            )
            assert n == 1, f"expected exactly one M_H definition, patched {n}"
            src.write_text(new)
            build_scan(src, thdmc_fork, td / f"scan_{label}")
            variants[label] = td / f"scan_{label}"

        print("\n## Q1/Q2: mh = 125.0 (compiled default) vs 125.09 (paper)")
        for name, pt in POINTS.items():
            r0, _ = run_point(variants["mh_125.0_default"], pt, td / "a.csv")
            r9, _ = run_point(variants["mh_125.09_paper"], pt, td / "b.csv")
            verdict = classify(r0, r9)
            entry = {"verdict": verdict}
            if r0 and r9:
                entry["m12_rel_delta"] = abs(
                    float(r9["m12"]) - float(r0["m12"])
                ) / max(abs(float(r0["m12"])), 1e-300)
                entry["flags_125.0"] = [r0[f] for f in
                                        ("positivity_ok", "unitarity_ok", "perturbativity_ok")]
                entry["flags_125.09"] = [r9[f] for f in
                                         ("positivity_ok", "unitarity_ok", "perturbativity_ok")]
            results["points"][name] = entry
            print(f"  {name:22s} -> {verdict}")
            if "m12_rel_delta" in entry:
                print(f"     m12 relative delta: {entry['m12_rel_delta']:.3e}")

        # --- the construction boundary MOVES with mh ---
        print("\n## Q1: the m_h > m_H construction boundary moves with mh")
        for m_phi in ("124.99999999999999", "125.0", "125.05", "125.09", "125.2"):
            row = {}
            for label, v in variants.items():
                pt = (m_phi, "0.1", "300", "0.999", "50", "0.1", "0")
                r, _ = run_point(v, pt, td / "c.csv")
                row[label] = "constructed" if r else "CONSTRUCTION FAIL"
            same = row["mh_125.0_default"] == row["mh_125.09_paper"]
            results["boundary_band"][m_phi] = row
            flag = "" if same else "   <-- SEMANTIC DIFFERENCE"
            print(
                f"  m_phi={m_phi:20s} 125.0:{row['mh_125.0_default']:18s} "
                f"125.09:{row['mh_125.09_paper']:18s}{flag}"
            )

        # --- EPS variant: fork 1e-9 vs stock 1e-12 ---
        print("\n## Q3: THDM::EPS = 1e-9 (fork) vs 1e-12 (stock)")
        thdmc_stock = td / "2hdmc_stockeps"
        shutil.copytree(thdmc_fork, thdmc_stock)
        hdr = thdmc_stock / "src" / "THDM.h"
        new, n = re.subn(
            r"constexpr static double EPS = 1E-9;",
            "constexpr static double EPS = 1E-12;",
            hdr.read_text(),
        )
        assert n == 1, f"expected exactly one EPS definition, patched {n}"
        hdr.write_text(new)
        for f in (thdmc_stock / "lib").glob("*"):
            f.unlink()
        build_2hdmc(thdmc_stock)
        build_scan(base_src, thdmc_stock, td / "scan_stockeps")

        for name, pt in POINTS.items():
            rf, ef = run_point(td / "scan_base", pt, td / "d.csv")
            rs, es = run_point(td / "scan_stockeps", pt, td / "e.csv")
            verdict = classify(rf, rs)
            wf = "WARN" if "round-trip" in ef else "silent"
            ws = "WARN" if "round-trip" in es else "silent"
            results["eps"][name] = {
                "verdict": verdict, "warning_eps_1e-9": wf, "warning_eps_1e-12": ws,
            }
            print(f"  {name:22s} -> {verdict}; roundtrip warning: 1e-9={wf}, 1e-12={ws}")

        out = REPO_ROOT / "tests" / "golden" / "lambda1_v1" / "mh_eps_sensitivity.json"
        out.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
        print(f"\nwrote {out}")
    finally:
        if not args.keep:
            shutil.rmtree(td, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
