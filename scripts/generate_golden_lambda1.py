#!/usr/bin/env python3
"""Regenerate the lambda1 golden characterization oracle.

Regeneration is a deliberate, reviewed act -- never a way to turn a red test
green. The generator runs every case TWICE and refuses to freeze output that is
not byte-identical across the two runs, so nondeterminism can never be baked
into the oracle.

Outputs (all under tests/golden/lambda1_v1/):
  expected.csv           one row per emitted CSV row, prefixed with case_id
  expected_markers.json  per-case stdout/stderr observations (cardinality, warnings)
  manifest.json          commit/dirty state, toolchain, dependency identity,
                         checksums, row counts, generation command

Usage:
  scripts/build_lambda1_characterization.sh
  python3 scripts/generate_golden_lambda1.py
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
GOLDEN_DIR = REPO_ROOT / "tests" / "golden" / "lambda1_v1"
CASES_PATH = GOLDEN_DIR / "cases.json"
BINARY = Path(
    os.environ.get("LAMBDA1_BINARY", REPO_ROOT / "build" / "lambda1_char" / "PhysLam1Scan")
)

# Markers emitted on stdout by PhysLam1Scan.
RE_ATTEMPTS = re.compile(r"^Total Attempts:\s*(\d+)\s*$", re.M)
RE_ROWS = re.compile(r"^Total CSV Rows:\s*(\d+)\s*$", re.M)
RE_TRIPLE = re.compile(r"^TRIPLE_OK_POINTS\s+(\d+)\s*$", re.M)
# Warnings emitted on stderr (by THDM, not by PhysLam1Scan).
RE_WARN_CONSTRUCT = re.compile(r"WARNING: Cannot set physical masses")
RE_WARN_ROUNDTRIP = re.compile(r"WARNING: set_param_phys_lam1 lambda_1 round-trip abs error")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_case(case: dict, out_csv: Path) -> tuple[str, dict]:
    """Run one case; return (raw csv text, observed markers)."""
    args = [str(BINARY)] + list(case["args"]) + [str(out_csv)]
    env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run(args, capture_output=True, text=True, env=env, cwd=REPO_ROOT)
    if proc.returncode != 0:
        raise SystemExit(
            f"case {case['case_id']}: binary exited {proc.returncode}\n{proc.stderr}"
        )
    csv_text = out_csv.read_text()
    csv_rows = max(len(csv_text.strip().splitlines()) - 1, 0)

    def _one(rx, hay, what):
        m = rx.search(hay)
        if not m:
            raise SystemExit(f"case {case['case_id']}: could not find {what} in stdout")
        return int(m.group(1))

    markers = {
        "total_attempts": _one(RE_ATTEMPTS, proc.stdout, "Total Attempts"),
        "total_csv_rows_reported": _one(RE_ROWS, proc.stdout, "Total CSV Rows"),
        "triple_ok_points": _one(RE_TRIPLE, proc.stdout, "TRIPLE_OK_POINTS"),
        "csv_rows_actual": csv_rows,
        "construction_warnings": len(RE_WARN_CONSTRUCT.findall(proc.stderr)),
        "roundtrip_warnings": len(RE_WARN_ROUNDTRIP.findall(proc.stderr)),
    }
    return csv_text, markers


def git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], capture_output=True, text=True, cwd=REPO_ROOT
    ).stdout.strip()


def main() -> int:
    if not BINARY.exists():
        raise SystemExit(
            f"binary not found: {BINARY}\nRun scripts/build_lambda1_characterization.sh first."
        )
    spec = json.loads(CASES_PATH.read_text())
    cases = spec["cases"]

    header: str | None = None
    out_lines: list[str] = []
    markers: dict[str, dict] = {}

    with tempfile.TemporaryDirectory() as td:
        for case in cases:
            cid = case["case_id"]
            # Determinism gate: run twice, require byte-identical results.
            text_a, mark_a = run_case(case, Path(td) / f"{cid}_a.csv")
            text_b, mark_b = run_case(case, Path(td) / f"{cid}_b.csv")
            if text_a != text_b:
                raise SystemExit(f"case {cid}: CSV not byte-identical across two runs")
            if mark_a != mark_b:
                raise SystemExit(f"case {cid}: markers differ across two runs")

            lines = text_a.strip().splitlines()
            if not lines:
                raise SystemExit(f"case {cid}: empty CSV (not even a header)")
            if header is None:
                header = lines[0]
            elif lines[0] != header:
                raise SystemExit(f"case {cid}: CSV header drift:\n{header}\n{lines[0]}")
            for row in lines[1:]:
                out_lines.append(f"{cid},{row}")
            markers[cid] = mark_a

    assert header is not None
    expected_csv = GOLDEN_DIR / "expected.csv"
    expected_csv.write_text("case_id," + header + "\n" + "\n".join(out_lines) + "\n")

    markers_path = GOLDEN_DIR / "expected_markers.json"
    markers_path.write_text(json.dumps(markers, indent=2, sort_keys=True) + "\n")

    gxx = subprocess.run(["g++", "--version"], capture_output=True, text=True).stdout
    status_short = git("status", "--short", "--untracked-files=no")
    manifest = {
        "schema_version": "lambda1_golden_v1",
        "purpose": (
            "Behavioral oracle for the lambda1-target construction path. Freezes CURRENT "
            "behavior; does not certify it as physically correct."
        ),
        "source": {
            "repo": "main_dihiggs",
            "commit": git("rev-parse", "HEAD"),
            "describe": git("describe", "--always", "--dirty"),
            "tracked_files_modified": bool(status_short),
            "git_status_short": status_short.splitlines(),
        },
        "binary": {
            "path": "build/lambda1_char/PhysLam1Scan",
            "built_by": "scripts/build_lambda1_characterization.sh",
            "link_mode": (
                "standalone against 2HDMC only; bypasses libdihiggs.a and does not link "
                "HiggsTools (PhysLam1Scan.cpp does not use it). See the script header for why "
                "the production make target cannot build from a clean checkout."
            ),
            "cxxflags": "-fopenmp -std=gnu++17 -Wall -Wextra -Wpedantic -O0 -g",
            "sources": ["dihiggs/src/PhysLam1Scan.cpp", "dihiggs/src/ParamUtils.cpp"],
        },
        "toolchain": {
            "compiler": gxx.splitlines()[0] if gxx else "unknown",
            "omp_num_threads": "1",
        },
        "dependency": {
            "2hdmc": "vendored fork at 2hdmc/ (adds THDM::set_param_phys_lam1; THDM::EPS=1e-9 vs stock 1e-12)",
            "2hdmc_tree_hash": git("rev-parse", "HEAD:2hdmc"),
            "2hdmc_src_tree_hash": git("rev-parse", "HEAD:2hdmc/src"),
            "physlam1scan_blob": git("rev-parse", "HEAD:dihiggs/src/PhysLam1Scan.cpp"),
            "paramutils_blob": git("rev-parse", "HEAD:dihiggs/src/ParamUtils.cpp"),
            "thdm_cpp_blob": git("rev-parse", "HEAD:2hdmc/src/THDM.cpp"),
            "thdm_h_blob": git("rev-parse", "HEAD:2hdmc/src/THDM.h"),
        },
        "compiled_constants": {
            "M_H": "125.0 (constexpr in dihiggs/src/PhysLam1Scan.cpp; NOT the paper's 125.09)",
            "THDM::EPS": "1e-9 (constexpr in 2hdmc/src/THDM.h; stock 2HDMC uses 1e-12)",
            "yukawas_type": "1 (hard-coded)",
            "mHp": "set equal to mA (mA_fixed is passed for both)",
        },
        "generation": {
            "command": "python3 scripts/generate_golden_lambda1.py",
            "each_case_run_twice_byte_identical": True,
        },
        "counts": {
            "cases": len(cases),
            "expected_csv_rows": len(out_lines),
            "csv_columns": len(header.split(",")) + 1,
        },
        "checksums": {
            "cases.json": sha256(CASES_PATH),
            "expected.csv": sha256(expected_csv),
            "expected_markers.json": sha256(markers_path),
        },
    }
    (GOLDEN_DIR / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )

    print(f"wrote {expected_csv} ({len(out_lines)} rows from {len(cases)} cases)")
    print(f"wrote {markers_path}")
    print(f"wrote {GOLDEN_DIR / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
