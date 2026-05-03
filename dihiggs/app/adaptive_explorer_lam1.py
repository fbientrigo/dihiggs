#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


def parse_lam1_range(raw: str) -> tuple[float, float]:
    parts = raw.split(":")
    if len(parts) != 2:
        raise ValueError(f"--lam1-range must be 'min:max', got: {raw}")
    return float(parts[0]), float(parts[1])


def parse_tb_value(raw: str) -> float:
    if "," in raw:
        raw = raw.split(",", 1)[0]
    return float(raw)


def run_physlam1scan(
    phys_exec: str,
    out_csv: Path,
    tb: float,
    lam1_value: float,
) -> tuple[int, str, str]:
    cmd = [
        phys_exec,
        "130",
        "130",
        "1",
        f"{lam1_value}",
        f"{lam1_value}",
        "1",
        "300",
        "0.999",
        f"{tb}",
        "0.1",
        "0.0",
        str(out_csv),
    ]

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    return proc.returncode, proc.stdout, proc.stderr


def main() -> int:
    parser = argparse.ArgumentParser(description="Adaptive explorer wrapper for PhysLam1Scan")
    parser.add_argument("--checkpoint-root", required=True)
    parser.add_argument("--tb-values", required=True)
    parser.add_argument("--lam1-range", required=True)
    parser.add_argument("--n-bins", required=True, type=int)
    parser.add_argument("--n-proposals", required=True, type=int)
    parser.add_argument("--iter", required=True, type=int)
    args = parser.parse_args()

    lam1_min, lam1_max = parse_lam1_range(args.lam1_range)
    tb_value = parse_tb_value(args.tb_values)

    checkpoint_root = Path(args.checkpoint_root)
    checkpoint_root.mkdir(parents=True, exist_ok=True)
    iter_dir = checkpoint_root / f"iter_{args.iter:04d}"
    iter_dir.mkdir(parents=True, exist_ok=True)

    proposal_index = 0
    lam1_value = lam1_min if args.n_proposals <= 1 else (lam1_min + lam1_max) * 0.5
    run_dir = iter_dir / f"run_{args.iter:04d}_{proposal_index:03d}"
    run_dir.mkdir(parents=True, exist_ok=True)
    out_csv = run_dir / "results.csv"

    phys_exec = str(Path(__file__).resolve().parent / "PhysLam1Scan")
    rc, stdout, stderr = run_physlam1scan(phys_exec, out_csv, tb_value, lam1_value)

    state = {
        "iter": args.iter,
        "proposals": [
            {
                "run_dir": str(run_dir),
                "bin_index": proposal_index,
                "lam1_min": lam1_value,
                "lam1_max": lam1_value,
                "tb": tb_value,
                "status": "success" if rc == 0 else "failed",
            }
        ],
    }
    (iter_dir / "adaptive_state.json").write_text(
        json.dumps(state, indent=2) + "\n", encoding="utf-8"
    )

    if stdout:
        sys.stdout.write(stdout)
    if stderr:
        sys.stderr.write(stderr)

    return rc


if __name__ == "__main__":
    raise SystemExit(main())
