"""
cli.py — Backward-compatible CLI entrypoint for the modular scan orchestrator.

Invocation examples
-------------------

Lambda1 scan (fully backward-compatible with legacy orchestrate_scans.py):
    python -m dihiggs.app.orchestrator \\
        --exec ./PhysScanWithFixings \\
        --campaign scan_999 \\
        --sin-ba 0.995 --lambda6 0.10 --lambda7 0.00 --mA 300 \\
        --tanbeta 10000,15000,20000 \\
        --lam1-min 0 --lam1-max 12 --n-lam1 666

M2 boundary scan (new):
    python -m dihiggs.app.orchestrator \\
        --engine m2 \\
        --exec ./Phys_M2BoundaryScan \\
        --campaign m2_scan_001 \\
        --sin-ba 1.0 --lambda6 0.001 --lambda7 0.0 --mA 500 \\
        --lambda1 1.0 \\
        --tanbeta 50.0 \\
        --axis-min 0 --axis-max 500000 --n-axis 50

Dry-run (prints commands, writes manifests, no C++ execution):
    python -m dihiggs.app.orchestrator --dry-run ...

Force overwrite:
    python -m dihiggs.app.orchestrator --force ...
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from dihiggs.app.orchestrator.engines.lambda1 import Lambda1Engine
from dihiggs.app.orchestrator.engines.m2 import M2Engine
from dihiggs.app.orchestrator.grid import ScanGrid
from dihiggs.app.orchestrator.io_utils import parse_csv_floats
from dihiggs.app.orchestrator.models import FixedParams
from dihiggs.app.orchestrator.runner import ScanRunner

# ---------------------------------------------------------------------------
# Defaults (same as legacy orchestrate_scans.py for backward compat)
# ---------------------------------------------------------------------------
_DEFAULT_MPHI_MIN = 130.0
_DEFAULT_MPHI_MAX = 290.0
_DEFAULT_N_MPHI = 15

# Lambda1 defaults (preserved from legacy)
_DEFAULT_LAM1_MIN = 0.0
_DEFAULT_LAM1_MAX = 12.0
_DEFAULT_N_LAM1 = 666

# M2 defaults
_DEFAULT_M2_MIN = 0.0
_DEFAULT_M2_MAX = 500_000.0
_DEFAULT_N_M2 = 50

_DEFAULT_MA = 300.0
_DEFAULT_SIN_BA = 0.995
_DEFAULT_L6 = 0.1
_DEFAULT_L7 = 0.0
_DEFAULT_TANBETA = "10000,15000,20000"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m dihiggs.app.orchestrator",
        description=(
            "Modular 2HDM scan orchestrator supporting PhysScanWithFixings "
            "(lambda1 axis) and Phys_M2BoundaryScan (M2 axis)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Engine selection
    p.add_argument(
        "--engine",
        type=str,
        default="lambda1",
        choices=["lambda1", "m2", "m2_tracker"],
        help=(
            "Physics engine: 'lambda1' → PhysScanWithFixings (second axis = lambda_1); "
            "'m2' → Phys_M2BoundaryScan (second axis = M^2); "
            "'m2_tracker' → Phys_M2BandTracker (second axis = M^2 tracked dynamically). "
        ),
    )

    # Executable
    p.add_argument(
        "--exec",
        dest="exec_path",
        type=str,
        default=None,
        help=(
            "Path to the compiled C++ binary. "
            "Defaults to './PhysScanWithFixings' for lambda1, "
            "'./Phys_M2BoundaryScan' for m2."
        ),
    )

    # Runtime
    p.add_argument(
        "--threads",
        type=int,
        default=None,
        help="OMP_NUM_THREADS. If omitted, inherits from environment.",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing CSV outputs.",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Create folders, write manifests, print commands — "
            "but do NOT execute C++."
        ),
    )
    p.add_argument(
        "--timeout",
        type=float,
        default=None,
        help="Per-task subprocess timeout in seconds.",
    )

    # Output organisation
    p.add_argument(
        "--outdir",
        type=str,
        default="./output",
        help="Root output directory.",
    )
    p.add_argument(
        "--lake-name",
        type=str,
        default="dihiggs_lake",
        help="Data-lake subdirectory name.",
    )
    p.add_argument(
        "--campaign",
        type=str,
        default="scan",
        help="Campaign label (e.g. 'scan_999').",
    )
    p.add_argument(
        "--run-name",
        type=str,
        default=None,
        help="Run identifier. Auto-generated from timestamp+host if omitted.",
    )

    # Resume
    p.add_argument(
        "--resume-scope",
        type=str,
        default="fixed",
        choices=["run", "fixed"],
        help=(
            "'fixed': also reuse CSVs from previous runs under the same "
            "fixed-param directory (matching grid signature). "
            "'run': only check current run."
        ),
    )
    p.add_argument(
        "--materialize",
        action="store_true",
        help="When reusing a previous-run CSV, copy it into the current run dir.",
    )

    # Fixed physics parameters
    p.add_argument("--mA", type=float, default=_DEFAULT_MA,
                   help="Fixed mA (CP-odd Higgs mass, also mHp) in GeV.")
    p.add_argument("--sin-ba", type=float, default=_DEFAULT_SIN_BA,
                   help="Fixed sin(b-a).")
    p.add_argument("--lambda6", type=float, default=_DEFAULT_L6,
                   help="Fixed lambda_6.")
    p.add_argument("--lambda7", type=float, default=_DEFAULT_L7,
                   help="Fixed lambda_7.")
    p.add_argument(
        "--lambda1",
        type=float,
        default=None,
        help=(
            "Fixed lambda_1. "
            "For --engine lambda1 this must be omitted (lambda1 is the scan axis)."
        ),
    )

    # m_phi grid (shared by both engines)
    p.add_argument("--mphi-min", type=float, default=_DEFAULT_MPHI_MIN,
                   help="m_phi scan min (GeV).")
    p.add_argument("--mphi-max", type=float, default=_DEFAULT_MPHI_MAX,
                   help="m_phi scan max (GeV).")
    p.add_argument("--n-mphi", type=int, default=_DEFAULT_N_MPHI,
                   help="m_phi grid points.")

    # Generic second axis (engine interprets the semantics)
    p.add_argument(
        "--axis-min",
        type=float,
        default=None,
        help=(
            "Second axis scan min. "
            "For lambda1 engine: lambda_1 min (dimensionless). "
            "For m2 engine: M^2 min (GeV^2)."
        ),
    )
    p.add_argument(
        "--axis-max",
        type=float,
        default=None,
        help="Second axis scan max.",
    )
    p.add_argument(
        "--n-axis",
        type=int,
        default=None,
        help="Number of second-axis grid points.",
    )

    # Legacy lambda1 aliases (backward compatible)
    p.add_argument("--lam1-min", "--lam1_min", dest="lam1_min",
                   type=float, default=None,
                   help="[lambda1 engine] Alias for --axis-min.")
    p.add_argument("--lam1-max", "--lam1_max", dest="lam1_max",
                   type=float, default=None,
                   help="[lambda1 engine] Alias for --axis-max.")
    p.add_argument("--n-lam1", "--N-lam1", dest="n_lam1",
                   type=int, default=None,
                   help="[lambda1 engine] Alias for --n-axis.")

    # tan(beta) list
    p.add_argument(
        "--tanbeta",
        type=str,
        default=_DEFAULT_TANBETA,
        help="Comma-separated tan(beta) values, e.g. '10000,15000,20000'.",
    )

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    # ---- select engine -------------------------------------------------------
    if args.engine == "lambda1":
        engine = Lambda1Engine()
    elif args.engine == "m2":
        engine = M2Engine()
    elif args.engine == "m2_tracker":
        from dihiggs.app.orchestrator.engines.m2_tracker import M2TrackerEngine
        engine = M2TrackerEngine()
    else:
        print(f"[ERROR] Unknown engine: {args.engine}", file=sys.stderr)
        return 2

    # ---- resolve executable --------------------------------------------------
    if args.exec_path:
        exec_path = Path(args.exec_path).expanduser().resolve()
    else:
        exec_path = (
            Path(".") / engine.executable_basename
        ).resolve()

    if not args.dry_run and not exec_path.exists():
        print(f"[ERROR] Executable not found: {exec_path}", file=sys.stderr)
        return 2

    # ---- resolve second axis -------------------------------------------------
    # Priority: --axis-{min,max}, then legacy --lam1-{min,max}
    axis_min = args.axis_min
    axis_max = args.axis_max
    n_axis = args.n_axis

    if args.engine == "lambda1":
        axis_min = axis_min if axis_min is not None else (
            args.lam1_min if args.lam1_min is not None else _DEFAULT_LAM1_MIN
        )
        axis_max = axis_max if axis_max is not None else (
            args.lam1_max if args.lam1_max is not None else _DEFAULT_LAM1_MAX
        )
        n_axis = n_axis if n_axis is not None else (
            args.n_lam1 if args.n_lam1 is not None else _DEFAULT_N_LAM1
        )
        # lambda1 must NOT be in fixed for this engine
        if args.lambda1 is not None:
            print(
                "[ERROR] --lambda1 is not allowed for --engine lambda1 "
                "(lambda1 is the scan axis). Omit --lambda1.",
                file=sys.stderr,
            )
            return 2
        fixed_lambda1 = None

    else:  # m2
        axis_min = axis_min if axis_min is not None else _DEFAULT_M2_MIN
        axis_max = axis_max if axis_max is not None else _DEFAULT_M2_MAX
        n_axis = n_axis if n_axis is not None else _DEFAULT_N_M2
        if args.lambda1 is not None:
            fixed_lambda1 = args.lambda1
        else:
            fixed_lambda1 = None

    # ---- build grid and fixed ------------------------------------------------
    grid = ScanGrid(
        mphi_min=args.mphi_min,
        mphi_max=args.mphi_max,
        n_mphi=args.n_mphi,
        axis_min=axis_min,
        axis_max=axis_max,
        n_axis=n_axis,
    )
    tanbeta_list = parse_csv_floats(args.tanbeta)
    if not tanbeta_list:
        print("[ERROR] --tanbeta list is empty.", file=sys.stderr)
        return 2

    # fixed_base.tan_beta is overridden per-task inside ScanRunner; 0.0 here.
    fixed_base = FixedParams(
        mA=args.mA,
        sin_ba=args.sin_ba,
        tan_beta=0.0,
        lambda6=args.lambda6,
        lambda7=args.lambda7,
        lambda1=fixed_lambda1,
    )

    # ---- construct and run ---------------------------------------------------
    runner = ScanRunner(
        engine=engine,
        executable=exec_path,
        grid=grid,
        fixed_base=fixed_base,
        tanbeta_list=tanbeta_list,
        campaign=args.campaign,
        outdir=Path(args.outdir).expanduser().resolve(),
        lake_name=args.lake_name,
        run_name=args.run_name,
        dry_run=args.dry_run,
        force=args.force,
        resume_scope=args.resume_scope,
        materialize=args.materialize,
        omp_threads=args.threads,
        timeout_s=args.timeout,
    )

    runner.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
