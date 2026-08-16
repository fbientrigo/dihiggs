"""
cli.py — CLI entrypoint for canonical and compatibility scan paths.

Invocation examples
-------------------

Canonical lambda1 v2 scan (row-preserving):
    python -m dihiggs.app.orchestrator \\
        --engine lambda1_legacy --exec ./PhysScanWithFixings \\
        --campaign scan_999 \\
        --sin-ba 0.995 --lambda6 0.10 --lambda7 0.00 --mA 300 \\
        --tanbeta 10000,15000,20000 \\
        --lam1-min 0 --lam1-max 12 --n-lam1 666

M2 boundary scan (new):
    python -m dihiggs.app.orchestrator \\
        --engine m2 \\
        --exec ./dihiggs/app/DihiggsPointV2Evaluator \\
        --campaign m2_scan_001 \\
        --sin-ba 1.0 --lambda6 0.001 --lambda7 0.0 --mA 500 \\
        --tanbeta 50.0 \\
        --axis-min 0 --axis-max 500000 --n-axis 50

Dry-run (prints commands, writes manifests, no C++ execution):
    python -m dihiggs.app.orchestrator --dry-run ...

Force overwrite:
    python -m dihiggs.app.orchestrator --force ...
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

from dihiggs.app.orchestrator import conventions
from dihiggs.app.orchestrator.engines.lambda1 import Lambda1Engine
from dihiggs.app.orchestrator.engines.m2 import M2Engine
from dihiggs.app.orchestrator.grid import ScanGrid
from dihiggs.app.orchestrator.io_utils import parse_csv_floats
from dihiggs.app.orchestrator.lambda1_v2 import (
    Lambda1Fixed,
    SCHEMA_VERSION,
    cartesian_rows,
    grid_values,
    run_lambda1_v2,
)
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

# CPU headroom: cores left free for the OS/sshd so the host never starves.
_CPU_HEADROOM = 2


def _safe_default_threads(reserve: int = _CPU_HEADROOM) -> int:
    """OMP thread count that always leaves CPU headroom for the OS/sshd.

    A scan that pins *every* core (OMP_NUM_THREADS == ncpu) starves sshd and
    the rest of the system and can freeze the whole box — this actually
    happened (jun 2026: a 12-thread run on a 12-core host stopped logging
    mid-run and had to be power-cycled). We therefore reserve a couple of
    cores by default so the machine stays interactively reachable no matter
    how long the job runs. Override with ``--threads N`` or ``--all-cores``.
    """
    ncpu = os.cpu_count() or 1
    return max(1, ncpu - reserve)


def resolve_omp_threads(
    all_cores: bool,
    threads: int | None,
    reserve: int = _CPU_HEADROOM,
) -> int | None:
    """Resolve the OMP_NUM_THREADS value with safe CPU-headroom defaults.

    Parameters mirror the CLI flags ``--all-cores`` and ``--threads``:

    - ``all_cores=True``  -> ``None`` (no cap; the engine may grab every core).
    - ``threads`` given   -> that explicit value (passed through unchanged).
    - otherwise           -> ``_safe_default_threads(reserve)`` so a couple of
      cores stay free for the OS/sshd.
    """
    if all_cores:
        return None
    if threads is not None:
        return threads
    return _safe_default_threads(reserve)
_DEFAULT_N_M2 = 50

_DEFAULT_MA = 300.0
_DEFAULT_SIN_BA = 0.995
_DEFAULT_L6 = 0.1
_DEFAULT_L7 = 0.0
_DEFAULT_TANBETA = "10000,15000,20000"


def _read_bronze_fixing(bronze_csv: str) -> dict:
    """
    Read (tan_beta, m_A, lambda6, lambda7, sin_ba) from the first data row of
    a chris/CalcLambda1ScanFixings bronze shard CSV.

    A bronze shard is one fixing tuple, so these columns are constant across
    all rows; the first row is authoritative.
    """
    with open(bronze_csv, newline="") as fh:
        reader = csv.DictReader(fh)
        row = next(reader, None)
    if row is None:
        raise ValueError(f"Bronze CSV has no data rows: {bronze_csv}")
    try:
        return {
            "tan_beta": float(row["tan_beta"]),
            "mA": float(row["m_A"]),
            "lambda6": float(row["lambda6"]),
            "lambda7": float(row["lambda7"]),
            "sin_ba": float(row["sin_ba"]),
        }
    except KeyError as exc:
        raise ValueError(
            f"Bronze CSV missing expected column {exc}: {bronze_csv}"
        ) from exc


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m dihiggs.app.orchestrator",
        description=(
            "Executable 2HDMC orchestrator: lambda1_v2 and m2 are canonical; "
            "lambda1_legacy and m2_tracker are compatibility/experimental paths."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Engine selection
    p.add_argument(
        "--engine",
        type=str,
        default="lambda1_v2",
        choices=["lambda1_v2", "lambda1_legacy", "lambda1", "m2", "m2_tracker", "gen_fixings"],
        help=(
            "Canonical 'lambda1_v2' → Lambda1EvaluatorV2 (explicit lambda1_target); "
            "legacy 'lambda1_legacy'/'lambda1' → PhysScanWithFixings; "
            "canonical 'm2' → DihiggsPointV2Evaluator (reconstructed lambda1); "
            "experimental 'm2_tracker' → bounded-pilot Phys_M2BandTracker; "
            "'gen_fixings' → GenScanWithFixings (Stage 2: calibrates/validates a "
            "chris/CalcLambda1ScanFixings bronze shard against 2HDMC). "
        ),
    )

    # Executable
    p.add_argument(
        "--exec",
        dest="exec_path",
        type=str,
        default=None,
        help=(
            "Path to the compiled C++ binary. Defaults to the selected executable "
            "under dihiggs/app for canonical v2 paths."
        ),
    )

    # Runtime
    p.add_argument(
        "--threads",
        type=int,
        default=None,
        help=(
            "OMP_NUM_THREADS for the C++ engine. Default: (logical CPUs - "
            f"{_CPU_HEADROOM}), leaving headroom so the host stays reachable. "
            "Pass an explicit value to override, or --all-cores to use every "
            "core (NOT recommended on a machine you also SSH into)."
        ),
    )
    p.add_argument(
        "--all-cores",
        action="store_true",
        help=(
            "Disable the automatic CPU-headroom cap and let the engine use "
            "every core (OMP_NUM_THREADS left to the environment). Can starve "
            "sshd and freeze the host; only use on a dedicated batch node."
        ),
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
    p.add_argument("--mh", type=float, default=conventions.M_H_GEV,
                   help="[m2] Explicit h1 mass in GeV. Defaults to the canonical "
                        "convention in conventions/physics_conventions.yaml "
                        "(%s GeV, %s)." % (conventions.M_H_GEV_TEXT, conventions.M_H_SOURCE))
    p.add_argument("--mHp", type=float, default=None,
                   help="[m2] Explicit charged-Higgs mass in GeV; defaults to mA.")
    p.add_argument("--yukawa-type", type=int, choices=(1, 2, 3, 4), default=1,
                   help="[m2] 2HDMC Yukawa type.")
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
            "Legacy fixed lambda_1 compatibility option. Rejected for canonical "
            "m2 because lambda1 is reconstructed output, not an input."
        ),
    )

    # m_phi grid (shared by both engines)
    p.add_argument("--mH-min", "--mphi-min", dest="mphi_min", type=float,
                   default=_DEFAULT_MPHI_MIN, help="mH scan minimum (GeV).")
    p.add_argument("--mH-max", "--mphi-max", dest="mphi_max", type=float,
                   default=_DEFAULT_MPHI_MAX, help="mH scan maximum (GeV).")
    p.add_argument("--n-mH", "--n-mphi", dest="n_mphi", type=int,
                   default=_DEFAULT_N_MPHI, help="Number of mH grid points.")

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

    # gen_fixings (Stage 2) options
    p.add_argument(
        "--bronze-csv",
        type=str,
        default=None,
        help=(
            "[gen_fixings engine] Path to a chris/CalcLambda1ScanFixings "
            "bronze shard CSV to calibrate/validate. Required for "
            "--engine gen_fixings."
        ),
    )
    p.add_argument(
        "--calibration-n",
        type=int,
        default=None,
        help=(
            "[gen_fixings engine] Number of calibration candidates (N) for "
            "the +/-frac random search. Omit to use the C++ binary's default (50)."
        ),
    )
    p.add_argument(
        "--calibration-frac",
        type=float,
        default=None,
        help=(
            "[gen_fixings engine] Fractional jitter half-width for the "
            "calibration random search (e.g. 0.10 = +/-10%%). Omit to use "
            "the C++ binary's default (0.10)."
        ),
    )
    p.add_argument(
        "--rng-seed",
        type=int,
        default=None,
        help=(
            "[gen_fixings engine] Base RNG seed for the calibration random "
            "search. Omit to use the C++ binary's default (0)."
        ),
    )

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.engine == "lambda1_v2":
        if args.lambda1 is not None:
            print(
                "[ERROR] --lambda1 is not accepted by lambda1_v2; "
                "lambda1_target is the grid axis.", file=sys.stderr,
            )
            return 2
        executable = (
            Path(args.exec_path).expanduser().resolve()
            if args.exec_path else Path("dihiggs/app/Lambda1EvaluatorV2").resolve()
        )
        if not args.dry_run and not executable.exists():
            print(f"[ERROR] Executable not found: {executable}", file=sys.stderr)
            return 2
        axis_min = args.axis_min if args.axis_min is not None else (
            args.lam1_min if args.lam1_min is not None else _DEFAULT_LAM1_MIN
        )
        axis_max = args.axis_max if args.axis_max is not None else (
            args.lam1_max if args.lam1_max is not None else _DEFAULT_LAM1_MAX
        )
        n_axis = args.n_axis if args.n_axis is not None else (
            args.n_lam1 if args.n_lam1 is not None else _DEFAULT_N_LAM1
        )
        fixed = Lambda1Fixed(
            args.mh, args.mA, args.mA if args.mHp is None else args.mHp,
            args.sin_ba, args.lambda6, args.lambda7,
        )
        rows = cartesian_rows(
            fixed=fixed,
            mH_values=grid_values(args.mphi_min, args.mphi_max, args.n_mphi),
            lambda1_values=grid_values(axis_min, axis_max, n_axis),
            tan_beta_values=parse_csv_floats(args.tanbeta),
        )
        try:
            manifest = run_lambda1_v2(
                executable=executable,
                rows=rows,
                outdir=Path(args.outdir).expanduser().resolve() / args.lake_name,
                campaign=args.campaign,
                run_name=args.run_name or "run",
                repo_root=Path.cwd(),
                dry_run=args.dry_run,
                force=args.force,
                timeout_s=args.timeout,
            )
        except (OSError, ValueError, RuntimeError) as error:
            print(f"[ERROR] {error}", file=sys.stderr)
            return 2
        print(f"[lambda1_v2] {manifest['status']} rows={len(rows)} schema={SCHEMA_VERSION}")
        return 0

    # ---- select engine -------------------------------------------------------
    if args.engine in ("lambda1", "lambda1_legacy"):
        engine = Lambda1Engine()
    elif args.engine == "m2":
        engine = M2Engine()
    elif args.engine == "m2_tracker":
        from dihiggs.app.orchestrator.engines.m2_tracker import M2TrackerEngine
        engine = M2TrackerEngine()
    elif args.engine == "gen_fixings":
        from dihiggs.app.orchestrator.engines.gen_fixings import GenFixingsEngine
        engine = GenFixingsEngine()
    else:
        print(f"[ERROR] Unknown engine: {args.engine}", file=sys.stderr)
        return 2

    # ---- resolve executable --------------------------------------------------
    if args.exec_path:
        exec_path = Path(args.exec_path).expanduser().resolve()
    else:
        exec_path = (Path("dihiggs/app") / engine.executable_basename).resolve()

    if not args.dry_run and not exec_path.exists():
        print(f"[ERROR] Executable not found: {exec_path}", file=sys.stderr)
        return 2

    # ---- resolve second axis -------------------------------------------------
    # Priority: --axis-{min,max}, then legacy --lam1-{min,max}
    axis_min = args.axis_min
    axis_max = args.axis_max
    n_axis = args.n_axis

    # mA/sin_ba/lambda6/lambda7/tan_beta: normally taken from CLI flags, but
    # for gen_fixings these are derived from the bronze shard itself (below)
    # to avoid silently writing mismatched provenance into the run_dir name.
    mA_val = args.mA
    sin_ba_val = args.sin_ba
    lambda6_val = args.lambda6
    lambda7_val = args.lambda7
    tanbeta_override = None

    if args.engine in ("lambda1", "lambda1_legacy"):
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
                f"[ERROR] --lambda1 is not allowed for --engine {args.engine} "
                "(lambda1 is the scan axis). Omit --lambda1.",
                file=sys.stderr,
            )
            return 2
        fixed_lambda1 = None

    elif args.engine == "gen_fixings":
        if not args.bronze_csv:
            print(
                "[ERROR] --bronze-csv is required for --engine gen_fixings.",
                file=sys.stderr,
            )
            return 2
        # No generated grid for this engine; ScanGrid is a required-but-unused
        # placeholder (see GenFixingsEngine.build_command docstring).
        axis_min = axis_min if axis_min is not None else 0.0
        axis_max = axis_max if axis_max is not None else 0.0
        n_axis = n_axis if n_axis is not None else 1
        fixed_lambda1 = args.lambda1

        # GenScanWithFixings reads (tan_beta, m_A, lambda6, lambda7, sin_ba)
        # per-row from the bronze shard itself, not from --mA/--sin-ba/etc.
        # Derive run-dir/grid-signature provenance from the shard so it can't
        # silently diverge from the data being processed.
        try:
            bronze_fixing = _read_bronze_fixing(args.bronze_csv)
        except (OSError, ValueError) as exc:
            print(f"[ERROR] {exc}", file=sys.stderr)
            return 2
        mA_val = bronze_fixing["mA"]
        sin_ba_val = bronze_fixing["sin_ba"]
        lambda6_val = bronze_fixing["lambda6"]
        lambda7_val = bronze_fixing["lambda7"]
        tanbeta_override = [bronze_fixing["tan_beta"]]
        print(
            f"[gen_fixings] Derived fixed params from bronze shard "
            f"{args.bronze_csv}: mA={mA_val}, sin_ba={sin_ba_val}, "
            f"lambda6={lambda6_val}, lambda7={lambda7_val}, "
            f"tan_beta={tanbeta_override[0]} "
            f"(--mA/--sin-ba/--lambda6/--lambda7/--tanbeta are ignored "
            f"for this engine)."
        )

    else:  # m2 or experimental m2_tracker
        axis_min = axis_min if axis_min is not None else _DEFAULT_M2_MIN
        axis_max = axis_max if axis_max is not None else _DEFAULT_M2_MAX
        n_axis = n_axis if n_axis is not None else _DEFAULT_N_M2
        if args.lambda1 is not None:
            label = "m2_tracker" if args.engine == "m2_tracker" else "m2"
            print(
                f"[ERROR] --lambda1 is rejected for --engine {label}: "
                "lambda1 is reconstructed output, not a fixed input.",
                file=sys.stderr,
            )
            return 2
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
    if tanbeta_override is not None:
        tanbeta_list = tanbeta_override
    else:
        tanbeta_list = parse_csv_floats(args.tanbeta)
    if not tanbeta_list:
        print("[ERROR] --tanbeta list is empty.", file=sys.stderr)
        return 2

    # fixed_base.tan_beta is overridden per-task inside ScanRunner; 0.0 here.
    fixed_base = FixedParams(
        mA=mA_val,
        sin_ba=sin_ba_val,
        tan_beta=0.0,
        lambda6=lambda6_val,
        lambda7=lambda7_val,
        lambda1=fixed_lambda1,
        bronze_shard_csv=args.bronze_csv if args.engine == "gen_fixings" else None,
        calibration_n=args.calibration_n if args.engine == "gen_fixings" else None,
        calibration_frac=args.calibration_frac if args.engine == "gen_fixings" else None,
        rng_seed=args.rng_seed if args.engine == "gen_fixings" else None,
        mh=args.mh if args.engine in ("m2", "m2_tracker") else None,
        mHp=(args.mA if args.mHp is None else args.mHp) if args.engine in ("m2", "m2_tracker") else None,
        yukawa_type=args.yukawa_type if args.engine == "m2" else None,
    )

    # ---- resolve OMP threads (safe-by-default CPU headroom) ------------------
    ncpu = os.cpu_count() or 1
    omp_threads = resolve_omp_threads(args.all_cores, args.threads)
    if args.all_cores:
        print(
            "[WARN] --all-cores: no CPU-headroom cap. This can starve sshd and "
            "freeze a host you also SSH into.",
            file=sys.stderr,
        )
    elif args.threads is not None:
        if omp_threads is not None and omp_threads >= ncpu:
            print(
                f"[WARN] --threads {omp_threads} >= {ncpu} logical CPUs: no "
                "headroom left for the OS/sshd; the host may become "
                "unreachable mid-run.",
                file=sys.stderr,
            )
    else:
        print(
            f"[INIT] OMP_NUM_THREADS={omp_threads} "
            f"(auto CPU headroom: {ncpu} logical CPUs - {_CPU_HEADROOM}). "
            "Override with --threads N or --all-cores."
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
        omp_threads=omp_threads,
        timeout_s=args.timeout,
    )

    runner.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
