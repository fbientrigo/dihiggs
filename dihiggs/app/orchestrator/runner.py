"""
runner.py — ScanRunner: the main scan execution loop.

ScanRunner orchestrates a list of tan(beta) values through a chosen engine,
handling:
  - dry-run (print commands, write manifests, no subprocess)
  - resume by exact grid signature (engine-aware)
  - cross-run CSV reuse (resume_scope=fixed + materialize)
  - force overwrite
  - OMP_NUM_THREADS injection
  - stdout/stderr capture and persistence on failure
  - per-task scan_meta.json
  - task_summary.jsonl append
  - run_manifest.json (initial write + final summary update)
  - Total Attempts and TRIPLE_OK_POINTS parsing
"""

from __future__ import annotations

import os
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from dihiggs.app.orchestrator.engines.base import EngineAdapter
from dihiggs.app.orchestrator.grid import ScanGrid, grid_signature
from dihiggs.app.orchestrator.io_utils import (
    append_jsonl,
    log_msg,
    utc_now_iso,
)
from dihiggs.app.orchestrator.layout import (
    build_run_dir,
    default_run_name,
    format_tanbeta_folder,
)
from dihiggs.app.orchestrator.manifest import (
    update_manifest_summary,
    write_initial_manifest,
    write_scan_meta,
)
from dihiggs.app.orchestrator.models import FixedParams, TaskResult, TaskSpec
from dihiggs.app.orchestrator.parsing import parse_cpp_stats
from dihiggs.app.orchestrator.resume import (
    find_previous_csv,
    load_done_signatures,
    materialize_csv,
    should_skip,
)


@dataclass
class ScanRunner:
    """
    Orchestrates a series of C++ scan tasks with full lifecycle management.

    Parameters
    ----------
    engine:
        Engine adapter (Lambda1Engine or M2Engine).
    executable:
        Path to the compiled C++ binary.
    grid:
        Scan grid (mphi + second axis).
    fixed_base:
        Fixed physics parameters common to all tasks.
        ``fixed_base.tan_beta`` is overridden per task from ``tanbeta_list``.
    tanbeta_list:
        List of tan(beta) values to iterate over.
    campaign:
        Campaign label (used in directory naming).
    outdir:
        Root output directory.
    lake_name:
        Data-lake subdirectory name.
    run_name:
        Optional run identifier; auto-generated if None.
    dry_run:
        If True, print commands and write manifests but do not execute C++.
    force:
        If True, overwrite existing CSV outputs.
    resume_scope:
        ``"run"`` — only reuse CSVs from the current run.
        ``"fixed"`` — also reuse CSVs from previous runs under the same
        fixed-param directory (matching grid signature).
    materialize:
        If True and reusing a previous-run CSV, copy it into the current
        run directory.
    omp_threads:
        Value for OMP_NUM_THREADS; None means inherit from environment.
    timeout_s:
        Per-task subprocess timeout in seconds; None means no limit.
    """

    engine: EngineAdapter
    executable: Path
    grid: ScanGrid
    fixed_base: FixedParams
    tanbeta_list: List[float]
    campaign: str
    outdir: Path
    lake_name: str = "dihiggs_lake"
    run_name: Optional[str] = None
    dry_run: bool = False
    force: bool = False
    resume_scope: str = "fixed"  # "run" | "fixed"
    materialize: bool = False
    omp_threads: Optional[int] = None
    timeout_s: Optional[float] = None

    # Internal state (set during run())
    _run_dir: Optional[Path] = field(default=None, repr=False, init=False)
    _log_path: Optional[Path] = field(default=None, repr=False, init=False)
    _jsonl_path: Optional[Path] = field(default=None, repr=False, init=False)
    _manifest_path: Optional[Path] = field(default=None, repr=False, init=False)

    # ---------------------------------------------------------------------------
    # Public entry point
    # ---------------------------------------------------------------------------

    def run(self) -> List[TaskResult]:
        """
        Execute the full scan campaign and return a list of TaskResults.

        Returns
        -------
        List[TaskResult]
            One record per tan(beta) value processed.
        """
        # ---- resolve run name and directories --------------------------------
        run_name = self.run_name or default_run_name()
        # Build a temporary fixed with the first tb for run_dir naming;
        # the fixed_dir name uses only non-tb fixed params (mA, sin_ba,
        # lambda6, lambda7), so any tb value is fine here.
        _tmp_fixed = FixedParams(
            mA=self.fixed_base.mA,
            sin_ba=self.fixed_base.sin_ba,
            tan_beta=self.tanbeta_list[0] if self.tanbeta_list else 0.0,
            lambda6=self.fixed_base.lambda6,
            lambda7=self.fixed_base.lambda7,
            lambda1=self.fixed_base.lambda1,
        )
        run_dir = build_run_dir(
            outdir=self.outdir,
            lake_name=self.lake_name,
            campaign=self.campaign,
            fixed=_tmp_fixed,
            run_name=run_name,
        )
        run_dir.mkdir(parents=True, exist_ok=True)

        self._run_dir = run_dir
        self._log_path = run_dir / "orchestrator.log"
        self._jsonl_path = run_dir / "task_summary.jsonl"
        self._manifest_path = run_dir / "run_manifest.json"

        log = self._log_path

        self._log("[INIT] Starting scan orchestration.")
        self._log(f"[INIT] engine={self.engine.engine_name}  "
                  f"exec={self.executable}  dry_run={self.dry_run}  "
                  f"force={self.force}  resume_scope={self.resume_scope}")
        self._log(f"[INIT] run_dir={run_dir}")

        # ---- manifest (initial write) ----------------------------------------
        axis_meta = self.engine.axis_metadata()
        # grid_sig uses a representative fixed (tb doesn't affect the sig for
        # the run-level manifest; per-task sigs include the actual tb).
        manifest = write_initial_manifest(
            manifest_path=self._manifest_path,
            run_dir=run_dir,
            log_path=self._log_path,
            jsonl_path=self._jsonl_path,
            campaign=self.campaign,
            run_name=run_name,
            engine_name=self.engine.engine_name,
            axis_metadata=axis_meta,
            executable=self.executable,
            grid=self.grid,
            fixed=_tmp_fixed,
            tanbeta_list=self.tanbeta_list,
            omp_threads=self.omp_threads,
            dry_run=self.dry_run,
            force=self.force,
            resume_scope=self.resume_scope,
            materialize=self.materialize,
            outdir=self.outdir,
            lake_name=self.lake_name,
        )
        self._log(f"[MANI] Wrote manifest: {self._manifest_path}")

        # ---- load existing done signatures for resume ------------------------
        done_sigs: Set[str] = load_done_signatures(self._jsonl_path)

        # ---- subprocess environment ------------------------------------------
        env = os.environ.copy()
        if self.omp_threads is not None:
            env["OMP_NUM_THREADS"] = str(self.omp_threads)

        # ---- scan loop -------------------------------------------------------
        results: List[TaskResult] = []
        fixed_dir = run_dir.parent  # campaign/fixed_.../
        grand_t0 = time.time()
        total_tasks = len(self.tanbeta_list)
        tasks_ok = 0
        tasks_skipped = 0
        tasks_failed = 0
        total_attempts_sum = 0

        for idx, tb in enumerate(self.tanbeta_list, start=1):
            folder_name, tb_tag = format_tanbeta_folder(tb)
            tb_dir = run_dir / folder_name
            tb_dir.mkdir(parents=True, exist_ok=True)

            output_csv = tb_dir / f"scan_tb_{tb_tag}.csv"
            meta_path = tb_dir / "scan_meta.json"

            # Build per-task fixed (inject tb)
            task_fixed = FixedParams(
                mA=self.fixed_base.mA,
                sin_ba=self.fixed_base.sin_ba,
                tan_beta=tb,
                lambda6=self.fixed_base.lambda6,
                lambda7=self.fixed_base.lambda7,
                lambda1=self.fixed_base.lambda1,
            )
            sig = grid_signature(self.engine.engine_name, self.grid, task_fixed)

            self._log(
                f"[{idx:03d}/{total_tasks:03d}] tb={tb}  sig={sig}  "
                f"-> {output_csv.name}"
            )

            # ---- resume check (same-run) ------------------------------------
            skip, reason = should_skip(
                output_csv=output_csv,
                grid_sig=sig,
                force=self.force,
                done_signatures=done_sigs,
            )
            if skip:
                tasks_skipped += 1
                self._log(f"[SKIP] tb={tb}  reason={reason}")
                result = TaskResult(
                    task_id=f"tb_{tb_tag}",
                    status="skip",
                    grid_sig=sig,
                    output_csv=output_csv,
                    extra={"reason": reason, "tanbeta": tb},
                )
                append_jsonl(self._jsonl_path, self._task_summary_record(
                    result, idx, total_tasks, tb, tb_tag,
                ))
                results.append(result)
                continue

            # ---- resume check (cross-run, resume_scope=fixed) ---------------
            if not self.force and self.resume_scope == "fixed":
                prev_csv = find_previous_csv(
                    fixed_dir=fixed_dir,
                    current_run_dir=run_dir,
                    tb_tag=tb_tag,
                    desired_grid_sig=sig,
                )
                if prev_csv is not None:
                    tasks_skipped += 1
                    self._log(
                        f"[REUSE] tb={tb}  found_in_previous_run={prev_csv}"
                    )
                    extra: Dict[str, Any] = {
                        "reason": "found_in_previous_run",
                        "previous_csv": str(prev_csv),
                        "tanbeta": tb,
                    }
                    if self.materialize:
                        ok = materialize_csv(prev_csv, output_csv)
                        extra["materialized"] = ok
                        if ok:
                            self._log(f"[LINK] materialized -> {output_csv}")
                        else:
                            self._log("[WARN] materialize failed; continuing.")
                    result = TaskResult(
                        task_id=f"tb_{tb_tag}",
                        status="skip",
                        grid_sig=sig,
                        output_csv=output_csv,
                        extra=extra,
                    )
                    append_jsonl(self._jsonl_path, self._task_summary_record(
                        result, idx, total_tasks, tb, tb_tag,
                    ))
                    results.append(result)
                    continue

            # ---- build command -----------------------------------------------
            cmd = self.engine.build_command(
                executable=self.executable,
                grid=self.grid,
                fixed=task_fixed,
                output_csv=output_csv,
            )

            # ---- dry-run: write manifests, skip subprocess ------------------
            if self.dry_run:
                self._log(f"[DRY ] cmd = {' '.join(cmd)}")
                write_scan_meta(
                    meta_path,
                    event="dry_run",
                    task_id=f"tb_{tb_tag}",
                    grid_sig=sig,
                    grid=self.grid,
                    fixed=task_fixed,
                    axis_metadata=axis_meta,
                    engine_name=self.engine.engine_name,
                    cmd=cmd,
                    output_csv=output_csv,
                )
                result = TaskResult(
                    task_id=f"tb_{tb_tag}",
                    status="dry_run",
                    grid_sig=sig,
                    output_csv=output_csv,
                    extra={"cmd": cmd, "tanbeta": tb},
                )
                append_jsonl(self._jsonl_path, self._task_summary_record(
                    result, idx, total_tasks, tb, tb_tag,
                ))
                results.append(result)
                continue

            # ---- execute subprocess -----------------------------------------
            t0 = time.time()
            result = self._execute(
                cmd=cmd,
                env=env,
                tb=tb,
                tb_tag=tb_tag,
                sig=sig,
                output_csv=output_csv,
                meta_path=meta_path,
                tb_dir=tb_dir,
                task_fixed=task_fixed,
                axis_meta=axis_meta,
            )
            result.elapsed_s = time.time() - t0

            if result.status == "done":
                tasks_ok += 1
                if result.total_attempts is not None:
                    total_attempts_sum += result.total_attempts
                done_sigs.add(sig)
            else:
                tasks_failed += 1

            self._log(
                f"[{result.status.upper()}] tb={tb}  "
                f"elapsed={result.elapsed_s:.2f}s  "
                f"attempts={result.total_attempts}  "
                f"tripleOK={result.triple_ok_points}"
            )
            append_jsonl(self._jsonl_path, self._task_summary_record(
                result, idx, total_tasks, tb, tb_tag,
            ))
            results.append(result)

        # ---- final summary ---------------------------------------------------
        wall = time.time() - grand_t0
        self._log("-" * 60)
        self._log(
            f"[FIN ] ok={tasks_ok}  skipped={tasks_skipped}  "
            f"failed={tasks_failed}  wall={wall:.1f}s"
        )
        if total_attempts_sum > 0 and wall > 0:
            pph = int(total_attempts_sum / (wall / 3600.0))
            self._log(f"[PERF] throughput={pph} attempts/hour")

        update_manifest_summary(
            manifest_path=self._manifest_path,
            manifest=manifest,
            tasks_total=total_tasks,
            tasks_ok=tasks_ok,
            tasks_skipped=tasks_skipped,
            tasks_failed=tasks_failed,
            total_attempts=total_attempts_sum,
            wall_time_sec=wall,
        )
        self._log(f"[MANI] Updated manifest summary: {self._manifest_path}")

        return results

    # ---------------------------------------------------------------------------
    # Private helpers
    # ---------------------------------------------------------------------------

    def _execute(
        self,
        *,
        cmd: List[str],
        env: Dict[str, str],
        tb: float,
        tb_tag: str,
        sig: str,
        output_csv: Path,
        meta_path: Path,
        tb_dir: Path,
        task_fixed: FixedParams,
        axis_meta: Dict[str, Any],
    ) -> TaskResult:
        """Run the subprocess and return a TaskResult."""
        try:
            proc = subprocess.run(
                cmd,
                env=env,
                capture_output=True,
                text=True,
                timeout=self.timeout_s,
            )
            attempts, triple_ok = parse_cpp_stats(proc.stdout)

            if proc.returncode == 0:
                write_scan_meta(
                    meta_path,
                    event="done",
                    task_id=f"tb_{tb_tag}",
                    grid_sig=sig,
                    grid=self.grid,
                    fixed=task_fixed,
                    axis_metadata=axis_meta,
                    engine_name=self.engine.engine_name,
                    cmd=cmd,
                    output_csv=output_csv,
                    total_attempts=attempts,
                    triple_ok_points=triple_ok,
                    returncode=proc.returncode,
                )
                return TaskResult(
                    task_id=f"tb_{tb_tag}",
                    status="done",
                    grid_sig=sig,
                    output_csv=output_csv,
                    total_attempts=attempts,
                    triple_ok_points=triple_ok,
                    returncode=proc.returncode,
                    extra={"tanbeta": tb},
                )
            else:
                # Save stdout/stderr for debugging
                stdout_path = tb_dir / "stdout.log"
                stderr_path = tb_dir / "stderr.log"
                stdout_path.write_text(proc.stdout or "", encoding="utf-8")
                stderr_path.write_text(proc.stderr or "", encoding="utf-8")

                write_scan_meta(
                    meta_path,
                    event="fail",
                    task_id=f"tb_{tb_tag}",
                    grid_sig=sig,
                    grid=self.grid,
                    fixed=task_fixed,
                    axis_metadata=axis_meta,
                    engine_name=self.engine.engine_name,
                    cmd=cmd,
                    output_csv=output_csv,
                    returncode=proc.returncode,
                    stdout_path=stdout_path,
                    stderr_path=stderr_path,
                )
                return TaskResult(
                    task_id=f"tb_{tb_tag}",
                    status="fail",
                    grid_sig=sig,
                    output_csv=output_csv,
                    returncode=proc.returncode,
                    stdout_path=stdout_path,
                    stderr_path=stderr_path,
                    extra={"tanbeta": tb},
                )

        except Exception as exc:
            exc_repr = repr(exc)
            write_scan_meta(
                meta_path,
                event="crash",
                task_id=f"tb_{tb_tag}",
                grid_sig=sig,
                grid=self.grid,
                fixed=task_fixed,
                axis_metadata=axis_meta,
                engine_name=self.engine.engine_name,
                cmd=cmd,
                output_csv=output_csv,
                exception=exc_repr,
            )
            return TaskResult(
                task_id=f"tb_{tb_tag}",
                status="crash",
                grid_sig=sig,
                output_csv=output_csv,
                error_message=exc_repr,
                extra={"tanbeta": tb, "cmd": cmd},
            )

    def _task_summary_record(
        self,
        result: TaskResult,
        idx: int,
        total: int,
        tb: float,
        tb_tag: str,
    ) -> Dict[str, Any]:
        """Build the JSONL record for task_summary.jsonl."""
        rec = result.as_dict()
        rec.update({
            "utc": utc_now_iso(),
            "index": idx,
            "total": total,
            "tanbeta": tb,
            "tb_tag": tb_tag,
            "engine_name": self.engine.engine_name,
        })
        return rec

    def _log(self, msg: str) -> None:
        """Log to stdout and orchestrator.log."""
        log_msg(msg, self._log_path)
