"""
manifest.py — Write and update run_manifest.json.

run_manifest.json captures all configuration for a scan run so that:
- Results are fully reproducible.
- Downstream tools can re-hydrate the run configuration.
- Cross-run checkpoint matching can validate grid semantics.

The manifest is written BEFORE the scan loop (with summary={}) and
updated AFTER the loop completes with final task statistics.
"""

from __future__ import annotations

import sys
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from dihiggs.app.orchestrator.grid import ScanGrid, grid_signature
from dihiggs.app.orchestrator.io_utils import (
    detect_git_info,
    safe_write_json,
    utc_now_iso,
)
from dihiggs.app.orchestrator.layout import host_metadata
from dihiggs.app.orchestrator.models import FixedParams


def write_initial_manifest(
    *,
    manifest_path: Path,
    run_dir: Path,
    log_path: Path,
    jsonl_path: Path,
    campaign: str,
    run_name: str,
    engine_name: str,
    axis_metadata: Dict[str, Any],
    executable: Path,
    grid: ScanGrid,
    fixed: FixedParams,
    tanbeta_list: List[float],
    omp_threads: Optional[int],
    dry_run: bool,
    force: bool,
    resume_scope: str,
    materialize: bool,
    outdir: Path,
    lake_name: str,
) -> Dict[str, Any]:
    """
    Write the initial run_manifest.json to *manifest_path*.

    The manifest contains everything needed to reproduce or audit the run.
    The ``summary`` field is written as an empty dict and filled in later
    by :func:`update_manifest_summary`.

    Returns
    -------
    dict
        The manifest payload (useful for the caller to retain a reference).
    """
    git = detect_git_info(Path.cwd())
    sig = grid_signature(engine_name, grid, fixed)

    manifest: Dict[str, Any] = {
        "schema_version": "orchestrator.v2",
        "created_utc": utc_now_iso(),
        "campaign": campaign,
        "run_name": run_name,
        "engine_name": engine_name,
        "engine": engine_name,
        "paths": {
            "outdir": str(outdir),
            "lake_name": lake_name,
            "run_dir": str(run_dir),
            "log": str(log_path),
            "task_summary_jsonl": str(jsonl_path),
        },
        "host": host_metadata(),
        "git": git,
        "runtime": {
            "executable": str(executable),
            "engine_name": engine_name,
            "dry_run": dry_run,
            "force": force,
            "omp_num_threads": omp_threads,
            "resume_scope": resume_scope,
            "materialize": materialize,
        },
        # Axis semantics — authoritative; downstream must read this, not infer
        # units from CSV column names.
        "axis_metadata": axis_metadata,
        "fixed_params": fixed.as_dict(),
        "scan_grid": {
            "mphi_min": grid.mphi_min,
            "mphi_max": grid.mphi_max,
            "n_mphi": grid.n_mphi,
            "axis_min": grid.axis_min,
            "axis_max": grid.axis_max,
            "n_axis": grid.n_axis,
            "grid_signature": sig,
        },
        "tanbeta_list": tanbeta_list,
        "summary": {},  # filled by update_manifest_summary()
    }
    if engine_name == "m2":
        manifest.update({
            "point_schema_version": "dihiggs.point.v2",
            "mass_convention": {
                "mh_GeV": 125.13 if fixed.mh is None else fixed.mh,
                "source": "PDG 2026 Higgs listing",
                "source_url": "https://pdg.lbl.gov/encoder_listings/s126.pdf",
            },
            "twohdmc_provenance": {
                "api": "THDM::set_param_phys",
                "repository_commit": git.get("commit", "unknown"),
            },
            "acceptance_definitions": axis_metadata["acceptance"],
        })
    safe_write_json(manifest_path, manifest)
    return manifest


def update_manifest_summary(
    manifest_path: Path,
    manifest: Dict[str, Any],
    *,
    tasks_total: int,
    tasks_ok: int,
    tasks_skipped: int,
    tasks_failed: int,
    total_attempts: int,
    wall_time_sec: float,
) -> None:
    """
    Fill in the ``summary`` block and rewrite *manifest_path*.

    Called once at the end of the scan loop.
    """
    manifest["summary"] = {
        "finished_utc": utc_now_iso(),
        "wall_time_sec": wall_time_sec,
        "tasks_total": tasks_total,
        "tasks_ok": tasks_ok,
        "tasks_skipped": tasks_skipped,
        "tasks_failed": tasks_failed,
        "total_attempts_parsed": total_attempts,
        "completion_status": "complete" if tasks_failed == 0 else "failed",
    }
    if manifest["engine_name"] == "m2":
        outputs = []
        for meta_path in sorted(Path(manifest["paths"]["run_dir"]).rglob("scan_meta.json")):
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            if meta.get("event") == "done":
                outputs.append({key: meta[key] for key in (
                    "output_csv", "output_sha256", "output_row_count", "output_columns", "command"
                )})
        manifest["outputs"] = outputs
    safe_write_json(manifest_path, manifest)


def write_scan_meta(
    meta_path: Path,
    *,
    event: str,
    task_id: str,
    grid_sig: str,
    grid: ScanGrid,
    fixed: FixedParams,
    axis_metadata: Dict[str, Any],
    engine_name: str,
    cmd: List[str],
    output_csv: Path,
    total_attempts: Optional[int] = None,
    triple_ok_points: Optional[str] = None,
    returncode: Optional[int] = None,
    stdout_path: Optional[Path] = None,
    stderr_path: Optional[Path] = None,
    exception: Optional[str] = None,
) -> None:
    """
    Write per-task scan_meta.json.

    Only ``event="done"`` records are eligible for checkpoint reuse;
    the resume logic checks ``meta["event"] == "done"`` before accepting.

    Parameters
    ----------
    meta_path:
        Destination path (e.g. ``run_dir/tb_10000/scan_meta.json``).
    event:
        ``"done"``, ``"fail"``, or ``"crash"``.
    axis_metadata:
        Engine-specific axis semantics — included so that every scan_meta.json
        is self-describing regarding units and axis meaning.
    """
    payload: Dict[str, Any] = {
        "schema_version": "scan_meta.v2",
        "event": event,
        "utc": utc_now_iso(),
        "task_id": task_id,
        "engine_name": engine_name,
        "grid_signature": grid_sig,
        # Axis semantics — the authoritative source for downstream analysis
        "axis_metadata": axis_metadata,
        "scan_grid": {
            "mphi_min": grid.mphi_min,
            "mphi_max": grid.mphi_max,
            "n_mphi": grid.n_mphi,
            "axis_min": grid.axis_min,
            "axis_max": grid.axis_max,
            "n_axis": grid.n_axis,
            "grid_signature": grid_sig,
        },
        "fixed_params": fixed.as_dict(),
        "command": cmd,
        "output_csv": str(output_csv),
        "total_attempts": total_attempts,
        "triple_ok_points": triple_ok_points,
    }
    if returncode is not None:
        payload["returncode"] = returncode
    if stdout_path is not None:
        payload["stdout_path"] = str(stdout_path)
    if stderr_path is not None:
        payload["stderr_path"] = str(stderr_path)
    if exception is not None:
        payload["exception"] = exception
    if event == "done" and output_csv.is_file():
        digest = hashlib.sha256()
        with output_csv.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        payload["output_sha256"] = digest.hexdigest()
        with output_csv.open(newline="", encoding="utf-8") as handle:
            rows = csv.reader(handle)
            payload["output_columns"] = next(rows, [])
            payload["output_row_count"] = sum(1 for _ in rows)

    safe_write_json(meta_path, payload)
