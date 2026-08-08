"""
resume.py — Resume / skip logic and cross-run CSV reuse.

Resume semantics
----------------
A task is eligible to be skipped (not re-run) when:
  1. ``--force`` is not set, AND
  2. Either:
     a. The expected output CSV exists in the current run directory and its
        per-task ``scan_meta.json`` records ``event="done"`` with the exact
        grid signature, OR
     b. ``--resume-scope=fixed`` is set and a non-empty CSV from a previous run
        under the same fixed-param directory exists with a matching per-task
        ``scan_meta.json`` carrying ``event="done"`` and the exact grid
        signature.

A historical ``done`` record without its output artifact is not resumable.
Likewise, run-level manifests are not sufficient evidence that an individual
CSV completed successfully.

Grid signature matching
-----------------------
The grid signature (from grid.grid_signature()) encodes the engine name,
grid numeric ranges, and fixed params.  Two tasks with different engines
but identical numeric ranges will have different signatures, preventing
accidental reuse of lambda1 CSVs for M2 scans.

Cross-run matching only reuses scan_meta.json records with event="done".
Partial/failed runs are never reused.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Optional, Set, Tuple

from dihiggs.app.orchestrator.io_utils import load_json_best_effort


# ---------------------------------------------------------------------------
# Load completed signatures from an existing task_summary.jsonl
# ---------------------------------------------------------------------------

def load_done_signatures(jsonl_path: Path) -> Set[str]:
    """
    Read a task_summary.jsonl and return the set of grid signatures
    whose tasks completed with status ``"done"``.

    These records are useful audit information, but a signature alone is not
    sufficient to skip execution: the expected output artifact and its
    per-task ``scan_meta.json`` must still be present and valid.

    Parameters
    ----------
    jsonl_path:
        Path to the task_summary.jsonl file from a previous (or current) run.

    Returns
    -------
    Set[str]
        Grid signatures that have been recorded as successfully completed.
    """
    done: Set[str] = set()
    if not jsonl_path.exists():
        return done
    try:
        with jsonl_path.open("r", encoding="utf-8") as fh:
            import json
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except Exception:
                    continue
                if rec.get("status") == "done" or rec.get("event") == "done":
                    sig = rec.get("grid_signature")
                    if sig:
                        done.add(sig)
    except Exception:
        pass
    return done


# ---------------------------------------------------------------------------
# Per-task resumability check
# ---------------------------------------------------------------------------

def should_skip(
    *,
    output_csv: Path,
    grid_sig: str,
    force: bool,
    done_signatures: Set[str],
) -> Tuple[bool, str]:
    """
    Determine whether a task should be skipped.

    Parameters
    ----------
    output_csv:
        Expected output CSV path for this task.
    grid_sig:
        Grid signature for this task.
    force:
        If True, always run (never skip).
    done_signatures:
        Historical in-run completion records.  Kept for API compatibility and
        auditing, but never trusted without the concrete output artifact and
        per-task done metadata.

    Returns
    -------
    (skip, reason):
        ``skip=True`` if the task should be skipped.
        ``reason`` is a human-readable string.
    """
    if force:
        return False, ""

    # A non-empty CSV may only be reused if its recorded per-task grid
    # signature matches the current effective grid + fixed parameters and the
    # task metadata explicitly says event="done".  If the artifact is absent,
    # an old JSONL done record is not enough: rerun rather than silently
    # producing a skipped task with no output.
    if output_csv.exists() and output_csv.stat().st_size > 0:
        existing_sig = _grid_sig_from_meta(output_csv.parent)
        if existing_sig is None:
            return False, "stale_csv_meta_missing_or_not_done"
        if existing_sig != grid_sig:
            return False, "stale_csv_grid_signature_mismatch"
        return True, "output_csv_exists"

    return False, ""


# ---------------------------------------------------------------------------
# Cross-run CSV reuse (resume_scope=fixed)
# ---------------------------------------------------------------------------

def _grid_sig_from_meta(tb_dir: Path) -> Optional[str]:
    """Read grid signature from scan_meta.json; only accept event=done."""
    meta = load_json_best_effort(tb_dir / "scan_meta.json")
    if not meta:
        return None
    if meta.get("event") != "done":
        return None
    sig = meta.get("grid_signature")
    return sig.strip() if isinstance(sig, str) and sig.strip() else None


def find_previous_csv(
    fixed_dir: Path,
    current_run_dir: Path,
    tb_tag: str,
    desired_grid_sig: str,
) -> Optional[Path]:
    """
    Search all previous runs under *fixed_dir* for a CSV matching
    the given *tb_tag* and *desired_grid_sig*.

    Only returns non-empty CSVs whose associated per-task scan_meta.json
    carries the exact grid signature and event=done.  A run-level manifest is
    intentionally insufficient because it cannot prove that a particular task
    finished successfully.

    Parameters
    ----------
    fixed_dir:
        The fixed-param directory (parent of the current run dir).
    current_run_dir:
        Excluded from the search (avoid reusing from the current run).
    tb_tag:
        tan(beta) tag, e.g. ``"10000"``.
    desired_grid_sig:
        Grid signature that the found CSV must match.

    Returns
    -------
    Path or None
        Path to the best matching CSV (most recently modified), or None.
    """
    pattern = f"*/tb_*/scan_tb_{tb_tag}.csv"
    candidates = list(fixed_dir.glob(pattern))

    best: Optional[Tuple[float, Path]] = None

    for csv_path in candidates:
        try:
            # Skip anything inside the current run directory
            if csv_path.is_relative_to(current_run_dir):
                continue

            st = csv_path.stat()
            if st.st_size <= 0:
                continue

            tb_dir = csv_path.parent
            sig = _grid_sig_from_meta(tb_dir)
            if sig != desired_grid_sig:
                continue

            mtime = st.st_mtime
            if best is None or mtime > best[0]:
                best = (mtime, csv_path)

        except Exception:
            continue

    return best[1] if best else None


def materialize_csv(src: Path, dst: Path) -> bool:
    """
    Copy *src* CSV into *dst* (materialize from a previous run).

    Returns True on success, False on any error.
    """
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        return True
    except Exception:
        return False
