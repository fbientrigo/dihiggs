"""
layout.py — Data-lake directory structure for scan outputs.

Folder hierarchy
----------------
<outdir>/<lake_name>/
    campaign=<campaign>/
        fixed_sinba=<...>_l6=<...>_l7=<...>_mA=<...>/   ← "fixed_dir"
            <run_name>/                                    ← "run_dir"
                run_manifest.json
                orchestrator.log
                task_summary.jsonl
                tb_<NNNNN>/                               ← per-tanbeta
                    scan_tb_<tag>.csv                     ← C++ output CSV
                    scan_meta.json                        ← per-task metadata
                    stdout.log   (on failure only)
                    stderr.log   (on failure only)

The fixed_dir name encodes the fixed physics parameters for deduplication
across campaigns.  It does NOT encode lambda1 or M2 axis parameters because
those belong to the grid signature.

Engine name is NOT in the fixed_dir because the same fixed physics can
legitimately be probed with both engines; the grid_signature differentiates.
"""

from __future__ import annotations

import os
import platform
import socket
from pathlib import Path
from typing import Tuple

from dihiggs.app.orchestrator.io_utils import (
    format_float_tag,
    sanitize_for_path,
    utc_now_compact,
)
from dihiggs.app.orchestrator.models import FixedParams


# ---------------------------------------------------------------------------
# Run-directory construction
# ---------------------------------------------------------------------------

def build_run_dir(
    outdir: Path,
    lake_name: str,
    campaign: str,
    fixed: FixedParams,
    run_name: str,
) -> Path:
    """
    Construct the full run directory path.

    Parameters
    ----------
    outdir:
        Root output directory (CERNBox, local NFS, etc.).
    lake_name:
        Subdirectory used as the data-lake root (e.g. ``"dihiggs_lake"``).
    campaign:
        Human label for the scan campaign (e.g. ``"scan_999"``).
    fixed:
        Fixed physics parameters; encoded into the fixed-param dir name.
    run_name:
        Unique run identifier (timestamp + host, or user-specified).

    Returns
    -------
    Path
        Full path: ``outdir/lake_name/campaign=.../fixed_.../run_name/``
    """
    lake_root = outdir / lake_name
    campaign_dir = lake_root / f"campaign={sanitize_for_path(campaign)}"

    # lambda6 / lambda7 are tiny couplings (~1e-3); 4 decimals aliases
    # distinct values (e.g. 0.0019683 and 0.0020273 both -> 0.0020), so use
    # higher precision for these two axes to keep run-dir identity unique.
    fixed_dir_name = (
        f"fixed_sinba={sanitize_for_path(format_float_tag(fixed.sin_ba, 4))}"
        f"_l6={sanitize_for_path(format_float_tag(fixed.lambda6, 7))}"
        f"_l7={sanitize_for_path(format_float_tag(fixed.lambda7, 7))}"
        f"_mA={sanitize_for_path(format_float_tag(fixed.mA, 1))}"
    )
    fixed_dir = campaign_dir / fixed_dir_name
    return fixed_dir / sanitize_for_path(run_name)


def default_run_name(git_short: str | None = None) -> str:
    """
    Generate a default run name from UTC timestamp, hostname, and git sha.

    Parameters
    ----------
    git_short:
        Short git SHA if available; ``"nogit"`` used otherwise.

    Returns
    -------
    str
        Example: ``"run_20260207T182031Z_host=myhost_git=abc1234"``
    """
    host = socket.gethostname()
    sha = git_short or "nogit"
    return f"run_{utc_now_compact()}_host={host}_git={sha}"


# ---------------------------------------------------------------------------
# Per-tanbeta folder naming  (preserved from legacy orchestrator)
# ---------------------------------------------------------------------------

def format_tanbeta_folder(tb: float) -> Tuple[str, str]:
    """
    Return ``(folder_name, file_tag)`` for a tan(beta) value.

    Integer values get zero-padded to at least 5 digits:
        10000 → ``("tb_10000", "10000")``

    Non-integer or negative values use a sanitised decimal representation:
        12345.6 → ``("tb_12345p6", "12345p6")``

    Parameters
    ----------
    tb:
        tan(beta) value.

    Returns
    -------
    folder_name:
        Directory name (e.g. ``"tb_10000"``).
    tag:
        File-name component (e.g. ``"10000"``).
    """
    if abs(tb - round(tb)) < 1e-12 and tb >= 0:
        tb_int = int(round(tb))
        folder = f"tb_{tb_int:05d}"
        tag = str(tb_int)
        return folder, tag

    raw = f"{tb:.6g}"
    safe = sanitize_for_path(raw)
    return f"tb_{safe}", safe


# ---------------------------------------------------------------------------
# Host metadata
# ---------------------------------------------------------------------------

def host_metadata() -> dict:
    """Return a dict of host-level metadata for inclusion in manifests."""
    import sys as _sys
    return {
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python": _sys.version.split()[0],
        "user": os.environ.get("USER") or os.environ.get("USERNAME"),
    }
