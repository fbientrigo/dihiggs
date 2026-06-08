"""
io_utils.py — Shared I/O helpers for the scan orchestrator.

Provides:
- safe_write_json   : atomic JSON write (tmp + rename)
- append_jsonl      : append one JSON record to a JSONL file
- load_json_best_effort : load JSON; return None on any error
- utc_now_iso       : ISO-8601 UTC timestamp string
- utc_now_compact   : compact UTC timestamp for folder naming
- detect_git_info   : best-effort git commit + dirty status
- log_msg           : timestamp-prefixed log to stdout + file

These utilities contain NO physics logic and have NO dependencies
beyond the Python standard library.
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional


# ---------------------------------------------------------------------------
# Timestamps
# ---------------------------------------------------------------------------

def utc_now_iso() -> str:
    """Return current UTC time as ISO-8601 string (seconds precision)."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def utc_now_compact() -> str:
    """Return compact UTC timestamp suitable for directory names, e.g. 20260207T182031Z."""
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


# ---------------------------------------------------------------------------
# JSON / JSONL helpers
# ---------------------------------------------------------------------------

def safe_write_json(path: Path, payload: Dict[str, Any]) -> None:
    """
    Atomically write *payload* as JSON to *path*.

    Uses a temporary sibling file + rename so that readers never see a
    partially-written file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=False)
        fh.write("\n")
    tmp.replace(path)


def append_jsonl(path: Path, record: Dict[str, Any]) -> None:
    """
    Append one JSON record as a line to the JSONL file at *path*.

    Creates the file (and parent directories) if they do not exist.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(record) + "\n")


def load_json_best_effort(path: Path) -> Optional[Dict[str, Any]]:
    """
    Load and parse a JSON file.

    Returns *None* if the file does not exist, cannot be read, or is
    not valid JSON.  Never raises.
    """
    try:
        if not path.exists():
            return None
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def log_msg(msg: str, log_file: Optional[Path] = None) -> None:
    """
    Write a timestamped message to *stdout* and, optionally, to *log_file*.

    Parameters
    ----------
    msg:
        Message text (without leading timestamp).
    log_file:
        If given, append the timestamped line to this file.
    """
    ts = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    full = f"{ts} {msg}"
    print(full, flush=True)
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with log_file.open("a", encoding="utf-8") as fh:
            fh.write(full + "\n")


# ---------------------------------------------------------------------------
# Git metadata
# ---------------------------------------------------------------------------

def detect_git_info(start_dir: Path) -> Dict[str, Optional[str]]:
    """
    Best-effort git metadata: commit hash, short hash, dirty flag, repo root.

    Never raises; returns ``None`` values when git is unavailable.
    """
    def _run(*args: str) -> Optional[str]:
        try:
            r = subprocess.run(
                ["git", *args],
                cwd=str(start_dir),
                capture_output=True,
                text=True,
                check=False,
            )
            return r.stdout.strip() if r.returncode == 0 else None
        except Exception:
            return None

    repo_root = _run("rev-parse", "--show-toplevel")
    commit = _run("rev-parse", "HEAD")
    dirty = _run("status", "--porcelain")

    return {
        "repo_root": repo_root,
        "commit": commit,
        "commit_short": commit[:7] if commit else None,
        "is_dirty": None if dirty is None else ("yes" if dirty else "no"),
    }


# ---------------------------------------------------------------------------
# String helpers
# ---------------------------------------------------------------------------

def sanitize_for_path(s: str) -> str:
    """
    Replace characters that are unsafe in filesystem paths.

    ``.`` → ``p``, ``-`` → ``m``, everything else non-alnum/_/= → ``_``.
    """
    import re
    s = s.replace(".", "p").replace("-", "m")
    s = re.sub(r"[^A-Za-z0-9_=]+", "_", s)
    return s.strip("_")


def format_float_tag(x: float, ndp: int = 4) -> str:
    """Format a float for stable, deterministic use in path components."""
    return f"{x:.{ndp}f}"


def parse_csv_floats(s: str) -> List[float]:
    """
    Parse a comma-separated list of floats from string *s*.

    Examples
    --------
    ``"100,200,500"`` → ``[100.0, 200.0, 500.0]``
    """
    if not s.strip():
        return []
    return [float(p.strip()) for p in s.split(",") if p.strip()]
