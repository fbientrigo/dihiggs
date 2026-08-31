from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git(root: Path, *args: str) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(root), *args], text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def repository_identity(root: Path) -> dict[str, Any]:
    status = _git(root, "status", "--porcelain")
    return {
        "repository_root": str(root),
        "repository_commit": _git(root, "rev-parse", "HEAD"),
        "repository_dirty": bool(status),
        "repository_status_sha256": hashlib.sha256(status.encode()).hexdigest(),
    }


def evaluator_identity(root: Path, executable: Path) -> dict[str, Any]:
    return {
        "executable": str(executable),
        "evaluator_executable_sha256": sha256_file(executable),
        "evaluator_source_sha256": sha256_file(root / "dihiggs/src/DihiggsPointV2Evaluator.cpp"),
        "twohdmc_library_sha256": sha256_file(root / "2hdmc/lib/lib2HDMC.a"),
        "dependency_identity": "vendored-2HDMC-1.8.0-plus-local-patch",
    }

