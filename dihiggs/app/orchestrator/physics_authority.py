"""Read the canonical v3 physics convention for new point production.

This module deliberately has one narrow responsibility: expose the active
SM-like Higgs mass and enough immutable provenance for a point manifest. It
does not select campaign choices or reinterpret historical artifacts.
"""

from __future__ import annotations

import hashlib
import subprocess
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Optional


_REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
_CONVENTIONS_PATH = _REPOSITORY_ROOT / "conventions" / "physics_conventions.yaml"


@lru_cache(maxsize=1)
def _load_contract() -> Dict[str, Any]:
    """Load the canonical contract or fail before producing a new point."""
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover - dependency is declared
        raise RuntimeError("PyYAML is required to load the physics authority") from exc

    try:
        contract = yaml.safe_load(_CONVENTIONS_PATH.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"cannot load physics authority: {_CONVENTIONS_PATH}") from exc

    if not isinstance(contract, dict):
        raise RuntimeError("physics authority must be a YAML mapping")
    if contract.get("schema_version") != "physics_conventions_v3":
        raise RuntimeError("unsupported physics authority schema")
    return contract


def active_mass_gev() -> float:
    """Return the active new-production mass, read from the v3 contract."""
    active = _load_contract()["mass_conventions"]["active"]
    if active.get("status") != "active_new_production":
        raise RuntimeError("physics authority has no active production convention")
    return float(active["value_GeV"])


def _authority_commit() -> str:
    """Return the commit containing the authority file when Git is available."""
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=_REPOSITORY_ROOT,
            capture_output=True,
            check=False,
            text=True,
        )
    except OSError:
        return "unknown"
    return completed.stdout.strip() if completed.returncode == 0 else "unknown"


def convention_provenance(
    *, source_commit: Optional[str], requested_mass_gev: Optional[float] = None
) -> Dict[str, Any]:
    """Return the v3 provenance block required on a new point handoff."""
    contract = _load_contract()
    active = contract["mass_conventions"]["active"]
    active_mass = float(active["value_GeV"])
    if requested_mass_gev is not None and requested_mass_gev != active_mass:
        raise RuntimeError(
            "new production must use the active v3 Higgs-mass convention; "
            "historical replay requires an explicit historical interface"
        )
    authority = contract["authority"]
    digest = hashlib.sha256(_CONVENTIONS_PATH.read_bytes()).hexdigest()
    return {
        "convention_id": active["id"],
        "m_h_GeV": active_mass,
        # `dihiggs.point.v2` already exposed this legacy spelling. Preserve it
        # until the point schema is versioned, while new consumers use m_h_GeV.
        "mh_GeV": active_mass,
        "m_h_GeV_text": active["value_GeV"],
        "source": "PDG 2026 Higgs-boson listing",
        "source_url": active["external_source"]["url"],
        "schema_version": contract["schema_version"],
        "source_repository": authority["repository"],
        "source_commit": source_commit or _authority_commit(),
        "source_path": authority["path"],
        "source_sha256": digest,
    }
