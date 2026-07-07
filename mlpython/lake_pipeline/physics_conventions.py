"""Shared physics conventions for the lake pipeline.

Single in-repo source of truth for constants that used to be copy-pasted into
individual lake_pipeline scripts (e.g. HBAR_C_GEV_MM). Values come from the
ecosystem-wide ``conventions/physics_conventions.yaml`` at the repo root, which
is kept byte-identical with the copies in dihiggs_boundary and
dihiggs_hep_cross so ctau derivations agree across all three repos.

The pinned literals below are the fallback when PyYAML or the file is
unavailable; ``tests/test_physics_conventions.py`` asserts the loaded values
equal these literals so the two can never silently disagree.
"""

import os

# Pinned fallback (must equal conventions/physics_conventions.yaml).
_HBAR_C_GEV_MM_PINNED = 1.973269804e-13  # c*tau[mm] = hbar*c / Gamma[GeV]
_C_MM_PER_NS_PINNED = 299.792458  # c*tau[mm] = C_MM_PER_NS * tau[ns]

# mlpython/lake_pipeline/physics_conventions.py -> repo root is three levels up.
_REPO_ROOT = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)
CONVENTIONS_PATH = os.path.join(_REPO_ROOT, "conventions", "physics_conventions.yaml")


def _load():
    try:
        import yaml
    except ImportError:
        return {}
    try:
        with open(CONVENTIONS_PATH) as fh:
            return yaml.safe_load(fh) or {}
    except OSError:
        return {}


_CONV = _load()

HBAR_C_GEV_MM = float(_CONV.get("hbar_c_gev_mm", _HBAR_C_GEV_MM_PINNED))
C_MM_PER_NS = float(_CONV.get("c_mm_per_ns", _C_MM_PER_NS_PINNED))


def ctau_mm_from_width_gev(total_width_gev):
    """c*tau in mm from a total width in GeV (>0). Mirrors the identical
    formula in dihiggs_hep_cross recast_math and the C++ evaluator in
    dihiggs_boundary."""
    return HBAR_C_GEV_MM / total_width_gev
