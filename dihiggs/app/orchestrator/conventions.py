"""Canonical physics conventions for the orchestrator.

Single in-repo accessor for the constants that the orchestrator used to
copy-paste as literals (``mh = 125.13`` appeared in five separate files, and
``engines/m2.py`` disagreed with ``manifest.py`` as a result). Values come from
the ecosystem-wide ``conventions/physics_conventions.yaml`` at the repo root,
which is kept byte-identical with the copies in dihiggs_boundary and
dihiggs_hep_cross.

The pinned literals below are the fallback when PyYAML or the file is
unavailable; ``tests/test_orchestrator/test_conventions.py`` asserts the loaded
values equal these literals so the two can never silently disagree.

Mirrors ``mlpython/lake_pipeline/physics_conventions.py`` and
``dihiggs_hep_cross/src/llp_recast/constants.py``.
"""

import os

# Pinned fallback (must equal conventions/physics_conventions.yaml).
# Decimal string, not float: the float64 repr of 125.20 is 1.25200000000000003e+02,
# which is not byte-equal to the 1.25200000000000000e+02 form SLHA/UFO cards need.
_M_H_GEV_TEXT_PINNED = "125.20"
_SOURCE_PINNED = "PDG 2026 Higgs listing"
_SOURCE_URL_PINNED = "https://pdg.lbl.gov/encoder_listings/s126.pdf"

# dihiggs/app/orchestrator/conventions.py -> repo root is four levels up.
_REPO_ROOT = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
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
_SM_LIKE = _CONV.get("sm_like_higgs") or {}

# Canonical SM-like Higgs mass for ACTIVE physical-point production.
M_H_GEV_TEXT = str(_SM_LIKE.get("m_h_GeV", _M_H_GEV_TEXT_PINNED))
M_H_GEV = float(M_H_GEV_TEXT)

M_H_SOURCE = str(_SM_LIKE.get("source", _SOURCE_PINNED))
M_H_SOURCE_URL = str(_SM_LIKE.get("source_url", _SOURCE_URL_PINNED))


def mass_convention_block(mh=None):
    """Provenance block recorded in run manifests.

    ``mh`` overrides the canonical value (an explicit ``--mh``); when it differs
    from the convention the block says so, so a manifest can never claim the
    canonical convention while the evaluator ran on something else.
    """
    value = M_H_GEV if mh is None else float(mh)
    block = {
        "mh_GeV": value,
        "source": M_H_SOURCE,
        "source_url": M_H_SOURCE_URL,
    }
    if value != M_H_GEV:
        block["source"] = "explicit --mh override"
        block["canonical_mh_GeV"] = M_H_GEV
        block["canonical_source"] = M_H_SOURCE
    return block
