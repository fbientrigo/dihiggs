"""
grid.py — Scan-grid descriptor and engine-aware grid signatures.

The ``ScanGrid`` dataclass describes the two-dimensional parameter grid:
- mphi axis  (m_phi in GeV) — always the same across engines
- axis       (lambda1 or M2/m12_sq depending on engine) — second dimension

The ``grid_signature()`` function hashes the grid + fixed params + engine name
into a 16-character hex string.  Including the engine name is **required** to
prevent resume collisions between lambda1 and M2 scans that happen to use
the same numeric ranges.

Semantics note
--------------
``axis_min`` / ``axis_max`` are deliberately generic.  Their physical meaning
depends on which engine is selected:

    engine=lambda1  →  axis is lambda_1       (dimensionless, range ~[0, 12])
    engine=m2       →  axis is M^2 = m12_sq/(sin(beta)cos(beta)) (GeV^2)

The ``ScanGrid`` itself stores no unit information; the EngineAdapter's
``axis_metadata()`` method is the authoritative unit source.  Never infer
units from column names or numeric ranges alone.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from dihiggs.app.orchestrator.models import FixedParams


@dataclass(frozen=True)
class ScanGrid:
    """
    Two-dimensional scan grid.

    Attributes
    ----------
    mphi_min, mphi_max:
        m_phi scan range in GeV.
    n_mphi:
        Number of m_phi grid points.
    axis_min, axis_max:
        Second axis scan range.  Units depend on the engine:
        - lambda1 engine: dimensionless (lambda_1)
        - M2 engine: GeV^2 (M^2 = m12_sq/(sin(beta)cos(beta)))
    n_axis:
        Number of second-axis grid points.
    """

    mphi_min: float
    mphi_max: float
    n_mphi: int
    axis_min: float
    axis_max: float
    n_axis: int


def grid_signature(engine_name: str, grid: "ScanGrid", fixed: "FixedParams") -> str:
    """
    Compute a 16-character hex signature for a (engine, grid, fixed) triple.

    The engine name is included to prevent resume collisions between
    lambda1 and M2 scans that share the same numeric ranges.

    Parameters
    ----------
    engine_name:
        Canonical engine identifier (e.g. ``"lambda1"`` or ``"m2"``).
    grid:
        Scan grid specification.
    fixed:
        Fixed physics parameters.  ``fixed.lambda1`` is included when
        present (relevant for M2 scans where lambda1 is a fixed constant).

    Returns
    -------
    str
        16-character lowercase hex digest.
    """
    # Build a deterministic dict; sort_keys=True ensures stable ordering.
    payload: dict = {
        "engine": engine_name,
        "mphi_min": grid.mphi_min,
        "mphi_max": grid.mphi_max,
        "n_mphi": grid.n_mphi,
        "axis_min": grid.axis_min,
        "axis_max": grid.axis_max,
        "n_axis": grid.n_axis,
        "mA": fixed.mA,
        "sin_ba": fixed.sin_ba,
        "tan_beta": fixed.tan_beta,
        "lambda6": fixed.lambda6,
        "lambda7": fixed.lambda7,
    }
    # Include lambda1 only when it is a fixed constant (M2 engine)
    if fixed.lambda1 is not None:
        payload["lambda1_fixed"] = fixed.lambda1
    if fixed.mh is not None:
        payload["mh"] = fixed.mh
    if fixed.mHp is not None:
        payload["mHp"] = fixed.mHp
    if fixed.yukawa_type is not None:
        payload["yukawa_type"] = fixed.yukawa_type

    # gen_fixings: fold in bronze shard identity + calibration config so
    # distinct shards/configs don't collide on the dedup key.
    if fixed.bronze_shard_csv is not None:
        payload["bronze_shard_csv"] = fixed.bronze_shard_csv
    if fixed.calibration_n is not None:
        payload["calibration_n"] = fixed.calibration_n
    if fixed.calibration_frac is not None:
        payload["calibration_frac"] = fixed.calibration_frac

    blob = json.dumps(payload, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]
