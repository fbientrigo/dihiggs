"""Three-tier bounds model for the discovery campaign.

Tier 1  GLOBAL_PHYSICS_BOUNDS   immutable physics/safety limits. Never widened by
                                any agent, lead, or master action.
Tier 2  Envelope                master-controlled adaptive search region. Must be
                                contained in tier 1. May be re-issued per epoch.
Tier 3  worker-local            proposal neighbourhoods around a parent; enforced
                                by the perturbation helpers, always clipped to the
                                active tier-2 envelope.

Search coordinates are (mH_GeV, mA_GeV, tan_beta, X, Q), NOT raw (M2, lambda6).
This is deliberate: master calibration established that the frozen evaluator's
theory gates respond to the invariant Q = (mH^2 - M2) * tan_beta^2 rather than to
M2 itself, and that positivity responds to X = lambda6 * tan_beta through the
measured relation lambda1 = m_h^2/v^2 - 1.5 * X.  Sampling in (Q, X) is therefore
well conditioned where sampling in (M2, lambda6) is not.

Derived, exactly as in the frozen substrate:
    lambda6 = X / tan_beta
    M2_GeV2 = mH_GeV**2 - Q / tan_beta**2
    mHp_GeV = mA_GeV
"""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from typing import Any

# ---------------------------------------------------------------- tier 1
# Immutable. Justification for each limit is recorded in
# docs/DISCOVERY_CALIBRATION_EVIDENCE.md.
GLOBAL_PHYSICS_BOUNDS: dict[str, tuple[float, float]] = {
    # scientific target masses for the 150/200/250 family continuation
    "mH_GeV": (150.0, 250.0),
    # mA = mHp. Must exceed mH or H2 -> A Z opens and the state is prompt.
    "mA_GeV": (150.0, 1200.0),
    # LLP regime needs cot(beta)-suppressed fermionic widths; below ~1e3 the
    # state is prompt for any Q. Upper limit is a numerical-safety ceiling.
    "tan_beta": (1.0e3, 1.0e8),
    # X = lambda6 * tan_beta. Master calibration measured lambda1 = m_h^2/v^2
    # - 1.5 * X exactly, so positivity requires X < ~0.172. That is a GATE
    # RESULT, not a bound: the ceiling here is deliberately permissive (it spans
    # the whole frozen v1 contract range, X in [200, 1e5]) so that the campaign
    # can re-measure and record the positivity wall rather than inherit it.
    "X": (0.0, 1.0e5),
    # Q = (mH^2 - M2) * tan_beta^2. Unitarity/perturbativity degrade above
    # |Q| ~ 3e5 in the calibrated slice. Again permissive on purpose: the frozen
    # theory gates decide validity, never this bound.
    "Q": (-1.0e15, 1.0e15),
}

# Hard structural requirement, independent of any envelope.
MIN_MA_MINUS_MH_GEV = 1.0
# lambda6 must stay resolvable in double precision after division by tan_beta.
MIN_LAMBDA6 = 1.0e-14


class BoundsError(ValueError):
    pass


@dataclass(frozen=True)
class Envelope:
    """Tier-2 master-controlled adaptive search envelope."""

    envelope_id: str
    rationale: str
    mH_GeV: tuple[float, float]
    mA_GeV: tuple[float, float]
    tan_beta: tuple[float, float]
    X: tuple[float, float]
    Q: tuple[float, float]

    def __post_init__(self) -> None:
        for name in ("mH_GeV", "mA_GeV", "tan_beta", "X", "Q"):
            low, high = getattr(self, name)
            glow, ghigh = GLOBAL_PHYSICS_BOUNDS[name]
            if not (math.isfinite(low) and math.isfinite(high)):
                raise BoundsError(f"nonfinite_envelope:{name}")
            if low > high:
                raise BoundsError(f"inverted_envelope:{name}")
            if low < glow or high > ghigh:
                raise BoundsError(
                    f"envelope_escapes_global_physics_bounds:{name}:"
                    f"[{low},{high}] not in [{glow},{ghigh}]"
                )

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    def span(self, name: str) -> float:
        low, high = getattr(self, name)
        return high - low


def derive(mH_GeV: float, mA_GeV: float, tan_beta: float, X: float, Q: float) -> dict[str, float]:
    """Map search coordinates to the frozen evaluator's physical parameters."""
    lambda6 = X / tan_beta
    M2 = mH_GeV * mH_GeV - Q / (tan_beta * tan_beta)
    return {
        "mH_GeV": float(mH_GeV),
        "mA_GeV": float(mA_GeV),
        "M2_GeV2": float(M2),
        "tan_beta": float(tan_beta),
        "lambda6": float(lambda6),
    }


def check_global(mH_GeV: float, mA_GeV: float, tan_beta: float, X: float, Q: float) -> None:
    values = {"mH_GeV": mH_GeV, "mA_GeV": mA_GeV, "tan_beta": tan_beta, "X": X, "Q": Q}
    for name, value in values.items():
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(float(value)):
            raise BoundsError(f"nonfinite_or_wrong_type:{name}")
        low, high = GLOBAL_PHYSICS_BOUNDS[name]
        if not low <= float(value) <= high:
            raise BoundsError(f"outside_global_physics_bounds:{name}:{value}")
    if mA_GeV - mH_GeV < MIN_MA_MINUS_MH_GEV:
        raise BoundsError(f"invalid_hierarchy:mA_minus_mH_below_{MIN_MA_MINUS_MH_GEV}")
    if X > 0.0 and X / tan_beta < MIN_LAMBDA6:
        raise BoundsError("lambda6_underflow_risk")


def clip_to_envelope(name: str, value: float, envelope: Envelope) -> float:
    low, high = getattr(envelope, name)
    return min(max(float(value), low), high)
