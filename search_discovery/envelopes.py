"""Tier-2 master-controlled adaptive envelopes.

E0 is deliberately WIDER than the region master calibration found to be valid.
It spans the positivity boundary (X ~ 0.172) and the unitarity/perturbativity
boundary (|Q| ~ 3e5) from both sides, so the calibration wave *measures* those
boundaries with the frozen evaluator instead of inheriting them as assumptions.
Only the master may change the active envelope; workers read it.
"""
from __future__ import annotations

import json
from pathlib import Path

from .bounds import Envelope

ENVELOPES: dict[str, Envelope] = {
    "E0": Envelope(
        envelope_id="E0",
        rationale=("calibration wave: spans both known theory boundaries from both sides "
                   "so the frozen evaluator maps them; covers the full 150-250 GeV target "
                   "mass range and the tan_beta decades that reach ctau ~ 1e2-1e4 mm"),
        mH_GeV=(150.0, 250.0),
        mA_GeV=(155.0, 900.0),
        tan_beta=(1.0e4, 5.0e7),
        X=(0.0, 0.30),
        Q=(-1.0e7, 1.0e7),
    ),
    "E1_mixed_exploit": Envelope(
        envelope_id="E1_mixed_exploit",
        rationale=("exploitation envelope around the calibrated theory-valid mixed basin; "
                   "tan_beta narrowed to the decade that brackets ctau_200 = 500-1000 mm"),
        mH_GeV=(150.0, 250.0),
        mA_GeV=(210.0, 700.0),
        tan_beta=(2.0e6, 1.2e7),
        X=(0.0, 0.17),
        Q=(-2.0e5, 2.0e5),
    ),
    "E2_photonic_hunt": Envelope(
        envelope_id="E2_photonic_hunt",
        rationale=("photonic-family hunt: pushes |Q| across the unitarity boundary and "
                   "down to low tan_beta, where the charged-Higgs-loop photon modes can "
                   "dominate. Deliberately includes theory-invalid ground so the frontier "
                   "between photon dominance and unitarity is mapped, not assumed"),
        mH_GeV=(150.0, 250.0),
        mA_GeV=(155.0, 1200.0),
        tan_beta=(1.0e3, 1.0e8),
        X=(0.0, 0.20),
        Q=(-1.0e9, 1.0e9),
    ),
}
DEFAULT_ENVELOPE_ID = "E0"


def save_envelope(run_dir: Path, envelope: Envelope) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "active_envelope.json").write_text(
        json.dumps(envelope.as_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def active_envelope(run_dir: Path) -> Envelope:
    path = Path(run_dir) / "active_envelope.json"
    if not path.exists():
        return ENVELOPES[DEFAULT_ENVELOPE_ID]
    document = json.loads(path.read_text(encoding="utf-8"))
    return Envelope(
        envelope_id=document["envelope_id"], rationale=document["rationale"],
        mH_GeV=tuple(document["mH_GeV"]), mA_GeV=tuple(document["mA_GeV"]),
        tan_beta=tuple(document["tan_beta"]), X=tuple(document["X"]), Q=tuple(document["Q"]),
    )
