#!/usr/bin/env python3
"""Generate resumable, non-certified T2 arithmetic and trajectory artifacts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from decimal import Decimal, localcontext
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
DIHIGGS = REPO
sys.path.insert(0, str(DIHIGGS))
from dihiggs.app.orchestrator.aligned_lambda7 import (  # noqa: E402
    build_aligned_lambda7_family, xstar_lambda6_tanbeta,
)

CAMPAIGN = Path(__file__).resolve().parents[1]
BASELINE = CAMPAIGN / "baseline"
OUT = CAMPAIGN / "t2_lightweight"


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source_hash() -> str:
    digest = hashlib.sha256()
    for path in (DIHIGGS / "dihiggs/app/orchestrator/aligned_lambda7.py", Path(__file__)):
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _inputs(t, L, mphi, mh="125.13"):
    return dict(tan_beta=str(t), lambda1_ref=str(L), mh_GeV=mh,
                mphi_GeV=str(mphi), mA_GeV=str(mphi + 50),
                mHp_GeV=str(mphi + 50), v_GeV="246")


def _screen(point):
    vals = point.values
    max_lam = max(abs(Decimal(vals[f"lambda{i}"])) for i in range(1, 8))
    max_z = max(abs(Decimal(vals[f"Z{i}"])) for i in range(1, 8))
    return "PASS" if max(max_lam, max_z) < Decimal("12.566370614359172953850573533118") else "FAIL"


def _point_row(slice_name, t, L, mphi, trajectory, selector, *, mA=None, mHp=None):
    args = _inputs(t, L, mphi)
    if mA is not None:
        args["mA_GeV"] = str(mA)
    if mHp is not None:
        args["mHp_GeV"] = str(mHp)
    if selector == "xstar":
        X = xstar_lambda6_tanbeta(**{k: args[k] for k in
            ("tan_beta", "lambda1_ref", "mh_GeV", "mphi_GeV", "v_GeV")})
        kwargs = {"X_lambda6_tanbeta": X}
        x_value = X
    elif selector == "fixed_x":
        x_value = "1"
        kwargs = {"X_lambda6_tanbeta": x_value}
    elif selector == "zero_x":
        x_value = "0"
        kwargs = {"X_lambda6_tanbeta": x_value}
    else:
        lambda6 = "1e-10"
        kwargs = {"lambda6": lambda6}
        # Fixed lambda6 is a distinct trajectory; retain its evaluated X.
        x_value = str(Decimal(lambda6) * Decimal(str(t)))
    point = build_aligned_lambda7_family(**args, **kwargs, precision=80)
    row = {
        "point_id": hashlib.sha256(f"{slice_name}|{t}|{L}|{mphi}|{trajectory}".encode()).hexdigest()[:20],
        "slice": slice_name, "tan_beta": str(t), "lambda1_ref": str(L),
        "mh_GeV": "125.13", "mphi_GeV": str(mphi), "mA_GeV": args["mA_GeV"],
        "mHp_GeV": args["mHp_GeV"], "v_GeV": "246", "trajectory": trajectory,
        "evaluated_X": x_value, "lambda6": point.values["lambda6"],
        **{f"lambda{i}": point.values[f"lambda{i}"] for i in range(1, 8)},
        **{f"Z{i}": point.values[f"Z{i}"] for i in range(1, 8)},
        "g_hphiphi_GeV": point.values["g_hphiphi_GeV"],
        "quartic_screen_4pi": _screen(point),
        "spectrum_reconstruction": point.constraints["spectrum_reconstruction"]["status"],
        "alignment": point.constraints["alignment"]["status"],
        "stationarity": point.constraints["stationarity"]["status"],
        "local_minimum_tree": point.constraints["local_minimum_tree"]["status"],
        "full_tree_unitarity": point.constraints["full_tree_unitarity"]["status"],
        "BFB_general": point.constraints["BFB_general"]["status"],
        "global_EW_minimum": point.constraints["global_EW_minimum"]["status"],
        "widths": point.constraints["widths"]["status"],
        "branching_fractions": point.constraints["branching_fractions"]["status"],
        "ctau_mm": point.constraints["ctau_mm"]["status"],
        "displaced_photon_selection": point.constraints["displaced_photon_selection"]["status"],
        "candidate_validity": "NOT_TESTED",
        "point_class": "screened; not certified",
        "rejected_failure_reason": "" if _screen(point) == "PASS" else "quartic_magnitude_screen_4pi",
    }
    return row


def generate_points(path: Path):
    """Append missing rows by stable point id; safe to resume after interruption."""
    rows = []
    for t in (100, 1000, 10000, 100000, 1000000):
        for mphi in (150, 200, 250):
            for L in (1, 1.2):
                for selector, label in (("zero_x", "fixed_X_0"),
                                        ("fixed_x", "fixed_X_1"), ("xstar", "Xstar")):
                    rows.append(_point_row("illustrative_90", t, L, mphi, label, selector))
    for historical_mass in (450, 500):
        for t in (100, 1000, 10000, 100000, 1000000):
            for mphi in (150, 200, 250):
                for L in (1, 1.2):
                    for selector, label in (("zero_x", "fixed_X_0"),
                                            ("fixed_x", "fixed_X_1"), ("xstar", "Xstar")):
                        row = _point_row(f"T1_repeat_mA_{historical_mass}", t, L,
                                         mphi, label, selector, mA=historical_mass,
                                         mHp=historical_mass)
                        row["historical_anchor_status"] = (
                            "separate legacy anchor validated; this T1 trajectory row is not the anchor"
                        )
                        rows.append(row)
    # A distinct fixed-lambda6=1e-10 exercise across the illustrative mass grid.
    for t in (100, 1000, 10000, 100000, 1000000):
        for mphi in (150, 200, 250):
            for L in (1, 1.2):
                rows.append(_point_row("fixed_lambda6_1e-10", t, L, mphi,
                                       "fixed_lambda6_1e-10", "fixed_lambda6"))
    source_hash = _source_hash()
    expected_ids = {row["point_id"] for row in rows}
    grid_hash = hashlib.sha256("\n".join(sorted(expected_ids)).encode()).hexdigest()
    expected = {}
    for row in rows:
        row["source_sha256"] = source_hash
        row["grid_sha256"] = grid_hash
        expected[row["point_id"]] = row
    existing = {}
    if path.exists():
        with path.open(encoding="utf-8") as f:
            for line_no, line in enumerate(f, 1):
                if not line.strip():
                    continue
                try:
                    record = json.loads(line)
                    point_id = record["point_id"]
                except (json.JSONDecodeError, KeyError) as exc:
                    raise ValueError(f"invalid checkpoint record at {path}:{line_no}") from exc
                if point_id in existing:
                    raise ValueError(f"duplicate checkpoint point_id {point_id} at line {line_no}")
                if point_id not in expected:
                    raise ValueError(f"checkpoint point_id {point_id} is not in the declared grid")
                if record != expected[point_id]:
                    raise ValueError(f"checkpoint record conflicts with current source/grid: {point_id}")
                existing[point_id] = record
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        for row in rows:
            if row["point_id"] not in existing:
                f.write(json.dumps(row, sort_keys=True) + "\n")
                f.flush()
                existing[row["point_id"]] = row
    return len(existing)


def run_t1_verifier():
    verifier = BASELINE / "verify_lambda7_T1.py"
    with tempfile.TemporaryDirectory(prefix="t2-t1-verify-") as temp:
        tempdir = Path(temp)
        shutil.copy2(verifier, tempdir / verifier.name)
        done = subprocess.run([sys.executable, str(tempdir / verifier.name)],
                              text=True, capture_output=True, check=True)
        report = json.loads((tempdir / "T1_verification.json").read_text())
        if report["exact_symbolic_checks_passed"] != 34 or report["numeric_controls"] != 60:
            raise RuntimeError("T1 verifier did not reproduce 34 checks and 60 controls")
        shutil.copy2(tempdir / "T1_arithmetic_controls.csv", OUT / "T1_arithmetic_controls.csv")
        shutil.copy2(tempdir / "T1_verification.json", OUT / "T1_verification.json")
        (OUT / "T1_verifier_stdout.txt").write_text(done.stdout)
        return report


def precision_audit():
    pars = _inputs(1_000_000, "1.2", 150)
    point = build_aligned_lambda7_family(**pars, X_lambda6_tanbeta="1", precision=120)
    t = Decimal(pars["tan_beta"])
    mh, mp, v = (Decimal(pars[x]) for x in ("mh_GeV", "mphi_GeV", "v_GeV"))
    L = Decimal(pars["lambda1_ref"])
    l6 = Decimal(point.values["lambda6"])
    scale = max(v * v, mh * mh, mp * mp)
    out = {"scope": "arithmetic precision audit; actual calculator routes are audited separately in T4",
           "case": pars | {"X_lambda6_tanbeta": "1"}, "dimensioned_residual_scale_GeV2": str(scale),
           "constructed_lambda1_ref": str(L), "rows": []}
    reference = {key: point.values[key] for key in
                 [*(f"lambda{i}" for i in range(1,8)), *(f"Z{i}" for i in range(1,8)),
                  "m12_sq_GeV2", "g_hphiphi_GeV", "CphiHpHm_potential_GeV",
                  "tadpole_max_abs_normalized", "even_hessian_max_abs_normalized",
                  "CP_odd_hessian_max_abs_normalized", "charged_hessian_max_abs_normalized"]}
    for digits in (50, 80, 120):
        with localcontext() as ctx:
            ctx.prec = digits
            precision_point = build_aligned_lambda7_family(**pars, X_lambda6_tanbeta="1", precision=digits)
            stable_l7 = Decimal(precision_point.values["lambda7"])
            # Reconstruct lambda1 from the stable family quartics via the inverse map.
            ap_inverse = (mh * mh / (v * v) + (mp * mp - mh * mh / 2) * t * t / (v * v)
                          - Decimal("1.5") * l6 * t + stable_l7 * t**3 / 2)
            float_l7 = Decimal.from_float(float(stable_l7))
            serialized_l7 = Decimal(str(float(stable_l7)))
            naive = (float(mh)**2 / float(v)**2 + (float(mp)**2 - float(mh)**2 / 2)
                     * float(t)**2 / float(v)**2 - 1.5 * float(l6) * float(t)
                     + 0.5 * float(stable_l7) * float(t)**3)
            out["rows"].append({
                "digits": digits, "stable_constructed_lambda1": str(L),
                "forward_constructor_precision_digits": digits,
                "forward_values": {key: precision_point.values[key] for key in reference},
                "distance_from_120_digit_reference": {key: str(abs(Decimal(precision_point.values[key])-Decimal(reference[key])))
                    for key in reference},
                "arbitrary_precision_inverse_lambda1": str(ap_inverse),
                "absolute_residual": str(abs(ap_inverse - L)),
                "normalized_residual": str(abs(ap_inverse - L) / max(Decimal(1), abs(ap_inverse))),
                "double_inverse_lambda1": repr(naive),
                "double_inverse_abs_residual": repr(abs(naive - float(L))),
                "lambda7_exact_to_double_loss": str(abs(float_l7 - stable_l7)),
                "lambda7_float_lexeme_roundtrip_loss": str(abs(serialized_l7 - stable_l7)),
                "lambda7_perturbation_sensitivity": {
                    str(eps): str(abs(Decimal(eps) * t**3 / 2))
                    for eps in ("1e-24", "1e-20", "1e-16")},
                "measured_double_inverse_error_budget": "1.3e-5 absolute lambda1 for this arithmetic case; not a backend budget",
                "refinement_status": "INCONCLUSIVE; see T4 actual calculator route results",
            })
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-verifier", action="store_true")
    parser.add_argument("--points-file", type=Path, default=None,
                        help="checkpoint JSONL path (default: campaign t2_lightweight/trajectory_points.jsonl)")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rows_path = args.points_file or (OUT / "trajectory_points.jsonl")
    count = generate_points(rows_path)
    report = {"campaign": CAMPAIGN.name, "status": "LIGHTWEIGHT_ONLY",
              "points_written_or_replayed": count,
              "points_file": str(rows_path), "checkpoint_source_sha256": _source_hash(),
              "point_counts": {"illustrative_90": 90, "T1_repeat_mA_450_slice": 90,
                               "T1_repeat_mA_500_slice": 90, "fixed_lambda6_1e-10": 30},
              "screen_semantics": "trajectory JSONL has arithmetic and quartic magnitude screen only; certified BFB/unitarity/global classifications are in physical_constraints_report.jsonl",
              "heavy_execution": "BLOCKED: lab-health and safe-run unavailable; no full pilot calculator sweep. A bounded nine-route T4 anchor probe was compiled and run separately.",
              "inputs": {"t_values": [100, 1000, 10000, 100000, 1000000],
                         "lambda1_ref": [1, 1.2], "illustrative_mphi_GeV": [150, 200, 250],
                         "illustrative_mA_mHp": "mphi+50", "trajectory_X": [0, 1, "Xstar"],
                         "historical_mass_slices_mA_mHp_GeV": [450, 500],
                         "validated_historical_anchors": {
                             "450": "verification/historical/historical450.json (mH=150, tb=300000, M2=22499.999999500335, lambda6=1e-10, lambda7=0)",
                             "500": "verification/historical/historical500.json (mH=200, tb=10000, M2=39999.9995713761, lambda6=1e-10, lambda7=0)"},
                         "fixed_lambda6_investigation": "1e-10"},
              "precision_audit": precision_audit(),
              "source_sha256": {"contract": _digest(BASELINE / "CODEX_CONTRACT_T1.md"),
                                "verifier": _digest(BASELINE / "verify_lambda7_T1.py")}}
    if not args.skip_verifier:
        report["T1_verifier"] = run_t1_verifier()
    (OUT / "t2_lightweight_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"output": str(OUT), "rows": count,
                      "T1_verifier": report.get("T1_verifier", {}).get("exact_symbolic_checks_passed", "skipped")}, indent=2))


if __name__ == "__main__":
    main()
