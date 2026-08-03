#!/usr/bin/env python3
"""Check the H2 construction coordinate and double-representation stability.

Only the exact m12_2 value produced by set_param_phys_lam1 and its immediately
adjacent representable doubles are compared. Wider offsets change reconstructed
lambda1 and are physical sensitivity tests, not numerical-error estimates.
"""
import csv
import json
import math
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CPP_SRC = ROOT / "benchmarks/check_H2scan_mH150_tb300000.cpp"
CPP_BIN = ROOT / "benchmarks/.check_H2scan_mH150_tb300000_bin"
CANONICAL_CSV = ROOT / "benchmarks/first_h2_bounded_scan.csv"
OUT_CSV = ROOT / "benchmarks/H2scan_mH150_tb300000_numerical_check.csv"
OUT_MD = ROOT / "benchmarks/H2scan_mH150_tb300000_numerical_check.md"
POINT_ID = "H2scan_mH150_tb300000"

MH, MH2, MA, MHP = 125.13, 150.0, 450.0, 450.0
SBA, L1, L6, L7, TB, GF = 1.0, 1.0, 1e-10, 0.0, 300000.0, 1.16637e-5
FIELDS = ("total_width_gev", "ctau_mm", "br_bb", "br_gammagamma", "br_Zgamma", "br_tautau", "br_gg")
REPLAY_FIELDS = ("total_width_gev", "ctau_mm", "br_bb", "br_gammagamma", "br_Zgamma")
CENTER_REPRO_REL_TOL = 1e-8
LAMBDA1_REPLAY_ABS_TOL = 1e-6


def faithful_m12_2() -> float:
    beta = math.atan(TB)
    cb, tb = math.cos(beta), math.tan(beta)
    alpha = -math.asin(SBA) + beta
    sa, ca = math.sin(alpha), math.cos(alpha)
    v2 = 1.0 / (math.sqrt(2.0) * GF)
    bracket = L1 + 1.5 * L6 * tb - 0.5 * L7 * tb**3
    return (MH2**2 * ca**2 + MH**2 * sa**2 - v2 * cb**2 * bracket) / tb


def sensitivity() -> float:
    beta = math.atan(TB)
    cb, tb = math.cos(beta), math.tan(beta)
    return tb / ((1.0 / (math.sqrt(2.0) * GF)) * cb**2)


def construction_M2(m12_2: float) -> float:
    beta = math.atan(TB)
    return m12_2 / (math.sin(beta) * math.cos(beta))


def canonical_row() -> dict[str, str]:
    with CANONICAL_CSV.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row["point_id"] == POINT_ID:
                return row
    raise SystemExit(f"{POINT_ID} not found in {CANONICAL_CSV}")


def compile_checker() -> None:
    subprocess.run([
        "g++", "-I", str(ROOT / "2hdmc/src"), "-std=gnu++17", "-O0", "-g",
        str(CPP_SRC), "-o", str(CPP_BIN), "-L", str(ROOT / "2hdmc/lib"),
        "-l2HDMC", "-lgsl", "-lgslcblas", "-lm",
    ], cwd=ROOT, check=True)


def evaluate(m12_2: float) -> dict[str, str]:
    result = subprocess.run([
        str(CPP_BIN), str(MH), str(MH2), str(MA), str(MHP), str(SBA),
        str(L6), str(L7), f"{m12_2:.17e}", str(TB),
    ], cwd=ROOT, check=True, capture_output=True, text=True)
    values = result.stdout.strip().split(",")
    if len(values) % 2:
        raise RuntimeError(f"unexpected checker output: {result.stdout!r}")
    return dict(zip(values[0::2], values[1::2]))


def relative_spread(values: list[float]) -> float:
    return (max(values) - min(values)) / max(max(abs(v) for v in values), 1e-300)


def finite_value(probe: dict[str, str], field: str) -> bool:
    try:
        return math.isfinite(float(probe[field]))
    except (KeyError, TypeError, ValueError):
        return False


def replay_gate(row: dict[str, str], replay: dict[str, str]) -> dict[str, object]:
    target = float(row["lambda1_target"])
    reconstructed = float(replay["lambda1_reconstructed"])
    lambda1_residual = abs(reconstructed - target)
    observable_reproduction = {
        key: abs(float(replay[key]) - float(row[key])) / max(abs(float(row[key])), 1e-300)
        if finite_value(replay, key) and finite_value(row, key) else math.inf
        for key in REPLAY_FIELDS
    }
    max_field = max(observable_reproduction, key=observable_reproduction.get)
    max_difference = observable_reproduction[max_field]
    passed = (
        replay.get("construction_ok") == "1"
        and replay.get("theory_ok") == "1"
        and math.isfinite(lambda1_residual)
        and lambda1_residual <= LAMBDA1_REPLAY_ABS_TOL
        and all(
            math.isfinite(value) and value <= CENTER_REPRO_REL_TOL
            for value in observable_reproduction.values()
        )
    )
    return {
        "lambda1_target": target,
        "lambda1_reconstructed": reconstructed,
        "lambda1_abs_residual": lambda1_residual,
        "lambda1_abs_residual_tolerance": LAMBDA1_REPLAY_ABS_TOL,
        "observable_reproduction": observable_reproduction,
        "maximum_observable_reproduction_relative_difference": max_difference,
        "maximum_observable_reproduction_field": max_field,
        "passed": passed,
    }


def main() -> None:
    row = canonical_row()
    compile_checker()
    center = faithful_m12_2()
    points = [
        ("previous_float", -1, math.nextafter(center, -math.inf)),
        ("center", 0, center),
        ("next_float", 1, math.nextafter(center, math.inf)),
    ]
    probes = []
    for name, ulps, value in points:
        result = evaluate(value)
        result.update(probe=name, ulp_offset=str(ulps), m12_2_gev2=f"{value:.17e}")
        probes.append(result)

    center_result = probes[1]
    m12_sq_construction = center
    m12_sq_roundtrip = float(center_result["m12_sq_reconstructed_gev2"])
    m12_sq_roundtrip_delta = m12_sq_roundtrip - m12_sq_construction
    m12_sq_roundtrip_relative_difference = abs(m12_sq_roundtrip_delta) / abs(m12_sq_construction)
    replay = replay_gate(row, center_result)
    soft_scale_export_status = (
        "VALIDATED_BY_SET_PARAM_PHYS_REPLAY"
        if replay["passed"] else "REJECTED_BY_SET_PARAM_PHYS_REPLAY"
    )
    reproduction = {}
    for key in FIELDS:
        if finite_value(center_result, key) and finite_value(row, key):
            reproduction[key] = abs(float(center_result[key]) - float(row[key])) / max(abs(float(row[key])), 1e-300)
        else:
            reproduction[key] = math.inf
    max_reproduction_field = max(reproduction, key=reproduction.get)
    max_reproduction = reproduction[max_reproduction_field]
    spreads = {
        key: relative_spread([float(probe[key]) for probe in probes])
        if all(finite_value(probe, key) for probe in probes) else math.nan
        for key in FIELDS
    }
    finite_spreads = {key: value for key, value in spreads.items() if math.isfinite(value)}
    max_field = max(finite_spreads, key=finite_spreads.get) if finite_spreads else ""
    max_spread = finite_spreads[max_field] if max_field else math.nan
    all_theory_valid = all(probe.get("theory_ok") == "1" for probe in probes)
    all_values_finite = all(
        finite_value(probe, key) for probe in [*probes, row] for key in FIELDS
    )
    all_relevant_values_physically_valid = all(
        float(probe["total_width_gev"]) > 0.0
        and float(probe["ctau_mm"]) > 0.0
        and all(0.0 <= float(probe[key]) <= 1.0 for key in FIELDS[2:])
        for probe in [*probes, row]
    ) if all_values_finite else False
    center_reproduces = all(math.isfinite(value) and value <= CENTER_REPRO_REL_TOL for value in reproduction.values())
    if not all_theory_valid:
        classification = "NON_USABLE_PROBE_THEORY_INVALID"
    elif not all_values_finite or not all_relevant_values_physically_valid:
        classification = "NONFINITE_OR_PHYSICALLY_UNUSABLE"
    elif not center_reproduces:
        classification = "CENTER_DOES_NOT_REPRODUCE_CANONICAL_ROW"
    elif max_spread > 0.05:
        classification = "NUMERICALLY_UNRESOLVED"
    elif max_spread > 0.02:
        classification = "USABLE_WITH_DECLARED_NUMERICAL_SYSTEMATIC"
    else:
        classification = "STABLE_AT_DOUBLE_REPRESENTATION_SCALE"

    ulp = math.ulp(center)
    slope = sensitivity()
    with OUT_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["m12_2_center_gev2", f"{center:.17e}"])
        writer.writerow(["m12_sq_construction_GeV2", f"{m12_sq_construction:.17e}"])
        writer.writerow(["M2_construction_GeV2", f"{construction_M2(m12_sq_construction):.17e}"])
        writer.writerow(["m12_sq_roundtrip_reconstructed_GeV2", f"{m12_sq_roundtrip:.17e}"])
        writer.writerow(["m12_sq_roundtrip_delta_GeV2", f"{m12_sq_roundtrip_delta:.17e}"])
        writer.writerow(["m12_sq_roundtrip_relative_difference", f"{m12_sq_roundtrip_relative_difference:.17e}"])
        writer.writerow(["soft_scale_export_status", soft_scale_export_status])
        writer.writerow(["lambda1_replay_target", f"{replay['lambda1_target']:.17e}"])
        writer.writerow(["lambda1_replay_reconstructed", f"{replay['lambda1_reconstructed']:.17e}"])
        writer.writerow(["lambda1_replay_abs_residual", f"{replay['lambda1_abs_residual']:.17e}"])
        writer.writerow(["lambda1_replay_abs_residual_tolerance", f"{replay['lambda1_abs_residual_tolerance']:.17e}"])
        writer.writerow(["maximum_replay_observable_relative_difference", f"{replay['maximum_observable_reproduction_relative_difference']:.17e}"])
        writer.writerow(["maximum_replay_observable_field", replay["maximum_observable_reproduction_field"]])
        writer.writerow(["m12_2_ulp_gev2", f"{ulp:.17e}"])
        writer.writerow(["rounding_bound_half_ulp_gev2", f"{0.5 * ulp:.17e}"])
        writer.writerow(["abs_dlambda1_dm12_2_per_gev2", f"{slope:.17e}"])
        writer.writerow(["lambda1_half_ulp_error_bound", f"{slope * 0.5 * ulp:.17e}"])
        writer.writerow(["center_reproduction_tolerance", f"{CENTER_REPRO_REL_TOL:.17e}"])
        writer.writerow(["maximum_center_reproduction_relative_difference", f"{max_reproduction:.17e}"])
        writer.writerow(["maximum_center_reproduction_field", max_reproduction_field])
        writer.writerow(["all_values_finite", "yes" if all_values_finite else "no"])
        writer.writerow(["all_relevant_values_physically_valid", "yes" if all_relevant_values_physically_valid else "no"])
        writer.writerow(["all_adjacent_theory_valid", "yes" if all_theory_valid else "no"])
        writer.writerow(["adjacent_float_maximum_spread", f"{max_spread:.17e}"])
        writer.writerow(["adjacent_float_maximum_spread_field", max_field])
        writer.writerow([])
        writer.writerow(["probe", "ulp_offset", "m12_2_gev2", "theory_ok", "lambda1_reconstructed", *FIELDS])
        for probe in probes:
            writer.writerow([probe["probe"], probe["ulp_offset"], probe["m12_2_gev2"], probe.get("theory_ok"), probe.get("lambda1_reconstructed"), *(probe[key] for key in FIELDS)])
        writer.writerow([])
        writer.writerow(["field", "center_reproduction_relative_difference", "adjacent_float_relative_spread"])
        for key in FIELDS:
            writer.writerow([key, f"{reproduction[key]:.17e}", f"{spreads[key]:.17e}"])
        writer.writerow(["classification", classification])

    lines = [
        f"# Numerical representation check: `{POINT_ID}`", "",
        "This check uses only the exact `m12_2` double and its immediately adjacent representable values.",
        "The previous `1e-12 GeV^2` offsets changed reconstructed `lambda1` by order one and therefore compared different physical models; their channel classifications are withdrawn.", "",
        f"- construction `m12_sq`: `{m12_sq_construction:.17e}` GeV^2",
        f"- construction `M2`: `{construction_M2(m12_sq_construction):.17e}` GeV^2",
        f"- round-trip reconstructed `m12_sq`: `{m12_sq_roundtrip:.17e}` GeV^2",
        f"- round-trip delta: `{m12_sq_roundtrip_delta:.17e}` GeV^2",
        f"- round-trip relative difference: `{m12_sq_roundtrip_relative_difference:.17e}`",
        "- the round-trip reconstructed value is diagnostic only and must not be used as the downstream physical construction input.", "",
        "## Soft-scale export gate", "",
        f"- `soft_scale_export_status`: `{soft_scale_export_status}`",
        f"- replay `lambda1_target`: `{replay['lambda1_target']:.17e}`",
        f"- replay `lambda1_reconstructed`: `{replay['lambda1_reconstructed']:.17e}`",
        f"- replay `lambda1` absolute residual: `{replay['lambda1_abs_residual']:.17e}`",
        f"- declared replay `lambda1` absolute residual tolerance: `{replay['lambda1_abs_residual_tolerance']:.17e}`",
        f"- maximum replay observable relative difference: `{replay['maximum_observable_reproduction_relative_difference']:.17e}` in `{replay['maximum_observable_reproduction_field']}`",
        "- replay gate fields: `total_width_gev`, `ctau_mm`, `br_bb`, `br_gammagamma`, `br_Zgamma`, and `lambda1_reconstructed`.", "",
        f"- center: `{center:.17e}` GeV^2",
        f"- one ULP: `{ulp:.17e}` GeV^2",
        f"- propagated half-ULP `lambda1` bound: `{slope * 0.5 * ulp:.17e}`", "",
        f"- center reproduction tolerance: `{CENTER_REPRO_REL_TOL:.17e}`",
        f"- maximum center reproduction difference: `{max_reproduction:.17e}` in `{max_reproduction_field}`",
        f"- all values finite: `{'yes' if all_values_finite else 'no'}`",
        f"- all relevant values physically valid: `{'yes' if all_relevant_values_physically_valid else 'no'}`",
        f"- adjacent-float maximum spread: `{max_spread:.17e}` in `{max_field}`",
        f"- all adjacent probes theory-valid: `{'yes' if all_theory_valid else 'no'}`", "",
        "| Probe | ULP | theory_ok | lambda1 | total_width_gev | ctau_mm | br_bb | br_gammagamma | br_Zgamma |",
        "|---|---:|:-:|---:|---:|---:|---:|---:|---:|",
    ]
    for probe in probes:
        lines.append(f"| {probe['probe']} | {probe['ulp_offset']} | {probe.get('theory_ok')} | {probe.get('lambda1_reconstructed')} | {probe.get('total_width_gev')} | {probe.get('ctau_mm')} | {probe.get('br_bb')} | {probe.get('br_gammagamma')} | {probe.get('br_Zgamma')} |")
    lines += ["", "| Field | Adjacent-float relative spread |", "|---|---:|"]
    lines += [f"| {key} | {spreads[key]:.6%} |" for key in FIELDS]
    lines += ["", f"## Classification: `{classification}`", "", f"Maximum spread: `{max_spread:.6%}` in `{max_field}`.", "", "This result addresses double-representation uncertainty only; it does not validate scan provenance, production normalization, detector acceptance, or publication readiness."]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({
        "classification": classification,
        "max_spread": max_spread,
        "max_field": max_field,
        "all_adjacent_theory_valid": all_theory_valid,
        "soft_scale_export_status": soft_scale_export_status,
        "m12_sq_construction_GeV2": m12_sq_construction,
        "M2_construction_GeV2": construction_M2(m12_sq_construction),
        "m12_sq_roundtrip_reconstructed_GeV2": m12_sq_roundtrip,
        "lambda1_reconstructed": replay["lambda1_reconstructed"],
        "maximum_replay_observable_relative_difference": replay["maximum_observable_reproduction_relative_difference"],
    }, indent=2))


if __name__ == "__main__":
    main()
