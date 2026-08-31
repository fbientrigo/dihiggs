#!/usr/bin/env python3
"""Basis check: for a handful of accepted trajectory points, re-express the
same physical point via the second canonical evaluator (Lambda1EvaluatorV2,
lambda1-input basis) and confirm physical quantities agree with the primary
M2-input evaluator (DihiggsPointV2Evaluator) to numerical precision.

Does not modify either evaluator. Read-only comparison.
"""
import csv
import math
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
BINARY_B = REPO_ROOT / "dihiggs" / "app" / "Lambda1EvaluatorV2"

sys.path.insert(0, str(REPO_ROOT / "scripts"))
import h2_m2_continuation as m  # noqa: E402


def run_lambda1(row, tag, workdir):
    in_csv = workdir / f"basis_{tag}_in.csv"
    out_csv = workdir / f"basis_{tag}_out.csv"
    # Lambda1EvaluatorV2 compares the header line via std::getline with no
    # CR-stripping, so this must be written with bare "\n" -- csv.writer's
    # default "\r\n" line terminator makes the header comparison fail.
    with open(in_csv, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("point_id,mh_gev,mH_gev,mA_gev,mHp_gev,sin_beta_minus_alpha,"
                  "tan_beta,lambda1_target,lambda6_input,lambda7_input\n")
        fh.write(",".join([tag, row["mh_GeV"], row["m_H2_GeV"], row["m_A_GeV"],
                            row["m_Hp_GeV"], row["sin_beta_minus_alpha"],
                            row["tan_beta"], row["lambda1"], row["lambda6"],
                            row["lambda7"]]) + "\n")
    subprocess.run([str(BINARY_B), str(in_csv), str(out_csv)], check=True,
                    capture_output=True, text=True, cwd=REPO_ROOT)
    with open(out_csv, newline="", encoding="utf-8") as fh:
        return next(csv.DictReader(fh))


def reldiff(a, b):
    a, b = float(a), float(b)
    scale = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / scale


def check_point(row, tag, workdir):
    b = run_lambda1(row, tag, workdir)
    beta = math.atan(float(row["tan_beta"]))
    denom = math.sin(beta) * math.cos(beta)
    M2_from_B = float(b["m12_sq_reconstructed_gev2"]) / denom

    comparisons = {
        "m12_sq_GeV2": reldiff(row["m12_sq_GeV2"], b["m12_sq_reconstructed_gev2"]),
        "M2_GeV2": reldiff(row["M2_GeV2"], M2_from_B),
        "tan_beta_reconstructed": reldiff(row["tan_beta"], b["tan_beta_reconstructed"]),
        "total_width_GeV": reldiff(row["total_width_GeV"], b["total_width_gev"]),
        "ctau_mm": reldiff(row["ctau_physical_mm"], b["ctau_mm"]),
        "br_bb": reldiff(row["BR_bb"], b["br_bb"]),
        "br_tautau": reldiff(row["BR_tautau"], b["br_tautau"]),
        "br_WW": reldiff(row["BR_WW"], b["br_WW"]),
        "br_ZZ": reldiff(row["BR_ZZ"], b["br_ZZ"]),
        "br_gg": reldiff(row["BR_gg"], b["br_gg"]),
        "br_gammagamma": reldiff(row["BR_gammagamma"], b["br_gammagamma"]),
        "br_hh": reldiff(row["BR_hh"], b["br_hh"]),
    }
    theory_ok_b = b["theory_ok"] == "1"
    return {
        "point_tag": tag,
        "m_H2_GeV": row["m_H2_GeV"],
        "theory_ok_A": row["theory_ok"] == m.VALID_TOKEN,
        "theory_ok_B": theory_ok_b,
        **{f"reldiff_{k}": f"{v:.3e}" for k, v in comparisons.items()},
        "max_reldiff": f"{max(comparisons.values()):.3e}",
        "pass_1e-8": max(comparisons.values()) < 1e-8,
    }


def main():
    slices = [
        ("mA300", REPO_ROOT / "results" / "h2_m2_continuation" / "trajectory_mA300.csv"),
        ("mA350", REPO_ROOT / "results" / "h2_m2_continuation" / "trajectory_mA350.csv"),
    ]
    workdir = REPO_ROOT / "results" / "h2_m2_continuation" / "_scratch"
    workdir.mkdir(parents=True, exist_ok=True)
    out_rows = []
    for slice_name, path in slices:
        if not path.exists():
            print(f"skip {slice_name}: {path} not found")
            continue
        rows = list(csv.DictReader(open(path)))
        picks = [("first", rows[0]), ("middle", rows[len(rows) // 2]),
                 ("highest", rows[-1])]
        for label, row in picks:
            tag = f"{slice_name}_{label}"
            result = check_point(row, tag, workdir)
            result["slice"] = slice_name
            result["label"] = label
            out_rows.append(result)
            print(tag, "mH2=", row["m_H2_GeV"], "max_reldiff=",
                  result["max_reldiff"], "pass=", result["pass_1e-8"])

    if out_rows:
        out_csv = REPO_ROOT / "results" / "h2_m2_continuation" / "basis_check.csv"
        fieldnames = ["slice", "label"] + [k for k in out_rows[0].keys()
                                            if k not in ("slice", "label")]
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            for r in out_rows:
                w.writerow(r)
        print("wrote", out_csv)


if __name__ == "__main__":
    main()
