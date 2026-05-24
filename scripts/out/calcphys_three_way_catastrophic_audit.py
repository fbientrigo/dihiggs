#!/usr/bin/env python3
import argparse
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

BRIEF_POINTS = [
    ("A", 128.006913407821),
    ("B", 128.453840782123),
    ("C", 129.258310055866),
    ("D", 129.973393854749),
    ("E", 131.046019553073),
]


def parse_decay35(path: Path):
    lines = path.read_text().splitlines()
    idx = None
    for i, ln in enumerate(lines):
        if ln.startswith("DECAY  35"):
            idx = i
            break
    if idx is None:
        return None
    tot = float(lines[idx].split()[2])
    br = {}
    j = idx + 1
    while j < len(lines) and not lines[j].startswith("DECAY"):
        s = lines[j].strip()
        if s and not s.startswith("#"):
            p = s.split()
            if len(p) >= 4:
                b = float(p[0])
                i1 = int(p[2])
                i2 = int(p[3])
                br[tuple(sorted((i1, i2)))] = b
        j += 1

    def pw(key):
        return tot * br.get(key, 0.0)

    return {
        "total_width": tot,
        "width_gaga": pw((22, 22)),
        "width_Zga": pw((22, 23)),
        "width_bb": pw((-5, 5)),
        "width_tautau": pw((-15, 15)),
        "width_gg": pw((21, 21)),
        "br_gaga": br.get((22, 22), math.nan),
    }


def safe_ratio(a: float, b: float) -> float:
    if b == 0:
        return math.inf
    return a / b


def abs_log10_ratio(r: float) -> float:
    if r > 0 and math.isfinite(r):
        return abs(math.log10(r))
    return math.inf


def run_cmd(cmd, env=None):
    cp = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    if cp.returncode != 0:
        raise RuntimeError(f"Command failed ({cp.returncode}): {' '.join(cmd)}\nSTDOUT:\n{cp.stdout}\nSTDERR:\n{cp.stderr}")
    return cp


def _format_md_cell(v) -> str:
    if pd.isna(v):
        return "nan"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return f"{v:.17g}"
    return str(v)


def dataframe_to_markdown_resilient(df: pd.DataFrame) -> str:
    try:
        return df.to_markdown(index=False)
    except Exception as exc:
        print(
            f"[warn] pandas.to_markdown unavailable; falling back to manual markdown table ({type(exc).__name__}: {exc})",
            file=sys.stderr,
        )

    headers = [str(c) for c in df.columns]
    head = "| " + " | ".join(headers) + " |"
    sep = "| " + " | ".join(["---"] * len(headers)) + " |"
    body = []
    for _, row in df.iterrows():
        vals = [_format_md_cell(row[c]).replace("|", "\\|") for c in df.columns]
        body.append("| " + " | ".join(vals) + " |")
    return "\n".join([head, sep] + body)


def run_calcphys(calcphys: Path, row: pd.Series, m12_value: float, out_file: Path):
    cmd = [
        str(calcphys),
        "125",
        f"{float(row['m_phi']):.17g}",
        f"{float(row['mA']):.17g}",
        f"{float(row['mA']):.17g}",
        f"{float(row['sin_ba']):.17g}",
        f"{float(row['lambda6']):.17g}",
        f"{float(row['lambda7']):.17g}",
        f"{float(m12_value):.17g}",
        f"{float(row['tan_beta']):.17g}",
        "1",
        str(out_file),
    ]
    run_cmd(cmd)
    parsed = parse_decay35(out_file)
    if parsed is None:
        raise RuntimeError(f"Failed parsing DECAY 35 from {out_file}")
    return parsed


def select_points(parity_df: pd.DataFrame, mode: str) -> pd.DataFrame:
    targets = BRIEF_POINTS if mode == "A_to_E" else BRIEF_POINTS[:1]
    rows = []
    for label, mphi in targets:
        subset = parity_df[(parity_df["m_phi"].sub(mphi).abs() < 1e-9)]
        if subset.empty:
            raise RuntimeError(f"No old parity rows found for {label} m_phi={mphi}")
        # prefer lam1=1e-12 representative, fallback nearest
        picked = subset.iloc[(subset["lam1"].sub(1e-12).abs()).argsort()[:1]].copy()
        picked.insert(0, "point_id", label)
        rows.append(picked)
    out = pd.concat(rows, ignore_index=True)
    return out


def recompute_new(scan_exe: Path, out_csv: Path, row: pd.Series, env: dict):
    cmd = [
        str(scan_exe),
        f"{float(row['m_phi']):.15g}", f"{float(row['m_phi']):.15g}", "1",
        f"{float(row['lam1']):.15g}", f"{float(row['lam1']):.15g}", "1",
        f"{float(row['mA']):.15g}", f"{float(row['sin_ba']):.15g}", f"{float(row['tan_beta']):.15g}",
        f"{float(row['lambda6']):.15g}", f"{float(row['lambda7']):.15g}",
        str(out_csv),
    ]
    run_cmd(cmd, env=env)
    df = pd.read_csv(out_csv)
    if df.empty:
        raise RuntimeError("PhysScanWithFixings emitted empty CSV")
    required = ["m12_2_used", "m12_2_gen_after_set", "delta_m12_2_gen_minus_used", "replay_semantics_version", "total_width", "width_gaga", "width_Zga"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"Schema failure: recomputed CSV missing required columns: {missing}")
    return df.iloc[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["pointA", "A_to_E"], required=True)
    ap.add_argument("--parity-csv", default="scripts/out/autoresearch_calcphys_eda_full/seed_000/parity_g180/calcphys_parity_rows.csv")
    ap.add_argument("--old-scan-csv", default="scripts/out/autoresearch_calcphys_eda_full/seed_000/scan_g180.csv")
    ap.add_argument("--scan-exe", default="./dihiggs/app/PhysScanWithFixings")
    ap.add_argument("--calcphys", default="./2hdmc/CalcPhys")
    ap.add_argument("--outdir", default="scripts/out/three_way_catastrophic_audit")
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    parity_df = pd.read_csv(args.parity_csv)
    old_scan_df = pd.read_csv(args.old_scan_csv)

    mode_norm = "A_to_E" if args.mode == "A_to_E" else "pointA"
    selected = select_points(parity_df, "A_to_E" if mode_norm == "A_to_E" else "pointA")

    # Attach full old-scan parameters via idx
    if "idx" not in selected.columns:
        raise RuntimeError("old parity CSV has no idx column")
    old_scan_df = old_scan_df.reset_index().rename(columns={"index": "idx"})
    merged = selected.merge(old_scan_df, on="idx", how="left", suffixes=("_oldparity", "_oldscan"))
    if merged[["mA", "tan_beta", "sin_ba", "lambda6", "lambda7"]].isna().any().any():
        raise RuntimeError("Failed to recover full replay inputs from old scan rows")

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"

    rows = []
    run_cmd_lines = []
    t0 = time.time()
    for _, r in merged.iterrows():
        point_id = str(r["point_id"])
        one_csv = outdir / f"recomputed_{point_id}.csv"
        decay_out_old = outdir / f"calcphys_oldm12_{point_id}.out"
        decay_out_new = outdir / f"calcphys_newm12_{point_id}.out"

        run_cmd_lines.append(
            " ".join([
                str(Path(args.scan_exe).resolve()),
                f"{float(r['m_phi_oldparity']):.15g}", f"{float(r['m_phi_oldparity']):.15g}", "1",
                f"{float(r['lam1_oldparity']):.15g}", f"{float(r['lam1_oldparity']):.15g}", "1",
                f"{float(r['mA']):.15g}", f"{float(r['sin_ba']):.15g}", f"{float(r['tan_beta']):.15g}",
                f"{float(r['lambda6']):.15g}", f"{float(r['lambda7']):.15g}",
                str(one_csv),
            ])
        )

        new_row = recompute_new(Path(args.scan_exe).resolve(), one_csv, pd.Series({
            "m_phi": r["m_phi_oldparity"], "lam1": r["lam1_oldparity"], "mA": r["mA"], "sin_ba": r["sin_ba"],
            "tan_beta": r["tan_beta"], "lambda6": r["lambda6"], "lambda7": r["lambda7"]
        }), env)

        calc_in = pd.Series({
            "m_phi": r["m_phi_oldparity"], "mA": r["mA"], "sin_ba": r["sin_ba"], "tan_beta": r["tan_beta"],
            "lambda6": r["lambda6"], "lambda7": r["lambda7"]
        })
        gt_old = run_calcphys(Path(args.calcphys).resolve(), calc_in, float(r["m12_oldparity"]), decay_out_old)
        gt_new = run_calcphys(Path(args.calcphys).resolve(), calc_in, float(new_row["m12_2_used"]), decay_out_new)

        old_total = float(r["total_width_scan"])
        new_total = float(new_row["total_width"])
        gt_old_total = float(gt_old["total_width"])
        gt_new_total = float(gt_new["total_width"])
        r_old_gt = safe_ratio(old_total, gt_old_total)
        r_new_gt = safe_ratio(new_total, gt_new_total)
        r_old_new = safe_ratio(old_total, new_total)

        old_gaga = float(r.get("width_gaga_scan", math.nan))
        new_gaga = float(new_row.get("width_gaga", math.nan))
        gt_old_gaga = float(gt_old.get("width_gaga", math.nan))
        gt_new_gaga = float(gt_new.get("width_gaga", math.nan))

        old_zga = float(r.get("width_Zga_scan", math.nan))
        new_zga = float(new_row.get("width_Zga", math.nan))
        gt_old_zga = float(gt_old.get("width_Zga", math.nan))
        gt_new_zga = float(gt_new.get("width_Zga", math.nan))

        row_out = {
            "point_id": point_id,
            "old_source_file": str(Path(args.parity_csv).resolve()),
            "old_row_index": int(r["idx"]),
            "m_phi": float(r["m_phi_oldparity"]),
            "lam1": float(r["lam1_oldparity"]),
            "tan_beta": float(r["tan_beta"]),
            "sin_ba": float(r["sin_ba"]),
            "lambda6": float(r["lambda6"]),
            "lambda7": float(r["lambda7"]),
            "mA": float(r["mA"]),
            "old_m12_or_m12_2": float(r["m12_oldparity"]),
            "new_m12_2_used": float(new_row["m12_2_used"]),
            "new_m12_2_gen_after_set": float(new_row["m12_2_gen_after_set"]),
            "delta_m12_2_gen_minus_used": float(new_row["delta_m12_2_gen_minus_used"]),
            "replay_semantics_version": str(new_row["replay_semantics_version"]),
            "old_total_width": old_total,
            "new_total_width": new_total,
            "calcphys_total_width": gt_new_total,
            "calcphys_total_width_old_m12": gt_old_total,
            "calcphys_total_width_new_m12_2_used": gt_new_total,
            "R_old_vs_calcphys_total_width": r_old_gt,
            "R_new_vs_calcphys_total_width": r_new_gt,
            "R_old_vs_new_total_width": r_old_new,
            "log10_abs_old_vs_calcphys": abs_log10_ratio(r_old_gt),
            "log10_abs_new_vs_calcphys": abs_log10_ratio(r_new_gt),
            "catastrophic_old_vs_calcphys": abs_log10_ratio(r_old_gt) >= 3.0,
            "catastrophic_new_vs_calcphys": abs_log10_ratio(r_new_gt) >= 3.0,
            "improvement_orders_total_width": abs_log10_ratio(r_old_gt) - abs_log10_ratio(r_new_gt),
            "old_width_gaga": old_gaga,
            "new_width_gaga": new_gaga,
            "calcphys_width_gaga": gt_new_gaga,
            "calcphys_width_gaga_old_m12": gt_old_gaga,
            "calcphys_width_gaga_new_m12_2_used": gt_new_gaga,
            "R_old_vs_calcphys_gaga": safe_ratio(old_gaga, gt_old_gaga),
            "R_new_vs_calcphys_gaga": safe_ratio(new_gaga, gt_new_gaga),
            "log10_abs_old_vs_calcphys_gaga": abs_log10_ratio(safe_ratio(old_gaga, gt_old_gaga)),
            "log10_abs_new_vs_calcphys_gaga": abs_log10_ratio(safe_ratio(new_gaga, gt_new_gaga)),
            "old_width_Zga": old_zga,
            "new_width_Zga": new_zga,
            "calcphys_width_Zga": gt_new_zga,
            "calcphys_width_Zga_old_m12": gt_old_zga,
            "calcphys_width_Zga_new_m12_2_used": gt_new_zga,
            "R_old_vs_calcphys_Zga": safe_ratio(old_zga, gt_old_zga),
            "R_new_vs_calcphys_Zga": safe_ratio(new_zga, gt_new_zga),
            "log10_abs_old_vs_calcphys_Zga": abs_log10_ratio(safe_ratio(old_zga, gt_old_zga)),
            "log10_abs_new_vs_calcphys_Zga": abs_log10_ratio(safe_ratio(new_zga, gt_new_zga)),
        }
        rows.append(row_out)

    out_df = pd.DataFrame(rows)
    csv_name = "A_to_E_three_way.csv" if mode_norm == "A_to_E" else "pointA_three_way.csv"
    md_name = "A_to_E_three_way.md" if mode_norm == "A_to_E" else "pointA_three_way.md"
    out_csv = outdir / csv_name
    out_md = outdir / md_name
    out_df.to_csv(out_csv, index=False)

    lines = [
        f"# Three-way catastrophic audit ({mode_norm})",
        "",
        "Comparison definition:",
        "1) old pre-fix parity row values vs CalcPhys",
        "2) newly recomputed PhysScanWithFixings (current fixed code path) vs CalcPhys",
        "3) old pre-fix values vs newly recomputed values",
        "",
        f"Method for new recomputation: direct PhysScanWithFixings single-point invocation (m_phi and lam1 one-point grid) with recovered old replay inputs mA/sin_ba/tan_beta/lambda6/lambda7; OMP_NUM_THREADS=1.",
        "",
        f"Rows audited: {len(out_df)}",
        f"Elapsed seconds: {time.time()-t0:.3f}",
        "",
        dataframe_to_markdown_resilient(out_df[[
            "point_id", "old_row_index", "m_phi", "lam1", "old_total_width", "new_total_width", "calcphys_total_width",
            "R_old_vs_calcphys_total_width", "R_new_vs_calcphys_total_width", "R_old_vs_new_total_width",
            "log10_abs_old_vs_calcphys", "log10_abs_new_vs_calcphys",
            "catastrophic_old_vs_calcphys", "catastrophic_new_vs_calcphys",
            "new_m12_2_used", "new_m12_2_gen_after_set", "delta_m12_2_gen_minus_used", "replay_semantics_version"
        ]]),
        "",
    ]
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    state = {
        "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "mode": mode_norm,
        "rows_audited": int(len(out_df)),
        "catastrophic_old_vs_calcphys": int(out_df["catastrophic_old_vs_calcphys"].sum()),
        "catastrophic_new_vs_calcphys": int(out_df["catastrophic_new_vs_calcphys"].sum()),
        "max_log10_abs_old_vs_calcphys": float(out_df["log10_abs_old_vs_calcphys"].max()),
        "max_log10_abs_new_vs_calcphys": float(out_df["log10_abs_new_vs_calcphys"].max()),
        "method_new_recompute": "PhysScanWithFixings single-point direct invocation",
        "inputs": {
            "parity_csv": str(Path(args.parity_csv).resolve()),
            "old_scan_csv": str(Path(args.old_scan_csv).resolve()),
        },
        "outputs": {
            "csv": str(out_csv),
            "md": str(out_md),
        },
    }
    (outdir / "state.json").write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")

    run_sh = outdir / "run_commands.sh"
    run_sh.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        "python3 scripts/out/calcphys_three_way_catastrophic_audit.py --mode pointA\n"
        "python3 scripts/out/calcphys_three_way_catastrophic_audit.py --mode A_to_E\n",
        encoding="utf-8",
    )
    os.chmod(run_sh, 0o755)

    print(f"Wrote: {out_csv}")
    print(f"Wrote: {out_md}")
    print(f"Wrote: {outdir / 'state.json'}")
    print(f"Wrote: {run_sh}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
