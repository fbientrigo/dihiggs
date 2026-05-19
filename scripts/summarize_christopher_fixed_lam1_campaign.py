#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path
from typing import Dict, List

DEFAULT_CAMPAIGN = "christopher_fixed_lam1_2026apr"
DEFAULT_LAKE_ROOT = Path("/mnt/c/Users/Asus/cern_db/dihiggs_lake")
EPS = 1e-9

REQUIRED_COLUMNS = {
    "m_phi",
    "lam1",
    "lambda6",
    "tan_beta",
    "positivity_ok",
    "unitarity_ok",
    "perturbativity_ok",
    "total_width",
    "br_gaga",
}


def _f(v: str) -> float:
    return float(v)


def _b(v: str) -> bool:
    vv = str(v).strip().lower()
    if vv in {"1", "1.0", "true", "t", "yes", "y"}:
        return True
    if vv in {"0", "0.0", "false", "f", "no", "n"}:
        return False
    return bool(float(vv))


def _constant(vals: List[float], tol: float = EPS) -> bool:
    if not vals:
        return False
    a = vals[0]
    return all(abs(x - a) <= tol for x in vals)


def _find_meta_for_csv(csv_path: Path) -> Dict:
    meta_path = csv_path.parent / "scan_meta.json"
    if not meta_path.exists():
        return {}
    try:
        return json.loads(meta_path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def summarize_csv(csv_path: Path, expected_n_mphi: int) -> Dict:
    row: Dict = {
        "csv_path": str(csv_path),
        "lambda1": "",
        "lambda6": "",
        "tan_beta": "",
        "n_rows": 0,
        "n_triple_ok": 0,
        "fraction_triple_ok": "",
        "min_m_phi": "",
        "max_m_phi": "",
        "median_total_width": "",
        "min_total_width": "",
        "max_total_width": "",
        "median_br_gaga": "",
        "min_br_gaga": "",
        "max_br_gaga": "",
        "status": "ERROR",
    }
    try:
        with csv_path.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                row["status"] = "EMPTY"
                return row
            missing = REQUIRED_COLUMNS - set(reader.fieldnames)
            if missing:
                row["status"] = "MISSING_COLUMNS"
                row["missing_columns"] = ",".join(sorted(missing))
                return row
            rows = list(reader)

        if not rows:
            row["status"] = "EMPTY"
            return row

        lam1_vals = [_f(r["lam1"]) for r in rows]
        l6_vals = [_f(r["lambda6"]) for r in rows]
        tb_vals = [_f(r["tan_beta"]) for r in rows]
        mphi_vals = [_f(r["m_phi"]) for r in rows]
        total_width_vals = [_f(r["total_width"]) for r in rows]
        br_gaga_vals = [_f(r["br_gaga"]) for r in rows]
        triple_ok = [
            _b(r["positivity_ok"]) and _b(r["unitarity_ok"]) and _b(r["perturbativity_ok"])
            for r in rows
        ]

        meta = _find_meta_for_csv(csv_path)
        row["n_rows"] = len(rows)
        row["n_triple_ok"] = int(sum(1 for x in triple_ok if x))
        row["fraction_triple_ok"] = row["n_triple_ok"] / len(rows)

        row["lambda1"] = lam1_vals[0] if _constant(lam1_vals) else (meta.get("grid", {}).get("lam1_min", ""))
        row["lambda6"] = l6_vals[0] if _constant(l6_vals) else (meta.get("fixed_params", {}).get("lambda6", ""))
        row["tan_beta"] = tb_vals[0] if _constant(tb_vals) else meta.get("tanbeta", "")

        row["min_m_phi"] = min(mphi_vals)
        row["max_m_phi"] = max(mphi_vals)

        row["median_total_width"] = statistics.median(total_width_vals)
        row["min_total_width"] = min(total_width_vals)
        row["max_total_width"] = max(total_width_vals)

        row["median_br_gaga"] = statistics.median(br_gaga_vals)
        row["min_br_gaga"] = min(br_gaga_vals)
        row["max_br_gaga"] = max(br_gaga_vals)

        constants_ok = _constant(lam1_vals) and _constant(l6_vals) and _constant(tb_vals)
        span_ok = abs(row["min_m_phi"] - min(mphi_vals)) <= 1e-9 and abs(row["max_m_phi"] - max(mphi_vals)) <= 1e-9

        if not constants_ok or not span_ok:
            row["status"] = "ERROR"
        elif len(rows) < expected_n_mphi:
            row["status"] = "LOW_YIELD"
        elif len(rows) == expected_n_mphi:
            row["status"] = "OK"
        else:
            row["status"] = "PARTIAL"
        return row
    except Exception:
        row["status"] = "ERROR"
        return row


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Summarize fixed-lambda1 Christopher campaign outputs")
    p.add_argument("--campaign", type=str, default=DEFAULT_CAMPAIGN)
    p.add_argument("--lake-root", type=str, default=str(DEFAULT_LAKE_ROOT))
    p.add_argument("--n-mphi", type=int, default=400)
    return p


def main() -> int:
    args = build_parser().parse_args()
    campaign = args.campaign
    lake_campaign = Path(args.lake_root) / f"campaign={campaign}"

    repo_root = Path(__file__).resolve().parents[1]
    outdir = repo_root / "scripts" / "out" / campaign
    outdir.mkdir(parents=True, exist_ok=True)
    out_csv = outdir / "summary.csv"
    out_md = outdir / "summary.md"

    csv_paths = sorted(lake_campaign.glob("fixed_*/*/tb_*/scan_tb_*.csv"))
    rows = [summarize_csv(p, expected_n_mphi=args.n_mphi) for p in csv_paths]

    fieldnames = [
        "csv_path",
        "lambda1",
        "lambda6",
        "tan_beta",
        "n_rows",
        "n_triple_ok",
        "fraction_triple_ok",
        "min_m_phi",
        "max_m_phi",
        "median_total_width",
        "min_total_width",
        "max_total_width",
        "median_br_gaga",
        "min_br_gaga",
        "max_br_gaga",
        "status",
    ]

    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})

    status_counts: Dict[str, int] = {}
    for r in rows:
        status_counts[r["status"]] = status_counts.get(r["status"], 0) + 1

    md_lines = [
        f"# Summary for campaign={campaign}",
        "",
        f"- Lake campaign path: `{lake_campaign}`",
        f"- CSV files found: {len(rows)}",
        f"- Expected n_mphi per curve: {args.n_mphi}",
        "",
        "## Status counts",
    ]
    for k in sorted(status_counts):
        md_lines.append(f"- {k}: {status_counts[k]}")

    md_lines.extend([
        "",
        "## Notes",
        "- LOW_YIELD means fewer than requested m_phi rows (accepted for internal filtering behavior).",
        "- MISSING_COLUMNS means required contract columns were absent.",
    ])

    out_md.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(f"Wrote {out_csv}")
    print(f"Wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
