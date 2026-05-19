#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REQUIRED_COLUMNS = [
    "m_phi", "mA", "alpha", "beta", "lambda6", "lambda7", "m12", "sin_ba", "tan_beta",
    "positivity_ok", "unitarity_ok", "perturbativity_ok",
    "width_bb", "width_tautau", "width_WW", "width_ZZ",
    "width_gaga", "width_Zga", "width_gg", "width_hh",
    "total_width", "br_gaga",
    "lam1", "computed_lam1", "lam2", "computed_lam2", "lam3", "lam4", "lam5",
]
PARTIAL_WIDTH_COLUMNS = ["width_bb","width_tautau","width_WW","width_ZZ","width_gaga","width_Zga","width_gg","width_hh"]
OK_COLUMNS = ["positivity_ok","unitarity_ok","perturbativity_ok"]

@dataclass
class FileAudit:
    campaign: str
    file_rel: str
    rows: int
    selected_rows_meta: int | None
    triple_ok_meta: int | None
    triple_ok_csv: int
    missing_required_cols: list[str]
    nonfinite_total_width: int
    negative_total_width: int
    nonfinite_partial_widths: int
    negative_partial_widths: int
    br_missing_when_total_positive: int
    br_out_of_range: int
    br_inconsistent: int
    known_partial_sum_gt_total: int
    lam1_roundtrip_large: int
    lam2_roundtrip_large: int
    meta_event: str
    meta_row_match: bool
    status: str
    notes: list[str]

def parse_args():
    p = argparse.ArgumentParser(description="Audit recomputed campaigns.")
    p.add_argument("--recomputed-root", type=Path, default=Path("/mnt/c/Users/Asus/cern_db/dihiggs_lake/recomputed"))
    p.add_argument("--campaign", action="append", default=[])
    p.add_argument("--ignore-campaign", action="append", default=[])
    p.add_argument("--output-dir", type=Path, default=Path("artifacts/recompute_audit"))
    p.add_argument("--max-files-per-campaign", type=int, default=None)
    p.add_argument("--neg-width-tol", type=float, default=1e-18)
    p.add_argument("--br-upper-eps", type=float, default=1e-8)
    p.add_argument("--rtol-br", type=float, default=1e-8)
    p.add_argument("--atol-br", type=float, default=1e-14)
    p.add_argument("--roundtrip-abs-threshold", type=float, default=1e-8)
    return p.parse_args()

def discover_campaigns(root: Path, campaigns, ignore):
    if campaigns:
        out = [root / c for c in campaigns if c not in ignore]
    else:
        out = sorted([p for p in root.iterdir() if p.is_dir() and p.name.startswith("campaign=") and p.name not in ignore])
    return [p for p in out if p.exists() and p.is_dir()]

def find_scan_csvs(campaign_dir: Path):
    return sorted(campaign_dir.glob("**/tb_*/scan_tb_*.csv"))

def safe_int_from_meta(payload, key):
    if key not in payload:
        return None
    try:
        return int(payload[key])
    except Exception:
        return None

def load_meta(csv_path: Path):
    meta_path = csv_path.parent / "scan_meta.json"
    if not meta_path.exists():
        return {}
    try:
        return json.loads(meta_path.read_text(encoding="utf-8"))
    except Exception:
        return {}

def num(df, col):
    return pd.to_numeric(df[col], errors="coerce")

def compute_file_audit(csv_path: Path, campaign_name: str, args):
    df = pd.read_csv(csv_path)
    meta = load_meta(csv_path)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    notes = []
    if missing:
        return FileAudit(campaign_name, str(csv_path.relative_to(args.recomputed_root)), len(df),
                         safe_int_from_meta(meta, "source_row_count_selected"),
                         safe_int_from_meta(meta, "triple_ok_points"), 0, missing,
                         0,0,0,0,0,0,0,0,0,0,str(meta.get("event","")),False,"BLOCK",["missing required columns"])

    total = num(df, "total_width")
    width_gaga = num(df, "width_gaga")
    br_gaga = num(df, "br_gaga")
    lam1 = num(df, "lam1")
    computed_lam1 = num(df, "computed_lam1")
    lam2 = num(df, "lam2")
    computed_lam2 = num(df, "computed_lam2")
    partials = [num(df, c) for c in PARTIAL_WIDTH_COLUMNS]
    partial_matrix = np.column_stack([s.to_numpy() for s in partials])

    ok_flags = [(num(df, c) == 1.0) for c in OK_COLUMNS]
    triple_ok_csv = int((ok_flags[0] & ok_flags[1] & ok_flags[2]).sum())

    nonfinite_total_width = int((~np.isfinite(total)).sum())
    negative_total_width = int((np.isfinite(total) & (total < -args.neg_width_tol)).sum())
    nonfinite_partial_widths = int((~np.isfinite(partial_matrix)).any(axis=1).sum())
    negative_partial_widths = int((np.isfinite(partial_matrix) & (partial_matrix < -args.neg_width_tol)).any(axis=1).sum())

    total_pos = np.isfinite(total) & (total > 0.0)
    br_missing_when_total_positive = int((total_pos & ~np.isfinite(br_gaga)).sum())
    br_out_of_range = int((np.isfinite(br_gaga) & ((br_gaga < -args.br_upper_eps) | (br_gaga > 1.0 + args.br_upper_eps))).sum())

    br_expected = width_gaga / total
    br_comparable = total_pos & np.isfinite(width_gaga) & np.isfinite(br_gaga) & (width_gaga >= 0.0)
    br_bad_mask = br_comparable & ~np.isclose(br_gaga, br_expected, rtol=args.rtol_br, atol=args.atol_br)
    br_inconsistent = int(br_bad_mask.sum())

    partial_sum = np.nansum(partial_matrix, axis=1)
    known_partial_sum_gt_total = int((np.isfinite(total) & np.isfinite(partial_sum) & (partial_sum > total + max(args.atol_br, 1e-14))).sum())

    lam1_roundtrip_large = int((np.isfinite(lam1) & np.isfinite(computed_lam1) & (np.abs(lam1 - computed_lam1) > args.roundtrip_abs_threshold)).sum())
    lam2_roundtrip_large = int((np.isfinite(lam2) & np.isfinite(computed_lam2) & (np.abs(lam2 - computed_lam2) > args.roundtrip_abs_threshold)).sum())

    selected_rows_meta = safe_int_from_meta(meta, "source_row_count_selected")
    triple_ok_meta = safe_int_from_meta(meta, "triple_ok_points")
    meta_event = str(meta.get("event", ""))
    meta_row_match = (selected_rows_meta is None) or (selected_rows_meta == len(df))

    status = "READY"
    if meta_event not in {"done", "skip_empty_selection", ""}:
        status = "BLOCK"; notes.append(f"unexpected meta event: {meta_event!r}")
    if not meta_row_match:
        status = "BLOCK"; notes.append("meta selected row count mismatch")
    if nonfinite_total_width or negative_total_width or nonfinite_partial_widths or negative_partial_widths:
        status = "BLOCK"; notes.append("non-finite or negative widths")
    if br_missing_when_total_positive or br_out_of_range or br_inconsistent:
        status = "BLOCK"; notes.append("BR integrity failure")
    if known_partial_sum_gt_total and status != "BLOCK":
        status = "SUSPECT"; notes.append("known partial sum > total width")
    if triple_ok_meta is not None and triple_ok_meta != triple_ok_csv and status != "BLOCK":
        status = "SUSPECT"; notes.append("meta triple_ok differs from CSV triple_ok")

    return FileAudit(
        campaign=campaign_name,
        file_rel=str(csv_path.relative_to(args.recomputed_root)),
        rows=len(df),
        selected_rows_meta=selected_rows_meta,
        triple_ok_meta=triple_ok_meta,
        triple_ok_csv=triple_ok_csv,
        missing_required_cols=missing,
        nonfinite_total_width=nonfinite_total_width,
        negative_total_width=negative_total_width,
        nonfinite_partial_widths=nonfinite_partial_widths,
        negative_partial_widths=negative_partial_widths,
        br_missing_when_total_positive=br_missing_when_total_positive,
        br_out_of_range=br_out_of_range,
        br_inconsistent=br_inconsistent,
        known_partial_sum_gt_total=known_partial_sum_gt_total,
        lam1_roundtrip_large=lam1_roundtrip_large,
        lam2_roundtrip_large=lam2_roundtrip_large,
        meta_event=meta_event,
        meta_row_match=meta_row_match,
        status=status,
        notes=notes,
    )

def aggregate_campaign(file_audits):
    if not file_audits:
        return {"campaign_status":"EMPTY","n_files":0}
    campaign_status = "READY"
    if any(f.status == "BLOCK" for f in file_audits):
        campaign_status = "BLOCK"
    elif any(f.status == "SUSPECT" for f in file_audits):
        campaign_status = "SUSPECT"
    return {
        "campaign_status": campaign_status,
        "n_files": len(file_audits),
        "rows_total": int(sum(f.rows for f in file_audits)),
        "triple_ok_csv_total": int(sum(f.triple_ok_csv for f in file_audits)),
        "br_missing_when_total_positive_total": int(sum(f.br_missing_when_total_positive for f in file_audits)),
        "br_out_of_range_total": int(sum(f.br_out_of_range for f in file_audits)),
        "br_inconsistent_total": int(sum(f.br_inconsistent for f in file_audits)),
        "known_partial_sum_gt_total_total": int(sum(f.known_partial_sum_gt_total for f in file_audits)),
        "negative_total_width_total": int(sum(f.negative_total_width for f in file_audits)),
        "negative_partial_widths_total": int(sum(f.negative_partial_widths for f in file_audits)),
        "lam1_roundtrip_large_total": int(sum(f.lam1_roundtrip_large for f in file_audits)),
        "lam2_roundtrip_large_total": int(sum(f.lam2_roundtrip_large for f in file_audits)),
    }

def save_histograms(campaign_out: Path, csv_paths):
    total_vals, br_vals, width_gaga_vals = [], [], []
    for csv_path in csv_paths:
        try:
            df = pd.read_csv(csv_path, usecols=["total_width","br_gaga","width_gaga"])
        except Exception:
            continue
        total = pd.to_numeric(df["total_width"], errors="coerce")
        br = pd.to_numeric(df["br_gaga"], errors="coerce")
        wg = pd.to_numeric(df["width_gaga"], errors="coerce")
        total_vals.append(total[np.isfinite(total) & (total > 0.0)].to_numpy())
        br_vals.append(br[np.isfinite(br)].to_numpy())
        width_gaga_vals.append(wg[np.isfinite(wg) & (wg > 0.0)].to_numpy())
    total_vals = np.concatenate(total_vals) if total_vals else np.array([])
    br_vals = np.concatenate(br_vals) if br_vals else np.array([])
    width_gaga_vals = np.concatenate(width_gaga_vals) if width_gaga_vals else np.array([])

    if total_vals.size:
        plt.figure(figsize=(8,5)); plt.hist(np.log10(total_vals), bins=80)
        plt.xlabel("log10(total_width)"); plt.ylabel("count"); plt.title("Histogram of total_width")
        plt.tight_layout(); plt.savefig(campaign_out / "hist_log10_total_width.png", dpi=160); plt.close()
    if br_vals.size:
        plt.figure(figsize=(8,5)); plt.hist(br_vals, bins=80)
        plt.xlabel("br_gaga"); plt.ylabel("count"); plt.title("Histogram of br_gaga")
        plt.tight_layout(); plt.savefig(campaign_out / "hist_br_gaga.png", dpi=160); plt.close()
    if width_gaga_vals.size:
        plt.figure(figsize=(8,5)); plt.hist(np.log10(width_gaga_vals), bins=80)
        plt.xlabel("log10(width_gaga)"); plt.ylabel("count"); plt.title("Histogram of width_gaga")
        plt.tight_layout(); plt.savefig(campaign_out / "hist_log10_width_gaga.png", dpi=160); plt.close()

def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    campaigns = discover_campaigns(args.recomputed_root, args.campaign, set(args.ignore_campaign))
    all_summary = {}
    for campaign_dir in campaigns:
        csv_paths = find_scan_csvs(campaign_dir)
        if args.max_files_per_campaign is not None:
            csv_paths = csv_paths[:args.max_files_per_campaign]
        file_audits = [compute_file_audit(p, campaign_dir.name, args) for p in csv_paths]
        agg = aggregate_campaign(file_audits)
        campaign_out = args.output_dir / campaign_dir.name
        campaign_out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([asdict(f) for f in file_audits]).to_csv(campaign_out / "file_audit.csv", index=False)
        (campaign_out / "summary.json").write_text(json.dumps(agg, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        save_histograms(campaign_out, csv_paths)
        all_summary[campaign_dir.name] = agg
        print(f"{campaign_dir.name}: {agg['campaign_status']} | files={agg['n_files']} rows={agg.get('rows_total',0)} "
              f"br_missing={agg.get('br_missing_when_total_positive_total',0)} "
              f"br_inconsistent={agg.get('br_inconsistent_total',0)}")
    (args.output_dir / "all_campaigns_summary.json").write_text(json.dumps(all_summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    rows = []
    for name, agg in all_summary.items():
        row = {"campaign": name}; row.update(agg); rows.append(row)
    if rows:
        pd.DataFrame(rows).sort_values(["campaign_status", "campaign"]).to_csv(args.output_dir / "all_campaigns_summary.csv", index=False)
    print(f"\nWrote audit artifacts to: {args.output_dir}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
