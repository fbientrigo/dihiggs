#!/usr/bin/env python3
from __future__ import annotations

import argparse
from bisect import bisect_left
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

DEFAULT_MERGED_CSV = Path("scripts/out/refactor_colab_compare/real_run/colab_vs_2hdmc_merged_comparison.csv")
DEFAULT_FALLBACK_MERGED_CSV = Path("scripts/out/current_baseline/colab_vs_2hdmc_merged_comparison.csv")
DEFAULT_CAMPAIGN_ROOT = Path("/mnt/c/Users/Asus/cern_db/dihiggs_lake/campaign=christopher_fixed_lam1_2026apr")
DEFAULT_SUMMARY_CSV = Path("scripts/out/christopher_fixed_lam1_2026apr/summary.csv")
DEFAULT_OUT_DIR = Path("scripts/out/christopher_fixed_lam1_2026apr/comparison")

PARAM_KEY_COLS = ["lambda1", "lambda6", "tan_beta"]
OPTIONAL_CHANNELS = ["width_cc"]
CHANNEL_COLS = [
    "width_bb",
    "width_cc",
    "width_tautau",
    "width_gg",
    "width_gaga",
    "width_Zga",
    "total_width",
    "br_gaga",
]


def _to_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def _safe_rel_diff(a: pd.Series, b: pd.Series, eps: float = 1e-30) -> pd.Series:
    a = _to_num(a)
    b = _to_num(b)
    d = b.abs().clip(lower=eps)
    out = (a - b).abs() / d
    return out.where(b.abs() > eps, np.nan)


def _choose_first(df: pd.DataFrame, candidates: Iterable[str], required: bool = True) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise ValueError(f"Missing required column candidates: {list(candidates)}")
    return None


def _normalize_merged(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    c_mphi = _choose_first(out, ["m_phi", "mH", "mh"])
    c_l1 = _choose_first(out, ["lambda1", "lam1", "lam1_g", "lambda1_input"])
    c_l6 = _choose_first(out, ["lambda6"])
    c_tb = _choose_first(out, ["tan_beta", "tb"])

    rename = {
        c_mphi: "m_phi",
        c_l1: "lambda1",
        c_l6: "lambda6",
        c_tb: "tan_beta",
    }
    out = out.rename(columns=rename)

    for c in ["m_phi", "lambda1", "lambda6", "tan_beta"]:
        out[c] = _to_num(out[c])

    if "triple_ok" in out.columns:
        out["triple_ok_colab_recomputed"] = out["triple_ok"].astype(bool)
    else:
        p = _choose_first(out, ["positivity_ok"], required=False)
        u = _choose_first(out, ["unitarity_ok"], required=False)
        t = _choose_first(out, ["perturbativity_ok"], required=False)
        if p and u and t:
            out["triple_ok_colab_recomputed"] = out[p].astype(bool) & out[u].astype(bool) & out[t].astype(bool)
        else:
            raise ValueError("Merged CSV has neither triple_ok nor full positivity/unitarity/perturbativity columns")

    map_cols = {
        "width_bb": ["width_bb"],
        "width_cc": ["width_cc"],
        "width_tautau": ["width_tautau"],
        "width_gg": ["width_gg"],
        "width_gaga": ["width_gaga"],
        "width_Zga": ["width_Zga"],
        "total_width": ["total_width"],
        "br_gaga": ["br_gaga_2hdmc", "br_gaga"],
    }
    for k, cands in map_cols.items():
        c = _choose_first(out, cands, required=False)
        if c is not None and c != k:
            out = out.rename(columns={c: k})

    return out


def _load_fixed_campaign(campaign_root: Path) -> pd.DataFrame:
    files = sorted(campaign_root.glob("**/scan_tb_*.csv"))
    if not files:
        raise FileNotFoundError(f"No scan_tb_*.csv files found under {campaign_root}")
    parts = []
    for f in files:
        d = pd.read_csv(f)
        d["__source_file"] = str(f)
        parts.append(d)
    out = pd.concat(parts, ignore_index=True)

    c_mphi = _choose_first(out, ["m_phi", "mH", "mh"])
    c_l1 = _choose_first(out, ["lambda1", "lam1"])
    c_l6 = _choose_first(out, ["lambda6"])
    c_tb = _choose_first(out, ["tan_beta", "tb"])

    out = out.rename(columns={c_mphi: "m_phi", c_l1: "lambda1", c_l6: "lambda6", c_tb: "tan_beta"})
    for c in ["m_phi", "lambda1", "lambda6", "tan_beta"]:
        out[c] = _to_num(out[c])

    p = _choose_first(out, ["positivity_ok"])
    u = _choose_first(out, ["unitarity_ok"])
    t = _choose_first(out, ["perturbativity_ok"])
    out["triple_ok_fixed_campaign"] = out[p].astype(bool) & out[u].astype(bool) & out[t].astype(bool)
    return out


def _group_key(df: pd.DataFrame, param_tol: float) -> pd.Series:
    vals = [np.round(_to_num(df[c]) / param_tol).astype("Int64").astype(str) for c in PARAM_KEY_COLS]
    return vals[0] + "|" + vals[1] + "|" + vals[2]


def match_points(colab_df: pd.DataFrame, fixed_df: pd.DataFrame, mphi_tol: float, param_tol: float) -> pd.DataFrame:
    c = colab_df.copy().reset_index(drop=True)
    f = fixed_df.copy().reset_index(drop=True)
    c["_key"] = _group_key(c, param_tol)
    f["_key"] = _group_key(f, param_tol)

    rows = []
    for k in sorted(set(c["_key"]).intersection(set(f["_key"]))):
        cg = c[c["_key"] == k].copy().sort_values("m_phi").reset_index(drop=True)
        fg = f[f["_key"] == k].copy().sort_values("m_phi").reset_index(drop=True)
        mvals = fg["m_phi"].to_numpy()
        used = np.zeros(len(fg), dtype=bool)

        for _, crow in cg.iterrows():
            x = float(crow["m_phi"])
            i = bisect_left(mvals, x)
            cand = []
            if i < len(mvals):
                cand.append(i)
            if i - 1 >= 0:
                cand.append(i - 1)
            best_j = None
            best_d = None
            for j in cand:
                if used[j]:
                    continue
                d = abs(float(mvals[j]) - x)
                if d <= mphi_tol and (best_d is None or d < best_d):
                    best_j = j
                    best_d = d
            if best_j is None:
                continue
            used[best_j] = True
            frow = fg.iloc[best_j]

            row = {
                "lambda1": crow["lambda1"],
                "lambda6": crow["lambda6"],
                "tan_beta": crow["tan_beta"],
                "m_phi_colab": crow["m_phi"],
                "m_phi_fixed": frow["m_phi"],
                "m_phi_abs_diff": abs(float(crow["m_phi"]) - float(frow["m_phi"])),
                "triple_ok_colab_recomputed": bool(crow["triple_ok_colab_recomputed"]),
                "triple_ok_fixed_campaign": bool(frow["triple_ok_fixed_campaign"]),
            }
            for ch in CHANNEL_COLS:
                ccol = ch
                fcol = ch
                row[f"{ch}_colab"] = crow[ccol] if ccol in cg.columns else np.nan
                row[f"{ch}_fixed"] = frow[fcol] if fcol in fg.columns else np.nan
                if ch.startswith("width_") or ch == "total_width":
                    row[f"rel_diff_{ch}"] = _safe_rel_diff(pd.Series([row[f"{ch}_colab"]]), pd.Series([row[f"{ch}_fixed"]])).iloc[0]
                if ch == "br_gaga":
                    row["abs_diff_br_gaga"] = abs(_to_num(pd.Series([row["br_gaga_colab"]])).iloc[0] - _to_num(pd.Series([row["br_gaga_fixed"]])).iloc[0])
            rows.append(row)

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["both_ok"] = out["triple_ok_colab_recomputed"] & out["triple_ok_fixed_campaign"]
    out["only_colab_ok"] = out["triple_ok_colab_recomputed"] & (~out["triple_ok_fixed_campaign"])
    out["only_fixed_ok"] = (~out["triple_ok_colab_recomputed"]) & out["triple_ok_fixed_campaign"]
    out["both_fail"] = (~out["triple_ok_colab_recomputed"]) & (~out["triple_ok_fixed_campaign"])
    return out


def _fraction(num: float, den: float) -> float:
    return float(num) / float(den) if den else float("nan")


def build_group_summary(colab_df: pd.DataFrame, fixed_df: pd.DataFrame, matched_df: pd.DataFrame) -> pd.DataFrame:
    cols = PARAM_KEY_COLS
    colab_n = colab_df.groupby(cols, dropna=False).size().rename("n_colab")
    fixed_n = fixed_df.groupby(cols, dropna=False).size().rename("n_fixed")

    summary = pd.concat([colab_n, fixed_n], axis=1).fillna(0).reset_index()
    summary["n_colab"] = summary["n_colab"].astype(int)
    summary["n_fixed"] = summary["n_fixed"].astype(int)

    if matched_df.empty:
        for c in ["n_matched", "n_triple_ok_colab", "n_triple_ok_fixed", "n_both_ok", "n_only_colab_ok", "n_only_fixed_ok", "n_both_fail"]:
            summary[c] = 0
        return summary

    gb = matched_df.groupby(cols, dropna=False)
    agg = gb.agg(
        n_matched=("m_phi_colab", "size"),
        n_triple_ok_colab=("triple_ok_colab_recomputed", "sum"),
        n_triple_ok_fixed=("triple_ok_fixed_campaign", "sum"),
        n_both_ok=("both_ok", "sum"),
        n_only_colab_ok=("only_colab_ok", "sum"),
        n_only_fixed_ok=("only_fixed_ok", "sum"),
        n_both_fail=("both_fail", "sum"),
        median_rel_diff_width_gg=("rel_diff_width_gg", "median"),
        median_rel_diff_width_gaga=("rel_diff_width_gaga", "median"),
        median_rel_diff_width_Zga=("rel_diff_width_Zga", "median"),
        median_rel_diff_total_width=("rel_diff_total_width", "median"),
        median_abs_diff_br_gaga=("abs_diff_br_gaga", "median"),
        max_rel_diff_width_gg=("rel_diff_width_gg", "max"),
        max_rel_diff_width_gaga=("rel_diff_width_gaga", "max"),
        max_rel_diff_total_width=("rel_diff_total_width", "max"),
    ).reset_index()

    summary = summary.merge(agg, on=cols, how="left")
    for c in ["n_matched", "n_triple_ok_colab", "n_triple_ok_fixed", "n_both_ok", "n_only_colab_ok", "n_only_fixed_ok", "n_both_fail"]:
        summary[c] = summary[c].fillna(0).astype(int)

    summary["fraction_triple_ok_colab"] = summary.apply(lambda r: _fraction(r["n_triple_ok_colab"], r["n_matched"]), axis=1)
    summary["fraction_triple_ok_fixed"] = summary.apply(lambda r: _fraction(r["n_triple_ok_fixed"], r["n_matched"]), axis=1)
    summary["fraction_mask_agreement"] = summary.apply(lambda r: _fraction(r["n_both_ok"] + r["n_both_fail"], r["n_matched"]), axis=1)

    return summary.sort_values(cols).reset_index(drop=True)


def build_mask_confusion_by_group(matched_df: pd.DataFrame) -> pd.DataFrame:
    if matched_df.empty:
        return pd.DataFrame(columns=PARAM_KEY_COLS + ["both_ok", "only_colab_ok", "only_fixed_ok", "both_fail"])
    out = matched_df.groupby(PARAM_KEY_COLS, dropna=False).agg(
        both_ok=("both_ok", "sum"),
        only_colab_ok=("only_colab_ok", "sum"),
        only_fixed_ok=("only_fixed_ok", "sum"),
        both_fail=("both_fail", "sum"),
    ).reset_index()
    return out.sort_values(PARAM_KEY_COLS).reset_index(drop=True)


def build_channel_diff_by_group(matched_df: pd.DataFrame) -> pd.DataFrame:
    if matched_df.empty:
        return pd.DataFrame(columns=PARAM_KEY_COLS)
    gb = matched_df.groupby(PARAM_KEY_COLS, dropna=False)
    out = gb.agg(
        median_rel_diff_width_gg=("rel_diff_width_gg", "median"),
        median_rel_diff_width_gaga=("rel_diff_width_gaga", "median"),
        median_rel_diff_width_Zga=("rel_diff_width_Zga", "median"),
        median_rel_diff_total_width=("rel_diff_total_width", "median"),
        median_abs_diff_br_gaga=("abs_diff_br_gaga", "median"),
        max_rel_diff_width_gg=("rel_diff_width_gg", "max"),
        max_rel_diff_width_gaga=("rel_diff_width_gaga", "max"),
        max_rel_diff_total_width=("rel_diff_total_width", "max"),
    ).reset_index()
    return out.sort_values(PARAM_KEY_COLS).reset_index(drop=True)


def _to_markdown(df: pd.DataFrame, title: str) -> str:
    if df.empty:
        return f"## {title}\n\n(empty)\n"
    try:
        body = df.to_markdown(index=False)
    except Exception:
        cols = list(df.columns)
        header = "| " + " | ".join(cols) + " |"
        sep = "| " + " | ".join(["---"] * len(cols)) + " |"
        lines = []
        for _, r in df.iterrows():
            vals = ["" if pd.isna(r[c]) else str(r[c]) for c in cols]
            lines.append("| " + " | ".join(vals) + " |")
        body = "\n".join([header, sep] + lines)
    return f"## {title}\n\n{body}\n"


def _render_report(out_dir: Path, summary: pd.DataFrame, matched_df: pd.DataFrame, mask_df: pd.DataFrame, merged_csv_path: Path, campaign_root: Path, summary_csv: Path) -> None:
    n_colab = int(summary["n_colab"].sum()) if not summary.empty else 0
    n_fixed = int(summary["n_fixed"].sum()) if not summary.empty else 0
    n_matched = int(summary["n_matched"].sum()) if "n_matched" in summary.columns else 0
    unmatched = (n_colab + n_fixed - 2 * n_matched)

    main_lines = [
        "# Christopher fixed-lambda1 comparison report",
        "",
        f"- merged_csv: {merged_csv_path}",
        f"- campaign_root: {campaign_root}",
        f"- campaign_summary_csv: {summary_csv}",
        f"- n_colab_total: {n_colab}",
        f"- n_fixed_total: {n_fixed}",
        f"- n_matched: {n_matched}",
        f"- n_unmatched_total_proxy: {unmatched}",
    ]
    if not mask_df.empty:
        sums = mask_df[["both_ok", "only_colab_ok", "only_fixed_ok", "both_fail"]].sum()
        main_lines.extend([
            "",
            "## Physical mask totals (matched rows)",
            f"- both_ok: {int(sums['both_ok'])}",
            f"- only_colab_ok: {int(sums['only_colab_ok'])}",
            f"- only_fixed_ok: {int(sums['only_fixed_ok'])}",
            f"- both_fail: {int(sums['both_fail'])}",
        ])
    (out_dir / "comparison_report.md").write_text("\n".join(main_lines) + "\n", encoding="utf-8")

    summary_md = _to_markdown(summary, "Summary by group")
    (out_dir / "comparison_summary_by_group.md").write_text(summary_md, encoding="utf-8")


def make_plots(summary: pd.DataFrame, matched_df: pd.DataFrame, fig_dir: Path) -> None:
    import matplotlib.pyplot as plt

    fig_dir.mkdir(parents=True, exist_ok=True)

    # 1) triple_ok_fraction_colab_vs_fixed
    if not summary.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        order = summary.sort_values(["lambda1", "lambda6"]).reset_index(drop=True)
        x = np.arange(len(order))
        w = 0.38
        ax.bar(x - w / 2, order["fraction_triple_ok_colab"], width=w, label="colab_recomputed")
        ax.bar(x + w / 2, order["fraction_triple_ok_fixed"], width=w, label="fixed_campaign")
        labels = [f"l1={r.lambda1:g}\nl6={r.lambda6:g}" for r in order.itertuples()]
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=0)
        ax.set_ylabel("fraction triple_ok")
        ax.set_title("Triple-OK fraction: Colab recomputed vs fixed campaign")
        ax.legend()
        fig.tight_layout()
        fig.savefig(fig_dir / "triple_ok_fraction_colab_vs_fixed.png", dpi=170)
        plt.close(fig)

    # 2) mask_agreement_by_group
    if not summary.empty:
        fig, ax = plt.subplots(figsize=(10, 5))
        order = summary.sort_values(["lambda1", "lambda6"]).reset_index(drop=True)
        x = np.arange(len(order))
        w = 0.2
        stacks = ["n_both_ok", "n_only_colab_ok", "n_only_fixed_ok", "n_both_fail"]
        labels = ["both_ok", "only_colab_ok", "only_fixed_ok", "both_fail"]
        for i, (c, lbl) in enumerate(zip(stacks, labels)):
            ax.bar(x + (i - 1.5) * w, order[c], width=w, label=lbl)
        ax.set_xticks(x)
        ax.set_xticklabels([f"l1={r.lambda1:g}\nl6={r.lambda6:g}" for r in order.itertuples()])
        ax.set_ylabel("count")
        ax.set_title("Mask agreement by group")
        ax.legend(ncol=2)
        fig.tight_layout()
        fig.savefig(fig_dir / "mask_agreement_by_group.png", dpi=170)
        plt.close(fig)

    def _curve_plot(ycol: str, out_name: str, logy: bool = True) -> None:
        if matched_df.empty or ycol + "_colab" not in matched_df.columns:
            return
        groups = matched_df.groupby(PARAM_KEY_COLS, dropna=False)
        n = len(groups)
        fig, axes = plt.subplots(n, 1, figsize=(10, max(3, 2.5 * n)), sharex=True)
        if n == 1:
            axes = [axes]
        for ax, (k, g) in zip(axes, groups):
            g = g.sort_values("m_phi_colab")
            ax.plot(g["m_phi_colab"], g[f"{ycol}_colab"], label="colab_recomputed", marker="o", ms=2)
            ax.plot(g["m_phi_fixed"], g[f"{ycol}_fixed"], label="fixed_campaign", marker="x", ms=2)
            ax.set_title(f"lambda1={k[0]:g}, lambda6={k[1]:g}, tan_beta={k[2]:g}")
            ax.set_ylabel(ycol)
            if logy:
                vals = pd.concat([_to_num(g[f"{ycol}_colab"]), _to_num(g[f"{ycol}_fixed"])])
                if (vals > 0).all():
                    ax.set_yscale("log")
            ax.grid(alpha=0.25)
        axes[-1].set_xlabel("m_phi")
        axes[0].legend(loc="best")
        fig.tight_layout()
        fig.savefig(fig_dir / out_name, dpi=170)
        plt.close(fig)

    _curve_plot("width_gg", "width_gg_curve_comparison.png", logy=True)
    _curve_plot("width_gaga", "width_gaga_curve_comparison.png", logy=True)
    _curve_plot("br_gaga", "br_gaga_curve_comparison.png", logy=False)
    _curve_plot("total_width", "total_width_curve_comparison.png", logy=True)

    # 7) rel_diff_by_channel_grouped
    if not summary.empty:
        fig, ax = plt.subplots(figsize=(10, 5))
        order = summary.sort_values(["lambda1", "lambda6"]).reset_index(drop=True)
        channels = [
            ("median_rel_diff_width_gg", "gg"),
            ("median_rel_diff_width_gaga", "gaga"),
            ("median_rel_diff_width_Zga", "Zga"),
            ("median_rel_diff_total_width", "total"),
        ]
        x = np.arange(len(order))
        w = 0.18
        for i, (c, lbl) in enumerate(channels):
            ax.bar(x + (i - 1.5) * w, order[c], width=w, label=lbl)
        ax.set_xticks(x)
        ax.set_xticklabels([f"l1={r.lambda1:g}\nl6={r.lambda6:g}" for r in order.itertuples()])
        ax.set_ylabel("median relative diff")
        ax.set_title("Median relative differences by channel and group")
        ax.legend()
        fig.tight_layout()
        fig.savefig(fig_dir / "rel_diff_by_channel_grouped.png", dpi=170)
        plt.close(fig)


def _pick_merged_path(cli_path: Path) -> Path:
    if cli_path.exists():
        return cli_path
    if cli_path == DEFAULT_MERGED_CSV and DEFAULT_FALLBACK_MERGED_CSV.exists():
        return DEFAULT_FALLBACK_MERGED_CSV
    raise FileNotFoundError(f"Merged CSV not found: {cli_path}")


def cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare Christopher Colab/recomputed artifact vs fixed-lambda1 campaign")
    p.add_argument("--merged-csv", type=Path, default=DEFAULT_MERGED_CSV)
    p.add_argument("--campaign-root", type=Path, default=DEFAULT_CAMPAIGN_ROOT)
    p.add_argument("--summary-csv", type=Path, default=DEFAULT_SUMMARY_CSV)
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    p.add_argument("--mphi-tol", type=float, default=1e-4)
    p.add_argument("--param-tol", type=float, default=1e-8)
    p.add_argument("--make-plots", dest="make_plots", action="store_true")
    p.add_argument("--no-plots", dest="make_plots", action="store_false")
    p.set_defaults(make_plots=True)
    return p.parse_args()


def main() -> int:
    args = cli()
    merged_path = _pick_merged_path(args.merged_csv)

    merged_df = _normalize_merged(pd.read_csv(merged_path))
    fixed_df = _load_fixed_campaign(args.campaign_root)

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    matched_df = match_points(merged_df, fixed_df, mphi_tol=args.mphi_tol, param_tol=args.param_tol)
    summary = build_group_summary(merged_df, fixed_df, matched_df)
    mask_df = build_mask_confusion_by_group(matched_df)
    ch_df = build_channel_diff_by_group(matched_df)

    summary.to_csv(out_dir / "comparison_summary_by_group.csv", index=False)
    matched_df.to_csv(out_dir / "matched_points.csv", index=False)
    mask_df.to_csv(out_dir / "mask_confusion_by_group.csv", index=False)
    ch_df.to_csv(out_dir / "channel_diff_by_group.csv", index=False)

    _render_report(out_dir, summary, matched_df, mask_df, merged_path, args.campaign_root, args.summary_csv)

    if args.make_plots:
        make_plots(summary, matched_df, out_dir / "figures")

    n_colab = len(merged_df)
    n_fixed = len(fixed_df)
    n_matched = len(matched_df)
    print(f"merged_csv={merged_path}")
    print(f"n_colab={n_colab} n_fixed={n_fixed} n_matched={n_matched} unmatched_proxy={n_colab + n_fixed - 2*n_matched}")
    print(f"wrote: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
