#!/usr/bin/env python3
"""
compare_chris_vs_2hdmc_lam1eq1.py
=================================

Detailed BR + ctau comparison between Christopher's analytic Stage-1 calculator
(chris/CalcLambda1ScanFixings, carried into the silver CSV as the ``chris_*``
columns) and the Stage-2 2HDMC engine (GenScanWithFixings / GenScanPointEvaluator,
the ``width_*`` / ``total_width`` columns), at the working point lambda_1 = 1.

The two sides live in the *same* silver row (the gen_fixings engine evaluates each
generic-basis variation candidate with 2HDMC and stores the matching analytic
``chris_*`` values + precomputed ``delta_*`` / ``ratio_*`` columns), so NO
point-matching is needed -- unlike scripts/compare_christopher_fixed_lam1.py,
from which only the report *shape* is borrowed.

Each Stage-1 base point (one m_phi at fixed mA, tan_beta, lambda6) spawns a cloud
of ``variation_idx = 0..N-1`` generic-basis neighbours; ``variation_idx == 0`` is
the unperturbed reconstruction (calibration_score ~ 0, masses == Christopher's
input). We report the baseline curve *and* the cloud spread, grouped by mA.

Outputs (under --outdir):
  - ratio_vs_mphi_<channel>.png   ratio (2HDMC / chris) vs m_phi, per mA group
  - br_gaga_vs_mphi.png           BR(H->gaga): 2HDMC vs chris, per mA group
  - ratio_hist.png                histogram of per-channel ratios (cluster at 1?)
  - conditioning.png              calibration_score and recovered lambda_1
  - comparison_summary.csv        per (mA, channel) median/mean/std ratio
  - comparison_report.md          human-readable summary

All physics rows are filtered to triple-OK (positivity & unitarity &
perturbativity) AND stability_ok before any ratio statistics are computed.
"""
from __future__ import annotations

import argparse
import csv
import math
from array import array
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# hbar*c in GeV*mm  (6.582119569e-25 GeV*s * 2.99792458e11 mm/s); matches the
# constant used in GenScanPointEvaluator.cpp and mlpython/2603.
HBAR_C_GEV_MM = 6.582119569e-25 * 2.99792458e11

# In the GENERIC basis the physical masses are *outputs* of (lambda_i, m12^2),
# so a +/-10% coupling jitter does NOT preserve masses: many cloud points drift
# (and some flip alignment to sba ~ -1). calibration_score = sum of squared
# relative errors of the recovered spectrum vs Christopher's input masses
# (mA, mH, mh=125, mHp=mA, sba=1). We treat a cloud point as "mass-preserving"
# (i.e. genuinely "around" Christopher's point) only when score < MASS_SCORE_TOL.
# variation_idx==0 (the unperturbed reconstruction) always has score ~ 0.
MASS_SCORE_TOL = 1e-2          # ~ <4.5% rms per-observable; excludes sba flips
MASS_SCORE_TOL_TIGHT = 1e-3    # ~ <1.4% rms per-observable (reporting only)

# Twice the top mass: above this, H->t tbar opens in 2HDMC but is ABSENT from
# Christopher's 6-channel analytic total, so total-width-derived quantities
# (ctau, br) diverge there even though the per-channel widths still agree.
TWO_MTOP_GEV = 2 * 173.0

DEFAULT_SILVER = Path(
    "data/lam1eq1_var10_2026jun/lake/dihiggs_lake/"
    "campaign=lam1eq1_var10_2026jun/"
    "fixed_sinba=1p0000_l6=0p0001_l7=0p0000_mA=300p0/run01/"
    "tb_01000/scan_tb_1000.csv"
)
DEFAULT_OUTDIR = Path("scripts/out/lam1eq1_var10_2026jun/comparison")

# Channels chris implements (2HDMC width column -> chris width column).
CHANNELS = [
    ("gaga", "width_gaga", "chris_width_gaga"),
    ("Zga", "width_Zga", "chris_width_Zga"),
    ("gg", "width_gg", "chris_width_gg"),
    ("bb", "width_bb", "chris_width_bb"),
    ("tautau", "width_tautau", "chris_width_tautau"),
]

# Only these columns are needed by the report. Loading just these (as packed
# array('d')) keeps memory ~O(N_rows * len(NEEDED) * 8 bytes) -- ~1 GB even at
# 6M rows -- instead of materialising a list of N_rows dicts (which OOMs).
NEEDED_COLUMNS = [
    "positivity_ok", "unitarity_ok", "perturbativity_ok", "stability_ok",
    "calibration_score", "variation_idx", "mA_target", "mH_target",
    "width_bb", "width_tautau", "width_gg", "width_gaga", "width_Zga",
    "chris_width_bb", "chris_width_tautau", "chris_width_gg",
    "chris_width_gaga", "chris_width_Zga", "chris_ctau_mm",
    "ratio_ctau_mm", "br_gaga", "lam1", "delta_m12_2_gen_minus_used",
]


def load_silver(path: Path, needed: list[str] | None = NEEDED_COLUMNS) -> dict[str, np.ndarray]:
    """Stream the silver CSV into a dict of numpy float columns.

    Single pass with csv.reader, accumulating each wanted column into a packed
    array('d') (8 bytes/value, no per-row Python objects). ``needed`` restricts
    which columns are kept (None = all); missing names are silently skipped so
    optional columns (e.g. delta_m12_2_gen_minus_used) don't error.
    """
    with path.open(newline="") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit(f"Empty CSV: {path}")
        index = {name: i for i, name in enumerate(header)}
        wanted = [n for n in (needed or header) if n in index]
        if not wanted:
            raise SystemExit(f"None of the needed columns present in {path}")
        buffers = {n: array("d") for n in wanted}
        positions = [(n, index[n]) for n in wanted]
        for fields in reader:
            for name, pos in positions:
                try:
                    buffers[name].append(float(fields[pos]))
                except (ValueError, IndexError):
                    buffers[name].append(math.nan)
    cols: dict[str, np.ndarray] = {}
    for name in wanted:
        cols[name] = np.array(buffers[name], dtype=np.float64)  # owned, writable copy
        buffers[name] = array("d")  # free the source buffer
    if not cols or next(iter(cols.values())).size == 0:
        raise SystemExit(f"No data rows in {path}")
    return cols


def common5_ctau_ratio(c: dict[str, np.ndarray]) -> np.ndarray:
    """Apples-to-apples ctau ratio using ONLY the 5 channels both sides compute
    (bb, tautau, gg, gaga, Zga) -- removes the H->tt / WW / ZZ / ... coverage
    gap that contaminates the full-total ctau ratio. ratio = ctau_2hdmc/ctau_chris
    = total_chris5 / total_2hdmc5."""
    ch = ["bb", "tautau", "gg", "gaga", "Zga"]
    tot_2hdmc = sum(c[f"width_{x}"] for x in ch)
    tot_chris = sum(c[f"chris_width_{x}"] for x in ch)
    return safe_ratio(tot_chris, tot_2hdmc)


def mass_preserving_mask(c: dict[str, np.ndarray], tol: float = MASS_SCORE_TOL) -> np.ndarray:
    """Cloud points that actually reproduce Christopher's masses (score < tol).
    The baseline (variation_idx==0) always qualifies (score ~ 0)."""
    return c["calibration_score"] < tol


def triple_ok_stable_mask(c: dict[str, np.ndarray]) -> np.ndarray:
    return (
        (c["positivity_ok"] > 0.5)
        & (c["unitarity_ok"] > 0.5)
        & (c["perturbativity_ok"] > 0.5)
        & (c["stability_ok"] > 0.5)
    )


def safe_ratio(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    out = np.full_like(num, np.nan, dtype=float)
    ok = np.isfinite(num) & np.isfinite(den) & (den != 0.0)
    out[ok] = num[ok] / den[ok]
    return out


def group_indices_by_mA(c: dict[str, np.ndarray], mask: np.ndarray) -> dict[float, np.ndarray]:
    groups: dict[float, list[int]] = defaultdict(list)
    mA = c["mA_target"]
    for i in np.where(mask)[0]:
        groups[round(float(mA[i]), 3)].append(int(i))
    return {k: np.array(v) for k, v in sorted(groups.items())}


def plot_ratio_vs_mphi(c, mask, groups, outdir: Path) -> None:
    mphi = c["mH_target"]
    vidx = c["variation_idx"]
    for label, wcol, ccol in CHANNELS:
        ratio = safe_ratio(c[wcol], c[ccol])
        fig, axes = plt.subplots(1, len(groups), figsize=(5 * len(groups), 4), squeeze=False)
        for ax, (mA, idx) in zip(axes[0], groups.items()):
            base = idx[vidx[idx] == 0]
            cloud = idx[vidx[idx] != 0]
            # cloud as faint scatter, baseline as line
            ax.scatter(mphi[cloud], ratio[cloud], s=8, alpha=0.25, color="C0", label="cloud (idx>0)")
            order = np.argsort(mphi[base])
            ax.plot(mphi[base][order], ratio[base][order], "-o", ms=3, color="C3", label="baseline (idx=0)")
            ax.axhline(1.0, color="k", lw=0.8, ls="--")
            ax.set_title(f"mA = {mA:g} GeV")
            ax.set_xlabel("m_H  (m_phi) [GeV]")
            ax.set_ylabel(f"ratio 2HDMC / chris  [{label}]")
            ax.legend(fontsize=7)
        fig.suptitle(f"Width ratio (2HDMC / chris analytic) -- {label} -- lambda_1=1")
        fig.tight_layout()
        fig.savefig(outdir / f"ratio_vs_mphi_{label}.png", dpi=120)
        plt.close(fig)

    # ctau ratio (precomputed column)
    fig, axes = plt.subplots(1, len(groups), figsize=(5 * len(groups), 4), squeeze=False)
    for ax, (mA, idx) in zip(axes[0], groups.items()):
        base = idx[vidx[idx] == 0]
        cloud = idx[vidx[idx] != 0]
        r = c["ratio_ctau_mm"]
        ax.scatter(mphi[cloud], r[cloud], s=8, alpha=0.25, color="C0", label="cloud (idx>0)")
        order = np.argsort(mphi[base])
        ax.plot(mphi[base][order], r[base][order], "-o", ms=3, color="C3", label="baseline (idx=0)")
        ax.axhline(1.0, color="k", lw=0.8, ls="--")
        ax.set_title(f"mA = {mA:g} GeV")
        ax.set_xlabel("m_H  (m_phi) [GeV]")
        ax.set_ylabel("ratio ctau 2HDMC / chris")
        ax.legend(fontsize=7)
    fig.suptitle("ctau ratio (2HDMC / chris analytic) -- lambda_1=1")
    fig.tight_layout()
    fig.savefig(outdir / "ratio_vs_mphi_ctau.png", dpi=120)
    plt.close(fig)


def plot_br_gaga(c, groups, outdir: Path) -> None:
    """BR(H->gaga): 2HDMC (br_gaga) vs chris (chris_width_gaga / chris_total_width)."""
    mphi = c["mH_target"]
    vidx = c["variation_idx"]
    chris_tot = safe_ratio(np.full_like(c["chris_ctau_mm"], HBAR_C_GEV_MM), c["chris_ctau_mm"])
    chris_br_gaga = safe_ratio(c["chris_width_gaga"], chris_tot)
    fig, axes = plt.subplots(1, len(groups), figsize=(5 * len(groups), 4), squeeze=False)
    for ax, (mA, idx) in zip(axes[0], groups.items()):
        base = idx[vidx[idx] == 0]
        order = np.argsort(mphi[base])
        ax.plot(mphi[base][order], c["br_gaga"][base][order], "-o", ms=3, color="C0", label="2HDMC br_gaga")
        ax.plot(mphi[base][order], chris_br_gaga[base][order], "-s", ms=3, color="C3", label="chris br_gaga")
        ax.set_yscale("log")
        ax.set_title(f"mA = {mA:g} GeV")
        ax.set_xlabel("m_H  (m_phi) [GeV]")
        ax.set_ylabel("BR(H -> gamma gamma)")
        ax.legend(fontsize=7)
    fig.suptitle("BR(H->gaga): 2HDMC vs chris analytic (baseline idx=0) -- lambda_1=1")
    fig.tight_layout()
    fig.savefig(outdir / "br_gaga_vs_mphi.png", dpi=120)
    plt.close(fig)


def plot_ratio_hist(c, mask, outdir: Path) -> None:
    # Restrict to triple-OK + stable AND mass-preserving so sba-flipped /
    # mass-drifted cloud points don't pollute the spread.
    clean = mask & mass_preserving_mask(c)
    fig, axes = plt.subplots(2, 4, figsize=(17, 7))
    items = [(lbl, safe_ratio(c[w], c[cc])) for lbl, w, cc in CHANNELS]
    items.append(("ctau (full total)", c["ratio_ctau_mm"]))
    items.append(("ctau (common 5ch)", common5_ctau_ratio(c)))
    for ax, (lbl, r) in zip(axes.ravel(), items):
        vals = r[clean]
        vals = vals[np.isfinite(vals) & (vals > 0)]
        if vals.size:
            ax.hist(np.log10(vals), bins=60, color="C0", alpha=0.8)
            ax.axvline(0.0, color="k", ls="--", lw=0.8)  # ratio == 1
            med = float(np.median(vals))
            ax.set_title(f"{lbl}: median={med:.3g}  (n={vals.size})")
        ax.set_xlabel(f"log10(ratio 2HDMC/chris) [{lbl}]")
    for ax in axes.ravel()[len(items):]:
        ax.axis("off")
    fig.suptitle("Ratio distributions (triple-OK + stable + mass-preserving) -- 0 = agreement")
    fig.tight_layout()
    fig.savefig(outdir / "ratio_hist.png", dpi=120)
    plt.close(fig)


def plot_conditioning(c, outdir: Path) -> None:
    vidx = c["variation_idx"]
    base = vidx == 0
    cloud = vidx != 0
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    # calibration_score vs m_phi
    ax = axes[0]
    ax.scatter(c["mH_target"][cloud], np.log10(np.clip(c["calibration_score"][cloud], 1e-30, None)),
               s=6, alpha=0.2, color="C0", label="cloud")
    ax.scatter(c["mH_target"][base], np.log10(np.clip(c["calibration_score"][base], 1e-30, None)),
               s=10, color="C3", label="baseline")
    ax.set_xlabel("m_H [GeV]")
    ax.set_ylabel("log10(calibration_score)")
    ax.set_title("mass-match score (0 = exact)")
    ax.legend(fontsize=7)
    # recovered lambda_1
    ax = axes[1]
    if "lam1" in c:
        ax.hist(c["lam1"][np.isfinite(c["lam1"])], bins=60, color="C0", alpha=0.8)
        ax.axvline(1.0, color="k", ls="--", lw=0.8)
    ax.set_xlabel("recovered lam1 (2HDMC get_param_gen)")
    ax.set_title("lambda_1 (target = 1)")
    # m12_2 generation delta if present
    ax = axes[2]
    key = "delta_m12_2_gen_minus_used"
    if key in c:
        d = c[key][np.isfinite(c[key])]
        ax.hist(d, bins=60, color="C0", alpha=0.8)
        ax.set_xlabel(key)
        ax.set_title("generic-basis m12^2 replay delta")
    else:
        ax.axis("off")
    fig.suptitle("Conditioning / calibration diagnostics -- lambda_1=1")
    fig.tight_layout()
    fig.savefig(outdir / "conditioning.png", dpi=120)
    plt.close(fig)


def write_summary(c, mask, groups, outdir: Path) -> Path:
    vidx = c["variation_idx"]
    massok = mass_preserving_mask(c)
    rows = []
    channels = [(lbl, safe_ratio(c[w], c[cc])) for lbl, w, cc in CHANNELS]
    channels.append(("ctau_full", c["ratio_ctau_mm"]))
    channels.append(("ctau_common5", common5_ctau_ratio(c)))
    for mA, idx in groups.items():
        idx_base = idx[vidx[idx] == 0]
        idx_cloud = idx[massok[idx]]  # baseline + mass-preserving cloud points
        for lbl, ratio in channels:
            # 'baseline' = unperturbed point; 'cloud_massok' = mass-preserving cloud
            for scope, ii in (("baseline", idx_base), ("cloud_massok", idx_cloud)):
                vals = ratio[ii]
                vals = vals[np.isfinite(vals)]
                if vals.size == 0:
                    rows.append((mA, lbl, scope, 0, math.nan, math.nan, math.nan))
                    continue
                rows.append((mA, lbl, scope, int(vals.size),
                             float(np.median(vals)), float(np.mean(vals)), float(np.std(vals))))
    out_csv = outdir / "comparison_summary.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mA", "channel", "scope", "n", "ratio_median", "ratio_mean", "ratio_std"])
        for r in rows:
            w.writerow(r)
    return out_csv


def write_report(c, mask, groups, outdir: Path, silver: Path) -> None:
    vidx = c["variation_idx"]
    n_total = c["variation_idx"].size
    n_ok = int(mask.sum())
    n_cloud = int((vidx > 0).sum())
    n_cloud_massok = int(((vidx > 0) & mass_preserving_mask(c, MASS_SCORE_TOL)).sum())
    n_cloud_tight = int(((vidx > 0) & mass_preserving_mask(c, MASS_SCORE_TOL_TIGHT)).sum())
    n_base = int((vidx == 0).sum())
    # survivors per base point (out of N variations)
    surv_counts = []
    base_keys = {}
    for i in range(n_total):
        key = (round(float(c["mA_target"][i]), 3), round(float(c["mH_target"][i]), 4))
        base_keys.setdefault(key, []).append(i)
    for key, idxs in base_keys.items():
        surv_counts.append(sum(1 for i in idxs if c["calibration_score"][i] < MASS_SCORE_TOL))
    surv_mean = float(np.mean(surv_counts)) if surv_counts else float("nan")
    n_per = max((len(v) for v in base_keys.values()), default=0)

    lines = [
        "# Chris (analytic) vs 2HDMC BR/ctau comparison -- lambda_1 = 1",
        "",
        f"- Silver CSV: `{silver}`",
        f"- Total silver rows (base x variations): **{n_total}**  "
        f"(base points: {n_base}, cloud points: {n_cloud})",
        f"- Triple-OK + stable rows: **{n_ok}**",
        f"- mA groups: {', '.join(f'{m:g}' for m in groups)}",
        "",
        "Both sides live in the same row (no point-matching). `ratio = 2HDMC / chris`;",
        "1.0 = perfect agreement. `variation_idx==0` is the unperturbed baseline",
        "(masses == Christopher's input, calibration_score ~ 0).",
        "",
        "## Mass preservation of the +/-10% cloud (generic-basis caveat)",
        "",
        "In the generic basis the masses are *outputs* of (lambda_i, m12^2), so a",
        "+/-10% coupling jitter does not hold masses fixed; many cloud points drift",
        "and some flip alignment (sba -> -1). Treat only cloud points with",
        f"`calibration_score < {MASS_SCORE_TOL:g}` as genuinely \"around\" Christopher's masses.",
        "",
        f"- cloud points reproducing masses (score < {MASS_SCORE_TOL:g}): "
        f"**{n_cloud_massok} / {n_cloud}**",
        f"- tighter (score < {MASS_SCORE_TOL_TIGHT:g}): **{n_cloud_tight} / {n_cloud}**",
        f"- mean survivors per base point (of {n_per}): **{surv_mean:.2f}**",
        "",
        "Cloud statistics below use only mass-preserving points (baseline always qualifies).",
        "",
        "## Median ratio (2HDMC / chris)",
        "",
        "| mA | channel | n_base | median (baseline) | n_cloud_ok | median (cloud, mass-ok) |",
        "|----|---------|--------|-------------------|------------|--------------------------|",
    ]
    massok = mass_preserving_mask(c)
    channels = [(lbl, safe_ratio(c[w], c[cc])) for lbl, w, cc in CHANNELS]
    channels.append(("ctau_full", c["ratio_ctau_mm"]))
    channels.append(("ctau_common5", common5_ctau_ratio(c)))
    for mA, idx in groups.items():
        base = idx[vidx[idx] == 0]
        cloud = idx[massok[idx]]
        for lbl, ratio in channels:
            vb = ratio[base][np.isfinite(ratio[base])]
            vc = ratio[cloud][np.isfinite(ratio[cloud])]
            mb = float(np.median(vb)) if vb.size else math.nan
            mc = float(np.median(vc)) if vc.size else math.nan
            lines.append(f"| {mA:g} | {lbl} | {vb.size} | {mb:.4g} | {vc.size} | {mc:.4g} |")
    lines += [
        "",
        "## Interpretation",
        "",
        "- **Per-channel partial widths** (bb, tautau, gg, gaga, Zga) are the clean,",
        "  apples-to-apples comparison: tree channels (bb, tautau) agree to ~1e-3;",
        "  loop channels (gaga, Zga, gg) agree to a few percent (analytic vs full 2HDMC).",
        f"- **Full-total ctau / br diverge for any point with m_H > ~2 m_top ({TWO_MTOP_GEV:.0f} GeV)**",
        "  because 2HDMC's total width then includes H->t tbar, a channel ABSENT from",
        "  Christopher's 6-channel analytic total. This is channel coverage, NOT a",
        "  calculation error. The effect grows with mA (more of the [130,mA] m_H range",
        "  lies above 2 m_top): the baseline median ctau_full stays ~1 up to mA=450,",
        "  sags at mA=500, and collapses at mA=600. Use the `ctau_common5` ratio",
        "  (built from the 5 shared channels on both sides) for an apples-to-apples",
        "  lifetime comparison -- it stays near 1 across all mA.",
        "- Cloud loop-channel spread also reflects a one-sided effect: 2HDMC sees the",
        "  jittered lambda_1 while the bronze `chris_*` widths were computed once at",
        "  lambda_1=1. The baseline (idx=0) is free of this.",
        "- Not run: the 7 `02_results/chris_hybrid` reference points (they use lambda_1=0,",
        "  so they are not directly comparable to this lambda_1=1 campaign).",
        "",
        "## Figures",
        "- `ratio_vs_mphi_<channel>.png` -- width ratio vs m_H per mA (baseline line + cloud)",
        "- `ratio_vs_mphi_ctau.png` -- full-total ctau ratio vs m_H per mA",
        "- `br_gaga_vs_mphi.png` -- BR(H->gaga) 2HDMC vs chris",
        "- `ratio_hist.png` -- ratio distributions (triple-OK + stable + mass-preserving)",
        "- `conditioning.png` -- calibration_score, recovered lambda_1, m12^2 replay delta",
        "",
    ]
    (outdir / "comparison_report.md").write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--silver-csv", type=Path, default=DEFAULT_SILVER)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = ap.parse_args()

    if not args.silver_csv.exists():
        raise SystemExit(f"silver CSV not found: {args.silver_csv}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    c = load_silver(args.silver_csv)
    mask = triple_ok_stable_mask(c)
    groups = group_indices_by_mA(c, mask)
    if not groups:
        raise SystemExit("No triple-OK + stable rows to compare.")

    plot_ratio_vs_mphi(c, mask, groups, args.outdir)
    plot_br_gaga(c, groups, args.outdir)
    plot_ratio_hist(c, mask, args.outdir)
    plot_conditioning(c, args.outdir)
    out_csv = write_summary(c, mask, groups, args.outdir)
    write_report(c, mask, groups, args.outdir, args.silver_csv)

    print(f"[OK] wrote comparison artifacts to {args.outdir}")
    print(f"     summary: {out_csv}")
    print(f"     report:  {args.outdir / 'comparison_report.md'}")


if __name__ == "__main__":
    main()
