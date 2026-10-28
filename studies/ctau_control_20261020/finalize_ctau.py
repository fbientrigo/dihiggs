"""Ranges, selected points, figure, LaTeX table and snippet from raw/scan_points.jsonl + raw/edges.jsonl only (no backend).

  python3 scripts/finalize_ctau.py
Validated segment = contiguous run of grid points around the reference whose status is PASS, where a failure of
ONLY the strict October gate C_V = -v Z7 (relative 1e-9) is accepted when the deviation is below 10 t eps
(double-precision round-off of a cancellation-dominated Z7, the tolerance form of the October kappa_top gate);
claimed only up to T_CAP (above it the double-precision backend's c tau t^-2 drifts by > 1e-6 and its mass
read-back by > 1e-4 GeV).  Both the strict and the round-off-tolerant ranges are written.
"""
import os
import csv, json, math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get("DIHIGGS_ROOT") or Path(__file__).resolve().parents[3])
H = ROOT / "reports/ctau_control_20261020"  # output package (lake), not code
EPS = 2.220446049250313e-16
T_CAP = 1e10
FAM = ("hadronic", "mixed", "photonic")
COL = {"hadronic": "#4C72B0", "mixed": "#C44E52", "photonic": "#55A868"}
REF = {"hadronic": "A6", "mixed": "M250p", "photonic": "P250n"}
rows = [json.loads(l) for l in open(H / "raw/scan_points.jsonl")]
edges = json.load(open(H / "raw/edges.json"))


def ok_roundoff(r):
    if r["status"] == "PASS":
        return True
    return (r["status"] == "INCONCLUSIVE" and r["failed_gates"] == "charged_loop_coupling_CV_eq_minus_vZ7"
            and r["CV_vs_minus_vZ7_rel"] < 10 * float(r["tan_beta"]) * EPS)


def segment(grid, t0, good):
    i0 = min(range(len(grid)), key=lambda i: abs(math.log(float(grid[i]["tan_beta"]) / t0)))
    lo = hi = i0
    while lo > 0 and good(grid[lo - 1]):
        lo -= 1
    while hi < len(grid) - 1 and good(grid[hi + 1]):
        hi += 1
    return grid[lo], grid[hi]


ranges, grids = {}, {}
for f in FAM:
    ref = next(r for r in rows if r["family"] == f and r["experiment"] == "reference")
    grid = sorted((r for r in rows if r["family"] == f and r["experiment"] == "A_grid"), key=lambda r: float(r["tan_beta"]))
    grids[f] = grid
    t0 = float(ref["tan_beta"])
    s_lo, s_hi = segment(grid, t0, lambda r: r["status"] == "PASS")
    v_lo, v_hi = segment(grid, t0, lambda r: ok_roundoff(r) and float(r["tan_beta"]) <= T_CAP * 1.0001)
    seg = [r for r in grid if float(v_lo["tan_beta"]) <= float(r["tan_beta"]) <= float(v_hi["tan_beta"])]
    b1 = json.load(open(H / "raw/summary.json"))[f]["B1"]
    slopes = [math.log(b["ctau_mm"] / a["ctau_mm"]) / math.log(float(b["tan_beta"]) / float(a["tan_beta"]))
              for a, b in zip(seg, seg[1:])]
    ranges[f] = dict(
        reference=REF[f], X=ref["X"], t_ref=ref["tan_beta"], ctau_ref_mm=ref["ctau_mm"], B_gamma_ref=ref["B_gamma"],
        tested_t=[grid[0]["tan_beta"], grid[-1]["tan_beta"]], tested_ctau_mm=[grid[0].get("ctau_mm"), grid[-1].get("ctau_mm")],
        n_grid=len(grid), n_grid_pass_strict=sum(r["status"] == "PASS" for r in grid),
        strict_segment_t=[s_lo["tan_beta"], s_hi["tan_beta"]], strict_segment_ctau_mm=[s_lo["ctau_mm"], s_hi["ctau_mm"]],
        validated_t=[v_lo["tan_beta"], v_hi["tan_beta"]], validated_ctau_mm=[v_lo["ctau_mm"], v_hi["ctau_mm"]],
        validated_n_points=len(seg), validated_n_CV_roundoff_only=sum(r["status"] != "PASS" for r in seg),
        slope_grid_min=min(slopes), slope_grid_max=max(slopes), slope_local_B1=[b["dln_ctau_dln_t"] for b in b1],
        max_abs_dB_gamma=max(abs(r["dB_gamma"]) for r in seg), max_abs_dBR=max(r["max_abs_dBR"] for r in seg),
        max_rel_dBR=max(r["max_rel_dBR"] for r in seg),
        class_along_segment=sorted({r["family_class"] for r in seg}),
        lower_edge_reason=("BFB_general (exact sufficient certificate) INCONCLUSIVE below t=" + edges[f"{f}_BFB_PASS_from"]["tan_beta"]
                           + ", exact negative witness (FAIL) below t=" + edges[f"{f}_BFB_INCONCLUSIVE_from"]["tan_beta"]
                           if f"{f}_BFB_PASS_from" in edges else
                           "lambda6 = 0: no BFB edge; the October gate 'lambda7 within x2 of its leading term -2d/t' fails at t = 1 "
                           "(lambda7 = +0.076, still nonzero); all 8 tree tests pass down to t = 1"),
        upper_edge_reason=f"not a physics edge: claim capped at t = {T_CAP:.0e} (double-precision backend); all gates except the "
                          "strict C_V round-off gate pass up to t = 1e13")
    for name, seg2 in (("restricted_0.1_1e4mm", [r for r in seg if 0.05 <= r["ctau_mm"] <= 2e4]),):
        ranges[f][name] = dict(max_abs_dB_gamma=max(abs(r["dB_gamma"]) for r in seg2), max_rel_dBR=max(r["max_rel_dBR"] for r in seg2),
                               slope_min=min(math.log(b["ctau_mm"] / a["ctau_mm"]) / math.log(float(b["tan_beta"]) / float(a["tan_beta"]))
                                             for a, b in zip(seg2, seg2[1:])),
                               slope_max=max(math.log(b["ctau_mm"] / a["ctau_mm"]) / math.log(float(b["tan_beta"]) / float(a["tan_beta"]))
                                             for a, b in zip(seg2, seg2[1:])))
json.dump(ranges, open(H / "raw/ranges.json", "w"), indent=1)

# ---------------- selected points: lifetime targets of each family + reference ----------------
COLS = ["family", "role", "target_ctau_mm", "ctau_mm", "tan_beta", "X", "lambda6", "lambda1", "lambda2", "lambda3", "lambda4", "lambda5",
        "lambda7", "m12_sq_GeV2", "Z7", "Y2_GeV2", "g_hphiphi_GeV", "Gamma_GeV", "Gamma_bb_GeV", "Gamma_gluglu_GeV", "Gamma_Zgamma_GeV",
        "Gamma_gammagamma_GeV", "BR_bb", "BR_gluglu", "BR_Zgamma", "BR_gammagamma", "BR_tautau", "BR_cc", "BR_ttstar", "BR_other_terms",
        "B_h", "B_gamma", "dB_gamma", "max_abs_dBR", "family_class", "status", "failed_gates", "max_abs_quartic", "backend_mass_max_abs_GeV",
        "lambda1_from_lambda7_residual"] + [f"tree_{k}" for k in ("spectrum_reconstruction", "stationarity", "alignment", "local_minimum_tree",
                                                                  "BFB_general", "global_EW_minimum", "full_tree_unitarity", "quartic_screen_4pi")]
sel = []
for f in FAM:
    ref = next(r for r in rows if r["family"] == f and r["experiment"] == "reference")
    sel.append(dict(ref, role=f"reference ({REF[f]}, 2026-10-19)", target_ctau_mm="", dB_gamma=0.0, max_abs_dBR=0.0))
    for r in rows:
        if r["family"] == f and r["experiment"] == "A_target" and 0.1 <= r["target_ctau_mm"] <= 1e4 and r["target_ctau_mm"] != 100:
            sel.append(dict(r, role="lifetime target, fixed X"))
with open(H / "selected_points.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, COLS, extrasaction="ignore")
    w.writeheader()
    w.writerows(sel)
# input file for the independent (no project code) potential check of 2026-10-19
with open(H / "raw/selected_for_independent_check/validated_benchmarks.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, ["candidate_id", "mphi_GeV", "tan_beta", "m12_sq_GeV2"] + [f"lambda{i}" for i in range(1, 8)])
    w.writeheader()
    for i, r in enumerate(sel):
        w.writerow(dict(candidate_id=f"{r['family'][0]}{i}_ctau{r['ctau_mm']:.3g}", mphi_GeV="250", tan_beta=r["tan_beta"],
                        m12_sq_GeV2=r["m12_sq_GeV2"], **{f"lambda{k}": r[f"lambda{k}"] for k in range(1, 8)}))

# ---------------- figure ----------------
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7})
fig, (a, b) = plt.subplots(2, 1, figsize=(5.6, 5.4), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1], hspace=0.06))
xlo, xhi = 1e2, 3e10
a.axhspan(0.1, 1e4, color="0.93", lw=0, zorder=0)
a.text(1.3e2, 1.6e3, r"requested range $0.1$–$10^4$ mm", fontsize=7.5, color="0.35")
a.axhline(100, color="0.55", lw=0.6, ls=(0, (4, 3)), zorder=1)
a.text(1.3e2, 125, "100 mm: reference benchmark", fontsize=7.5, color="0.35")
for f in FAM:
    rg = ranges[f]
    tlo, thi = float(rg["validated_t"][0]), float(rg["validated_t"][1])
    g = [r for r in grids[f] if "ctau_mm" in r and xlo <= float(r["tan_beta"]) <= xhi]
    seg = [r for r in g if tlo <= float(r["tan_beta"]) <= thi]
    badl = [r for r in g if float(r["tan_beta"]) < tlo]
    X = float(rg["X"])
    lab = {"hadronic": r"hadronic, $X=0$", "mixed": rf"mixed, $X={X:+.0f}$", "photonic": rf"photonic, $X={X:+.0f}$"}[f]
    for ax, key in ((a, "ctau_mm"), (b, "B_gamma")):
        ax.plot([float(r["tan_beta"]) for r in seg], [r[key] for r in seg], color=COL[f], lw=1.6, label=lab if ax is a else None, zorder=3)
        if badl:
            ax.plot([float(r["tan_beta"]) for r in badl + seg[:1]], [r[key] for r in badl + seg[:1]], color=COL[f], lw=0.9, ls=":", zorder=2)
            ax.scatter([float(r["tan_beta"]) for r in badl], [r[key] for r in badl], marker="x", s=9, lw=0.6, color=COL[f], zorder=2)
        ax.scatter([float(rg["t_ref"])], [rg["ctau_ref_mm"] if key == "ctau_mm" else rg["B_gamma_ref"]], marker="*", s=70,
                   color=COL[f], edgecolor="k", lw=0.5, zorder=5)
a.set_xscale("log"); a.set_yscale("log")
a.set_xlim(xlo, xhi); a.set_ylim(1e-8, 1e9)
a.set_ylabel(r"$c\tau_\phi$ [mm]")
a.plot([], [], color="0.4", lw=0.9, ls=":", marker="x", ms=3, label="tested, fails BFB (or INCONCLUSIVE)")
a.scatter([], [], marker="*", s=60, color="w", edgecolor="k", lw=0.5, label="2026-10-19 benchmark")
a.legend(loc="lower right", fontsize=7.2, frameon=False, bbox_to_anchor=(1.0, 0.07))
a.text(0.99, 0.03, r"$m_\phi=250$ GeV, $Y_2=0$, $\lambda_7\neq0$; fixed $X=\lambda_6\tan\beta$", transform=a.transAxes,
       ha="right", va="bottom", fontsize=7.5)
for y in (0.5, 0.95):
    b.axhline(y, color="0.6", lw=0.5, ls=(0, (2, 2)))
b.text(1.25e10, 0.95, "0.95", fontsize=6.8, color="0.4", va="center", bbox=dict(fc="w", ec="none", pad=0.5))
b.text(1.25e10, 0.5, "0.5", fontsize=6.8, color="0.4", va="center", bbox=dict(fc="w", ec="none", pad=0.5))
b.set_ylim(-0.03, 1.08)
b.set_ylabel(r"$B_\gamma=$BR$(\gamma\gamma)+$BR$(Z\gamma)$", fontsize=8)
b.set_xlabel(r"$\tan\beta$")
for ax in (a, b):
    ax.tick_params(which="both", direction="in", top=True, right=True)
fig.savefig(H / "figs/ctau_vs_tanb_250.pdf", bbox_inches="tight")
fig.savefig(H / "figs/ctau_vs_tanb_250.png", dpi=220, bbox_inches="tight")

# ---------------- LaTeX table (mixed 250 GeV) ----------------
def sci(x, n=3):
    m, e = f"{x:.{n - 1}e}".split("e")
    e = int(e)
    return m if e == 0 else rf"{m}\times10^{{{e}}}"


def ctau_tex(x):
    return f"{x:.6g}" if 1e-3 <= x < 1e6 else sci(x, 4)


mx = sorted([r for r in rows if r["family"] == "mixed" and r["experiment"] == "A_target" and 0.1 <= r["target_ctau_mm"] <= 1e4],
            key=lambda r: r["target_ctau_mm"])
rgm = ranges["mixed"]
edge_lo = next(r for r in grids["mixed"] if r["tan_beta"] == rgm["validated_t"][0])
edge_hi = next(r for r in grids["mixed"] if r["tan_beta"] == rgm["validated_t"][1])


def line(tg, r):
    st = "PASS" if r["status"] == "PASS" else r["status"]
    d = r["dB_gamma"]
    dtex = "0" if d == 0 else sci(d, 2)
    return (f"{tg} & ${ctau_tex(r['ctau_mm'])}$ & ${sci(float(r['tan_beta']), 6 if float(r['tan_beta']) < 1e9 else 3)}$ & "
            f"${sci(float(r['lambda6']), 4)}$ & {float(r['X']):.5f} & {r['B_gamma']:.6f} & ${dtex}$ & {st} \\\\")


body = [line(f"${ctau_tex(r['target_ctau_mm'])}$", r) for r in mx]
tab = r"""%% generated by reports/ctau_control_20261020/scripts/finalize_ctau.py -- do not edit by hand
\begin{table}[htbp]
\centering
\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{rrrrrrrl}
\hline
target $c\tau$ [mm] & $c\tau_\phi$ [mm] & $\tan\beta$ & $\lambda_6=X/\tan\beta$ & $X$ & $B_\gamma$ & $\Delta B_\gamma$ & status \\
\hline
""" + "\n".join(body) + r"""
\hline
\multicolumn{8}{l}{\emph{limits of the validated segment on this trajectory}}\\
""" + line("lower edge (BFB)", edge_lo) + "\n" + line(r"cap (numerics)", edge_hi) + r"""
\hline
\end{tabular}
\caption{Lifetime control along the fixed-$X$ trajectory of the mixed benchmark M250p ($m_\phi=250$~GeV, $Y_2=0$,
$X=\lambda_6\tan\beta$ held fixed). Only $\tan\beta$ is changed; $\lambda_6=X/\tan\beta$ exactly (rounded here for display),
and $\lambda_{2\ldots5}$, $m_{12}^2$ and the nonzero generic-basis $\lambda_7$ are reconstructed at 50 digits for every point.
$\lambda_7$ is not shown: at $\tan\beta\sim10^{7}$ a rounded $\lambda_7$ is not a usable input
($\partial\lambda_1/\partial\lambda_7=\tan^3\!\beta/2$); full-precision values are in \texttt{selected\_points.csv} and \texttt{scan\_points.csv}.
$\Delta B_\gamma$ is relative to the 100~mm reference. Status: tree-level theory and numerical validation only
(BFB, unitarity, perturbativity, local and sufficient global-minimum certificate); collider constraints not tested.
The lower edge is set by bounded-from-below ($|Z_7|\simeq1.1$); the upper row is where the claim is capped because the
double-precision decay backend degrades, not a physical limit.}
\label{tab:ctau-control-mixed250}
\end{table}
"""
open(H / "tex/table_ctau_mixed250.tex", "w").write(tab)
print(json.dumps(ranges, indent=1))
