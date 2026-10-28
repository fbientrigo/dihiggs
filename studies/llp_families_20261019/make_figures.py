"""Diagnostic figure from raw/scan_ledger.csv + validated_benchmarks.csv only:  python3 scripts/make_figures.py
Left: B_h and B_gamma vs tan(beta) on the ctau = 100 mm surface (lambda6 solved); right: |lambda6| required."""
import os
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get("DIHIGGS_ROOT") or Path(__file__).resolve().parents[3])
H = ROOT / "reports/llp_families_20261019"  # output package (lake), not code
rows = [r for r in csv.DictReader(open(H / "raw/scan_ledger.csv")) if r["phase"] in ("coarse", "final_rounded", "hadronic_solve")]
vb = list(csv.DictReader(open(H / "validated_benchmarks.csv")))
col = {"200": "#4C72B0", "250": "#DD8452", "300": "#55A868"}
fig, (a, b) = plt.subplots(1, 2, figsize=(9.2, 3.6))
for m in ("200", "250", "300"):
    for sign, ls in (("-", "-"), ("+", ":")):
        pts = sorted((float(r["tan_beta"]), float(r["B_h"]), float(r["B_gamma"]), abs(float(r["lambda6"])))
                     for r in rows if r["mphi_GeV"] == m and r["sign"] == sign and r["phase"] != "hadronic_solve")
        if not pts:
            continue
        t, bh, bg, l6 = zip(*pts)
        lab = rf"$m_\phi={m}$ GeV" if sign == "-" else None
        a.plot(t, bg, ls, color=col[m], lw=1.2, label=lab)
        a.plot(t, bh, ls, color=col[m], lw=0.8, alpha=0.6)
        b.plot(t, l6, ls, color=col[m], lw=1.2, label=(lab if sign == "-" else None))
for v in vb:
    a.plot(float(v["tan_beta"]), float(v["B_gamma"]), "o", color=col[v["mphi_GeV"]], ms=4, mec="k", mew=0.4)
for y, txt in ((0.75, r"hadronic: $B_h\geq0.75$"), (0.95, r"photonic: $B_\gamma\geq0.95$"), (0.5, r"mixed: $B_\gamma\geq0.5$")):
    a.axhline(y, color="0.6", lw=0.6, ls="--"); a.text(1.5e6, y + 0.012, txt, fontsize=6.5, color="0.35")
a.set_xscale("log"); a.set_xlim(1.4e6, 1e8); a.set_xlabel(r"$\tan\beta$"); a.set_ylabel(r"$B_\gamma$ (thick), $B_h$ (thin)")
a.set_title(r"$c\tau_\phi=100$ mm, $Y_2=0$; solid $\lambda_6<0$, dotted $\lambda_6>0$", fontsize=8)
a.legend(fontsize=7, frameon=False, loc="center right")
b.set_xscale("log"); b.set_xlim(1.4e6, 1e8); b.set_yscale("log"); b.set_xlabel(r"$\tan\beta$"); b.set_ylabel(r"$|\lambda_6|$ solving $c\tau_\phi=100$ mm")
b.legend(fontsize=7, frameon=False)
fig.tight_layout()
(H / "figs").mkdir(exist_ok=True)
fig.savefig(H / "figs/composition_vs_tanb.pdf"); fig.savefig(H / "figs/composition_vs_tanb.png", dpi=150)
print("written figs/composition_vs_tanb.{pdf,png}")
