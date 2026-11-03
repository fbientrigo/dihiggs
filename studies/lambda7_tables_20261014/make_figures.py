"""Two figures from the CSVs of this package only:  python3 scripts/make_figures.py"""
import os
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get("DIHIGGS_ROOT") or Path(__file__).resolve().parents[3])
HERE = ROOT / "reports/lambda7_tables_20261014"  # output package (lake), not code
OUT = HERE / "note/figs"; OUT.mkdir(parents=True, exist_ok=True)
bm = {r["id"]: r for r in csv.DictReader(open(HERE / "BENCHMARKS.csv"))}
dc = {r["id"]: r for r in csv.DictReader(open(HERE / "DECAYS.csv"))}
pr = {r["id"]: r for r in csv.DictReader(open(HERE / "PRODUCTION.csv"))}
ids = list(bm)

# 1) decay composition
chans = [("BR_bb", r"$b\bar b$", "#4C72B0"), ("BR_gg", "$gg$", "#55A868"), ("BR_gammagamma", r"$\gamma\gamma$", "#C44E52"),
         ("BR_Zgamma", r"$Z\gamma$", "#DD8452"), ("BR_other_sum", r"other ($c\bar c,\tau\tau,\mu\mu$)", "#8C8C8C")]
fig, ax = plt.subplots(figsize=(7.5, 3.6))
bottom = [0.0] * len(ids)
for key, lab, col in chans:
    vals = [float(dc[i][key]) for i in ids]
    ax.bar(ids, vals, bottom=bottom, label=lab, color=col, edgecolor="white", linewidth=0.5)
    bottom = [b + v for b, v in zip(bottom, vals)]
for k, i in enumerate(ids):
    ax.text(k, 1.02, f"{float(dc[i]['ctau_mm']):.2g} mm", ha="center", fontsize=7)
ax.set_ylim(0, 1.12); ax.set_ylabel("branching ratio")
ax.set_title(r"$m_\phi=200$ GeV, $m_A=m_{H^\pm}=300$ GeV; label = $c\tau_\phi$", fontsize=9)
ax.legend(ncol=5, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.1), frameon=False)
fig.tight_layout(); fig.savefig(OUT / "br_composition.pdf"); fig.savefig(OUT / "br_composition.png", dpi=150)

# 2) sigma vs g/g0 (direct integrations) with the g^2 line through the g0 integration
fig, ax = plt.subplots(figsize=(4.6, 3.4))
x = [float(bm[i]["delta_g_over_g0"]) + 1 for i in ids]
y = [float(pr[i]["sigma_LO_fb"]) for i in ids]
e = [float(pr[i]["sigma_MC_err_fb"]) for i in ids]
s0 = float(pr["A1"]["sigma_LO_fb"])
xs = [0.05 + k * 0.01 for k in range(150)]
ax.plot(xs, [s0 * t * t for t in xs], color="#8C8C8C", lw=1, label=r"$\sigma(g_0)\,(g/g_0)^2$ (reference only)")
ax.errorbar(x, y, yerr=e, fmt="o", color="#C44E52", ms=4, label="MG5 LO integrations")
for xi, yi, i in zip(x, y, ids):
    if i in ("A1", "B1", "B2", "B3"):
        ax.annotate(i, (xi, yi), textcoords="offset points", xytext=(5, -10), fontsize=7)
ax.set_yscale("log"); ax.set_xlabel(r"$g_{h\phi\phi}/g_0$"); ax.set_ylabel(r"$\sigma(gg\to h^*\to\phi\phi)$ [fb]")
ax.legend(fontsize=7, frameon=False); fig.tight_layout()
fig.savefig(OUT / "sigma_vs_g.pdf"); fig.savefig(OUT / "sigma_vs_g.png", dpi=150)
print("figures written to", OUT)
