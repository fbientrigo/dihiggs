#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2
import sys
import os

if len(sys.argv) < 2:
    print("Uso: python analyze_hbhs.py <csv_file>")
    sys.exit(1)

fname = sys.argv[1]
df = pd.read_csv(fname)

# Calcular chi2/NDOF y p-value
df["chi2_ndof"] = df["hs_chi2"] / df["hs_ndof"]
df["p_value"] = 1 - chi2.cdf(df["hs_chi2"], df["hs_ndof"])

print("\n=== Resumen general ===")
print(f"Total puntos: {len(df)}")
print(f"Puntos permitidos por HB: {df['hb_allowed'].sum()} / {len(df)}")

print("\n=== Ejemplos de mejores puntos (ordenados por chi2/ndof) ===")
print(df.sort_values("chi2_ndof").head(10)[
    ["m_phi","tan_beta","hs_chi2","hs_ndof","chi2_ndof","p_value","hb_allowed"]
])

# Crear carpeta para plots
plot_dir = "plots"
os.makedirs(plot_dir, exist_ok=True)

# Histograma de chi2/NDOF
plt.figure()
plt.hist(df["chi2_ndof"], bins=30, alpha=0.7, color="steelblue", edgecolor="black")
plt.xlabel(r"$\chi^2/\mathrm{NDOF}$")
plt.ylabel("Número de puntos")
plt.title(r"Distribución de $\chi^2$/NDOF en el scan")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, "hist_chi2_ndof.png"), dpi=150)
plt.close()

# Scatter m_phi vs chi2/NDOF
plt.figure()
sc = plt.scatter(df["m_phi"], df["chi2_ndof"], c=df["hb_allowed"], cmap="coolwarm", s=20)
plt.xlabel(r"$m_\phi$ [GeV]")
plt.ylabel(r"$\chi^2/\mathrm{NDOF}$")
plt.title(r"Ajuste HiggsSignals vs $m_\phi$")
plt.colorbar(sc, label="HB allowed (1=yes, 0=no)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, "scatter_mphi_chi2ndof.png"), dpi=150)
plt.close()

print(f"\nImágenes guardadas en la carpeta '{plot_dir}'")
