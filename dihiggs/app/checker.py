#!/usr/bin/env python3
import pandas as pd
import sys

if len(sys.argv) < 2:
    print("Uso: python check_ok.py <archivo.csv>")
    sys.exit(1)

file = sys.argv[1]
df = pd.read_csv(file)

# columnas de interés
conds = ["positivity_ok", "unitarity_ok", "perturbativity_ok"]

total = len(df)
if total == 0:
    print("CSV vacío")
    sys.exit(0)

def pct(mask):
    return 100.0 * mask.sum() / total

print(f"Total puntos: {total}\n")

# porcentajes individuales
for c in conds:
    print(f"{c:20s}: {pct(df[c] == 1):6.2f}%")

# intersecciones de pares
print("\nIntersecciones pares:")
for i in range(len(conds)):
    for j in range(i+1, len(conds)):
        ci, cj = conds[i], conds[j]
        mask = (df[ci] == 1) & (df[cj] == 1)
        print(f"{ci} ∧ {cj:15s}: {pct(mask):6.2f}%")

# intersección triple
mask_all = (df["positivity_ok"] == 1) & (df["unitarity_ok"] == 1) & (df["perturbativity_ok"] == 1)
print(f"\nIntersección triple: {pct(mask_all):6.2f}%")
print(f"\n\t N=",mask_all.sum())
