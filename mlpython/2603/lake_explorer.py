import os
import glob
import time
import gc
import warnings
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# =========================================================
# CONFIG
# =========================================================
DATA_LAKE_DIR = "/mnt/c/Users/Asus/cern_db/dihiggs_lake"
IMG_DIR = "/home/fabi/dihiggs/mlpython/paper_img/reunion_marzo_streaming"

TARGET_COLS = [
    "m_phi", "lam1", "tan_beta", "lambda6",
    "positivity_ok", "unitarity_ok", "perturbativity_ok",
    "br_gaga", "total_width"
]

# Tipos compactos para reducir RAM.
# float32 basta para plots/estadísticas rápidas de reunión.
DTYPE_MAP = {
    "m_phi": "float32",
    "lam1": "float32",
    "tan_beta": "float32",
    "lambda6": "float32",
    "positivity_ok": "Int8",
    "unitarity_ok": "Int8",
    "perturbativity_ok": "Int8",
    "br_gaga": "float32",
    "total_width": "float32",
}

# Muestra máxima guardada por archivo para plots
MAX_VALID_SAMPLE_PER_FILE = 100

# Tope global de puntos usados para plots
MAX_PLOT_POINTS = 350_000

# Para el scatter de BR vs width
MAX_SCATTER_POINTS = 120_000

# Progreso
PROGRESS_EVERY = 250

RNG = np.random.default_rng(42)


def format_time(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.2f} seg"
    elif seconds < 3600:
        return f"{seconds/60:.2f} min"
    else:
        return f"{seconds/3600:.2f} hrs"


def read_csv_fast(path: str) -> pd.DataFrame | None:
    """
    Lee un CSV cargando solo las columnas útiles.
    Devuelve None si falla.
    """
    try:
        df = pd.read_csv(
            path,
            usecols=lambda c: c in TARGET_COLS,
            dtype=DTYPE_MAP,
            engine="c",
        )
        if df.empty:
            return None

        # Asegurarse de que estén las columnas mínimas necesarias
        needed = {
            "m_phi", "lam1", "tan_beta", "lambda6",
            "positivity_ok", "unitarity_ok", "perturbativity_ok",
            "br_gaga", "total_width"
        }
        if not needed.issubset(df.columns):
            return None

        return df

    except Exception:
        return None


def compact_sample(valid_df: pd.DataFrame, nmax: int) -> pd.DataFrame:
    """
    Guarda solo una muestra pequeña por archivo para no explotar RAM.
    """
    n = len(valid_df)
    if n == 0:
        return valid_df.iloc[:0].copy()
    if n <= nmax:
        return valid_df.copy()
    seed = int(RNG.integers(0, 2**31 - 1))
    return valid_df.sample(n=nmax, random_state=seed).copy()


def save_summary_txt(path: str, text: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


# =========================================================
# DESCUBRIMIENTO
# =========================================================
print(f"[*] Escaneando el Data Lake en: {DATA_LAKE_DIR}")
all_csv_files = sorted(glob.glob(os.path.join(DATA_LAKE_DIR, "**/*.csv"), recursive=True))
total_files = len(all_csv_files)
print(f"[+] Se encontraron {total_files:,} archivos CSV.\n")

if total_files == 0:
    raise SystemExit("[-] No hay archivos CSV.")

os.makedirs(IMG_DIR, exist_ok=True)

# =========================================================
# EXTRACCIÓN STREAMING
# =========================================================
print("[*] Iniciando extracción STREAMING (sin cargar todo a RAM)...")
t0 = time.perf_counter()

total_rows = 0
valid_rows = 0
fail_pos = 0
fail_uni = 0
fail_pert = 0
skipped_files = 0

lambda6_counts = Counter()
tanb_counts = Counter()

sampled_valid_parts = []

for i, fpath in enumerate(all_csv_files, start=1):
    df = read_csv_fast(fpath)

    if df is None:
        skipped_files += 1
        continue

    n_rows = len(df)
    total_rows += n_rows

    # Fallos globales
    fail_pos += int((df["positivity_ok"] == 0).sum())
    fail_uni += int((df["unitarity_ok"] == 0).sum())
    fail_pert += int((df["perturbativity_ok"] == 0).sum())

    # Mask de validez física
    valid_mask = (
        (df["positivity_ok"] == 1) &
        (df["unitarity_ok"] == 1) &
        (df["perturbativity_ok"] == 1)
    )

    valid_df = df.loc[
        valid_mask,
        ["m_phi", "lam1", "tan_beta", "lambda6", "br_gaga", "total_width"]
    ].copy()

    n_valid = len(valid_df)
    valid_rows += n_valid

    if n_valid > 0:
        # Contadores exactos, redondeando para evitar ruido de floats
        lambda6_counts.update(valid_df["lambda6"].round(6).value_counts().to_dict())
        tanb_counts.update(valid_df["tan_beta"].round(6).value_counts().to_dict())

        # Guardamos solo una muestra pequeña para plots
        sampled_valid_parts.append(compact_sample(valid_df, MAX_VALID_SAMPLE_PER_FILE))

    del df, valid_df, valid_mask
    gc.collect()

    if (i % PROGRESS_EVERY == 0) or (i == total_files):
        elapsed = time.perf_counter() - t0
        rate = i / elapsed if elapsed > 0 else 0.0
        eta = (total_files - i) / rate if rate > 0 else float("inf")
        print(
            f"[{i:>5,}/{total_files:,}] "
            f"elapsed={format_time(elapsed)} | "
            f"ETA={format_time(eta)} | "
            f"rows={total_rows:,} | valid={valid_rows:,} | skipped={skipped_files:,}"
        )

elapsed_total = time.perf_counter() - t0
print(f"\n[+] Extracción completada en {format_time(elapsed_total)}.")

# =========================================================
# CONSOLIDAR SOLO LA MUESTRA PARA PLOTS
# =========================================================
if sampled_valid_parts:
    sampled_valid = pd.concat(sampled_valid_parts, ignore_index=True)
else:
    sampled_valid = pd.DataFrame(columns=["m_phi", "lam1", "tan_beta", "lambda6", "br_gaga", "total_width"])

del sampled_valid_parts
gc.collect()

if len(sampled_valid) > MAX_PLOT_POINTS:
    seed = int(RNG.integers(0, 2**31 - 1))
    sampled_valid = sampled_valid.sample(n=MAX_PLOT_POINTS, random_state=seed).copy()

# =========================================================
# ESTADÍSTICAS
# =========================================================
eff = (100.0 * valid_rows / total_rows) if total_rows > 0 else 0.0

summary = []
summary.append("=" * 60)
summary.append("ESTADÍSTICAS GLOBALES DEL DATA LAKE")
summary.append("=" * 60)
summary.append(f"Archivos CSV encontrados      : {total_files:,}")
summary.append(f"Archivos saltados             : {skipped_files:,}")
summary.append(f"Puntos totales computados     : {total_rows:,}")
summary.append(f"Puntos físicamente válidos    : {valid_rows:,}")
summary.append(f"Eficiencia de muestreo        : {eff:.4f}%")
summary.append("")
summary.append("Desglose de fallos:")
summary.append(f" - Fallaron positividad       : {fail_pos:,}")
summary.append(f" - Fallaron unitariedad       : {fail_uni:,}")
summary.append(f" - Fallaron perturbatividad   : {fail_pert:,}")
summary.append("")
summary.append("Top 5 lambda6 con más puntos válidos:")
for k, v in lambda6_counts.most_common(5):
    summary.append(f" - lambda6 = {k}: {v:,}")
summary.append("")
summary.append("Top 5 tan_beta con más puntos válidos:")
for k, v in tanb_counts.most_common(5):
    summary.append(f" - tan_beta = {k}: {v:,}")
summary.append("")
summary.append(f"Muestra retenida para plots   : {len(sampled_valid):,}")
summary.append(f"Tiempo total                  : {format_time(elapsed_total)}")
summary.append("=" * 60)

summary_text = "\n".join(summary)
print("\n" + summary_text + "\n")

save_summary_txt(os.path.join(IMG_DIR, "summary.txt"), summary_text)

# =========================================================
# PLOTS
# =========================================================
if len(sampled_valid) == 0:
    print("[-] No hay puntos válidos para graficar.")
    raise SystemExit(0)

print("[*] Generando gráficos...")

# ---------- Plot 1: hexbin lam1 vs m_phi ----------
plt.figure(figsize=(10, 6))
hb = plt.hexbin(
    sampled_valid["lam1"],
    sampled_valid["m_phi"],
    gridsize=40,
    bins="log",
    mincnt=1,
    cmap="YlGnBu",
)
plt.colorbar(hb, label="log10(N puntos de la muestra)")
plt.title(r"Distribución de puntos válidos (muestra): $\lambda_1$ vs $m_\phi$")
plt.xlabel(r"$\lambda_1$")
plt.ylabel(r"$m_\phi$ [GeV]")
plt.tight_layout()
out1 = os.path.join(IMG_DIR, "1_espacio_fase_hexbin.png")
plt.savefig(out1, dpi=250)
plt.close()
print(f" [+] Guardado: {out1}")

# ---------- Plot 2: violin/box simplificado para top lambda6 ----------
# Usamos solo los lambda6 más poblados para no hacer un monstruo ilegible.
top_lambda6 = [k for k, _ in lambda6_counts.most_common(10)]

vdf = sampled_valid.copy()
vdf["lambda6_round"] = vdf["lambda6"].round(6)
vdf = vdf[vdf["lambda6_round"].isin(top_lambda6)].copy()

groups = []
labels = []

for lv in sorted(vdf["lambda6_round"].dropna().unique()):
    arr = vdf.loc[vdf["lambda6_round"] == lv, "lam1"].dropna().to_numpy()
    if len(arr) > 4000:
        arr = RNG.choice(arr, size=4000, replace=False)
    if len(arr) >= 5:
        groups.append(arr)
        labels.append(f"{lv:g}")

if groups:
    plt.figure(figsize=(max(10, 0.9 * len(labels) + 2), 6))
    parts = plt.violinplot(groups, showmedians=True, showextrema=False)
    plt.xticks(np.arange(1, len(labels) + 1), labels, rotation=45)
    plt.title(r"Rango permitido de $\lambda_1$ condicionado por los $\lambda_6$ más poblados")
    plt.xlabel(r"$\lambda_6$")
    plt.ylabel(r"$\lambda_1$ (válidos, muestra)")
    plt.tight_layout()
    out2 = os.path.join(IMG_DIR, "2_lambda6_vs_lam1_violin_top10.png")
    plt.savefig(out2, dpi=250)
    plt.close()
    print(f" [+] Guardado: {out2}")
else:
    print(" [!] Plot 2 omitido: no hubo suficientes datos agrupables por lambda6.")

# ---------- Plot 3: total_width vs br_gaga ----------
df_br = sampled_valid[sampled_valid["br_gaga"] > 1e-15].copy()

if len(df_br) > MAX_SCATTER_POINTS:
    seed = int(RNG.integers(0, 2**31 - 1))
    df_br = df_br.sample(n=MAX_SCATTER_POINTS, random_state=seed).copy()

if len(df_br) > 0:
    plt.figure(figsize=(10, 6))
    plt.scatter(df_br["total_width"], df_br["br_gaga"], s=8, alpha=0.25)
    plt.xscale("log")
    plt.yscale("log")
    plt.axvline(x=1e-10, linestyle="--", color="red", label=r"Región long-lived ($\Gamma < 10^{-10}$ GeV)")
    plt.title(r"Correlación: anchura total vs BR($h \to \gamma\gamma$)")
    plt.xlabel(r"$\Gamma_{total}$ [GeV]")
    plt.ylabel(r"BR($h \to \gamma\gamma$)")
    plt.legend()
    plt.tight_layout()
    out3 = os.path.join(IMG_DIR, "3_long_lived_particles.png")
    plt.savefig(out3, dpi=250)
    plt.close()
    print(f" [+] Guardado: {out3}")
else:
    print(" [!] Plot 3 omitido: no hubo puntos con br_gaga > 1e-15.")

print("\n[🚀] Listo. Estadísticas y gráficos exportados.")
print(f"[*] Revisa: {IMG_DIR}")