import os
import glob
import time
import gc
import warnings
from collections import Counter

import numpy as np
import polars as pl
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# =========================================================
# CONFIG
# =========================================================
DATA_LAKE_DIR = "/mnt/c/Users/Asus/cern_db/dihiggs_lake"
IMG_DIR = "/home/fabi/dihiggs/mlpython/paper_img/reunion_marzo_polars"

TARGET_COLS = [
    "m_phi",
    "lam1",
    "tan_beta",
    "lambda6",
    "positivity_ok",
    "unitarity_ok",
    "perturbativity_ok",
    "br_gaga",
    "total_width",
]

SCHEMA_OVERRIDES = {
    "m_phi": pl.Float32,
    "lam1": pl.Float32,
    "tan_beta": pl.Float32,
    "lambda6": pl.Float32,
    "positivity_ok": pl.Int8,
    "unitarity_ok": pl.Int8,
    "perturbativity_ok": pl.Int8,
    "br_gaga": pl.Float32,
    "total_width": pl.Float32,
}

# batches por archivo
BATCH_SIZE = 50_000

# conservar solo una muestra pequeña para plots
MAX_SAMPLE_PER_BATCH = 1500
MAX_PLOT_POINTS = 180_000
MAX_SCATTER_POINTS = 120_000

PROGRESS_EVERY = 200

RNG = np.random.default_rng(42)


def format_time(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.2f} seg"
    elif seconds < 3600:
        return f"{seconds/60:.2f} min"
    else:
        return f"{seconds/3600:.2f} hrs"


def safe_read_batches(path: str):
    """
    Devuelve un BatchedCsvReader de Polars o None si falla.
    """
    try:
        reader = pl.read_csv_batched(
            source=path,
            columns=TARGET_COLS,
            schema_overrides=SCHEMA_OVERRIDES,
            infer_schema_length=100,
            batch_size=BATCH_SIZE,
            low_memory=True,
            rechunk=False,
            ignore_errors=True,
            try_parse_dates=False,
        )
        return reader
    except Exception:
        return None


def sample_df(df: pl.DataFrame, nmax: int) -> pl.DataFrame:
    """
    Devuelve una muestra aleatoria pequeña para plots.
    """
    if df.height <= nmax:
        return df
    return df.sample(n=nmax, shuffle=True, seed=int(RNG.integers(0, 2**31 - 1)))


def update_counter_from_grouped(counter: Counter, grouped_df: pl.DataFrame, key_col: str, count_col: str = "len") -> None:
    """
    Actualiza un Counter Python desde un group_by.len() de Polars.
    """
    if grouped_df.height == 0:
        return
    for key, cnt in grouped_df.select([key_col, count_col]).iter_rows():
        counter[key] += int(cnt)


# =========================================================
# DESCUBRIMIENTO
# =========================================================
print(f"[*] Escaneando el Data Lake en: {DATA_LAKE_DIR}")
all_csv_files = sorted(glob.glob(os.path.join(DATA_LAKE_DIR, "**/*.csv"), recursive=True))
total_files = len(all_csv_files)
print(f"[+] Se encontraron {total_files:,} archivos CSV.\n")

if total_files == 0:
    raise SystemExit("[-] No se encontraron CSVs.")

os.makedirs(IMG_DIR, exist_ok=True)

# =========================================================
# ACUMULADORES GLOBALES
# =========================================================
total_rows = 0
valid_rows = 0
fail_pos = 0
fail_uni = 0
fail_pert = 0
skipped_files = 0

lambda6_counts = Counter()
tanb_counts = Counter()

plot_samples = []

# =========================================================
# EXTRACCIÓN STREAMING POR BATCHES
# =========================================================
print("[*] Iniciando extracción con Polars (batched, memoria controlada)...")
t0 = time.perf_counter()

for i, fpath in enumerate(all_csv_files, start=1):
    reader = safe_read_batches(fpath)

    if reader is None:
        skipped_files += 1
        continue

    try:
        batches = reader.next_batches(8)  # pedir varios batches a la vez amortiza overhead
    except Exception:
        skipped_files += 1
        continue

    while batches:
        for batch in batches:
            if batch is None or batch.height == 0:
                continue

            total_rows += batch.height

            # Estadísticas de fallos sobre el batch completo
            stats = batch.select(
                [
                    (pl.col("positivity_ok") == 0).sum().alias("fail_pos"),
                    (pl.col("unitarity_ok") == 0).sum().alias("fail_uni"),
                    (pl.col("perturbativity_ok") == 0).sum().alias("fail_pert"),
                ]
            )

            fail_pos += int(stats[0, "fail_pos"])
            fail_uni += int(stats[0, "fail_uni"])
            fail_pert += int(stats[0, "fail_pert"])

            # Filtro físico
            valid = batch.filter(
                (pl.col("positivity_ok") == 1)
                & (pl.col("unitarity_ok") == 1)
                & (pl.col("perturbativity_ok") == 1)
            ).select(
                ["m_phi", "lam1", "tan_beta", "lambda6", "br_gaga", "total_width"]
            )

            n_valid = valid.height
            valid_rows += n_valid

            if n_valid > 0:
                # Conteos exactos de lambda6 y tan_beta redondeados
                # Redondear evita que el ruido de float convierta un mismo valor físico en muchas claves vecinas.
                lambda6_grouped = (
                    valid.with_columns(pl.col("lambda6").round(6))
                    .group_by("lambda6")
                    .len()
                )
                tanb_grouped = (
                    valid.with_columns(pl.col("tan_beta").round(6))
                    .group_by("tan_beta")
                    .len()
                )

                update_counter_from_grouped(lambda6_counts, lambda6_grouped, "lambda6")
                update_counter_from_grouped(tanb_counts, tanb_grouped, "tan_beta")

                # Retener solo muestra pequeña
                sampled = sample_df(valid, MAX_SAMPLE_PER_BATCH)
                plot_samples.append(sampled)

                del lambda6_grouped, tanb_grouped, sampled

            # Liberar todo lo innecesario
            del stats, valid, batch
            gc.collect()

        try:
            batches = reader.next_batches(8)
        except Exception:
            break

    del reader
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
if plot_samples:
    plot_df = pl.concat(plot_samples, rechunk=False)

    if plot_df.height > MAX_PLOT_POINTS:
        plot_df = plot_df.sample(
            n=MAX_PLOT_POINTS,
            shuffle=True,
            seed=int(RNG.integers(0, 2**31 - 1)),
        )
else:
    plot_df = pl.DataFrame(
        schema={
            "m_phi": pl.Float32,
            "lam1": pl.Float32,
            "tan_beta": pl.Float32,
            "lambda6": pl.Float32,
            "br_gaga": pl.Float32,
            "total_width": pl.Float32,
        }
    )

del plot_samples
gc.collect()

# =========================================================
# RESUMEN
# =========================================================
eff = (100.0 * valid_rows / total_rows) if total_rows > 0 else 0.0

summary_lines = [
    "=" * 60,
    "ESTADÍSTICAS GLOBALES DEL DATA LAKE",
    "=" * 60,
    f"Archivos CSV encontrados      : {total_files:,}",
    f"Archivos saltados             : {skipped_files:,}",
    f"Puntos totales computados     : {total_rows:,}",
    f"Puntos físicamente válidos    : {valid_rows:,}",
    f"Eficiencia de muestreo        : {eff:.4f}%",
    "",
    "Desglose de fallos:",
    f" - Fallaron positividad       : {fail_pos:,}",
    f" - Fallaron unitariedad       : {fail_uni:,}",
    f" - Fallaron perturbatividad   : {fail_pert:,}",
    "",
    "Top 5 lambda6 con más puntos válidos:",
]

for k, v in lambda6_counts.most_common(5):
    summary_lines.append(f" - lambda6 = {k}: {v:,}")

summary_lines.append("")
summary_lines.append("Top 5 tan_beta con más puntos válidos:")
for k, v in tanb_counts.most_common(5):
    summary_lines.append(f" - tan_beta = {k}: {v:,}")

summary_lines.extend(
    [
        "",
        f"Muestra retenida para plots   : {plot_df.height:,}",
        f"Tiempo total                  : {format_time(elapsed_total)}",
        "=" * 60,
    ]
)

summary_text = "\n".join(summary_lines)
print("\n" + summary_text + "\n")

with open(os.path.join(IMG_DIR, "summary.txt"), "w", encoding="utf-8") as f:
    f.write(summary_text)

# =========================================================
# PLOTS
# =========================================================
if plot_df.height == 0:
    print("[-] No hubo puntos válidos para graficar.")
    raise SystemExit(0)

# Convertimos solo la muestra final a pandas para matplotlib
pdf = plot_df.to_pandas()
del plot_df
gc.collect()

print("[*] Generando gráficos...")

# 1) Hexbin lam1 vs m_phi
plt.figure(figsize=(10, 6))
hb = plt.hexbin(
    pdf["lam1"],
    pdf["m_phi"],
    gridsize=40,
    bins="log",
    mincnt=1,
    cmap="YlGnBu",
)
plt.colorbar(hb, label="log10(N puntos en muestra)")
plt.title(r"Distribución de puntos válidos (muestra): $\lambda_1$ vs $m_\phi$")
plt.xlabel(r"$\lambda_1$")
plt.ylabel(r"$m_\phi$ [GeV]")
plt.tight_layout()
out1 = os.path.join(IMG_DIR, "1_espacio_fase_hexbin.png")
plt.savefig(out1, dpi=250)
plt.close()
print(f" [+] Guardado: {out1}")

# 2) Violin de lam1 para top 10 lambda6 más poblados
top_lambda6 = [k for k, _ in lambda6_counts.most_common(10)]
vdf = pdf[np.isin(np.round(pdf["lambda6"], 6), top_lambda6)].copy()

groups = []
labels = []

for lv in sorted(np.unique(np.round(vdf["lambda6"].dropna(), 6))):
    arr = vdf.loc[np.round(vdf["lambda6"], 6) == lv, "lam1"].dropna().to_numpy()
    if len(arr) > 4000:
        arr = RNG.choice(arr, size=4000, replace=False)
    if len(arr) >= 5:
        groups.append(arr)
        labels.append(f"{lv:g}")

if groups:
    plt.figure(figsize=(max(10, 0.9 * len(labels) + 2), 6))
    plt.violinplot(groups, showmedians=True, showextrema=False)
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
    print(" [!] Plot 2 omitido: no hubo suficientes grupos por lambda6.")

# 3) total_width vs br_gaga
df_br = pdf[pdf["br_gaga"] > 1e-15].copy()

if len(df_br) > MAX_SCATTER_POINTS:
    idx = RNG.choice(len(df_br), size=MAX_SCATTER_POINTS, replace=False)
    df_br = df_br.iloc[idx].copy()

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