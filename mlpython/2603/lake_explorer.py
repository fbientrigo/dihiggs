import pandas as pd
import glob
import os
import time
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import warnings

# Ignorar warnings menores de matplotlib
warnings.filterwarnings("ignore")

# --- 1. CONFIGURACIÓN DEL DATA LAKE ---
DATA_LAKE_DIR = "/mnt/c/Users/Asus/cern_db/dihiggs_lake"

# Columnas estrictamente necesarias para la física (ahorra ~70% de RAM)
TARGET_COLS = [
    'm_phi', 'lam1', 'tan_beta', 'lambda6', 
    'positivity_ok', 'unitarity_ok', 'perturbativity_ok', 
    'br_gaga', 'total_width'
]

def load_csv_batch(file_list):
    """Carga un lote de CSVs filtrando solo las columnas de interés físico."""
    dfs = []
    for f in file_list:
        try:
            # Usamos usecols pasándole una lambda para que sea robusto si alguna columna falta
            df = pd.read_csv(f, usecols=lambda c: c in TARGET_COLS)
            if not df.empty:
                dfs.append(df)
        except Exception:
            pass
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return pd.DataFrame()

def format_time(seconds):
    if seconds < 60: return f"{seconds:.2f} seg"
    elif seconds < 3600: return f"{seconds/60:.2f} min"
    else: return f"{seconds/3600:.2f} hrs"

# ==============================================================================
# FASE 1: DESCUBRIMIENTO Y PROFILING (ESTIMACIÓN LINEAL)
# ==============================================================================
print(f"[*] Escaneando el Data Lake en: {DATA_LAKE_DIR}")
all_csv_files = glob.glob(os.path.join(DATA_LAKE_DIR, "**/*.csv"), recursive=True)
total_files = len(all_csv_files)
print(f"[+] Se encontraron {total_files:,} archivos CSV en total.\n")

if total_files == 0:
    print("[-] No hay datos para analizar. Saliendo.")
    exit()

print("[-] Iniciando Profiling de I/O...")
# Prueba con 10 archivos
t0 = time.perf_counter()
df_10 = load_csv_batch(all_csv_files[:10])
t_10 = time.perf_counter() - t0
pts_10 = len(df_10)

# Prueba con 100 archivos (si hay suficientes)
sample_size = min(100, total_files)
t0 = time.perf_counter()
df_100 = load_csv_batch(all_csv_files[:sample_size])
t_100 = time.perf_counter() - t0
pts_100 = len(df_100)

# Estimación lineal O(N) basada en el sample más grande
estimated_total_time = (t_100 / sample_size) * total_files
estimated_total_pts = int((pts_100 / sample_size) * total_files)

print("\n" + "="*50)
print(" ⏱️ ESTIMACIÓN DE CARGA (Escalamiento Lineal)")
print("="*50)
print(f" -> Tiempo lectura 10 archivos  : {t_10:.4f} seg ({pts_10:,} puntos)")
print(f" -> Tiempo lectura {sample_size} archivos : {t_100:.4f} seg ({pts_100:,} puntos)")
print(f" -> ESTIMADO PARA {total_files:,} archivos : {format_time(estimated_total_time)}")
print(f" -> Puntos totales proyectados : ~{estimated_total_pts:,}")
print("="*50 + "\n")

# Breve pausa para que el usuario lea
time.sleep(3)

# ==============================================================================
# FASE 2: EXTRACCIÓN MASIVA DE DATOS
# ==============================================================================
print(f"[*] Ejecutando extracción masiva de datos (esto tomará ~{format_time(estimated_total_time)})...")
t_start_massive = time.perf_counter()

df = load_csv_batch(all_csv_files)

t_end_massive = time.perf_counter() - t_start_massive
print(f"[+] Carga completada en {format_time(t_end_massive)}. Rendimiento real: {len(df):,} puntos.")

# Calcular el flag "Triple OK"
df['is_valid'] = (df['positivity_ok'] == 1) & (df['unitarity_ok'] == 1) & (df['perturbativity_ok'] == 1)
df_valid = df[df['is_valid']].copy()
# ==============================================================================
# FASE 3: ESTADÍSTICAS DETALLADAS
# ==============================================================================
total_pts = len(df)
valid_pts = len(df_valid)

print("\n" + "="*50)
print(" 📊 ESTADÍSTICAS GLOBALES DEL DATA LAKE")
print("="*50)
print(f"Puntos Totales Computados : {total_pts:,}")
print(f"Puntos Físicamente Válidos: {valid_pts:,}")
print(f"Eficiencia de Muestreo    : {(valid_pts/total_pts)*100:.3f}%")

print("\n[+] Desglose de fallos en puntos descartados:")
print(f" -> Fallaron Positividad   : {(df['positivity_ok'] == 0).sum():,}")
print(f" -> Fallaron Unitariedad   : {(df['unitarity_ok'] == 0).sum():,}")
print(f" -> Fallaron Perturbatividad: {(df['perturbativity_ok'] == 0).sum():,}")

# --- FIX: Corrección del typo value_counts() ---
print("\n[+] Top 5 Valores de Lambda_6 con más puntos válidos:")
print(df_valid['lambda6'].value_counts().head(5))

print("\n[+] Top 5 Valores de tan(beta) con más puntos válidos:")
print(df_valid['tan_beta'].value_counts().head(5))
print("="*50 + "\n")

# Liberar memoria de los puntos inválidos que ya no usaremos
del df
import gc
gc.collect()

# ==============================================================================
# FASE 4: VISUALIZACIÓN Y GUARDADO DE ALTA CALIDAD (Para la Reunión)
# ==============================================================================
if valid_pts == 0:
    print("[-] No hay suficientes puntos válidos para graficar.")
    exit()

print("[*] Generando y guardando gráficos de alta resolución...")

# Asegurar que el directorio de imágenes del paper exista
IMG_DIR = "/home/fabi/dihiggs/mlpython/paper_img/reunion_marzo"
os.makedirs(IMG_DIR, exist_ok=True)

sns.set_theme(style="whitegrid", palette="muted")

# 1. Mapa de Densidad del Espacio de Fase (m_phi vs lam1)
plt.figure(figsize=(10, 6))
hb = plt.hexbin(df_valid['lam1'], df_valid['m_phi'], gridsize=40, cmap='YlGnBu', bins='log', mincnt=1)
cb = plt.colorbar(hb, label='log10(N puntos válidos)')
plt.title(r'Distribución de Puntos Válidos: $\lambda_1$ vs $m_\phi$')
plt.xlabel(r'$\lambda_1$')
plt.ylabel(r'$m_\phi$ [GeV]')
plt.tight_layout()
out_path1 = os.path.join(IMG_DIR, "1_espacio_fase_hexbin.png")
plt.savefig(out_path1, dpi=300)
print(f" [+] Guardado: {out_path1}")
plt.close()

# 2. Impacto de Lambda_6 sobre el espacio permitido de Lambda_1
plt.figure(figsize=(12, 6))
sns.violinplot(data=df_valid, x='lambda6', y='lam1', inner="quartile", palette="Set2")
plt.title(r"Rango permitido de $\lambda_1$ condicionado por $\lambda_6$")
plt.xlabel(r"$\lambda_6$")
plt.ylabel(r"$\lambda_1$ (Válidos)")
plt.xticks(rotation=45)
plt.tight_layout()
out_path2 = os.path.join(IMG_DIR, "2_simetria_z2_violin.png")
plt.savefig(out_path2, dpi=300)
print(f" [+] Guardado: {out_path2}")
plt.close()

# 3. Total Width vs BR_gaga (Búsqueda de Vértices Desplazados)
plt.figure(figsize=(10, 6))
df_br = df_valid[df_valid['br_gaga'] > 1e-15]
sns.scatterplot(data=df_br, x='total_width', y='br_gaga', alpha=0.3, s=10, edgecolor=None, color='indigo')
plt.xscale('log')
plt.yscale('log')
plt.title(r"Correlación: Anchura Total vs Branching Ratio a $\gamma\gamma$")
plt.xlabel(r"$\Gamma_{total}$ [GeV]")
plt.ylabel(r"BR($h \to \gamma\gamma$)")
plt.axvline(x=1e-10, color='red', linestyle='--', label=r'Región Long-Lived ($\Gamma < 10^{-10}$ GeV)')
plt.legend()
plt.tight_layout()
out_path3 = os.path.join(IMG_DIR, "3_long_lived_particles.png")
plt.savefig(out_path3, dpi=300)
print(f" [+] Guardado: {out_path3}")
plt.close()

print("\n[🚀] ¡Todos los gráficos han sido exportados exitosamente para la presentación!")