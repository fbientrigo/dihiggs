import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# --- Configuración ---
# Definir rutas de archivos relativas a la estructura del proyecto
base_path = "../app/out_l6_test/validation/"
files = {
    "Simetría Exacta (lambda6=0)": "case_A_l6_zero.csv",
    "Ruptura Leve (lambda6=0.01)": "case_B_l6_point01.csv",
    "Ruptura Fuerte (lambda6=0.1)": "case_B_l6_point1.csv"
}

# Constante de conversión: GeV^-1 a mm
# hbar * c ~ 1.973e-13 mm * GeV
HBAR_C_MM_GEV = 1.973269804e-13

# --- 1. Carga y Limpieza de Datos ---
dfs = {}
print("--- Carga de Datos ---")
for label, fname in files.items():
    try:
        fpath = os.path.join(base_path, fname)
        if not os.path.exists(fpath):
            print(f"[WARN] Archivo no encontrado: {fpath}")
            continue
            
        df_raw = pd.read_csv(fpath)
        
        # Filtrar puntos teóricamente válidos (Positivity, Unitarity, Perturbativity)
        # Asumimos que las columnas son flags 1/0
        df_valid = df_raw[
            (df_raw['positivity_ok'] == 1) & 
            (df_raw['unitarity_ok'] == 1) & 
            (df_raw['perturbativity_ok'] == 1)
        ].sort_values('m_phi')
        
        if not df_valid.empty:
            dfs[label] = df_valid
            print(f"[{label}] -> {len(df_valid)} puntos físicos válidos.")
        else:
            print(f"[{label}] -> 0 puntos válidos encontrados (o archivo vacío).")
            
    except Exception as e:
        print(f"[ERROR] Fallo al leer {fname}: {e}")

# --- 2. Visualización: Impacto en Total Width ---
if dfs:
    plt.figure(figsize=(12, 10))
    
    # Subplot 1: Total Width (GeV)
    plt.subplot(2, 1, 1)
    colors = ['#2ca02c', '#ff7f0e', '#d62728'] # Green, Orange, Red
    
    for (label, df), color in zip(dfs.items(), colors):
        plt.plot(df['m_phi'], df['total_width'], 'o-', 
                 label=label, color=color, markersize=4, alpha=0.7, linewidth=1.5)
    
    plt.yscale('log')
    plt.title(r'Impacto de la Ruptura de Simetría en $\Gamma_{tot}$ ($H \to hh$ dominance)', fontsize=14)
    plt.ylabel(r'Total Width $\Gamma_{H}$ [GeV]', fontsize=12)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(fontsize=10)

    # Subplot 2: Lifetime Estimado (c*tau)
    plt.subplot(2, 1, 2)
    
    for (label, df), color in zip(dfs.items(), colors):
        # Calcular lifetime: c*tau = (hbar*c) / Gamma
        # Evitar división por cero si width es 0
        width_safe = df['total_width'].replace(0, np.nan)
        ctau = HBAR_C_MM_GEV / width_safe
        
        plt.plot(df['m_phi'], ctau, 'o-', 
                 label=label, color=color, markersize=4, alpha=0.7, linewidth=1.5)

    plt.yscale('log')
    plt.xlabel(r'Masa del Higgs Pesado $m_{\phi}$ [GeV]', fontsize=12)
    plt.ylabel(r'Lifetime $c\tau$ [mm]', fontsize=12)
    plt.title(r'Colapso del Lifetime ($c\tau$)', fontsize=14)
    
    # Línea de referencia para detectores (1 mm es visible)
    plt.axhline(y=1.0, color='gray', linestyle=':', linewidth=2, label='Umbral de Detección (1 mm)')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(fontsize=10)
    
    plt.tight_layout()
    plt.savefig('validation_1203_plot.png')
    print("Gráfica guardada en validation_1203_plot.png")

    # --- 3. Reporte Estadístico ---
    print("\n--- Resumen Estadístico del Width (GeV) ---")
    stats = []
    for label, df in dfs.items():
        w = df['total_width']
        stats.append({
            "Escenario": label,
            "Min Width": f"{w.min():.4e}",
            "Max Width": f"{w.max():.4e}",
            "Mean Width": f"{w.mean():.4e}",
            "Dominant Channel": "H -> hh (Expected)" # Placeholder teórico
        })
    
    stat_df = pd.DataFrame(stats)
    print(stat_df.to_markdown(index=False))

else:
    print("No hay datos suficientes para graficar.")
