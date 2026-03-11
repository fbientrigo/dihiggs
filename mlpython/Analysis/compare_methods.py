import argparse
import pandas as pd
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    # --- 1. CONFIGURACIÓN MEDIANTE ARGPARSE ---
    parser = argparse.ArgumentParser(description="Compara el rendimiento de exploración: Bayesiano vs Clásico.")
    parser.add_argument(
        "-d", "--datalake-dir", 
        type=str, 
        default="/mnt/c/Users/Asus/cern_db/dihiggs_lake",
        help="Ruta al directorio raíz del Data Lake (por defecto asume entorno WSL local)."
    )
    args = parser.parse_args()

    DATA_LAKE_DIR = args.datalake_dir

    # Columnas de interés
    TARGET_COLS = [
        'm_phi', 'lam1', 'tan_beta', 'lambda6', 
        'positivity_ok', 'unitarity_ok', 'perturbativity_ok', 
        'br_gaga', 'total_width'
    ]

    print(f"[*] Escaneando Data Lake en: {DATA_LAKE_DIR}")
    print(f"[*] Buscando el caso de control: l6=0.1, mA=300...")

    # Buscamos TODOS los CSVs
    all_files = glob.glob(os.path.join(DATA_LAKE_DIR, "**/*.csv"), recursive=True)

    # Filtramos estrictamente por la nomenclatura de carpetas de tu C++
    target_files = [f for f in all_files if 'l6=0p1000' in f and 'mA=300p0' in f]

    print(f"[+] Encontrados {len(target_files)} archivos que coinciden con el caso de control.")

    if not target_files:
        print("[-] No se encontraron datos. Revisa la ruta especificada.")
        return

    # --- 2. EXTRACCIÓN Y CLASIFICACIÓN ---
    df_list = []
    for f in target_files:
        try:
            temp_df = pd.read_csv(f, usecols=lambda c: c in TARGET_COLS)
            if temp_df.empty: continue
                
            # Clasificador Heurístico de Origen
            path_lower = f.lower()
            if 'adaptive' in path_lower or 'exploracion' in path_lower or 'noche' in path_lower:
                temp_df['Metodo'] = 'Bayesiano (Adaptativo)'
            else:
                temp_df['Metodo'] = 'Clásico (Grid/Random)'
                
            df_list.append(temp_df)
        except Exception:
            pass

    df = pd.concat(df_list, ignore_index=True)
    df['is_valid'] = (df['positivity_ok'] == 1) & (df['unitarity_ok'] == 1) & (df['perturbativity_ok'] == 1)

    # --- 3. ANÁLISIS DE RENDIMIENTO (EFICIENCIA) ---
    stats = df.groupby('Metodo')['is_valid'].agg(['count', 'sum'])
    stats.rename(columns={'count': 'Puntos Simulados', 'sum': 'Puntos Válidos'}, inplace=True)
    stats['Eficiencia (%)'] = (stats['Puntos Válidos'] / stats['Puntos Simulados']) * 100

    print("\n" + "="*60)
    print(" 🏆 COMPARATIVA DE RENDIMIENTO COMPUTACIONAL")
    print("="*60)
    print(stats.to_string())
    print("="*60 + "\n")

    df_valid = df[df['is_valid']].copy()

    if df_valid.empty:
        print("[-] No hay puntos Triple-OK para graficar.")
        return

    # --- 4. GRAFICACIÓN COMPARATIVA DE ALTA CALIDAD ---
    # Usaremos rutas relativas al repositorio para que sea agnóstico al PC
    IMG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper_img", "reunion_marzo")
    os.makedirs(IMG_DIR, exist_ok=True)
    sns.set_theme(style="whitegrid", palette="Set1")

    # GRAFICO A: Eficiencia Computacional (Barplot)
    plt.figure(figsize=(8, 6))
    ax = sns.barplot(x=stats.index, y=stats['Eficiencia (%)'], palette=['#e74c3c', '#3498db'])
    # Se añade la 'r' antes de la string para escapar caracteres LaTeX y evitar SyntaxWarnings
    plt.title(r"Eficiencia de Búsqueda: Bayesiano vs Clásico" + "\n" + r"($\lambda_6 = 0.1, m_A = 300$ GeV)", fontsize=14)
    plt.ylabel("Eficiencia (% de puntos Triple-OK)", fontsize=12)
    plt.xlabel("Método de Muestreo", fontsize=12)
    for i, v in enumerate(stats['Eficiencia (%)']):
        ax.text(i, v + 0.5, f"{v:.2f}%", ha='center', fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(IMG_DIR, "A_eficiencia_comparativa.png"), dpi=300)
    plt.close()

    # GRAFICO B: Mapeo del Espacio de Fase Permitido
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df_valid, x='lam1', y='m_phi', hue='Metodo', 
                    style='Metodo', alpha=0.5, s=25, palette=['#e74c3c', '#3498db'])
    plt.title(r"Puntos Físicamente Válidos Descubiertos ($\lambda_1$ vs $m_\phi$)", fontsize=14)
    plt.xlabel(r"$\lambda_1$", fontsize=12)
    plt.ylabel(r"$m_\phi$ [GeV]", fontsize=12)
    plt.legend(title="Método")
    plt.tight_layout()
    plt.savefig(os.path.join(IMG_DIR, "B_espacio_fase_comparativo.png"), dpi=300)
    plt.close()

    # GRAFICO C: Calidad de la Física (Distribución de Branching Ratio)
    plt.figure(figsize=(10, 6))
    sns.kdeplot(data=df_valid[df_valid['br_gaga'] > 1e-10], x='br_gaga', hue='Metodo', 
                fill=True, common_norm=False, log_scale=True, palette=['#e74c3c', '#3498db'])
    plt.title(r"Densidad de Puntos por Branching Ratio ($h \to \gamma\gamma$)", fontsize=14)
    plt.xlabel(r"BR($h \to \gamma\gamma$)", fontsize=12)
    plt.ylabel("Densidad de Modelos", fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(IMG_DIR, "C_branching_ratio_comparativo.png"), dpi=300)
    plt.close()

    print(f"[🚀] Gráficos de comparativa guardados en: {IMG_DIR}")

if __name__ == "__main__":
    main()