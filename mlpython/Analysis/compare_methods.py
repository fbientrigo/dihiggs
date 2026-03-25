import argparse
import pandas as pd
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="Compara el rendimiento de exploración: Bayesiano vs Clásico.")
    parser.add_argument("-d", "--datalake-dir", type=str, default="/mnt/c/Users/Asus/cern_db/dihiggs_lake")
    args = parser.parse_args()

    DATA_LAKE_DIR = args.datalake_dir
    TARGET_COLS = ['m_phi', 'lam1', 'tan_beta', 'lambda6', 'positivity_ok', 'unitarity_ok', 'perturbativity_ok', 'br_gaga', 'total_width']

    print(f"\n[*] Escaneando Data Lake en: {DATA_LAKE_DIR}")
    print(f"[*] Buscando el caso de control: l6=0.1, mA=300...")

    all_files = glob.glob(os.path.join(DATA_LAKE_DIR, "**/*.csv"), recursive=True)
    target_files = [f for f in all_files if 'l6=0p1000' in f and 'mA=300p0' in f]

    print(f"[+] Encontrados {len(target_files)} archivos que coinciden con el caso de control.\n")

    if not target_files:
        print("[-] No se encontraron datos.")
        return

    # --- 1. EXTRACCIÓN Y TRANSPARENCIA ---
    df_list = []
    
    # Diccionarios para auditoría transparente
    audit_stats = {
        'Bayesiano (Adaptativo)': {'archivos': 0, 'simulados': 0, 'validos': 0, 'ejemplos': []},
        'Clásico (Grid/Random)': {'archivos': 0, 'simulados': 0, 'validos': 0, 'ejemplos': []}
    }

    print("🔍 AUDITORÍA DE ARCHIVOS (Leyendo datos...)\n" + "-"*60)
    for f in target_files:
        try:
            # Lógica de Clasificación
            path_lower = f.lower()
            if 'adaptive' in path_lower or 'exploracion' in path_lower or 'noche' in path_lower:
                metodo = 'Bayesiano (Adaptativo)'
            else:
                metodo = 'Clásico (Grid/Random)'
            
            temp_df = pd.read_csv(f, usecols=lambda c: c in TARGET_COLS)
            if temp_df.empty: continue
            
            # Calcular validez (Triple-OK)
            is_valid = (temp_df['positivity_ok'] == 1) & (temp_df['unitarity_ok'] == 1) & (temp_df['perturbativity_ok'] == 1)
            pts_simulados = len(temp_df)
            pts_validos = is_valid.sum()

            # Llenar dataframes
            temp_df['Metodo'] = metodo
            temp_df['is_valid'] = is_valid
            df_list.append(temp_df)

            # Llenar auditoría
            audit_stats[metodo]['archivos'] += 1
            audit_stats[metodo]['simulados'] += pts_simulados
            audit_stats[metodo]['validos'] += pts_validos
            
            # Guardar el nombre de la carpeta contenedora para transparencia (no el archivo completo)
            parent_dir = os.path.basename(os.path.dirname(os.path.dirname(f)))
            if parent_dir not in audit_stats[metodo]['ejemplos'] and len(audit_stats[metodo]['ejemplos']) < 5:
                audit_stats[metodo]['ejemplos'].append(parent_dir)

        except Exception as e:
            pass

    # Imprimir Reporte de Transparencia
    for m, stats in audit_stats.items():
        print(f"📌 {m.upper()}:")
        print(f"   - Archivos leídos: {stats['archivos']}")
        print(f"   - Puntos Simulados: {stats['simulados']:,}")
        print(f"   - Puntos Válidos encontrados: {stats['validos']:,}")
        print(f"   - Campañas detectadas (Muestra): {', '.join(stats['ejemplos'])}")
        print()

    df = pd.concat(df_list, ignore_index=True)

    # --- 2. ANÁLISIS DE RENDIMIENTO ---
    stats = df.groupby('Metodo')['is_valid'].agg(['count', 'sum'])
    stats.rename(columns={'count': 'Puntos Simulados', 'sum': 'Puntos Válidos'}, inplace=True)
    stats['Eficiencia (%)'] = (stats['Puntos Válidos'] / stats['Puntos Simulados']) * 100

    print("="*60)
    print(" 🏆 COMPARATIVA DE RENDIMIENTO COMPUTACIONAL")
    print("="*60)
    print(stats.to_string())
    print("="*60 + "\n")

    df_valid = df[df['is_valid']].copy()

    # --- FIX WARNINGS: Diccionario de colores explícito ---
    color_map = {'Bayesiano (Adaptativo)': '#e74c3c', 'Clásico (Grid/Random)': '#3498db'}

    # --- 3. GRAFICACIÓN COMPARATIVA DE ALTA CALIDAD ---
    IMG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper_img", "reunion_marzo")
    os.makedirs(IMG_DIR, exist_ok=True)
    sns.set_theme(style="whitegrid")

    # GRAFICO A: Eficiencia Computacional
    plt.figure(figsize=(8, 6))
    # FIX: Se añade hue=stats.index y legend=False para apagar el warning de Seaborn
    ax = sns.barplot(x=stats.index, y=stats['Eficiencia (%)'], hue=stats.index, palette=color_map, legend=False)
    plt.title(r"Eficiencia de Búsqueda: Bayesiano vs Clásico" + "\n" + r"($\lambda_6 = 0.1, m_A = 300$ GeV)", fontsize=14)
    plt.ylabel("Eficiencia (% de puntos Triple-OK)", fontsize=12)
    plt.xlabel("Método de Muestreo", fontsize=12)
    for i, v in enumerate(stats['Eficiencia (%)']):
        ax.text(i, v + 0.5, f"{v:.2f}%", ha='center', fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(IMG_DIR, "A_eficiencia_comparativa.png"), dpi=300)
    plt.close()

    if not df_valid.empty:
        # GRAFICO B: Mapeo del Espacio de Fase Permitido
        plt.figure(figsize=(10, 6))
        # FIX: palette recibe el diccionario exacto, rasterized=True previene que el plot sea muy pesado
        sns.scatterplot(data=df_valid, x='lam1', y='m_phi', hue='Metodo', style='Metodo', 
                        alpha=0.5, s=25, palette=color_map, rasterized=True)
        plt.title(r"Puntos Físicamente Válidos Descubiertos ($\lambda_1$ vs $m_\phi$)", fontsize=14)
        plt.xlabel(r"$\lambda_1$", fontsize=12)
        plt.ylabel(r"$m_\phi$ [GeV]", fontsize=12)
        plt.legend(title="Método")
        plt.tight_layout()
        plt.savefig(os.path.join(IMG_DIR, "B_espacio_fase_comparativo.png"), dpi=300)
        plt.close()

        # GRAFICO C: Calidad de la Física
        plt.figure(figsize=(10, 6))
        sns.kdeplot(data=df_valid[df_valid['br_gaga'] > 1e-10], x='br_gaga', hue='Metodo', 
                    fill=True, common_norm=False, log_scale=True, palette=color_map)
        plt.title(r"Densidad de Puntos por Branching Ratio ($h \to \gamma\gamma$)", fontsize=14)
        plt.xlabel(r"BR($h \to \gamma\gamma$)", fontsize=12)
        plt.ylabel("Densidad de Modelos", fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(IMG_DIR, "C_branching_ratio_comparativo.png"), dpi=300)
        plt.close()

    print(f"[🚀] Gráficos guardados en: {IMG_DIR}")

if __name__ == "__main__":
    main()