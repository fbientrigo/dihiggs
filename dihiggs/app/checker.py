#!/usr/bin/env python3
import pandas as pd
import sys
import numpy as np
import argparse

def main():
    # 1. Parsing de argumentos robusto (Mejor práctica que sys.argv directo)
    parser = argparse.ArgumentParser(description="Analizador de puntos físicos válidos en 2HDM.")
    parser.add_argument("file", help="Ruta al archivo CSV a analizar")
    parser.add_argument("--precision", type=int, default=4, help="Decimales para agrupar m_phi (default: 4)")
    args = parser.parse_args()

    # Columnas estrictamente necesarias (Optimización de Memoria HPC)
    # Solo cargamos lo que usamos. Si el CSV tiene 50 columnas, esto reduce el uso de RAM un 90%.
    required_cols = ["m_phi", "positivity_ok", "unitarity_ok", "perturbativity_ok"]
    
    try:
        # Engine 'c' es más rápido. usecols filtra en lectura, no en memoria.
        df = pd.read_csv(args.file, usecols=lambda c: c in required_cols, engine='c')
    except ValueError as e:
        print(f"Error: Al menos una de las columnas requeridas {required_cols} no está en el CSV.")
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: Archivo {args.file} no encontrado.")
        sys.exit(1)

    total = len(df)
    if total == 0:
        print("El CSV está vacío.")
        sys.exit(0)

    # Definir condiciones
    conds = ["positivity_ok", "unitarity_ok", "perturbativity_ok"]
    
    # Pre-calcular máscaras booleanas (vectorización)
    # Convertimos a bool explícito para acelerar operaciones lógicas bitwise
    masks = {c: df[c].astype(bool) for c in conds}
    
    # --- REPORTE GENERAL ---
    print(f"{'='*40}")
    print(f"ANÁLISIS DE: {args.file}")
    print(f"Total puntos escaneados: {total}")
    print(f"{'='*40}\n")

    def pct(count):
        return 100.0 * count / total

    # 1. Estadísticas individuales
    print("--- Filtros Individuales ---")
    for c in conds:
        n = masks[c].sum()
        print(f"{c:20s}: {pct(n):6.2f}% (n={n})")

    # 2. Intersección Triple (La condición física final)
    valid_mask = masks["positivity_ok"] & masks["unitarity_ok"] & masks["perturbativity_ok"]
    n_valid = valid_mask.sum()
    
    print(f"\n--- INTERSECCIÓN TRIPLE (Puntos Físicos) ---")
    print(f"Valid (ALL OK)      : {pct(n_valid):6.2f}%")
    print(f"Total N             : {n_valid}")

    # --- ANÁLISIS DE MASAS (Tu requerimiento específico) ---
    if n_valid > 0:
        print(f"\n{'-'*40}")
        print(f"DISTRIBUCIÓN DE m_phi EN PUNTOS VÁLIDOS")
        print(f"{'-'*40}")
        
        # Filtramos solo las filas válidas
        valid_df = df.loc[valid_mask].copy()
        
        # Redondeamos para agrupar flotantes cercanos (HPC Best Practice para grids)
        valid_df["m_phi_rounded"] = valid_df["m_phi"].round(args.precision)
        
        # Contamos ocurrencias por masa
        mass_counts = valid_df["m_phi_rounded"].value_counts().sort_index()
        
        print(f"{'m_phi (GeV)':<15} | {'N Puntos':<10} | {'% del Total Válido'}")
        print("-" * 45)
        
        for mass, count in mass_counts.items():
            percent_of_valid = (count / n_valid) * 100.0
            print(f"{mass:<15.4f} | {count:<10} | {percent_of_valid:6.2f}%")
            
        # Set python puro para copiar y pegar si es necesario
        print(f"\nSet de m_phi válidos (Python syntax):")
        print(set(mass_counts.index.tolist()))
    else:
        print("\n[!] No se encontraron puntos físicos válidos. No se puede analizar m_phi.")

if __name__ == "__main__":
    main()