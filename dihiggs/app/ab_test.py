import subprocess
import time
import os
import pandas as pd
from pathlib import Path

# --- CONFIGURACIÓN DEL TEST ---
# Rutas a los ejecutables compilados
EXE_A = "../dihiggs/app/PhysParamScan"       # Estrategia A: Linear Search (Old)
EXE_B = "../dihiggs/app/PhysScanWithFixings" # Estrategia B: Analytical Inversion (New)

# Parámetros físicos comunes para la comparación
MPHI_MIN = 150.0
MPHI_MAX = 300.0
N_MPHI   = 20
MA       = 300.0
SIN_BA   = 0.999
TAN_BETA = 1000.0  # Probar con valores altos (1e3, 1e4)
LAM6     = 0.1
LAM7     = 0.0

# Parámetros específicos de cada estrategia
# A: Busca m12 alrededor de una conjetura
N_M12_A  = 10000 

# B: Escanea Lambda1
LAM1_MIN = 0.0
LAM1_MAX = 4.0
N_LAM1_B = 200

OUTPUT_DIR = "ab_test_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def run_strategy_A():
    """Ejecuta PhysParamScan (Linear Search)"""
    output_file = f"{OUTPUT_DIR}/strategy_A_tb{int(TAN_BETA)}.csv"
    
    # Firma: output N_mphi tan_beta tol_exp N_m12 y_type --flags...
    cmd = [
        EXE_A,
        output_file,
        str(N_MPHI),
        str(TAN_BETA),
        "3",            # tol_exp
        str(N_M12_A),   # N_m12
        "1",            # y_type
        f"--mphi_min={MPHI_MIN}",
        f"--mphi_max={MPHI_MAX}",
        f"--mA={MA}",
        f"--sin_ba={SIN_BA}",
        f"--lambda6={LAM6}",
        f"--lambda7={LAM7}"
    ]
    
    print(f"\n[TEST A] Ejecutando Linear Search (N_m12={N_M12_A})...")
    start = time.time()
    subprocess.run(cmd, check=True)
    elapsed = time.time() - start
    
    # Análisis rápido de resultados
    try:
        df = pd.read_csv(output_file)
        valid_points = len(df)
    except:
        valid_points = 0
        
    return elapsed, valid_points, output_file

def run_strategy_B():
    """Ejecuta PhysScanWithFixings (Analytical Inversion)"""
    output_file = f"{OUTPUT_DIR}/strategy_B_tb{int(TAN_BETA)}.csv"
    
    # Firma: mphi_min mphi_max N_mphi lam1_min lam1_max N_lam1 mA sin_ba tan_beta l6 l7 output
    cmd = [
        EXE_B,
        str(MPHI_MIN),
        str(MPHI_MAX),
        str(N_MPHI),
        str(LAM1_MIN),
        str(LAM1_MAX),
        str(N_LAM1_B),
        str(MA),
        str(SIN_BA),
        str(TAN_BETA),
        str(LAM6),
        str(LAM7),
        output_file
    ]
    
    print(f"\n[TEST B] Ejecutando Analytical Inversion (N_lam1={N_LAM1_B})...")
    start = time.time()
    subprocess.run(cmd, check=True)
    elapsed = time.time() - start
    
    try:
        df = pd.read_csv(output_file)
        valid_points = len(df)
    except:
        valid_points = 0
        
    return elapsed, valid_points, output_file

def main():
    print(f"=== A/B TESTING: Parameter Scan Strategies ===")
    print(f"Tan Beta: {TAN_BETA}")
    print(f"m_phi range: [{MPHI_MIN}, {MPHI_MAX}] ({N_MPHI} steps)")
    
    # Run A
    try:
        time_a, valid_a, file_a = run_strategy_A()
    except Exception as e:
        print(f"Estrategia A falló: {e}")
        time_a, valid_a = 1e-9, 0

    # Run B
    try:
        time_b, valid_b, file_b = run_strategy_B()
    except Exception as e:
        print(f"Estrategia B falló: {e}")
        time_b, valid_b = 1e-9, 0

    # Reporte
    print("\n" + "="*60)
    print(f"{'METRIC':<20} | {'STRATEGY A (Linear)':<20} | {'STRATEGY B (Inversion)':<20}")
    print("-" * 60)
    print(f"{'Time (s)':<20} | {time_a:<20.4f} | {time_b:<20.4f}")
    print(f"{'Total Valid Points':<20} | {valid_a:<20} | {valid_b:<20}")
    
    speed_a = valid_a / time_a if time_a > 0 else 0
    speed_b = valid_b / time_b if time_b > 0 else 0
    
    print(f"{'Speed (pts/s)':<20} | {speed_a:<20.1f} | {speed_b:<20.1f}")
    
    if speed_a > 0:
        improvement = speed_b / speed_a
        print(f"\n>> Factor de Mejora (Speedup): {improvement:.2f}x")
    else:
        print("\n>> Estrategia A no produjo puntos válidos (Speedup infinito).")
        
    print("="*60)

if __name__ == "__main__":
    main()