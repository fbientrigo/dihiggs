#!/usr/bin/env python3
import argparse
import os
import subprocess
import time
import sys
from pathlib import Path
from datetime import datetime

# =============================================================================
# CONFIGURACIÓN DEL ESPACIO DE PARÁMETROS
# =============================================================================
# Rangos y constantes físicas
MPHI_MIN = 130.0
MPHI_MAX = 290.0
N_MPHI   = 15      # Puntos en m_phi

LAM1_MIN = 0.0
LAM1_MAX = 12.0
# N_LAM1 se define por argumento (default 666)

MA_FIXED = 300.0   # mA = mHp
SIN_BA   = 1.0     # Alineamiento exacto
LAMBDA_6 = 0.1
LAMBDA_7 = 0.0

# Lista de TanBeta a explorar
#TAN_BETA_LIST = [10, 50, 100, 200, 500, 1000, 1500, 2000, 5000] # previos ya heco en night
#TAN_BETA_LIST = [8000, 9000, 10000, 15000, 20000, 25000, 30000] # agregados
#TAN_BETA_LIST = [50000, 70000, 100000] # agregados
#TAN_BETA_LIST = [0.1, 1, 30,200000,500000, 1000000]
TAN_BETA_LIST = [35000,80000,150000]

# =============================================================================
# UTILIDADES
# =============================================================================

def log_message(msg, log_file):
    """Escribe un mensaje en consola y en el archivo de log con timestamp."""
    timestamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    full_msg = f"{timestamp} {msg}"
    print(full_msg)
    with open(log_file, "a") as f:
        f.write(full_msg + "\n")

def parse_args():
    parser = argparse.ArgumentParser(description="Orquestador de escaneos 2HDM (C++ Driver)")
    parser.add_argument("--N-lam1", type=int, default=666, 
                        help="Número de puntos para lambda_1 (Default: 666. Usa 10 para test).")
    parser.add_argument("--exec", type=str, default="./PhysScanWithFixings",
                        help="Ruta al ejecutable C++ compilado.")
    parser.add_argument("--outdir", type=str, default="scan_results_nightly",
                        help="Directorio raíz para los resultados.")
    parser.add_argument("--threads", type=str, default=None,
                        help="Número de hilos OMP (opcional, por defecto usa todos).")
    parser.add_argument("--force", action="store_true",
                        help="Sobrescribir archivos existentes si se encuentran.")
    return parser.parse_args()

# =============================================================================
# MAIN ORCHESTRATOR
# =============================================================================

def main():
    args = parse_args()
    
    # 1. Configuración de Entorno
    root_dir = Path(args.outdir)
    root_dir.mkdir(parents=True, exist_ok=True)
    log_path = root_dir / "orchestrator.log"
    
    executable = Path(args.exec).resolve()
    if not executable.exists():
        log_message(f"[ERROR] No se encuentra el ejecutable: {executable}", log_path)
        sys.exit(1)

    # Configurar OpenMP si se solicita
    env = os.environ.copy()
    if args.threads:
        env["OMP_NUM_THREADS"] = args.threads
        log_message(f"[INIT] Forzando OMP_NUM_THREADS = {args.threads}", log_path)

    log_message(f"[INIT] Iniciando orquestación de escaneo.", log_path)
    log_message(f"[CONF] m_phi: [{MPHI_MIN}, {MPHI_MAX}] (N={N_MPHI})", log_path)
    log_message(f"[CONF] lam_1: [{LAM1_MIN}, {LAM1_MAX}] (N={args.N_lam1})", log_path)
    log_message(f"[CONF] Constantes: mA={MA_FIXED}, sin(b-a)={SIN_BA}, l6={LAMBDA_6}, l7={LAMBDA_7}", log_path)
    log_message(f"[CONF] TanBeta List: {TAN_BETA_LIST}", log_path)

    total_tasks = len(TAN_BETA_LIST)
    
    # Variables estadísticas globales
    grand_total_start = time.time()
    total_points_processed = 0
    tasks_run_count = 0
    
    # 2. Bucle Principal (Resiliente)
    for idx, tb in enumerate(TAN_BETA_LIST):
        
        # Crear subdirectorio ordenado (ej. tb_01000)
        folder_name = f"tb_{int(tb):05d}" 
        current_dir = root_dir / folder_name
        current_dir.mkdir(exist_ok=True)
        
        output_csv = current_dir / f"scan_tb_{tb}.csv"
        
        # Chequeo de Resiliencia: ¿Ya existe?
        if output_csv.exists() and not args.force:
            log_message(f"[SKIP] ({idx+1}/{total_tasks}) TanBeta={tb} ya existe en {output_csv}", log_path)
            continue

        log_message(f"[RUN ] ({idx+1}/{total_tasks}) Procesando TanBeta={tb} ...", log_path)
        
        # 3. Construcción del Comando (Mapping EXACTO al main C++)
        # argv[1]=mphi_min, argv[2]=mphi_max, argv[3]=N_mphi
        # argv[4]=lam1_min, argv[5]=lam1_max, argv[6]=N_lam1
        # argv[7]=mA, argv[8]=sin_ba, argv[9]=tan_beta
        # argv[10]=l6, argv[11]=l7, argv[12]=output
        
        cmd = [
            str(executable),
            f"{MPHI_MIN:.4f}",
            f"{MPHI_MAX:.4f}",
            str(N_MPHI),
            f"{LAM1_MIN:.4f}",
            f"{LAM1_MAX:.4f}",
            str(args.N_lam1),
            f"{MA_FIXED:.4f}",
            f"{SIN_BA:.4f}",
            f"{tb:.4f}",        # El tan_beta actual del loop
            f"{LAMBDA_6:.4f}",
            f"{LAMBDA_7:.4f}",
            str(output_csv)
        ]

        # 4. Ejecución y Medición
        t0 = time.time()
        try:
            result = subprocess.run(cmd, env=env, capture_output=True, text=True)
            t1 = time.time()
            elapsed = t1 - t0
            
            if result.returncode == 0:
                tasks_run_count += 1
                
                # Parsing de salida C++ para estadísticas
                ok_points = "Unknown"
                attempts  = 0
                
                for line in result.stdout.splitlines():
                    if "TRIPLE_OK_POINTS" in line:
                        try:
                            ok_points = line.split()[-1]
                        except: pass
                    if "Total Attempts:" in line:
                        try:
                            attempts = int(line.split(":")[-1].strip())
                        except: pass
                
                total_points_processed += attempts
                
                log_message(f"[DONE] TanBeta={tb} completado en {elapsed:.2f}s. "
                            f"Intentos: {attempts}, TripleOK: {ok_points}", log_path)
            else:
                log_message(f"[FAIL] TanBeta={tb} falló con código {result.returncode}", log_path)
                log_message(f"[STDERR] {result.stderr}", log_path)

        except Exception as e:
            log_message(f"[CRASH] Excepción ejecutando TanBeta={tb}: {e}", log_path)

    # =========================================================================
    # ESTADÍSTICAS FINALES
    # =========================================================================
    grand_total_duration = time.time() - grand_total_start
    
    log_message("-" * 60, log_path)
    log_message("[FIN ] Orquestación finalizada.", log_path)
    log_message(f"       Tiempo Total Ejecución: {grand_total_duration:.2f} s ({grand_total_duration/60:.2f} min)", log_path)
    log_message(f"       Total Puntos Calculados: {total_points_processed}", log_path)
    
    if total_points_processed > 0 and grand_total_duration > 0:
        sec_per_point = grand_total_duration / total_points_processed
        points_per_hour = total_points_processed / (grand_total_duration / 3600.0)
        
        log_message(f"       Velocidad Promedio:      {sec_per_point:.5f} segundos/punto", log_path)
        log_message(f"       Throughput:              {int(points_per_hour)} puntos/hora", log_path)
    else:
        log_message("       No se procesaron puntos nuevos (¿todo skippeado?).", log_path)
        
    log_message("-" * 60, log_path)

if __name__ == "__main__":
    main()
