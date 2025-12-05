#!/bin/bash
# run_lifetime_validation.sh
# Objetivo: A/B Test para verificar sensibilidad del Lifetime (Width) a Lambda6
# Uso:
#   ./run_lifetime_validation.sh <exp_tan_beta> <N_LAM1>
# Ejemplos:
#   ./run_lifetime_validation.sh 2 1000   # tan_beta = 10^2,  N_LAM1 = 1000
#   ./run_lifetime_validation.sh 5 5000   # tan_beta = 10^5,  N_LAM1 = 5000

set -euo pipefail

# ============================
# 1) Parseo de argumentos
# ============================
if [ "$#" -ne 2 ]; then
    echo "Uso: $0 <exp_tan_beta> <N_LAM1>"
    echo "Ejemplo: $0 3 1000   # tan_beta = 10^3, N_LAM1 = 1000"
    exit 1
fi

EXP_TB="$1"
N_LAM1="$2"

# Comprobamos que EXP_TB sea entero no negativo
if ! [[ "$EXP_TB" =~ ^[0-9]+$ ]]; then
    echo "Error: <exp_tan_beta> debe ser un entero no negativo (0,1,2,3,...)"
    exit 1
fi

# Comprobamos que N_LAM1 sea entero positivo
if ! [[ "$N_LAM1" =~ ^[0-9]+$ ]] || [ "$N_LAM1" -le 0 ]; then
    echo "Error: <N_LAM1> debe ser un entero positivo (1,2,3,...)"
    exit 1
fi

# tan_beta = 10^EXP_TB
TAN_BETA=$((10**EXP_TB))

BINARY="./PhysScanWithFixings"
OUT_DIR="./out_l6_test/validation_tb1e${EXP_TB}"
mkdir -p "$OUT_DIR"

echo "[CFG] exp_tan_beta = ${EXP_TB}"
echo "[CFG] tan_beta     = 10^${EXP_TB} = ${TAN_BETA}"
echo "[CFG] N_LAM1       = ${N_LAM1}"
echo "[CFG] Output dir   = ${OUT_DIR}"

# ============================
# 2) Parámetros fijos
# ============================
MPHI_MIN=150
MPHI_MAX=300
N_MPHI=15

LAM1_MIN=0.0
LAM1_MAX=1.0

MA=300       # Heavy doublet degenerado (mA ~ mphi para simplificar)
SIN_BA=0.999 # Casi alineamiento
LAM7=0.0     # Apagamos lambda7 para aislar lambda6

# ============================
# 3) Caso A: Simetría Z2 Exacta (Lambda6 = 0.0)
# ============================
echo "[RUN] Caso A: Lambda6 = 0.0 (Baseline, tan_beta = ${TAN_BETA}, N_LAM1 = ${N_LAM1})"
"$BINARY" \
    "$MPHI_MIN" "$MPHI_MAX" "$N_MPHI" \
    "$LAM1_MIN" "$LAM1_MAX" "$N_LAM1" \
    "$MA" "$SIN_BA" "$TAN_BETA" \
    0.0 "$LAM7" \
    "${OUT_DIR}/case_A_l6_zero.csv"

# ============================
# 4) Caso B1: Ruptura Leve (Lambda6 = 0.01)
# ============================
echo "[RUN] Caso B1: Lambda6 = 0.01 (Symmetry Breaking)"
"$BINARY" \
    "$MPHI_MIN" "$MPHI_MAX" "$N_MPHI" \
    "$LAM1_MIN" "$LAM1_MAX" "$N_LAM1" \
    "$MA" "$SIN_BA" "$TAN_BETA" \
    0.01 "$LAM7" \
    "${OUT_DIR}/case_B_l6_point01.csv"

# ============================
# 5) Caso B2: Ruptura Fuerte (Lambda6 = 0.1)
# ============================
echo "[RUN] Caso B2: Lambda6 = 0.1 (Symmetry Breaking)"
"$BINARY" \
    "$MPHI_MIN" "$MPHI_MAX" "$N_MPHI" \
    "$LAM1_MIN" "$LAM1_MAX" "$N_LAM1" \
    "$MA" "$SIN_BA" "$TAN_BETA" \
    0.1 "$LAM7" \
    "${OUT_DIR}/case_B_l6_point1.csv"

echo "[DONE] Resultados generados en ${OUT_DIR}."
echo "       Archivos: case_A_l6_zero.csv, case_B_l6_point01.csv, case_B_l6_point1.csv"
