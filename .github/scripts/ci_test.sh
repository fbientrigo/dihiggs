#!/bin/bash
set -e

# Configuración Visual
BLUE='\033[0;34m'
NC='\033[0m'

# Determinar el root del proyecto (2 niveles arriba de donde está el script)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/../.." && pwd )"

cd "$PROJECT_ROOT"

echo -e "${BLUE}[CI] Working directory: $(pwd)${NC}"

# --- 1. Verificación de Entorno ---
if [ -z "$HB_DATASET" ]; then
    echo "❌ HB_DATASET not set. Is the DevContainer environment correct?"
    exit 1
fi

# --- 2. Compilación Core ---
echo -e "${BLUE}[CI] Compiling Core Library...${NC}"
cd dihiggs
make clean && make
cd "$PROJECT_ROOT"

# --- 3. Ejecución de Tests ---
echo -e "${BLUE}[CI] Running Unit Tests...${NC}"
cd tests/unit
make clean && make
./dihiggs_unit_tests -s