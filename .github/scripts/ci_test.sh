#!/bin/bash
# No usamos 'set -e' al inicio para poder capturar fallos específicos y reportar logs antes de salir
set -u # Error si usamos variables no definidas

# Configuración Visual
BLUE='\033[0;34m'
CYAN='\033[0;36m'
RED='\033[0;31m'
NC='\033[0m'

# Determinar el root del proyecto
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/../.." && pwd )"
cd "$PROJECT_ROOT"

echo -e "${BLUE}[CI-DEBUG] Root directory: $(pwd)${NC}"

# --- 1. Verificación de Entorno y Estructura ---
echo -e "${CYAN}[1/4] Verificando presencia de archivos fuente...${NC}"

# Inspección profunda: Esto nos dirá si los archivos sobrevivieron al Git Checkout
echo "Estructura detectada en dihiggs/src/:"
if [ -d "dihiggs/src" ]; then
    ls -la dihiggs/src/
else
    echo -e "${RED}❌ ERROR: El directorio dihiggs/src NO EXISTE en $(pwd)${NC}"
    exit 1
fi

if [ -z "${HB_DATASET:-}" ]; then
    echo -e "${RED}❌ HB_DATASET not set.${NC}"
    exit 1
fi

# --- 2. Compilación Core con Verbose Logs ---
echo -e "${CYAN}[2/4] Compilando Core Library...${NC}"
cd dihiggs

# Forzamos una limpieza y luego compilamos
make clean

# Verificamos si existe la carpeta obj/ antes de empezar
mkdir -p obj lib

echo "Ejecutando make para libdihiggs.a..."
# Si falla, imprimimos un diagnóstico antes de que el script termine
if ! make; then
    echo -e "${RED}❌ El proceso de compilación falló.${NC}"
    echo "Diagnóstico de Make:"
    echo " - Makefile actual: $(ls Makefile)"
    echo " - ¿Existe src/ParamUtils.cpp? $(ls -l src/ParamUtils.cpp 2>/dev/null || echo 'NO')"
    echo " - ¿Existe carpeta obj? $(ls -d obj 2>/dev/null || echo 'NO')"
    exit 1
fi
cd "$PROJECT_ROOT"

# --- 3. Compilación y Ejecución de Tests ---
echo -e "${CYAN}[3/4] Compilando Unit Tests...${NC}"
cd tests/unit

if ! make; then
    echo -e "${RED}❌ Falló la compilación de los tests.${NC}"
    exit 1
fi

echo -e "${CYAN}[4/4] Ejecutando Tests...${NC}"
if ! ./dihiggs_unit_tests -s; then
    echo -e "${RED}❌ Algunos tests fallaron o el binario no se ejecutó.${NC}"
    exit 1
fi

echo -e "${BLUE}✅ CI Pipeline finalizado con éxito.${NC}"