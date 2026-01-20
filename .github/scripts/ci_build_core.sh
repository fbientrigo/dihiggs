#!/bin/bash
set -u
BLUE='\033[0;34m'
NC='\033[0m'
PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." && pwd )"
cd "$PROJECT_ROOT"

echo -e "${BLUE}[BUILD] Compilando Core Library...${NC}"

# 1. Asegurar dependencias (Makefiles robustos ya se encargan, pero forzamos por seguridad)
if [ ! -d "dihiggs" ]; then echo "❌ No existe dihiggs dir"; exit 1; fi

cd dihiggs
make clean
# Compilamos la librería estática
make lib/libdihiggs.a
# Opcional: Compilar ya los tests si queremos, o dejarlos para la fase de test.
# Para HPC, mejor compilar todo lo pesado aquí.
cd ../tests/unit
make dihiggs_unit_tests

echo -e "${BLUE}✅ Build completado. Artefactos listos.${NC}"