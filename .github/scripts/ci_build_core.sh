#!/bin/bash
set -u
BLUE='\033[0;34m'
NC='\033[0m'
GREEN='\033[0;32m'
RED='\033[0;31m'

# Project root is one level up from .github/scripts/
PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." && pwd )"
cd "$PROJECT_ROOT"

echo -e "${BLUE}[BUILD] Starting core build pipeline...${NC}"

# Step 1: Ensure dependencies exist
if [ ! -d "dihiggs" ]; then
    echo -e "${RED}❌ No dihiggs directory${NC}"
    exit 1
fi
if [ ! -d "2hdmc" ]; then
    echo -e "${RED}❌ No 2hdmc directory${NC}"
    exit 1
fi

# Step 2: Build 2HDMC library
echo -e "${BLUE}[BUILD] Building 2HDMC library...${NC}"
cd 2hdmc
make clean 2>/dev/null || true
make lib/lib2HDMC.a
if [ ! -f "lib/lib2HDMC.a" ]; then
    echo -e "${RED}❌ 2HDMC library build failed${NC}"
    exit 1
fi

# Step 3: Verify set_param_phys_lam1 symbol exists (plan line 366)
echo -e "${BLUE}[BUILD] Verifying set_param_phys_lam1 symbol in lib2HDMC.a...${NC}"
if nm lib/lib2HDMC.a 2>/dev/null | grep -q "set_param_phys_lam1"; then
    echo -e "${GREEN}✅ set_param_phys_lam1 symbol verified${NC}"
else
    echo -e "${RED}❌ set_param_phys_lam1 symbol NOT found in lib2HDMC.a${NC}"
    echo "   This indicates the 2HDMC library was not rebuilt after adding the new API."
    echo "   Please rebuild 2HDMC before continuing."
    exit 1
fi

cd ..

# Step 4: Build dihiggs (includes PhysLam1Scan)
echo -e "${BLUE}[BUILD] Building dihiggs (includes PhysLam1Scan)...${NC}"
cd dihiggs
make clean 2>/dev/null || true
make
if [ ! -f "app/PhysLam1Scan" ]; then
    echo -e "${RED}❌ PhysLam1Scan binary not built${NC}"
    exit 1
fi
echo -e "${GREEN}✅ PhysLam1Scan built successfully${NC}"

cd ..

# Step 5: Build unit tests
echo -e "${BLUE}[BUILD] Building unit tests...${NC}"
cd tests/unit
make clean 2>/dev/null || true
make dihiggs_unit_tests
if [ ! -f "dihiggs_unit_tests" ]; then
    echo -e "${RED}❌ Unit tests binary not built${NC}"
    exit 1
fi
echo -e "${GREEN}✅ Unit tests built successfully${NC}"

cd ../..

echo -e "${BLUE}[BUILD] Verifying all artifacts exist...${NC}"

# Verify artifacts
ARTIFACTS_OK=true
if [ -f "dihiggs/app/PhysLam1Scan" ]; then
    echo -e "${GREEN}✅ dihiggs/app/PhysLam1Scan${NC}"
else
    echo -e "${RED}❌ dihiggs/app/PhysLam1Scan MISSING${NC}"
    ARTIFACTS_OK=false
fi

if [ -f "tests/unit/dihiggs_unit_tests" ]; then
    echo -e "${GREEN}✅ tests/unit/dihiggs_unit_tests${NC}"
else
    echo -e "${RED}❌ tests/unit/dihiggs_unit_tests MISSING${NC}"
    ARTIFACTS_OK=false
fi

if [ "$ARTIFACTS_OK" = true ]; then
    echo -e "${GREEN}✅ Build completed. All artifacts ready.${NC}"
    exit 0
else
    echo -e "${RED}❌ Build completed but some artifacts are missing${NC}"
    exit 1
fi
