#!/bin/bash
set -u
BLUE='\033[0;34m'
NC='\033[0m'
GREEN='\033[0;32m'
RED='\033[0;31m'

TYPE=${1:-unit}  # Default to 'unit', can be 'oracle'
PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." && pwd )"
cd "$PROJECT_ROOT"

echo -e "${BLUE}[TEST] Running Test Suite: $TYPE${NC}"

# Ensure binaries are executable
if [ "$TYPE" == "unit" ] || [ "$TYPE" == "oracle" ]; then
    if [ -f "tests/unit/dihiggs_unit_tests" ]; then
        chmod +x tests/unit/dihiggs_unit_tests
        echo -e "${BLUE}[TEST] dihiggs_unit_tests permissions set${NC}"
    else
        echo -e "${RED}❌ tests/unit/dihiggs_unit_tests not found${NC}"
        echo "   Did you run ci_build_core.sh first?"
        exit 1
    fi
fi

# PhysLam1Scan should also be executable (for manual testing)
if [ -f "dihiggs/app/PhysLam1Scan" ]; then
    chmod +x dihiggs/app/PhysLam1Scan
    echo -e "${BLUE}[TEST] dihiggs/app/PhysLam1Scan permissions set${NC}"
fi

if [ "$TYPE" == "unit" ]; then
    echo -e "${BLUE}[TEST] Running unit tests (all tags)...${NC}"
    ./tests/unit/dihiggs_unit_tests -s
    RC=$?

elif [ "$TYPE" == "oracle" ]; then
    echo -e "${BLUE}[TEST] Running oracle tests (tagged [oracle])...${NC}"
    echo -e "${BLUE}[TEST] Setting OMP_NUM_THREADS=1 for thread safety${NC}"

    # Execute oracle tests with single-threaded enforcement
    # This is critical because 2HDMC/GSL have global handlers that are not thread-safe
    export OMP_NUM_THREADS=1

    ./tests/unit/dihiggs_unit_tests "[oracle]" -s
    RC=$?

    if [ $RC -eq 0 ]; then
        echo -e "${GREEN}✅ Oracle tests passed${NC}"
    else
        echo -e "${RED}❌ Oracle tests failed (exit code: $RC)${NC}"
    fi

else
    echo -e "${RED}❌ Unknown test type: $TYPE${NC}"
    echo "   Usage: $0 [unit|oracle]"
    exit 1
fi

exit $RC
