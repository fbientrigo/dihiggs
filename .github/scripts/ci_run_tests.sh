#!/bin/bash
set -u
TYPE=${1:-unit} # Por defecto 'unit'
PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." && pwd )"
cd "$PROJECT_ROOT"

echo "Running Test Suite: $TYPE"

if [ "$TYPE" == "unit" ]; then
    cd tests/unit
    chmod +x dihiggs_unit_tests
    ./dihiggs_unit_tests -s
elif [ "$TYPE" == "oracle" ]; then
    echo "🚧 Oracle tests not implemented yet (Placeholder)"
    # Aquí iría: ./dihiggs_unit_tests "[oracle]"
else
    echo "Unknown test type"
    exit 1
fi