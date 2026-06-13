#!/usr/bin/env bash
# Build the Stage 1 OpenMP fast scanner (chris/CalcLambda1ScanFixings).
# Self-contained: no 2HDMC/GSL dependencies.
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

echo "== compile CalcLambda1Core.c =="
gcc -c CalcLambda1Core.c -std=c99 -O3 -Wall -o CalcLambda1Core.o

echo "== compile + link CalcLambda1ScanFixings =="
g++ CalcLambda1ScanFixings.cpp CalcLambda1Core.o -std=c++17 -O3 -Wall -Wextra -fopenmp -lm \
    -o CalcLambda1ScanFixings

echo "== done =="
