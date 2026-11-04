#!/usr/bin/env bash
set -euo pipefail

# Builds the PhysLam1Scan binary for the lambda1 golden characterization suite.
#
# PhysLam1Scan is the frozen legacy oracle characterized by
# tests/test_golden_lambda1.py; it is not part of `make -C dihiggs all`.
# It links standalone against 2HDMC only (no libdihiggs.a, no HiggsTools).
#
# COMPILER FLAGS are copied verbatim from dihiggs/Makefile's CXXFLAGS so that
# the characterized binary matches the production build settings:
#   -fopenmp -std=gnu++17 -Wall -Wextra -Wpedantic -O0 -g
#
# No production source file is modified by this script.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

OUT_DIR="${LAMBDA1_CHAR_OUT_DIR:-${REPO_ROOT}/build/lambda1_char}"
mkdir -p "${OUT_DIR}"

echo "[L1] Building vendored 2HDMC (lib/lib2HDMC.a)..."
make -C 2hdmc lib/lib2HDMC.a >/dev/null

echo "[L1] Building PhysLam1Scan (standalone, production CXXFLAGS)..."
g++ -I2hdmc/src -Idihiggs/src \
    -fopenmp -std=gnu++17 -Wall -Wextra -Wpedantic -O0 -g \
    dihiggs/src/PhysLam1Scan.cpp dihiggs/src/ParamUtils.cpp \
    -o "${OUT_DIR}/PhysLam1Scan" \
    -L2hdmc/lib -l2HDMC -lgsl -lgslcblas -lm

echo "[L1] OK: ${OUT_DIR}/PhysLam1Scan"
