#!/usr/bin/env bash
set -euo pipefail

# Builds the PhysLam1Scan binary for the lambda1 golden characterization suite.
#
# WHY THIS DOES NOT USE `make -C dihiggs app/PhysLam1Scan`
# -------------------------------------------------------
# The production Makefile path cannot build PhysLam1Scan from a clean checkout,
# for two independent pre-existing reasons (both documented in
# docs/characterization_lambda1.md as UNRESOLVED; neither is repaired here):
#
#   1. libdihiggs.a compiles COMMON_SRCS, which includes M2PointEvaluator.cpp,
#      which #includes M2PointEvaluator.hpp -- and dihiggs/include/*.hpp is NOT
#      tracked in git. `make -C dihiggs app/PhysLam1Scan` therefore fails with
#      "fatal error: M2PointEvaluator.hpp: No such file or directory".
#   2. LDFLAGS unconditionally links -lHiggsTools, which PhysLam1Scan.cpp never
#      uses (it includes only THDM.h, Constraints.h, DecayTable.h, ParamUtils.hpp).
#
# So we link PhysLam1Scan standalone against 2HDMC only. This is the same
# pattern the repository itself already uses for GenScanWithFixings, whose
# Makefile rule carries an equivalent comment ("standalone, bypasses
# libdihiggs.a ... deliberately does NOT link -lHiggsTools").
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
