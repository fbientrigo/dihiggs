#!/usr/bin/env bash
# Build and run the converter for Christopher's 7 benchmark points.
# Mirrors the 2hdmc/Makefile program rule:
#   g++ <src> -Isrc -std=c++11 -Wall -O3 -fopenmp -Llib -l2HDMC -lgsl -lgslcblas -lm
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THDMC="$REPO/2hdmc"
OUTDIR="$REPO/02_results/chris_hybrid"

mkdir -p "$OUTDIR"

echo "== git state =="
git -C "$REPO" rev-parse HEAD
git -C "$REPO" status --short --branch | head -5

echo "== ensure lib2HDMC.a is up to date (does not touch existing binaries) =="
make -C "$THDMC" lib

echo "== compile converter =="
COMPILE_CMD="g++ $REPO/chris/chris_points_to_hybrid.cpp -I$THDMC/src -std=c++11 -Wall -O3 -fopenmp -L$THDMC/lib -l2HDMC -lgsl -lgslcblas -lm -o $REPO/chris/chris_points_to_hybrid"
echo "$COMPILE_CMD"
$COMPILE_CMD

echo "== copy input for provenance =="
cp "$REPO/chris/chris_points.csv" "$OUTDIR/chris_points_input.csv"

echo "== run converter =="
RUN_CMD="$REPO/chris/chris_points_to_hybrid $REPO/chris/chris_points.csv $OUTDIR/chris_points_hybrid.csv"
echo "$RUN_CMD"
$RUN_CMD

echo "== done =="
