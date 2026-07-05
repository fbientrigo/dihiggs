#!/usr/bin/env bash
# Reproducibly plot Christopher's production run using existing mlpython/2603 tools.
#
# Usage:
#   bash mlpython/2603/run_chris_existing_plots.sh [INPUT_CSV [OUTPUT_ROOT]]
#
# Defaults:
#   INPUT_CSV   = data/lam1eq1_allchris_var10000_2026jun/chunked/silver_all.csv
#   OUTPUT_ROOT = data/lam1eq1_allchris_var10000_2026jun/chunked/plots_chris
#
# Smoke command:
#   bash mlpython/2603/run_chris_existing_plots.sh \
#     data/lam1eq1_allchris_var10000_2026jun/chunked/silver_all.csv \
#     data/lam1eq1_allchris_var10000_2026jun/chunked/plots_chris

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

INPUT_CSV="${1:-$REPO_ROOT/data/lam1eq1_allchris_var10000_2026jun/chunked/silver_all.csv}"
OUTPUT_ROOT="${2:-$REPO_ROOT/data/lam1eq1_allchris_var10000_2026jun/chunked/plots_chris}"

WRAPPER="$REPO_ROOT/scripts/run_chris_existing_plots.py"

if python3 -c "import polars" 2>/dev/null; then
    exec python3 "$WRAPPER" --input "$INPUT_CSV" --output-root "$OUTPUT_ROOT"
else
    # polars not on PATH — use uv to create a temporary env with the required packages
    exec uv run \
        --with "polars,pyarrow,pandas,matplotlib,numpy" \
        "$WRAPPER" --input "$INPUT_CSV" --output-root "$OUTPUT_ROOT"
fi
