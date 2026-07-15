#!/usr/bin/env bash
set -euo pipefail

# Smallest practical entry point for the lambda1 golden characterization suite:
# build the vendored 2HDMC fork + PhysLam1Scan, then run the golden tests with
# the binary REQUIRED (absence fails, it does not skip).
#
# See docs/characterization_lambda1.md. Runtime is dominated by the 2HDMC build
# (~15 s); the golden suite itself runs 9 single-point/small-grid cases.
#
# No HiggsTools, no HB/HS datasets, and no campaign scan are involved.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

bash scripts/build_lambda1_characterization.sh

DIHIGGS_REQUIRE_LAMBDA1_BINARY=1 python3 -m pytest -v tests/test_golden_lambda1.py
