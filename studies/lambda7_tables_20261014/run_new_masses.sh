#!/usr/bin/env bash
# Direct 13 TeV LO integrations for the newly covered production kernels.
# B03/B10 at 250 GeV share m_phi and g exactly; only their external final-state
# decay widths differ, which do not enter this stable-phi 2->2 production ME.
set -euo pipefail
s="$(dirname "$0")/run_mg5.sh"
"$s" A_m250_g0_s101 508.130081300813 8.163e-09 101 250
"$s" A_m250_g0_s107 508.130081300813 8.163e-09 107 250
"$s" A_m300_g0_s101 731.707317073171 1.9733e-15 101 300
"$s" A_m300_g0_s107 731.707317073171 1.9733e-15 107 300
"$s" B_m300_l7zero_s101 63.648442645227125 1.9663956600676296e-15 101 300
"$s" B_m300_l7zero_s107 63.648442645227125 1.9663956600676296e-15 107 300
