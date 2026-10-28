#!/usr/bin/env bash
# Direct LO integration of g g > H > h2 h2 (13 TeV, nn23lo1, dynamical scale -1) for one
# value of the physical coupling.  The UFO vertex is GC_90 = i*GHphiphi, so GHphiphi = -g
# gives the Feynman rule -i g of the canonical appendix.
#   usage: run_mg5.sh RUN_ID g_hphiphi_GeV width_phi_GeV seed [m_phi_GeV]
set -euo pipefail
id=$1; g=$2; wphi=$3; seed=$4; mphi=${5:-200}
here="${DIHIGGS_ROOT:?set DIHIGGS_ROOT (. dihiggs_hep_cross/scripts/env.sh)}/reports/lambda7_tables_20261014"  # output package (lake), not code
tmpl=/tmp/dihiggs-honly-replay-cards-20261007T0130/mphi200_physical_seed101   # prepared, never run
py311=${PY311:?set PY311 to a python3.11 with the MG5 deps}
dir="$here/mg5/runs/$id"
[[ -e "$dir" ]] && { echo "refusing to overwrite $dir" >&2; exit 2; }
mkdir -p "$here/mg5/runs" "$here/mg5/logs"
cp -a "$tmpl" "$dir"
pc="$dir/Cards/param_card.dat"; rc="$dir/Cards/run_card.dat"
sed -i -E "s/^( +3 ).*# GHphiphi.*/\1-$g # GHphiphi = -g_hphiphi (FR = -i g)/" "$pc"
sed -i -E "s/^( +25 ).*# MH.*/\1125.13 # MH/" "$pc"
sed -i -E "s/^( +9000006 ).*# Mh2.*/\1$mphi # Mh2/" "$pc"
sed -i -E "s/^DECAY +9000006 .*/DECAY   9000006 $wphi # 2HDMC total width of phi/" "$pc"
sed -i -E "s/^ *[0-9]+ += iseed.*/  $seed = iseed   ! rnd seed/" "$rc"
sed -i -E "s/^ *[A-Za-z]+ += use_syst.*/  False = use_syst/" "$rc"
sed -i -E "s/^ *[^ ]+[[:space:]]*= run_tag.*/  $id = run_tag/" "$rc"
printf 'run_mode = 2\nnb_core = 2\n' >> "$dir/Cards/me5_configuration.txt"
sed -i -E "s|^mg5_path = .*|mg5_path = $DIHIGGS_ROOT/experiments/trilinear_normalization_20261010T191851Z/engine/mg5|" "$dir/Cards/me5_configuration.txt"
grep -E "GHphiphi| MH$|DECAY +9000006|DECAY +25 " "$pc"; grep -E "iseed|nevents|pdlabel|lhaid|dynamical_scale|ebeam" "$rc"
export PYTHONPATH=/tmp/dihiggs-honly-replay-20261007T012700Z/python_deps
(cd "$dir" && "$py311" ./bin/generate_events -f "$id") > "$here/mg5/logs/$id.log" 2>&1
grep -E "Cross-section|Nb of events" "$here/mg5/logs/$id.log" | tail -2
