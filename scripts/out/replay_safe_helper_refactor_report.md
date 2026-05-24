# Replay-safe helper refactor report

## Branch
- `refactor/replay-safe-output-helper`

## Scope guardrails respected
- Physics formulas changed: **No**
- Scan semantics changed: **No**
- Model construction sequence changed: **No**
- DecayTable call path changed: **No**
- Legacy columns deleted: **No**
- Historical database recomputed: **No**
- Historical artifacts overwritten/deleted: **No**

## Files changed (this refactor)
- Modified:
  - `dihiggs/Makefile`
  - `dihiggs/src/PhysScanWithFixings.cpp`
- Added:
  - `dihiggs/src/ReplaySafeOutput.hpp`
  - `dihiggs/src/ReplaySafeOutput.cpp`
  - `docs/contracts/replay_safe_scan_output_contract.md`
  - `scripts/check_replay_safe_output_smoke.py`

## Additive replay-safe contract
Created:
- `docs/contracts/replay_safe_scan_output_contract.md`

Contract states:
- `m12` is legacy/replay-unsafe unless explicit semantics are present.
- `m12_2_used` is the canonical replay value.
- precision requirement: scientific format with 17 digits where possible.
- required replay-safe/provenance fields include:
  - `m12_2_used`
  - `m12_2_gen_after_set`
  - `delta_m12_2_gen_minus_used`
  - `replay_semantics_version`
  - `yukawa_type`
  - `higgs_state`
  - `model_api_path`
  - `calc_engine`
  - `git_commit`
  - `git_dirty`

## Helper implementation summary
Added lightweight helper:
- `dihiggs/src/ReplaySafeOutput.hpp`
- `dihiggs/src/ReplaySafeOutput.cpp`

Responsibilities implemented:
- replay semantics version constant (`pswf_m12_used_v1`)
- scientific double formatting at 17 digits
- replay metadata struct with required fields
- replay-safe metadata column names provider
- CSV metadata append serializer
- optional git provenance (`DIHIGGS_GIT_COMMIT`, `DIHIGGS_GIT_DIRTY`; fallback `unknown`)

## PhysScanWithFixings integration (first engine only)
- Integrated helper only in `PhysScanWithFixings.cpp`.
- Preserved legacy `m12` column.
- Added/centralized replay-safe metadata columns through helper.
- Kept `set_param_phys_lam1 -> get_param_gen -> set_param_phys` sequence unchanged.
- Kept `DecayTable` usage unchanged.

## Columns now emitted (header)
From smoke CSV:
- `m_phi,mA,alpha,beta,lambda6,lambda7,m12,`
- `m12_2_used,m12_2_gen_after_set,delta_m12_2_gen_minus_used,replay_semantics_version,`
- `yukawa_type,higgs_state,model_api_path,calc_engine,git_commit,git_dirty,`
- `sin_ba,tan_beta,positivity_ok,unitarity_ok,perturbativity_ok,`
- `width_bb,width_tautau,width_WW,width_ZZ,width_gaga,width_Zga,width_gg,width_hh,total_width,br_gaga,`
- `lam1,computed_lam1,lam2,computed_lam2,lam3,lam4,lam5`

## Commands run
1. Build check:
   - `make app/PhysScanWithFixings` (workdir `dihiggs/`)
2. Point A one-point smoke:
   - `OMP_NUM_THREADS=1 ./dihiggs/app/PhysScanWithFixings 128.00691340782123 128.00691340782123 1 1e-12 1e-12 1 620 1.0 829846.6999999998 0.0001 0.0 scripts/out/replay_safe_helper_refactor/pointA_smoke.csv`
3. Replay-safe metadata/precision smoke check:
   - `python scripts/check_replay_safe_output_smoke.py --csv scripts/out/replay_safe_helper_refactor/pointA_smoke.csv`
4. Point A non-catastrophic confirmation against known CalcPhys reference:
   - computed `R_new_vs_calcphys_ref = 1.000132477494591`
   - computed `abs_log10_R = 5.753043421747889e-05`
   - catastrophic threshold check (`abs_log10_R < 3`) => **True (non-catastrophic)**

## Smoke artifacts
- `scripts/out/replay_safe_helper_refactor/pointA_smoke.csv`
- `scripts/out/replay_safe_helper_refactor_report.md`

## Validation result summary
- Build: **PASS**
- CSV replay-safe columns present: **PASS**
- `replay_semantics_version` present: **PASS**
- `m12_2_used` finite: **PASS**
- `delta_m12_2_gen_minus_used` finite: **PASS**
- precision (scientific, >=17 sig digits on replay metadata fields): **PASS**
- `yukawa_type == 1` in PSWF path: **PASS**
- `higgs_state` recorded (`h2/DECAY35` convention): **PASS**
- Point A remains non-catastrophic: **PASS**

## Remaining work (before extending to other engines)
1. Integrate the same helper into other scan engines (additive only).
2. Keep engine-specific physics/model sequences untouched while adding metadata columns.
3. Add cross-engine replay-safe contract checks in CI for column presence/precision.
4. Optionally inject build-time git metadata in build system (`DIHIGGS_GIT_COMMIT`, `DIHIGGS_GIT_DIRTY`) for richer provenance.
5. Validate downstream consumers tolerate additive columns without schema breakage.
