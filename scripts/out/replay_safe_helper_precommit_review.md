# Replay-safe helper pre-commit review

## Git branch
- `refactor/replay-safe-output-helper`

## Git status summary
- Repository status at review time:
  - modified entries: 6
  - untracked entries: 36
  - total status entries: 42
- Scoped status for this refactor/review:
  - `M dihiggs/Makefile`
  - `M dihiggs/src/PhysScanWithFixings.cpp`
  - `?? dihiggs/src/ReplaySafeOutput.cpp`
  - `?? dihiggs/src/ReplaySafeOutput.hpp`
  - `?? docs/contracts/replay_safe_scan_output_contract.md`
  - `?? scripts/check_replay_safe_output_smoke.py`
  - `?? scripts/out/replay_safe_helper_refactor_report.md`
  - `?? scripts/out/replay_safe_helper_precommit_AE/` (review artifact dir)

## Files changed (scoped)
- `docs/contracts/replay_safe_scan_output_contract.md`
- `dihiggs/src/ReplaySafeOutput.hpp`
- `dihiggs/src/ReplaySafeOutput.cpp`
- `dihiggs/src/PhysScanWithFixings.cpp`
- `dihiggs/Makefile`
- `scripts/check_replay_safe_output_smoke.py`
- `scripts/out/replay_safe_helper_refactor_report.md`

## Review checklist (pass/fail)

| # | Check | Result | Evidence / note |
|---|---|---|---|
| 1 | ReplaySafeOutput helper is additive-only | PASS | New helper files only add metadata/formatting utilities. |
| 2 | Helper does not own/alter physics calculations | PASS | No THDM/Constraints/DecayTable physics logic in helper. |
| 3 | PhysScanWithFixings model construction sequence unchanged | PASS | Sequence remains `set_param_phys_lam1 -> get_param_gen -> set_param_phys (canonical) -> get_param_gen`. |
| 4 | DecayTable calls unchanged | PASS | Same `DecayTable tab(canonical_model)` and same width accessors (`hdd/hll/hvv/hgaga/hZga/hgg/hhh/gammatot`). |
| 5 | Legacy `m12` preserved | PASS | Header still includes `m12`; row still emits legacy slot before replay metadata block. |
| 6 | `m12_2_used` captured from exact value used for set_param_phys / decay | PASS | `m12_2_used` taken from `model.get_param_gen(...)` and passed directly to canonical `set_param_phys(...)`; DecayTable uses canonical model. |
| 7 | `m12_2_gen_after_set` captured after model setup | PASS | Captured from `canonical_model.get_param_gen(...)` after successful canonical setup. |
| 8 | `delta_m12_2_gen_minus_used` computed correctly | PASS | Helper computes `m12_2_gen_after_set - m12_2_used`. |
| 9 | `replay_semantics_version` always emitted for valid rows | PASS | `make_metadata(...)` always sets default `pswf_m12_used_v1`; metadata serialization is mandatory on each emitted row. |
| 10 | `yukawa_type` emitted and equals 1 for this path | PASS | PSWF path sets `canonical_model.set_yukawas_type(1)` and emits `canonical_model.get_yukawas_type()`. Smoke output shows 1. |
| 11 | `higgs_state` convention documented and emitted | PASS | Contract documents `h2/DECAY35`; helper default emits `h2/DECAY35`. |
| 12 | `calc_engine` and `model_api_path` emitted | PASS | Both fields in metadata columns and CSV serialization. |
| 13 | `git_commit`/`git_dirty` safe when env vars absent | PASS | Fallback to `"unknown"` if env vars missing. |
| 14 | Floating formatting scientific + 17 digits | PASS | `configure_scientific_17(...)` used for row stream; metadata doubles serialized with scientific 17-digit formatter. |
| 15 | Makefile change does not break partial/other app targets | PASS | `make app/PhysScanWithFixings` and `make app/PhysParamScan` succeed; adding helper to `COMMON_SRCS` links cleanly. |
| 16 | Smoke script checks required columns robustly | PASS | Validates required columns, non-empty semantics, finite replay fields, precision threshold, expected yukawa/higgs_state, non-empty provenance fields. |
| 17 | No historical artifacts overwritten/deleted | PASS | Refactor edits only code/docs/new scripts; review-time A-E run used isolated outdir `scripts/out/replay_safe_helper_precommit_AE`. No deletes performed. |

## Commands run
1. `git branch --show-current && git status --short --branch`
2. `git diff --name-only -- docs/contracts/replay_safe_scan_output_contract.md dihiggs/src/ReplaySafeOutput.hpp dihiggs/src/ReplaySafeOutput.cpp dihiggs/src/PhysScanWithFixings.cpp dihiggs/Makefile scripts/check_replay_safe_output_smoke.py scripts/out/replay_safe_helper_refactor_report.md`
3. `make app/PhysScanWithFixings` (workdir: `dihiggs/`)
4. `source ~/higgs_env_py312/bin/activate && python scripts/check_replay_safe_output_smoke.py --csv scripts/out/replay_safe_helper_refactor/pointA_smoke.csv`
5. Optional A-E runner attempt:
   - `source ~/higgs_env_py312/bin/activate && python scripts/out/calcphys_three_way_catastrophic_audit.py --mode A_to_E --outdir scripts/out/replay_safe_helper_precommit_AE`
   - Result: process completed recomputation artifacts/CSV, then failed at markdown generation due missing Python optional dependency `tabulate`.
6. `make app/PhysParamScan` (workdir: `dihiggs/`) for cross-target sanity.

## A-E result
- Existing A-E runner was available and executed point recomputation.
- Output CSV produced:
  - `scripts/out/replay_safe_helper_precommit_AE/A_to_E_three_way.csv`
- Parsed CSV summary:
  - rows: 5 (A,B,C,D,E)
  - `catastrophic_old_vs_calcphys`: 5
  - `catastrophic_new_vs_calcphys`: 0
  - `max_log10_abs_new_vs_calcphys`: `5.753035335351849e-05`
- Runner final markdown/state write failed because `tabulate` is not installed in `~/higgs_env_py312`.

## Remaining risks
1. Git provenance fields currently default to `unknown` unless env vars are set at runtime (`DIHIGGS_GIT_COMMIT`, `DIHIGGS_GIT_DIRTY`).
2. A-E audit script has a non-critical dependency on `tabulate` for markdown emission; CSV evidence is still produced.
3. Repository has many unrelated dirty/untracked files; commit should be path-scoped to intended refactor files.

## Recommendation
**READY_TO_COMMIT**

Rationale: all 17 scoped review checks passed; required short validations passed; optional A-E replay evidence is positive (new catastrophic=0) despite optional markdown dependency issue.
