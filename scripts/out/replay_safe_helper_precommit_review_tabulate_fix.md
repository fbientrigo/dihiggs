# Replay-safe helper pre-commit review: tabulate fallback fix

## 1) Scope
- Branch: `refactor/replay-safe-output-helper`
- Goal: make optional A-E markdown/state emission robust when `tabulate` is missing.
- Constraints respected:
  - no mandatory dependency added
  - no C++ physics code touched
  - no ReplaySafeOutput helper changes
  - no long audits run

## 2) Files changed
- `scripts/out/calcphys_three_way_catastrophic_audit.py`
  - Added resilient markdown rendering helper:
    - tries `DataFrame.to_markdown()`
    - falls back to manual markdown table if unavailable (e.g. missing `tabulate`)
  - Keeps CSV + state writing path unchanged semantically.
- New run artifacts:
  - `scripts/out/replay_safe_helper_precommit_AE_tabulate_fix/A_to_E_three_way.csv`
  - `scripts/out/replay_safe_helper_precommit_AE_tabulate_fix/A_to_E_three_way.md`
  - `scripts/out/replay_safe_helper_precommit_AE_tabulate_fix/state.json`
  - `scripts/out/replay_safe_helper_precommit_AE_tabulate_fix/run_commands.sh`

## 3) Commands run
1. Patch reporting script to add fallback markdown table.
2. Re-run only A-E precommit runner:
   - `source ~/higgs_env_py312/bin/activate && python scripts/out/calcphys_three_way_catastrophic_audit.py --mode A_to_E --outdir scripts/out/replay_safe_helper_precommit_AE_tabulate_fix`
3. Verify required artifacts and catastrophic counts:
   - existence checks for CSV/MD/state JSON
   - parse CSV and state JSON counts

## 4) Validation results
- Runner exit code: `0` (success)
- Observed warning (expected in env without tabulate):
  - `pandas.to_markdown unavailable; falling back to manual markdown table`
- Artifact checks:
  - A-E CSV exists: PASS
  - A-E markdown exists: PASS
  - A-E state JSON exists: PASS
- A-E catastrophic counts (from CSV):
  - old catastrophic = `5`
  - new catastrophic = `0`
- State JSON counters:
  - `catastrophic_old_vs_calcphys = 5`
  - `catastrophic_new_vs_calcphys = 0`

## 5) Risk check
- Change is additive and localized to reporting robustness.
- No physics formulas, scan semantics, or replay-safe C++ logic modified.
- Fallback only affects markdown rendering path when optional dependency is absent.

## 6) Recommendation
`READY_TO_COMMIT`
