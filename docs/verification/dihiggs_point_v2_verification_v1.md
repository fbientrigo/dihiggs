# DiHiggs Point v2 verification v1

Status: canonical core and bounded engineering pilot verified; production experimental gates are not integrated.

- Baseline: `05a6217a` (PR #58 final head supplied by the plan)
- Successor: uncommitted successor worktree
- Schema: `dihiggs.point.v2`
- Higgs mass: 125.13 GeV ([PDG 2026 listing](https://pdg.lbl.gov/encoder_listings/s126.pdf))
- Scope: engineering validation only; no campaign-level scientific conclusion

## Pilot results

| Case | Rows | construction_ok | theory_ok_v1 | SHA-256 |
|---|---:|---|---|---|
| L01_accepted_anchor | 1 | 1 | 1.00000000000000000e+00 | `43bc772bef1b68c94460d34e61504319936e80406b537e0a8f6b81398d88a3de` |
| L05_theory_rejected_anchor | 1 | 1 | 0.00000000000000000e+00 | `df20dd9195876436efacba4adfa95ed4cc53a90dbc372914b3abb4cd5a224f19` |
| L06_llp_anchor | 1 | 1 | 1.00000000000000000e+00 | `9bbc77fa5f7e2e5448ac235824423c0bbd7c901f5c5bde2212b650592ccfad97` |
| ordering_boundary | 2 | 0,1 | nan,0.00000000000000000e+00 | `0f3059e197dc9cd9a4c64cfbf20e982489ad3a4469acbee9191e11ef56b5df96` |
| construction_failure | 1 | 0 | nan | `631129a5714a59e5a9e1c014b81ad630324f3547027ec860b5475e25f5e5b58b` |

## Readiness dispositions

- D01-D04: verified by focused executable tests
- Q01: quarantined and non-gating
- RV01: verified for canonical producer
- RV02: verified only in the tested pilot domains; boundaries are not globally validated
- RV03: verified for canonical producer
- RV04: deferred as approved
- RV05: verified not production-integrated; experimental fields remain unevaluated

## Deferred

full scans, HB/HS datasets, recasting, training, plotting, maintained TeX claims.

## Pre-existing full-suite collection failures

- mlpython/lake_pipeline/test_nulls.py: missing optional polars dependency
- tests/test_recompute_readiness.py: missing scripts.recompute_readiness
- tests/test_run_quarantine.py: missing scripts.run_quarantine

Exact commands, field results, output paths, and checksums are in the adjacent JSON report.
