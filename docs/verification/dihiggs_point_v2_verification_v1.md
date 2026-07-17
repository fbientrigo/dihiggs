# DiHiggs Point v2 verification v1

Status: canonical core and bounded engineering pilot verified; production experimental gates are not integrated.

- Baseline: `05a6217a` (PR #58 final head supplied by the plan)
- Verified implementation commit: `2b277f416caec1503becd38225a996cf0a9b42d8`
- Schema: `dihiggs.point.v2`
- Higgs mass: 125.13 GeV ([PDG 2026 listing](https://pdg.lbl.gov/encoder_listings/s126.pdf))
- Scope: engineering validation only; no campaign-level scientific conclusion

## Pilot results

| Case | Rows | construction_ok | theory_ok_v1 | SHA-256 |
|---|---:|---|---|---|
| L01_accepted_anchor | 1 | 1 | 1.00000000000000000e+00 | `c7a334a80459d25e85260f1f8e7699e1af7aae7b5e099aee4ba184d39eb9d44f` |
| L05_theory_rejected_anchor | 1 | 1 | 0.00000000000000000e+00 | `978c66e162b4846635be420070ebaef377d4b4caf010253f6e8e1b187ba39f02` |
| L06_llp_anchor | 1 | 1 | 1.00000000000000000e+00 | `a3bb0a88677602b17a081a59338b9e4b185f06acf438b255aad32f885d7d4500` |
| ordering_boundary | 2 | 0,1 | nan,0.00000000000000000e+00 | `ddf6cbd24386b3b8fb2c188e211cacdcc1591bdd537fc5174c54d1c8b9a19343` |
| construction_failure | 1 | 0 | nan | `212cc6a0a45b65da76fdf4a2552c10560e1e2e568c206b1adf97d7dc91de113e` |

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
