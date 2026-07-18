# DiHiggs Point v2 verification v1

Status: canonical core and bounded engineering pilot verified; production experimental gates are not integrated.

- Baseline: `cbb5079c8ecc012395a525ccd1dd54f8481d5be9` (PR #58 final head)
- Verified implementation commit: `6bfad7662fd87750d838bf2fe0bd7ac00ee2326a`
- Schema: `dihiggs.point.v2`
- Higgs mass: 125.13 GeV ([PDG 2026 listing](https://pdg.lbl.gov/encoder_listings/s126.pdf))
- Scope: engineering validation only; no campaign-level scientific conclusion

## Pilot results

| Case | Rows | construction_ok | theory_ok_v1 | SHA-256 |
|---|---:|---|---|---|
| L01_accepted_anchor | 1 | 1 | 1.00000000000000000e+00 | `99536756dc07908ee932ba006527963f92e2dc1cbba5b524f5eb691020b37db8` |
| L05_theory_rejected_anchor | 1 | 1 | 0.00000000000000000e+00 | `de2581a2b98a4aabfb4bddb84b851473f331e2044e09c42be80d7cd15c088a2a` |
| L06_llp_anchor | 1 | 1 | 1.00000000000000000e+00 | `09438d2d7d0655625d1508911644e283a753d5b81aa6e8e7f7a684717a95ee49` |
| ordering_boundary | 2 | 0,1 | nan,0.00000000000000000e+00 | `c5aa0bdc84045588a54811df96de3147269ae5fe111435f1d5fc251a37fc65ef` |
| construction_failure | 1 | 0 | nan | `e5b6b47086f4e4203f92b3bd0127f59eca2166c83d83e411cef1cef8bf4b008f` |

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

## Yukawa correction before/after

The reported widths are selected channels; `width_unaccounted_GeV` records closure against the total.

| Anchor | Quantity | Before (zero-Yukawa replay) | After (verified Type I) |
|---|---|---:|---:|
| L01 | total_width_GeV | 3.29976630039622283e-06 | 5.44093329247643563e-06 |
| L01 | width_bb_GeV | 0.00000000000000000e+00 | 1.68923454963338956e-06 |
| L01 | width_tautau_GeV | 0.00000000000000000e+00 | 1.64557991401996309e-07 |
| L01 | br_gammagamma | 7.20599999999999973e-03 | 4.37072735717628467e-03 |
| L01 | ctau_mm | 5.98000000000000061e-08 | 3.62671199576841017e-08 |
| L06 | total_width_GeV | 1.10178135764866426e-17 | 5.80801813354687673e-11 |
| L06 | width_bb_GeV | 0.00000000000000000e+00 | 3.92309076675650660e-11 |
| L06 | width_tautau_GeV | 0.00000000000000000e+00 | 4.14221289016167297e-12 |
| L06 | br_gammagamma | 6.85127951799999968e-01 | 4.96465589482382474e-04 |
| L06 | ctau_mm | 1.79098129599999993e+04 | 3.39749249852109442e-03 |

## Pre-existing full-suite collection failures

- mlpython/lake_pipeline/test_nulls.py: missing optional polars dependency
- tests/test_recompute_readiness.py: missing scripts.recompute_readiness
- tests/test_run_quarantine.py: missing scripts.run_quarantine

Exact commands, field results, output paths, and checksums are in the adjacent JSON report.
