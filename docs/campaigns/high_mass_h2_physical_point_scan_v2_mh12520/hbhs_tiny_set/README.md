# HB/HS enrichment — tiny 3-point set

Proof that the `HBHS_BLOCKED` blocker is removed. **Three points only**; this is
not a rescan of the 208-point set.

## Before

All 208 accepted high-mass points were classified `HBHS_BLOCKED` — *attempted,
not skipped*. `dihiggs_boundary`'s only point producer hard-coded
`mh = 125.09` and `sin_ba = 1.0`, and took `M` (so `M2 >= 0`), while the
campaign runs at `sin_ba in {0.995, 0.999, 1.0}` with negative `M2` for half the
set. No campaign point could be expressed in its input contract.

## After

`src/evaluate_point.cpp` gained a canonical v2 input schema that reads
`mh, mHp, sin_ba, lambda7, M2` off the point. Result: `total=3 enriched=3
skipped=0 errors=0`.

| point | m_H2 | sin_ba | hb_allowed | hb_max_obsratio | limiting process | hs_chi2 | hs_delta_chi2 | exp_ok |
|---|---:|---:|---|---:|---|---:|---:|---|
| `canonical_anchor_150` | 150 | 1.000 | 0 | 23.36 | `LHC13 [ggH]->X1->Z(X2->[bb])` | 156.81 | 5.17 | 0 |
| `point_feb838bf5e5e421f` | 200 | 0.995 | 0 | 3.01 | `LHC13 [vbfH;H]>[ZZ]` | 166.65 | 15.01 | 0 |
| `point_61bf5e4c3200c037` | 250 | 0.995 | 0 | 4.75 | `LHC13 [ggH>ZZ]` | 165.88 | 14.24 | 0 |

`hs_chi2_sm_ref = 151.638` (159 observables).

**Physics reading**: all three are theory-valid but **HiggsBounds-excluded**,
each by a named LHC13 heavy-scalar search, with observed/expected ratios of
3–23. That is a real experimental verdict where previously there was none.
The 200 and 250 GeV points are the mission-chain points
(`point_011af49c15561f96`, `point_0960bd688a662562` in v1) recalculated at the
canonical convention.

## SM-reference sensitivity to the convention

`dhb.runner.SM_REFERENCE_MH` now follows the canonical convention. HiggsSignals'
own mass anchor is `125.25 +- 0.17` (`hsdataset-v1.1/h125/mass_LHC13_LHCComb_36.json`),
so the move is *toward* the data:

| mh (GeV) | SM-reference chi2 |
|---:|---:|
| 125.09 (old boundary constant) | 152.549 |
| 125.13 (old upstream convention) | 152.050 |
| **125.20 (canonical)** | **151.638** |
| 125.25 (HS anchor) | 151.552 |

Nothing in the HS dataset requires 125.09; this is a numerics shift, not a
compatibility break.

## Reproduce

    cd dihiggs_boundary
    bash scripts/build_evaluate_point.sh
    ./build/bin/evaluate_point tiny_points.csv evaluate_point.csv
    PYTHONPATH=python python3 -m dhb.enrich \
      --input evaluate_point.csv --output hbhs_enriched.csv \
      --config configs/theory_atlas_v0.yaml
