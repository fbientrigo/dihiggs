# High-mass H2 physical point search — report

Producer: `DihiggsPointV2Evaluator` (commit `2264ffe1b8d8f59bec69341d6aa1a917396ccb2b`,
post-PR#80, includes the `width_tt_GeV`/`br_tt` Gate A fix). Experiment class:
`PHYSICAL_POINT_SCAN` (tan_beta, M2, lambda6 allowed to vary in search of
theory-valid points, as distinct from the fixed-parameter `MASS_CONTROL`
pilots in `docs/pilots/high_mass_h2_v1/`). Model variant:
`PHYSICAL_DECAYS_NO_HEAVY_CASCADES`. `m_h` frozen at 125.13 GeV throughout,
per `HIGH_MASS_H2_CONTRACT.md` section 3.

Search script: `scripts/high_mass_physical_point_search.py`. Total evaluator
CLI invocations: 1080 (9 mass regions x 3 sin(beta-alpha) x 10 tan_beta x 4
lambda6 combinations), each internally gridding M2 over 51 points as a
fraction of `m_A^2` in `[-1, 3]`. Total individual points evaluated: 55,080.
Wall time: 458.9 s (6 concurrent worker threads).

## Regions searched

| Region | m_H2 (GeV) | Delta_heavy (GeV) | m_A = m_Hp (GeV) | Purpose |
|---|---:|---:|---:|---|
| R150_anchor | 150 | 300 | 450 | regression-adjacent anchor |
| R200_below_hh | 200 | 300 | 500 | below H2->hh threshold (250.26) |
| R250_near_hh_threshold | 250 | 300 | 550 | just below H2->hh threshold |
| R300_above_hh | 300 | 500 | 800 | above H2->hh threshold |
| R350_near_tt_threshold | 350 | 450 | 800 | near H2->tt threshold (345) |
| R400_above_tt | 400 | 400 | 800 | above H2->tt threshold |
| R500 | 500 | 500 | 1000 | mid-range |
| R600 | 600 | 600 | 1200 | mid-high |
| R800_near2000 | 800 | 1200 | 2000 | upper edge of both mass targets |

Fixed across all regions: Yukawa Type I, `lambda7 = 0.0` exactly. Varied per
attempt: `sin(beta-alpha) in {1.0, 0.999, 0.995}`, `tan_beta in {2, 5, 10,
30, 100, 300, 1e3, 1e4, 1e5, 3e5}`, `lambda6 in {0.0, 1e-10, 0.1, -0.1}`,
`M2` swept as a fraction of `m_A^2` in `[-1, 3]` (51 points).

## Result: theory-valid points found only up to ~250 GeV in this grid

**208 of 55,080 attempted points (0.38%) were accepted** — i.e. passed the
evaluator's own `theory_ok` (positivity, unitarity, perturbativity,
stability all true) AND have all four forbidden cascade flags
(`H2_to_AZ_open`, `H2_to_HpW_open`, `H2_to_AA_open`, `H2_to_HpHm_open`)
false. All 208 accepted points fall in the three lowest-mass regions:

| Region | Accepted | tan_beta values found | sin(beta-alpha) values found | lambda6 values found | M2 range (GeV^2) |
|---|---:|---|---|---|---|
| R150_anchor | 77 | {2, 5} | {0.995, 0.999, 1.0} | {-0.1, 0, 1e-10, 0.1} | [-8.9e4, 2.4e4] |
| R200_below_hh | 62 | {2, 5} | {0.995, 0.999, 1.0} | {-0.1, 0, 1e-10, 0.1} | [-7.0e4, 5.0e4] |
| R250_near_hh_threshold | 69 | {2, 5, 10} | {0.995, 0.999, 1.0} | {-0.1, 0, 1e-10, 0.1} | [-3.6e4, 6.1e4] |
| R300_above_hh through R800_near2000 (6 regions) | **0** | — | — | — | — |

**No theory-valid point was found at or above 300 GeV anywhere in this
search grid.** This is recorded as a genuine result, not a search failure to
paper over: every one of the 6120 attempted points per region above 250 GeV
was rejected, overwhelmingly at the positivity (vacuum-stability /
boundedness-from-below) stage, with unitarity a secondary rejector and
perturbativity a minor contributor (see `rejection_summary.csv` for exact
counts per region/stage/reason). This says the specific
`(tan_beta, sin(beta-alpha), lambda6, M2-as-fraction-of-mA^2)` grid explored
here does not contain a theory-valid point once `m_A = m_Hp` grows past
~550-800 GeV at these `Delta_heavy` values — not that no such point can
exist in the full 2HDM parameter space. A denser or differently-shaped grid
(e.g. exploring `lambda6` magnitudes beyond +-0.1, or `M2` fractions outside
`[-1,3]`) was out of this task's time-box; see "Next steps" below.

Interesting contrast with the pre-existing `MASS_CONTROL` pilot anchor (P0,
150 GeV, `docs/pilots/high_mass_h2_v1/pilot_points.csv`): that anchor used
near-decoupling (`tan_beta = 300000`), while this `PHYSICAL_POINT_SCAN`
found its accepted points instead cluster at **low** `tan_beta in {2, 5,
10}` with `sin(beta-alpha)` close to (but not necessarily exactly at)
alignment. These are two different, both theory-valid, corners of parameter
space for a similar mass point — a useful independent confirmation that the
evaluator is not simply reproducing one hand-tuned solution.

## Was H2->tt or H2->hh important?

**Not assessable from the accepted-point sample in this run**, because every
accepted point has `m_H2 <= 250 GeV`, which is below both the `H2->hh`
threshold (250.26 GeV, confirmed: `H2_to_hh_open = False` for all 208
accepted rows, `BR_hh = 0.0` for all of them) and far below the `H2->tt`
threshold (345 GeV; `BR_tt = 0.0` for all 208). The regions specifically
targeting these thresholds (`R300_above_hh`, `R350_near_tt_threshold`,
`R400_above_tt`) found zero theory-valid points in this grid, so this
question remains open pending a wider/denser search at those masses — this
is explicitly flagged as the most valuable next step (see below), not
glossed over.

For contrast, the earlier `MASS_CONTROL` pilot (not theory-tuned, so not
directly comparable, but informative about the evaluator's physics) *did*
show `BR_tt` rising through the threshold and `BR_hh` becoming visible above
250 GeV (see `docs/HIGH_MASS_H2_CONTRACT.md` section 7's pilot table:
`br_hh` ~0.29-0.33, `width_tt_GeV` rising from 0 to 2.03e-2 across
`P1`-`P5`), so both channels are known to be physically significant once
kinematically open — the open question is only whether a *theory-valid*
point also exists up there, which this particular grid did not find.

## Width closure and data quality

All 208 accepted points pass `explicit_width_sum + width_unaccounted ==
total_width` to `rel_tol=1e-9` (checked directly, zero violations). No
significant width is hidden in `width_unaccounted_GeV`.

## HB/HS enrichment: attempted, blocked (not skipped)

The `dihiggs_boundary` HB/HS pipeline was inspected read-only
(`python/dhb/runner.py`, `enrich.py`, `adapter.py`, `schema.py`,
`src/evaluate_point.cpp`, `docs/evaluate_point_contract.md`,
`tests/test_golden_evaluate_point.py`) to determine whether it could enrich
these 208 points. Its only tested, working point-producer
(`evaluate_point.cpp`, wrapped by `HbhsRunner`) hard-codes:

```cpp
constexpr double kMh = 125.09;   // this campaign requires 125.13, frozen (HIGH_MASS_H2_CONTRACT.md sec 3)
constexpr double kSinBa = 1.0;   // this campaign's accepted points include sin_ba in {0.995, 0.999}
```

(independently confirmed at `dihiggs_boundary/src/evaluate_point.cpp:21-22`
by this report's author, not merely quoted from the search agent). Feeding
this campaign's points through that pipeline would require either modifying
`dihiggs_boundary` (out of scope — read-only access per mission constraints)
or writing a new C++ program that reimplements `evaluate_point.cpp`'s 2HDMC
effective-coupling extraction with `mh=125.13` and general `sin(beta-alpha)`
— a nontrivial reimplementation, not a small conversion, and outside this
task's time-box.

Per the mission's explicit instruction ("do not block the entire mission if
HB/HS infrastructure itself has an unrelated operational issue"), all 208
accepted points are classified **`HBHS_BLOCKED`**, not silently dropped and
not miscategorized as `THEORY_VALID`-and-nothing-more (their evaluator-level
theory validity is still separately recorded and true) or as
experimentally-allowed (never claimed).

## Classification summary

```
HBHS_BLOCKED: 208   (all accepted points; theory-valid, HB/HS attempted and blocked)
```

No point in this run reached `THEORY_AND_HBHS_VALID` or plain
`HBHS_NOT_EVALUATED` (HB/HS was attempted for all of them, not skipped).

## Files in this directory

- `high_mass_valid_points.csv` — 208 accepted points, full physics record.
- `high_mass_attempted_points.csv` — all 55,080 attempted points (accepted and rejected).
- `rejection_summary.csv` — counts by mass_region x rejection_stage x rejection_reason.
- `mass_validity_map.png` / `.pdf` — diagnostic scatter (m_H2 vs m_A, colored by rejection stage) + acceptance-rate bar chart per region.
- `manifest.json` — machine-readable summary (search config hash, evaluator provenance, classification counts, HB/HS blocker detail).
- `checksums.sha256` — sha256 of every file in this directory.

## Next steps (not executed in this task, time-boxed out)

1. Widen/refine the search grid specifically for `m_H2 >= 300 GeV` — the
   positivity rejection dominance suggests either different `lambda6`
   magnitudes or a different `M2`-fraction range is needed as `m_A` grows;
   this was not explored here.
2. If a theory-valid point is found above 345 GeV, re-run the HB/HS blocker
   assessment — it may still be blocked by the same `mh`/`sin_ba` mismatch,
   but should be re-confirmed per-point rather than assumed.
3. Feed a representative subset of the 208 accepted (150-250 GeV) points
   into the `model_point_to_llp_recast_v2` contract and the generic
   FACTORIZED_G_ONLY UFO builder (parallel tasks in this mission) to
   validate the full downstream handoff, independent of whether the
   300+ GeV search is ever widened.
