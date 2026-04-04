# Comparison: baseline vs trig-identity rewrite (internal `THDM.cpp` rewrite)

## Goal

Evaluate whether the trig-identity rewrite in `set_param_phys_lam1` improves numerical behavior relative to the previous implementation.

## Same benchmark setup

- Deterministic validator: `CalcRoundTrip` (generic + near-alignment + invalid guard)
- Random stress validator: fixed seed, 400 random valid points + 3 invalid guard cases
- Same tolerance policy and observable set (`m12_2`, `tan_beta`, `lambda1..7`, `hmass(1..4)`).

## Baseline (previous implementation)

- Deterministic validator: PASS
- Random stress validator: PASS
- Global max abs diff: `2.6193447411060333e-09`
- Global max rel diff: `2.2463382265407876e-12`

## Current trig rewrite (new implementation)

### Deterministic validator (`CalcRoundTrip`)

```text
Case generic: FAIL
  m12_2: ref=20000 new=-57677.5
  tanb : ref=2 new=2
Case near-alignment: FAIL
  m12_2: ref=20000 new=-56239.7
  tanb : ref=2 new=2
Case invalid-input: PASS
```

Exit code: `1`

### Random stress validator

First failure (point 0):

```text
FAIL: Point 0 m12_2 mismatch: ref=109048.89997870798 new=149612707.64933252
abs_diff=149503658.74935383 rel_diff=0.99927112541647001
```

Exit code: `2`

## Regression magnitude

- Baseline global max abs diff: `2.6193447411060333e-09`
- Current first failure abs diff: `149503658.74935383`
- Regression factor (abs diff): `~5.707673999651724e+16` worse than baseline.

## Root cause analysis (line-level)

In `src/THDM.cpp`, the comment states correct identities:

- `sin(beta) = tan_beta / sqrt(1 + tan_beta^2)`
- `cos(beta) = 1 / sqrt(1 + tan_beta^2)`

But implementation currently assigns them swapped:

```cpp
double cb = tan_beta/sqrt(1.+tan_beta*tan_beta);  // this is sin(beta), not cos(beta)
double sb = 1./sqrt(1.+tan_beta*tan_beta);        // this is cos(beta), not sin(beta)
```

This inversion propagates to `sa`/`ca` reconstruction, causing wrong `ca2/sa2` and therefore wrong `m12_2` reconstruction.

## Conclusion

Current trig rewrite is a **strong regression** (fails both deterministic and stress tests).
The observed error is not random FP noise; it is a deterministic algebraic assignment bug (`sb`/`cb` swapped).

## After fixing the `sb`/`cb` swap (rerun)

Rerun with the same deterministic and random harnesses:

- Deterministic validator: PASS (all 3 cases), exit code `0`
- Random stress validator: PASS, exit code `0`
- New global max abs diff: `3.6088749766349792e-08`
- New global max rel diff: `1.0093685887093751e-11`

Comparison:

- vs broken trig rewrite: **massive recovery** (orders of magnitude better; failures removed)
- vs original baseline: precision is still worse
  - abs diff factor: `13.777777777777779x` worse
  - rel diff factor: `4.493395414740085x` worse

Interpretation: the correction fixed the algebraic bug and restored correctness, but the trig-only variant remains less precise than the original baseline implementation.

## Evidence files

- `.sisyphus/evidence/current-calcroundtrip-run.txt`
- `.sisyphus/evidence/current-calcroundtrip-run-exit.txt`
- `.sisyphus/evidence/current-metrics-run.txt`
- `.sisyphus/evidence/current-metrics-run-exit.txt`
- `.sisyphus/evidence/current-vs-baseline-comparison.json`
- `.sisyphus/evidence/current-trig-identity-check.txt`
- `.sisyphus/evidence/rerun-fix-calcroundtrip-run.txt`
- `.sisyphus/evidence/rerun-fix-calcroundtrip-run-exit.txt`
- `.sisyphus/evidence/rerun-fix-metrics-run.txt`
- `.sisyphus/evidence/rerun-fix-metrics-run-exit.txt`
- `.sisyphus/evidence/rerun-fix-vs-baseline.json`
