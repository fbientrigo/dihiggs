# CalcRoundTrip metrics (commit 2 support)

This file is intended to be committed with the validator work so reviewers can see quantitative outcomes without opening transient `.sisyphus` artifacts.

## Test scope

1. Deterministic validator (`src/CalcRoundTrip.cpp`):
   - generic point
   - near-alignment point
   - invalid-input guard case
2. Extended random validator (400 random valid points + 3 invalid guard points) for numerical stress testing.

## Deterministic results

Output:

```text
Case generic: PASS
Case near-alignment: PASS
Case invalid-input: PASS
All tests passed
```

Exit code: `0`

## Random stress results (400 points)

Output summary:

```text
PASS: 400 random round-trip points validated
PASS: invalid-input guards validated (3 cases)
Global max abs diff: 2.6193447411060333e-09 in m12_2 at point 66
Global max rel diff: 2.2463382265407876e-12 in hmass(1) at point 147
```

Per-observable maxima:

- `m12_2`: max_abs = `2.6193447411060333e-09`, max_rel = `7.3971567400955148e-15`
- `tan_beta`: max_abs = `0`, max_rel = `0`
- `lambda1`: max_abs = `4.6566128730773926e-10`, max_rel = `2.3350149538672452e-14`
- `lambda2`: max_abs = `1.1368683772161603e-13`, max_rel = `7.915361206189481e-16`
- `lambda3`: max_abs = `3.4106051316484809e-13`, max_rel = `5.6186457402196685e-15`
- `lambda4`: max_abs = `2.2737367544323206e-13`, max_rel = `6.1498595353138059e-15`
- `lambda5`: max_abs = `2.2737367544323206e-13`, max_rel = `4.6634475196139837e-15`
- `lambda6`: max_abs = `0`, max_rel = `0`
- `lambda7`: max_abs = `0`, max_rel = `0`
- `hmass(1)`: max_abs = `2.8079227831767639e-10`, max_rel = `2.2463382265407876e-12`
- `hmass(2)`: max_abs = `2.8194335754960775e-11`, max_rel = `2.42076070055429e-14`
- `hmass(3)`: max_abs = `1.6709122974134516e-10`, max_rel = `6.8135169195459164e-13`
- `hmass(4)`: max_abs = `1.6217427400988527e-10`, max_rel = `6.4182559634008434e-13`

## Source of observed discrepancies/errors

No physics consistency failures were found. The comparison itself passes all tolerances.

Observed "error-like" signals and root cause:

1. **`WARNING: Cannot set physical masses such that m_H < m_h`**
   - Source: intentionally injected invalid guard test case (`m_h > m_H`) in the stress validator.
   - Nature: expected behavior, not a defect.

2. **Non-zero numeric deltas (up to `~2.6e-9` absolute)**
   - Source: floating-point roundoff from inversion + trigonometric evaluation path (`asin/cos/tan`) and re-derivation of derived quantities.
   - Nature: expected double-precision behavior; relative errors remain tiny (`<= 2.25e-12`) and well within tolerance thresholds used by validators.

3. **`make build` target failure (outside validator logic)**
   - Source: this Makefile defines `all`, not `build`.
   - Nature: build-system target mismatch, not a physics/model issue.

## Evidence files (transient)

- `.sisyphus/evidence/rigorous-calcroundtrip-postbuild.txt`
- `.sisyphus/evidence/rigorous-calcroundtrip-postbuild-exit.txt`
- `.sisyphus/evidence/rigorous-metrics-run.txt`
- `.sisyphus/evidence/rigorous-metrics-build-exit.txt`
- `.sisyphus/evidence/rigorous-metrics-run-exit.txt`
- `.sisyphus/evidence/rigorous-random-after-relink.txt`
- `.sisyphus/evidence/rigorous-random-after-relink-exit.txt`
