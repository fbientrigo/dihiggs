# lam1 error analysis for `2026-04-12_lam1_errors_seed123456_n5000.csv`

## Dataset summary

- samples: `5000`
- warnings (`abs_error > 1e-12`): `115`
- warning rate: `0.023000`

## Absolute error statistics

- min: `0.00000000000000000e+00`
- median: `0.00000000000000000e+00`
- mean: `1.10596062615513800e-05`
- p90: `0.00000000000000000e+00`
- p99: `1.22070312500000000e-04`
- max: `3.90625000000000000e-03`

## Worst sample

- sample_index: `1974`
- attempt_index: `1974`
- lambda1_input: `2.25418284845568320e+13`
- lambda1_recomputed: `2.25418284845568359e+13`
- abs_error: `3.90625000000000000e-03`
- warning_flag: `1`
- point: `mh=122.53866342, mH=282.92991854, mA=223.47203363, mHp=391.64896657, sba=0.01351507, lambda6=0.73815877, lambda7=0.47223929, m12_2=499.58397535, tan_beta=46247.84737292`

## Interpretation

- This dataset measures the round-trip error introduced by `set_param_phys_lam1(...)` after reconstructing `m12_2` and delegating through `set_param_phys(...)`.
- `warning_flag=1` means the stored absolute difference exceeded the hard threshold `1e-12` used by `THDM::EPS`.
- Compare the max and warning rate to decide whether your preferred scan region needs a looser warning threshold or further numerical stabilization.
