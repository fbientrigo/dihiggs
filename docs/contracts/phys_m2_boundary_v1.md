# phys_m2_boundary_v1

## Purpose

`phys_m2_boundary_v1` defines a 2HDMC physical-parameter scan where the externally controlled second scan coordinate is `M2`, not `lambda1`.

The executable contract is:

```text
Phys_M2BoundaryScan \
  mphi_min mphi_max N_mphi \
  M2_min M2_max N_M2 \
  mA sin_ba tan_beta lambda6 lambda7 \
  output.csv
```

The fixed light Higgs mass is `mh = 125.0`. For this version, `mHp = mA`.

## Controlled Variables

The scan grid is controlled by:

- `m_phi`: linearly sampled from `mphi_min` to `mphi_max` with `N_mphi` points.
- `M2_input`: linearly sampled from `M2_min` to `M2_max` with `N_M2` points.

The fixed physical inputs are:

- `mA`
- `sin_ba`
- `tan_beta`
- `lambda6`
- `lambda7`

The model uses `model.set_yukawas_type(1)`.

## Derived Variables

For each point, the engine computes:

```text
beta = atan(tan_beta)
sin_beta = sin(beta)
cos_beta = cos(beta)
sin_beta_cos_beta = sin_beta * cos_beta
m12_sq_input = M2_input * sin_beta_cos_beta
```

`m12_sq_input` is the value supplied to:

```text
THDM::set_param_phys(mh, mH, mA, mHp, sin_ba, lambda6, lambda7, m12_sq_input, tan_beta)
```

## Observed Variables

`lambda1` is not an input coordinate in this contract. The generic-basis parameters are observed after `set_param_phys(...)` through `get_param_gen(...)` and are emitted as:

- `lam1_out`
- `lam2_out`
- `lam3_out`
- `lam4_out`
- `lam5_out`
- `lam6_out`
- `lam7_out`
- `m12_sq_out`
- `tanb_out`

Constraint and decay outputs include:

- `positivity_ok`
- `unitarity_ok`
- `perturbativity_ok`
- `stability_ok`
- `width_bb`
- `width_tautau`
- `width_WW`
- `width_ZZ`
- `width_gaga`
- `width_Zga`
- `width_gg`
- `width_hh`
- `total_width`
- `br_gaga`

## Why Lambda1 Is No Longer A Scan Coordinate

Older campaigns used `lambda1` as the controlled scan coordinate and reconstructed or canonicalized `m12^2` around that choice. This contract instead scans the boundary variable:

```text
M2 = m12_sq / (sin(beta) * cos(beta))
```

That makes `M2_input` the replayed coordinate and makes `m12_sq_input` a direct derived value. `lambda1` is therefore only a model output read back from 2HDMC, not a requested value and not a round-trip target.

## Precision Policy

The trigonometric quantities `sin(beta)`, `cos(beta)`, and `sin_beta_cos_beta` are computed in `long double`. The derived value is computed as:

```text
m12_sq_input_ld = M2_input_ld * sin_beta_ld * cos_beta_ld
```

The engine casts `m12_sq_input_ld` once to `double` immediately before passing it to `set_param_phys(...)`.

CSV output uses scientific notation with `setprecision(17)` for numeric values.

## Compatibility Warning

`phys_m2_boundary_v1` is not schema-compatible with old lambda1-based campaigns. In particular, it does not emit `computed_lam1`, and `lambda1` must not be interpreted as a controlled or requested coordinate. Downstream tooling that expects `lam1`/`computed_lam1` round-trip columns needs an explicit migration path before reading this contract.
