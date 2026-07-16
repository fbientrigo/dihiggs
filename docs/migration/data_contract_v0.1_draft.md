# Canonical lambda1 evaluator data contract — v2

Status: executable contract for `dihiggs/app/Lambda1EvaluatorV2`.
Machine-readable crosswalk: `field_crosswalk_v0.1.json`.

## Input

The input header is exact and mandatory:

```text
point_id,mh_gev,mH_gev,mA_gev,mHp_gev,sin_beta_minus_alpha,tan_beta,lambda1_target,lambda6_input,lambda7_input
```

All nine numeric model inputs are required finite Float64 values. The evaluator
passes them to `THDM::set_param_phys_lam1` in that order. `point_id` passes
through unchanged; a blank ID becomes `lambda1_<16-hex>`, an FNV-1a-64 digest
of the canonical round-trip Float64 strings. A malformed row hashes its exact
physical input line instead.

Every physical line after the header produces exactly one output row, including
blank, malformed, construction-failing, and theory-rejected lines. Input-open,
empty-file, header, output-open/write, and final-rename failures are fatal.

## Serialization

Each numeric input has an exact `*_raw` field and a companion parsed Float64
field. Raw fields preserve the submitted decimal lexeme byte-for-byte. Parsed
and computed doubles use scientific notation with
`std::numeric_limits<double>::max_digits10`; parsing a persisted value therefore
reconstructs the same `double`. Digits beyond Float64 remain provenance only,
because 2HDMC consumes `double`.

The exact ordered 64-field output header is the `output_fields` array in the
machine-readable crosswalk; a test compares it directly with evaluator output.

`schema_version` is `dihiggs.lambda1.v2`. Commit and dirty provenance come from
`DIHIGGS_GIT_COMMIT` and `DIHIGGS_GIT_DIRTY`, falling back to `unknown`.

## Executable semantics

- `construction_ok` is exactly the return of `set_param_phys_lam1`.
- `triple_ok = positivity_ok && unitarity_ok && perturbativity_ok`.
- `theory_ok = triple_ok`; width validity and residual warnings do not change it.
- `stability_ok` records `Constraints::check_stability()` independently in the
  row, although the current dependency aliases it to the positivity mechanism.
- `lambda1_abs_residual` and `lambda1_residual_warning` record 2HDMC validation;
  a warning does not reject the point.
- All nine partial widths and `total_width_gev` refer to H2 (2HDMC scalar 2).
- `width_ok` means the total width is finite and strictly positive.
- When `width_ok=1`, every BR is its partial width divided by that same total;
  there is no legacy `1e-15` guard. `ctau_mm = 1.973269804e-13 / total_width_gev`.
- Failed or malformed rows carry `nan` for unavailable Float64 results and zero
  for unavailable flags.

## Mechanically resolved migration fields

The executable now resolves row preservation, stable IDs, explicit h1/h2/h3/H±
masses, explicit `mh`, target versus reconstructed lambda1, correctly named
`m12_sq`, full-width precision, dimensionless BRs, millimetre lifetime, width
validity, residual visibility, and evaluator provenance. Legacy names remain in
the crosswalk only as migration aliases.

## Scientific questions deliberately unchanged

1. Whether positivity and stability should remain aliased or require distinct
   physics implementations.
2. Whether canonical paper acceptance should include stability, width validity,
   or a lambda1 residual threshold beyond executable `theory_ok`.
3. Whether a frozen paper campaign samples `M`, `M2`, `m12_sq`, or
   `lambda1_target`, and the measure used on that coordinate.
4. Which explicit `mh` value a future paper campaign adopts. V2 accepts 125.0
   and 125.09 as distinct inputs and preserves their boundary behavior; it does
   not silently choose between them.

Those decisions require scientific review, not a serialization migration.
