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

The output header contains exactly these 64 fields, in order:

1. `schema_version`
2. `evaluator_commit`
3. `evaluator_dirty`
4. `evaluator_api`
5. `point_id`
6. `mh_input_gev_raw`
7. `mh_input_gev`
8. `mH_input_gev_raw`
9. `mH_input_gev`
10. `mA_input_gev_raw`
11. `mA_input_gev`
12. `mHp_input_gev_raw`
13. `mHp_input_gev`
14. `sin_beta_minus_alpha_input_raw`
15. `sin_beta_minus_alpha_input`
16. `tan_beta_input_raw`
17. `tan_beta_input`
18. `lambda1_target_raw`
19. `lambda1_target`
20. `lambda6_input_raw`
21. `lambda6_input`
22. `lambda7_input_raw`
23. `lambda7_input`
24. `construction_ok`
25. `rejection_stage`
26. `rejection_reason`
27. `lambda1_reconstructed`
28. `lambda2_reconstructed`
29. `lambda3_reconstructed`
30. `lambda4_reconstructed`
31. `lambda5_reconstructed`
32. `lambda6_reconstructed`
33. `lambda7_reconstructed`
34. `lambda1_abs_residual`
35. `lambda1_residual_warning`
36. `m12_sq_reconstructed_gev2`
37. `tan_beta_reconstructed`
38. `positivity_ok`
39. `unitarity_ok`
40. `perturbativity_ok`
41. `stability_ok`
42. `triple_ok`
43. `theory_ok`
44. `width_bb_gev`
45. `width_cc_gev`
46. `width_tautau_gev`
47. `width_WW_gev`
48. `width_ZZ_gev`
49. `width_gammagamma_gev`
50. `width_Zgamma_gev`
51. `width_gg_gev`
52. `width_hh_gev`
53. `total_width_gev`
54. `br_bb`
55. `br_cc`
56. `br_tautau`
57. `br_WW`
58. `br_ZZ`
59. `br_gammagamma`
60. `br_Zgamma`
61. `br_gg`
62. `br_hh`
63. `width_ok`
64. `ctau_mm`

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
