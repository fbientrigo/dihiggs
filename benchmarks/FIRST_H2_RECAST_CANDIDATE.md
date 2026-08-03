# First H2 benchmark freeze: `H2scan_mH150_tb300000`

## Scoped verdict

The fresh 15-point campaign and corrected numerical check pass the scoped
handoff gate for the first `H2H2 -> 4b` downstream chain.

| Gate | Result |
|---|---|
| `channel_validity.bb_dvjets` | `VALID_FOR_FIRST_BB_RECAST` |
| `scan_provenance_status` | `VALIDATED` |
| `numerical_check_status` | `PASS` |

Numerical stability for `gammagamma` and `Zgamma` is recorded below as a
representation-scale classification only; it is not downstream authorization
for those channels. Production normalization remains pending MadGraph.

## Selected row

| Field | Value |
|---|---:|
| `point_id` | `H2scan_mH150_tb300000` |
| `m_H2_GeV` | `150.0` |
| `mA_GeV` | `450.0` |
| `mHp_GeV` | `450.0` |
| `tan_beta` | `300000.0` |
| `sin_beta_minus_alpha` | `1.0` |
| `lambda1_target` | `1.0` |
| `lambda6_input` | `1.00000000000000004e-10` |
| `lambda7_input` | `0.0` |
| `yukawa_type` | `Type I` |
| `construction_ok` / `theory_ok` / `width_ok` | `1 / 1 / 1` |
| `total_width_GeV` | `4.56118529862185007e-14` |
| `ctau_mm` | `4.32622152973311191` |
| `br_bb` | `0.756737485808578692` |
| `br_gammagamma` | `0.0002657929417704355` |
| `br_Zgamma` | `0.0000204915508605281493` |
| `lambda1_abs_residual` | `9.34182041500974947e-7` |
| `lambda1_residual_warning` | `1` |

## Soft scale export

All soft-scale values are in `GeV²` and were derived from the selected fresh
CSV row, not hardcoded.

| Field | Value |
|---|---:|
| `source_soft_scale_field` | `m12_sq_reconstructed_gev2` |
| `source_soft_scale_value` | `0.07499975253962432` |
| `m12_sq_GeV2` | `0.07499975253962432` |
| `M2_GeV2` | `22499.925761493243` |
| `soft_scale_relation` | `m12_sq = M2 * sin(beta) * cos(beta)` |
| `beta` | `atan(tan_beta)` |
| `soft_scale_consistency_relative_error` | `0.0` |

## Fresh scan provenance

| Field | Value |
|---|---|
| `scan_producer_commit` | `31adde89d831195adde927b364046723ba29e3fe` |
| `scan_manifest_schema` | `dihiggs.h2_bounded_scan.manifest.v2` |
| `scan_execution_mode` | `FRESH_FROM_EMPTY_OUTPUT` |
| `scan_resumed_rows` | `0` |
| `scan_row_count` | `15` |
| runner SHA-256 | `c8b8fdff3a18121d9b165b8d23fab4253bb5a9e9cf582fe04b96f21e7e46e859` |
| evaluator source SHA-256 | `5e3a81ff53e1778d1700824acbe20d8808aaefdd46a349ef5339d64ddf795eb4` |
| ReplaySafeOutput source SHA-256 | `77102e8db53b0fbf8bc88abe92bb10ce44d82420cd9a975c7f6c1dddca0615d4` |
| evaluator binary SHA-256 | `0685693d60b2bf55c6c0ab2ebe3f7fa9f8da2fb0dc3421cb2ed1fdc68a01b31a` |
| linked `lib2HDMC` | `2hdmc/lib/lib2HDMC.a` |
| linked `lib2HDMC` SHA-256 | `4dc3443da9e5c86f0b37b924b92f3b2f03d71bfb524e3d96a7d387630fd76742` |

The evaluator is statically linked to the rebuilt `lib2HDMC.a`; `ldd` showed
no `lib2HDMC.so` entry.

## Corrected numerical evidence

The checker uses the exact center and its immediately adjacent representable
double values only.

| Field | Value |
|---|---:|
| center reproduction tolerance | `1.00000000000000002e-08` |
| maximum center reproduction difference | `0.00000000000000000e+00` |
| field producing maximum center difference | `total_width_gev` |
| all values finite | yes |
| all relevant values physically valid | yes |
| all adjacent probes theory-valid | yes |
| adjacent-float maximum spread | `1.66564239607794041e-06` |
| field producing maximum adjacent spread | `br_Zgamma` |
| final classification | `STABLE_AT_DOUBLE_REPRESENTATION_SCALE` |

Numerical classifications: `bb_dvjets`, `gammagamma`, and `Zgamma` are each
`STABLE_AT_DOUBLE_REPRESENTATION_SCALE`. Only `bb_dvjets` receives downstream
validity authorization in this freeze.

## Lifetime and production status

`python3 benchmarks/verify_pilot_ctau_invariant.py` passed for all 15 fresh
benchmark rows and all checked pilot rows using the unchanged
`1.973269804e-13 GeV mm` constant and `1e-4` relative tolerance.

| Field | Value |
|---|---|
| `sigma_production_fb` | `null` |
| `sigma_status` | `PENDING_MADGRAPH` |
