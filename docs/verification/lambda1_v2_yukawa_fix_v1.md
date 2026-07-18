# Lambda1 v2 Yukawa initialization correction

- Baseline: `cbb5079`
- Corrected producer commit: `6bfad7662fd87750d838bf2fe0bd7ac00ee2326a`
- Type I is installed only after successful construction and verified with `get_yukawas_type()`.
- Width columns are selected H2 partial widths; `width_unaccounted_gev` is the closure diagnostic.

Output SHA-256: `59ac49bac86bf2c3922d008827ddadafd782b6173d3d2bc023c795cd8b2bf330`

| Case | construction | installed | theory | total width | ctau (mm) |
|---|---:|---:|---:|---:|---:|
| L01 | 1 | 1.00000000000000000e+00 | 1 | 5.44093366939228375e-06 | 3.62671174453115750e-08 |
| L05 | 1 | 1.00000000000000000e+00 | 0 | 5.43236111832007560e-06 | 3.63243488608544801e-08 |
| L06 | 1 | 1.00000000000000000e+00 | 1 | 5.80801225661750930e-11 | 3.39749593632779259e-03 |
| ordering | 1 | 1.00000000000000000e+00 | 1 | 4.02448750676615524e-06 | 4.90315798143849954e-08 |
| construction | 0 | nan | 0 | nan | nan |
