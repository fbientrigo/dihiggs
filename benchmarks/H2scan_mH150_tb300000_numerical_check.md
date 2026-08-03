# Numerical representation check: `H2scan_mH150_tb300000`

This check uses only the exact `m12_2` double and its immediately adjacent representable values.
The previous `1e-12 GeV^2` offsets changed reconstructed `lambda1` by order one and therefore compared different physical models; their channel classifications are withdrawn.

- construction `m12_sq`: `7.49999999996479594e-02` GeV^2
- construction `M2`: `2.24999999995003345e+04` GeV^2
- round-trip reconstructed `m12_sq`: `7.49997525396243203e-02` GeV^2
- round-trip delta: `-2.47460023639067828e-07` GeV^2
- round-trip relative difference: `3.29946698186972485e-06`
- the round-trip reconstructed value is diagnostic only and must not be used as the downstream physical construction input.

## Soft-scale export gate

- `soft_scale_export_status`: `VALIDATED_BY_SET_PARAM_PHYS_REPLAY`
- replay `lambda1_target`: `1.00000000000000000e+00`
- replay `lambda1_reconstructed`: `1.00000093418204150e+00`
- replay `lambda1` absolute residual: `9.34182041500974947e-07`
- declared replay `lambda1` absolute residual tolerance: `9.99999999999999955e-07`
- maximum replay observable relative difference: `0.00000000000000000e+00` in `total_width_gev`
- replay gate fields: `total_width_gev`, `ctau_mm`, `br_bb`, `br_gammagamma`, `br_Zgamma`, and `lambda1_reconstructed`.

- center: `7.49999999996479594e-02` GeV^2
- one ULP: `1.38777878078144568e-17` GeV^2
- propagated half-ULP `lambda1` bound: `3.09033350026551164e-06`

- center reproduction tolerance: `1.00000000000000002e-08`
- maximum center reproduction difference: `0.00000000000000000e+00` in `total_width_gev`
- all values finite: `yes`
- all relevant values physically valid: `yes`
- adjacent-float maximum spread: `1.66564239607794041e-06` in `br_Zgamma`
- all adjacent probes theory-valid: `yes`

| Probe | ULP | theory_ok | lambda1 | total_width_gev | ctau_mm | br_bb | br_gammagamma | br_Zgamma |
|---|---:|:-:|---:|---:|---:|---:|---:|---:|
| previous_float | -1 | 1 | 1.00001173568050961e+00 | 4.56118529774213195e-14 | 4.32622153056751202e+00 | 7.56737485954531053e-01 | 2.65792771709345886e-04 | 2.04915281061213073e-05 |
| center | 0 | 1 | 1.00000093418204150e+00 | 4.56118529862185007e-14 | 4.32622152973311191e+00 | 7.56737485808578692e-01 | 2.65792941770435500e-04 | 2.04915508605281493e-05 |
| next_float | 1 | 1 | 9.99995533432807449e-01 | 4.56118529906170692e-14 | 4.32622152931591319e+00 | 7.56737485735602400e-01 | 2.65793026801000338e-04 | 2.04915622377361324e-05 |

| Field | Adjacent-float relative spread |
|---|---:|
| total_width_gev | 0.000000% |
| ctau_mm | 0.000000% |
| br_bb | 0.000000% |
| br_gammagamma | 0.000096% |
| br_Zgamma | 0.000167% |
| br_tautau | 0.000000% |
| br_gg | 0.000000% |

## Classification: `STABLE_AT_DOUBLE_REPRESENTATION_SCALE`

Maximum spread: `0.000167%` in `br_Zgamma`.

This result addresses double-representation uncertainty only; it does not validate scan provenance, production normalization, detector acceptance, or publication readiness.
