# CalcLam1ScanFromLake official-lake benchmark

## Configuration

- lake_root: `/mnt/c/Users/Asus/cern_db/dihiggs_lake`
- binary: `/home/fabi/wt_dihiggs_exploratory/2hdmc/CalcLam1ScanFromLake`
- mode: `triple_ok`
- n_campaigns: `3`
- n_files_per_campaign: `1`
- samples_per_file: `500`
- seed: `123456`
- campaign_regex: `<none>`
- require_header_match: `True`

## Selected official files

- `campaign=physdisc_zoom_tb6p3e6_l60p00082` — `/mnt/c/Users/Asus/cern_db/dihiggs_lake/campaign=physdisc_zoom_tb6p3e6_l60p00082/fixed_sinba=1p0000_l6=0p0008_l7=0p0000_mA=300p0/physdisc_zoom_tb6p3e6_l60p00082_iter0001_bin095/tb_6300000/scan_tb_6300000.csv`
- `campaign=physdisc_nb_tb6p3e6_l60p00082` — `/mnt/c/Users/Asus/cern_db/dihiggs_lake/campaign=physdisc_nb_tb6p3e6_l60p00082/fixed_sinba=1p0000_l6=0p0008_l7=0p0000_mA=300p0/physdisc_nb_tb6p3e6_l60p00082_iter0001_bin127/tb_6300000/scan_tb_6300000.csv`
- `campaign=physdisc_nb_tb4p4e6_l60p00115` — `/mnt/c/Users/Asus/cern_db/dihiggs_lake/campaign=physdisc_nb_tb4p4e6_l60p00115/fixed_sinba=1p0000_l6=0p0011_l7=0p0000_mA=300p0/physdisc_nb_tb4p4e6_l60p00115_iter0002_bin125/tb_4400000/scan_tb_4400000.csv`

## Results

| campaign | returncode | eligible_rows | actual_samples_written | baseline_failures | probe_failures | max_abs_error | max_rel_error | max_scaled_error | elapsed_sec |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| campaign=physdisc_zoom_tb6p3e6_l60p00082 | 6 |  |  |  |  |  |  |  | 0.115 |
| campaign=physdisc_nb_tb6p3e6_l60p00082 | 6 |  |  |  |  |  |  |  | 0.057 |
| campaign=physdisc_nb_tb4p4e6_l60p00115 | 6 |  |  |  |  |  |  |  | 0.357 |

## Aggregate summary

- campaigns tested: **3**
- files tested: **3**
- total eligible rows: **0**
- total actual samples written: **0**
- total baseline failures: **0**
- total probe failures: **0**
- runs with non-zero return code: **3**
- max(max_abs_error): **0**
- max(max_rel_error): **0**
- max(max_scaled_error): **0**

## Interpretation

One or more validator runs returned non-zero exit codes, so this benchmark should be treated as an execution audit rather than a physics-validation result. Verify the binary path and CLI contract before drawing robustness conclusions.
