# CalcLam1ScanFromLake official-lake benchmark

## Configuration

- lake_root: `/mnt/c/Users/Asus/cern_db/dihiggs_lake`
- binary: `/home/fabi/wt_dihiggs_exploratory/2hdmc/CalcLam1ScanFromLake`
- mode: `triple_ok`
- n_campaigns: `3`
- n_files_per_campaign: `1`
- samples_per_file: `500`
- seed: `123456`
- campaign_regex: `scan_l6_tb_xhigh|sinba_0p9999_tb_hi`
- require_header_match: `True`

## Selected official files

- `campaign=scan_l6_tb_xhigh` — `/mnt/c/Users/Asus/cern_db/dihiggs_lake/campaign=scan_l6_tb_xhigh/fixed_sinba=1p0000_l6=0p0010_l7=0p0000_mA=300p0/branchcurve_l6_0p00100000_tb_10000p00000000_seed_try_01_global_seed/tb_10000/scan_tb_10000.csv`
- `campaign=sinba_0p9999_tb_hi` — `/mnt/c/Users/Asus/cern_db/dihiggs_lake/campaign=sinba_0p9999_tb_hi/fixed_sinba=0p9999_l6=0p1000_l7=0p0000_mA=300p0/run_20260126T214835Z_host=zephyrus_git=a718276/tb_20000/scan_tb_20000.csv`

## Results

| campaign | returncode | eligible_rows | actual_samples_written | baseline_failures | probe_failures | max_abs_error | max_rel_error | max_scaled_error | elapsed_sec |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| campaign=scan_l6_tb_xhigh | 0 | 13340 | 500 | 0 | 0 | 0 | 0 | 0 | 0.564 |
| campaign=sinba_0p9999_tb_hi | 0 | 17444 | 500 | 0 | 0 | 0 | 0 | 0 | 16.008 |

## Aggregate summary

- campaigns tested: **2**
- files tested: **2**
- total eligible rows: **30784**
- total actual samples written: **1000**
- total baseline failures: **0**
- total probe failures: **0**
- runs with non-zero return code: **0**
- max(max_abs_error): **0**
- max(max_rel_error): **0**
- max(max_scaled_error): **0**

## Interpretation

In this tested sample from multiple official campaigns, no baseline/probe failures were reported. These results support consistency for the exercised files and sample counts, without claiming full parameter-space coverage.
