# `orchestrate_scans.py`

CLI orchestrator for running 2HDM scans and persisting outputs in a resumable run layout.

## Usage

```bash
python dihiggs/app/orchestrate_scans.py [OPTIONS]
```

## Key Flags

- Runtime:
  - `--exec`: path to C++ executable (default `./PhysScanWithFixings`)
  - `--threads`: sets `OMP_NUM_THREADS`
  - `--dry-run`: create folders/manifest and print commands without executing
  - `--force`: rerun even when output CSV already exists
- Output layout:
  - `--outdir`: output root
  - `--lake-name`: subfolder under `--outdir` (default `dihiggs_lake`)
  - `--campaign`: campaign label
  - `--run-name`: optional explicit run folder name
- Resumability:
  - `--resume-scope run|fixed`:
    - `run`: only checks checkpoints in current `run_...` directory
    - `fixed`: also reuses checkpoints from previous runs under same fixed-parameter directory
  - `--materialize`: when reusing previous-run CSV, copy it into current run folder
- Fixed parameters:
  - `--mA`, `--sin-ba`, `--lambda6`, `--lambda7`
- Grid parameters:
  - `--mphi-min`, `--mphi-max`, `--n-mphi`
  - `--lam1-min`, `--lam1-max`, `--n-lam1` (also accepts legacy `--N-lam1`)
  - `--n-lam1-map`: per-`tanbeta` lambda1 point override
- Scan list:
  - `--tanbeta`: comma-separated values (for example `10000,15000,20000`)

## `--n-lam1-map` Format

Use comma-separated `tanbeta:n_lam1` pairs:

```text
--n-lam1-map "10000:200,15000:400"
```

- Keys are `tanbeta` values (float)
- Values are `n_lam1` (positive integer)
- Any `tanbeta` not listed falls back to `--n-lam1`

## Outputs

The orchestrator writes to:

```text
<outdir>/<lake-name>/campaign=<campaign>/fixed_sinba=..._l6=..._l7=..._mA=.../run_.../
```

Typical files in each run:

- `run_manifest.json`: run metadata, fixed parameters, grid, runtime settings, summary
- `orchestrator.log`: timestamped orchestrator logs
- `task_summary.jsonl`: one JSON record per tan(beta) task event
- `tb_.../scan_tb_<tag>.csv`: scan outputs
- `tb_.../scan_meta.json`: per-task checkpoint metadata (grid signature, status)

## Resuming and Cross-Run Reuse

When `--resume-scope fixed` is used, previous CSVs are reused only if the effective grid signature matches.

Cross-run lookup is constrained to paths matching:

```text
run_*/tb_*/scan_tb_<tb_tag>.csv
```

This `run_*` constraint is required by the current implementation (`find_previous_csv_matching_grid`).

If `--materialize` is set and a previous matching CSV is found, the file is copied into the current run folder; otherwise it is only referenced for skip logic.
