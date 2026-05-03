# HPC Operations and Reproducibility

This page captures practical guidance for running large scans on shared filesystems and reproducing adaptive runs.

## Threads and OMP bookkeeping

- Use `--threads N` with `orchestrate_scans.py` to set `OMP_NUM_THREADS=N` for scan subprocesses.
- If `--threads` is omitted, `OMP_NUM_THREADS` is not forced by the orchestrator.
- Thread settings are recorded in `run_manifest.json` under `runtime.omp_num_threads`.
- Adaptive checkpoints also record thread context in `adaptive_state.json` (`omp_num_threads`, plus orchestrator `threads`), and may include `omp_num_threads_manifest` per proposal when a run manifest is available.

## Filesystem behavior on HPC

- Runs generate many small files (`task_summary.jsonl`, per-`tb_*` CSVs and `scan_meta.json`, plus logs/manifests).
- On networked storage, metadata-heavy workloads can dominate runtime; avoid over-fragmenting work.
- Keep `--n-bins` modest in adaptive mode to reduce per-iteration directory and file churn.
- Prefer output locations optimized for many file creates/updates; avoid slow sync mounts for active writes when possible.

## Resumability and reuse semantics

- `--resume-scope run`: reuse only artifacts already present in the current `run_dir`.
- `--resume-scope fixed`: may reuse previous results under the same fixed-parameter directory when grid signatures match.
- Cross-run reuse requires matching `grid_signature` (`mphi` and `lam1` ranges with point counts, including effective `n_lam1`).
- `--materialize` copies reused previous-run CSVs into the current run directory; without it, skip records point to prior files but do not duplicate data locally.

### `run_*` reuse constraint

Cross-run discovery scans `run_*/tb_*/scan_tb_<tb_tag>.csv` under the fixed directory. Reuse from previous runs therefore depends on run directory names matching `run_*`.

## Reproducibility checklist

For each adaptive iteration (`iter_XXXX`), preserve all three checkpoint artifacts:

- `commands.sh`: canonical replay command sequence (`--replay` executes these lines).
- `proposals.csv`: proposal IDs, statuses, and exact command strings.
- `adaptive_state.json`: budgets, bin definitions, per-bin and per-`tb` outcomes, run mapping, and orchestrator context.

Also retain each referenced orchestrator `run_manifest.json` and `task_summary.jsonl` for full auditability.

## Operational footguns

- Passing custom `--run-name` values that do not start with `run_` can break future cross-run reuse discovery.
- Using `--resume-scope run` when you expected historical reuse can cause unnecessary recomputation.
- Forgetting `--materialize` can surprise downstream tools that expect reused CSVs physically present in the current `run_dir`.
- Large `--n-bins` with many `tanbeta` values can overload filesystem metadata operations before CPU becomes the bottleneck.
