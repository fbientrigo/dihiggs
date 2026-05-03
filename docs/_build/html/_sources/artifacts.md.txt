# Artifacts and Schemas

This page documents artifacts produced by scan orchestration and adaptive exploration.

## Orchestrator Run Artifacts

Under each run directory:

```text
<outdir>/<lake-name>/campaign=<campaign>/fixed_.../run_.../
```

Core files:

- `run_manifest.json`: run metadata and effective scan grid.
- `orchestrator.log`: timestamped execution log.
- `task_summary.jsonl`: per-task events (`done`, `skip`, `fail`, `crash`, `dry_run`).
- `tb_*/scan_tb_<tb_tag>.csv`: scan output for each `tanbeta`.
- `tb_*/scan_meta.json`: per-task metadata with `event`, `grid_signature`, and counters.

## Adaptive Iteration Ledger

`adaptive_explorer.py` writes one checkpoint directory per iteration:

```text
<checkpoint-root>/iter_XXXX/
```

Required files:

- `adaptive_state.json`
- `proposals.csv`
- `commands.sh`

### `adaptive_state.json` schema (high-level)

```text
{
  "iter_index": int,
  "n_iters": int,
  "lam1_min": float,
  "lam1_max": float,
  "n_bins": int,
  "n_coarse": int,
  "total_budget": int,
  "floor_points": int,
  "orchestrator": {...},
  "bins": [{"id","index","lam1_lo","lam1_hi"}],
  "budgets": {"lam1_bin_###": int},
  "budgets_by_tb_by_bin": {"<tb_tag>": {"lam1_bin_###": int}} | null,
  "trials_by_bin": {"lam1_bin_###": int},
  "successes_by_bin": {"lam1_bin_###": int},
  "trials_by_tb_by_bin": {"<tb_tag>": {"lam1_bin_###": int}},
  "successes_by_tb_by_bin": {"<tb_tag>": {"lam1_bin_###": int}},
  "proposals": [
    {
      "proposal_id": str,
      "status": "DONE",
      "command": str,
      "run_name": str,
      "run_dir": str | null,
      "n_lam1": int,
      "attempts_total": int?,
      "triple_ok_total": int?,
      "successes_by_tb": {"<tb_tag>": int}?,
      "trials_by_tb": {"<tb_tag>": int}?
    }
  ]
}
```

### `proposals.csv` schema

Header and status constraints are fixed:

```text
proposal_id,status,command
iter0001_bin000,DONE,"python ..."
```

- `status` is validated as `PENDING` or `DONE`.

### `commands.sh` schema

```bash
#!/usr/bin/env bash
set -euo pipefail
python .../orchestrate_scans.py ...
python .../orchestrate_scans.py ...
```

## Replay Behavior

- `python dihiggs/app/adaptive_explorer.py --replay <checkpoint_dir>`: executes command lines from `commands.sh` in order.
- `python dihiggs/app/adaptive_explorer.py --replay <checkpoint_dir> --list-commands`: prints `commands.sh` byte-identical.

## Strict Viability Semantics

- Success metric is strict viability (`triple_ok_points`) per bin and per `tb_tag`.
- For skipped tasks, recovery logic checks sibling `scan_meta.json` first (`triple_ok_points`).
- If metadata is unavailable, fallback is strict CSV recount with all three flags equal to 1.

## Iteration Semantics

- Iteration 1 (`iter_index == 1`) is coarse: each `lam1` bin uses `n_coarse`.
- Iteration 2 and later use adaptive budgets; when per-`tb_tag` data exists, explorer emits `--n-lam1-map` for orchestrator.
