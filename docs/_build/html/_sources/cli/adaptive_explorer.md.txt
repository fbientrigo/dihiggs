# `adaptive_explorer.py`

Adaptive CLI loop over `orchestrate_scans.py` that allocates `lam1` effort by observed viability.

## Usage

```bash
python dihiggs/app/adaptive_explorer.py [OPTIONS]
```

## Adaptive Loop

- Iteration 1 is coarse exploration: every `lam1` bin gets `--n-coarse` points.
- Iteration 2 and later switch to adaptive allocation from prior outcomes.
- If per-`tanbeta` history is available, the explorer uses per-`tb_tag` budgets and passes them to the orchestrator through `--n-lam1-map`.
- In adaptive iterations, the wrapper sends `--n-lam1 <floor_points>` and lets `--n-lam1-map` carry the per-`tb_tag` effective counts.

## Objective and Budgeting

- The objective is strict viability: maximize expected strict-valid points (`triple_ok_points`) rather than raw attempts.
- Adaptive scoring is Beta-Binomial over successes and trials per bin (and per `tb_tag` when available).
- Per-`tanbeta` budgeting uses `total_budget` as budget-per-`tb_tag` and enforces `floor_points` per bin.

## Skip Recovery and Metrics

- Adaptive accounting reads `task_summary.jsonl` from each produced run directory.
- For `skip` events, strict-valid recovery first checks sibling `scan_meta.json` (`triple_ok_points`).
- If metadata is missing, recovery falls back to counting strict-valid rows in the CSV (`positivity_ok=1`, `unitarity_ok=1`, `perturbativity_ok=1`).
- This ordering avoids undercounting when skips happen across resumed runs.

## Checkpoints and Replay

- `--checkpoint-root` is required; each iteration writes `iter_XXXX/`.
- Each iteration ledger includes:
  - `adaptive_state.json`: full iteration state, budgets, outcomes, and proposals
  - `proposals.csv`: proposal id/status/command ledger
  - `commands.sh`: replayable orchestrator commands
- Replay modes:
  - `--replay <iter_dir>` executes commands from `commands.sh`
  - `--replay <iter_dir> --list-commands` prints `commands.sh` byte-identical to stdout

## Smoke Example

```bash
python dihiggs/app/adaptive_explorer.py \
  --checkpoint-root /tmp/adaptive_smoke \
  --n-iters 2 \
  --n-bins 3 \
  --n-coarse 4 \
  --total-budget 12 \
  --floor-points 1 \
  --run-name-prefix run_smoke \
  --exec ./PhysScanWithFixings \
  --outdir /mnt/c/Users/Asus/cern_db \
  --orchestrator-dry-run
```
