# benchmark_search_daemon.py — operator README

Bounded, resumable background harvesting loop over the existing frozen LLP
discovery search substrate. It never reimplements physics, the ledger, the
archive, or the proposal strategies — it drives the existing
`search_substrate/` + `search_discovery/` code. See `IMPLEMENTATION_DECISION.md`
for the full architecture rationale and `RECONNAISSANCE_REPORT.md` for what
already existed before this daemon was added.

```
resume state -> propose -> dedup -> evaluate -> append evidence
  -> update archive -> checkpoint -> repeat
```

## What it does

Each cycle: loads/replays state, asks the configured proposal strategies
(`search_discovery/policy_adapter.py`) for a bounded batch of candidates,
deduplicates against everything already evaluated, evaluates the survivors
through the frozen evaluator (a worker pool, configurable size), appends every
result to the append-only ledger (successes *and* failures — failures are
data), updates the Top-5-per-family archive, checkpoints, prints a summary,
and repeats — bounded by `runtime`/`budgets`/`stopping` settings in the
config. It can run for hours or days unattended; no LLM is required for the
loop to operate.

## What it never does automatically

- Never writes into `runs/discovery_v1/`, `runs/discovery_photonic_v1/`, or
  anywhere under `deliverables/` — a path-allowlist guard
  (`search_discovery/daemon.py:guard_run_dir` / `guard_write_path`) refuses to
  start, or to write outside its own configured `run_dir`, including any
  subdirectory of a protected path. This is what keeps the daemon off the
  live, review-blocked campaign state.
- Never modifies `search_substrate/` or `search_discovery/`'s pre-existing
  files — `tests/test_frozen_substrate_unchanged.py` pins their hashes as a
  standing regression test.
- Never edits the frozen physics evaluator, theory gates, or Yukawa contract.
- Never promotes a candidate into a downstream benchmark snapshot. A
  discovered/archived candidate is evidence, not a promotion — see
  "Promotion protocol" below.
- Never rewrites or truncates a committed ledger row. The only two exceptions,
  both narrowly scoped and non-destructive to committed evidence: (1) on
  startup it truncates a torn *trailing* line left by a hard kill mid-append
  (never a completed record — `repair_torn_ledger_tail`), and (2)
  `rebuild_archive.py` resets the derived `discovery_archive.json` cache (never
  the append-only `discovery_archive_history.jsonl`) and only after the same
  path guard.
- Never runs a real physics evaluation in `--dry-run` mode.

## Start

```bash
cd /home/fabian/atlas_dihiggs/dihiggs
/home/fabian/atlas_dihiggs/.venv/bin/python benchmark_search_daemon.py --config configs/harvest_conservative.yaml
```

(There is no `dihiggs/.venv` — the one virtualenv for this workspace lives at
`/home/fabian/atlas_dihiggs/.venv`; every command below uses its absolute
path for that reason.)

For unattended background execution, wrap it with the repo's existing
resource-capped launcher rather than inventing a new one:

```bash
tmux new -A -s lab
./safe_run.sh --reserve 2 --mem-high-pct 70 --mem-max-pct 85 --name benchmark_harvest -- \
  /home/fabian/atlas_dihiggs/.venv/bin/python benchmark_search_daemon.py --config configs/harvest_conservative.yaml
```

## Resume

Safe to re-run identically after any stop (clean, killed, or crashed) — it
loads `runs/<run_dir>/daemon_checkpoint.json` for loop position (cycle count,
budgets consumed, policy/RNG state) and replays the ledger for candidate-level
dedup, so committed work is never lost or re-evaluated:

```bash
/home/fabian/atlas_dihiggs/.venv/bin/python benchmark_search_daemon.py --config configs/harvest_conservative.yaml --resume
```

## Dry run

Validates config, paths, schemas, envelope, and dedup logic, and prints the
planned cycle allocation — invokes the evaluator zero times
(`tests/test_daemon_process.py::test_dry_run_invokes_the_evaluator_zero_times`
proves this):

```bash
/home/fabian/atlas_dihiggs/.venv/bin/python benchmark_search_daemon.py --config configs/harvest_conservative.yaml --dry-run
```

## Stop

`Ctrl-C` (SIGINT) or `kill -TERM <pid>`. Both stop new candidate generation,
drain in-flight work up to `runtime.shutdown_grace_s`, write a final
checkpoint, and exit — SIGINT exits `130`, SIGTERM exits `143`. If a worker is
still stuck past the grace period, the daemon force-kills it rather than
hanging (see `search_discovery/daemon.py:_shutdown_pool` and
`_kill_process_group_except_self`) — verified by repeatedly killing a live
daemon subprocess mid-cycle and confirming a clean, bounded exit with no
leftover worker/`resource_tracker` processes.

The daemon isolates itself into its own process group at startup precisely so
this cleanup sweep can reliably reach every descendant (workers, and the
independent `multiprocessing.resource_tracker` helper, neither of which are
always fully accounted for by `ProcessPoolExecutor`'s own bookkeeping) without
depending on private internals. The one thing this cannot cover, by the
nature of Unix signals, is `kill -9`/SIGKILL sent directly to the daemon's own
PID — that bypasses our handler entirely, so no cleanup code runs and worker
processes could be orphaned. For unattended background operation, run under
`safe_run.sh`'s systemd `--user` scope (see "Start" above): a systemd scope
stop kills by cgroup membership, not process group, so it still reliably
cleans up the whole descendant tree even in that case.

## Inspect

```bash
# loop position / budgets consumed
cat runs/harvest_v1/daemon_checkpoint.json

# full diagnostics (valid fraction, dup rate, gate-failure breakdown, basins,
# best-by-family) via the pre-existing stats command, replayed from the ledger
/home/fabian/atlas_dihiggs/.venv/bin/python -m search_discovery checkpoint --run-dir runs/harvest_v1

# raw evidence
wc -l runs/harvest_v1/ledger.jsonl
tail -n1 runs/harvest_v1/discovery_archive.json
```

## Recover / rebuild derived state

The ledger is the only authoritative source; the archive is a derived cache.
To reconstruct it independently of the live `discovery_archive.json` (T4):

```bash
/home/fabian/atlas_dihiggs/.venv/bin/python -m search_discovery.rebuild_archive --run-dir runs/harvest_v1 \
  --target-dir /tmp/rebuilt_check
diff <(/home/fabian/atlas_dihiggs/.venv/bin/python -m json.tool runs/harvest_v1/discovery_archive.json) \
     <(/home/fabian/atlas_dihiggs/.venv/bin/python -m json.tool /tmp/rebuilt_check/discovery_archive.json)
```

Omit `--target-dir` to rebuild in place (still guarded against protected
run-dirs; never deletes `discovery_archive_history.jsonl`).

## Adding a new proposal policy

See the docstring at the top of `search_discovery/policy_adapter.py`. In
short: implement either `BatchStrategy.propose(context, budget, state) ->
ProposeResult` (the daemon core evaluates and calls `Ledger.append` for you,
then your `tell(proposals, results, state)`), or
`SelfContainedStrategy.run(context, budget, evaluator, ledger, state)` if your
strategy needs its own ask-evaluate-tell loop per step (like `climb`/
`continue_lineage`). Register it in `BATCH_REGISTRY` or
`SELF_CONTAINED_REGISTRY`, then reference it by name in a config's
`policy.allocation`. A strategy must never call `Ledger.append` itself in the
batch path, and must never mark a candidate valid — only the frozen evaluator
decides that.

## Promotion protocol (not automated)

A Top-5 archive entry becoming a new frozen downstream benchmark is a
separate, human-triggered decision, never performed by the daemon:

1. Candidate exists in the ledger with complete evidence and provenance.
2. Review the archive entry (`discovery_archive.json`) and its ledger
   evidence.
3. Explicit human/scientific promotion decision.
4. Write a new versioned artifact under `deliverables/` (e.g.
   `selected_candidates_v2.json`) — never overwrite an existing `*_v1`/`*_vN`
   file in place.
5. Downstream workflows are repointed explicitly at the new version.

## Tests

```bash
cd /home/fabian/atlas_dihiggs/dihiggs
/home/fabian/atlas_dihiggs/.venv/bin/python -m pytest tests/ -q
```
