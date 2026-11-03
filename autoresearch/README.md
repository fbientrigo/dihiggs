# Autoresearch DiHiggs Harness

## Current path: budgeted loop over dihiggs_boundary (start here)

`python -m autoresearch.run_supervisor <config.json> --cpu-hours H --disk-gb G --wall-hours W` runs a campaign. Each round
the policy chain (UCB1 -> abstain -> `agy` -> rules, every hop journaled) picks one key from the registered action menu
(`harness/boundary_actions.py`: strategy x tan_beta slice x regime, lambda7 != 0). The action is
`lake explored --export` -> `dhb.<strategy> --explored` -> `lake commit`; boundary computes every number, autoresearch only
chooses. `<outdir>/journal.db` is an append-only SQLite journal (write-ahead intents, integer budget; kill -9 reconciles
and never relaunches). The budget is never exceeded: the next action's upper bound is reserved before it starts. Only
`survivor_v2`, `mh` and `pso` actions are runnable (`adaptive_m2_mesh`, `profiled_m2` cannot read `--explored` yet and are
masked). survivor_v2 is a fixed mesh: once the lake holds its quota a repeat brings 0 rows. `mh` (Metropolis-Hastings chains
on eps + 1[theory_ok]) and `pso` (swarm on theory_ok x distance to the nearest evaluated point) from `dhb.stochastic_search`
evaluate 2048 fresh continuous points per attempt, seeded from the attempt id. An arm whose last 2 attempts brought 0 fresh
rows is masked; when every arm is, the campaign ends `SATURATED` (hop `stop:saturated`, no agy call). After `lake verify`
confirms an action's commit, its `out/` data is deleted locally (manifests, logs and boundary config stay; disk is charged
before). Config: `configs/boundary_lambda7_search.json` (24 arms); smoke: `configs/boundary_lambda7_smoke.json`. `python -m autoresearch.tui` shows budget remaining, per-action
runs and coverage (Overview) from the same `journal.db`.

**Smoke (`make autoresearch-smoke`, ~40 s, status TESTED — not VALIDATED).** Unit tests, then a real run with a throwaway
local-dir lake (`LAKE_REMOTE`/`LAKE_CACHE` in a temp dir; the real lake is untouched): compiles the evaluator, runs
boundary (lambda7 = 0.01, tan_beta = 1e7, 32 points), commits to the lake; a second campaign over the same lake computes
0 fresh rows (32 seeded) and adds no commit; a daemon `kill -9`'d mid-action is restarted and ends with the attempt closed,
32 rows and exactly 1 commit; finally the TUI is mounted headless (Textual pilot) and its Overview text is printed.
Verified 2026-11-02 from a fresh `git clone --recursive` of atlas_dihiggs with no env vars: workspace `make setup && make smoke`
on Debian/Python 3.13.5 and in CI ([Fresh-clone smoke](https://github.com/fbientrigo/atlas_dihiggs/actions/workflows/smoke.yml),
ubuntu/Python 3.13, run 41531514883): 215 unit tests pass, all smoke checks `ok`, `autoresearch-smoke TESTED`.

Trimmed in #18: the pre-journal modules nothing reached (bounded_adaptive_search, campaign_state, run_health,
safe_automation_layer, mvp_graph, proposal_registry, ...; net -14.6k lines with their tests) and the eight loose guides were removed; they remain in
git history.

The older explorer-based harness (`cmd` arms wrapping `adaptive_explorer.py` / `branch_continuation_explorer_v2.py`,
checkpoint parsers, axis contract, `dihiggs_explorers*.json` configs) was removed with the legacy lambda1 subsystem. Only
`action` arms run today. Reproduce the old harness with `git checkout legacy-lambda1-final`.
