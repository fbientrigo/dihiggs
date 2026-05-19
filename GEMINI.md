# Gemini Operator Guide (DiHiggs Autoresearch)

## 1) Purpose
Run bounded, evidence-gated DiHiggs/2HDMC sensitivity probes to improve ctau_m safely and incrementally.

## 2) Current architecture
- Backend execution: `dihiggs/app/orchestrate_scans.py`
- Safety gates: `autoresearch/harness/safe_automation_layer.py` (v1.1 fail-closed)
- Adaptive runner: `autoresearch/harness/bounded_adaptive_search.py`
- Decision helpers: `autoresearch/harness/decision_layer.py`
- Worker loop wrapper: `scripts/gemini_worker_loop.py`
- Progress artifacts: `progress.jsonl`, `progress.md`, `attempts.jsonl`, `commands.jsonl`

## 2.1) Roles
Gemini Worker:
- read latest `progress.md`/`progress.jsonl` and summaries
- propose or run bounded contracts only
- execute only approved harness paths via worker loop
- summarize each iteration and stop on configured conditions

Hermes Supervisor:
- audit worker outputs every configured cadence
- review failures/drift/invalid JSON
- validate any code changes and repair bugs
- approve envelope widening or semantic changes

Autoresearch runtime:
- enforces envelope bounds
- enforces Safe Automation gates
- enforces triple_ok-only physics interpretation
- writes stop reports and blocks invalid evidence

## 3) Safe operating rules
- Never run broad scans.
- Never enable unrestricted proposal generation.
- Never edit C++ or physics/scoring/triple_ok semantics.
- Keep within declared envelopes.
- Exclude failed attempts from physics conclusions.
- Use triple_ok-only rows for physics claims.

## 4) Current best state
- Best ctau_m: `5.390644846600711e-4 m`
- Best point approx: `mphi=200, mA=500, tan_beta=126904, lambda6~0.0019, lambda1=1.0`
- Observed leverage: high mA + low mphi improved ctau.

## 5) Run bounded adaptive searches
Example (plan):
`source ~/higgs_env_py312/bin/activate && python autoresearch/harness/bounded_adaptive_search.py --contract <contract.json> --output-dir scripts/out/<run> --plan-only`

Example (execute bounded contract):
`source ~/higgs_env_py312/bin/activate && python autoresearch/harness/bounded_adaptive_search.py --contract <contract.json> --output-dir scripts/out/<run> --execute`

## 6) Inspect outputs quickly
- `cat scripts/out/<run>/progress.md`
- `tail -n 20 scripts/out/<run>/progress.jsonl`
- `cat scripts/out/<run>/adaptive_state.md`
- `cat scripts/out/<run>/decisions.md`

## 7) What NOT to do
- No autonomous exploration outside envelope.
- No all-row metrics as physics claims.
- No silent gate bypass.
- No edits to C++ or triple_ok semantics.

## 8) Reporting format (every subcampaign)
1. Main result.
2. Gates status.
3. triple_ok rate.
4. Best ctau_m + best point.
5. Budget remaining.
6. Recommendation + reason.

## 9) Near-term objective
Within declared envelope: improve best ctau_m by **>=2x** over `5.390644846600711e-4 m`.
Stop immediately on gate failures that require scope expansion.

## 10) Exploration policy for Gemini templates
- Do not repeat the same exploitation template more than 2 times without improvement.
- If improvement over the last 2 iterations is <10%, run `sensitivity_probe`.
- Preferred sensitivity variable order:
  1) `tan_beta`
  2) `lambda6`
  3) `mA`
  4) `lambda1` only when explicitly enabled
- Always keep `lambda7=0.0` and `sin_ba=1.0`.

## 11) Active supervised cycle contract (cycle01)
- Execute exactly 3 iterations unless a gate fails.
- Iteration 1: `run_template` with `template_name=sensitivity_probe`, `variable=tan_beta`, `relative_steps=[-0.10,0.0,0.10]`.
- Iteration 2: `run_template` with `template_name=sensitivity_probe`, `variable=lambda6`, `relative_steps=[-0.10,0.0,0.10]`.
- Iteration 3 decision:
  - if either prior probe best_gain_fraction >= 0.05, exploit best_direction using bounded local contract;
  - else run `sensitivity_probe` on `mA` with same relative steps.
- Keep fixed controls: `lambda7=0.0`, `sin_ba=1.0`, `lambda1=1.0`, `n_lam1=1`.
- Anchor near current best low-mphi region unless progress evidence says otherwise.
- Stop immediately if any probe fails validation or contamination guard triggers.