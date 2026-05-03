# Decisions — autoresearch-dihiggs-explorers

Architectural choices and rationale recorded during implementation.

---

- Defined the Task 1 mode contract as `dihiggs-explorers` in code and config, but left actual execution wiring unimplemented on purpose because runner work belongs to later plan tasks.
- Kept canonical `cell_id` logic separate via `canonical_cell_id()` and made `encode_cell_id()` delegate to it, preserving old call sites while making determinism explicit.
- Added a minimal README schema note because the broader README still documented an older config shape; the small snippet removes Task 1 ambiguity without pulling in Task 2+ details.
# Decisions — autoresearch-dihiggs-explorers

- Task 2: kept `parse_adaptive_checkpoint(...)` as backward-compatible wrapper for current harness call sites, but added `parse_adaptive_checkpoint_inventory(...)` so Task 5 can consume failure metadata without reparsing or changing black-box explorer boundaries.
- Task 2: `AdaptiveDiscovery.tb` and `AdaptiveDiscovery.lam1_raw` remain compatibility aliases (`bin_index` and `lam1_min`) until runner-side event emission is updated in a later task.
- Task 3: kept `parse_branch_checkpoint(...)` as backward-compatible wrapper over `parse_branch_checkpoint_inventory(...)`, but made branch rows carry explicit `lam1_raw_value` so callers can bin legacy `params.lam1` even when `result.selected` is absent.
- Task 3: when branch JSONL has a harness-style event envelope, parser accepts only `ATTEMPT_COMPLETED`; direct branch artifact rows without `event_type` still parse for black-box explorer compatibility.
- Task 4: kept legacy preflight helpers intact for older tests/callers, but switched `DiHiggsRunner` and `dihiggs-explorers --preflight-only` to the new PhysScan-focused preflight path so end-to-end runtime gating gets stricter without rewriting unrelated harness behavior.
- Task 5: rewrote `autoresearch/harness/dihiggs_runner.py` around first-arm execution only, matching the MVP scope exactly; bandit update/reward logic stays out of this runner until later plan work.
- Task 5: runner emits state-compatible `ATTEMPT_CREATED`/`ATTEMPT_STARTED`/`ATTEMPT_EVALUATED`/terminal events so existing `state.py` derivation and reporting tooling can consume DIHiggs mode attempts without special-case schema handling.
- Task 6: restored `autoresearch/benchmarks/score.py` as a dual-contract module instead of forcing one API shape; benchmark event scoring and older run-dir CSV tests stay alive simultaneously.
- Task 6: duplicate events count toward reuse via either `is_duplicate=true` or legacy `eval_status="SKIP_REUSED"`, but are excluded from yield, coverage, diversity, and collapse aggregation in the new multi-axis path.
- Task 7: kept DiHiggs bandit stats in `autoresearch/opencode_logs/<campaign_id>/state/dihiggs_mode_state.json` rather than generic scheduler state, because this task only needs per-arm `n`/`mean_reward` persistence and should not pre-empt Task 8 resume semantics.
- Task 7: arm reward is campaign-level composite score delta from `compute_metrics(events_path, config)` and `compute_composite(metrics)["score"]`, remapped from `[-1, 1]` into `[0, 1]` so bandit credit reflects net campaign improvement.
- Test isolation fix: `_base_config(...)` in the DiHiggs mode fixture tests now generates a unique `campaign_id` per call, so persisted `dihiggs_mode_state.json` under `autoresearch/opencode_logs/` cannot leak warm-start arm counts across test cases.
