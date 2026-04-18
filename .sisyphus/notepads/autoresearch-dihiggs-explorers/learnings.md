# Learnings — autoresearch-dihiggs-explorers

Conventions, patterns, and wisdom accumulated during implementation.

---

- Task 1 contract anchored in `autoresearch/harness/dihiggs_axis_contract.py`: mode name constant, canonical `cell_id`, and required axes live there.
- Canonical `cell_id` stays backward-compatible by always starting `tb=<...>|bin=<lam1_bin>` and sorting only extra keys for determinism.
- Missing axes semantics already fit sparse aggregation needs in `autoresearch/benchmarks/score.py`: absent axis values are skipped, not converted into synthetic `NA` bins.
- `autoresearch/configs/dihiggs_explorers.json` must follow the authoritative Appendix example from the plan, not the older README-era sample.
- Adaptive checkpoint parsing should separate full proposal inventory from evaluable discoveries: `parse_adaptive_checkpoint_inventory(...)` keeps DONE/failed proposal metadata, while `parse_adaptive_checkpoint(...)` filters to stable, evaluable `run_dir` entries for current harness consumers.
- Robust adaptive iter discovery should accept numeric `iter_*` directories generically and sort by parsed integer suffix instead of hard-coding four digits.
- Branch continuation parsing needs dual JSONL modes: bare artifact rows with direct `step_label/result` fields, and enveloped harness rows where only `event_type == "ATTEMPT_COMPLETED"` should feed inventory/discovery parsing.
- Branch `lam1_raw` should be stored explicitly on parsed discoveries so missing `selected` blocks do not silently collapse to `0.0`; legacy `params.lam1` remains valid fallback metadata.
- Branch continuation parsing now follows same split: `parse_branch_checkpoint_inventory(...)` keeps per-line metadata from `events.jsonl`, while `parse_branch_checkpoint(...)` filters to evaluable discoveries with stable `run_dir`, `tanbeta`, and selected `lam1` for runner compatibility.
- Branch event parser should accept both current raw attempt JSONL rows and older `payload`/`params`-wrapped rows so existing branch parser callers/tests keep working during migration.
- Task 4 preflight now has a dedicated PhysScanWithFixings gate in `autoresearch/harness/dihiggs_preflight.py`: required dataset envs are enforced only when the selected config actually needs them, `libgsl` is required, and missing `gfortran` is reported as optional/skip so mocked unit environments can still pass.
- `autoresearch/harness/__main__.py` handles `--mode dihiggs-explorers --preflight-only` before any scheduler launch side effects, while `scheduler.run_campaign(..., preflight_only=True)` returns a minimal summary for unit coverage and future mode wiring.
- Task 4: `run_physscanwithfixings_preflight(...)` is stricter than legacy `run_all_preflight_checks(...)`: it hard-fails on missing `libgsl`, `HB_DATASET`, `HS_DATASET`, or non-executable `PhysScanWithFixings`, and embeds exact README snippets for remediation.
- Task 5 MVP runner wiring can safely reuse `DiHiggsRunner` from `scheduler.run_campaign(...)` as long as scheduler owns the strict preflight first and passes it through; this keeps stable `attempt_id`, per-slice adaptive emission, parser reuse, and resume-safe dedupe without importing `dihiggs/app/*`.
- Duplicate artifact policy needs two layers: identical `campaign_id + arm_id + fingerprint + slice_key` suppresses re-emission entirely, while same artifact fingerprint on a new slice emits `ATTEMPT_EVALUATED` with `eval_status=SKIP_REUSED`, `successes=0`, `trials=0`, and `is_duplicate=true`.
- The effective dihiggs mode cap should use the tighter of `limits.max_new_run_dirs_per_round` and `limits.max_new_run_dirs_per_arm_call`; fixture tests should set both high enough when exercising multi-discovery scenarios.
- Task 5 runner MVP uses checkpoint inventory APIs, not the filtered `parse_*_checkpoint(...)` helpers, because missing `run_dir` rows must still become failure events (`CRASH`) instead of disappearing during parsing.
- Adaptive mode emission is per tan(beta) slice: one evaluated attempt per `successes_by_tb`/`trials_by_tb` key, all sharing the same artifact fingerprint when they come from the same physical `run_dir`.
- Duplicate handling needs two layers: stable `attempt_id = sha256(campaign_id + arm_id + fingerprint + slice_key)` for idempotent replay, plus `is_duplicate=true` when a new slice references an already-seen artifact fingerprint from an earlier round.
- Task 6 score path now defaults enabled multi-axis metrics to physical axes `tb`, `lam1_bin`, `mphi_bin` and default pairs `tb|lam1_bin`, `lam1_bin|mphi_bin` when config omits explicit specs; metrics still prefer `payload.axes_binned` and fall back to `cell_id` parsing.
- Task 6 collapse pressure in `compute_metrics(...)` must use `metrics.multi_axis.collapse_axes` (default `tb+lam1_bin`) rather than raw `cell_id`, so different `mphi_bin` values can still count as the same collapsed slice.
- Task 6 scoring needs `compute_metrics` overload preservation: legacy 3-arg run-dir CSV extraction still serves `autoresearch/tests/test_score.py`, while 2-arg event-log scoring serves benchmark campaigns and `run_benchmarks.py`.
- Task 6 multi-axis normalization belongs in the events-derived scoring path only; keeping `compute_diversity(...)` on raw entropy preserves older helper tests like `test_multiaxis_metrics_missing_axes` while `compute_metrics(events_path, config)` applies entropy/log(domain) per plan.
- Multi-axis fallback parsing should trust `payload.axes_binned` first, then decode `cell_id` segments like `bin=<n>` -> `lam1_bin` so mixed old/new fixtures score identically.
- Task 7 arm selection now lives in `DiHiggsRunner`: if `adaptation.warm_start_each_arm` is true, arms are pulled in config order until each has `n>0`; after that, selection uses `mean_reward + ucb_c * sqrt(log(total_rounds) / n)` from persisted state.
- Task 7 tests stay deterministic by patching `run_single_round()` and campaign score snapshots directly, which isolates UCB math and mean-reward persistence without requiring real explorer subprocesses or checkpoint fixtures.
- DiHiggs fixture tests that write real runner state must not hardcode `campaign_id`, because bandit state persists outside pytest/unittest tempdirs in `autoresearch/opencode_logs/<campaign_id>/state/dihiggs_mode_state.json`; UUID-suffixed campaign IDs eliminate cross-test arm-selection drift.
- Scheduler entrypoint for `dihiggs-explorers` must call `DiHiggsRunner.run(...)`, not `run_single_round(...)`, or `runtime.max_rounds` and persisted UCB updates never take effect from CLI/scheduler flows.
- Negative composite deltas should be tested explicitly: patching `_compute_campaign_score` from `0.9 -> 0.3` proves exact normalization `reward = (delta + 1) / 2 = 0.2` and catches accidental sign/clamp regressions.
