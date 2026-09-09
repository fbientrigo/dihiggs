# Issue #82 reconnaissance

Scope: audit only. No search was launched and no evaluator, benchmark, or downstream snapshot was changed.

| Capability | Existing implementation | Evidence/path | Reusable? | Gap |
|---|---|---|---|---|
| Deterministic evaluator | Canonical evaluator and guarded orchestrator | `dihiggs/app/orchestrate_scans.py`, `dihiggs/app/orchestrator/runner.py` | Yes | Daemon must call this path only. |
| Bounded proposal execution | Adaptive and branch explorers plus bounded search | `dihiggs/app/adaptive_explorer.py`, `dihiggs/app/branch_continuation_explorer_v2.py`, `autoresearch/harness/bounded_adaptive_search.py` | Partly | Not a unified candidate-level daemon interface. |
| Proposal validation and lineage | Event-sourced proposal registry, contract validation, parent checks | `autoresearch/harness/proposal_registry.py` | Yes | Proposal ID is not a canonical physical-candidate ID. |
| Policy selection | UCB arm selection, adaptive/branch explorers, goal proposer | `autoresearch/harness/dihiggs_adaptation.py`, `autoresearch/harness/mvp_goal_proposer.py` | Yes | No minimal common `ask/tell` boundary separating policy from execution. |
| Resume/checkpoints | Grid-signature resume, per-task metadata, explorer checkpoints | `dihiggs/app/orchestrator/resume.py`, `dihiggs/app/adaptive_checkpoint.py` | Yes | Mutable state writes are not atomic; no daemon crash boundary. |
| Evidence ledger | JSONL attempt and proposal events | `autoresearch/harness/dihiggs_runner.py`, `autoresearch/harness/proposal_registry.py` | Partly | Attempt events have no single-writer lock; candidate and attempt identities are conflated by run context. |
| Reconciliation/dedup | SQLite upsert, identity contract, watermarks, graph deltas | `autoresearch/harness/mvp_upsert_pipeline.py`, `run_identity_contract.py`, `reconcile_watermark.py`, `mvp_graph.py` | Yes | Reconciliation is post-hoc; it cannot prevent concurrent physical duplicate evaluation. |
| Archive/frontier | Coverage state, graph, and frontier-aware goal proposal | `coverage_contract_bridge.py`, `mvp_graph.py`, `mvp_goal_proposer.py` | Partly | No documented rebuildable daemon archive/frontier derived solely from authoritative evidence. |
| Supervisor/runtime budgets | Round loop, timeout, max rounds/duration, preflight | `autoresearch/harness/campaign_supervisor.py`, `dihiggs_runner.py` | Partly | Serial runner; no worker pool, SIGINT/SIGTERM checkpointing, or cross-process lock. |
| Dry-run | Orchestrator and adapter dry-run paths | `dihiggs/app/orchestrate_scans.py`, `autoresearch/harness/orchestrator_adapter.py` | Yes | Daemon-level dry-run must guarantee zero evaluator invocations. |
| Failure preservation | Failure task metadata, run-health and quarantine paths | `dihiggs/app/orchestrator/runner.py`, `autoresearch/harness/run_health.py`, `autonomy_scheduler.py` | Yes | Need one normalized daemon event taxonomy. |
| Provenance/manifests | Campaign manifest, rerun manifest, external-tool lock | `campaign_manifest.py`, `rerun_manifest.py`, `external_tools.lock.yaml` | Yes | Need daemon state/evidence schema and documented ownership. |
| Frozen downstream snapshots | Frozen benchmark artifacts are checked in | `benchmarks/FIRST_H2_RECAST_CANDIDATE.json`, `docs/pilots/` | Yes | No regression proving automatic discovery never writes them. |
| Conservative operations config | Smoke and explorer configs | `autoresearch/configs/smoke_test_minimal.json`, `autoresearch/configs/dihiggs_explorers*.json` | Partly | Existing smoke config is not a checked-in, no-LLM, conservative harvesting config. |
| LLM controller | No mandatory runtime dependency found | N/A | Yes—omit | Optional controller is not needed for the first implementation. |

## Evidence from the audit

- `CampaignSupervisor` already enforces preflight, round/time bounds and persists a campaign state, but writes that state directly and executes rounds serially.
- `proposal_registry.append_event` fsyncs an append-only JSONL record, while `DiHiggsRunner.emit_attempt_event` uses an unlocked append. Neither is a cross-process single-writer commit path.
- The SQLite upsert pipeline performs transactional derived-state updates and its identity contract is useful for reconciliation, but its attempt identity includes a run fingerprint. It cannot reserve a canonical candidate before evaluation.
- Existing test coverage exercises resume, stale-output rejection, proposal replay, deduplication during reconciliation, dry-run, checkpoints, bounded search and supervisor behavior. The focused audit suite passed: `97 passed, 43 skipped`.
- No existing implementation handles `SIGINT`/`SIGTERM` for this lifecycle, persists an atomic daemon checkpoint, or proves a concurrent worker pool cannot duplicate an evaluation or corrupt the ledger.
- No existing test protects `table_snapshot_v1` or the current frozen downstream artifacts from automated discovery writes.

## MINIMAL_IMPLEMENTATION_PLAN

1. Add one small harvesting entrypoint which composes the existing proposal validation and canonical orchestrator; do not add an evaluator path, distributed service, LLM runtime, or optimizer framework.
2. Define a versioned JSONL evidence record with separate deterministic `candidate_id` and execution `attempt_id`; reserve/deduplicate candidates under one local lock before dispatch.
3. Use a small local worker pool with bounded workers, evaluation timeout, total/cycle budgets, and a single evidence writer. Persist its mutable checkpoint atomically with `os.replace`.
4. Rebuild the daemon archive/frontier from evidence using the existing coverage/graph helpers; retain SQLite only as a derived cache where it already helps.
5. Reuse the existing dry-run adapter and add focused tests for kill/resume, concurrent dedup, append-only evidence, archive rebuild, malformed proposals, zero-evaluator dry-run, and frozen-snapshot non-mutation.
6. Add a conservative checked-in config and concise operations README; run only the required bounded `INFRASTRUCTURE_SMOKE_ONLY` exercise after implementation review.

## Lead review gate

The required reconnaissance is complete. The plan intentionally excludes a new physics campaign, changes to the canonical evaluator, snapshot promotion, a dashboard, distributed infrastructure, and an LLM controller. Implementation must not begin until this report is reviewed.
