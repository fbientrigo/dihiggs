# Autoresearch: run DiHiggs adaptive/branch explorers (black-box) with multi-axis diversity + composite scoring

## TL;DR
> **Summary**: Extend `python -m autoresearch.harness` with a new mode that runs the DiHiggs explorers as subprocesses, ingests their checkpoint outputs to extract `run_dir`s, evaluates artifacts into harness events, and adapts which explorer+hparams to run next to maximize the composite score while enforcing multi-axis exploration diversity.
> **Deliverables**:
> - New harness mode (CLI) under `autoresearch` that executes `dihiggs/app/adaptive_explorer.py` and `dihiggs/app/branch_continuation_explorer_v2.py` as black boxes
> - Robust parsers for adaptive/branch checkpoints → `run_dir` inventory
> - Multi-axis coverage/diversity support in `autoresearch/benchmarks/score.py` (backwards compatible)
> - Real end-to-end runnable path: build `dihiggs/app/PhysScanWithFixings`, run a minimal campaign, produce reports
> - Unit tests (unittest) + gated real E2E smoke
> **Effort**: XL
> **Parallel**: YES — 4 waves
> **Critical Path**: Multi-axis metrics contract → checkpoint parsers → new harness mode runner → real build preflight + E2E smoke

## Context
### Original Request
- Provide an autoresearch/auto-learning system that can run and adapt the DiHiggs search.
- Integrate `dihiggs/app/` adaptive explorer and branch adaptive explorer (branch continuation) which run with hyperparameters.
- Autoresearch must execute those commands and adapt strategy to get more “physics points” while keeping exploration variety.

### Interview Summary (decisions locked)
- Deliverable: **Python CLI** (`python -m autoresearch.harness ...`).
- Objective: **composite score** from `autoresearch/benchmarks/score.py` (`compute_metrics` + `compute_composite`).
- Integration: **subprocess black-box** (execute DiHiggs explorer CLIs; parse outputs/checkpoints).
- Variety: **multi-axis** (not just lam1; include additional axes such as tanbeta/mphi/track/strategy).
- Must run **real end-to-end** (compile and use `dihiggs/app/PhysScanWithFixings`).

### Metis Review (gaps addressed)
- Metis agent was not available (tool timeouts/errors). Mitigation applied:
  - Added explicit, decision-complete contracts for: event schema, checkpoint parsing, axis binning, and failure mapping.
  - Added strong guardrails + extensive verification tasks (unit + gated real E2E) to catch brittle integration.

### Oracle Review (architecture decisions)
- Treat DiHiggs explorers as **strictly opaque** producers of checkpoints; harness does only: choose arm → run subprocess (timeouts, logs) → discover `run_dir`s → evaluate → emit **idempotent** events.
- Multi-axis diversity must be **sparse**:
  - Track coverage/diversity on **per-axis bins** plus a small fixed set of **pairwise** interactions (avoid Cartesian `total_cells`).
- Idempotence/resume is a first-class requirement:
  - Persist processed fingerprints + emitted attempt_ids and never re-emit.
  - Cap discoveries per subprocess call.

### Repo Constraint: required harness interpreter
- `autoresearch/harness/__main__.py` calls `require_higgs_env()` which enforces running the harness with `~/higgs_env_py312/bin/python` (see `autoresearch/harness/python_env.py`).
- All new mode behavior must remain compatible with this constraint (do not remove the guard; update docs and tests accordingly).

## Work Objectives
### Core Objective
Create a new autoresearch harness mode that:
1) Runs DiHiggs explorers as subprocesses with config-defined hyperparameter bundles,
2) Extracts all produced `run_dir`s from explorer checkpoints,
3) Evaluates each `run_dir` with existing harness evaluator/scoring,
4) Logs harness-compatible `events.jsonl`,
5) Computes the composite score and adapts the next explorer/hparam bundle to maximize it while maintaining multi-axis diversity.

### Design Principles (non-negotiable)
- **Black-box integration**: explorers invoked only via subprocess argv lists.
- **Idempotent events**: each logical attempt has a stable `attempt_id`; runner must be safe to restart.
- **Sparse diversity**: per-axis + selected pairwise bins; no full grid enumeration.
- **Fail loud, fail actionable**: missing prerequisites must surface as clear preflight errors.

### Deliverables
- `python -m autoresearch.harness --mode dihiggs-explorers --config <json> ...` runnable.
- New config file: `autoresearch/configs/dihiggs_explorers.json` (schema defined below).
- New mode implementation module(s) under `autoresearch/harness/`.
- New parsers for:
  - Adaptive explorer checkpoints: `checkpoint_root/iter_XXXX/adaptive_state.json` → proposals[].run_dir (+ per-tb successes/trials)
  - Branch continuation artifacts: `checkpoint_root/<track_id>/events.jsonl` → result.run_dir (+ selected BranchState)
- Backwards-compatible multi-axis coverage/diversity in `autoresearch/benchmarks/score.py`.
- Tests:
  - Unit tests for parsing, event emission, metrics
  - Gated real E2E smoke that runs a minimal scan with `PhysScanWithFixings`.

### Definition of Done (agent-verifiable)
- [ ] Build succeeds and binary exists: `ls -l dihiggs/app/PhysScanWithFixings`
- [ ] Harness new mode runs and produces:
  - `autoresearch/opencode_logs/<campaign_id>/state/events.jsonl`
  - `autoresearch/opencode_logs/<campaign_id>/reports/summary.json`
- [ ] `python -m unittest discover -s autoresearch/tests -v` passes
- [ ] Real E2E smoke (gated) passes and produces at least one real `run_dir` with `run_manifest.json` + `task_summary.jsonl`

### Must Have
- Zero coupling via imports to DiHiggs explorers (strict subprocess calls).
- Robust parsing of checkpoints even with partial/failed runs.
- Deterministic, canonical `cell_id` encoding + axis extraction rules.
- Duplicate handling: duplicates contribute to penalty (counted as reuse).

### Must NOT Have (guardrails)
- No “rewrite the explorers”; treat `dihiggs/app/adaptive_explorer.py` and `branch_continuation_explorer_v2.py` as black boxes.
- No large-scale refactors across `dihiggs/` unrelated to running.
- No hardcoding user-specific paths in code (must be config-driven).

## Verification Strategy
> ZERO HUMAN INTERVENTION — all verification is agent-executed.
- Test decision: **tests-after** using existing **unittest** (repo already uses `python -m unittest discover`).
- CI parity checks (local): build C++ + run unit binary + run python unittests.
- Evidence policy: each task writes evidence to `.sisyphus/evidence/task-{N}-{slug}.txt` (or .log/.json as appropriate).

## Execution Strategy
### Parallel Execution Waves
Wave 1 (contracts + foundations)
- Metrics/axes contract + config schema
- Checkpoint parsers (adaptive + branch)
- Preflight/build checks (PhysScanWithFixings + datasets)

Wave 2 (new harness mode MVP)
- Implement mode runner: subprocess execution + parsing + evaluate_run_dir + event emission
- Backwards-compatible multi-axis compute_metrics updates

Wave 3 (adaptation + robustness)
- Strategy selection (bandit) + reward definition tied to composite score deltas
- Resume/idempotence + failure mapping + duplicate penalty wiring

Wave 4 (real end-to-end + docs + gated smoke)
- Minimal real E2E smoke + reporting
- Update README(s) with exact commands

### Dependency Matrix (high-level)
- Parsers + axis schema → event emission → compute_metrics update → adaptation loop → E2E smoke

### Agent Dispatch Summary
- Wave 1: 3 tasks (unspecified-high / deep)
- Wave 2: 2 tasks (deep)
- Wave 3: 2 tasks (deep)
- Wave 4: 2 tasks (unspecified-high)

## TODOs
> Implementation + Test = ONE task.
> EVERY task includes QA scenarios.

- [x] 1. Define config + axis contract for multi-axis diversity

  **What to do**:
  - Define the **new mode name**: `dihiggs-explorers`.
  - Define a JSON config schema for this mode (create `autoresearch/configs/dihiggs_explorers.json`). Use the **exact example** in “Appendix: Config Example” as the initial default.
  - Decide and encode the **axis contract** (axes are optional per attempt; missing axes are skipped for that axis metric).
    - Required physical axes (used for coverage + diversity):
      - `tb` (categorical; from config)
      - `lam1_bin` (int bin; domain size `search.lam1.n_bins`)
      - `mphi_bin` (int bin; optional; domain size `search.mphi.n_bins`)
    - Meta axes (diversity only, low weight):
      - `explorer` (categorical: `adaptive|branch`)
      - `strategy` (categorical: arm IDs)
    - Pairwise interactions (diversity only):
      - (`tb`,`lam1_bin`) and (`lam1_bin`,`mphi_bin`)
    - Additional metadata (not normalized; for debugging only): `track_id`, `lambda6`, `window_label`, `expansion_index`.
  - Define **payload encoding** for axes (preferred over string parsing):
    - `payload.axes_binned` (dict) contains binned/categorical axis values.
    - `payload.axes_raw` (dict) optionally contains raw floats (mphi, lam1).
  - Define **canonical cell_id format** (stable key order, pipe-separated) for backwards compatibility and human readability:
    - Always include `tb=<...>` and `bin=<lam1_bin>`.
    - May include extra segments, but metrics must primarily use `payload.axes_binned` when present.
  - Add config knobs for metrics weights and floors:
    - `metrics.weights` (use existing composite weights by default)
    - `metrics.floors` (lower default floors for coverage/diversity in multi-axis mode)
    - `metrics.multi_axis.coverage_axes`, `metrics.multi_axis.diversity_axes`, `metrics.multi_axis.diversity_pairs` (with weights).
  - Write a short schema doc block (in code docstring + README snippet) to remove ambiguity.

  **Must NOT do**:
  - Do not change any DiHiggs explorer CLI interfaces.

  **Recommended Agent Profile**:
  - Category: `deep` — Reason: requires consistent cross-module contracts and back-compat.

  **Parallelization**: Can Parallel: YES | Wave 1 | Blocks: [2,4] | Blocked By: []

  **References**:
  - Event wrapper and required payload keys: `autoresearch/harness/scheduler.py:80-87,433-448,534-548`
  - compute_metrics required payload keys: `autoresearch/benchmarks/score.py:131-171`
  - Current cell_id format: `autoresearch/harness/scheduler.py:132-164`

  **Acceptance Criteria**:
  - [ ] A concrete JSON example for `autoresearch/configs/dihiggs_explorers.json` exists in the codebase.
  - [ ] A deterministic `cell_id` construction rule is implemented and unit-tested.

  **QA Scenarios**:
  ```
  Scenario: Load config and build canonical cell_id
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_dihiggs_cell_id -v
    Expected: tests pass; cell_id matches exact expected string
    Evidence: .sisyphus/evidence/task-1-config-axis-contract.txt

  Scenario: Missing axis values are skipped (not counted as 'NA')
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_multiaxis_metrics_missing_axes -v
    Expected: metrics computed without exceptions; axis-specific coverage ignores missing axis
    Evidence: .sisyphus/evidence/task-1-missing-axes.txt
  ```

  **Commit**: YES | Message: `feat(autoresearch): define dihiggs explorers config + axis contract` | Files: [`autoresearch/configs/dihiggs_explorers.json`, tests + docs-in-code]

- [x] 2. Implement parsers: adaptive_explorer checkpoint → run_dir inventory

  **What to do**:
  - Add a parser module under `autoresearch/harness/` (e.g., `dihiggs_parsers.py`) with:
    - `parse_adaptive_checkpoint(iter_dir) -> list[AdaptiveRun]`
      - Load `adaptive_state.json`.
      - Iterate `proposals[]` and extract:
        - proposal_id, bin_index, lam1_min/max, run_dir, elapsed_sec, returncode
        - per-tb successes/trials from `successes_by_tb`/`trials_by_tb` when present
      - Return only proposals with `status == "DONE"` and `run_dir != null` for evaluation; keep failures for event emission.
    - `find_latest_iter_dirs(checkpoint_root)` and robust sorting by iter_XXXX.
  - Add unit tests using fixture JSON embedded in test file (no reliance on real checkpoints).

  **Must NOT do**:
  - Do not depend on importing `dihiggs/app/adaptive_explorer.py`.

  **Recommended Agent Profile**:
  - Category: `unspecified-high` — Reason: parsing + fixtures + edge cases.

  **Parallelization**: Can Parallel: YES | Wave 1 | Blocks: [5] | Blocked By: [1]

  **References**:
  - Checkpoint writer schema: `dihiggs/app/adaptive_explorer.py:495-538` (top-level keys)
  - Proposal keys: `dihiggs/app/adaptive_explorer.py:452-493`
  - Checkpoint location: `dihiggs/app/adaptive_explorer.py:_iter_dir 121-123`
  - run_dir marker regex: `dihiggs/app/adaptive_artifacts.py:15,241-255`

  **Acceptance Criteria**:
  - [ ] Parser returns correct run_dir list and per-tb metrics for fixture state.
  - [ ] Parser gracefully handles proposals with `run_dir=null`.

  **QA Scenarios**:
  ```
  Scenario: Parse adaptive_state.json fixture
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_parse_adaptive_checkpoint -v
    Expected: parsed proposals count, run_dir strings, and tb splits match expected
    Evidence: .sisyphus/evidence/task-2-parse-adaptive.txt

  Scenario: Corrupt JSON line / missing keys
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_parse_adaptive_checkpoint_errors -v
    Expected: parser skips/raises with clear error message (per design) and test asserts it
    Evidence: .sisyphus/evidence/task-2-parse-adaptive-errors.txt
  ```

  **Commit**: YES | Message: `feat(autoresearch): parse adaptive explorer checkpoints` | Files: [new parser + tests]

- [x] 3. Implement parsers: branch_continuation artifacts → run_dir inventory + metadata

  **What to do**:
  - Extend the same parser module with:
    - `parse_branch_track_events(track_dir) -> list[BranchAttempt]` reading `<track_dir>/events.jsonl`.
    - Extract from each JSON line:
      - track_id (dir name), step_label, point.tanbeta, point.lambda6
      - result.status, result.command, result.run_dir, result.csv_paths, result.elapsed_sec
      - selected state when present: selected.m_phi, selected.lam1, selected.csv_path, selected.row_index
  - Provide logic to discover all track dirs under a branch checkpoint_root.
  - Add unit tests using fixture JSONL lines embedded in tests.

  **Recommended Agent Profile**:
  - Category: `unspecified-high`

  **Parallelization**: Can Parallel: YES | Wave 1 | Blocks: [5] | Blocked By: [1]

  **References**:
  - events.jsonl write & schema: `dihiggs/app/branch_continuation_explorer_v2.py:1129-1176`
  - AttemptResult fields: `dihiggs/app/branch_continuation_explorer_v2.py:120-131`
  - run_dir parse usage: `dihiggs/app/branch_continuation_explorer_v2.py:996-1024`
  - CSV collection logic: `dihiggs/app/branch_continuation_explorer_v2.py:444-453`

  **Acceptance Criteria**:
  - [ ] Parser extracts run_dir and selected m_phi/lam1 when available.
  - [ ] Parser handles missing run_dir (orchestrator_failed) without crashing.

  **QA Scenarios**:
  ```
  Scenario: Parse branch events.jsonl fixture
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_parse_branch_events -v
    Expected: parsed attempts include track_id, run_dir, and selected state fields
    Evidence: .sisyphus/evidence/task-3-parse-branch.txt

  Scenario: Missing selected field
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_parse_branch_events_missing_selected -v
    Expected: parser returns attempt with selected=None; no exception
    Evidence: .sisyphus/evidence/task-3-parse-branch-missing-selected.txt
  ```

  **Commit**: YES | Message: `feat(autoresearch): parse branch continuation artifacts` | Files: [parser + tests]

- [x] 4. Add preflight: real end-to-end build & runtime checks for PhysScanWithFixings

  **What to do**:
  - Implement a preflight function invoked by the new mode before launching explorers:
    - Verify required system deps are present (best-effort): `libgsl` available; optional `gfortran`.
    - Verify env vars `HB_DATASET` and `HS_DATASET` are set and directories exist (if required by selected strategy).
    - Verify binary exists + executable: `dihiggs/app/PhysScanWithFixings`.
    - Verify build instructions are printed when missing (exact commands).
  - Add an optional `--preflight-only` flag to the new mode to run checks and exit 0/1.
  - Unit test: preflight should fail with clear error if binary missing (mock filesystem via temp dirs).

  **Recommended Agent Profile**:
  - Category: `deep` — Reason: cross-cutting runtime usability.

  **Parallelization**: Can Parallel: YES | Wave 1 | Blocks: [5,9] | Blocked By: [1]

  **References**:
  - Build steps:
    - `README.md:92-99` (2HDMC build)
    - `README.md:115-121` (HiggsTools build)
    - `README.md:146-148` (HB_DATASET/HS_DATASET)
    - `README.md:154-158` (dihiggs make)
  - Binary target location: `dihiggs/Makefile:33-43` and `dihiggs/Makefile:63-75`
  - CI env var check: `.github/scripts/ci_test.sh:30-33`

  **Acceptance Criteria**:
  - [ ] `python -m autoresearch.harness --mode dihiggs-explorers --preflight-only --config ...` returns non-zero with actionable message if missing prerequisites.

  **QA Scenarios**:
  ```
  Scenario: Preflight fails when PhysScanWithFixings missing
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_dihiggs_preflight_missing_binary -v
    Expected: clear error includes expected build commands and target path
    Evidence: .sisyphus/evidence/task-4-preflight-missing-binary.txt

  Scenario: Preflight passes with mocked binary and datasets
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_dihiggs_preflight_pass -v
    Expected: exit success
    Evidence: .sisyphus/evidence/task-4-preflight-pass.txt
  ```

  **Commit**: YES | Message: `feat(autoresearch): add dihiggs preflight checks` | Files: [mode runner + tests]

- [x] 5. Implement new harness mode runner (MVP): run subprocess explorers, parse checkpoints, evaluate run_dirs, emit events

  **What to do**:
  - Add a new mode handler called by `autoresearch/harness/__main__.py` when `--mode dihiggs-explorers`.
  - Implement the runner to:
    1) Load config.
    2) Run preflight.
    3) Execute exactly ONE strategy bundle (for MVP: the first arm in config) as subprocess.
    4) Parse checkpoints using tasks 2 & 3 parsers.
    5) Build a stable internal record per discovered run_dir:
       - `DiscoveredRun{arm_id, explorer, run_dir, discovered_at_utc, source_ref}` where source_ref points to (iter_dir, proposal_id) or (track_id, step_label, attempt_index).
    6) For each discovered `run_dir` (up to `limits.max_new_run_dirs_per_round`):
       - Call `autoresearch/harness/evaluator.evaluate_run_dir(run_dir)`.
       - Fingerprint artifacts using `autoresearch/harness/scoring.fingerprint_artifacts(eval_result.artifact_ref)`.
       - Duplicate policy (decision complete):
         - If fingerprint seen in prior rounds: mark all events for this run_dir as duplicate → emit `eval_status="SKIP_REUSED"`, `successes=0`, `trials=0`, `is_duplicate=true`.
         - If fingerprint new: process normally and only add fingerprint to seen-set **after** all slice events for this run_dir are emitted.
       - Emit events with **stable attempt_id**:
         - `attempt_id = sha256(campaign_id + arm_id + fingerprint + slice_key)`.
         - slice_key is `"tb=<tag>"` for adaptive per-tb slices; `"tb=<tag>|track=<track_id>|step=<step_label>"` for branch.
         - Ensure runner never emits the same attempt_id twice (resume safe).
       - Emit `ATTEMPT_CREATED`, `ATTEMPT_STARTED`, `ATTEMPT_EVALUATED` (and `ATTEMPT_SCORED/FINISHED` if utility computed).
    7) Write reports via existing reporting pipeline.
  - For adaptive explorer outputs, emit **per-tb slice** events:
    - For each proposal and each tb_tag in `successes_by_tb`, emit one `ATTEMPT_EVALUATED` with that split (successes/trials), sharing the same run_dir + fingerprint.
  - Define failure mapping:
    - If a run_dir is missing for an attempted orchestrator run → emit ATTEMPT_EVALUATED with `eval_status="CRASH"`, successes=0, trials=0, elapsed_sec from attempt metadata.
    - If run_dir exists but missing `run_manifest.json` or `task_summary.jsonl` → emit `eval_status="BROKEN_ARTIFACT"`.

  **Must NOT do**:
  - No imports from `dihiggs/app/*`.

  **Recommended Agent Profile**:
  - Category: `deep`

  **Parallelization**: Can Parallel: NO | Wave 2 | Blocks: [6,7,8,9] | Blocked By: [1,2,3,4]

  **References**:
  - Harness event wrapper schema: `autoresearch/harness/scheduler.py:80-87`
  - append_event JSONL writer: `autoresearch/harness/state.py:32-40`
  - compute_metrics payload requirements: `autoresearch/benchmarks/score.py:131-171`
  - Adaptive proposal run_dir key: `dihiggs/app/adaptive_explorer.py:466-470,452-493`
  - Branch attempt run_dir key: `dihiggs/app/branch_continuation_explorer_v2.py:1167-1176`

  **Acceptance Criteria**:
  - [ ] Running the new mode in a synthetic test harness (fixtures) produces a valid events.jsonl that `compute_metrics` can parse.
  - [ ] `events.jsonl` includes ATTEMPT_EVALUATED payload keys required by `compute_metrics`.

  **QA Scenarios**:
  ```
  Scenario: MVP mode runs against fixture checkpoints
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_dihiggs_mode_fixtures -v
    Expected: events.jsonl written; compute_metrics returns without errors
    Evidence: .sisyphus/evidence/task-5-mode-mvp-fixtures.txt

  Scenario: Failure mapping when run_dir missing
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_dihiggs_mode_missing_run_dir -v
    Expected: ATTEMPT_EVALUATED emitted with eval_status=CRASH, successes=0, trials=0
    Evidence: .sisyphus/evidence/task-5-mode-missing-run-dir.txt
  ```

  **Commit**: YES | Message: `feat(autoresearch): add dihiggs-explorers harness mode (MVP)` | Files: [new mode + tests]

- [x] 6. Extend compute_metrics for multi-axis coverage/diversity (backwards compatible)

  **What to do**:
  - Update `autoresearch/benchmarks/score.py` to support optional `config["metrics"]["multi_axis"]`.
    - If absent or disabled: preserve existing behavior exactly (tb/bin only).
    - If enabled: compute **sparse multi-axis** metrics from `ATTEMPT_EVALUATED` payloads using this decision-complete formula:
      - Prefer `payload.axes_binned` when present; fallback to parsing `cell_id` segments.
      - Physical per-axis coverage: tb, lam1_bin, mphi_bin.
      - Diversity: weighted mixture of per-axis entropies plus selected pairwise entropies:
        - Pairs: (tb, lam1_bin) and (lam1_bin, mphi_bin).
      - Missing axes are skipped for that axis/pair (do not count as "NA").
      - Duplicates (`payload.is_duplicate==true`) are excluded from coverage/diversity/yield aggregation.
      - Coverage aggregation:
        - For each configured coverage axis i with domain size Di and weight wi:
          - coverage_i = unique_bins_i / Di
        - coverage_norm = sum_i(wi * coverage_i) / sum_i(wi) (after skipping axes with no observations)
      - Diversity aggregation:
        - For each configured diversity axis j:
          - diversity_j = entropy(counts_j) / log(Dj)
        - For each configured diversity pair k with domain size Dk = Da*Db:
          - diversity_pair_k = entropy(counts_pairs_k) / log(Dk)
        - diversity_norm = weighted mean over axes+pairs (after skipping missing)
      - Collapse metric (max_cell_share) must use `metrics.multi_axis.collapse_axes` (default: tb+lam1_bin) rather than full cell_id.
    - Output extra diagnostics in metrics dict:
      - `coverage_by_axis`, `diversity_by_axis`, `diversity_by_pair`, `n_evaluated_unique`.
  - Duplicate policy in metrics:
    - If `payload.is_duplicate==true`, count as reuse (increase reuse_ratio) and exclude from yield/coverage/diversity.
  - Add unit tests:
    - Back-compat: existing fixtures still produce identical coverage/diversity for old-style `tb=...|bin=...`.
    - Multi-axis: known events produce expected axis coverages and diversity values.

  **Recommended Agent Profile**:
  - Category: `deep`

  **Parallelization**: Can Parallel: NO | Wave 2 | Blocks: [7] | Blocked By: [1,5]

  **References**:
  - Existing parsing: `autoresearch/benchmarks/score.py:131-171,236-250`
  - Current bin extraction: `autoresearch/benchmarks/score.py:_cell_bin_index 79-89`

  **Acceptance Criteria**:
  - [ ] Old benchmarks still pass their tests (no behavior change without new config).
  - [ ] Multi-axis metrics are computed and exposed in the returned metrics dict.

  **QA Scenarios**:
  ```
  Scenario: Backwards compatibility
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_benchmarks_score_backcompat -v
    Expected: metrics values match golden numbers
    Evidence: .sisyphus/evidence/task-6-metrics-backcompat.txt

  Scenario: Multi-axis metrics
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_benchmarks_score_multiaxis -v
    Expected: axis coverages/diversities match expected; duplicates excluded
    Evidence: .sisyphus/evidence/task-6-metrics-multiaxis.txt
  ```

  **Commit**: YES | Message: `feat(metrics): support multi-axis diversity/coverage + duplicate penalty` | Files: [score.py + tests]

- [x] 7. Implement adaptation loop (bandit) driven by composite score deltas

  **What to do**:
  - In the new mode runner, implement an outer loop for `max_rounds`:
    - Choose a strategy bundle (arm) using UCB1 over rewards in [0,1].
    - Reward definition per round:
      - Compute score_before = compute_composite(compute_metrics(events_path,...)).score
      - Run arm → emit new events
      - Compute score_after similarly
      - delta = score_after - score_before (range [-1,1])
      - reward = (delta + 1) / 2
    - Track per-arm (n, mean_reward) in an in-campaign state file.
  - Add config knobs:
    - `adaptation.ucb_c` (exploration constant)
    - `adaptation.warm_start_each_arm` (bool; default true)
  - Add unit tests with fake metrics to ensure arm selection is deterministic.

  **Recommended Agent Profile**:
  - Category: `deep`

  **Parallelization**: Can Parallel: NO | Wave 3 | Blocks: [8,9] | Blocked By: [5,6]

  **Acceptance Criteria**:
  - [ ] Given deterministic synthetic scores, the runner selects the expected arm sequence.

  **QA Scenarios**:
  ```
  Scenario: UCB chooses untried arm then best mean
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_dihiggs_ucb_selection -v
    Expected: arm selection order matches expected
    Evidence: .sisyphus/evidence/task-7-ucb.txt
  ```

  **Commit**: YES | Message: `feat(autoresearch): add bandit adaptation for dihiggs-explorers mode` | Files: [runner + tests]

- [ ] 8. Add resume/idempotence + robust subprocess logging

  **What to do**:
  - Persist a campaign-local state file (JSON) under `autoresearch/opencode_logs/<campaign_id>/state/dihiggs_mode_state.json` containing:
    - completed_round_ids
    - arm_stats
    - seen_fingerprints
    - mapping run_dir -> attempt_id(s)
    - emitted_attempt_ids
    - per-arm checkpoint cursors (last iter_XXXX processed; last track_ids processed)
  - On startup with `--resume`, load state + replay events.jsonl to reconstruct counters.
  - Ensure subprocess stdout/stderr are always captured to files under `attempt_logs_root`.
  - Add checkpoint parsing robustness:
    - Only parse checkpoint directories that are “stable” (mtime unchanged for N seconds) when resuming after a crash.
  - Add a hard cap to avoid infinite loops when explorers produce no run_dirs:
    - `max_empty_rounds` config; after N empty rounds, abort with non-zero.

  **Recommended Agent Profile**:
  - Category: `deep`

  **Parallelization**: Can Parallel: NO | Wave 3 | Blocked By: [5,7] | Blocks: [9]

  **Acceptance Criteria**:
  - [ ] If runner is interrupted mid-campaign, `--resume` continues without duplicating events.

  **QA Scenarios**:
  ```
  Scenario: Resume does not duplicate attempt_ids
    Tool: Bash
    Steps: python -m unittest autoresearch.tests.test_dihiggs_resume_idempotent -v
    Expected: events.jsonl length unchanged after resume; new attempts appended only once
    Evidence: .sisyphus/evidence/task-8-resume.txt
  ```

  **Commit**: YES | Message: `feat(autoresearch): add resume/idempotence for dihiggs-explorers mode` | Files: [runner + tests]

- [ ] 9. Real end-to-end smoke: build binary + run minimal campaign using explorers

  **What to do**:
  - Add a gated integration test (unittest) controlled by env var `AUTORESEARCH_ENABLE_DIHIGGS_E2E=1` that:
    1) Asserts `dihiggs/app/PhysScanWithFixings` exists (or runs preflight and skips with message).
    2) Runs `python -m autoresearch.harness --mode dihiggs-explorers --config autoresearch/configs/dihiggs_explorers.json --iters 1 --attempts 1` (mapping: iters=max_rounds, attempts=max_new_run_dirs_per_round).
    3) Verifies at least one `run_dir` exists with `run_manifest.json` + `task_summary.jsonl`.
    4) Verifies reports exist.
  - Ensure the default config chooses **tiny budgets** and small grids to keep runtime bounded.

  **Recommended Agent Profile**:
  - Category: `unspecified-high` — Reason: environment + heavy deps + gating.

  **Parallelization**: Can Parallel: NO | Wave 4 | Blocked By: [4,5,7,8] | Blocks: [10]

  **References**:
  - Real build sequence: `README.md:92-99,115-121,154-158`
  - Orchestrator-run expectation: explorers rely on `dihiggs/app/orchestrate_scans.py` and stdout marker `[PATH] run_dir = ...` (`dihiggs/app/adaptive_artifacts.py:15`).

  **Acceptance Criteria**:
  - [ ] With env var enabled and binary present, the integration test passes.

  **QA Scenarios**:
  ```
  Scenario: Run gated real E2E
    Tool: Bash
    Steps: AUTORESEARCH_ENABLE_DIHIGGS_E2E=1 python -m unittest autoresearch.tests.test_dihiggs_e2e_real -v
    Expected: test passes; prints run_dir; artifacts exist
    Evidence: .sisyphus/evidence/task-9-e2e-real.txt
  ```

  **Commit**: YES | Message: `test(autoresearch): add gated real E2E smoke for dihiggs-explorers mode` | Files: [test + config tweaks]

- [ ] 10. Documentation: exact run commands + troubleshooting

  **What to do**:
  - Update `autoresearch/README.md` with:
    - Build prerequisites (GSL, gfortran) and exact build commands for 2hdmc/higgstools/dihiggs.
    - Dataset env vars HB_DATASET/HS_DATASET.
    - Example invocation for new mode with real end-to-end.
    - How to add new strategy bundles.
    - Common failure modes (missing run_dir marker, missing datasets, permission issues on outdir).

  **Recommended Agent Profile**:
  - Category: `writing`

  **Parallelization**: Can Parallel: YES | Wave 4 | Blocked By: [9] | Blocks: []

  **Acceptance Criteria**:
  - [ ] README includes a copy-paste run sequence that results in a campaign directory with reports.

  **QA Scenarios**:
  ```
  Scenario: Docs include runnable commands
    Tool: Bash
    Steps:
      - ~/higgs_env_py312/bin/python -c "from pathlib import Path; t=Path('autoresearch/README.md').read_text(); assert '--mode dihiggs-explorers' in t"
      - ~/higgs_env_py312/bin/python -c "import json; json.load(open('autoresearch/configs/dihiggs_explorers.json'))"
      - ~/higgs_env_py312/bin/python -m autoresearch.harness --help
    Expected: all commands exit 0; README references the new mode; config JSON is valid
    Evidence: .sisyphus/evidence/task-10-docs-runbook.txt
  ```

  **Commit**: YES | Message: `docs(autoresearch): add dihiggs-explorers runbook` | Files: [`autoresearch/README.md`]

## Final Verification Wave (MANDATORY)
> Run 4 review agents in PARALLEL. ALL must approve. Present consolidated results to user and wait for explicit “okay”.

- [ ] F1. Plan Compliance Audit — oracle

  **What to do**:
  - Run an oracle review against the implemented changes and this plan.

  **QA Scenarios**:
  ```
  Scenario: Oracle approves plan compliance
    Tool: task(subagent_type="oracle")
    Steps: Prompt oracle with: "Review implementation against .sisyphus/plans/autoresearch-dihiggs-explorers.md. Return APPROVE or REQUEST_CHANGES with critical issues only."
    Expected: APPROVE, or REQUEST_CHANGES is addressed and review re-run until APPROVE
    Evidence: .sisyphus/evidence/f1-oracle-plan-compliance.txt
  ```

- [ ] F2. Code Quality Review — unspecified-high

  **What to do**:
  - Run an independent review focusing on subprocess safety, parsing robustness, and idempotence.

  **QA Scenarios**:
  ```
  Scenario: Reviewer finds no critical issues
    Tool: task(category="unspecified-high")
    Steps: Prompt reviewer with: "Review changes for robustness (timeouts, partial checkpoints, duplicates, resume). Return APPROVE or REQUEST_CHANGES."
    Expected: APPROVE, or REQUEST_CHANGES is addressed and review re-run
    Evidence: .sisyphus/evidence/f2-code-quality-review.txt
  ```

- [ ] F3. Real Manual QA — unspecified-high

  **What to do**:
  - Execute the real minimal run (not just unit tests) and verify artifacts/reports.

  **QA Scenarios**:
  ```
  Scenario: Minimal real run succeeds
    Tool: Bash
    Steps:
      - ~/higgs_env_py312/bin/python -m autoresearch.harness --mode dihiggs-explorers --config autoresearch/configs/dihiggs_explorers.json --iters 1 --attempts 1 --preflight-only
      - AUTORESEARCH_ENABLE_DIHIGGS_E2E=1 ~/higgs_env_py312/bin/python -m unittest autoresearch.tests.test_dihiggs_e2e_real -v
    Expected: preflight exits 0; gated E2E test passes
    Evidence: .sisyphus/evidence/f3-real-manual-qa.txt
  ```

- [ ] F4. Scope Fidelity Check — deep

  **What to do**:
  - Verify no scope creep and that required deliverables are present.

  **QA Scenarios**:
  ```
  Scenario: Scope fidelity approved
    Tool: task(category="deep")
    Steps: Prompt reviewer with: "Check implemented changes match plan scope (no extra refactors, no dihiggs explorer rewrites). Return APPROVE or REQUEST_CHANGES."
    Expected: APPROVE, or REQUEST_CHANGES is addressed and re-run
    Evidence: .sisyphus/evidence/f4-scope-fidelity.txt
  ```

## Commit Strategy
- Commit 1: config + axis contract + cell_id utilities + unit tests
- Commit 2: checkpoint parsers + fixtures tests
- Commit 3: new harness mode MVP (single strategy) + tests
- Commit 4: compute_metrics multi-axis + duplicate penalty + tests
- Commit 5: adaptation (bandit) + resume/idempotence + tests
- Commit 6: gated real E2E smoke + docs

## Success Criteria
- Running the following produces a non-empty campaign with a composite score and reports:
  - `~/higgs_env_py312/bin/python -m autoresearch.harness --mode dihiggs-explorers --config autoresearch/configs/dihiggs_explorers.json --iters 5 --attempts 10`
- Multi-axis metrics are visible in the reported summary (coverage_norm/diversity_norm derived from multiple axes, not only lam1).
- Real binary build and minimal real scan are verified by the gated E2E test.

## Appendix: Config Example (authoritative initial default)
> This JSON must be copied into `autoresearch/configs/dihiggs_explorers.json` (with paths adjusted). The runner must support `{placeholders}` exactly as specified.

```json
{
  "campaign_id": "dihiggs_explorers_dev",
  "paths": {
    "repo_root": "/home/fabi/wt_dihiggs_exploratory",
    "outdir": "/tmp/dihiggs_runs",
    "lake_name": "dihiggs_lake"
  },
  "runtime": {
    "python_exe": "{sys_executable}",
    "threads": 4,
    "timeout_sec": 7200,
    "max_empty_rounds": 2
  },
  "dihiggs": {
    "phys_exec": "/home/fabi/wt_dihiggs_exploratory/dihiggs/app/PhysScanWithFixings",
    "hb_dataset_env": "HB_DATASET",
    "hs_dataset_env": "HS_DATASET"
  },
  "search": {
    "tb_values": [10000, 15000, 20000],
    "lam1": {"min": 0.0, "max": 12.0, "n_bins": 12},
    "mphi": {"min": 120.0, "max": 300.0, "n_bins": 10}
  },
  "limits": {
    "max_new_run_dirs_per_round": 10,
    "max_new_run_dirs_per_arm_call": 50
  },
  "metrics": {
    "floors": {"yield_norm": 0.05, "coverage_norm": 0.05, "diversity_norm": 0.05},
    "weights": {
      "yield_norm": 0.40,
      "coverage_norm": 0.20,
      "diversity_norm": 0.15,
      "efficiency_norm": 0.10,
      "reuse_penalty": 0.10,
      "failure_penalty": 0.05
    },
    "multi_axis": {
      "enabled": true,
      "collapse_axes": ["tb", "lam1_bin"],
      "coverage_axes": [
        {"name": "tb", "kind": "categorical", "domain": [10000, 15000, 20000], "weight": 0.40},
        {"name": "lam1_bin", "kind": "int", "domain_size": 12, "weight": 0.30},
        {"name": "mphi_bin", "kind": "int", "domain_size": 10, "weight": 0.30}
      ],
      "diversity_axes": [
        {"name": "tb", "weight": 0.20},
        {"name": "lam1_bin", "weight": 0.20},
        {"name": "mphi_bin", "weight": 0.20},
        {"name": "explorer", "kind": "categorical", "domain": ["adaptive", "branch"], "weight": 0.05},
        {"name": "strategy", "kind": "categorical", "domain": ["adaptive-v1", "branch-v1"], "weight": 0.05}
      ],
      "diversity_pairs": [
        {"axes": ["tb", "lam1_bin"], "weight": 0.15},
        {"axes": ["lam1_bin", "mphi_bin"], "weight": 0.15}
      ]
    }
  },
  "adaptation": {
    "ucb_c": 1.0,
    "warm_start_each_arm": true
  },
  "arms": [
    {
      "id": "adaptive-v1",
      "explorer": "adaptive",
      "timeout_sec": 7200,
      "cmd": [
        "{python}",
        "dihiggs/app/adaptive_explorer.py",
        "--lam1-min", "0.0",
        "--lam1-max", "12.0",
        "--n-bins", "12",
        "--seed", "0",
        "--n-coarse", "5",
        "--total-budget", "30",
        "--floor-points", "1",
        "--n-iters", "1",
        "--checkpoint-root", "{checkpoint_root}",
        "--exec", "{phys_exec}",
        "--outdir", "{outdir}",
        "--lake-name", "{lake_name}",
        "--campaign", "{campaign}",
        "--threads", "{threads}"
      ],
      "env": {"OMP_NUM_THREADS": "4"}
    },
    {
      "id": "branch-v1",
      "explorer": "branch",
      "timeout_sec": 7200,
      "cmd": [
        "{python}",
        "dihiggs/app/branch_continuation_explorer_v2.py",
        "--mode", "level-curve",
        "--tb0", "10000",
        "--lambda6-0", "0.10",
        "--n-up", "1",
        "--n-down", "1",
        "--exec", "{phys_exec}",
        "--outdir", "{outdir}",
        "--lake-name", "{lake_name}",
        "--campaign", "{campaign}",
        "--checkpoint-root", "{checkpoint_root}",
        "--local-mphi-half-span", "8.0",
        "--local-lam1-half-span", "0.5",
        "--n-mphi-local", "5",
        "--n-lam1-local", "10",
        "--max-expansions", "1",
        "--expansion-factor", "2.0",
        "--seed-policy", "center",
        "--threads", "{threads}"
      ],
      "env": {"OMP_NUM_THREADS": "4"}
    }
  ]
}
```

**Placeholder substitution contract (must implement)**
- `{sys_executable}`: resolved at runtime to `sys.executable` when loading config.
- `{python}`: `runtime.python_exe` (after resolving `{sys_executable}`) unless overridden per-arm.
- `{checkpoint_root}`: a unique directory per round+arm under campaign logs, e.g. `.../experiments/round_0003/arm_adaptive-v1/checkpoints`.
- `{campaign}`: a deterministic sub-campaign name (campaign_id + round + arm_id).
- `{phys_exec}`, `{outdir}`, `{lake_name}`, `{threads}`: from config/runtime.
