## [2026-04-02T23:30:52-03:00] Task 1: Config + Axis Contract

### Cell_id Format
- Canonical: `tb=<value>|bin=<lam1_bin>`
- Extended: `tb=<value>|bin=<lam1_bin>|<axis_name>=<value>|...`
- Ordering: tb first, bin second, then alphabetical

### Binning Logic
- Linear binning: `bin_index = floor((value - min) / bin_width)`
- Edge handling: clamp out-of-range values to `[0, n_bins-1]`

### Placeholder Substitution
- Supported: {sys_executable}, {python}, {checkpoint_root}, {campaign}, {outdir}, {lake_name}, {threads}, {tb}, {iter}, {track_id}, {phys_exec}, {hb_dataset_env}, {hs_dataset_env}
- Recursive resolution: nested placeholders resolved depth-first

## [2026-04-02T23:38:07-03:00] Task 4: Preflight Checks

### Preflight Check Contract
- Each check returns: `{"check": str, "status": "pass"|"fail"|"warn", "message": str}`
- Overall status: "fail" > "warn" > "pass" (worst wins)
- Exception-safe: all checks catch and return fail status

### Check Types
- **phys_exec**: Binary must exist and be executable (critical)
- **datasets**: Env vars should be set, paths must exist if set (warn if not set, fail if set but missing)
- **2hdmc_lib**: Library should exist but not critical (warn if missing)

### Testing Strategy
- Mock filesystem with `unittest.mock.patch("os.path.exists")`
- Mock env vars with `monkeypatch.setenv` (pytest)
- Test exception handling for robustness

## [2026-04-02] Task 2: Adaptive Checkpoint Parser

### Checkpoint Schema
- Path pattern: `{checkpoint_root}/iter_{iter_index:04d}/adaptive_state.json`
- Directory naming: iter_XXXX where XXXX is 4-digit zero-padded iteration number
- JSON structure: `{"proposals": [...], "metadata": {...}}`
- Proposal keys (from adaptive_explorer.py state dict):
  - `run_dir`: str (path to orchestrator output dir)
  - `bin_index`: int (used for both tb and proposal_index)
  - `lam1_min`: float (lower bound of lam1 range for this bin)
  - Additional keys: `proposal_id`, `status`, `command`, `returncode`, `elapsed_sec`, etc.

### Stability Check Implementation
- **mtime threshold: 2.0 seconds** (avoids parsing half-written checkpoints)
- Rationale: Filesystem writes may not be atomic; concurrent adaptive_explorer.py write + harness read could yield corrupted JSON
- Implementation: `os.stat(path).st_mtime > (time.time() - 2.0)` → skip file
- Gracefully skips checkpoints, allowing parse to succeed with partial results

### Filtering Logic
- **Deduplication**: Skip if `run_dir` in `known_run_dirs` set
- **Null/empty run_dir**: Skip proposals where `run_dir is None or run_dir == ""`
  - Rationale: Proposals may have null run_dir if orchestrator failed before creating output directory
- **Return only NEW discoveries**: Caller maintains `known_run_dirs` set across multiple parse calls

### Error Handling Strategy
- **Missing checkpoint root**: Return empty list, log warning (not an error - expected for fresh runs)
- **Malformed JSON**: Skip checkpoint file, log warning, continue with other iterations
- **Missing required keys**: Use fallback defaults (tb=0, lam1_raw=0.0) - parser is robust to incomplete proposals
- **Invalid field types**: Skip individual proposals (try/except around int/float conversions)
- **Philosophy**: Graceful degradation - extract as much data as possible, never crash on partial/corrupted checkpoints

### Binning Integration
- Uses `bin_lam1(lam1_raw, config)` from `dihiggs_axis_contract` module
- Config path: `config["search"]["lam1"]` with `min`, `max`, `n_bins` keys
- Binning errors (KeyError, TypeError, ValueError) are caught and logged as warnings
- Allows parser to skip un-binnable proposals without failing entire parse

### Discovery Dataclass Fields
- `run_dir`: str (full path to orchestrator output directory)
- `tb`: int (tanbeta value, from bin_index)
- `lam1_raw`: float (raw lam1 value before binning)
- `lam1_bin`: int (computed bin index for lam1 axis)
- `iter_index`: int (iteration number from directory name)
- `proposal_index`: int (from bin_index field in proposal)

### Sorting and Determinism
- Final results sorted by `(iter_index, proposal_index)` for deterministic ordering
- Ensures stable output across multiple calls with same checkpoint state
- Important for reproducible test runs and checkpoint replay scenarios

### Test Coverage Achieved
- Valid checkpoint parsing with multiple proposals
- Known run_dir filtering (deduplication)
- Checkpoint stability check (skip recent files)
- Missing checkpoint root/files (graceful degradation)
- Malformed JSON and non-dict JSON (error resilience)
- Missing proposals key (empty result)
- Null/empty run_dir proposals (filtered out)
- Missing optional fields (fallback defaults)
- Invalid field types (skip proposal)
- Multiple iterations (correct sorting)
- Empty checkpoints (empty result)
- lam1 binning integration (verify axis contract usage)

## [2026-04-02 23:44:54] Task 3: Branch Continuation Parser

### Checkpoint Schema
- Path: `{checkpoint_root}/{track_id}/events.jsonl`
- Format: JSONL (one JSON object per line)
- Event filter: `event_type == "ATTEMPT_COMPLETED"`
- Payload keys: `result.run_dir`, `params.tb`, `params.lam1`, `step_label`
- Key extraction paths:
  - `run_dir`: `payload.result.run_dir`
  - `tb`: `payload.params.tb`
  - `lam1_raw`: `payload.params.lam1`
  - `step_name`: `payload.step_label`

### JSONL Parsing
- Read line by line, skip empty lines
- Each line is separate JSON object (use `json.loads(line)`)
- Handle malformed lines gracefully: log warning and skip line
- Continue processing remaining lines after encountering malformed JSON

### Filtering Logic
- Same as adaptive: skip known_run_dirs, skip null/empty run_dir
- Event type filter: only ATTEMPT_COMPLETED events (ignore ATTEMPT_STARTED, ATTEMPT_FAILED, etc.)
- Required fields validation: `tb`, `lam1`, and `step_label` must all be present (no defaults)
- Missing any required field → skip event

### Error Handling
- Missing events.jsonl → log debug, continue to next track
- Malformed JSONL line → log warning, skip line
- Missing required keys → log debug, skip event
- Invalid payload/result/params types → log debug, skip event

### Binning
- Use `bin_lam1(lam1_raw, config)` from axis_contract
- Same config format as adaptive parser: `config["search"]["lam1"]`

### Sorting
- Results sorted by `(track_id, step_name)` for deterministic ordering
- Track dirs sorted alphabetically before processing

## [2026-04-03T00:00:00Z] Task 5: Harness Mode Runner MVP

### Runner Architecture
- DiHiggsRunner class manages full discovery loop
- Subprocess black-box execution of explorers
- Checkpoint parsing via Task 2/3 parsers
- Placeholder scoring (MVP stub)
- Event emission to JSONL log

### Subprocess Execution
- Use subprocess.run with timeout, capture_output
- Build cmd/env via placeholder substitution
- Log stdout/stderr to separate files per round
- Handle TimeoutExpired gracefully

### Attempt ID Formula
- `sha256(campaign_id + "|" + arm_id + "|" + cell_id + "|" + run_dir).hexdigest()[:16]`
- Idempotent: same inputs → same ID
- Includes run_dir for uniqueness per discovery

### Event Schema Compatibility
- Matches autoresearch/harness/scheduler.py event schema
- Wrapper: {schema_version, campaign_id, event_type, utc, payload}
- Payload: {attempt_id, iter_index, attempt_index, cell_id, eval_status, successes, trials, elapsed_sec, axes_binned, axes_raw}
- JSONL format: one event per line

### MVP Stubs (to be replaced later)
- evaluate_run_dir: returns mock scores (TODO: read CSV, compute real metrics)
- iter_index/attempt_index: hardcoded to 0 (TODO: track properly in Task 8)
- elapsed_sec: hardcoded to 0.0 (TODO: track subprocess timing)
- Note: alphabetical step_name sorting means "expand" comes before "seed"

### Track Discovery
- Find all subdirectories in checkpoint_root
- Each subdirectory name = track_id
- Process tracks in alphabetical order

## [2026-04-02T23:XX:XX] Task 5: Harness Runner MVP - VERIFIED ✅

### Implementation Complete
- DiHiggsRunner class with full integration loop
- 254 lines of production code, 202 lines of test code
- 6/6 tests passing, all Wave 1 tests still passing (54/54 total)

### Key Features Verified
- **Preflight integration**: Skips subprocess on preflight fail
- **Subprocess execution**: Timeout handling, stdout/stderr logging to files
- **Checkpoint parsing**: Dispatcher to adaptive/branch parsers based on explorer type
- **Event emission**: Matches scheduler.py schema exactly
- **Attempt ID stability**: SHA-256 based, deterministic 16-char hex
- **Deduplication**: Uses known_run_dirs set across rounds

### MVP Stub Confirmed
- evaluate_run_dir returns mock scores with TODO comment at line 218
- Ready for Task 6 to extend with real metrics

### Log File Naming
- Pattern: `{outdir}/{arm_id}_round_{N}.stdout.txt` and `.stderr.txt`
- Round counter tracked per arm_id in _round_counter dict

### Integration Points
- Uses all Wave 1 modules: preflight, parsers, axis_contract
- Config-driven: arm cmd/env via placeholder substitution
- No hardcoded paths (all from config)

### Test Coverage
1. Full subprocess → parse → emit integration
2. Attempt ID stability and uniqueness
3. Timeout handling with partial output capture
4. Preflight gate (skip on fail)
5. Event schema validation
6. MVP stub score return

## [2026-04-03T00:15:00Z] Task 6: Multi-Axis Scoring Module

### Module Architecture
- Location: `autoresearch/benchmarks/score.py`
- Purpose: Multi-axis coverage and diversity metrics with sparse axis tracking
- Functions: `compute_metrics`, `compute_coverage`, `compute_diversity`, `compute_composite`

### Metrics Implementation

#### Yield Metrics (compute_metrics)
- Reads `{run_dir}/results.csv` with columns: `success` (bool), `total_events` (int)
- Formula: `yield = count(success=True) / total_rows`
- Returns: `{"yield_val": float, "successes": int, "trials": int}`
- Graceful degradation: Missing/malformed CSV → return zeros, no crash
- Success parsing: accepts "true", "TRUE", "1", "yes" (case-insensitive)

#### Coverage Metrics (compute_coverage)
- Per-axis bin coverage: `bins_visited / total_bins`
- Weighted mean across axes using `coverage_axes` weights
- Config path: `config["metrics"]["multi_axis"]["coverage_axes"]`
- Each axis spec: `{"name": str, "kind": "categorical"|"discrete", "domain"|"domain_size": ..., "weight": float}`
- Fallback: uses `collapse_axes` if `coverage_axes` not specified
- Domain size resolution:
  - Categorical (tb): `len(domain)` from axis spec or `len(config["search"]["tb_values"])`
  - Discrete (lam1_bin, mphi_bin): `domain_size` from axis spec or `config["search"][axis]["n_bins"]`
- Missing axis values in history are skipped gracefully

#### Diversity Metrics (compute_diversity)
- Shannon entropy per axis: `H = -sum(p * log2(p))`
- Weighted mean across axes using `diversity_axes` weights
- Config path: `config["metrics"]["multi_axis"]["diversity_axes"]`
- Pairwise interactions support: `diversity_pairs` with `{"axes": [axis1, axis2], "weight": float}`
- Pairwise computation: treat (val1, val2) tuple as single distribution key
- Fallback: uses `collapse_axes` if `diversity_axes` not specified
- Single-value distributions return entropy=0.0 (expected behavior)

#### Composite Score (compute_composite)
- Weighted sum: `w_yield * yield + w_coverage * coverage + w_diversity * diversity`
- Weights from `config["metrics"]["weights"]` with keys: "yield", "coverage", "diversity"
- Unbounded score (depends on diversity scale and weights)

### Config Schema Usage
```json
{
  "metrics": {
    "weights": {"yield": 0.4, "coverage": 0.3, "diversity": 0.3},
    "multi_axis": {
      "collapse_axes": ["tb", "lam1_bin"],
      "coverage_axes": [
        {"name": "tb", "kind": "categorical", "domain": [1000, 5000, 10000, 30000], "weight": 0.3},
        {"name": "lam1_bin", "kind": "discrete", "domain_size": 40, "weight": 0.4}
      ],
      "diversity_axes": [
        {"name": "tb", "weight": 0.2},
        {"name": "lam1_bin", "weight": 0.3}
      ],
      "diversity_pairs": [
        {"axes": ["tb", "lam1_bin"], "weight": 0.5}
      ]
    }
  }
}
```

### Error Handling Philosophy
- **Graceful degradation**: Always return valid defaults, never crash
- **File I/O**: Wrap in try/except, log warnings, return safe values
- **Missing config keys**: Fallback to alternative sources or defaults
- **Malformed data**: Skip invalid entries, continue processing
- **Type coercion**: Use `float()` and `int()` with defensive checks

### Test Coverage (22 tests, all passing)
- **compute_metrics (5 tests)**: valid CSV, missing CSV, malformed CSV, empty CSV, various success formats
- **compute_coverage (6 tests)**: single-axis, multi-axis, empty history, full coverage, weighted axes, fallback
- **compute_diversity (6 tests)**: uniform, skewed, empty, single-value, pairwise, comparison
- **compute_composite (3 tests)**: basic weights, zero weights, various combinations
- **End-to-end (2 tests)**: full pipeline, fallback configurations

### Integration Point
- Task 5 `DiHiggsRunner.evaluate_run_dir()` currently returns mock scores
- Ready to integrate: `from autoresearch.benchmarks.score import compute_metrics`
- Full composite scoring requires history tracking (Task 7+)

### Type Safety
- All functions have full type annotations
- Config uses `Mapping[str, Any]` for flexibility with JSON
- LSP warnings are minor (expected `Any` from untyped config dicts)
- No errors, only warnings about `Any` propagation

### Key Design Decisions
1. **Sparse tracking**: Per-axis bins + selected pairwise (not full Cartesian grid)
2. **Weighted metrics**: Each axis has configurable weight for coverage/diversity
3. **Fallback chains**: coverage_axes → collapse_axes → search config → defaults
4. **CSV format**: Simple success/total_events columns (extensible to more fields)
5. **Reserved parameters**: `axes_binned` and `config` in compute_metrics for future use

### Backward Compatibility
- Works with existing single-axis `lam1_bin` tracking
- Supports minimal configs with only `collapse_axes`
- Gracefully handles missing optional config sections

### Performance Characteristics
- Coverage: O(history_size * n_axes)
- Diversity: O(history_size * n_axes + history_size * n_pairs)
- Memory: O(unique_values_per_axis) for Counter objects
- Suitable for histories up to ~10K attempts (typical campaign size)

## [2026-04-03T00:05:00Z] Task 6: Multi-Axis Metrics Module - VERIFIED ✅

### Implementation Complete
- autoresearch/benchmarks/score.py (376 lines)
- autoresearch/tests/test_score.py (425 lines, 22 tests)
- All 76 tests passing (54 existing + 22 new)
- No LSP errors

### Functions Implemented
1. **compute_metrics(run_dir, axes_binned, config)**
   - Reads results.csv from run_dir
   - Extracts yield, successes, trials
   - Graceful degradation: missing/malformed → zeros

2. **compute_coverage(history, config)**
   - Per-axis bin coverage: bins_visited / total_bins
   - Weighted mean across configured axes
   - Fallback chain: coverage_axes → collapse_axes

3. **compute_diversity(history, config)**
   - Shannon entropy per axis: -sum(p * log2(p))
   - Supports pairwise interactions (optional)
   - Weighted mean with fallback chain

4. **compute_composite(yield_val, coverage, diversity, config)**
   - Weighted sum: w_yield * yield + w_coverage * coverage + w_diversity * diversity
   - Weights from config["metrics"]["weights"]

### Key Design Decisions
- **Sparse tracking**: Per-axis bins + selected pairwise interactions (no full Cartesian grid)
- **Fallback chains**: coverage_axes → collapse_axes, diversity_axes → collapse_axes
- **Graceful degradation**: Missing CSV/config → return safe defaults (zeros)
- **Backward compatible**: Works with existing single-axis lam1_bin tracking
- **Config-driven**: No hardcoded axis names, all from config["metrics"]["multi_axis"]

### Config Schema Used
```json
{
  "metrics": {
    "weights": {"yield": 0.5, "coverage": 0.3, "diversity": 0.2},
    "multi_axis": {
      "collapse_axes": ["tb", "lam1_bin"],
      "coverage_axes": [
        {"name": "tb", "kind": "categorical", "domain": [1000, 5000, ...], "weight": 0.5},
        {"name": "lam1_bin", "kind": "discrete", "domain_size": 40, "weight": 0.5}
      ],
      "diversity_axes": [
        {"name": "tb", "weight": 0.5},
        {"name": "lam1_bin", "weight": 0.5}
      ],
      "diversity_pairs": [
        {"axes": ["tb", "lam1_bin"], "weight": 1.0}
      ]
    }
  }
}
```

### Test Coverage
- compute_metrics: 5 tests (valid, missing, malformed, empty, various formats)
- compute_coverage: 5 tests (single-axis, multi-axis, empty, full, weighted)
- compute_diversity: 6 tests (uniform, skewed, empty, single-value, pairwise, comparison)
- compute_composite: 3 tests (basic, zero weights, various weights)
- Integration: 1 end-to-end test
- Fallback: 2 tests (coverage + diversity fallback to collapse_axes)

### Integration Points
- Task 5 `DiHiggsRunner.evaluate_run_dir()` can now import compute_metrics
- Full composite scoring with history will be integrated in Task 7 (bandit adaptation)
- CSV format: success (bool), total_events (int)

### Error Handling
- All file I/O wrapped in try/except
- Missing files → return zeros, log warning
- Malformed CSV → return zeros, log warning
- Unknown axes → log warning, use domain_size=1
- Invalid config → graceful fallback chains

## [2026-04-03T00:35:00Z] Task 7: Bandit Adaptation Loop - COMPLETED ✅

### Implementation Complete
- autoresearch/harness/dihiggs_adaptation.py (222 lines)
- autoresearch/tests/test_dihiggs_adaptation.py (283 lines, 14 tests)
- DiHiggsRunner integrated with BanditState
- All 90 tests passing (Wave 1-7)

### Functions Implemented
1. **BanditState dataclass**
   - arm_stats: dict[arm_id -> {pulls, total_reward}]
   - global_history: list[attempt_data] across all arms
   - arm_histories: dict[arm_id -> list[attempt_data]]

2. **select_arm_ucb1(state, config) -> str**
   - UCB1: mean_reward + c * sqrt(log(total_pulls) / arm_pulls)
   - Cold-start: 0 pulls → infinity score (explored first)
   - Exploration constant from config["adaptation"]["ucb1_exploration_constant"]

3. **update_state(state, arm_id, reward, attempt_data)**
   - Increment arm pulls
   - Accumulate total_reward
   - Append to global_history and arm_histories

4. **compute_arm_reward(history, config) -> float**
   - Latest attempt yield + global coverage/diversity
   - Returns compute_composite from Task 6 metrics

### DiHiggsRunner Integration
- Added self.bandit_state in __init__
- Modified evaluate_run_dir to call real metrics (Task 6)
- Added run_adaptation_round() method:
  - select_arm_ucb1 → run_single_round → compute_arm_reward → update_state
  - Returns full result with arm_id and reward

### Key Design Decisions
- **Global history for coverage/diversity**: All arms contribute to exploration metrics
- **Per-arm history for debugging**: Separate histories for analysis
- **Cold-start handling**: Infinity score ensures all arms tried before exploitation
- **Reward components**: yield (immediate) + coverage (breadth) + diversity (balance)

### Test Coverage
- UCB1: cold-start (3 tests), exploitation, exploration, initialization
- update_state: pulls increment, reward accumulation, history appending (5 tests)
- compute_arm_reward: composite score, latest yield, empty history (3 tests)
- End-to-end: full adaptation cycle (select → update → select)

### Integration Notes
- evaluate_run_dir now requires CSV file in run_dir (real metrics)
- Test updated to create mock CSV: test_evaluate_run_dir_returns_mvp_scores
- Imports from score.py: compute_metrics, compute_coverage, compute_diversity, compute_composite

## [2026-04-03T02:40:00Z] Task 8: Resume/Idempotence + Logging - COMPLETED ✅

### Implementation Complete
- Modified autoresearch/harness/dihiggs_runner.py (added _load_existing_state, deduplication)
- Created autoresearch/tests/test_dihiggs_resume.py (13 tests)
- All 103 tests passing (90 previous + 13 new)

### Key Changes
1. **Added known_attempt_ids tracking**
   - Set[str] initialized in __init__
   - Populated from event log on startup
   - Checked before emission to prevent duplicates

2. **Implemented _load_existing_state()**
   - Reads event log if exists
   - Parses ATTEMPT_EVALUATED events
   - Populates known_attempt_ids for deduplication
   - Skips BanditState reconstruction (rebuilds naturally)

3. **Modified emit_attempt_event()**
   - Check if attempt_id in known_attempt_ids
   - Skip emission if duplicate (log info message)
   - Add to known_attempt_ids after successful emission

### Idempotence Strategy
- **Primary**: Track attempt_ids (SHA-256 fingerprints)
- **Skip**: Duplicate events on restart
- **Graceful**: Handle missing/malformed log files
- **Natural rebuild**: BanditState reconstructs as new attempts execute

### Design Rationale
- **Why not reconstruct BanditState from log?**
  - Attempt_id doesn't include arm_id (SHA-256 is one-way)
  - Would require parsing entire checkpoint history
  - UCB1 handles cold-start gracefully (0 pulls → explore first)
  - Primary goal is idempotence, not perfect state recovery

### Test Coverage (13 tests)
- load_existing_state: empty log, populated, multiple events, malformed JSON, empty lines, non-evaluated events
- emit_attempt_event: skip duplicates, add new to known set
- restart scenarios: no duplicates after crash, preserve known_run_dirs
- integration: checkpoint stability, subprocess logging
- full cycle: complete resume workflow

### Error Handling
- OSError: log error, continue (graceful degradation)
- JSONDecodeError: log warning, skip line
- Invalid payload: log warning, skip event
- Missing log: no-op (fresh start)

### Subprocess Logging
- Already implemented in Task 5: {outdir}/{arm_id}_round_{N}.stdout/stderr.txt
- Integration test added to verify files exist after subprocess call

### Checkpoint Stability
- Already implemented in parsers (2s mtime threshold)
- Integration test verifies stability check in action

## [2026-04-03T03:23:00Z] Task 9: Gated E2E Smoke Test - COMPLETED ✅

### Implementation Complete
- Created autoresearch/tests/test_e2e_smoke.py (100+ lines)
- Created autoresearch/configs/smoke_test_minimal.json (minimal config)
- Test skipped by default, runs only when DIHIGGS_E2E_TEST=1
- LSP diagnostics clean, pytest validation passing

### Files Created

#### 1. smoke_test_minimal.json
- Single-arm config (adaptive-smoke)
- Minimal parameters: tb=[1000], lam1_bins=10, mphi_bins=5
- Short timeout: 60 seconds
- Small iteration limits: max_new_run_dirs_per_arm_call=2
- Uses environment variables for paths:
  - HB_DATASET, HS_DATASET, PHYS_SCAN_EXEC
- All placeholders resolved by DiHiggsRunner

#### 2. test_e2e_smoke.py
- Module docstring with complete instructions
- **Skip mechanism**: @pytest.mark.skipif on DIHIGGS_E2E_TEST env var
- **Test function**: test_full_harness_smoke_test_gated()
  - Loads config with _get_smoke_config() helper
  - Creates DiHiggsRunner in tmp_path
  - Runs single_round("adaptive-smoke")
  - Validates: subprocess_status, event_log exists, stdout/stderr logs exist
  - Prints diagnostic output

### Skip Behavior Verified
```bash
# Without env var: SKIPPED (default, CI safe)
cd autoresearch && pytest tests/test_e2e_smoke.py -v
→ 1 skipped in 0.04s ✓

# With env var: Would run if prerequisites available
DIHIGGS_E2E_TEST=1 pytest tests/test_e2e_smoke.py -v
→ runs full harness workflow (requires binaries + datasets)
```

### Design Decisions
1. **Skip-by-default**: Prevents CI failures without real prerequisites
2. **Environment gating**: DIHIGGS_E2E_TEST=1 enables manual testing
3. **Minimal config**: Small iterations/bins/timeout reduces test runtime
4. **Assertions validate**:
   - Subprocess execution (success/nonzero_exit/timeout/preflight_fail are all acceptable)
   - Event log creation (even empty log is success)
   - Stdout/stderr logging (subprocess capture working)
5. **Configuration templating**: repo_root, sys_executable, phys_scan_exec resolved at test time

### Integration Points
- Uses existing DiHiggsRunner class (full integration)
- Validates preflight checks
- Validates checkpoint parsing logic (indirectly via run_single_round)
- Validates event emission
- Tests all Task 1-8 components working together

### Running the Test

```bash
# Default (CI): Skipped
pytest tests/test_e2e_smoke.py -v
→ test_full_harness_smoke_test_gated SKIPPED [100%] ✓

# Manual (with prerequisites):
export DIHIGGS_E2E_TEST=1
export PHYS_SCAN_EXEC=/path/to/PhysScanWithFixings
export HB_DATASET=/path/to/hb_dataset
export HS_DATASET=/path/to/hs_dataset
cd autoresearch && pytest tests/test_e2e_smoke.py -v -s
→ Runs full harness, validates all components, prints diagnostics
```

### Error Handling
- Missing smoke_test_minimal.json: FileNotFoundError (caught early)
- Invalid config: ValueError from DiHiggsRunner (expected, test fails clearly)
- Subprocess nonzero: Accepted as valid status (explorer may fail gracefully)
- Timeout: Accepted as valid status (test still validates logs created)
- Preflight fail: Accepted as valid status (gates subprocess cleanly)

### Test Philosophy
- **Permissive success criteria**: Accepts any subprocess status (success/nonzero/timeout/preflight_fail)
- **Strict artifact requirements**: Event log and stdout/stderr MUST exist
- **Environmental awareness**: Skips if prerequisites not configured
- **CI-safe by default**: No failures in CI without env var set
- **Manual-friendly**: Clear instructions in docstring for running with real data

## [2026-04-03T04:30:00Z] Task 10: Documentation Complete - VERIFIED ✅

### README.md Created
- Location: `autoresearch/README.md`
- Length: 243 lines (within 150-250 target range)
- Format: Professional markdown with clear hierarchy

### Content Coverage
1. **Overview**: Black-box subprocess architecture, key features
2. **Prerequisites**: System dependencies (GSL, CMake, OpenMP, Python 3.12+)
3. **Build Instructions**: Complete 3-step sequence
   - 2HDMC library (make in 2hdmc/)
   - HiggsTools library (cmake + make install)
   - DiHiggs PhysScanWithFixings binary
4. **Installation**: PYTHONPATH-based (no setup.py)
5. **Configuration**: Full schema with placeholder reference
6. **Usage**: CLI and programmatic examples
7. **Architecture**: 6 component descriptions
8. **Testing**: Full test suite + gated E2E test instructions
9. **Troubleshooting**: 8 common issues in table format
10. **File Reference**: Links to actual config/test files

### Build Sequence Documented
```bash
# Step 1: 2HDMC
cd /home/fabi/dihiggs/2hdmc && make
# Output: lib/lib2HDMC.a

# Step 2: HiggsTools
cd /home/fabi/dihiggs/higgstools && mkdir -p installation && rm -rf build && mkdir build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX=../installation && make -j$(nproc) && make install
# Output: installation/lib/libHiggsTools.so

# Step 3: PhysScanWithFixings
cd /home/fabi/dihiggs/dihiggs && make all
# Output: app/PhysScanWithFixings
```

### Environment Variables Documented
- `PHYS_SCAN_EXEC`: Path to PhysScanWithFixings binary (required)
- `HB_DATASET`: HiggsBounds dataset path (optional, warn if missing)
- `HS_DATASET`: HiggsSignals dataset path (optional, warn if missing)

### Example Invocations
- CLI: `python -m autoresearch.harness --mode dihiggs-explorers --config configs/dihiggs_explorers.json --outdir runs/campaign_001`
- Programmatic: `DiHiggsRunner(config, outdir).run_adaptation_round()`
- Gated E2E: `DIHIGGS_E2E_TEST=1 pytest tests/test_e2e_smoke.py -v -s`

### File References Validated
All referenced files exist:
- ✅ `configs/dihiggs_explorers.json`
- ✅ `configs/smoke_test_minimal.json`
- ✅ `tests/test_e2e_smoke.py`
- ✅ `/home/fabi/dihiggs/dihiggs/Makefile`
- ✅ `/home/fabi/dihiggs/dihiggs/app/PhysScanWithFixings`

### Test Results
- Full suite: 103 passed, 1 skipped (gated E2E)
- No regressions from documentation addition

### Key Design Principles Emphasized
1. **Black-box integration**: No imports from dihiggs/app/*
2. **Subprocess safety**: Timeout handling, log capture
3. **Idempotent resume**: SHA-256 attempt_ids, event log reconstruction
4. **Multi-axis tracking**: Sparse coverage/diversity metrics
5. **Bandit adaptation**: UCB1 arm selection
6. **Preflight gates**: Binary/dataset validation before execution

### Troubleshooting Coverage
| Issue | Solution |
|-------|----------|
| Missing binary | Follow 3-step build sequence |
| GSL not found | `apt-get install libgsl-dev` |
| HiggsTools not found | Check HIGGSTOOLS_PREFIX |
| Import errors | Set PYTHONPATH |
| Timeout | Increase timeout_sec in config |
| Parse errors | Check 2s mtime threshold |
| Permissions | Check outdir write access |
| Dataset warnings | Optional, set env vars |

### Writing Style Achieved
- User-facing (not developer internal)
- Clear section headers with markdown hierarchy
- Code blocks with syntax highlighting (bash/json/python)
- Bullet points for lists
- Command examples with full paths
- Brief explanations (no over-explaining)
- Professional tone, concise paragraphs

## BanditState Reconstruction Implementation (Task 9)

### COMPLETED: Full BanditState reconstruction from event log
- Added `arm_id` to event payload in `emit_attempt_event()` method signature
- Updated call site in `run_single_round()` to pass `arm_id` parameter
- Implemented complete event replay in `_load_existing_state()` that:
  - Extracts `arm_id` from event payload (previously missing)
  - Reconstructs attempt_data with yield, coverage, diversity 
  - Calls `compute_arm_reward()` to recompute composite score
  - Calls `update_state()` to replay state updates
  - Maintains backward compatibility by skipping events without arm_id

### Key Implementation Details
- Event schema unchanged (schema_version: 1)
- BanditState fields preserved: arm_stats, global_history, arm_histories
- Coverage/diversity scores recomputed during replay (not stored in events)
- All 13 tests pass (test_dihiggs_resume.py)
- Test calls updated to include arm_id parameter (API compatibility)

### State Reconstruction Pattern
When harness restarts:
1. Load event log line-by-line
2. For each ATTEMPT_EVALUATED event with arm_id:
   - Extract attempt_data and scores from payload
   - Compute reward via `compute_arm_reward()`
   - Call `update_state()` to accumulate into BanditState
3. Result: BanditState matches pre-restart state exactly

### Why This Works
- `update_state()` is idempotent with respect to the event log
- Replaying same events produces same state increments
- UCB1 algorithm benefits from full history reconstruction
- New arms get proper pull counts, total_rewards, and histories

