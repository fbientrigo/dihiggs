# Autoresearch DiHiggs Explorers Harness

A Python harness for orchestrating adaptive parameter space exploration using the DiHiggs explorers. Wraps `adaptive_explorer.py` and `branch_continuation_explorer_v2.py` as black-box subprocesses with checkpoint parsing, multi-axis metrics, and bandit-based arm selection.

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration Reference](#configuration-reference)
- [Usage Examples](#usage-examples)
- [CLI Reference](#cli-reference)
- [Architecture Deep-Dive](#architecture-deep-dive)
- [Metrics Explained](#metrics-explained)
- [Troubleshooting](#troubleshooting)
- [Lam1 Backend](#lam1-backend)
- [Development](#development)

## Overview

The autoresearch harness automates parameter space exploration for DiHiggs physics analysis. It treats explorer scripts as black-box subprocesses, parses their checkpoint files, computes multi-dimensional metrics, and uses the UCB1 bandit algorithm to dynamically select the best exploration strategy.

### Key Features

**Black-Box Integration**
Explorers run as isolated subprocesses. No imports from `dihiggs/app/`. The harness communicates through checkpoint files and command-line arguments.

**UCB1 Bandit Adaptation**
The Upper Confidence Bound (UCB1) algorithm balances exploration (trying under-sampled arms) with exploitation (favoring high-performing arms). Arms with zero pulls always get priority for cold-start handling.

**Multi-Axis Diversity**
Tracks exploration across multiple parameter axes (tan(beta), lambda1, m_phi) with weighted coverage and Shannon entropy-based diversity metrics.

**Idempotent Resume**
Campaigns can be safely restarted. The harness reconstructs full state from the event log, skipping duplicate attempts and continuing from where it left off.

**Structured Event Logging**
All attempts are logged as JSONL events with schema versioning, enabling downstream analysis and reproducibility.

### Architecture Diagram

```
                    +------------------+
                    |   Config JSON    |
                    +--------+---------+
                             |
                             v
+--------------------------------------------------+
|                    Preflight                      |
|  - Check PhysScanWithFixings binary              |
|  - Verify HB_DATASET / HS_DATASET env vars       |
|  - Validate 2HDMC library                        |
+--------+-----------------------------------------+
         |
         v
+--------------------------------------------------+
|              DiHiggsRunner (Orchestrator)         |
|                                                   |
|  +----------------+    +---------------------+   |
|  |  Arm Selection |<---|   BanditState       |   |
|  |  (UCB1)        |    |   (UCB1 stats)      |   |
|  +--------+-------+    +---------------------+   |
|           |                                       |
|           v                                       |
|  +----------------+    +---------------------+   |
|  |  Subprocess    |--->|  Checkpoint Parser  |   |
|  |  Execution     |    |  (adaptive/branch)  |   |
|  +--------+-------+    +----------+----------+   |
|           |                       |               |
|           v                       v               |
|  +----------------+    +---------------------+   |
|  |  Metrics       |--->|  Event Emitter      |   |
|  |  (score.py)    |    |  (JSONL append)     |   |
|  +--------+-------+    +----------+----------+   |
|           |                       |               |
+-----------|-----------------------|---------------+
            |                       |
            v                       v
    +---------------+      +------------------+
    |  run_dirs/    |      |  events.jsonl    |
    |  (PhysScan    |      |  (attempt log)   |
    |   outputs)    |      +------------------+
    +---------------+
```

## Prerequisites

### Python Version

Python 3.12 or higher is required. The harness uses modern type hints and dataclass features.

```bash
python --version  # Should show 3.12.x or higher
```

### System Dependencies

```bash
# Ubuntu/Debian
sudo apt-get install libgsl-dev cmake

# macOS
brew install gsl cmake
```

OpenMP (usually bundled with GCC) is also required for parallel execution.

### Required Binaries

The harness requires a compiled `PhysScanWithFixings` binary from the DiHiggs project.

### Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `PHYS_SCAN_EXEC` | Yes | Path to PhysScanWithFixings binary |
| `HB_DATASET` | No | Path to HiggsBounds dataset (warns if missing) |
| `HS_DATASET` | No | Path to HiggsSignals dataset (warns if missing) |

Set these in your shell:

```bash
export PHYS_SCAN_EXEC=/path/to/dihiggs/app/PhysScanWithFixings
export HB_DATASET=/path/to/higgsbounds/dataset
export HS_DATASET=/path/to/higgssignals/dataset
```

### 2HDMC Library

The 2HDMC library must be compiled before running explorers:

```bash
cd /path/to/2hdmc
make
```

Expected output: `lib/lib2HDMC.a` (or `.so`/`.dylib` on some platforms)

## Installation

### Step 1: Build Dependencies

Build the three required libraries in sequence:

#### 1. Build 2HDMC Library

```bash
cd /path/to/2hdmc
make
```

Output: `lib/lib2HDMC.a`

#### 2. Build HiggsTools Library

```bash
cd /path/to/higgstools
mkdir -p installation
rm -rf build && mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../installation
make -j$(nproc) && make install
```

Output: `installation/lib/libHiggsTools.so`

#### 3. Build DiHiggs PhysScanWithFixings

```bash
cd /path/to/dihiggs
make all
```

Output: `app/PhysScanWithFixings`

### Step 2: Configure Python Environment

No setup.py is required. Use PYTHONPATH:

```bash
export PYTHONPATH="/path/to/2hdmc:$PYTHONPATH"
```

Add to your `.bashrc` or `.zshrc` for persistence.

### Step 3: Verify Installation

```bash
cd /path/to/2hdmc
python -c "from autoresearch.harness.dihiggs_runner import DiHiggsRunner; print('Import OK')"
```

## Configuration Reference

Configuration is specified via JSON files. Two example configs are provided:

- `configs/dihiggs_explorers.json` - Full production configuration
- `configs/smoke_test_minimal.json` - Minimal test configuration

### Top-Level Schema

```json
{
  "campaign_id": "unique-campaign-name",
  "paths": { ... },
  "runtime": { ... },
  "dihiggs": { ... },
  "search": { ... },
  "limits": { ... },
  "metrics": { ... },
  "adaptation": { ... },
  "arms": [ ... ]
}
```

### Field Reference

#### `campaign_id`

Unique identifier for this campaign. Used in event logs and attempt IDs.

| Property | Type | Description |
|----------|------|-------------|
| `campaign_id` | string | Unique campaign name (e.g., "dihiggs-explorers-dev-2026-04-01") |

#### `paths`

Directory and file path configuration.

| Property | Type | Description |
|----------|------|-------------|
| `repo_root` | string | Root directory of the dihiggs repository |
| `outdir` | string | Output directory for campaign results |
| `lake_name` | string | Name of the event log file (default: "events.jsonl") |

#### `runtime`

Execution environment settings.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `python_exe` | string | "{sys_executable}" | Python interpreter path (placeholder) |
| `threads` | integer | 4 | Number of OpenMP threads for PhysScan |
| `timeout_sec` | integer | 3600 | Subprocess timeout in seconds |
| `max_empty_rounds` | integer | 3 | Consecutive empty rounds before stopping |

#### `dihiggs`

DiHiggs-specific configuration.

| Property | Type | Description |
|----------|------|-------------|
| `phys_exec` | string | Path to PhysScanWithFixings binary (supports `{repo_root}` placeholder) |
| `hb_dataset_env` | string | Environment variable name for HiggsBounds dataset path |
| `hs_dataset_env` | string | Environment variable name for HiggsSignals dataset path |

#### `search`

Parameter space search configuration.

| Property | Type | Description |
|----------|------|-------------|
| `tb_values` | array[int] | List of tan(beta) values to explore (e.g., [1000, 5000, 10000, 30000]) |
| `lam1.min` | float | Minimum lambda1 value |
| `lam1.max` | float | Maximum lambda1 value |
| `lam1.n_bins` | integer | Number of bins for lambda1 discretization |
| `mphi.min` | float | Minimum m_phi value |
| `mphi.max` | float | Maximum m_phi value |
| `mphi.n_bins` | integer | Number of bins for m_phi discretization |

#### `limits`

Execution limits to prevent runaway exploration.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `max_new_run_dirs_per_round` | integer | 50 | Maximum new run directories per adaptation round |
| `max_new_run_dirs_per_arm_call` | integer | 20 | Maximum new run directories per single arm invocation |

#### `metrics`

Scoring and metric computation configuration.

**Floors** (minimum values for normalization):

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `floors.yield_norm` | float | 0.01 | Minimum yield for normalization |
| `floors.coverage_norm` | float | 0.0 | Minimum coverage for normalization |
| `floors.diversity_norm` | float | 0.0 | Minimum diversity for normalization |

**Weights** (for composite score):

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `weights.yield` | float | 0.4 | Weight for yield component |
| `weights.coverage` | float | 0.3 | Weight for coverage component |
| `weights.diversity` | float | 0.3 | Weight for diversity component |

**Multi-axis configuration**:

| Property | Type | Description |
|----------|------|-------------|
| `multi_axis.enabled` | boolean | Enable multi-axis tracking |
| `multi_axis.collapse_axes` | array[string] | Axes to collapse for cell ID generation |
| `multi_axis.coverage_axes` | array[object] | Per-axis coverage configuration |
| `multi_axis.diversity_axes` | array[object] | Per-axis diversity configuration |
| `multi_axis.diversity_pairs` | array[object] | Pairwise diversity interactions |

Coverage axis specification:

```json
{
  "name": "tb",
  "kind": "categorical",
  "domain": [1000, 5000, 10000, 30000],
  "weight": 0.3
}
```

Or for discrete axes:

```json
{
  "name": "lam1_bin",
  "kind": "discrete",
  "domain_size": 40,
  "weight": 0.4
}
```

Diversity axis specification:

```json
{
  "name": "tb",
  "weight": 0.2
}
```

Diversity pair specification (for joint entropy):

```json
{
  "axes": ["tb", "lam1_bin"],
  "weight": 0.5
}
```

#### `adaptation`

UCB1 bandit algorithm configuration.

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `ucb_c` | float | 1.414 | UCB1 exploration constant (sqrt(2) is common) |
| `warm_start_each_arm` | boolean | true | Whether to warm-start arms with previous state |

The UCB1 formula is:

```
score = mean_reward + c * sqrt(log(total_pulls) / arm_pulls)
```

Higher `ucb_c` values favor exploration; lower values favor exploitation.

#### `arms`

Array of explorer configurations. Each arm represents a different exploration strategy.

Arm specification:

| Property | Type | Required | Description |
|----------|------|----------|-------------|
| `id` | string | Yes | Unique arm identifier |
| `explorer` | string | Yes | Explorer type: "adaptive" or "branch" |
| `timeout_sec` | integer | No | Arm-specific timeout (overrides runtime.timeout_sec) |
| `cmd` | array[string] | Yes | Command template with placeholders |
| `env` | object | No | Environment variables with placeholders |

**Adaptive explorer example**:

```json
{
  "id": "adaptive-v1",
  "explorer": "adaptive",
  "timeout_sec": 7200,
  "cmd": [
    "{python}",
    "{repo_root}/dihiggs/app/adaptive_explorer.py",
    "--checkpoint-root", "{checkpoint_root}/adaptive-v1",
    "--tb-values", "{tb_values}",
    "--lam1-range", "{lam1_min}", "{lam1_max}",
    "--n-bins", "{lam1_n_bins}",
    "--n-proposals", "10",
    "--iter", "{iter}"
  ],
  "env": {
    "HB_DATASET": "{hb_dataset}",
    "HS_DATASET": "{hs_dataset}",
    "OMP_NUM_THREADS": "{threads}"
  }
}
```

**Branch explorer example**:

```json
{
  "id": "branch-v1",
  "explorer": "branch",
  "timeout_sec": 7200,
  "cmd": [
    "{python}",
    "{repo_root}/dihiggs/app/branch_continuation_explorer_v2.py",
    "--checkpoint-root", "{checkpoint_root}/branch-v1",
    "--track-id", "{track_id}",
    "--tb", "{tb}",
    "--n-steps", "5"
  ],
  "env": {
    "HB_DATASET": "{hb_dataset}",
    "HS_DATASET": "{hs_dataset}",
    "OMP_NUM_THREADS": "{threads}"
  }
}
```

### Placeholders

The following placeholders are available in `cmd` and `env` fields:

| Placeholder | Description |
|-------------|-------------|
| `{sys_executable}` | Path to current Python interpreter |
| `{python}` | Alias for `{sys_executable}` |
| `{repo_root}` | Repository root path from config |
| `{outdir}` | Campaign output directory |
| `{checkpoint_root}` | Checkpoint directory (outdir/checkpoints) |
| `{campaign}` | Campaign ID |
| `{threads}` | Thread count from runtime config |
| `{tb_values}` | Comma-separated tan(beta) values |
| `{tb}` | First tan(beta) value |
| `{track_id}` | Generated track ID for branch explorer |
| `{iter}` | Current iteration number for this arm |
| `{lam1_min}` | Lambda1 minimum from search config |
| `{lam1_max}` | Lambda1 maximum from search config |
| `{lam1_n_bins}` | Lambda1 bin count |
| `{hb_dataset}` | HiggsBounds dataset path from env |
| `{hs_dataset}` | HiggsSignals dataset path from env |

## Usage Examples

### Example 1: Basic Run (Single Arm)

Run a single iteration with one adaptive explorer arm:

```python
from autoresearch.harness.dihiggs_runner import DiHiggsRunner
import json

# Load configuration
with open("configs/smoke_test_minimal.json") as f:
    config = json.load(f)

# Create runner
runner = DiHiggsRunner(config=config, outdir="runs/basic_test")

# Run single round with specific arm
result = runner.run_single_round("adaptive-smoke")

print(f"Discoveries: {result['discoveries']}")
print(f"Events emitted: {result['events_emitted']}")
print(f"Status: {result['subprocess_status']}")
```

### Example 2: Multi-Arm Exploration (UCB1 Selection)

Let the UCB1 bandit algorithm select between multiple arms:

```python
from autoresearch.harness.dihiggs_runner import DiHiggsRunner
import json

with open("configs/dihiggs_explorers.json") as f:
    config = json.load(f)

runner = DiHiggsRunner(config=config, outdir="runs/ucb1_campaign")

# Run 10 adaptation rounds
for i in range(10):
    result = runner.run_adaptation_round()
    print(f"Round {i+1}: Selected {result['arm_id']}, "
          f"discoveries={result['discoveries']}")
```

The UCB1 algorithm will:
1. First explore each arm once (cold-start with infinity score)
2. Then balance exploration of under-sampled arms with exploitation of high-reward arms
3. Log all attempts to `events.jsonl` for analysis

### Example 3: Resume After Crash

The harness is fully idempotent. Restarting a campaign safely continues from the last state:

```python
from autoresearch.harness.dihiggs_runner import DiHiggsRunner
import json

# Same config, same outdir as before
with open("configs/dihiggs_explorers.json") as f:
    config = json.load(f)

# This will:
# 1. Load existing events.jsonl
# 2. Reconstruct BanditState from ATTEMPT_EVALUATED events
# 3. Skip any attempts already in known_attempt_ids
runner = DiHiggsRunner(config=config, outdir="runs/ucb1_campaign")

# Continue from where we left off
for i in range(10):
    result = runner.run_adaptation_round()
    print(f"Resumed round: Selected {result['arm_id']}")
```

Key points for resume:
- Use the same `outdir` path
- The same `campaign_id` must be in the config
- Event log (`events.jsonl`) must be intact
- Checkpoints are preserved in `outdir/checkpoints/`

### Example 4: Custom Scoring Weights

Emphasize diversity over raw yield by adjusting metric weights:

```json
{
  "campaign_id": "diversity-focused-campaign",
  "metrics": {
    "weights": {
      "yield": 0.2,
      "coverage": 0.3,
      "diversity": 0.5
    }
  },
  "adaptation": {
    "ucb_c": 2.0
  }
}
```

With these weights:
- Diversity contributes 50% to the composite score
- The higher `ucb_c` (2.0 vs 1.414) increases exploration
- The bandit will favor arms that explore under-represented parameter regions

## CLI Reference

The harness is designed for programmatic use. Import and use the `DiHiggsRunner` class directly.

### Basic Usage Pattern

```python
from autoresearch.harness.dihiggs_runner import DiHiggsRunner
import json

# Load configuration
with open("path/to/config.json") as f:
    config = json.load(f)

# Initialize runner
runner = DiHiggsRunner(
    config=config,
    outdir="./results"
)

# Run specific arm
result = runner.run_single_round(arm_id="adaptive-v1")

# Or run with automatic arm selection
result = runner.run_adaptation_round()
```

### DiHiggsRunner API

#### Constructor

```python
DiHiggsRunner(config: dict[str, object], outdir: str)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `config` | dict | Full configuration dictionary |
| `outdir` | str | Output directory for campaign artifacts |

#### Methods

**run_single_round(arm_id, update_bandit=False)**

Execute a single round with a specific arm.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `arm_id` | str | required | Arm identifier from config |
| `update_bandit` | bool | False | Whether to update bandit state with results |

Returns: `dict` with keys:
- `discoveries`: Number of new run directories found
- `events_emitted`: Number of events written to log
- `preflight`: Preflight check results
- `subprocess_status`: "success", "timeout", "nonzero_exit", or "skipped_preflight_fail"

**run_adaptation_round()**

Automatically select and execute an arm using UCB1.

Returns: `dict` with keys:
- `arm_id`: Selected arm identifier
- `discoveries`: Number of new run directories
- `events_emitted`: Number of events written
- `preflight`: Preflight check results
- `subprocess_status`: Subprocess execution status

### Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `PHYS_SCAN_EXEC` | Yes | Path to PhysScanWithFixings binary |
| `HB_DATASET` | No | HiggsBounds dataset directory |
| `HS_DATASET` | No | HiggsSignals dataset directory |
| `PYTHONPATH` | Yes | Must include path to 2hdmc directory |

## Architecture Deep-Dive

### Component Overview

| Component | File | Responsibility |
|-----------|------|----------------|
| Preflight | `dihiggs_preflight.py` | Validates binaries, datasets, libraries before execution |
| Parsers | `dihiggs_parsers.py` | Reads adaptive JSON and branch JSONL checkpoints |
| Metrics | `benchmarks/score.py` | Computes yield, coverage, diversity, composite scores |
| Adaptation | `dihiggs_adaptation.py` | BanditState, UCB1 selection, reward computation |
| Axis Contract | `dihiggs_axis_contract.py` | Cell encoding/decoding, placeholder substitution |
| Runner | `dihiggs_runner.py` | Orchestrates the full pipeline |

### Data Flow

```
1. CONFIG LOADING
   Load JSON config with all arm definitions

2. PREFLIGHT CHECKS
   - Verify PhysScanWithFixings exists and is executable
   - Check HB_DATASET / HS_DATASET environment variables
   - Validate 2HDMC library presence

3. STATE RECONSTRUCTION (if resuming)
   - Read existing events.jsonl
   - Parse ATTEMPT_EVALUATED events
   - Rebuild known_attempt_ids set
   - Reconstruct BanditState from event history

4. ARM SELECTION (UCB1)
   - Compute UCB1 score for each arm
   - Cold-start arms (0 pulls) get infinity score
   - Select arm with highest score

5. SUBPROCESS EXECUTION
   - Substitute placeholders in cmd and env
   - Run explorer as subprocess with timeout
   - Capture stdout/stderr to files

6. CHECKPOINT PARSING
   - Parse adaptive_state.json or events.jsonl
   - Extract new run directories
   - Filter out already-known runs

7. METRICS COMPUTATION
   - Read results.csv from each run_dir
   - Compute yield (success rate)
   - Compute coverage (bins visited / total bins)
   - Compute diversity (Shannon entropy)
   - Compute composite (weighted sum)

8. EVENT EMISSION
   - Generate attempt_id (SHA-256 fingerprint)
   - Skip if attempt_id already known
   - Append ATTEMPT_EVALUATED event to events.jsonl

9. BANDIT UPDATE
   - Compute arm reward from attempt data
   - Update arm_stats (pulls, total_reward)
   - Append to global_history and arm_histories

10. ITERATION
    Return to step 4 for next round
```

### Event Log Schema

Events are appended to `events.jsonl` as JSON lines.

**Schema Version**: 1

**ATTEMPT_EVALUATED Event**:

```json
{
  "schema_version": 1,
  "campaign_id": "dihiggs-explorers-dev-2026-04-01",
  "event_type": "ATTEMPT_EVALUATED",
  "utc": "2026-04-01T12:34:56.789012Z",
  "payload": {
    "arm_id": "adaptive-v1",
    "attempt_id": "a1b2c3d4e5f67890",
    "iter_index": 0,
    "attempt_index": 0,
    "cell_id": "tb=1000|bin=15",
    "eval_status": "success",
    "successes": 42.0,
    "trials": 100,
    "elapsed_sec": 0.0,
    "axes_binned": {
      "tb": 1000,
      "lam1_bin": 15
    },
    "axes_raw": {
      "tb": 1000,
      "lam1": 3.14159
    }
  }
}
```

| Field | Description |
|-------|-------------|
| `schema_version` | Event schema version (currently 1) |
| `campaign_id` | Campaign identifier |
| `event_type` | Type of event (ATTEMPT_EVALUATED) |
| `utc` | ISO 8601 timestamp in UTC |
| `payload.arm_id` | Arm that generated this attempt |
| `payload.attempt_id` | Unique fingerprint (SHA-256, first 16 chars) |
| `payload.cell_id` | Encoded cell identifier (tb|bin format) |
| `payload.eval_status` | Evaluation status (success, failed, etc.) |
| `payload.successes` | Number of successful events |
| `payload.trials` | Total number of trials |
| `payload.axes_binned` | Discretized axis values |
| `payload.axes_raw` | Original (non-binned) axis values |

### BanditState Reconstruction

On startup, the harness reconstructs `BanditState` from the event log:

```python
# Pseudocode of reconstruction logic
for event in events_jsonl:
    if event.event_type == "ATTEMPT_EVALUATED":
        # Add to known attempts (deduplication)
        known_attempt_ids.add(event.payload.attempt_id)
        
        # Reconstruct attempt data
        attempt_data = {
            "cell_id": event.payload.cell_id,
            "axes_binned": event.payload.axes_binned,
            "axes_raw": event.payload.axes_raw,
            "yield": event.payload.successes,
        }
        
        # Compute reward and update state
        reward = compute_arm_reward([attempt_data], config)
        update_state(bandit_state, event.payload.arm_id, reward, attempt_data)
```

This ensures:
- No duplicate events are emitted
- Bandit statistics are accurate after restart
- Global history is preserved for coverage/diversity

### Attempt ID Generation

Attempt IDs are stable fingerprints computed as:

```python
fingerprint = f"{campaign_id}|{arm_id}|{cell_id}|{run_dir}"
attempt_id = hashlib.sha256(fingerprint.encode("utf-8")).hexdigest()[:16]
```

This guarantees idempotence: the same campaign, arm, cell, and run directory always produce the same attempt ID.

## Metrics Explained

### Yield

**Definition**: Raw success count from PhysScan execution.

**Computation**: Read `results.csv` from the run directory, count rows where `success=True`, divide by total rows.

**Range**: [0.0, 1.0] (0% to 100% success rate)

**Interpretation**: Higher yield means the parameter point produced more valid physics events. This is the primary objective but can lead to exploitation of known-good regions.

### Coverage

**Definition**: Fraction of axis bins that have been explored.

**Computation**:
```
per_axis_coverage = bins_visited / total_bins
overall_coverage = weighted_mean(per_axis_coverages)
```

**Range**: [0.0, 1.0]

**Interpretation**: Higher coverage means broader exploration of the parameter space. Prevents getting stuck in local optima.

### Diversity

**Definition**: Entropy-based measure of distribution uniformity across axes.

**Computation**: Shannon entropy for each axis:
```
H = -sum(p * log2(p)) for all p > 0
where p = frequency(value) / total_attempts
```

Pairwise interactions can also be configured for joint entropy.

**Range**: [0.0, unbounded] (higher = more diverse)

**Interpretation**: Higher diversity means attempts are spread evenly across axis values. Prevents clustering in popular regions.

### Composite Score

**Definition**: Weighted combination of yield, coverage, and diversity.

**Formula**:
```
composite = w_yield * yield + w_coverage * coverage + w_diversity * diversity
```

Default weights: yield=0.4, coverage=0.3, diversity=0.3

**Usage**: This is the reward signal for the UCB1 bandit algorithm. Arms are selected to maximize expected composite score.

### Metric Configuration Examples

**Yield-focused** (exploit known-good regions):
```json
{
  "metrics": {
    "weights": {
      "yield": 0.7,
      "coverage": 0.2,
      "diversity": 0.1
    }
  }
}
```

**Exploration-focused** (discover new regions):
```json
{
  "metrics": {
    "weights": {
      "yield": 0.2,
      "coverage": 0.4,
      "diversity": 0.4
    }
  }
}
```

**Balanced** (default):
```json
{
  "metrics": {
    "weights": {
      "yield": 0.4,
      "coverage": 0.3,
      "diversity": 0.3
    }
  }
}
```

## Troubleshooting

### Preflight Failures

**Issue**: `PhysScanWithFixings not found`

**Solution**:
```bash
# Verify binary exists
ls -la $PHYS_SCAN_EXEC

# If missing, rebuild
cd /path/to/dihiggs
make all

# Set environment variable
export PHYS_SCAN_EXEC=/path/to/dihiggs/app/PhysScanWithFixings
```

**Issue**: `Dataset env var(s) not set`

**Solution**:
```bash
# Optional but recommended
export HB_DATASET=/path/to/higgsbounds/dataset
export HS_DATASET=/path/to/higgssignals/dataset

# Or disable warnings by setting to empty
export HB_DATASET=""
```

**Issue**: `2HDMC library not found`

**Solution**:
```bash
cd /path/to/2hdmc
make clean && make
ls -la lib/lib2HDMC.*
```

### Checkpoint Parsing Errors

**Issue**: `Skipping unstable checkpoint (recently modified)`

**Cause**: Checkpoint file was modified within the last 2 seconds.

**Solution**: Wait a moment and retry. The stability threshold prevents reading half-written files.

**Issue**: `Failed to load checkpoint ... JSONDecodeError`

**Cause**: Corrupted or incomplete checkpoint file.

**Solution**:
```bash
# Check checkpoint integrity
cat /path/to/checkpoint/adaptive_state.json | python -m json.tool

# If corrupted, remove and restart (explorers will regenerate)
rm -rf /path/to/checkpoint/iter_*/adaptive_state.json
```

### Timeout Handling

**Issue**: `Timeout while executing arm ...`

**Cause**: Explorer subprocess exceeded `timeout_sec`.

**Solutions**:
1. Increase timeout in config:
   ```json
   {
     "runtime": {
       "timeout_sec": 7200
     }
   }
   ```
2. Reduce work per arm call:
   ```json
   {
     "limits": {
       "max_new_run_dirs_per_arm_call": 10
     }
   }
   ```
3. Check if PhysScan is hanging (check stdout/stderr files in outdir)

### Resume Not Working

**Issue**: Campaign restarts from beginning instead of resuming.

**Checklist**:
1. Verify `outdir` path is identical to previous run
2. Check that `events.jsonl` exists and is readable:
   ```bash
   ls -la $OUTDIR/events.jsonl
   wc -l $OUTDIR/events.jsonl
   ```
3. Verify `campaign_id` in config matches previous run
4. Check for JSON parse errors in event log:
   ```bash
   python -c "
   import json
   with open('events.jsonl') as f:
       for i, line in enumerate(f, 1):
           try:
               json.loads(line)
           except json.JSONDecodeError as e:
               print(f'Line {i}: {e}')
   "
   ```

### Low Discovery Rate

**Issue**: Arms run but produce zero discoveries.

**Diagnostics**:
```bash
# Check explorer stdout/stderr
cat runs/campaign_001/adaptive-v1_round_0.stdout.txt
cat runs/campaign_001/adaptive-v1_round_0.stderr.txt

# Verify PhysScan works standalone
$PHYS_SCAN_EXEC --help

# Check preflight status
python -c "
from autoresearch.harness.dihiggs_preflight import run_all_preflight_checks
import json
with open('configs/dihiggs_explorers.json') as f:
    config = json.load(f)
result = run_all_preflight_checks(config)
print(json.dumps(result, indent=2))
"
```

### Import Errors

**Issue**: `ModuleNotFoundError: No module named 'autoresearch'`

**Solution**:
```bash
export PYTHONPATH="/path/to/2hdmc:$PYTHONPATH"

# Verify
python -c "from autoresearch.harness.dihiggs_runner import DiHiggsRunner; print('OK')"
```

## Lam1 Backend

### Overview

The lam1 backend uses the `PhysLam1Scan` executable (instead of `PhysScanWithFixings`) to perform single-parameter scans across λ₁ values. This backend is designed for explorers targeting lambda1 parameter exploration.

Configuration: `configs/dihiggs_explorers_lam1.json`

### Configuration Selection

To use the lam1 backend, load the lam1-specific configuration:

```python
from autoresearch.harness.dihiggs_runner import DiHiggsRunner
import json

with open('configs/dihiggs_explorers_lam1.json') as f:
    config = json.load(f)

runner = DiHiggsRunner(config=config, outdir='runs/lam1_campaign')
```

Key differences from the base configuration:
- Uses `{repo_root}/dihiggs/app/PhysLam1Scan` as the physical scan executable
- Lambda1 scans range from -20.0 to 20.0 with 40 bins
- m_phi is fixed at 130 GeV (single bin)
- Four tan(beta) values: 1000, 5000, 10000, 30000
- Single adaptive arm: `adaptive-lam1-v1`

### Expected Artifacts

After running the lam1 backend, the campaign produces several artifacts:

#### 1. Event Log

**File**: `{outdir}/events.jsonl`

Structured event log containing all attempt records. Each line is a JSON object with schema version 1.

#### 2. Checkpoints

**Directory**: `{outdir}/checkpoints/`

Contains subdirectories for each explorer arm (e.g., `checkpoints/adaptive-lam1-v1/`).

**State File**: `{outdir}/checkpoints/adaptive-lam1-v1/iter_*/adaptive_state.json`

Stores the current state of the adaptive explorer, including:
- Current iteration index
- Best discovered parameters
- Coverage of the λ₁ parameter space
- Historical yield and constraint satisfaction metrics

#### 3. Results CSV

**File**: `{outdir}/run_dirs/*/results.csv`

The `PhysLam1Scan` executable outputs results.csv for each run directory. Schema matches the base `PhysScanWithFixings` output with 29 columns (see main README for full schema).

Key columns for lam1 backend:
- `m_phi`: Fixed at 130.0 GeV
- `lam1`: The scanned λ₁ value
- `computed_lam1`: Verification that computed λ₁ matches input
- `positivity_ok`, `unitarity_ok`, `perturbativity_ok`: Physics constraint satisfaction
- `total_width`: Total decay width (GeV)
- `br_gaga`: Branching ratio to photons

### Resume Behavior

The lam1 backend supports safe campaign resumption. To resume a previous campaign:

```python
from autoresearch.harness.dihiggs_runner import DiHiggsRunner
import json

# Use SAME outdir and campaign_id
with open('configs/dihiggs_explorers_lam1.json') as f:
    config = json.load(f)

runner = DiHiggsRunner(config=config, outdir='runs/lam1_campaign')

# Continue from last state
for i in range(10):
    result = runner.run_adaptation_round()
    print(f'Round {i+1}: {result["arm_id"]}, discoveries={result["discoveries"]}')
```

Resume guarantees:
- All previous ATTEMPT_EVALUATED events are reconstructed from `events.jsonl`
- Duplicate attempts are skipped (based on SHA-256 fingerprint)
- BanditState is fully recovered, including arm statistics
- Coverage and diversity metrics reflect full history
- Campaign continues from last round number

### Troubleshooting

#### PhysLam1Scan Not Found

**Error**: `PhysLam1Scan executable not found or not executable`

**Solution**:
```bash
# Verify binary exists
ls -la $PHYS_SCAN_EXEC  # or check {repo_root}/dihiggs/app/PhysLam1Scan

# If missing, rebuild dihiggs
cd /path/to/dihiggs
make all

# Set environment variable if using custom path
export PHYS_SCAN_EXEC=/path/to/custom/PhysLam1Scan
```

#### Missing HiggsBounds/HiggsSignals Datasets

**Error**: `HB_DATASET environment variable not set` or `HS_DATASET environment variable not set`

**Impact**: Preflight warning but execution continues if binary is present.

**Solution**:
```bash
# Set dataset paths (optional but recommended)
export HB_DATASET=/path/to/higgsbounds/dataset
export HS_DATASET=/path/to/higgssignals/dataset
```

#### OpenMP Thread Safety

**Error**: Segmentation fault or hang during parallel execution

**Cause**: PhysLam1Scan is single-threaded and may not be thread-safe.

**Solution**:
```bash
# Always set OMP_NUM_THREADS=1
export OMP_NUM_THREADS=1
```

This is already set in the config `env` section with `"OMP_NUM_THREADS": "{threads}"` (defaults to 4 from runtime config). For safety, the harness always enforces `OMP_NUM_THREADS=1` when executing PhysLam1Scan.

#### Low or Zero Discoveries

**Symptom**: Arm runs but finds no new parameter points

**Diagnostics**:
```bash
# Check adapter output
cat {outdir}/adaptive-lam1-v1_round_0.stdout.txt
cat {outdir}/adaptive-lam1-v1_round_0.stderr.txt

# Verify PhysLam1Scan works standalone
{repo_root}/dihiggs/app/PhysLam1Scan 130 130 1 -20 -20 1 300 0.999 50 0.1 0.0 /tmp/test.csv

# Check event log for errors
tail -n 5 {outdir}/events.jsonl | python -m json.tool
```

**Common causes**:
- Adaptive explorer exhausted λ₁ range in previous runs (check coverage in adaptive_state.json)
- Parameter points outside valid physics region (check positivity, unitarity, perturbativity flags)
- Mismatch between config and actual parameter ranges
## Development

### Running Tests

The test suite includes unit tests for all components:

```bash
cd /path/to/2hdmc/autoresearch
pytest tests/ -v
```

Test coverage:
- `test_dihiggs_runner.py` - Runner orchestration
- `test_dihiggs_adaptation.py` - UCB1 bandit logic
- `test_dihiggs_preflight.py` - Preflight checks
- `test_dihiggs_parsers_adaptive.py` - Adaptive checkpoint parsing
- `test_dihiggs_parsers_branch.py` - Branch checkpoint parsing
- `test_dihiggs_axis_contract.py` - Cell encoding/decoding
- `test_score.py` - Metrics computation
- `test_dihiggs_resume.py` - State reconstruction
- `test_e2e_smoke.py` - End-to-end smoke test (gated)

**Gated E2E Test** (requires real binaries):

```bash
export DIHIGGS_E2E_TEST=1
export PHYS_SCAN_EXEC=/path/to/dihiggs/app/PhysScanWithFixings
pytest tests/test_e2e_smoke.py -v -s
```

### Adding New Explorers

To add a new explorer type:

1. **Create parser** in `harness/dihiggs_parsers.py`:
   ```python
   @dataclass
   class MyExplorerDiscovery:
       run_dir: str
       tb: int
       lam1_raw: float
       lam1_bin: int
       # ... additional fields
   
   def parse_my_explorer_checkpoint(
       checkpoint_root: str, config: Mapping[str, object], known_run_dirs: set[str]
   ) -> list[MyExplorerDiscovery]:
       # Implementation
       pass
   ```

2. **Register in runner** in `harness/dihiggs_runner.py`:
   ```python
   elif explorer == "my_explorer":
       discoveries = parse_my_explorer_checkpoint(checkpoint_root, self.config, self.known_run_dirs)
   ```

3. **Add arm config** in your JSON config:
   ```json
   {
     "id": "my-explorer-v1",
     "explorer": "my_explorer",
     "cmd": ["{python}", "{repo_root}/dihiggs/app/my_explorer.py", ...]
   }
   ```

### Extending the Axis Contract

To add new axes for tracking:

1. **Update AxisContract** in `harness/dihiggs_axis_contract.py`:
   ```python
   @dataclass(frozen=True)
   class AxisContract:
       physical_axes: tuple[str, ...] = ("tb", "lam1_bin", "mphi_bin", "my_axis")
   ```

2. **Add binning function**:
   ```python
   def bin_my_axis(raw_value: float, config: Mapping[str, object]) -> int:
       search = _expect_mapping(config["search"], "search")
       my_axis_cfg = _expect_mapping(search["my_axis"], "search.my_axis")
       return _linear_bin(float(raw_value), my_axis_cfg, "my_axis")
   ```

3. **Update config schema**:
   ```json
   {
     "search": {
       "my_axis": {
         "min": 0.0,
         "max": 100.0,
         "n_bins": 20
       }
     }
   }
   ```

4. **Update coverage/diversity axes** in metrics configuration.

### File Reference

| File | Purpose |
|------|---------|
| `configs/dihiggs_explorers.json` | Production configuration template |
| `configs/smoke_test_minimal.json` | Minimal test configuration |
| `harness/dihiggs_runner.py` | Main orchestration logic |
| `harness/dihiggs_preflight.py` | Pre-execution validation |
| `harness/dihiggs_parsers.py` | Checkpoint file parsing |
| `harness/dihiggs_adaptation.py` | UCB1 bandit algorithm |
| `harness/dihiggs_axis_contract.py` | Cell encoding and placeholders |
| `benchmarks/score.py` | Metrics computation |
| `tests/test_*.py` | Unit and integration tests |
