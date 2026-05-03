# Autoresearch UX Improvements Summary

## What Changed

Three major UX improvements have been added to make autoresearch more user-friendly and autonomous:

### 1. ✅ Interactive Config Builder
### 2. ✅ Time-Budget Continuous Exploration
### 3. ✅ Results Aggregation and Visualization

---

## 1. Interactive Config Builder

**Problem**: Manually editing JSON configs is tedious and error-prone.

**Solution**: CLI wizard for generating campaign configurations interactively.

### Quick Start

```bash
# Launch interactive wizard
python autoresearch/harness/interactive_config.py --wizard

# Use a template as starting point
python autoresearch/harness/interactive_config.py --template lambda1

# Edit existing config
python autoresearch/harness/interactive_config.py --edit configs/my_config.json

# List available templates
python autoresearch/harness/interactive_config.py --list-templates
```

### Features

- **Templates**: Pre-configured starting points (lambda1 exploration, m_phi scans, etc.)
- **Guided prompts**: Step-by-step configuration with sensible defaults
- **Physics parameters**: Interactive editing of λ₁ ranges, tan(β), m_φ, etc.
- **Time budgets**: Built-in support for duration-based campaigns
- **Validation**: Automatic config structure generation

### Example Session

```
=== Campaign Settings ===

Max rounds (0 = run until convergence) [0]: 20
Enable convergence detection? (y/n) [y]: y
Use time budget (instead of/in addition to max_rounds)? (y/n) [n]: y
Maximum campaign duration (hours) [8.0]: 12.0
  ✓ Campaign will stop after 12.0 hours

=== Physics Parameters ===

Lambda1 minimum [0]: 0
Lambda1 maximum [12]: 12
Lambda1 bins (resolution) [2000]: 2000
tan(beta) [1000]: 1000
Proposals per round [100]: 200

✅ Configuration saved to: configs/my-campaign.json
```

---

## 2. Time-Budget Continuous Exploration

**Problem**: Campaign stops after fixed rounds, not truly autonomous.

**Solution**: Duration-based and deadline-based stopping conditions.

### Quick Start

**Duration-based (run for X hours)**:
```json
{
  "supervisor": {
    "max_rounds": 0,
    "max_duration_hours": 8.0,
    "enable_convergence": true
  }
}
```

**Deadline-based (stop at specific time)**:
```json
{
  "supervisor": {
    "max_rounds": 0,
    "stop_at_timestamp": 1735689600.0,
    "enable_convergence": true
  }
}
```

**Unlimited (run until convergence)**:
```json
{
  "supervisor": {
    "max_rounds": 0,
    "enable_convergence": true
  }
}
```

### Features

- **Time tracking**: Elapsed time tracked in `campaign_state.json`
- **Cumulative timing**: Across campaign restarts/resumes
- **Priority system**: Time budget > round limit > convergence
- **Monitoring integration**: `live_monitor.py` shows elapsed time and ETA

### Monitoring

```bash
# Check elapsed time
python autoresearch/harness/live_monitor.py runs/my-campaign

# Output:
# State: RUNNING
# Rounds: 14 / ∞
# Elapsed: 4.1 hours (of 8.0 hour budget, 51%)
# ETA: ~3.9 hours remaining
```

### Stop Reasons

Campaign can complete with different reasons:
- `time_budget_exceeded (max_duration_hours=8.0, elapsed=8.00)`
- `time_budget_exceeded (stop_at_timestamp=1735689600.0)`
- `round_limit_reached` (if `max_rounds` > 0)
- `converged` (if convergence detector signals completion)

---

## 3. Results Aggregation and Visualization

**Problem**: User can't see the 2000 physics results from completed campaigns.

**Solution**: Automated results aggregation, analysis, and visualization tool.

### Quick Start

```bash
# Generate summary + plots
python autoresearch/harness/aggregate_results.py runs/my-campaign --plots

# Save aggregated CSV + plots
python autoresearch/harness/aggregate_results.py runs/my-campaign --output results.csv --plots
```

### Features

**Summary Statistics**:
- Total evaluations, iterations, viable points
- Constraint pass rates (positivity, unitarity, perturbativity)
- Parameter coverage (min/max/mean for all parameters)
- Physics results (branching ratios, decay widths)

**Visualization Plots** (saved to `{campaign_dir}/analysis_plots/`):
1. **lambda1_coverage.png**: λ₁ exploration over time
2. **constraint_evolution.png**: Constraint satisfaction by iteration
3. **parameter_space.png**: λ₁ vs m_φ colored by viability
4. **lambda1_histogram.png**: λ₁ sampling distribution

**Aggregated CSV**:
- All 2000+ results in single file
- Full physics data (parameters, constraints, widths, branching ratios)
- Metadata (iteration, run_id, file paths)
- Ready for Jupyter notebooks, ML training, custom analysis

### Example Output

```
============================================================
CAMPAIGN RESULTS SUMMARY
============================================================

📊 Evaluation Statistics:
  Total evaluations: 2,000
  Unique iterations: 20
  Viable points (all constraints): 680

✓ Constraint Pass Rates:
  Positivity:      99.0%
  Unitarity:       100.0%
  Perturbativity:  35.0%
  ALL constraints: 34.0%

📈 Parameter Coverage:
  lam1        : [    0.000,    12.000]  (mean:     6.000)
  m_phi       : [  130.000,   130.000]  (mean:   130.000)
  tan_beta    : [ 1000.000,  1000.000]  (mean:  1000.000)

📊 All plots saved to: runs/my-campaign/analysis_plots/
✅ Analysis complete!
```

---

## Complete Workflow

### 1. Create Campaign Config (Interactive)

```bash
python autoresearch/harness/interactive_config.py --wizard
```

Follow prompts to configure:
- Max rounds or time budget
- Thread count and timeouts
- Physics parameters (λ₁ range, tan(β), etc.)
- Proposals per round

Saves to: `configs/my-campaign.json`

### 2. Run Campaign

```bash
python autoresearch/run_supervisor.py configs/my-campaign.json
```

**Autonomous behavior**:
- UCB1 bandit selects best explorer arms
- Bayesian budget allocation focuses on promising regions
- Autoscaling adjusts resources based on convergence
- Runs until time budget or convergence

### 3. Monitor Progress

```bash
# Real-time monitoring with time tracking
python autoresearch/harness/live_monitor.py --watch runs/my-campaign
```

Shows:
- Current round and state
- Elapsed time / time budget (if configured)
- Discoveries per round
- Convergence status

### 4. Analyze Results

```bash
# Aggregate and visualize
python autoresearch/harness/aggregate_results.py runs/my-campaign --plots --output results.csv
```

Generates:
- Printed summary statistics
- 4 visualization plots
- Aggregated CSV for further analysis

### 5. Iterate

Based on results:
- Identify promising regions
- Create refined campaign configs
- Launch follow-up explorations

---

## File Locations

### New Tools

- `autoresearch/harness/interactive_config.py` - Interactive config wizard
- `autoresearch/harness/aggregate_results.py` - Results analysis tool

### Modified Files

- `autoresearch/harness/campaign_supervisor.py` - Added time-budget support

### Documentation

- `autoresearch/TIME_BUDGET_GUIDE.md` - Complete time-budget feature guide
- `autoresearch/RESULTS_ANALYSIS_GUIDE.md` - Results visualization guide
- `autoresearch/UX_IMPROVEMENTS_SUMMARY.md` - This file

---

## Key Differences from Before

| Feature | Before | After |
|---------|--------|-------|
| **Config creation** | Manual JSON editing | Interactive CLI wizard |
| **Stopping condition** | Fixed rounds only | Rounds, time budget, or convergence |
| **Autonomy** | Stops after max_rounds | Runs continuously until budget/convergence |
| **Results visibility** | Hidden in 2000 directories | Aggregated summary + plots |
| **Monitoring** | Abstract metrics only | Real physics results visible |
| **UX** | Technical, manual | User-friendly, automated |

---

## Belief Updating (Already Implemented)

**Note**: The system already implements sophisticated belief updating:

1. **UCB1 Bandit Algorithm**: Balances exploration (trying under-sampled arms) with exploitation (favoring high-performing arms)
2. **Bayesian Budget Allocation**: Beta-binomial posterior for allocating points across bins
3. **Composite Reward**: Combines yield, coverage, and diversity for intelligent arm selection
4. **Autoscaling**: Dynamically adjusts resources based on convergence signals

**Time-budget mode works seamlessly with these existing mechanisms.** The campaign continues updating beliefs and refining exploration until the time limit is reached.

---

## Next Steps

1. **Try the interactive wizard**:
   ```bash
   python autoresearch/harness/interactive_config.py --wizard
   ```

2. **Run a time-budgeted campaign**:
   ```bash
   python autoresearch/run_supervisor.py configs/my-campaign.json
   ```

3. **Visualize your results**:
   ```bash
   python autoresearch/harness/aggregate_results.py runs/my-campaign --plots
   ```

4. **Explore the guides**:
   - `TIME_BUDGET_GUIDE.md` - Time-budget feature details
   - `RESULTS_ANALYSIS_GUIDE.md` - Results analysis workflows
   - `MONITORING_AND_CONTROL.md` - Campaign monitoring (already existed)

---

## Questions?

- **Interactive config**: See `TIME_BUDGET_GUIDE.md` section "Interactive Configuration"
- **Time budgets**: See `TIME_BUDGET_GUIDE.md` for full details
- **Results analysis**: See `RESULTS_ANALYSIS_GUIDE.md` for examples
- **Monitoring**: See `MONITORING_AND_CONTROL.md` for live monitoring

All tools include `--help` for usage information:
```bash
python autoresearch/harness/interactive_config.py --help
python autoresearch/harness/aggregate_results.py --help
python autoresearch/harness/live_monitor.py --help
```
