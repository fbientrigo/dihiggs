# Autoresearch Results Analysis Guide

## Quick Start

After completing a campaign, visualize and analyze your physics results:

```bash
# Aggregate all results and generate plots
python autoresearch/harness/aggregate_results.py runs/my-campaign --plots

# Save aggregated CSV for further analysis
python autoresearch/harness/aggregate_results.py runs/my-campaign --output results.csv --plots
```

## What You Get

### 1. Summary Statistics

The tool prints a comprehensive summary:

- **Evaluation Statistics**: Total points explored, unique iterations, viable points
- **Constraint Pass Rates**: Positivity, unitarity, perturbativity, and combined success rates
- **Parameter Coverage**: Min/max/mean for all explored parameters (λ₁, m_φ, tan(β), etc.)
- **Physics Results**: Branching ratios, decay widths, and other observables

Example output:

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

============================================================
```

### 2. Visualization Plots

Four plots are automatically generated in `{campaign_dir}/analysis_plots/`:

#### `lambda1_coverage.png`
Shows how λ₁ values were explored over time. Points are scattered by evaluation index.

**Use case**: Verify uniform coverage across the λ₁ range, identify exploration patterns.

#### `constraint_evolution.png`
Tracks constraint satisfaction rates over iterations.

**Use case**: Check if the adaptive algorithm is learning to avoid invalid regions, monitor exploration efficiency.

#### `parameter_space.png`
2D scatter plot of λ₁ vs m_φ, color-coded by constraint satisfaction (green=viable, red=invalid).

**Use case**: Visualize which regions of parameter space are viable, identify interesting boundaries.

#### `lambda1_histogram.png`
Distribution of λ₁ values explored.

**Use case**: Verify sampling strategy, check for clustering or gaps in coverage.

### 3. Aggregated CSV

The `--output` flag saves all 2000+ results into a single CSV file:

```bash
python autoresearch/harness/aggregate_results.py runs/my-campaign --output all_results.csv
```

**Columns include**:
- All physics parameters (m_phi, lambda1-7, tan_beta, sin_ba, etc.)
- Constraint flags (positivity_ok, unitarity_ok, perturbativity_ok)
- Decay widths (width_bb, width_tautau, width_WW, width_ZZ, width_gaga, etc.)
- Branching ratios (br_gaga, etc.)
- Metadata (iteration, run_id, file_path)

**Use cases**:
- Import into Jupyter notebooks for custom analysis
- Filter viable points for further study
- Train ML models on the results
- Identify promising regions for follow-up campaigns

## Advanced Analysis

### Filter Viable Points

```python
import pandas as pd

df = pd.read_csv("all_results.csv")

# Filter only viable points
viable = df[
    (df['positivity_ok'] == 1) &
    (df['unitarity_ok'] == 1) &
    (df['perturbativity_ok'] == 1)
]

print(f"Found {len(viable)} viable points out of {len(df)}")

# Find highest branching ratio
best = viable.loc[viable['br_gaga'].idxmax()]
print(f"Best BR(h→γγ): {best['br_gaga']:.6f} at λ₁={best['lam1']:.3f}")
```

### Analyze Exploration Efficiency

```python
import pandas as pd

df = pd.read_csv("all_results.csv")

# Group by iteration and compute success rate
by_iter = df.groupby('iteration').agg({
    'positivity_ok': 'mean',
    'unitarity_ok': 'mean',
    'perturbativity_ok': 'mean'
})

print("Constraint satisfaction evolution:")
print(by_iter)
```

### Custom Visualization

```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("all_results.csv")

# Create custom plot: λ₁ vs total decay width
fig, ax = plt.subplots(figsize=(10, 6))

viable = df[df['perturbativity_ok'] == 1]
invalid = df[df['perturbativity_ok'] == 0]

ax.scatter(invalid['lam1'], invalid['total_width'], alpha=0.3, s=10, c='red', label='Invalid')
ax.scatter(viable['lam1'], viable['total_width'], alpha=0.5, s=15, c='green', label='Viable')

ax.set_xlabel('λ₁')
ax.set_ylabel('Total Width (GeV)')
ax.set_title('λ₁ vs Total Higgs Width')
ax.legend()
plt.show()
```

## Troubleshooting

### No results.csv files found

**Problem**: `ERROR: No results.csv files found in campaign checkpoints`

**Solution**:
- Verify campaign actually completed: check `campaign_state.json` for `"campaign_state": "COMPLETED"`
- Check checkpoint directories exist: `ls runs/my-campaign/checkpoints/`
- Verify run directories were created: `find runs/my-campaign -name "results.csv" | wc -l`

### matplotlib/seaborn not installed

**Problem**: `WARNING: matplotlib/seaborn not installed. Skipping plots.`

**Solution**:
```bash
pip install matplotlib seaborn
```

Or run without plots:
```bash
python autoresearch/harness/aggregate_results.py runs/my-campaign --output results.csv
```

### Memory issues with large campaigns

**Problem**: Script crashes on campaigns with 10,000+ evaluations

**Solution**: Process in batches or use pandas chunking:
```python
# Process results in chunks
for chunk in pd.read_csv("all_results.csv", chunksize=1000):
    # Analyze each chunk
    print(chunk['lam1'].mean())
```

## Integration with Existing Tools

### Use with live_monitor.py

```bash
# Monitor campaign while it runs
python autoresearch/harness/live_monitor.py --watch runs/my-campaign

# When complete, analyze results
python autoresearch/harness/aggregate_results.py runs/my-campaign --plots
```

### Use with interactive_config.py

```bash
# Create new campaign config
python autoresearch/harness/interactive_config.py --wizard

# Run campaign
python autoresearch/run_supervisor.py configs/my-new-campaign.json

# Analyze results
python autoresearch/harness/aggregate_results.py runs/my-new-campaign --plots
```

## Next Steps

After analyzing your results:

1. **Identify promising regions**: Look for parameter combinations with high constraint satisfaction and interesting physics (e.g., high BR(h→γγ))

2. **Refine exploration**: Create a follow-up campaign focused on promising regions:
   ```bash
   python autoresearch/harness/interactive_config.py --template lambda1
   # Adjust λ₁ range to focus on interesting region
   ```

3. **Export for publication**: Use the aggregated CSV and plots for papers, presentations, or further ML analysis

4. **Continuous exploration**: Use time-budget mode for extended campaigns:
   ```json
   {
     "supervisor": {
       "max_rounds": 0,
       "max_duration_hours": 48.0,
       "enable_convergence": true
     }
   }
   ```
