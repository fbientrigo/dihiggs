# Time-Budget Continuous Exploration

## Overview

Autoresearch now supports **time-based stopping conditions** for continuous, autonomous exploration campaigns. Instead of specifying a fixed number of rounds, you can configure campaigns to run for a specific duration or until a deadline.

## Key Features

### 1. Duration-Based Stopping

Run campaigns for a fixed number of hours:

```json
{
  "supervisor": {
    "max_rounds": 0,
    "max_duration_hours": 8.0,
    "enable_convergence": true
  }
}
```

**Behavior**:
- Campaign runs for exactly 8 hours (or until convergence, whichever comes first)
- `max_rounds: 0` means no round limit
- Adaptive algorithm continues exploring until time budget is exhausted
- Each round's elapsed time is tracked in `campaign_state.json`

### 2. Deadline-Based Stopping

Run until a specific timestamp:

```json
{
  "supervisor": {
    "max_rounds": 0,
    "stop_at_timestamp": 1735689600.0,
    "enable_convergence": true
  }
}
```

**Behavior**:
- Campaign stops at the specified Unix timestamp
- Useful for scheduled runs (e.g., "stop at 5pm")
- Can be combined with `max_duration_hours` (whichever comes first wins)

### 3. Unlimited Exploration with Convergence

Run indefinitely until the algorithm converges:

```json
{
  "supervisor": {
    "max_rounds": 0,
    "enable_convergence": true
  }
}
```

**Behavior**:
- No time or round limits
- Campaign stops only when convergence detector signals completion
- Ideal for "explore until we've covered the space" scenarios

## Configuration Options

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `max_rounds` | int | 0 | Hard cap on adaptation rounds. Use `0` for no limit. |
| `max_duration_hours` | float | None | Maximum campaign duration in hours. |
| `stop_at_timestamp` | float | None | Unix timestamp to stop at. |
| `enable_convergence` | bool | true | Whether to stop on convergence. |

**Priority order** (first condition met triggers stop):
1. Exception/failure
2. Time budget exceeded (`max_duration_hours` or `stop_at_timestamp`)
3. Round limit reached (`max_rounds` > 0)
4. Convergence detected (`enable_convergence: true`)

## Interactive Configuration

Use the interactive wizard to set time budgets:

```bash
python autoresearch/harness/interactive_config.py --wizard
```

The wizard will prompt:

```
=== Campaign Settings ===

Max rounds (0 = run until convergence) [0]: 0
Enable convergence detection? (y/n) [y]: y
Use time budget (instead of/in addition to max_rounds)? (y/n) [n]: y
Maximum campaign duration (hours) [8.0]: 12.0
  ✓ Campaign will stop after 12.0 hours
```

## Monitoring Time-Budget Campaigns

### Check Elapsed Time

```bash
# View current state
python autoresearch/harness/live_monitor.py runs/my-campaign

# Output shows:
# Elapsed time: 3.2 hours (of 8.0 hour budget)
# Estimated rounds remaining: ~15 (based on avg round time)
```

### Real-Time Progress

```bash
# Watch mode with time tracking
python autoresearch/harness/live_monitor.py --watch runs/my-campaign

# Updates every 10 seconds:
# Round 14/∞ | Elapsed: 4.1h / 8.0h (51%) | ETA: 3.9h
```

### Campaign State

The `campaign_state.json` file tracks timing:

```json
{
  "campaign_state": "RUNNING",
  "campaign_start_time": 1735689600.0,
  "elapsed_hours": 4.125,
  "round_index": 14,
  "stop_reason": null
}
```

When time budget is exceeded:

```json
{
  "campaign_state": "COMPLETED",
  "elapsed_hours": 8.001,
  "stop_reason": "time_budget_exceeded (max_duration_hours=8.0, elapsed=8.00)"
}
```

## Use Cases

### 1. Overnight Exploration

```json
{
  "supervisor": {
    "max_rounds": 0,
    "max_duration_hours": 12.0,
    "enable_convergence": true
  }
}
```

**Scenario**: Start a campaign before leaving work, let it run overnight, and review results in the morning.

### 2. Continuous Research Assistant

```json
{
  "supervisor": {
    "max_rounds": 0,
    "enable_convergence": true,
    "checkpoint_interval": 1
  }
}
```

**Scenario**: Run indefinitely on a dedicated server, exploring parameter space autonomously until convergence.

### 3. Scheduled Campaigns

```python
import time
import json

# Stop at 5pm today
stop_time = time.mktime(time.strptime("2025-04-07 17:00:00", "%Y-%m-%d %H:%M:%S"))

config = {
    "supervisor": {
        "max_rounds": 0,
        "stop_at_timestamp": stop_time,
        "enable_convergence": false
    }
}

with open("configs/scheduled_campaign.json", "w") as f:
    json.dump(config, f, indent=2)
```

**Scenario**: Ensure campaign completes before a deadline (e.g., before daily backup, before presentation).

### 4. Budget-Constrained Exploration

```json
{
  "supervisor": {
    "max_rounds": 100,
    "max_duration_hours": 6.0,
    "enable_convergence": true
  }
}
```

**Scenario**: Explore as much as possible within both time and iteration constraints. Stops at whichever limit is reached first.

## Best Practices

### 1. Combine Time Budget with Convergence

Always enable convergence detection even with time budgets:

```json
{
  "supervisor": {
    "max_duration_hours": 8.0,
    "enable_convergence": true
  }
}
```

**Why**: Campaign may converge before time budget is exhausted, saving resources.

### 2. Use Checkpointing

Ensure frequent checkpoints for long campaigns:

```json
{
  "supervisor": {
    "checkpoint_interval": 1,
    "max_duration_hours": 24.0
  }
}
```

**Why**: If campaign is interrupted, state can be recovered from checkpoints.

### 3. Monitor Resource Usage

For multi-hour campaigns, monitor system resources:

```bash
# Check CPU/memory usage
htop

# Check disk usage (results accumulate)
du -sh runs/my-campaign/

# Monitor campaign status
watch -n 30 'python autoresearch/harness/live_monitor.py runs/my-campaign'
```

### 4. Set Realistic Time Budgets

Estimate round duration before setting time budgets:

```bash
# Run a short test campaign
python autoresearch/run_supervisor.py configs/test.json

# Check round timing in logs
grep "round_index" runs/test/campaign_state.json

# Calculate average round time
# If 10 rounds in 2 hours = 12 minutes/round
# For 100 rounds: set max_duration_hours = 20.0 (with buffer)
```

## Troubleshooting

### Campaign stops immediately

**Problem**: `stop_reason: "time_budget_exceeded (max_duration_hours=0.0, elapsed=0.00)"`

**Solution**: Check config - `max_duration_hours: 0.0` means zero budget. Use `null` or remove field for no limit:
```json
{
  "supervisor": {
    "max_duration_hours": null
  }
}
```

### Time budget not respected

**Problem**: Campaign runs beyond specified duration

**Check**:
1. Verify config was loaded: `cat runs/my-campaign/campaign_state.json | jq .elapsed_hours`
2. Check if exception occurred: `cat runs/my-campaign/campaign_state.json | jq .last_error`
3. Ensure campaign started with updated supervisor code

### Resume behavior with time budget

**Problem**: Resumed campaign ignores original time budget

**Solution**: Time budget is cumulative across resumes:
- `campaign_start_time` is preserved from original run
- `elapsed_hours` accumulates across all sessions
- To reset timer, delete `campaign_state.json` before resuming (⚠️ loses other state)

## Migration from Fixed Rounds

### Old Configuration

```json
{
  "supervisor": {
    "max_rounds": 50,
    "enable_convergence": false
  }
}
```

**Behavior**: Exactly 50 rounds, no early stopping.

### New Configuration (Time-Based)

```json
{
  "supervisor": {
    "max_rounds": 0,
    "max_duration_hours": 10.0,
    "enable_convergence": true
  }
}
```

**Behavior**: Run for 10 hours OR until convergence, whichever comes first. More adaptive, resource-efficient.

## Advanced: Dynamic Time Budgets

Calculate time budgets programmatically:

```python
import json
import time

def create_timed_config(base_config_path, duration_hours, output_path):
    """Create a config with a specific time budget."""
    with open(base_config_path) as f:
        config = json.load(f)
    
    config["supervisor"]["max_rounds"] = 0
    config["supervisor"]["max_duration_hours"] = duration_hours
    config["supervisor"]["enable_convergence"] = True
    
    with open(output_path, "w") as f:
        json.dump(config, f, indent=2)
    
    print(f"✓ Created config with {duration_hours}h budget: {output_path}")

# Example: Create configs for different time budgets
create_timed_config("configs/lambda1.json", 2.0, "configs/lambda1_2h.json")
create_timed_config("configs/lambda1.json", 8.0, "configs/lambda1_8h.json")
create_timed_config("configs/lambda1.json", 24.0, "configs/lambda1_24h.json")
```

## Integration with Existing Features

### Time Budget + Autoscaling

```json
{
  "supervisor": {
    "max_duration_hours": 12.0,
    "enable_convergence": true
  },
  "autoscaling": {
    "enabled": true,
    "scale_threads": true,
    "scale_timeout": true
  }
}
```

**Behavior**: Campaign runs for 12 hours, automatically adjusting thread count and timeouts based on convergence signals.

### Time Budget + Alerts

```json
{
  "supervisor": {
    "max_duration_hours": 8.0
  },
  "alerts": {
    "channels": ["stderr", "alerts.jsonl"],
    "severity_levels": ["warning", "error"]
  }
}
```

**Behavior**: Alerts are emitted throughout the campaign, and a final alert signals time budget exhaustion.

## Next Steps

1. **Try it out**: Create a time-budgeted config with the interactive wizard
2. **Monitor progress**: Use `live_monitor.py --watch` to track elapsed time
3. **Analyze results**: Use `aggregate_results.py` to visualize exploration efficiency
4. **Optimize**: Adjust time budgets based on observed round duration and convergence behavior

## Belief Updating (Already Implemented)

**Note**: The autoresearch system already implements belief updating via:
- **UCB1 bandit algorithm**: Balances exploration/exploitation across different explorer arms
- **Bayesian budget allocation**: Beta-binomial posterior for allocating points across bins
- **Composite reward**: Combines yield, coverage, and diversity for arm selection

Time-budget mode works seamlessly with these existing adaptive mechanisms. The campaign continues updating beliefs and refining exploration until the time limit is reached.
