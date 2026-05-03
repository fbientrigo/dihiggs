# Autoresearch Campaign Monitoring & Control Guide

Complete reference for checking campaign status, verifying progress, stopping/restarting jobs, and troubleshooting.

---

## 🆕 Live Monitor (Recommended)

**New automated tool for real-time monitoring with stale detection.**

```bash
# Single status check
python -m autoresearch.harness.live_monitor /home/fabi/dihiggs/runs/<campaign-name>

# Watch mode (auto-refresh every 10s)
python -m autoresearch.harness.live_monitor /home/fabi/dihiggs/runs/<campaign-name> --watch

# Compact mode (single line, great for tmux status bar)
python -m autoresearch.harness.live_monitor /home/fabi/dihiggs/runs/<campaign-name> --watch --compact

# Custom refresh interval
python -m autoresearch.harness.live_monitor /home/fabi/dihiggs/runs/<campaign-name> --watch --interval 5
```

**What it shows:**
- ✅ Campaign state (RUNNING/COMPLETED/CONVERGED)
- 🚨 Automatic stale detection
- Process status (PID, CPU, memory)
- Progress (events, rounds, run directories)
- Last activity timestamp
- Current metrics (yield, coverage, diversity)
- Convergence state and reason

**Stale detection logic:**
- State says RUNNING but no supervisor process → **STALE**
- No file activity in >5 minutes while process running → **STALE**
- Terminal states (COMPLETED/CONVERGED) → Not stale

## Quick Reference Commands

```bash
# Check if supervisor is running
ps aux | grep run_supervisor | grep -v grep

# Monitor real-time logs
tail -f /home/fabi/dihiggs/runs/<campaign-name>/supervisor.log

# Count completed evaluations
wc -l /home/fabi/dihiggs/runs/<campaign-name>/events.jsonl

# Count run directories created
ls -1d /home/fabi/dihiggs/runs/<campaign-name>/checkpoints/*/iter_*/run_* 2>/dev/null | wc -l

# Check last activity timestamp
ls -lht /home/fabi/dihiggs/runs/<campaign-name>/checkpoints/*/iter_*/run_* | head -5

# Stop running campaign
pkill -f "run_supervisor.*<campaign-name>"
```

---

## 1. Is My Campaign Actually Running?

### Check Process Status

```bash
# List all supervisor processes
ps aux | grep run_supervisor | grep -v grep

# Example output if running:
# fabi  12345  ... python -m autoresearch.harness.run_supervisor ...

# No output = no campaign running
```

**What to look for:**
- Process exists with `run_supervisor` in command line
- Check the config path in the command to identify which campaign
- Note the PID (process ID) for stopping later

### Check File Activity

```bash
# See most recently modified run directories
ls -lht /home/fabi/dihiggs/runs/<campaign-name>/checkpoints/*/iter_*/run_* | head -10

# Example output:
# drwxr-xr-x ... Apr  6 23:53 run_0042/
# drwxr-xr-x ... Apr  6 23:52 run_0041/
```

**Interpretation:**
- **Timestamps within last few minutes** → campaign is active
- **Timestamps hours/days old** → campaign is stale/stopped
- **No new files appearing** → campaign may be stuck or completed

---

## 2. Verify Progress (Not Stale)

### Count Total Evaluations

```bash
# Count ATTEMPT_EVALUATED events
grep -c '"event_type": "ATTEMPT_EVALUATED"' /home/fabi/dihiggs/runs/<campaign-name>/events.jsonl

# Or just count all events
wc -l /home/fabi/dihiggs/runs/<campaign-name>/events.jsonl
```

**Expected behavior:**
- Count should **increase** over time as campaign runs
- If count is static for 10+ minutes → campaign likely stopped/converged

### Count Run Directories

```bash
# Count all run_* directories
find /home/fabi/dihiggs/runs/<campaign-name>/checkpoints -type d -name "run_*" | wc -l

# Or per iteration
ls -1d /home/fabi/dihiggs/runs/<campaign-name>/checkpoints/*/iter_* | while read iter; do
  echo "$iter: $(ls -1d $iter/run_* 2>/dev/null | wc -l) runs"
done
```

**What to expect:**
- If `--n-proposals=100` and `max_rounds=20` → expect up to 2000 run directories
- If number stopped growing → check campaign state (see below)

### Check Campaign State

```bash
# Authoritative terminal state
jq -r '.campaign_state' /home/fabi/dihiggs/runs/<campaign-name>/campaign_state.json

# Possible values:
# - "RUNNING" → still active
# - "CONVERGED" → stopped due to convergence criteria
# - "COMPLETED" → stopped due to max_rounds limit
# - "FAILED" → encountered error
```

**Important distinction:**
- `campaign_state.json` → **authoritative** (trust this)
- `campaign_status.json` → last snapshot (may say RUNNING even after terminal)

### Check Current Round

```bash
# See which round is active/completed
jq -r '.round_number' /home/fabi/dihiggs/runs/<campaign-name>/campaign_state.json

# Compare to max_rounds in config
jq -r '.campaign.max_rounds' /home/fabi/dihiggs/autoresearch/configs/<your-config>.json
```

---

## 3. Real-Time Monitoring

### Tail Supervisor Logs

```bash
# Follow supervisor activity in real-time
tail -f /home/fabi/dihiggs/runs/<campaign-name>/supervisor.log

# What to look for:
# - "Starting round N" → new round began
# - "Executing arm: adaptive-smoke" → arm is running
# - "Round N complete" → round finished
# - "Campaign converged" or "Max rounds reached" → campaign stopping
```

### Watch Event Stream

```bash
# See new events as they're emitted
tail -f /home/fabi/dihiggs/runs/<campaign-name>/events.jsonl | jq -r '.event_type'

# Or with details
tail -f /home/fabi/dihiggs/runs/<campaign-name>/events.jsonl | jq '{event: .event_type, arm: .arm_id, round: .round_number}'
```

### Watch Run Directory Creation

```bash
# Count run directories every 10 seconds
watch -n 10 "find /home/fabi/dihiggs/runs/<campaign-name>/checkpoints -type d -name 'run_*' | wc -l"
```

### Monitor PhysLam1Scan Processes

```bash
# See active physics simulation processes
ps aux | grep PhysLam1Scan | grep -v grep

# Count how many running in parallel
ps aux | grep PhysLam1Scan | grep -v grep | wc -l
```

**Note:** If you see multiple `PhysLam1Scan` processes, that's normal - the supervisor runs proposals in parallel.

---

## 4. Stopping Campaigns

### Safe Stop (Graceful)

```bash
# Find the supervisor process
ps aux | grep "run_supervisor.*<campaign-name>" | grep -v grep

# Kill by pattern (safer if only one campaign running)
pkill -f "run_supervisor.*<campaign-name>"

# Or kill by PID
kill <PID>
```

**What happens:**
- Python receives SIGTERM signal
- May allow current operations to complete (depends on signal handling)
- State files (`campaign_state.json`, `events.jsonl`) remain intact
- Can resume later if desired

### Force Stop (Immediate)

```bash
# Force kill if graceful stop doesn't work
pkill -9 -f "run_supervisor.*<campaign-name>"

# Or by PID
kill -9 <PID>
```

**Use when:**
- Graceful stop didn't work after 30 seconds
- Campaign is hung/frozen
- You need immediate termination

**Risk:** May leave partial writes, but state files are usually consistent.

### Kill All Physics Processes (Emergency)

```bash
# Stop all PhysLam1Scan processes
pkill -f PhysLam1Scan

# Force kill if needed
pkill -9 -f PhysLam1Scan
```

**When to use:**
- Supervisor is killed but physics processes still running
- Runaway processes consuming resources
- Need clean slate before restart

### Verify Stop

```bash
# Confirm no supervisor running
ps aux | grep run_supervisor | grep -v grep

# Confirm no physics processes
ps aux | grep PhysLam1Scan | grep -v grep

# Both should return no results
```

---

## 5. Restarting & Resuming

### Resume Existing Campaign

**If you want to continue from where it left off:**

```bash
# Check current state
jq -r '.campaign_state' /home/fabi/dihiggs/runs/<campaign-name>/campaign_state.json

# If RUNNING or want to continue from terminal state:
python -m autoresearch.harness.run_supervisor \
  /home/fabi/dihiggs/autoresearch/configs/<your-config>.json
```

**Important:**
- If `campaign_state` is already **CONVERGED** or **COMPLETED**, supervisor will **exit immediately**
- To override terminal state, you need to manually edit `campaign_state.json` (not recommended)
- Better approach: start fresh campaign (see below)

### Start Fresh Campaign

**If you want to start over:**

```bash
# Option 1: Use new output directory
# Edit config: change "outdir" to new path
vim /home/fabi/dihiggs/autoresearch/configs/<your-config>.json
# Change: "outdir": "/home/fabi/dihiggs/runs/<NEW-campaign-name>"

# Option 2: Delete old campaign data
rm -rf /home/fabi/dihiggs/runs/<campaign-name>/*

# Then start supervisor
python -m autoresearch.harness.run_supervisor \
  /home/fabi/dihiggs/autoresearch/configs/<your-config>.json
```

**Recommendation:** Always use new `outdir` for fresh runs - safer and preserves old results for comparison.

### Restart After Config Change

**If you modified config and want to restart:**

```bash
# 1. Stop existing campaign
pkill -f "run_supervisor.*<campaign-name>"

# 2. Verify stop
ps aux | grep run_supervisor | grep -v grep

# 3. Start with modified config
python -m autoresearch.harness.run_supervisor \
  /home/fabi/dihiggs/autoresearch/configs/<your-config>.json
```

**Config changes that require restart:**
- `max_rounds`, `enable_convergence`, `--n-proposals`
- Any `arms[]` or `metrics` settings
- Convergence thresholds

**Config changes that DON'T require restart:**
- Comments, formatting
- Fields not read by supervisor

---

## 6. Troubleshooting Common Issues

### "Campaign converged too quickly"

**Symptoms:**
- Only ran a few rounds (5-10) before stopping
- Far fewer evaluations than expected
- Metrics look flat/unchanging

**Diagnosis:**
```bash
# Check convergence settings
jq '.campaign.convergence' /home/fabi/dihiggs/autoresearch/configs/<your-config>.json

# Check actual rounds completed
jq -r '.round_number' /home/fabi/dihiggs/runs/<campaign-name>/campaign_state.json
```

**Fix:**
- Set `"enable_convergence": false` to disable convergence detection
- Or increase `plateau_rounds_threshold` (e.g., 50)
- Or tighten `metric_delta_threshold` (e.g., 0.001)
- Use fresh `outdir` and restart

### "Only 1 point evaluated per round"

**Symptoms:**
- `events.jsonl` has very few lines
- All events show same lam1 value
- Run directories count << `--n-proposals` × rounds

**Diagnosis:**
```bash
# Check how many proposals per round
jq -r '.arms[].cmd' /home/fabi/dihiggs/autoresearch/configs/<your-config>.json | grep -o 'n-proposals [0-9]*'

# Check run directory count
find /home/fabi/dihiggs/runs/<campaign-name>/checkpoints -type d -name "run_*" | wc -l
```

**Fix:**
- This was the bug we fixed in `adaptive_explorer_lam1.py`
- Ensure you're using the updated version
- Set `--n-proposals 100` (or desired count) in config
- Set `limits.max_new_run_dirs_per_round` and `max_new_run_dirs_per_arm_call` high enough (e.g., 200)

### "Campaign says RUNNING but no progress"

**Symptoms:**
- `campaign_status.json` shows "RUNNING"
- But no new files, no process, timestamps old

**Diagnosis:**
```bash
# Check authoritative state
jq -r '.campaign_state' /home/fabi/dihiggs/runs/<campaign-name>/campaign_state.json

# Check process
ps aux | grep run_supervisor | grep -v grep

# Check last file activity
ls -lht /home/fabi/dihiggs/runs/<campaign-name>/checkpoints/*/iter_*/run_* | head -5
```

**Explanation:**
- `campaign_status.json` is last **snapshot**, not live state
- `campaign_state.json` is **authoritative**
- If `campaign_state` is terminal and no process → campaign finished

**Fix:**
- Trust `campaign_state.json`, ignore `campaign_status.json`
- Start new campaign with new `outdir`

### "Supervisor starts but exits immediately"

**Symptoms:**
- Run supervisor command
- Process exits within seconds
- No error message

**Diagnosis:**
```bash
# Check campaign state
jq -r '.campaign_state' /home/fabi/dihiggs/runs/<campaign-name>/campaign_state.json

# Check supervisor log
tail -50 /home/fabi/dihiggs/runs/<campaign-name>/supervisor.log
```

**Likely cause:**
- Campaign already in terminal state (CONVERGED/COMPLETED)
- Supervisor sees terminal state, exits early (lines 155-156 in `campaign_supervisor.py`)

**Fix:**
- Use new `outdir` for fresh run
- Or delete old state: `rm /home/fabi/dihiggs/runs/<campaign-name>/campaign_state.json`

### "PhysLam1Scan processes still running after stop"

**Symptoms:**
- Killed supervisor
- But `ps aux | grep PhysLam1Scan` shows processes

**Fix:**
```bash
# Kill all physics processes
pkill -f PhysLam1Scan

# Force if needed
pkill -9 -f PhysLam1Scan
```

**Why this happens:**
- Supervisor spawns subprocesses
- Killing supervisor doesn't auto-kill children
- Need explicit cleanup

---

## 7. Best Practices

### Before Starting Campaign
- [ ] Choose unique `outdir` name (include date/parameters)
- [ ] Set `max_rounds` limit (avoid infinite runs)
- [ ] Review `--n-proposals` and `max_new_run_dirs_per_*` limits
- [ ] Decide on convergence: enable or disable
- [ ] Run preflight check: `python -m autoresearch.harness.run_supervisor <config> --preflight-only`

### During Campaign
- [ ] Monitor logs: `tail -f supervisor.log`
- [ ] Check progress every ~30 min: count run directories
- [ ] Verify metrics make sense: `jq '.metrics' campaign_status.json`
- [ ] Watch for premature convergence

### After Campaign
- [ ] Check terminal state: `campaign_state.json`
- [ ] Verify expected evaluation count
- [ ] Archive results if successful
- [ ] Clean up failed/test runs: `rm -rf runs/<test-campaign>`

### When Stopping Campaign
- [ ] Use graceful stop first: `pkill -f run_supervisor`
- [ ] Wait 30 seconds before force kill
- [ ] Check for orphaned PhysLam1Scan processes
- [ ] Verify state files are intact

### When Restarting
- [ ] Always use new `outdir` for fresh campaigns
- [ ] Don't reuse terminal state directories
- [ ] Review config changes before restart
- [ ] Confirm old processes are dead

---

## 8. File Reference

### Campaign Output Structure

```
/home/fabi/dihiggs/runs/<campaign-name>/
├── campaign_state.json          # Authoritative terminal state (TRUST THIS)
├── campaign_status.json         # Last snapshot (may be stale)
├── events.jsonl                 # All ATTEMPT_EVALUATED events
├── supervisor_events.jsonl      # Supervisor round events
├── supervisor.log               # Human-readable log
├── snapshot_history/            # Historical snapshots
│   ├── snapshot_r00.json
│   ├── snapshot_r01.json
│   └── ...
└── checkpoints/
    └── adaptive-smoke/          # Arm name
        ├── iter_00/             # Round 0
        │   ├── adaptive_state.json
        │   ├── run_00000/
        │   │   ├── results.csv  # PhysLam1Scan output
        │   │   └── stdout.txt
        │   ├── run_00001/
        │   └── ...
        ├── iter_01/             # Round 1
        └── ...
```

### Key Files Explained

| File | Purpose | When Updated | Trust Level |
|------|---------|--------------|-------------|
| `campaign_state.json` | Authoritative state | On terminal transition | **TRUST** |
| `campaign_status.json` | Dashboard snapshot | Before terminal transition | Stale after terminal |
| `events.jsonl` | Event log | On every ATTEMPT_EVALUATED | Accurate |
| `supervisor.log` | Human log | Real-time | Good for debugging |
| `adaptive_state.json` | Bandit state per round | After each arm execution | Round-specific |
| `results.csv` | Physics output | After each PhysLam1Scan | Per-run data |

---

## 9. Example Workflows

### Workflow: Start and Monitor New Campaign

```bash
# 1. Prepare config (use unique outdir)
vim /home/fabi/dihiggs/autoresearch/configs/my_campaign.json
# Set: "outdir": "/home/fabi/dihiggs/runs/lam1-test-$(date +%Y%m%d-%H%M)"

# 2. Run preflight check
python -m autoresearch.harness.run_supervisor \
  /home/fabi/dihiggs/autoresearch/configs/my_campaign.json \
  --preflight-only

# 3. Start campaign
python -m autoresearch.harness.run_supervisor \
  /home/fabi/dihiggs/autoresearch/configs/my_campaign.json \
  > /tmp/supervisor_stdout.log 2>&1 &

# 4. Note the PID
echo $! > /tmp/supervisor.pid

# 5. Monitor in real-time
tail -f /home/fabi/dihiggs/runs/lam1-test-*/supervisor.log

# 6. Check progress periodically
watch -n 30 "find /home/fabi/dihiggs/runs/lam1-test-*/checkpoints -type d -name 'run_*' | wc -l"
```

### Workflow: Stop and Restart with New Config

```bash
# 1. Stop existing campaign
pkill -f "run_supervisor.*my_campaign"

# 2. Verify stop
ps aux | grep run_supervisor | grep -v grep

# 3. Edit config
vim /home/fabi/dihiggs/autoresearch/configs/my_campaign.json
# Change settings (e.g., increase --n-proposals)
# Change "outdir" to new path

# 4. Restart
python -m autoresearch.harness.run_supervisor \
  /home/fabi/dihiggs/autoresearch/configs/my_campaign.json \
  > /tmp/supervisor_stdout.log 2>&1 &
```

### Workflow: Investigate Stalled Campaign

```bash
# 1. Check process
ps aux | grep run_supervisor | grep -v grep
# Result: no process found

# 2. Check authoritative state
jq -r '.campaign_state' /home/fabi/dihiggs/runs/<campaign-name>/campaign_state.json
# Result: "COMPLETED"

# 3. Check round number
jq -r '.round_number' /home/fabi/dihiggs/runs/<campaign-name>/campaign_state.json
# Result: 19 (reached max_rounds=20)

# 4. Check evaluation count
wc -l /home/fabi/dihiggs/runs/<campaign-name>/events.jsonl
# Result: 2000 events

# Conclusion: Campaign completed successfully, not stalled
```

---

## 10. Advanced: Understanding State Transitions

### Normal Campaign Lifecycle

```
INITIALIZE → campaign_state="RUNNING"
    ↓
ROUND 0 → execute arms → parse checkpoints → emit events
    ↓
ROUND 1 → ...
    ↓
ROUND N → check convergence OR max_rounds
    ↓
TERMINAL STATE → campaign_state="CONVERGED" or "COMPLETED"
    ↓
PERSIST STATE → write campaign_state.json
    ↓
EXIT
```

### Resume Behavior

```
START → load campaign_state.json
    ↓
IF campaign_state IN [CONVERGED, COMPLETED, FAILED]
    ↓
    EXIT IMMEDIATELY (lines 155-156)
    ↓
ELSE
    ↓
    CONTINUE FROM round_number
```

**Key insight:** Reusing same `outdir` means reusing state. If terminal, supervisor exits.

---

## Questions?

If you encounter issues not covered here:

1. **Check logs first**: `tail -100 supervisor.log`
2. **Verify state files**: `jq . campaign_state.json`
3. **Count artifacts**: run directories, events, processes
4. **Check timestamps**: `ls -lht checkpoints/*/iter_*/run_*`

Most issues fall into:
- Campaign already terminal (use new outdir)
- Convergence too aggressive (disable or tune)
- Proposals not respected (check `--n-proposals` and limits)
- Orphaned processes (pkill cleanup)

---

**Last Updated:** Apr 7, 2026  
**Autoresearch Version:** lam1 backend with fixed n-proposals loop
