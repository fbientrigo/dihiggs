## Lam1 config: what actually changes execution

This backend does **not** interpret `search.lam1.n_bins` as "run 2000 lambda1 points".

### What the current lam1 backend really does

- `search.lam1.min` / `search.lam1.max`
  - define the lambda1 range passed to `adaptive_explorer_lam1.py`
- `search.lam1.n_bins`
  - is passed through to the wrapper CLI
  - is used for post-hoc `lam1_bin` bucketing in the harness
  - does **not** make the backend execute that many physical lambda1 values
- `arms[].cmd ... --n-proposals 1`
  - makes the current wrapper choose exactly one lambda1 value: `lam1_min`
  - with your config, that means every round proposes `lambda1 = 0.0`
- `metrics.multi_axis.coverage_axes[].domain_size`
  - controls the lam1 denominator used in reported coverage
  - if this says `24`, coverage is computed with 24 lam1 buckets even if `search.lam1.n_bins` is 2000

### What this means for your current config

With:

```json
"search": {
  "lam1": { "min": 0.0, "max": 12.0, "n_bins": 2000 }
},
"arms": [{
  "cmd": ["...", "--n-proposals", "1", "..."]
}]
```

the current lam1 backend does **not** scan 2000 lambda1 points.
It runs one single-point `PhysLam1Scan` call per round, and with `--n-proposals 1` that point is repeatedly `lambda1 = 0.0`.

### Safe config changes that do help

1. Use a **fresh** `paths.outdir` for every new experiment.
2. If you want status metrics to match the intended lam1 granularity, set:

```json
"metrics": {
  "multi_axis": {
    "coverage_axes": [
      { "name": "tb", "kind": "categorical", "domain": [1000], "weight": 0.5 },
      { "name": "lam1_bin", "kind": "discrete", "domain_size": 2000, "weight": 0.5 }
    ]
  }
}
```

3. While debugging this backend, disable convergence so the campaign does not stop early on flat metrics:

```json
"supervisor": {
  "enable_convergence": false,
  "max_rounds": 20
}
```

### Config changes that do **not** do what they seem to do

- Increasing `search.lam1.n_bins` does **not** increase the number of physical lambda1 evaluations.
- Changing `search.mphi` does **not** make the current lam1 wrapper scan that full mphi range.
- Increasing `--n-proposals` in the current wrapper does **not** create that many proposals; it only changes which single lambda1 value gets chosen.

### How to verify what actually ran

- `events.jsonl`: each `ATTEMPT_EVALUATED` line is one emitted attempt
- `checkpoints/.../adaptive_state.json`: shows the actual `lam1_min` / `lam1_max` used that round
- `campaign_state.json`: authoritative terminal state
- `campaign_status.json`: last dashboard snapshot, which can still say `RUNNING` after terminal convergence

### Bottom line

If you want **2000 actual lambda1 evaluations**, that requires a code change in `dihiggs/app/adaptive_explorer_lam1.py` (or a different backend), not just a JSON edit.
