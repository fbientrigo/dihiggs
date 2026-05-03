# Quick run: supported lam1 autoresearch setup

This setup is prepared for the **currently supported** lam1 backend slice:

- `tan_beta = 1000`
- `lambda_1 in [0, 12]`
- `m_phi = 130` fixed
- `m_A = m_H+ = 300` fixed by current backend
- `sin_ba = 0.999` fixed by current backend
- `lambda_6 = 0.1` fixed by current backend
- `lambda_7 = 0.0` fixed by current backend

## One-time environment

From repo root:

```bash
cd /home/fabi/dihiggs
export PYTHONPATH="$PWD:$PYTHONPATH"
export HB_DATASET="/home/fabi/dihiggs/datasets/HBDataset"
export HS_DATASET="/home/fabi/dihiggs/datasets/HSDataset"
```

## Optional preflight check

```bash
python autoresearch/run_supervisor.py autoresearch/configs/dihiggs_explorers_lam1_tb1000_supported.json --preflight-only
```

## Run the campaign

```bash
python autoresearch/run_supervisor.py autoresearch/configs/dihiggs_explorers_lam1_tb1000_supported.json --status
```

## Main output directory

```text
/home/fabi/dihiggs/runs/lam1-tb1000-mphi130-supported
```

## Files to inspect while it runs

- `campaign_state.json`
- `campaign_status.json`
- `campaign_status.html`
- `events.jsonl`
- `alerts.jsonl`
- `telemetry.db`
- `run_dirs/`

## Signs it is working

- `round_index` increases in `campaign_status.json`
- `events.jsonl` keeps growing
- `run_dirs/` gets new subdirectories
- `campaign_status.html` is updated
- `events_emitted` is nonzero in at least some rounds

## Important note

This prepared setup does **not** yet match your full requested physics slice. It uses the current lam1 backend exactly as wired in this repo today, without changing physics logic.
