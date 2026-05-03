# Handoff Runbook: Generate New Runs From Existing Datalake Proposals

This handoff explains how to operate the current MVP flow so new exploration decisions come from the existing datalake state.

## Goal

Use the datalake at:

`/mnt/c/Users/Asus/cern_db/dihiggs_lake`

to:

1. scan and normalize existing run artifacts,
2. reconcile coverage/discovery/graph state,
3. generate bounded next-run proposals,
4. validate behavior with repeatable tests.

---

## What Is Implemented (MVP)

- Datalake scanner + manifest: `autoresearch/harness/datalake_manifest.py`
- Registry/discovery/coverage upsert pipeline: `autoresearch/harness/mvp_upsert_pipeline.py`
- Graph delta updates: `autoresearch/harness/mvp_graph.py`
- Proposal engine: `autoresearch/harness/mvp_goal_proposer.py`
- Fixed-interval autonomy scheduler: `autoresearch/harness/autonomy_scheduler.py`
- Dry-run dispatch adapter: `autoresearch/harness/orchestrator_adapter.py`
- Supervisor runner entrypoint: `autoresearch/run_supervisor.py`

> Current dispatch behavior is dry-run/enqueue-style by default (safe mode). You will get structured proposals + audit trail without force-launching external orchestration.

---

## 1) One-Time Environment Setup

From repo root:

```bash
cd /home/fabi/dihiggs
export PYTHONPATH="$PWD:$PYTHONPATH"
export HB_DATASET="/home/fabi/dihiggs/datasets/HBDataset"
export HS_DATASET="/home/fabi/dihiggs/datasets/HSDataset"
```

---

## 2) Verify Config + Preflight

Use your campaign config (example shown):

```bash
python autoresearch/run_supervisor.py autoresearch/configs/dihiggs_explorers_lam1.json --preflight-only
```

Expected: JSON output with preflight checks and no blocking failure.

---

## 3) Scan the Existing Datalake

Run the scanner directly:

```bash
python -m autoresearch.harness.datalake_manifest /mnt/c/Users/Asus/cern_db/dihiggs_lake --pretty
```

Expected:

- `status` is `success` or `success_with_quarantine`
- non-zero `campaign_count` and `manifest_row_count`

If status is `root_not_found`, fix the mount/path first.

---

## 4) Run One Scheduler Tick From Datalake State (Proposal Generation)

This command performs one full reconcile → proposal → dispatch(dry-run) tick and writes a machine-readable scheduler snapshot.

```bash
python - <<'PY'
from pathlib import Path
import json

from autoresearch.harness.datalake_manifest import scan_datalake
from autoresearch.harness.telemetry_store import init_db
from autoresearch.harness.autonomy_scheduler import FixedIntervalAutonomyScheduler

repo = Path('/home/fabi/dihiggs')
cfg_path = repo / 'autoresearch/configs/dihiggs_explorers_lam1.json'
cfg = json.loads(cfg_path.read_text(encoding='utf-8'))

lake_root = Path('/mnt/c/Users/Asus/cern_db/dihiggs_lake')
report = scan_datalake(lake_root)

if report['status'] == 'root_not_found':
    raise SystemExit('Datalake root not found: check mount/path')

outdir = Path(str(cfg['paths']['outdir']))
outdir.mkdir(parents=True, exist_ok=True)

manifest_snapshot = outdir / 'datalake_manifest_snapshot.json'
manifest_snapshot.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n', encoding='utf-8')

conn = init_db(outdir / 'telemetry.db')
scheduler = FixedIntervalAutonomyScheduler(
    conn,
    config=cfg,
    source_file=manifest_snapshot,
    state_path=outdir / 'autonomy_scheduler_state.json',
)

result = scheduler.tick(scanner_rows=report['manifest'])
print(json.dumps({
    'status': result['status'],
    'tick_index': result['tick_index'],
    'proposal_count': result['proposal_count'],
    'dispatch_count': result['dispatch_count'],
    'quarantine': result['quarantine'],
    'state_path': result['state_path'],
}, indent=2, sort_keys=True))
PY
```

Inspect:

- `<outdir>/autonomy_scheduler_state.json`
  - `stage_telemetry`
  - `proposal_audit`
  - `quarantine`
  - `proposals`

These are the proposed new spaces produced from current datalake-derived state.

---

## 5) Execute Campaign Rounds (Operational Run)

Run supervisor loop with status print:

```bash
python autoresearch/run_supervisor.py autoresearch/configs/dihiggs_explorers_lam1.json --status
```

Main outputs (from `paths.outdir`):

- `campaign_state.json`
- `campaign_status.json`
- `events.jsonl`
- `supervisor_events.jsonl`
- `telemetry.db`

Use this when you want the campaign to actually keep operating round-by-round.

---

## 6) How to Test It Yourself (Recommended Order)

### A. Targeted autonomy/proposal tests

```bash
cd /home/fabi/dihiggs/autoresearch
pytest tests/test_autonomy_scheduler.py tests/test_mvp_goal_proposer.py tests/test_supervisor_telemetry.py -q
```

### B. End-to-end supervisor/autonomy dry-run checks

```bash
pytest tests/test_supervisor_e2e.py tests/test_task11_trigger.py -q
```

### C. Full regression suite

```bash
pytest tests -q
```

Expected current baseline:

- full suite passes (historically: `221 passed, 2 skipped`)

---

## 7) Success Criteria Checklist

- [ ] Datalake scanner returns `success` or `success_with_quarantine`
- [ ] Scheduler tick returns non-error status and writes `autonomy_scheduler_state.json`
- [ ] `proposal_audit` exists and contains per-goal decision outcomes
- [ ] `quarantine` block is explicit (`clear` or escalated status)
- [ ] Supervisor run updates campaign artifacts in outdir
- [ ] Tests pass in the order above

---

## 8) If Proposals Are Not Appearing

1. Confirm datalake path exists and is readable.
2. Check scanner report counts (`campaign_count`, `manifest_row_count`).
3. Check scheduler snapshot `quarantine` and `stage_telemetry` for blocked reasons.
4. Verify proposal budget knobs in config:
   - `autonomy_scheduler.max_proposals`
   - `autonomy_scheduler.max_dispatches`
5. Re-run targeted tests in section 6A.

---

## Notes for Handoff Consumer

- Use this runbook as the operational baseline for datalake-driven proposal generation.
- Keep dry-run dispatch enabled until you explicitly wire/validate a non-dry orchestrator adapter.
