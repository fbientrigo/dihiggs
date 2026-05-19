from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace('+00:00', 'Z')


DEFAULT_HEARTBEAT = {
    'run_id': None,
    'status': 'unknown',
    'pid': None,
    'updated_at_utc': None,
    'started_at_utc': None,
    'current_block': None,
    'total_blocks': None,
    'requested_attempted_points': None,
    'planned_points': 0,
    'attempted_points': 0,
    'raw_csv_rows': 0,
    'accepted_csv_rows': 0,
    'accepted_physics_points': 0,
    'learning_points': 0,
    'accepted_attempted_ratio': None,
    'learning_accepted_ratio': None,
    'best_ctau_m': None,
    'best_hyperparams': {},
    'last_csv_mtime': None,
    'last_event': None,
    'hard_stop_reason': None,
    'contract_satisfied': None,
    'evidence_paths': {},
}


def build_heartbeat(overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    hb = dict(DEFAULT_HEARTBEAT)
    if overrides:
        hb.update(overrides)
    hb['updated_at_utc'] = now_iso()
    return hb


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    os.replace(tmp, path)


def write_heartbeat(path: Path, current: dict[str, Any], updates: dict[str, Any] | None = None) -> dict[str, Any]:
    hb = dict(DEFAULT_HEARTBEAT)
    hb.update(current or {})
    if updates:
        hb.update(updates)
    hb['updated_at_utc'] = now_iso()
    atomic_write_json(path, hb)
    return hb


def read_heartbeat(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding='utf-8'))
    except Exception:
        return None
