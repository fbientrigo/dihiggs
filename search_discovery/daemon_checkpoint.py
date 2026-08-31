"""Durable daemon-level loop checkpoint.

Distinct from `cli.cmd_checkpoint`, which is a read-only stats report derived
by replaying the ledger and never writes a file. This module is the daemon's
own resumable position: cycle count, budget consumed, policy allocation
state, RNG state. Candidate-level dedup/resume still comes entirely from
replaying the ledger (see `daemon._replay_seen_ids`) -- this file only tracks
where the outer loop was, so there is exactly one source of truth for "was
this candidate evaluated" (the ledger) and exactly one for "where is the loop"
(this file).

Same atomic write discipline as the rest of the substrate: write a temp file,
`os.replace` into place.
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "dihiggs.llp_discovery.daemon_checkpoint.v1"


def path_for(run_dir: Path) -> Path:
    return Path(run_dir) / "daemon_checkpoint.json"


def load(run_dir: Path) -> dict[str, Any] | None:
    p = path_for(run_dir)
    if not p.exists():
        return None
    return json.loads(p.read_text(encoding="utf-8"))


def fresh(config_digest: str) -> dict[str, Any]:
    now = time.time()
    return {
        "schema_version": SCHEMA_VERSION,
        "cycle": 0,
        "total_evaluations": 0,
        "total_proposed": 0,
        "total_duplicates_dropped": 0,
        "config_digest": config_digest,
        "policy_state": {},
        "consecutive_no_improvement": 0,
        "consecutive_evaluator_failures": 0,
        "families_validated": 0,
        "started_utc": now,
        "last_checkpoint_utc": now,
        "stopped": False,
        "stop_reason": None,
    }


def save(run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    state = dict(state)
    state["last_checkpoint_utc"] = time.time()
    p = path_for(run_dir)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(state, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(tmp, p)
    return state
