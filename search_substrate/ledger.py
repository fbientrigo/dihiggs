from __future__ import annotations

import json
import fcntl
import os
import time
from pathlib import Path
from typing import Any, Iterator


class Ledger:
    """Append-only event ledger. Existing events are never rewritten."""

    def __init__(self, run_dir: Path):
        self.run_dir = run_dir
        self.path = run_dir / "ledger.jsonl"
        self.run_dir.mkdir(parents=True, exist_ok=True)

    def append(self, event: dict[str, Any]) -> dict[str, Any]:
        record = dict(event)
        record.setdefault("timestamp_unix", time.time())
        record.setdefault("event_id", f"event_{time.time_ns()}")
        with self.path.open("a", encoding="utf-8") as handle:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                handle.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
                handle.flush()
                os.fsync(handle.fileno())
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        return record

    def events(self) -> Iterator[dict[str, Any]]:
        if not self.path.exists():
            return
        with self.path.open(encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    yield json.loads(line)

    def evaluations(self) -> dict[str, dict[str, Any]]:
        latest: dict[str, dict[str, Any]] = {}
        for event in self.events():
            candidate_id = event.get("candidate_id")
            if candidate_id and event.get("event") == "EVALUATION_TERMINATED":
                latest[candidate_id] = event
        return latest

    def summary(self) -> dict[str, Any]:
        counts: dict[str, int] = {}
        for event in self.events():
            key = str(event.get("lifecycle", event.get("event", "UNKNOWN")))
            counts[key] = counts.get(key, 0) + 1
        return {"ledger": str(self.path), "event_counts": counts, "evaluations": len(self.evaluations())}
