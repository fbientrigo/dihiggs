"""Append-only SQLite journal for the daemon: decisions, write-ahead intents and integer budget charges.

Smallest version of research-policy-engine ADR-019/022 that autoresearch needs. Rows are never updated or deleted
(triggers abort it); event ids are deterministic (`intent:a0003`), so re-appending the same id is a no-op. An attempt
reserves its upper-bound cost in `intent:` before anything launches and is charged its actual cost in `closed:`; status
and ledger are derived from events, never stored, so a kill -9 can never leave them half-written.
Costs are integers: cpu_core_s (threads x wall of the boundary process), wall_s, disk_mb.
"""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path

RESOURCES = ("cpu_core_s", "wall_s", "disk_mb")


class BudgetExhausted(Exception):
    """The next reservation does not fit: a clean stop, not a failure."""


class Journal:
    def __init__(self, path: str | Path):
        self.db = sqlite3.connect(str(path), isolation_level=None)
        self.db.execute("PRAGMA journal_mode=WAL")
        self.db.execute("PRAGMA synchronous=FULL")
        self.db.execute("CREATE TABLE IF NOT EXISTS events(seq INTEGER PRIMARY KEY AUTOINCREMENT, event_id TEXT NOT NULL UNIQUE,"
                        " kind TEXT NOT NULL, utc TEXT NOT NULL DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ','now')), payload TEXT NOT NULL)")
        for op in ("UPDATE", "DELETE"):
            self.db.execute(f"CREATE TRIGGER IF NOT EXISTS events_no_{op.lower()} BEFORE {op} ON events "
                            "BEGIN SELECT RAISE(ABORT, 'journal is append-only'); END")

    def append(self, event_id: str, kind: str, payload: dict) -> dict:
        """Idempotent: an existing event_id returns its stored payload (the first write wins)."""
        row = self.db.execute("SELECT payload FROM events WHERE event_id=?", (event_id,)).fetchone()
        if row:
            return json.loads(row[0])
        self.db.execute("INSERT INTO events(event_id, kind, payload) VALUES(?,?,?)", (event_id, kind, json.dumps(payload, sort_keys=True)))
        return payload

    def get(self, event_id: str) -> dict | None:
        row = self.db.execute("SELECT payload FROM events WHERE event_id=?", (event_id,)).fetchone()
        return json.loads(row[0]) if row else None

    def events(self, kind: str | None = None) -> list[tuple[str, str, dict]]:
        sql, args = "SELECT event_id, kind, payload FROM events", ()
        if kind:
            sql, args = sql + " WHERE kind=?", (kind,)
        return [(e, k, json.loads(p)) for e, k, p in self.db.execute(sql + " ORDER BY seq", args)]

    def open_attempts(self) -> list[str]:
        """Attempts with an intent and no close: what a restart must reconcile, never relaunch."""
        return [eid.split(":", 1)[1] for eid, _, _ in self.events("intent") if not self.get("closed:" + eid.split(":", 1)[1])]

    def ledger(self) -> dict[str, int]:
        """Charged + still-reserved per resource: closed attempts count their actual, open ones their reservation."""
        out = dict.fromkeys(RESOURCES, 0)
        for eid, _, intent in self.events("intent"):
            closed = self.get("closed:" + eid.split(":", 1)[1])
            for r, amount in (closed["actual"] if closed else intent["reserved"]).items():
                out[r] += amount
        return out

    def reserve(self, attempt: str, reserved: dict[str, int], limits: dict[str, int], payload: dict) -> dict:
        """Write-ahead intent with its reservation, refused (BudgetExhausted) if it would cross any limit."""
        if self.get("intent:" + attempt):
            return self.get("intent:" + attempt)
        used = self.ledger()
        over = [r for r in RESOURCES if r in limits and used[r] + reserved[r] > limits[r]]
        if over:
            raise BudgetExhausted(", ".join(f"{r}: used {used[r]} + reserve {reserved[r]} > limit {limits[r]}" for r in over))
        return self.append("intent:" + attempt, "intent", {**payload, "reserved": reserved})

    def close(self, attempt: str, actual: dict[str, int], payload: dict) -> dict:
        intent = self.get("intent:" + attempt)
        if intent is None:
            raise ValueError(f"close without intent: {attempt}")
        return self.append("closed:" + attempt, "closed", {**payload, "actual": actual, "reserved": intent["reserved"]})
