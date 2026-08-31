"""Code-owned discovery archive: Top-5 diverse families per provisional family.

Separate from search_substrate.archive, which stays fail-closed and untouched.
Callers can never supply archive state, a family label, or a score; they submit
a completed cross-mass family record and the code decides.
"""
from __future__ import annotations

import fcntl
import json
import math
import os
from pathlib import Path
from typing import Any

from .bounds import Envelope
from .helpers import to_unit_vector

FAMILIES = ("mixed", "photonic")
TOP_K = 5
NEAR_DUPLICATE_RADIUS = 0.05


def family_distance(left: dict[str, Any], right: dict[str, Any], envelope: Envelope) -> float:
    a = to_unit_vector(left["coordinates"], envelope)
    b = to_unit_vector(right["coordinates"], envelope)
    return float(math.sqrt(sum((a - b) ** 2) / len(a)))


class DiscoveryArchive:
    def __init__(self, run_dir: Path):
        self.run_dir = Path(run_dir)
        self.path = self.run_dir / "discovery_archive.json"
        self.history_path = self.run_dir / "discovery_archive_history.jsonl"
        self.run_dir.mkdir(parents=True, exist_ok=True)

    def load(self) -> dict[str, list[dict[str, Any]]]:
        if not self.path.exists():
            return {family: [] for family in FAMILIES}
        document = json.loads(self.path.read_text(encoding="utf-8"))
        return {family: list(document.get(family, [])) for family in FAMILIES}

    def consider(self, record: dict[str, Any], envelope: Envelope) -> dict[str, Any]:
        family = record.get("family")
        if family not in FAMILIES:
            return {"promoted": False, "reason": f"unknown_family:{family}"}
        if not record.get("cross_mass_valid"):
            return {"promoted": False, "reason": "cross_mass_validation_failed"}
        if not record.get("same_X"):
            return {"promoted": False, "reason": "X_not_fixed_across_masses"}
        lock_path = self.path.with_suffix(".lock")
        with lock_path.open("a", encoding="utf-8") as lock:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
            try:
                archive = self.load()
                entries = [e for e in archive[family] if e.get("family_id") != record.get("family_id")]
                score = float(record.get("score", 0.0))
                for existing in list(entries):
                    if family_distance(existing["anchor"], record["anchor"], envelope) <= NEAR_DUPLICATE_RADIUS:
                        if float(existing.get("score", 0.0)) >= score:
                            return {"promoted": False, "reason": "near_duplicate_is_better_or_equal",
                                    "incumbent_family_id": existing.get("family_id")}
                        entries.remove(existing)
                entries.append(record)
                entries.sort(key=lambda item: (-float(item.get("score", 0.0)), str(item.get("family_id", ""))))
                evicted = entries[TOP_K:]
                archive[family] = entries[:TOP_K]
                # A record can be appended, sort last on the tie-break, and be
                # trimmed in this same call. Reporting that as promoted=True is a
                # false admission, so detect it explicitly.
                admitted = any(e.get("family_id") == record.get("family_id")
                               for e in archive[family])
                temporary = self.path.with_suffix(".json.tmp")
                temporary.write_text(json.dumps(archive, indent=2, sort_keys=True, allow_nan=False) + "\n",
                                     encoding="utf-8")
                os.replace(temporary, self.path)
                # full history is never deleted, even for evicted candidates
                with self.history_path.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps({"event": "ARCHIVE_CONSIDER", "family": family,
                                             "family_id": record.get("family_id"), "score": score,
                                             "promoted": True,
                                             "evicted": [e.get("family_id") for e in evicted]},
                                            sort_keys=True) + "\n")
            finally:
                fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        if not admitted:
            return {"promoted": False, "reason": "trimmed_by_tie_break_in_same_call",
                    "note": "score tied with incumbents; ordering fell to the family_id hash",
                    "evicted_family_ids": [e.get("family_id") for e in evicted]}
        return {"promoted": True, "reason": "promoted",
                "evicted_family_ids": [e.get("family_id") for e in evicted
                                       if e.get("family_id") != record.get("family_id")]}
