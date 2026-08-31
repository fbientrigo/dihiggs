from __future__ import annotations

import json
import fcntl
import math
import os
from pathlib import Path
from typing import Any

from .contract import BOUNDS


FAMILIES = ("mixed", "photonic")
TOP_K = 5
NEAR_DUPLICATE_RADIUS = 0.05


def _distance(left: dict[str, Any], right: dict[str, Any]) -> float:
    names = ("mH_GeV", "mA_GeV", "M2_GeV2", "tan_beta", "lambda6")
    values = []
    for name in names:
        low, high = BOUNDS[name]
        a = left["parameters"][name]
        b = right["parameters"][name]
        values.append(((a - b) / (high - low)) ** 2)
    return math.sqrt(sum(values) / len(values))


class Archive:
    """Code-owned derived Top-K archive; callers cannot supply archive state."""

    def __init__(self, run_dir: Path):
        self.run_dir = run_dir
        self.path = run_dir / "archive.json"
        self.run_dir.mkdir(parents=True, exist_ok=True)

    def load(self) -> dict[str, list[dict[str, Any]]]:
        if not self.path.exists():
            return {family: [] for family in FAMILIES}
        document = json.loads(self.path.read_text(encoding="utf-8"))
        return {family: list(document.get(family, [])) for family in FAMILIES}

    def consider(self, family_record: dict[str, Any]) -> dict[str, Any]:
        family = family_record.get("family")
        if family not in FAMILIES:
            return {"promoted": False, "reason": "unknown_family"}
        if family_record.get("status") != "PROMOTABLE" or not family_record.get("validity_gate"):
            return {"promoted": False, "reason": "not_promotable"}
        lock_path = self.path.with_suffix(".lock")
        with lock_path.open("a", encoding="utf-8") as lock:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
            try:
                archive = self.load()
                entries = archive[family]
                candidate_id = family_record.get("family_id")
                duplicate = next((entry for entry in entries if entry.get("family_id") == candidate_id), None)
                if duplicate is not None and float(duplicate.get("total_score", 0.0)) >= float(family_record.get("total_score", 0.0)):
                    return {"promoted": False, "reason": "existing_family_is_better_or_equal"}
                entries = [entry for entry in entries if entry.get("family_id") != candidate_id]
                for entry in list(entries):
                    if _distance(entry["anchor"], family_record["anchor"]) <= NEAR_DUPLICATE_RADIUS:
                        if float(entry.get("total_score", 0.0)) >= float(family_record.get("total_score", 0.0)):
                            return {"promoted": False, "reason": "near_duplicate_is_better_or_equal"}
                        entries.remove(entry)
                entries.append(family_record)
                entries.sort(key=lambda item: (-float(item.get("total_score", 0.0)), str(item.get("family_id", ""))))
                evicted = entries[TOP_K:]
                archive[family] = entries[:TOP_K]
                temporary = self.path.with_suffix(".json.tmp")
                temporary.write_text(json.dumps(archive, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
                os.replace(temporary, self.path)
            finally:
                fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        return {"promoted": True, "reason": "promoted", "evicted_family_ids": [item.get("family_id") for item in evicted]}
