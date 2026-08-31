"""Rebuild the Top-5-per-family discovery archive by replaying the ledger.

`discovery_archive.json` is a derived cache written by `DiscoveryArchive.consider`
(archive.py); the ledger is authoritative. This proves and exercises that the
cache is reconstructible: replay every `FAMILY_VALIDATION` event, in ledger
order, through the same `consider()` admission logic, and the result should
match the live archive exactly (mission acceptance test T4).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from search_substrate.ledger import Ledger

from .archive import FAMILIES, DiscoveryArchive
from .daemon import guard_run_dir
from .envelopes import active_envelope

# Keys Ledger.append adds on top of the original record; strip them so a
# replayed record is byte-identical to what cmd_validate_family/daemon
# originally passed into archive.consider().
_LEDGER_META_KEYS = {"event", "lifecycle", "timestamp_unix", "event_id"}


def rebuild(run_dir: Path, target_dir: Path | None = None) -> dict[str, list[dict[str, Any]]]:
    """Replay ledger FAMILY_VALIDATION events into a fresh archive.

    Pass `target_dir` (different from `run_dir`) to rebuild into a scratch
    location without touching the live `discovery_archive.json`.
    """
    run_dir = Path(run_dir)
    envelope = active_envelope(run_dir)
    ledger = Ledger(run_dir)
    destination = Path(target_dir) if target_dir is not None else run_dir
    # Refuse to write into the live, review-blocked campaign directories or
    # deliverables/ (reuses daemon.py's guard, not reimplemented here).
    # Reading `run_dir`'s ledger is unaffected -- that is read-only replay,
    # not a write, and reconstructing from a live/protected ledger is exactly
    # the reference use case (mission T4); only the write *destination* is
    # guarded. Point --target-dir elsewhere to rebuild from a protected
    # run_dir's ledger without touching it.
    guard_run_dir(destination)
    fresh_archive = DiscoveryArchive(destination)
    if fresh_archive.path.exists():
        fresh_archive.path.unlink()  # discovery_archive.json is a derived cache, safe to reset
    # discovery_archive_history.jsonl is append-only eviction evidence
    # (archive.py: "never deleted, even for evicted candidates") and is never
    # unlinked by rebuild, even when destination == run_dir; replay below can
    # only ever append further ARCHIVE_CONSIDER rows on top of it via the
    # same DiscoveryArchive.consider() used by the live campaign.
    for event in ledger.events():
        if event.get("event") != "FAMILY_VALIDATION":
            continue
        if not (event.get("cross_mass_valid") and event.get("same_X")):
            continue
        if event.get("family") not in FAMILIES:
            continue
        record = {k: v for k, v in event.items() if k not in _LEDGER_META_KEYS}
        fresh_archive.consider(record, envelope)
    return fresh_archive.load()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Rebuild the discovery archive from the ledger alone")
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--target-dir", default=None,
                        help="write the rebuilt archive here instead of overwriting the live one")
    args = parser.parse_args(argv)
    target = Path(args.target_dir) if args.target_dir else None
    result = rebuild(Path(args.run_dir), target)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
