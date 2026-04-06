"""
Checkpoint parser for adaptive_explorer.py discoveries.

Extracts discovered run_dir paths from checkpoint files, mapping them to axes
for metrics computation.
"""

from __future__ import annotations

import json
import logging
import os
import time
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from autoresearch.harness.dihiggs_axis_contract import bin_lam1


logger = logging.getLogger(__name__)


@dataclass
class AdaptiveDiscovery:
    """Represents a single discovered run from adaptive explorer checkpoints."""

    run_dir: str
    tb: int
    lam1_raw: float
    lam1_bin: int
    iter_index: int
    proposal_index: int


@dataclass
class BranchDiscovery:
    """Represents a single discovered run from branch continuation checkpoints."""

    run_dir: str
    tb: int
    lam1_raw: float
    lam1_bin: int
    track_id: str
    step_name: str


def parse_adaptive_checkpoint(
    checkpoint_root: str, config: Mapping[str, object], known_run_dirs: set[str]
) -> list[AdaptiveDiscovery]:
    """
    Parse adaptive explorer checkpoints to extract NEW run directories.

    Scans checkpoint_root/iter_XXXX/adaptive_state.json files, extracts proposals
    with run_dir, bins lam1, and returns only NEW discoveries (not in known_run_dirs).

    Args:
        checkpoint_root: Root directory containing iter_XXXX subdirectories
        config: Configuration dict with search.lam1 binning parameters
        known_run_dirs: Set of already-known run directories to filter out

    Returns:
        List of AdaptiveDiscovery instances sorted by (iter_index, proposal_index)

    Notes:
        - Only parses checkpoints stable for ≥2 seconds (avoids half-written files)
        - Handles missing/malformed checkpoints gracefully (logs warning, continues)
        - Filters out proposals with null/empty run_dir or already in known_run_dirs
    """
    checkpoint_path = Path(checkpoint_root)
    if not checkpoint_path.exists():
        logger.warning(f"Checkpoint root does not exist: {checkpoint_root}")
        return []

    discoveries: list[AdaptiveDiscovery] = []
    now = time.time()
    stability_threshold = 2.0  # seconds

    # Find all iter_XXXX directories
    iter_dirs: list[tuple[int, Path]] = []
    for entry in checkpoint_path.iterdir():
        if not entry.is_dir():
            continue
        # Match pattern: iter_XXXX where XXXX is 4-digit number
        name = entry.name
        if name.startswith("iter_") and len(name) == 9:
            try:
                iter_num = int(name[5:])
                iter_dirs.append((iter_num, entry))
            except ValueError:
                continue

    # Sort by iteration number
    iter_dirs.sort(key=lambda x: x[0])

    for iter_index, iter_dir in iter_dirs:
        state_file = iter_dir / "adaptive_state.json"
        if not state_file.exists():
            logger.debug(f"Checkpoint file not found: {state_file}")
            continue

        # Check stability: skip if modified within last 2 seconds
        try:
            mtime = os.stat(state_file).st_mtime
            if mtime > (now - stability_threshold):
                logger.debug(
                    f"Skipping unstable checkpoint (recently modified): {state_file}"
                )
                continue
        except OSError as exc:
            logger.warning(f"Failed to stat checkpoint file {state_file}: {exc}")
            continue

        # Load and parse checkpoint
        try:
            with open(state_file, "r", encoding="utf-8") as f:
                state: dict[str, Any] = json.load(f)
        except (OSError, json.JSONDecodeError) as exc:
            logger.warning(f"Failed to load checkpoint {state_file}: {exc}")
            continue

        if not isinstance(state, dict):
            logger.warning(f"Checkpoint {state_file} is not a JSON object")
            continue

        # Extract proposals
        proposals = state.get("proposals", [])
        if not isinstance(proposals, list):
            logger.warning(
                f"Checkpoint {state_file} has invalid proposals (not a list)"
            )
            continue

        for proposal_data in proposals:
            if not isinstance(proposal_data, dict):
                logger.debug(f"Skipping non-dict proposal in {state_file}")
                continue

            # Extract required fields
            run_dir = proposal_data.get("run_dir")
            if run_dir is None or run_dir == "":
                # Proposal without run_dir (e.g., failed orchestrator call)
                continue

            run_dir_str = str(run_dir)
            if run_dir_str in known_run_dirs:
                # Already processed
                continue

            # Extract other fields with fallback defaults
            try:
                tb = int(proposal_data.get("bin_index", 0))
            except (TypeError, ValueError):
                logger.debug(
                    f"Skipping proposal in {state_file}: invalid bin_index (tb)"
                )
                continue

            try:
                lam1_raw = float(proposal_data.get("lam1_min", 0.0))
            except (TypeError, ValueError):
                logger.debug(
                    f"Skipping proposal in {state_file}: invalid lam1_min"
                )
                continue

            try:
                proposal_index = int(proposal_data.get("bin_index", 0))
            except (TypeError, ValueError):
                logger.debug(
                    f"Skipping proposal in {state_file}: invalid bin_index (proposal_index)"
                )
                continue

            # Bin lam1 using axis contract
            try:
                lam1_bin = bin_lam1(lam1_raw, config)
            except (KeyError, TypeError, ValueError) as exc:
                logger.warning(
                    f"Failed to bin lam1={lam1_raw} for proposal in {state_file}: {exc}"
                )
                continue

            # Create discovery
            discovery = AdaptiveDiscovery(
                run_dir=run_dir_str,
                tb=tb,
                lam1_raw=lam1_raw,
                lam1_bin=lam1_bin,
                iter_index=iter_index,
                proposal_index=proposal_index,
            )
            discoveries.append(discovery)

    # Sort by (iter_index, proposal_index) for deterministic ordering
    discoveries.sort(key=lambda d: (d.iter_index, d.proposal_index))

    return discoveries


def parse_branch_checkpoint(
    checkpoint_root: str, config: Mapping[str, object], known_run_dirs: set[str]
) -> list[BranchDiscovery]:
    """
    Parse branch continuation explorer checkpoints to extract NEW run directories.

    Scans checkpoint_root/{track_id}/events.jsonl files, extracts ATTEMPT_COMPLETED
    events with run_dir, bins lam1, and returns only NEW discoveries (not in known_run_dirs).

    Args:
        checkpoint_root: Root directory containing track_id subdirectories
        config: Configuration dict with search.lam1 binning parameters
        known_run_dirs: Set of already-known run directories to filter out

    Returns:
        List of BranchDiscovery instances sorted by (track_id, step_name)

    Notes:
        - JSONL format: one JSON object per line
        - Only processes events with event_type == "ATTEMPT_COMPLETED"
        - Extracts: run_dir from payload.result.run_dir, tb from payload.params.tb,
          lam1 from payload.params.lam1, step_name from payload.step_label
        - Handles missing/malformed checkpoints gracefully (logs warning, continues)
        - Filters out events with null/empty run_dir or already in known_run_dirs
    """
    checkpoint_path = Path(checkpoint_root)
    if not checkpoint_path.exists():
        logger.warning(f"Checkpoint root does not exist: {checkpoint_root}")
        return []

    discoveries: list[BranchDiscovery] = []
    now = time.time()
    stability_threshold = 2.0  # seconds

    # Find all track_id subdirectories, sort alphabetically
    track_dirs: list[tuple[str, Path]] = []
    for entry in checkpoint_path.iterdir():
        if entry.is_dir():
            track_dirs.append((entry.name, entry))

    track_dirs.sort(key=lambda x: x[0])

    for track_id, track_dir in track_dirs:
        events_jsonl = track_dir / "events.jsonl"
        if not events_jsonl.exists():
            logger.debug(f"Events file not found: {events_jsonl}")
            continue

        # Check stability: skip if modified within last 2 seconds
        try:
            mtime = os.stat(events_jsonl).st_mtime
            if mtime > (now - stability_threshold):
                logger.debug(
                    f"Skipping unstable events file (recently modified): {events_jsonl}"
                )
                continue
        except OSError as exc:
            logger.warning(f"Failed to stat events file {events_jsonl}: {exc}")
            continue

        # Read JSONL line by line
        try:
            with open(events_jsonl, "r", encoding="utf-8") as f:
                for line_num, line in enumerate(f, start=1):
                    line = line.strip()
                    if not line:
                        continue

                    try:
                        event: dict[str, Any] = json.loads(line)
                    except json.JSONDecodeError as exc:
                        logger.warning(
                            f"Malformed JSON at {events_jsonl}:{line_num}: {exc}"
                        )
                        continue

                    if not isinstance(event, dict):
                        logger.debug(
                            f"Skipping non-dict event at {events_jsonl}:{line_num}"
                        )
                        continue

                    # Filter: only ATTEMPT_COMPLETED events
                    event_type = event.get("event_type")
                    if event_type != "ATTEMPT_COMPLETED":
                        continue

                    # Extract payload
                    payload = event.get("payload")
                    if not isinstance(payload, dict):
                        logger.debug(
                            f"Skipping event with invalid payload at {events_jsonl}:{line_num}"
                        )
                        continue

                    # Extract result.run_dir
                    result = payload.get("result")
                    if not isinstance(result, dict):
                        logger.debug(
                            f"Skipping event with invalid result at {events_jsonl}:{line_num}"
                        )
                        continue

                    run_dir = result.get("run_dir")
                    if run_dir is None or run_dir == "":
                        continue

                    run_dir_str = str(run_dir)
                    if run_dir_str in known_run_dirs:
                        continue

                    # Extract params
                    params = payload.get("params")
                    if not isinstance(params, dict):
                        logger.debug(
                            f"Skipping event with invalid params at {events_jsonl}:{line_num}"
                        )
                        continue

                    # Extract tb (required)
                    if "tb" not in params:
                        logger.debug(
                            f"Skipping event at {events_jsonl}:{line_num}: missing tb"
                        )
                        continue
                    try:
                        tb = int(params["tb"])
                    except (TypeError, ValueError):
                        logger.debug(
                            f"Skipping event at {events_jsonl}:{line_num}: invalid tb"
                        )
                        continue

                    # Extract lam1 (required)
                    if "lam1" not in params:
                        logger.debug(
                            f"Skipping event at {events_jsonl}:{line_num}: missing lam1"
                        )
                        continue
                    try:
                        lam1_raw = float(params["lam1"])
                    except (TypeError, ValueError):
                        logger.debug(
                            f"Skipping event at {events_jsonl}:{line_num}: invalid lam1"
                        )
                        continue

                    # Extract step_name (required)
                    step_name = payload.get("step_label")
                    if step_name is None:
                        logger.debug(
                            f"Skipping event at {events_jsonl}:{line_num}: missing step_label"
                        )
                        continue
                    step_name_str = str(step_name)

                    # Bin lam1
                    try:
                        lam1_bin = bin_lam1(lam1_raw, config)
                    except (KeyError, TypeError, ValueError) as exc:
                        logger.warning(
                            f"Failed to bin lam1={lam1_raw} at {events_jsonl}:{line_num}: {exc}"
                        )
                        continue

                    # Create discovery
                    discovery = BranchDiscovery(
                        run_dir=run_dir_str,
                        tb=tb,
                        lam1_raw=lam1_raw,
                        lam1_bin=lam1_bin,
                        track_id=track_id,
                        step_name=step_name_str,
                    )
                    discoveries.append(discovery)

        except OSError as exc:
            logger.warning(f"Failed to read events file {events_jsonl}: {exc}")
            continue

    # Sort by (track_id, step_name) for deterministic ordering
    discoveries.sort(key=lambda d: (d.track_id, d.step_name))

    return discoveries
