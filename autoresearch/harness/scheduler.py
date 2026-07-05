from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path

from autoresearch.harness.dihiggs_axis_contract import DIHIGGS_EXPLORERS_MODE, bin_lam1, encode_cell_id


def _mapping(value: object) -> Mapping[str, object]:
    return value if isinstance(value, Mapping) else {}


def _int_mapping(value: object) -> dict[int, int]:
    if not isinstance(value, Mapping):
        return {}
    out: dict[int, int] = {}
    for key, raw in value.items():
        try:
            out[int(str(key))] = int(raw)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            continue
    return out


def _stable_id(*parts: object, length: int = 16) -> str:
    data = "|".join(str(part) for part in parts)
    return hashlib.sha256(data.encode("utf-8")).hexdigest()[:length]


def _load_adaptive_proposals(checkpoint_root: Path) -> list[tuple[Path, int, Mapping[str, object]]]:
    proposals: list[tuple[Path, int, Mapping[str, object]]] = []
    if not checkpoint_root.exists():
        return proposals

    for iter_dir in sorted(path for path in checkpoint_root.iterdir() if path.is_dir()):
        state_file = iter_dir / "adaptive_state.json"
        if not state_file.exists():
            continue
        try:
            state = json.loads(state_file.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        raw_proposals = state.get("proposals") if isinstance(state, Mapping) else None
        if not isinstance(raw_proposals, list):
            continue
        for proposal_index, proposal in enumerate(raw_proposals):
            if isinstance(proposal, Mapping):
                proposals.append((state_file, proposal_index, proposal))
    return proposals


def _artifact_status(run_dir: str | None) -> str:
    if not run_dir:
        return "CRASH"
    run_path = Path(run_dir)
    if not run_path.exists():
        return "CRASH"
    if not (run_path / "results.csv").exists() and not (run_path / "task_summary.jsonl").exists():
        return "BROKEN_ARTIFACT"
    return "DONE"


def _event_payload(
    *,
    campaign_id: str,
    arm_id: str,
    state_file: Path,
    proposal_index: int,
    proposal: Mapping[str, object],
    tb: int,
    successes: int,
    trials: int,
    duplicate_run: bool,
    config: Mapping[str, object],
) -> dict[str, object]:
    lam1_raw = float(proposal.get("lam1_min", 0.0) or 0.0)
    axes_binned = {"tb": tb, "lam1_bin": bin_lam1(lam1_raw, config)}
    axes_raw = {"tb": tb, "lam1": lam1_raw}
    cell_id = encode_cell_id(axes_binned)
    run_dir = proposal.get("run_dir")
    run_dir_str = str(run_dir) if run_dir not in (None, "") else ""

    status = "SKIP_REUSED" if duplicate_run else _artifact_status(run_dir_str or None)
    event_successes = 0 if status in {"SKIP_REUSED", "CRASH", "BROKEN_ARTIFACT"} else successes
    event_trials = 0 if status in {"SKIP_REUSED", "CRASH", "BROKEN_ARTIFACT"} else trials
    source_ref = {
        "kind": "adaptive_state",
        "path": str(state_file),
        "proposal_index": proposal_index,
        "proposal_id": str(proposal.get("proposal_id", proposal_index)),
    }
    fingerprint = _stable_id(campaign_id, arm_id, cell_id, run_dir_str, source_ref["proposal_id"], status, tb, length=64)

    payload: dict[str, object] = {
        "arm_id": arm_id,
        "attempt_id": _stable_id(campaign_id, arm_id, cell_id, run_dir_str, source_ref["proposal_id"], tb),
        "iter_index": 0,
        "attempt_index": proposal_index,
        "cell_id": cell_id,
        "eval_status": status,
        "successes": event_successes,
        "trials": event_trials,
        "elapsed_sec": 0.0,
        "axes_binned": axes_binned,
        "axes_raw": axes_raw,
        "source_ref": source_ref,
        "fingerprint_sha256": fingerprint,
        "is_duplicate": duplicate_run,
    }
    if run_dir_str:
        payload["run_dir"] = run_dir_str
    return payload


def run_campaign(
    *,
    config: dict[str, object],
    mode: str,
    iters: int,
    attempts_per_iter: int,
) -> dict[str, object]:
    if mode != DIHIGGS_EXPLORERS_MODE:
        raise ValueError(f"Unsupported campaign mode: {mode}")

    paths = _mapping(config.get("paths"))
    outdir = Path(str(paths.get("outdir", "autoresearch-out")))
    lake_name = str(paths.get("lake_name", "events.jsonl"))
    outdir.mkdir(parents=True, exist_ok=True)
    reports_dir = outdir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    events_path = outdir / lake_name

    arms = config.get("arms")
    if not isinstance(arms, list) or not arms or not isinstance(arms[0], Mapping):
        raise ValueError("config.arms must contain at least one arm")
    arm_id = str(arms[0].get("id", "adaptive-v1"))

    checkpoint_root = outdir / "checkpoints" / arm_id
    proposals = _load_adaptive_proposals(checkpoint_root)
    seen_run_dirs: set[str] = set()
    events: list[dict[str, object]] = []

    for state_file, proposal_index, proposal in proposals:
        successes_by_tb = _int_mapping(proposal.get("successes_by_tb"))
        trials_by_tb = _int_mapping(proposal.get("trials_by_tb"))
        tb_values = sorted(set(successes_by_tb) | set(trials_by_tb))
        if not tb_values:
            tb_values = [int(proposal.get("tb", proposal.get("bin_index", 0)) or 0)]

        run_dir_obj = proposal.get("run_dir")
        run_dir = str(run_dir_obj) if run_dir_obj not in (None, "") else ""
        duplicate_run = bool(run_dir and run_dir in seen_run_dirs)
        for tb in tb_values:
            payload = _event_payload(
                campaign_id=str(config.get("campaign_id", "")),
                arm_id=arm_id,
                state_file=state_file,
                proposal_index=proposal_index,
                proposal=proposal,
                tb=tb,
                successes=successes_by_tb.get(tb, 0),
                trials=trials_by_tb.get(tb, 0),
                duplicate_run=duplicate_run,
                config=config,
            )
            events.append(
                {
                    "schema_version": 1,
                    "campaign_id": config.get("campaign_id"),
                    "event_type": "ATTEMPT_EVALUATED",
                    "utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
                    "payload": payload,
                }
            )
        if run_dir:
            seen_run_dirs.add(run_dir)

    events_path.write_text("".join(json.dumps(event) + "\n" for event in events), encoding="utf-8")
    summary = {
        "status": "pass",
        "arm_id": arm_id,
        "discoveries": len(proposals),
        "attempts": len(events),
        "events_path": str(events_path),
    }
    (reports_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary
