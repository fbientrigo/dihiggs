#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.jev_policy import build_decision_event, build_jev_request, choose_with_jev


WORKSPACE = ROOT.parent
RUNS_ROOT = WORKSPACE / "runs"
PREFERRED_STATE = RUNS_ROOT / "jev_harness" / "worker_state.json"
FALLBACK_STATE = RUNS_ROOT / "gemini_harness" / "worker_state.json"
DEFAULT_EVENTS = RUNS_ROOT / "jev_harness" / "decision_events.jsonl"


def default_state_path() -> Path:
    """Prefer external workspace state for jev_harness; fall back to gemini_harness if present."""
    if PREFERRED_STATE.exists():
        return PREFERRED_STATE
    if FALLBACK_STATE.exists():
        return FALLBACK_STATE
    return PREFERRED_STATE


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def load_object(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return obj


def baseline_state(config: dict[str, Any]) -> dict[str, Any]:
    # *_mm is authoritative when present; *_m is always derived (mm/1000), never
    # read as stored, so a stale/mismatched *_m can never silently win (factor-1000 guard).
    if "baseline_best_ctau_mm" in config:
        base_mm = float(config["baseline_best_ctau_mm"])
        base_m = base_mm / 1000.0
    else:
        base_m = float(config["baseline_best_ctau_m"])
        base_mm = base_m * 1000.0
    return {
        "status": "jev_baseline",
        "best_ctau_mm": base_mm,
        "best_ctau_m": base_m,
        "executed_points_total": 0,
        "iterations_completed": 0,
        "no_improvement_count": 0,
        "rejected_points_total": 0,
        "contaminated_summary_count": 0,
        "history": [],
        "last_iteration_summary_for_prompt": None,
    }


def append_jsonl(path: Path, obj: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(obj, sort_keys=True) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description="Make one bounded DiHiggs decision with Jev via Vercel AI Gateway")
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--state", type=Path, default=None)
    ap.add_argument("--events", type=Path, default=DEFAULT_EVENTS)
    ap.add_argument("--request-only", action="store_true", help="Print the System One request without calling Jev")
    args = ap.parse_args()

    cfg = load_object(args.config)
    state_path = args.state if args.state is not None else default_state_path()
    state = load_object(state_path) if state_path.exists() else baseline_state(cfg)

    if args.request_only:
        print(json.dumps(build_jev_request(cfg, state), indent=2, sort_keys=True))
        return 0

    decision = choose_with_jev(cfg, state)
    iteration = int(state.get("iterations_completed", 0) or 0) + 1
    event = build_decision_event(
        decision,
        timestamp=now_iso(),
        iteration=iteration,
        iteration_id=f"jev_decision_{iteration:03d}",
        cycle_id=str(state.get("cycle_id", "standalone_jev")),
    )
    try:
        from dihiggs.app.orchestrator.physics_authority import convention_provenance
        event["mass_convention"] = convention_provenance(source_commit=None)
    except Exception:
        pass
    append_jsonl(args.events, event)

    output = {
        "selected": decision.selected,
        "worker_action": decision.worker_action,
        "model": decision.raw_response.get("model", "typesafe-ai/jev"),
        "answers": decision.raw_response.get("answers"),
        "usage": decision.raw_response.get("usage"),
        "event_written_to": str(args.events),
    }
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
