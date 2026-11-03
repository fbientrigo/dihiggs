#!/usr/bin/env python3
"""ScientificState -> Jev decision -> allowlisted template -> canonical V2 plan -> Lambda1EvaluatorV2.

Canonical vertical workflow:
1. Jev policy selects an allowlisted template based on compact ScientificState.
2. Canonical adapter decomposes template into fixed (mA, lambda6) subruns for Lambda1EvaluatorV2.
3. Plan is materialized under RUNS_ROOT / "jev_harness" (never in the source tree).
4. If --execute is specified: subruns execute via Lambda1EvaluatorV2, producing dihiggs.lambda1.v2
   CSVs; outcomes are extracted (ctau_mm, valid rows, best point); worker state is updated.
5. If --execute is omitted (default): plan-only mode verifies row counts and creates input CSVs.
"""
from __future__ import annotations

import argparse
import json
import sys
import uuid
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.jev_adapter import (
    CANONICAL_ENGINE,
    CANONICAL_PRODUCER,
    CANONICAL_SCHEMA,
    CanonicalPlan,
    execute_canonical_plan,
    materialize_canonical_plan,
    now_iso,
    plan_template_v2,
    update_scientific_state,
)
from scripts.jev_decide_once import (
    DEFAULT_EVENTS,
    RUNS_ROOT,
    append_jsonl,
    baseline_state,
    default_state_path,
    load_object,
)
from scripts.jev_policy import (
    JevDecision,
    build_decision_event,
    build_jev_request,
    choose_with_gemini,
    choose_with_jev,
    response_to_worker_action,
)

HARNESS_DIR = RUNS_ROOT / "jev_harness"


def materialize_and_plancheck(
    decision_action: dict[str, Any],
    cfg: dict[str, Any],
    *,
    iteration: int,
    cycle_id: str,
    iteration_id: str,
    output_root: Path | None = None,
) -> dict[str, Any]:
    """Materialize the Jev-selected template into canonical V2 subruns (plan-only)."""
    harness_dir = output_root or HARNESS_DIR
    plan = plan_template_v2(
        template_name=str(decision_action["template_name"]),
        cfg=cfg,
        iteration=iteration,
        cycle_id=cycle_id,
        iteration_id=iteration_id,
        overrides=decision_action.get("template_overrides") or {},
        output_root=harness_dir,
    )

    plan_dir = harness_dir / f"plancheck_{cycle_id}_{iteration_id}"
    summary = materialize_canonical_plan(plan, plan_dir, dry_run=True)
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run one Jev decision through canonical Lambda1EvaluatorV2 planning or execution"
    )
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--state", type=Path, default=None)
    ap.add_argument("--events", type=Path, default=DEFAULT_EVENTS)
    ap.add_argument("--execute", action="store_true", help="Execute canonical Lambda1EvaluatorV2 subruns")
    ap.add_argument("--exec-path", type=Path, default=None, help="Override Lambda1EvaluatorV2 binary path")
    ap.add_argument("--timeout", type=float, default=None, help="Per-subrun timeout in seconds")
    ap.add_argument("--force", action="store_true", help="Force overwrite existing run directory")
    ap.add_argument(
        "--decision-fixture",
        type=Path,
        default=None,
        help="Load a saved raw Jev response from this JSON file instead of calling Jev. "
        "Makes no network request; the fixture is parsed through the same "
        "response_to_worker_action path a live call would use.",
    )
    ap.add_argument(
        "--policy",
        choices=["jev", "gemini"],
        default="jev",
        help="Decision backend when --decision-fixture is not given. 'gemini' uses the "
        "local `gemini` CLI (free) answering the same payload Jev would get; 'jev' "
        "makes a real paid Jev call. Mutually exclusive with --decision-fixture.",
    )
    args = ap.parse_args()
    if args.decision_fixture is not None and args.policy != "jev":
        ap.error("--policy is mutually exclusive with --decision-fixture")

    cfg = load_object(args.config)
    state_path = args.state if args.state is not None else default_state_path()
    state = load_object(state_path) if state_path.exists() else baseline_state(cfg)

    if args.decision_fixture is not None:
        raw_response = load_object(args.decision_fixture)
        decision = JevDecision(
            worker_action=response_to_worker_action(raw_response),
            request_payload=build_jev_request(cfg, state),
            raw_response=raw_response,
        )
    elif args.policy == "gemini":
        decision = choose_with_gemini(cfg, state)
    else:
        decision = choose_with_jev(cfg, state)
    iteration = int(state.get("iterations_completed", 0) or 0) + 1
    cycle_id = str(state.get("cycle_id", "standalone_jev"))
    iteration_id = f"jev_iter{iteration:03d}_{uuid.uuid4().hex[:6]}"

    event = build_decision_event(
        decision,
        timestamp=now_iso(),
        iteration=iteration,
        iteration_id=iteration_id,
        cycle_id=cycle_id,
    )
    if args.decision_fixture is not None:
        event["backend"] = "offline_fixture"
        event["decision_fixture_path"] = str(args.decision_fixture)
    else:
        event["backend"] = "gemini_cli" if args.policy == "gemini" else "jev"
    try:
        from dihiggs.app.orchestrator.physics_authority import convention_provenance
        event["mass_convention"] = convention_provenance(source_commit=None)
    except Exception:
        pass

    if decision.worker_action.get("action") == "stop":
        event["outcome"] = {"status": "stopped_by_jev"}
    elif args.execute:
        try:
            plan = plan_template_v2(
                template_name=str(decision.worker_action["template_name"]),
                cfg=cfg,
                iteration=iteration,
                cycle_id=cycle_id,
                iteration_id=iteration_id,
                overrides=decision.worker_action.get("template_overrides") or {},
                output_root=HARNESS_DIR,
            )
            run_dir = HARNESS_DIR / f"run_{cycle_id}_{iteration_id}"
            execution_res = execute_canonical_plan(
                plan=plan,
                plan_dir=run_dir,
                executable=args.exec_path,
                force=args.force,
                timeout_s=args.timeout,
            )
            event["outcome"] = execution_res["outcome"]

            updated_state = update_scientific_state(
                current_state=state,
                config=cfg,
                outcome=execution_res["outcome"],
                iteration=iteration,
                cycle_id=cycle_id,
                template_name=plan.template_name,
            )
            state_path.parent.mkdir(parents=True, exist_ok=True)
            state_path.write_text(json.dumps(updated_state, indent=2) + "\n", encoding="utf-8")
        except Exception as exc:
            event["outcome"] = {"status": "execution_failed", "error": str(exc)}
    else:
        # Default plan-only mode
        try:
            outcome = materialize_and_plancheck(
                decision.worker_action,
                cfg,
                iteration=iteration,
                cycle_id=cycle_id,
                iteration_id=iteration_id,
                output_root=HARNESS_DIR,
            )
            event["outcome"] = outcome
        except Exception as exc:
            event["outcome"] = {"status": "materialization_failed", "error": str(exc)}

    append_jsonl(args.events, event)
    print(
        json.dumps(
            {
                "selected": decision.selected,
                "outcome": event["outcome"],
                "event_written_to": str(args.events),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
