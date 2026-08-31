"""Bounded, resumable background harvesting loop over the frozen discovery substrate.

    resume state -> propose -> dedup -> evaluate -> append evidence
    -> update archive -> checkpoint -> repeat

This module never reimplements physics, the ledger, the archive, or the
proposal strategies -- it imports and drives them. See
IMPLEMENTATION_DECISION.md for the reuse/extension map.

Concurrency model: a worker-pool process evaluates one candidate each
(the frozen evaluator subprocess call is the unit of work); every
`Ledger.append` happens back in this single daemon process after a worker's
result is collected, so there is exactly one ledger writer regardless of
`--workers`. `climb`/`continue_lineage` strategies are the one documented
exception -- they are pre-existing self-contained ask-evaluate-tell loops
that already call `Ledger.append` internally, and the daemon runs them
synchronously in its own process (never dispatched to the pool), so they
never overlap with a worker-pool write either. See policy_adapter.py.
"""
from __future__ import annotations

import concurrent.futures as cf
import hashlib
import json
import multiprocessing
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import yaml

from search_substrate.ledger import Ledger

from . import daemon_checkpoint as checkpoint_mod
from . import policy_adapter as pa
from .archive import DiscoveryArchive
from .cli import _seen_ids as _replay_seen_ids
from .envelopes import ENVELOPES, active_envelope, save_envelope
from .evaluator import DiscoveryEvaluator
from .family import validate_family
from .helpers import deduplicate as _dedupe
from .objective import lifetime_distance

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVALUATOR = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"

# Live campaign state under independent review (RECONNAISSANCE_REPORT.md /
# IMPLEMENTATION_DECISION.md). The daemon must never write here.
PROTECTED_RUN_DIR_NAMES = {"discovery_v1", "discovery_photonic_v1"}
PROTECTED_PATH_NAMES = {"deliverables"}

SIGNAL_EXIT_CODES = {signal.SIGINT: 130, signal.SIGTERM: 143}


class DaemonConfigError(ValueError):
    pass


# --------------------------------------------------------------- guard rails
def guard_run_dir(run_dir: Path) -> None:
    """Refuse to operate on the live, review-blocked campaign directories or
    deliverables/, including any subdirectory of them -- path containment
    against every ancestor, not just an exact match on run_dir's own name."""
    resolved = Path(run_dir).resolve()
    for node in (resolved, *resolved.parents):
        if node.name in PROTECTED_RUN_DIR_NAMES:
            raise DaemonConfigError(f"refusing_protected_run_dir:{resolved}_under:{node}")
    for node in (resolved, *resolved.parents):
        if node.name in PROTECTED_PATH_NAMES:
            raise DaemonConfigError(f"refusing_protected_path:{resolved}_under:{node}")


def guard_write_path(path: Path, run_dir: Path) -> None:
    run_dir_resolved = Path(run_dir).resolve()
    resolved = Path(path).resolve()
    try:
        resolved.relative_to(run_dir_resolved)
    except ValueError as error:
        raise DaemonConfigError(
            f"daemon_write_outside_run_dir:{resolved}:not_under:{run_dir_resolved}"
        ) from error


# -------------------------------------------------------------------- config
DEFAULTS: dict[str, Any] = {
    "campaign": {"name": "unnamed_harvest"},
    "run_dir": "runs/harvest_v1",
    "evaluator": str(DEFAULT_EVALUATOR),
    "envelope_id": "E1_mixed_exploit",
    "runtime": {
        "workers": 2, "max_evaluations_per_cycle": 20,
        "checkpoint_every": 1, "summary_every": 1, "shutdown_grace_s": 30,
    },
    "budgets": {
        "max_total_evaluations": None, "max_cycle_walltime_s": 900, "max_total_walltime_s": None,
    },
    "policy": {
        "allocation": {"explore": 1.0}, "family_validate_per_cycle": 1, "seed": 20260825,
    },
    "stopping": {"patience_cycles": 20, "max_consecutive_evaluator_failures": 5},
    "safety": {"downstream_snapshots_read_only": True},
}


def _merge(base: dict[str, Any], override: dict[str, Any] | None) -> dict[str, Any]:
    out = dict(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _merge(out[key], value)
        else:
            out[key] = value
    return out


def load_config(path: Path) -> dict[str, Any]:
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    config = _merge(DEFAULTS, raw)
    if config["envelope_id"] not in ENVELOPES:
        raise DaemonConfigError(f"unknown_envelope_id:{config['envelope_id']}")
    pa.resolve_allocation(config["policy"]["allocation"])  # validates strategy names, raises if bad
    return config


def _config_digest(config: dict[str, Any]) -> str:
    material = {"envelope_id": config["envelope_id"], "policy": config["policy"]}
    return hashlib.sha256(json.dumps(material, sort_keys=True, default=str).encode()).hexdigest()[:16]


# ---------------------------------------------------------- crash recovery
def repair_torn_ledger_tail(ledger_path: Path) -> bool:
    """Recover from a hard kill mid-append.

    Ledger.append (frozen, search_substrate/ledger.py) always writes
    `json.dumps(record) + "\\n"` then flush+fsync as a single write, so a
    fully committed row always ends with a newline. If the file's last byte
    is not a newline, the very last write was torn by a hard kill/crash and
    was never a successfully committed record -- truncate back to the last
    complete line. Every complete line (every record that finished writing)
    is left byte-for-byte untouched; this only discards bytes that were never
    valid append-only evidence to begin with. Returns True if a repair
    happened.
    """
    if not ledger_path.exists():
        return False
    raw = ledger_path.read_bytes()
    if not raw or raw.endswith(b"\n"):
        return False
    last_newline = raw.rfind(b"\n")
    tail = raw[last_newline + 1:]
    try:
        if tail.strip():
            json.loads(tail.decode("utf-8"))
        return False
    except (UnicodeDecodeError, json.JSONDecodeError):
        with ledger_path.open("r+b") as handle:
            handle.truncate(last_newline + 1)
        return True


# -------------------------------------------------------------- diagnostics
def _best_by_family(ledger: Ledger) -> dict[str, Any]:
    best: dict[str, Any] = {"mixed": None, "photonic": None}
    for event in ledger.events():
        if event.get("event") != "EVALUATION" or not event.get("validity_gate"):
            continue
        family = (event.get("provisional_family") or {}).get("family")
        if family not in ("mixed", "photonic"):
            continue
        distance = lifetime_distance(event.get("ctau_mm"))
        current = best[family]
        if current is None or distance < current["_distance"]:
            best[family] = {"candidate_id": event.get("candidate_id"),
                            "coordinates": event.get("coordinates"), "_distance": distance}
    return {family: ({k: v for k, v in row.items() if k != "_distance"} if row else None)
            for family, row in best.items()}


def _print_full_diagnostics(run_dir: Path) -> None:
    """Reuse cli.cmd_checkpoint's exact diagnostic computation (valid_fraction,
    duplicate_rate, gate_failure_counts, basin diversity, ...) via the existing
    CLI entry point, at summary cadence, rather than reimplementing it."""
    proc = subprocess.run(
        [sys.executable, "-m", "search_discovery", "checkpoint", "--run-dir", str(run_dir)],
        cwd=ROOT, capture_output=True, text=True, check=False,
    )
    if proc.returncode != 0:
        print(json.dumps({"diagnostics_error": proc.stderr.strip()[-2000:]}))
        return
    print(proc.stdout.strip())


# ---------------------------------------------------------------- evaluation
def _evaluate_worker(root: str, executable: str, proposal: dict[str, Any],
                      attempt_dir: str, envelope) -> dict[str, Any]:
    ev = DiscoveryEvaluator(Path(root), Path(executable))
    return ev.evaluate(proposal, Path(attempt_dir), envelope)


def _wait_batch(futures: list[cf.Future], stop_flag: dict[str, Any], grace_s: float,
                cycle_deadline: float | None) -> list[dict[str, Any]]:
    """Collect worker results, bounded by whichever deadline applies.

    During graceful shutdown (`stop_flag["stop"]`), a future that is still
    running after `grace_s` is presumed cancellable and reported as a FAILED
    shutdown-timeout row (existing behavior).

    Otherwise, `cycle_deadline` (derived from `budgets.max_cycle_walltime_s`)
    bounds how long this cycle will wait on outstanding futures at all. A
    future that has not finished by then may still be genuinely running (a
    slow, not hung, evaluation) -- fabricating a FAILED ledger row for it
    would misrepresent evidence that might still arrive, so it is simply left
    unresolved (not appended, not retried this cycle) and the cycle proceeds
    to checkpoint/next cycle instead of hanging indefinitely.
    """
    results = []
    shutdown_deadline = time.time() + grace_s
    for future in futures:
        if stop_flag["stop"]:
            remaining = max(0.0, shutdown_deadline - time.time())
        elif cycle_deadline is not None:
            remaining = max(0.0, cycle_deadline - time.time())
        else:
            remaining = None
        try:
            results.append(future.result(timeout=remaining))
        except cf.TimeoutError:
            future.cancel()
            if stop_flag["stop"]:
                results.append({"status": "FAILED", "failure_stage": "shutdown_timeout",
                                "failure_reason": "cancelled_during_graceful_shutdown"})
            else:
                break  # cycle walltime exceeded; leave remaining futures unresolved this cycle
        except Exception as error:  # worker process crash, pickling error, etc.
            results.append({"status": "FAILED", "failure_stage": "worker_crash", "failure_reason": repr(error)})
    return results


def _new_ledger_candidate_ids_since(ledger_path: Path, start_offset: int) -> set[str]:
    """Candidate ids committed to the ledger after `start_offset`.

    climb/continue_lineage (policy_adapter.py's "self_contained" strategies)
    call `Ledger.append` internally rather than going through the daemon's
    own batch dedup path, so without this the daemon's live `seen_ids` set
    never learns about candidates they just committed and a later cycle can
    re-propose/re-evaluate the same seed as if it were new.
    """
    if not ledger_path.exists():
        return set()
    ids: set[str] = set()
    with ledger_path.open("rb") as handle:
        handle.seek(start_offset)
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                continue  # a torn trailing line from this same in-flight append; picked up next cycle
            if event.get("event") == "EVALUATION" and event.get("candidate_id"):
                ids.add(event["candidate_id"])
    return ids


# --------------------------------------------------------------------- cycle
def _run_cycle(config: dict[str, Any], run_dir: Path, envelope, ledger: Ledger, archive: DiscoveryArchive,
               evaluator_path: Path, pool: cf.ProcessPoolExecutor, state: dict[str, Any],
               seen_ids: set[str], stop_flag: dict[str, Any]) -> dict[str, Any]:
    cycle = state["cycle"]
    cycle_budget = config["runtime"]["max_evaluations_per_cycle"]
    max_total = config["budgets"]["max_total_evaluations"]
    if max_total is not None:
        cycle_budget = max(0, min(cycle_budget, max_total - state["total_evaluations"]))

    allocation = pa.resolve_allocation(config["policy"]["allocation"])
    batch_alloc = [(n, w) for n, w in allocation if n in pa.BATCH_REGISTRY]
    sc_alloc = [(n, w) for n, w in allocation if n in pa.SELF_CONTAINED_REGISTRY]

    context = pa.ProposeContext(run_dir=run_dir, envelope=envelope, worker_id=config["campaign"]["name"],
                                generation=cycle, rng_seed=config["policy"]["seed"] + cycle,
                                best_by_family=_best_by_family(ledger))

    cycle_evaluated = cycle_proposed = cycle_duplicates = cycle_failures = 0
    improvements = 0
    strategy_reports: dict[str, Any] = {}
    batch_results_for_family_scan: list[dict[str, Any]] = []
    grace_s = config["runtime"]["shutdown_grace_s"]
    max_cycle_walltime_s = config["budgets"]["max_cycle_walltime_s"]
    cycle_deadline = time.time() + max_cycle_walltime_s if max_cycle_walltime_s is not None else None

    if batch_alloc and cycle_budget > 0:
        for name, budget in pa.split_budget(batch_alloc, cycle_budget).items():
            if budget <= 0 or stop_flag["stop"]:
                continue
            strategy = pa.BATCH_REGISTRY[name]
            propose_result = strategy.propose(context, budget, state["policy_state"].get(name, {}))
            cycle_proposed += len(propose_result.proposals)
            kept, dropped = _dedupe(propose_result.proposals, seen_ids, envelope)
            cycle_duplicates += dropped

            futures = []
            for proposal in kept:
                attempt_dir = run_dir / "attempts" / f"c{cycle:05d}_{name}_{proposal['worker_id']}_{proposal['proposal_id'][:50]}"
                futures.append(pool.submit(_evaluate_worker, str(ROOT), str(evaluator_path),
                                           proposal, str(attempt_dir), envelope))
            results = _wait_batch(futures, stop_flag, grace_s, cycle_deadline)

            for proposal, result in zip(kept, results):
                ledger.append({"event": "EVALUATION", "lifecycle": result.get("status", "FAILED"), **result})
                cycle_evaluated += 1
                cid = result.get("candidate_id")
                if cid:
                    seen_ids.add(cid)
                if result.get("status") == "TERMINATED":
                    state["consecutive_evaluator_failures"] = 0
                    if result.get("validity_gate"):
                        improvements += 1
                        batch_results_for_family_scan.append(result)
                elif result.get("failure_stage") == "evaluator_execution":
                    cycle_failures += 1
                    state["consecutive_evaluator_failures"] += 1
                else:
                    state["consecutive_evaluator_failures"] = 0

            new_state = strategy.tell(kept, results, propose_result.state)
            state["policy_state"][name] = new_state
            strategy_reports[name] = {"proposed": len(propose_result.proposals), "kept": len(kept),
                                      "duplicates": dropped, "evaluated": len(results)}

    if sc_alloc and cycle_budget > 0 and not stop_flag["stop"]:
        evaluator = DiscoveryEvaluator(ROOT, evaluator_path)
        for name, budget in pa.split_budget(sc_alloc, cycle_budget).items():
            if budget <= 0 or stop_flag["stop"]:
                continue
            strategy = pa.SELF_CONTAINED_REGISTRY[name]
            ledger_offset_before = ledger.path.stat().st_size if ledger.path.exists() else 0
            out = strategy.run(context, budget, evaluator, ledger, state["policy_state"].get(name, {}))
            seen_ids.update(_new_ledger_candidate_ids_since(ledger.path, ledger_offset_before))
            state["policy_state"][name] = out.pop("state", state["policy_state"].get(name, {}))
            cycle_evaluated += out.get("evaluations", 0)
            if out.get("ok") and not out.get("skipped"):
                improvements += 1
            strategy_reports[name] = out

    families_validated_this_cycle = 0
    family_budget = config["policy"]["family_validate_per_cycle"]
    if family_budget > 0 and not stop_flag["stop"]:
        # Any valid TERMINATED candidate is a usable anchor: mirror cli.cmd_validate_family's
        # own convention of forcing mH_GeV=200.0 on whatever coordinates were found, then let
        # the frozen evaluator decide cross-mass/window validity. Requiring an exact mH==200.0
        # match here would (almost) never fire, since the strategies sample mH continuously.
        anchors = list(batch_results_for_family_scan)
        anchors.sort(key=lambda r: 0.0 if r.get("ctau_mm") and 500.0 <= r["ctau_mm"] <= 1000.0 else 1.0)
        evaluator = DiscoveryEvaluator(ROOT, evaluator_path)
        for result in anchors[:family_budget]:
            coords = dict(result["coordinates"])
            coords["mH_GeV"] = 200.0
            record = validate_family(coords, evaluator, run_dir, envelope,
                                     config["campaign"]["name"], "daemon_family_validation", ledger)
            ledger.append({"event": "FAMILY_VALIDATION", "lifecycle": "TERMINATED", **record})
            decision = archive.consider(record, envelope)
            ledger.append({"event": "ARCHIVE_DECISION", "lifecycle": "TERMINATED",
                           "family_id": record["family_id"], **decision})
            families_validated_this_cycle += 1
            state["families_validated"] += 1
            cycle_evaluated += 3
            if decision.get("promoted"):
                improvements += 1

    state["total_evaluations"] += cycle_evaluated
    state["total_proposed"] += cycle_proposed
    state["total_duplicates_dropped"] += cycle_duplicates
    state["consecutive_no_improvement"] = 0 if improvements > 0 else state["consecutive_no_improvement"] + 1

    return {
        "proposed": cycle_proposed, "duplicates_dropped": cycle_duplicates,
        "evaluated": cycle_evaluated, "evaluator_execution_failures": cycle_failures,
        "improvements": improvements, "families_validated": families_validated_this_cycle,
        "strategies": strategy_reports, "total_evaluations": state["total_evaluations"],
    }


def _shutdown_pool(pool: cf.ProcessPoolExecutor, bounded: bool) -> None:
    """Release the worker pool.

    On a normal (non-signal) exit, shut down the ordinary way (wait for
    in-flight work). On a signal-triggered shutdown, `_wait_batch` has
    already spent up to `shutdown_grace_s` waiting for in-flight futures;
    anything still running by the time we get here is treated as hung, not
    merely slow. `ProcessPoolExecutor.shutdown(wait=True)` (the default,
    and what a bare `with` block's __exit__ does) would otherwise block the
    daemon process indefinitely on one stuck evaluator subprocess, defeating
    `shutdown_grace_s` entirely -- so cancel outstanding work and kill the
    worker processes outright instead.

    `pool._processes` is a private implementation detail, and the pool's own
    manager thread can still be in the middle of spawning a replacement
    worker (to keep the pool at `max_workers`) at the exact moment we
    snapshot it -- a worker that comes into existence a moment later would
    be missed by a single pass and would keep holding this process's
    inherited stdout/stderr file descriptors open, which can wedge an
    external supervisor reading this process's output via a pipe even after
    this process itself has exited. Poll a few times over a short window
    instead of trusting one snapshot.
    """
    if not bounded:
        pool.shutdown(wait=True)
        return
    pool.shutdown(wait=False, cancel_futures=True)
    for _ in range(10):
        processes = getattr(pool, "_processes", None) or {}
        alive = [proc for proc in processes.values() if proc.is_alive()]
        if not alive:
            break
        for proc in alive:
            proc.kill()
        time.sleep(0.1)


def _kill_process_group_except_self(pgid: int) -> None:
    """SIGKILL every live process in `pgid` except the caller.

    `pool._processes` is not the whole descendant tree: the pool's manager
    thread can still be spawning a replacement worker as we shut down, and
    `ProcessPoolExecutor` also lazily starts an independent
    `multiprocessing.resource_tracker` process that is never tracked in
    `_processes` at all. Either can outlive `_shutdown_pool` above and keep
    holding this process's inherited stdout/stderr pipes open, wedging an
    external supervisor's `Popen(...).communicate()` even after this process
    has exited. Sweeping by process-group membership (set up via
    `os.setpgrp()` earlier in `run()`) catches all of it without depending
    on any executor-private bookkeeping. Best-effort: /proc access or a
    lost race on a given pid are not fatal.
    """
    self_pid = os.getpid()
    try:
        candidates = os.listdir("/proc")
    except OSError:
        return
    for name in candidates:
        if not name.isdigit():
            continue
        pid = int(name)
        if pid == self_pid:
            continue
        try:
            if os.getpgid(pid) != pgid:
                continue
            os.kill(pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError, OSError):
            continue


# --------------------------------------------------------------------- loop
def run(config: dict[str, Any], run_dir: Path, resume: bool, workers_override: int | None,
        dry_run: bool, stop_after_cycles: int | None = None) -> int:
    guard_run_dir(run_dir)
    run_dir = Path(run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    envelope = ENVELOPES[config["envelope_id"]]
    save_envelope(run_dir, envelope)
    ledger = Ledger(run_dir)
    archive = DiscoveryArchive(run_dir)
    guard_write_path(ledger.path, run_dir)
    guard_write_path(archive.path, run_dir)
    if repair_torn_ledger_tail(ledger.path):
        print(f"RECOVERED: truncated a torn trailing ledger line in {ledger.path} "
              f"(hard kill mid-append; no complete record was discarded)", file=sys.stderr)

    digest = _config_digest(config)
    state = checkpoint_mod.load(run_dir) if resume else None
    if state is None:
        state = checkpoint_mod.fresh(digest)
    elif state.get("config_digest") != digest:
        print(f"WARNING: resuming {run_dir} with a changed policy/envelope config "
              f"(checkpoint digest {state.get('config_digest')} != current {digest})", file=sys.stderr)

    seen_ids = _replay_seen_ids(ledger)

    if dry_run:
        return _run_dry(config, run_dir, envelope, state, seen_ids, ledger)

    workers = workers_override if workers_override is not None else config["runtime"]["workers"]
    evaluator_path = Path(config["evaluator"])
    if not evaluator_path.is_absolute():
        evaluator_path = ROOT / evaluator_path

    stop_flag = {"stop": False, "signal": None}

    def _handle_signal(signum, _frame):
        stop_flag["stop"] = True
        stop_flag["signal"] = signum

    previous_handlers = {
        signal.SIGINT: signal.signal(signal.SIGINT, _handle_signal),
        signal.SIGTERM: signal.signal(signal.SIGTERM, _handle_signal),
    }

    start_time = time.time()
    exit_code = 0
    # New, isolated process group for this run: every worker and every
    # multiprocessing-internal helper process (notably resource_tracker,
    # which ProcessPoolExecutor starts lazily and independently of
    # pool._processes) inherits it. This is what makes the bounded-shutdown
    # sweep below able to reliably clean up the whole descendant tree by
    # group membership rather than by walking a private executor attribute.
    # Best-effort: if the sandbox disallows changing the process group (or
    # we're already a group leader), the per-process kill loop in
    # _shutdown_pool is still a reasonable fallback.
    try:
        os.setpgrp()
    except OSError:
        pass
    # "spawn", not the platform-default "fork": a forked worker inherits
    # copies of this process's own stdout/stderr file descriptors, which can
    # keep an external reader (e.g. a test or supervisor doing
    # Popen(..., stdout=PIPE).communicate()) blocked past EOF even after this
    # process has genuinely exited, if a worker outlives a forced shutdown by
    # even a moment. spawn() starts each worker as a fresh interpreter that
    # does not blanket-inherit unrelated fds.
    pool = cf.ProcessPoolExecutor(max_workers=max(1, workers),
                                  mp_context=multiprocessing.get_context("spawn"))
    try:
        while True:
            stop_reason = _stop_reason(config, state, start_time, stop_after_cycles, stop_flag)
            if stop_reason is not None:
                state["stop_reason"] = stop_reason
                if stop_flag["stop"]:
                    exit_code = SIGNAL_EXIT_CODES.get(stop_flag["signal"], 1)
                break
            cycle_report = _run_cycle(config, run_dir, envelope, ledger, archive, evaluator_path,
                                      pool, state, seen_ids, stop_flag)
            state["cycle"] += 1
            if state["cycle"] % max(1, config["runtime"]["checkpoint_every"]) == 0:
                state = checkpoint_mod.save(run_dir, state)
            if state["cycle"] % max(1, config["runtime"]["summary_every"]) == 0:
                print(json.dumps({"cycle": state["cycle"], **cycle_report}, sort_keys=True, default=str))
    finally:
        signal.signal(signal.SIGINT, previous_handlers[signal.SIGINT])
        signal.signal(signal.SIGTERM, previous_handlers[signal.SIGTERM])
        _shutdown_pool(pool, bounded=stop_flag["stop"])

    state["stopped"] = True
    state = checkpoint_mod.save(run_dir, state)
    print(json.dumps({
        "shutdown_summary": True, "cycle": state["cycle"],
        "total_evaluations": state["total_evaluations"], "stop_reason": state.get("stop_reason"),
        "exit_code": exit_code, "wall_time_s": round(time.time() - start_time, 3),
    }, sort_keys=True))
    if stop_flag["stop"]:
        # All essential state (ledger row, checkpoint) is already durably
        # written above. Two things can otherwise keep this process's exit
        # from being externally observable promptly: (1) a worker or the
        # independent resource_tracker helper outliving _shutdown_pool and
        # still holding this process's inherited stdout/stderr pipes open,
        # wedging a supervisor's Popen(...).communicate() past EOF even
        # though this process itself has exited -- swept by process-group
        # membership just below; (2) concurrent.futures.process's own atexit
        # hook joining every ProcessPoolExecutor manager thread ever created,
        # which can hang on a pool torn down by force rather than its own
        # shutdown protocol -- avoided by skipping Python's normal
        # interpreter shutdown entirely via os._exit.
        sys.stdout.flush()
        sys.stderr.flush()
        try:
            _kill_process_group_except_self(os.getpgrp())
        except OSError:
            pass
        os._exit(exit_code)
    return exit_code


def _stop_reason(config: dict[str, Any], state: dict[str, Any], start_time: float,
                 stop_after_cycles: int | None, stop_flag: dict[str, Any]) -> str | None:
    if stop_flag["stop"]:
        return f"signal:{stop_flag['signal']}"
    max_total = config["budgets"]["max_total_evaluations"]
    if max_total is not None and state["total_evaluations"] >= max_total:
        return "max_total_evaluations_reached"
    max_wall = config["budgets"]["max_total_walltime_s"]
    if max_wall is not None and (time.time() - start_time) >= max_wall:
        return "max_total_walltime_reached"
    if state["consecutive_no_improvement"] >= config["stopping"]["patience_cycles"]:
        return "patience_exhausted"
    if state["consecutive_evaluator_failures"] >= config["stopping"]["max_consecutive_evaluator_failures"]:
        return "too_many_consecutive_evaluator_failures"
    if stop_after_cycles is not None and state["cycle"] >= stop_after_cycles:
        return "stop_after_cycles_reached"
    return None


def _run_dry(config: dict[str, Any], run_dir: Path, envelope, state: dict[str, Any],
            seen_ids: set[str], ledger: Ledger) -> int:
    evaluator_path = Path(config["evaluator"])
    if not evaluator_path.is_absolute():
        evaluator_path = ROOT / evaluator_path
    allocation = pa.resolve_allocation(config["policy"]["allocation"])
    shares = pa.split_budget(allocation, config["runtime"]["max_evaluations_per_cycle"])
    context = pa.ProposeContext(run_dir=run_dir, envelope=envelope, worker_id=config["campaign"]["name"],
                                generation=state["cycle"], rng_seed=config["policy"]["seed"],
                                best_by_family=_best_by_family(ledger))
    plan: dict[str, Any] = {}
    for name, budget in shares.items():
        if name in pa.BATCH_REGISTRY:
            strategy = pa.BATCH_REGISTRY[name]
            result = strategy.propose(context, budget, state["policy_state"].get(name, {}))
            kept, dropped = _dedupe(result.proposals, seen_ids, envelope)
            plan[name] = {"kind": "batch", "budget": budget, "proposed": len(result.proposals),
                          "duplicates_would_drop": dropped, "would_evaluate": len(kept)}
        else:
            plan[name] = {"kind": "self_contained", "budget": budget,
                          "would_run": budget > 0 and _best_by_family(ledger).get("mixed") is not None}
    report = {
        "dry_run": True, "evaluator_subprocess_invocations": 0,
        "run_dir": str(run_dir), "envelope_id": envelope.envelope_id,
        "evaluator_binary_exists": evaluator_path.exists(),
        "resume_checkpoint_found": Path(checkpoint_mod.path_for(run_dir)).exists(),
        "seen_candidate_ids_from_ledger": len(seen_ids),
        "resolved_allocation": {name: round(weight, 4) for name, weight in allocation},
        "cycle_plan": plan, "runtime": config["runtime"], "budgets": config["budgets"],
        "stopping": config["stopping"], "safety": config["safety"],
    }
    print(json.dumps(report, indent=2, sort_keys=True, default=str))
    return 0
