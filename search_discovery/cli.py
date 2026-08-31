"""Worker-facing discovery CLI.

Every command that proposes points also deduplicates, evaluates through the
frozen binary, appends to the append-only ledger, and offers results to the
archive. Steps 3-6 of the campaign loop therefore cannot be skipped by an agent.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np

from search_substrate.ledger import Ledger

from .archive import DiscoveryArchive
from .bounds import GLOBAL_PHYSICS_BOUNDS, Envelope
from .envelopes import ENVELOPES, active_envelope, save_envelope
from .evaluator import DiscoveryEvaluator
from .family import validate_family
from .helpers import CMAES, COORD_ORDER, deduplicate, perturb, random_candidates, sobol_candidates, to_unit_vector
from .objective import objective, max_photon_objective

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVALUATOR = ROOT / "dihiggs/app/DihiggsPointV2Evaluator"


def _emit(document: Any) -> None:
    print(json.dumps(document, indent=2, sort_keys=True, allow_nan=False, default=str))


def _seen_ids(ledger: Ledger) -> set[str]:
    out = set()
    for event in ledger.events():
        cid = event.get("candidate_id")
        if cid and event.get("event") == "EVALUATION":
            out.add(cid)
    return out


def _compact(result: dict[str, Any], target: str | None = None) -> dict[str, Any]:
    if result.get("status") != "TERMINATED":
        return {"candidate_id": result.get("candidate_id"), "status": result.get("status"),
                "failure_stage": result.get("failure_stage"), "failure_reason": result.get("failure_reason")}
    obs = result.get("observables") or {}
    prov = result.get("provisional_family") or {}
    row = {
        "candidate_id": result.get("candidate_id"),
        "valid": result.get("validity_gate"),
        "failed_gates": sorted(n for n, ok in (result.get("gates") or {}).items() if not ok),
        "ctau_mm": result.get("ctau_mm"),
        "provisional_family": prov.get("family"),
        "photon_fraction": prov.get("photon_fraction"),
        "fermionic_fraction": prov.get("fermionic_fraction"),
        "br_bb": obs.get("br_bb"), "br_gg": obs.get("br_gg"),
        "br_gammagamma": obs.get("br_gammagamma"), "br_Zgamma": obs.get("br_Zgamma"),
        "total_width_GeV": obs.get("total_width_GeV"),
        "lambda1": obs.get("lambda1_reconstructed"),
        "coordinates": result.get("coordinates"),
    }
    if target == "__maxphoton__":
        row["objective"] = max_photon_objective(result)
    elif target:
        row["objective"] = objective(target, result)
    return row


def _run_batch(proposals, evaluator, ledger, run_dir, envelope, archive_target=None):
    seen = _seen_ids(ledger)
    kept, dropped = deduplicate(proposals, seen, envelope)
    results = []
    for index, proposal in enumerate(kept):
        attempt = run_dir / "attempts" / f"{proposal['worker_id']}_{proposal['proposal_id'][:60]}"
        result = evaluator.evaluate(proposal, attempt, envelope)
        ledger.append({"event": "EVALUATION", "lifecycle": result.get("status", "FAILED"), **result})
        results.append(result)
    return results, dropped


def cmd_envelope(args):
    if args.list:
        _emit({"global_physics_bounds": GLOBAL_PHYSICS_BOUNDS,
               "envelopes": {k: v.as_dict() for k, v in ENVELOPES.items()}})
        return
    envelope = active_envelope(Path(args.run_dir)) if args.run_dir else ENVELOPES["E0"]
    _emit({"global_physics_bounds": GLOBAL_PHYSICS_BOUNDS, "active_envelope": envelope.as_dict()})


def cmd_set_envelope(args):
    envelope = ENVELOPES[args.envelope_id]
    save_envelope(Path(args.run_dir), envelope)
    _emit({"active_envelope": envelope.as_dict()})


def cmd_explore(args):
    run_dir = Path(args.run_dir).resolve()
    envelope = active_envelope(run_dir)
    ledger, archive = Ledger(run_dir), DiscoveryArchive(run_dir)
    evaluator = DiscoveryEvaluator(ROOT, Path(args.evaluator))
    maker = sobol_candidates if args.mode == "sobol" else random_candidates
    proposals = maker(args.count, args.seed, envelope, args.worker_id, args.generation)
    results, dropped = _run_batch(proposals, evaluator, ledger, run_dir, envelope)
    valid = [r for r in results if r.get("validity_gate")]
    _emit({"worker_id": args.worker_id, "mode": args.mode, "seed": args.seed,
           "envelope_id": envelope.envelope_id,
           "proposed": len(proposals), "duplicates_dropped": dropped, "evaluated": len(results),
           "valid": len(valid), "valid_fraction": round(len(valid) / max(1, len(results)), 4),
           "results": [_compact(r) for r in results]})


def cmd_refine(args):
    run_dir = Path(args.run_dir).resolve()
    envelope = active_envelope(run_dir)
    ledger, archive = Ledger(run_dir), DiscoveryArchive(run_dir)
    evaluator = DiscoveryEvaluator(ROOT, Path(args.evaluator))
    parent = None
    for event in ledger.events():
        if event.get("candidate_id") == args.parent and event.get("coordinates"):
            parent = event
    if parent is None:
        _emit({"error": "parent_candidate_not_found_in_ledger", "parent": args.parent})
        sys.exit(2)
    proposals = perturb(parent["coordinates"], args.count, args.seed, args.radius,
                        envelope, args.worker_id, args.parent, args.generation)
    results, dropped = _run_batch(proposals, evaluator, ledger, run_dir, envelope)
    target = args.target
    parent_obj = objective(target, parent) if target else None
    rows = [_compact(r, target) for r in results]
    improved = [r for r in rows if r.get("objective", {}).get("value", 1e9)
                < (parent_obj or {}).get("value", 1e9)]
    _emit({"worker_id": args.worker_id, "parent": args.parent, "radius": args.radius,
           "parent_objective": parent_obj, "duplicates_dropped": dropped,
           "evaluated": len(results), "improved_count": len(improved),
           "results": rows})


def cmd_optimize(args):
    run_dir = Path(args.run_dir).resolve()
    envelope = active_envelope(run_dir)
    ledger = Ledger(run_dir)
    evaluator = DiscoveryEvaluator(ROOT, Path(args.evaluator))
    rng = np.random.default_rng(args.seed)
    x0 = rng.random(len(COORD_ORDER)) if args.x0 is None else np.array(json.loads(args.x0), dtype=float)
    if args.start_candidate:
        for event in ledger.events():
            if event.get("candidate_id") == args.start_candidate and event.get("coordinates"):
                x0 = to_unit_vector(event["coordinates"], envelope)
    es = CMAES(x0, args.sigma0, args.seed, args.popsize)
    history, best = [], None
    for _ in range(args.generations):
        proposals = es.proposals(envelope, args.worker_id, args.seed)
        results, dropped = _run_batch(proposals, evaluator, ledger, run_dir, envelope)
        fitnesses, vectors = [], []
        for proposal, result in zip(proposals, results):
            score = (max_photon_objective(result, args.require_window)
                     if args.objective == "max_photon_valid" else objective(args.target, result))
            fitnesses.append(score["value"])
            vectors.append(to_unit_vector(
                {k: proposal["coordinates"][k] for k in COORD_ORDER}, envelope))
            if best is None or score["value"] < best["objective"]["value"]:
                best = {"candidate_id": result.get("candidate_id"), "objective": score,
                        "compact": _compact(result, "__maxphoton__"
                                            if args.objective == "max_photon_valid" else args.target)}
        es.tell(vectors, fitnesses)
        history.append({"generation": es.generation, "sigma": round(es.sigma, 6),
                        "best_in_gen": round(float(min(fitnesses)), 6),
                        "median": round(float(np.median(fitnesses)), 6),
                        "duplicates_dropped": dropped})
    stagnated = (len(history) >= 3
                 and abs(history[-1]["best_in_gen"] - history[-3]["best_in_gen"]) < 1e-6)
    _emit({"worker_id": args.worker_id, "target_family": args.target, "seed": args.seed,
           "generations": args.generations, "sigma_final": es.sigma,
           "stagnated": stagnated, "history": history, "best": best})


def cmd_validate_family(args):
    run_dir = Path(args.run_dir).resolve()
    envelope = active_envelope(run_dir)
    ledger, archive = Ledger(run_dir), DiscoveryArchive(run_dir)
    evaluator = DiscoveryEvaluator(ROOT, Path(args.evaluator))
    anchor = None
    for event in ledger.events():
        if event.get("candidate_id") == args.anchor_candidate and event.get("coordinates"):
            anchor = event
    if anchor is None:
        _emit({"error": "anchor_candidate_not_found", "anchor": args.anchor_candidate})
        sys.exit(2)
    coords = dict(anchor["coordinates"])
    if abs(coords["mH_GeV"] - 200.0) > 1e-9:
        coords["mH_GeV"] = 200.0  # continuation is defined from the 200 GeV anchor
    record = validate_family(coords, evaluator, run_dir, envelope,
                             args.worker_id, args.strategy_id, ledger)
    ledger.append({"event": "FAMILY_VALIDATION", "lifecycle": "TERMINATED", **record})
    decision = archive.consider(record, envelope)
    ledger.append({"event": "ARCHIVE_DECISION", "lifecycle": "TERMINATED",
                   "family_id": record["family_id"], **decision})
    record["archive_decision"] = decision
    record.pop("members", None)
    _emit(record)


def cmd_candidates(args):
    """List evaluated candidates, best-first, for lineage selection."""
    run_dir = Path(args.run_dir).resolve()
    rows = []
    for event in Ledger(run_dir).events():
        if event.get("event") != "EVALUATION" or event.get("status") != "TERMINATED":
            continue
        if args.valid_only and not event.get("validity_gate"):
            continue
        prov = event.get("provisional_family") or {}
        if args.family and prov.get("family") != args.family:
            continue
        score = objective(args.target, event) if args.target else None
        rows.append({"candidate_id": event.get("candidate_id"),
                     "strategy_id": event.get("strategy_id"),
                     "worker_id": event.get("worker_id"),
                     "valid": event.get("validity_gate"),
                     "ctau_mm": event.get("ctau_mm"),
                     "provisional_family": prov.get("family"),
                     "photon_fraction": prov.get("photon_fraction"),
                     "objective": (score or {}).get("value"),
                     "coordinates": event.get("coordinates")})
    key = (lambda r: (r["objective"] if r["objective"] is not None else 1e9)) if args.target \
        else (lambda r: -(r["ctau_mm"] or 0.0))
    rows.sort(key=key)
    _emit({"count": len(rows), "candidates": rows[: args.limit]})


def _eval_coords(coords, evaluator, ledger, run_dir, worker_id, strategy_id, tag):
    import hashlib
    proposal = {"proposal_id": f"{tag}_" + hashlib.sha256(
                    repr(sorted(coords.items())).encode()).hexdigest()[:12],
                "strategy_id": strategy_id, "worker_id": worker_id, "generation": 0,
                "parent_ids": [], "random_seed": None,
                "rationale": f"explicit coordinate evaluation ({tag})", "coordinates": coords}
    result = evaluator.evaluate(proposal, run_dir / "attempts" / proposal["proposal_id"], None)
    ledger.append({"event": "EVALUATION", "lifecycle": result.get("status", "FAILED"), **result})
    return result


def cmd_evaluate_point(args):
    """Evaluate one explicit coordinate. Global physics bounds still apply."""
    run_dir = Path(args.run_dir).resolve()
    ledger = Ledger(run_dir)
    evaluator = DiscoveryEvaluator(ROOT, Path(args.evaluator))
    coords = {"mH_GeV": args.mH, "mA_GeV": args.mA, "tan_beta": args.tan_beta,
              "X": args.X, "Q": args.Q}
    result = _eval_coords(coords, evaluator, ledger, run_dir, args.worker_id,
                          "explicit_point", "point")
    _emit(_compact(result, args.target))


def cmd_reanchor(args):
    """Re-anchor a candidate to mH=200 and solve tan_beta for ctau_200 in target.

    Uses the measured scaling ctau ~ tan_beta^2 at fixed Q as the update step,
    but every iterate is CONFIRMED by the frozen evaluator; nothing is assumed.
    """
    run_dir = Path(args.run_dir).resolve()
    ledger = Ledger(run_dir)
    evaluator = DiscoveryEvaluator(ROOT, Path(args.evaluator))
    source = None
    for event in ledger.events():
        if event.get("candidate_id") == args.candidate and event.get("coordinates"):
            source = event
    if source is None:
        _emit({"error": "candidate_not_found", "candidate": args.candidate}); sys.exit(2)
    coords = dict(source["coordinates"])
    s = coords["mA_GeV"] ** 2 - coords["mH_GeV"] ** 2      # preserve S
    coords["mH_GeV"] = 200.0
    coords["mA_GeV"] = math.sqrt(200.0 ** 2 + s)
    trace = []
    for iteration in range(args.max_iterations):
        result = _eval_coords(coords, evaluator, ledger, run_dir, args.worker_id,
                              "reanchor", f"reanchor{iteration}")
        ctau = result.get("ctau_mm")
        valid = bool(result.get("validity_gate"))
        trace.append({"iteration": iteration, "tan_beta": coords["tan_beta"],
                      "ctau_mm": ctau, "valid": valid,
                      "candidate_id": result.get("candidate_id"),
                      "failed_gates": sorted(n for n, ok in (result.get("gates") or {}).items() if not ok)})
        if ctau and 500.0 <= ctau <= 1000.0 and valid:
            _emit({"candidate": args.candidate, "converged": True, "iterations": iteration + 1,
                   "final_coordinates": coords, "trace": trace,
                   "final": _compact(result, args.target)})
            return
        if not ctau or ctau <= 0:
            break
        scale = math.sqrt(700.0 / ctau)
        new_tb = coords["tan_beta"] * scale
        lo, hi = GLOBAL_PHYSICS_BOUNDS["tan_beta"]
        coords["tan_beta"] = min(max(new_tb, lo), hi)
        if abs(coords["tan_beta"] - trace[-1]["tan_beta"]) / trace[-1]["tan_beta"] < 1e-9:
            break
    _emit({"candidate": args.candidate, "converged": False,
           "iterations": len(trace), "final_coordinates": coords, "trace": trace})


def cmd_reselect_archive(args):
    """Deterministic diversity-aware Top-K re-selection over ALL recorded families.

    Rationale: every cross-mass-valid, in-window family scores exactly 1.0, so the
    plain score ordering leaves the Top-K decided by a sha256 tie-break rather
    than by physics. The campaign brief requires the surviving five to span
    distinct basins. This command re-derives the archive from the ledger's
    complete, already-evaluated family history using max-min dispersion:

        seed  = highest score, lexicographically smallest family_id on ties
        step  = add the family maximising the MINIMUM normalised distance to the
                already-selected set, among those of the highest remaining score

    It runs no physics and creates no candidates. It is deterministic and
    reproducible from the ledger alone.
    """
    run_dir = Path(args.run_dir).resolve()
    envelope = active_envelope(run_dir)
    ledger, archive = Ledger(run_dir), DiscoveryArchive(run_dir)
    from .archive import FAMILIES, TOP_K, family_distance
    pool: dict[str, dict[str, Any]] = {}
    for event in ledger.events():
        if event.get("event") != "FAMILY_VALIDATION":
            continue
        if not (event.get("cross_mass_valid") and event.get("same_X")):
            continue
        if event.get("family") not in FAMILIES:
            continue
        pool[event["family_id"]] = event
    selected: dict[str, list[dict[str, Any]]] = {}
    report: dict[str, Any] = {}
    for family in FAMILIES:
        pending = [record for record in pool.values() if record["family"] == family]
        pending.sort(key=lambda r: (-float(r.get("score", 0.0)), str(r["family_id"])))
        chosen: list[dict[str, Any]] = []
        while pending and len(chosen) < TOP_K:
            if not chosen:
                chosen.append(pending.pop(0))
                continue
            best_score = max(float(r.get("score", 0.0)) for r in pending)
            tier = [r for r in pending if float(r.get("score", 0.0)) == best_score]
            def min_distance(record):
                return min(family_distance(record["anchor"], c["anchor"], envelope) for c in chosen)
            tier.sort(key=lambda r: (-min_distance(r), str(r["family_id"])))
            pick = tier[0]
            pending.remove(pick)
            chosen.append(pick)
        selected[family] = chosen
        distances = [family_distance(a["anchor"], b["anchor"], envelope)
                     for i, a in enumerate(chosen) for b in chosen[i + 1:]]
        report[family] = {
            "pool_size": len([r for r in pool.values() if r["family"] == family]),
            "selected": len(chosen),
            "min_pairwise_distance": round(min(distances), 6) if distances else None,
            "max_pairwise_distance": round(max(distances), 6) if distances else None,
            "entries": [{"family_id": r["family_id"], "score": r["score"],
                         "ctau_200_mm": r["ctau_200_mm"],
                         "ctau_by_mass_mm": r.get("ctau_by_mass_mm"),
                         "in_target_interval": r.get("in_target_interval"),
                         "coordinates": r["anchor"]["coordinates"]} for r in chosen],
        }
    stripped = {f: [{k: v for k, v in r.items() if k != "members"} for r in rows]
                for f, rows in selected.items()}
    temporary = archive.path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(stripped, indent=2, sort_keys=True, allow_nan=False) + "\n",
                         encoding="utf-8")
    os.replace(temporary, archive.path)
    ledger.append({"event": "ARCHIVE_RESELECTED", "lifecycle": "TERMINATED",
                   "method": "max_min_dispersion_v1",
                   "selected_family_ids": {f: [r["family_id"] for r in rows]
                                           for f, rows in selected.items()}})
    _emit({"method": "max_min_dispersion_v1", "envelope_id": envelope.envelope_id, **report})


def cmd_climb(args):
    """Constrained ascent on photon fraction from a validated mixed benchmark."""
    run_dir = Path(args.run_dir).resolve()
    envelope = active_envelope(run_dir)
    ledger = Ledger(run_dir)
    evaluator = DiscoveryEvaluator(ROOT, Path(args.evaluator))
    start = None
    for event in ledger.events():
        if event.get("candidate_id") == args.parent and event.get("coordinates"):
            start = event
    if start is None:
        _emit({"error": "parent_not_found", "parent": args.parent}); sys.exit(2)
    from .climb import ascend
    out = ascend(dict(start["coordinates"]), evaluator, ledger, run_dir, envelope,
                 args.worker_id, args.max_steps, args.step0)
    if out.get("ok"):
        out["trace"] = [t for t in out["trace"]]
    _emit(out)


def cmd_continue_lineage(args):
    """Constrained ascent on pf200 from a validated mixed family, full-family checked."""
    run_dir = Path(args.run_dir).resolve()
    envelope = active_envelope(run_dir)
    ledger = Ledger(run_dir)
    evaluator = DiscoveryEvaluator(ROOT, Path(args.evaluator))
    seed_coords = None
    if args.seed_family:
        for event in ledger.events():
            if event.get("event") == "FAMILY_VALIDATION" and event.get("family_id") == args.seed_family:
                seed_coords = dict(event["anchor"]["coordinates"])
        if seed_coords is None:
            for event in ledger.events():
                if event.get("event") == "ARCHIVE_RESELECTED":
                    pass  # reselection doesn't carry coordinates; family_id lookup above is authoritative
    elif args.seed_candidate:
        for event in ledger.events():
            if event.get("candidate_id") == args.seed_candidate and event.get("coordinates"):
                seed_coords = dict(event["coordinates"])
    if seed_coords is None:
        _emit({"error": "seed_not_found", "seed_family": args.seed_family,
               "seed_candidate": args.seed_candidate})
        sys.exit(2)
    from .continuation import ascend_lineage
    out = ascend_lineage(seed_coords, evaluator, ledger, run_dir, envelope, args.worker_id,
                         args.lineage_id, args.max_generations, args.step0, args.min_step,
                         args.patience, not args.no_branch)
    _emit(out)


def cmd_checkpoint(args):
    """Compact campaign state for a master allocation decision. No narrative."""
    run_dir = Path(args.run_dir).resolve()
    envelope = active_envelope(run_dir)
    ledger, archive = Ledger(run_dir), DiscoveryArchive(run_dir)
    from .objective import lifetime_distance
    by_strategy: dict[str, dict[str, Any]] = {}
    seen_ids: set[str] = set()
    duplicates = 0
    valid_vectors: list[tuple[np.ndarray, str]] = []
    best = {"mixed": None, "photonic": None}
    counts = {"attempts": 0, "valid": 0, "mixed": 0, "photonic": 0, "in_target": 0}
    gate_failures: dict[str, int] = {}
    for event in ledger.events():
        if event.get("event") != "EVALUATION":
            continue
        counts["attempts"] += 1
        cid = event.get("candidate_id")
        if cid in seen_ids:
            duplicates += 1
        elif cid:
            seen_ids.add(cid)
        sid = str(event.get("strategy_id", "?"))
        slot = by_strategy.setdefault(sid, {"attempts": 0, "valid": 0, "mixed": 0,
                                            "photonic": 0, "in_target": 0})
        slot["attempts"] += 1
        for name, ok in (event.get("gates") or {}).items():
            if not ok:
                gate_failures[name] = gate_failures.get(name, 0) + 1
        if not event.get("validity_gate"):
            continue
        counts["valid"] += 1; slot["valid"] += 1
        fam = (event.get("provisional_family") or {}).get("family")
        if fam in ("mixed", "photonic"):
            counts[fam] += 1; slot[fam] += 1
            distance = lifetime_distance(event.get("ctau_mm"))
            if best[fam] is None or distance < best[fam]["lifetime_distance"]:
                best[fam] = {"candidate_id": cid, "ctau_mm": event.get("ctau_mm"),
                             "lifetime_distance": distance,
                             "coordinates": event.get("coordinates")}
        ctau = event.get("ctau_mm")
        if ctau and 500.0 <= ctau <= 1000.0:
            counts["in_target"] += 1; slot["in_target"] += 1
        if event.get("coordinates"):
            try:
                valid_vectors.append((to_unit_vector(event["coordinates"], envelope), str(fam)))
            except Exception:
                pass
    # distinct basins: greedy clustering of valid points in unit space
    basins: list[dict[str, Any]] = []
    for vector, fam in valid_vectors:
        for basin in basins:
            if float(np.sqrt(((vector - basin["centre"]) ** 2).mean())) <= args.basin_radius:
                basin["n"] += 1
                basin["centre"] = basin["centre"] + (vector - basin["centre"]) / basin["n"]
                break
        else:
            basins.append({"centre": vector, "n": 1, "family": fam})
    current = archive.load()
    _emit({
        "envelope_id": envelope.envelope_id,
        "attempts_total": counts["attempts"],
        "valid_total": counts["valid"],
        "valid_fraction": round(counts["valid"] / max(1, counts["attempts"]), 4),
        "mixed_valid_found": counts["mixed"],
        "photonic_valid_found": counts["photonic"],
        "valid_in_target_ctau": counts["in_target"],
        "duplicate_rate": round(duplicates / max(1, counts["attempts"]), 4),
        "gate_failure_counts": dict(sorted(gate_failures.items(), key=lambda kv: -kv[1])),
        "best_by_family": best,
        "distinct_basins": len(basins),
        "basin_sizes": sorted((b["n"] for b in basins), reverse=True)[:10],
        "by_strategy": by_strategy,
        "archive_counts": {f: len(v) for f, v in current.items()},
    })


def cmd_summary(args):
    run_dir = Path(args.run_dir).resolve()
    ledger, archive = Ledger(run_dir), DiscoveryArchive(run_dir)
    by_strategy: dict[str, dict[str, Any]] = {}
    total = valid = 0
    families: dict[str, int] = {}
    best_ctau_distance = None
    in_target = 0
    for event in ledger.events():
        if event.get("event") != "EVALUATION":
            continue
        total += 1
        sid = str(event.get("strategy_id", "?"))
        slot = by_strategy.setdefault(sid, {"attempts": 0, "valid": 0, "families": {}})
        slot["attempts"] += 1
        if event.get("validity_gate"):
            valid += 1
            slot["valid"] += 1
            fam = (event.get("provisional_family") or {}).get("family")
            if fam:
                families[fam] = families.get(fam, 0) + 1
                slot["families"][fam] = slot["families"].get(fam, 0) + 1
            ctau = event.get("ctau_mm")
            if ctau and 500.0 <= ctau <= 1000.0:
                in_target += 1
            from .objective import lifetime_distance
            distance = lifetime_distance(ctau)
            if best_ctau_distance is None or distance < best_ctau_distance:
                best_ctau_distance = distance
    current = archive.load()
    _emit({"run_dir": str(run_dir), "total_evaluations": total, "valid": valid,
           "valid_fraction": round(valid / max(1, total), 4),
           "valid_by_provisional_family": families,
           "valid_in_target_ctau_interval": in_target,
           "best_anchor_lifetime_distance": best_ctau_distance,
           "by_strategy": by_strategy,
           "archive_counts": {f: len(v) for f, v in current.items()},
           "archive": {f: [{"family_id": e["family_id"], "score": e["score"],
                            "ctau_200_mm": e["ctau_200_mm"], "X": e["X"], "Q": e["Q"],
                            "tan_beta": e["tan_beta"],
                            "coordinates": e["anchor"]["coordinates"]} for e in v]
                       for f, v in current.items()}})


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="search_discovery", description="LLP discovery search layer")
    sub = parser.add_subparsers(dest="command", required=True)
    ev = lambda p: p.add_argument("--evaluator", default=str(DEFAULT_EVALUATOR))

    p = sub.add_parser("envelope"); p.add_argument("--run-dir"); p.add_argument("--list", action="store_true")
    p.set_defaults(func=cmd_envelope)
    p = sub.add_parser("set-envelope"); p.add_argument("--run-dir", required=True)
    p.add_argument("--envelope-id", required=True); p.set_defaults(func=cmd_set_envelope)

    p = sub.add_parser("explore"); p.add_argument("--run-dir", required=True)
    p.add_argument("--worker-id", required=True); p.add_argument("--seed", type=int, required=True)
    p.add_argument("--count", type=int, default=32); p.add_argument("--generation", type=int, default=0)
    p.add_argument("--mode", choices=("sobol", "random"), default="sobol"); ev(p)
    p.set_defaults(func=cmd_explore)

    p = sub.add_parser("refine"); p.add_argument("--run-dir", required=True)
    p.add_argument("--worker-id", required=True); p.add_argument("--parent", required=True)
    p.add_argument("--seed", type=int, required=True); p.add_argument("--count", type=int, default=16)
    p.add_argument("--radius", type=float, default=0.05); p.add_argument("--generation", type=int, default=1)
    p.add_argument("--target", choices=("mixed", "photonic")); ev(p)
    p.set_defaults(func=cmd_refine)

    p = sub.add_parser("optimize"); p.add_argument("--run-dir", required=True)
    p.add_argument("--worker-id", required=True); p.add_argument("--seed", type=int, required=True)
    p.add_argument("--target", choices=("mixed", "photonic"), required=True)
    p.add_argument("--objective", choices=("target", "max_photon_valid"), default="target")
    p.add_argument("--require-window", action="store_true")
    p.add_argument("--generations", type=int, default=10); p.add_argument("--sigma0", type=float, default=0.25)
    p.add_argument("--popsize", type=int); p.add_argument("--x0"); p.add_argument("--start-candidate"); ev(p)
    p.set_defaults(func=cmd_optimize)

    p = sub.add_parser("validate-family"); p.add_argument("--run-dir", required=True)
    p.add_argument("--anchor-candidate", required=True); p.add_argument("--worker-id", default="master")
    p.add_argument("--strategy-id", default="cross_mass_validation"); ev(p)
    p.set_defaults(func=cmd_validate_family)

    p = sub.add_parser("evaluate-point"); p.add_argument("--run-dir", required=True)
    p.add_argument("--worker-id", default="master"); p.add_argument("--mH", type=float, required=True)
    p.add_argument("--mA", type=float, required=True); p.add_argument("--tan-beta", type=float, required=True)
    p.add_argument("--X", type=float, required=True); p.add_argument("--Q", type=float, required=True)
    p.add_argument("--target", choices=("mixed", "photonic")); ev(p)
    p.set_defaults(func=cmd_evaluate_point)

    p = sub.add_parser("reanchor"); p.add_argument("--run-dir", required=True)
    p.add_argument("--candidate", required=True); p.add_argument("--worker-id", default="master")
    p.add_argument("--max-iterations", type=int, default=8)
    p.add_argument("--target", choices=("mixed", "photonic"), default="mixed"); ev(p)
    p.set_defaults(func=cmd_reanchor)

    p = sub.add_parser("reselect-archive"); p.add_argument("--run-dir", required=True)
    p.set_defaults(func=cmd_reselect_archive)

    p = sub.add_parser("climb"); p.add_argument("--run-dir", required=True)
    p.add_argument("--parent", required=True); p.add_argument("--worker-id", default="climb")
    p.add_argument("--max-steps", type=int, default=40); p.add_argument("--step0", type=float, default=0.08); ev(p)
    p.set_defaults(func=cmd_climb)

    p = sub.add_parser("continue-lineage"); p.add_argument("--run-dir", required=True)
    p.add_argument("--lineage-id", required=True); p.add_argument("--worker-id", required=True)
    p.add_argument("--seed-family"); p.add_argument("--seed-candidate")
    p.add_argument("--max-generations", type=int, default=30); p.add_argument("--step0", type=float, default=0.10)
    p.add_argument("--min-step", type=float, default=0.006); p.add_argument("--patience", type=int, default=6)
    p.add_argument("--no-branch", action="store_true"); ev(p)
    p.set_defaults(func=cmd_continue_lineage)

    p = sub.add_parser("checkpoint"); p.add_argument("--run-dir", required=True)
    p.add_argument("--basin-radius", type=float, default=0.08); p.set_defaults(func=cmd_checkpoint)

    p = sub.add_parser("candidates"); p.add_argument("--run-dir", required=True)
    p.add_argument("--limit", type=int, default=20); p.add_argument("--valid-only", action="store_true")
    p.add_argument("--family", choices=("mixed", "photonic"))
    p.add_argument("--target", choices=("mixed", "photonic")); p.set_defaults(func=cmd_candidates)

    p = sub.add_parser("summary"); p.add_argument("--run-dir", required=True); p.set_defaults(func=cmd_summary)

    args = parser.parse_args(argv)
    args.func(args)
    return 0
