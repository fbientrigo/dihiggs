from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

from .archive import Archive
from .contract import CONTRACT_PATH, ROOT, SCHEMA_VERSION, ContractError, normalize_proposal, load_contract
from .evaluator import CanonicalEvaluator, family_proposals
from .helpers import deduplicate, perturb, random_candidates
from .ledger import Ledger
from .score import score_components


def _load(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _emit(document: Any) -> None:
    print(json.dumps(document, indent=2, sort_keys=True, allow_nan=False))


def _record_evaluation(ledger: Ledger, evaluator: CanonicalEvaluator, proposal: dict[str, Any], run_dir: Path, ordinal: int) -> dict[str, Any]:
    try:
        normalized = normalize_proposal(proposal)
    except ContractError as error:
        ledger.append({"event": "PROPOSAL", "lifecycle": "PROPOSED", "proposal_id": proposal.get("proposal_id"), "source_candidate": proposal})
        failure = evaluator.evaluate(proposal, run_dir / "attempts" / f"{ordinal:06d}_malformed")
        ledger.append({"event": "EVALUATION_FAILED", "lifecycle": "FAILED", **failure})
        return failure
    ledger.append({"event": "PROPOSAL", "lifecycle": "PROPOSED", "candidate_id": normalized["candidate_id"], "proposal_id": normalized["proposal_id"], "source_candidate": normalized})
    ledger.append({"event": "EVALUATION_STARTED", "lifecycle": "EVALUATING", "candidate_id": normalized["candidate_id"], "proposal_id": normalized["proposal_id"]})
    result = evaluator.evaluate(proposal, run_dir / "attempts" / f"{ordinal:06d}_{normalized['candidate_id']}")
    ledger.append({"event": "EVALUATION_TERMINATED" if result.get("status") == "TERMINATED" else "EVALUATION_FAILED", "lifecycle": result.get("status", "FAILED"), **result})
    return result


def _family_record(results: list[dict[str, Any]], proposals: list[dict[str, Any]]) -> dict[str, Any]:
    normalized = [normalize_proposal(item) for item in proposals]
    ids = [item["candidate_id"] for item in normalized]
    family_id = "family_" + hashlib.sha256("|".join(ids).encode()).hexdigest()[:16]
    families = {result.get("classification", {}).get("family") for result in results}
    classification = results[1].get("classification", {}) if len(results) > 1 else {}
    valid = len(results) == 3 and all(result.get("status") == "TERMINATED" and result.get("validity_gate") for result in results)
    same_x = len({round(float(item["derived"]["X"]), 12) for item in normalized}) == 1
    all_masses = [item["parameters"]["mH_GeV"] for item in normalized] == [150.0, 200.0, 250.0]
    cross = 1.0 if valid and same_x and all_masses else 0.0
    anchor = results[1] if len(results) > 1 else results[0]
    total = score_components(valid=valid, classification=classification, ctau_mm=anchor.get("ctau_mm"), cross_mass_consistency=cross)
    promotable = valid and same_x and all_masses and classification.get("status") == "CLASSIFIED" and len(families) == 1
    return {
        "schema_version": "dihiggs.llp.family_evaluation.v1", "family_id": family_id,
        "family": classification.get("family"), "status": "PROMOTABLE" if promotable else ("BLOCKED_SCIENTIFIC_DECISION" if classification.get("status") == "BLOCKED_SCIENTIFIC_DECISION" else "REJECTED"),
        "validity_gate": valid, "same_X": same_x, "required_masses": all_masses, "X": normalized[1]["derived"]["X"] if len(normalized) > 1 else None,
        "cross_mass_consistency": cross, "score_components": total, "total_score": total["total_score"],
        "anchor": normalized[1] if len(normalized) > 1 else normalized[0], "members": results,
        "provenance": {"source_lineage": [item.get("proposal_id") for item in normalized]},
    }


def cmd_evaluate(args: argparse.Namespace) -> None:
    run_dir = Path(args.run_dir).resolve()
    evaluator = CanonicalEvaluator(ROOT, Path(args.evaluator))
    ledger = Ledger(run_dir)
    result = _record_evaluation(ledger, evaluator, _load(Path(args.candidate)), run_dir, int(time.time_ns() % 1000000))
    _emit(result)


def cmd_validate_family(args: argparse.Namespace) -> None:
    run_dir = Path(args.run_dir).resolve()
    evaluator = CanonicalEvaluator(ROOT, Path(args.evaluator))
    ledger = Ledger(run_dir)
    anchor = _load(Path(args.family))
    proposals = family_proposals(anchor)
    results = [_record_evaluation(ledger, evaluator, proposal, run_dir, index) for index, proposal in enumerate(proposals)]
    record = _family_record(results, proposals)
    ledger.append({"event": "FAMILY_EVALUATION", "lifecycle": "TERMINATED", **record})
    promotion = Archive(run_dir).consider(record)
    record["archive_decision"] = promotion
    _emit(record)


def cmd_request(args: argparse.Namespace) -> None:
    proposals = random_candidates(args.count, args.seed, args.worker_id)
    if args.evaluated_ids:
        evaluated = set(_load(Path(args.evaluated_ids)))
        proposals = deduplicate(proposals, evaluated)
    _emit({"schema_version": SCHEMA_VERSION, "contract": str(CONTRACT_PATH), "proposals": proposals})


def cmd_perturb(args: argparse.Namespace) -> None:
    parent = _load(Path(args.parent))
    _emit({"schema_version": SCHEMA_VERSION, "proposals": perturb(parent, args.count, args.seed, args.radius, args.worker_id)})


def cmd_summary(args: argparse.Namespace) -> None:
    run_dir = Path(args.run_dir).resolve()
    _emit({"contract": str(CONTRACT_PATH), "contract_version": SCHEMA_VERSION, **Ledger(run_dir).summary(), "archive": Archive(run_dir).load()})


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Deterministic DiHiggs LLP benchmark-search substrate")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("evaluate", help="evaluate one normalized proposal")
    p.add_argument("--candidate", required=True); p.add_argument("--run-dir", required=True); p.add_argument("--evaluator", default=str(ROOT / "dihiggs/app/DihiggsPointV2Evaluator")); p.set_defaults(func=cmd_evaluate)
    p = sub.add_parser("validate-family", help="derive and evaluate 150/200/250 GeV Q,S family")
    p.add_argument("--family", required=True); p.add_argument("--run-dir", required=True); p.add_argument("--evaluator", default=str(ROOT / "dihiggs/app/DihiggsPointV2Evaluator")); p.set_defaults(func=cmd_validate_family)
    p = sub.add_parser("request-candidates", help="deterministic random proposals")
    p.add_argument("--count", type=int, default=1); p.add_argument("--seed", type=int, required=True); p.add_argument("--worker-id", default="orca-worker"); p.add_argument("--evaluated-ids"); p.set_defaults(func=cmd_request)
    p = sub.add_parser("perturb", help="deterministic local perturbations")
    p.add_argument("--parent", required=True); p.add_argument("--count", type=int, default=1); p.add_argument("--seed", type=int, required=True); p.add_argument("--radius", type=float, default=0.05); p.add_argument("--worker-id", default="orca-worker"); p.set_defaults(func=cmd_perturb)
    p = sub.add_parser("summary", help="inspect derived run summary")
    p.add_argument("--run-dir", required=True); p.set_defaults(func=cmd_summary)
    p = sub.add_parser("contract", help="print frozen contract")
    p.set_defaults(func=lambda args: _emit(load_contract()))
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except Exception as error:
        print(json.dumps({"status": "FAILED", "error": str(error)}), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
