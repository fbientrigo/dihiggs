from __future__ import annotations

import argparse
import json
from pathlib import Path

from .cli import _record_evaluation
from .contract import ROOT
from .evaluator import CanonicalEvaluator
from .helpers import random_candidates
from .ledger import Ledger
from .archive import Archive


def main() -> int:
    parser = argparse.ArgumentParser(description="small deterministic substrate smoke campaign; never an LLM search")
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--seed", type=int, default=20260825)
    parser.add_argument("--evaluator", default=str(ROOT / "dihiggs/app/DihiggsPointV2Evaluator"))
    args = parser.parse_args()
    run_dir = Path(args.run_dir).resolve()
    evaluator = CanonicalEvaluator(ROOT, Path(args.evaluator))
    ledger = Ledger(run_dir)
    proposals = random_candidates(2, args.seed, "dry-campaign")
    results = [_record_evaluation(ledger, evaluator, item, run_dir, index) for index, item in enumerate(proposals)]
    report = {
        "schema_version": "dihiggs.llp.dry_campaign.v1", "seed": args.seed, "count": len(results),
        "results": [{"candidate_id": item.get("candidate_id"), "status": item.get("status"), "validity_gate": item.get("validity_gate"), "classification": item.get("classification"), "ctau_mm": item.get("ctau_mm"), "failure_reason": item.get("failure_reason")} for item in results],
        "summary": ledger.summary(), "archive": Archive(run_dir).load(),
    }
    (run_dir / "dry_campaign.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

