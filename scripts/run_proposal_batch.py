#!/usr/bin/env python3
"""CLI front end for dihiggs.app.orchestrator.proposal_batch (issue #72).

    python scripts/run_proposal_batch.py \\
        --proposals proposals.csv --outdir out --campaign-id my_campaign \\
        --run-id my_run

See docs/contracts/adaptive_proposal_batch_v1.yaml for the input contract.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dihiggs.app.orchestrator.proposal_batch import run_proposal_batch  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proposals", required=True, type=Path, help="proposals.csv path")
    parser.add_argument("--outdir", required=True, type=Path, help="output root directory")
    parser.add_argument("--campaign-id", required=True)
    parser.add_argument("--run-id", required=True, dest="run_name")
    parser.add_argument("--exec", dest="executable", type=Path,
                         default=ROOT / "dihiggs/app/DihiggsPointV2Evaluator",
                         help="DihiggsPointV2Evaluator binary (default: dihiggs/app/DihiggsPointV2Evaluator)")
    parser.add_argument("--timeout-s", type=float, default=None)
    parser.add_argument("--keep-raw", action="store_true",
                         help="keep each proposal's raw single-row evaluator CSV alongside attempts.csv")
    args = parser.parse_args()

    try:
        manifest = run_proposal_batch(
            executable=args.executable, proposals_csv=args.proposals, outdir=args.outdir,
            campaign_id=args.campaign_id, run_name=args.run_name, repo_root=ROOT,
            timeout_s=args.timeout_s, keep_raw=args.keep_raw,
        )
    except RuntimeError as error:
        print(f"[proposal_batch] FAILED: {error}", file=sys.stderr)
        return 1

    print(f"[proposal_batch] {manifest['status']} counts={manifest['counts']}")
    return 0 if manifest["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
