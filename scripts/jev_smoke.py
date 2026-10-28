#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.jev_policy import call_jev


def main() -> int:
    # Real Jev request, but tiny and side-effect free. This verifies the exact path
    # used by the DiHiggs policy rather than relying on a model-list endpoint.
    payload = {
        "state": {"purpose": "credential and connectivity smoke test", "expected": "service reachable"},
        "questions": {
            "reachable": {
                "type": "noul",
                "instructions": "Is this state clearly a connectivity smoke test?",
            }
        },
    }
    result = call_jev(payload, retries=0)
    print(json.dumps({
        "ok": True,
        "transport": "local_proxy" if os.environ.get("JEV_PROXY_URL") else "vercel_direct",
        "model": result.get("model", "typesafe-ai/jev"),
        "answers": result.get("answers"),
        "usage": result.get("usage"),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
