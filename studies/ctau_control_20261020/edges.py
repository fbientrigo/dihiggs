"""Locate the lower validity edge of the fixed-X trajectories (BFB_general transitions) by bisection in ln tan beta.

  <venv>/bin/python -I scripts/edges.py   -> raw/edges.jsonl, raw/edges.json
"""
import json, math, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ctl

rows, out = [], {}
order = {"FAIL": 0, "INCONCLUSIVE": 1, "PASS": 2}
for fam, ref in ctl.references().items():
    X = ref["X"]
    if float(X) == 0:
        continue
    for level in ("INCONCLUSIVE", "PASS"):  # smallest t whose BFB status is at least `level`
        lo, hi = math.log(10.0), math.log(float(ref["tan_beta"]))
        for _ in range(28):
            mid = (lo + hi) / 2
            r = ctl.evaluate(X, math.exp(mid), fam)
            r.update(experiment="edge_bisection", note=f"BFB>={level}")
            rows.append(r)
            lo, hi = (mid, hi) if order[r["tree_BFB_general"]] < order[level] else (lo, mid)
        r = ctl.evaluate(X, math.exp(hi), fam)
        out[f"{fam}_BFB_{level}_from"] = dict(tan_beta=r["tan_beta"], lambda6=r["lambda6"], Z7=r["Z7"], ctau_mm=r["ctau_mm"],
                                              BFB=r["tree_BFB_general"], status=r["status"])
        print(fam, level, out[f"{fam}_BFB_{level}_from"], flush=True)
with open(ctl.HERE / "raw/edges.jsonl", "w") as f:
    for r in rows:
        f.write(json.dumps(r, default=str) + "\n")
json.dump(out, open(ctl.HERE / "raw/edges.json", "w"), indent=1)
