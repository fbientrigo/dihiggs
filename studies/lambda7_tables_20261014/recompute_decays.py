"""Re-evaluate the selected lambda7 benchmarks with the frozen 2HDMC backend of the
2026-10-11 consolidation (imported read-only, nothing there is modified) and write
raw/decays_recomputed.json.

  /tmp/dihiggs-t2-t4-venv/bin/python -I scripts/recompute_decays.py

Each point is rebuilt from its stored inputs (tanb, lambda1, lambda6, rho, mA=mH+), so the
lambda7 / Z7 / M^2 / g columns come from the constructor, not from the CSV. The lambda7=0
control (B02_L7ZERO) is the C2 control of the consolidation: lambda1, lambda6, spectrum fixed,
rho solved so that lambda7=0 (g changes).
"""
import os
import csv, json, sys
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(os.environ.get("DIHIGGS_ROOT") or Path(__file__).resolve().parents[3])
HERE = ROOT / "reports/lambda7_tables_20261014"  # output package (lake), not code
CONS = ROOT / "dihiggs/campaigns/lambda7_consolidation/20261011T1530Z"
sys.path.insert(0, str(CONS))
import consolidate as C  # noqa: E402  (pulls stage3/study_d from the frozen campaign)

BM = {r["benchmark_id"]: r for r in csv.DictReader(open(CONS / "BENCHMARKS.csv"))}
SELECT = ["B02", "B05", "B09", "B12", "B17", "B13", "B14",
          "B03", "B10", "B18"]


def rho_lambda7_zero(b):
    with localcontext() as c:
        c.prec = 60
        t, v0, mh, mp, L1, l6 = (Decimal(b[k]) for k in ("tan_beta", "v_GeV", "mh_GeV", "mphi_GeV", "lambda1", "lambda6"))
        a = mh * mh / v0 / v0
        return str(1 - (mh * mh / 2 + v0 * v0 * (3 * l6 / (2 * t) + (L1 - a) / (t * t))) / (mp * mp))


def t1_residual(fp):
    """|lambda7(T1 identity) - lambda7(constructor)|, absolute."""
    i, v = fp.inputs, fp.values
    with localcontext() as c:
        c.prec = 60
        t, v0, mh, mp = (Decimal(i[k]) for k in ("tan_beta", "v_GeV", "mh_GeV", "mphi_GeV"))
        rho = Decimal(v.get("rho", "0")); a = mh * mh / v0 / v0
        dd = (mp * mp * (1 - rho) - mh * mh / 2) / v0 / v0
        l7 = -2 * dd / t + 3 * Decimal(v["lambda6"]) / t**2 + 2 * (Decimal(i["lambda1_ref"]) - a) / t**3
        return float(abs(l7 - Decimal(v["lambda7"])))


def one(pid, src, rho, study):
    b = BM[src]
    row = C.row_of(study, b["mphi_GeV"], b["tan_beta"], lam6=b["lambda6"], L1=b["lambda1"], Delta=b["Delta_GeV"], rho=rho)
    fp, cons, st, d = C.evaluate(row)
    v, i = fp.values, fp.inputs
    return dict(id=pid, source=src, inputs=i, values={k: v[k] for k in v if not isinstance(v[k], (dict, list))},
                constraints={k: cons[k]["status"] for k in C.T3.ORDER}, tree_status=st,
                t1_lambda7_abs_residual=t1_residual(fp),
                width_total_GeV=d["width_total_phi_GeV"], ctau_mm=d["ctau_mm"], widths_GeV=d["width_channels_GeV"],
                physical_state_ordering=d.get("physical_state_ordering"))


out = [one(s, s, BM[s]["rho"], BM[s]["study"]) for s in SELECT]
out.append(one("B02_L7ZERO", "B02", rho_lambda7_zero(BM["B02"]), "D"))
out.append(one("B17_L7ZERO", "B17", rho_lambda7_zero(BM["B17"]), "D"))
out.append(one("B18_L7ZERO", "B18", rho_lambda7_zero(BM["B18"]), "D"))
for r in out:
    print(r["id"], r["tree_status"], "l7=%.3e" % float(r["values"]["lambda7"]), "T1res=%.1e" % r["t1_lambda7_abs_residual"],
          "g=%.4f" % float(r["values"]["g_hphiphi_GeV"]), "ctau=%.4g" % r["ctau_mm"])
(HERE / "raw").mkdir(exist_ok=True)
json.dump(out, open(HERE / "raw/decays_recomputed.json", "w"), indent=1, default=str)
