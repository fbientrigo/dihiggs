"""Shared evaluator for the c tau control study (read-only reuse of the 2026-10-19 machinery).

A point is (m_phi, X, tan beta); lambda6 = X / tan beta is formed in 60-digit Decimal arithmetic, then the
2026-10-19 constructor (core.evaluate(full=True): aligned-lambda7 T1 constructor at 50 digits, 80-digit and
closed-form cross-routes, 8 tree tests), the 2HDMC analytic bridge, the native C++ decay inventory and the
independent LO loop widths are run for EVERY point.  Nothing in llp_families_20261019 is written: the
evaluation cache of core is redirected to this directory and starts empty.
"""
from __future__ import annotations

import os
import csv, math, sys
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(os.environ.get("DIHIGGS_ROOT") or Path(__file__).resolve().parents[3])
HERE = ROOT / "reports/ctau_control_20261020"  # output package (lake), not code
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "llp_families_20261019"))
import core  # noqa: E402
import finalize as F  # noqa: E402

core.CACHE = HERE / "raw/eval_cache.jsonl"  # never the October cache
M = 250
REF_IDS = {"hadronic": "A6", "mixed": "M250p", "photonic": "P250n"}
# gates of finalize.gates that are about the c tau = 100 mm benchmark choice or the family label, not validity
NOT_VALIDITY = {"ctau_80_120_mm", "composition_class"}
TREE = ("spectrum_reconstruction", "stationarity", "alignment", "local_minimum_tree", "BFB_general",
        "global_EW_minimum", "full_tree_unitarity", "quartic_screen_4pi")


def references():
    """Full-precision 250 GeV references from the October validated_benchmarks.csv."""
    rows = {r["candidate_id"]: r for r in csv.DictReader(open(OCT / "validated_benchmarks.csv"))}
    return {f: rows[cid] for f, cid in REF_IDS.items()}


def l6_of(X, t):
    with localcontext() as c:
        c.prec = 60
        x, tt = Decimal(X), Decimal(t)
        return "0" if x == 0 else str(x / tt)


def tstr(t):
    """tan beta as an exact decimal string with 6 significant digits."""
    return f"{float(f'{t:.5e}'):.0f}" if t >= 1e5 else f"{float(f'{t:.5e}'):.6g}"


def evaluate(X, t, family=None):
    t = t if isinstance(t, str) else tstr(t)
    l6 = l6_of(X, t)
    r = F.full_eval(M, t, l6)
    r["accounting"] = r.get("accounting") or core.accounting(r)
    g = F.gates(r, family)
    val = {k: v for k, v in g.items() if k not in NOT_VALIDITY}
    with localcontext() as c:
        c.prec = 60
        x_dev = abs(Decimal(r["X_lambda6_tanbeta"]) - Decimal(X))
    val["X_preserved_1e-40"] = (x_dev < Decimal("1e-40"), str(x_dev))
    tree = {k: r[f"test_{k}"] for k in TREE}
    failed = [k for k, v in val.items() if not v[0]]
    if any(s == "FAIL" for s in tree.values()):
        status = "FAIL"
    elif failed:
        status = "INCONCLUSIVE"
    else:
        status = "PASS"
    acc = r["accounting"]
    out = dict(family=family, mphi_GeV=M, X=X, tan_beta=t, lambda6=l6, lambda7=r["lambda7"],
               lambda1=core.L1, lambda2=r["lambda2"], lambda3=r["lambda3"], lambda4=r["lambda4"], lambda5=r["lambda5"],
               Z7=r["Z7"], m12_sq_GeV2=r["m12_sq_GeV2"], Y2_GeV2=r["Y2_GeV2"], g_hphiphi_GeV=r["g_hphiphi_GeV"],
               Gamma_GeV=r["Gamma_GeV"], ctau_mm=r["ctau_mm"],
               Gamma_bb_GeV=r["width_bb_GeV"], Gamma_gluglu_GeV=r["width_gg_GeV"],
               Gamma_Zgamma_GeV=r["width_Zgamma_GeV"], Gamma_gammagamma_GeV=r["width_diphoton_GeV"],
               BR_bb=r["BR_bb"], BR_gluglu=r["BR_gg"], BR_Zgamma=r["BR_Zgamma"], BR_gammagamma=r["BR_diphoton"],
               BR_tautau=acc["BR"].get("tautau"), BR_cc=acc["BR"].get("cc"), BR_ttstar=acc["BR"].get("ttstar"),
               BR_other_terms=acc["BR_other_terms"], BR_inventory_sum_minus_1=sum(acc["BR"].values()) + acc["BR_other_terms"] - 1,
               B_h=r["B_h"], B_gamma=r["B_gamma"], family_class=r["family_class"],
               max_abs_quartic=max(abs(float(r[k])) for k in ("lambda2", "lambda3", "lambda4", "lambda5", "lambda7", "Z7")),
               lambda1_from_lambda7_residual=r["lambda1_from_lambda7_abs_residual"],
               lambda7_identity_rel_residual=r["lambda7_identity_rel_residual"],
               backend_mass_max_abs_GeV=r["backend_mass_max_abs_GeV"], backend_sba_minus_1=r["backend_sba"] - 1,
               ratio_gaga_backend_over_LO=r["ratio_gaga_backend_over_LO"], ratio_Zga_backend_over_LO_os=r["ratio_Zga_backend_over_LO_os"],
               CV_vs_minus_vZ7_rel=r["backend_CphiHpHm_vs_minus_vZ7_rel"],
               **{f"tree_{k}": v for k, v in tree.items()},
               status=status, failed_gates=";".join(failed),
               failed_gate_details=";".join(f"{k}={val[k][1]}" for k in failed))
    return out


def dev(rec, ref):
    """Composition drift w.r.t. the reference record."""
    keys = ("BR_bb", "BR_gluglu", "BR_Zgamma", "BR_gammagamma")
    return dict(dB_gamma=rec["B_gamma"] - ref["B_gamma"], dB_h=rec["B_h"] - ref["B_h"],
                max_abs_dBR=max(abs(rec[k] - ref[k]) for k in keys),
                max_rel_dBR=max(abs(rec[k] / ref[k] - 1) for k in keys))
