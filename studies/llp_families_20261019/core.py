"""Point evaluator for the Y2=0, lambda7!=0 LLP families (read-only use of the frozen backend).

Primary route : stage3.build(study F)  = aligned_lambda7 T1 constructor (50 digits) + 8 tree tests
Route 2       : study_d.build_generalized_family(rho=0) at 80 digits (closed form) and
                study_d.linear_route (numerical inversion of the Z1..Z7 map, 120 digits)
Identity      : lambda7 = -2d/t + 3 lambda6/t^2 + 2(L-a)/t^3 in mpmath (100 digits)
Decays        : 2HDMC analytic-bridge (C++, double), same argv as the October campaign.

Run with the compute venv:  <venv>/bin/python -I scripts/<driver>.py
"""
from __future__ import annotations

import os
import json, math, subprocess, sys
from decimal import Decimal, localcontext
from pathlib import Path

import mpmath as mp

ROOT = Path(os.environ.get("DIHIGGS_ROOT") or Path(__file__).resolve().parents[3])
CONS = ROOT / "dihiggs/campaigns/lambda7_consolidation/20261011T1530Z"
sys.path.insert(0, str(CONS))
import consolidate as C  # noqa: E402  frozen; imports stage3 / study_d
S, T3 = C.S, C.T3

HERE = ROOT / "reports/llp_families_20261019"  # output package (lake), not code
ACCOUNTING = ROOT / "dihiggs/campaigns/t2_t4_20261008T043858Z/decay_accounting/DecayAccountingAudit"
HBARC = 1.973269804e-13  # GeV mm
MH, V, L1, DELTA = "125.13", "246", "1.2", "100"
CACHE = HERE / "raw/eval_cache.jsonl"
_cache: dict | None = None


def classify(Bh, Bg):
    if Bh >= 0.75:
        return "hadronic"
    if Bg >= 0.95 and Bh <= 0.03:
        return "photonic"
    if 0.5 <= Bg < 0.95 and Bh >= 0.05:
        return "mixed"
    return "transitional"


def identity_checks(m, t, l6, lam7_str, lam1_str=L1):
    with mp.workdps(100):
        t_, l6_, L = mp.mpf(t), mp.mpf(l6), mp.mpf(lam1_str)
        a = mp.mpf(MH) ** 2 / mp.mpf(V) ** 2
        d = (mp.mpf(m) ** 2 - mp.mpf(MH) ** 2 / 2) / mp.mpf(V) ** 2
        l7_id = -2 * d / t_ + 3 * l6_ / t_**2 + 2 * (L - a) / t_**3
        l7 = mp.mpf(lam7_str)
        l1_back = a + d * t_**2 - mp.mpf(3) / 2 * l6_ * t_ + l7 * t_**3 / 2
        return dict(lambda7_identity=mp.nstr(l7_id, 30),
                    lambda7_identity_rel_residual=float(abs(l7 - l7_id) / abs(l7_id)) if l7_id else float("inf"),
                    lambda1_from_lambda7_abs_residual=float(abs(l1_back - L)),
                    lambda7_leading_term=mp.nstr(-2 * d / t_, 20))


def bridge(fp, qcd=True):
    return S.bridge_decays(fp, qcd_on=qcd)


def evaluate(m, t, l6, full=False):
    """m [GeV] int/str, t str, l6 str. Returns a flat record (cached on (m,t,l6,full))."""
    global _cache
    key = f"{m}|{t}|{l6}|{int(full)}"
    if _cache is None:
        _cache = {}
        if CACHE.exists():
            for line in open(CACHE):
                r = json.loads(line); _cache[r["key"]] = r
    if key in _cache:
        return _cache[key]
    row = C.row_of("F", m, t, lam6=l6, L1=L1, Delta=DELTA)
    fp, cons, st = T3.build(row, 50)
    v, i = fp.values, fp.inputs
    rec = dict(key=key, mphi_GeV=float(m), tan_beta=str(t), lambda6=str(l6), tree_status=st,
               **{f"test_{k}": cons[k]["status"] for k in T3.ORDER})
    for k in ("lambda1", "lambda2", "lambda3", "lambda4", "lambda5", "lambda7", "Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7",
              "M_sq_GeV2", "Y2_GeV2", "m12_sq_GeV2", "m11_sq_GeV2", "m22_sq_GeV2", "g_hphiphi_GeV",
              "X_lambda6_tanbeta", "CphiHpHm_potential_GeV", "phi_mass_reconstructed_GeV", "h_mass_reconstructed_GeV",
              "A_mass_reconstructed_GeV", "Hp_mass_reconstructed_GeV"):
        rec[k] = str(v[k])
    with localcontext() as c:
        c.prec = 60
        rec["Y2_is_exact_zero"] = Decimal(rec["Y2_GeV2"]) == 0
        rec["M2_minus_mh2_over_2"] = str(Decimal(rec["M_sq_GeV2"]) - Decimal(MH) ** 2 / 2)
        rec["g_minus_2m2_over_v"] = str(Decimal(rec["g_hphiphi_GeV"]) - 2 * Decimal(m) ** 2 / Decimal(V))
    rec.update(identity_checks(m, t, l6, rec["lambda7"]))
    d = bridge(fp)
    tot = d["width_total_phi_GeV"]
    rec.update(Gamma_GeV=tot, ctau_mm=HBARC / tot,
               **{f"BR_{k}": d["width_channels_GeV"][k] / tot for k in C.CH},
               BR_reported_sum=d["BR_sum_reported"],
               backend_mphi_GeV=d["m_H_sorted_backend_GeV"], backend_mh_GeV=d["m_h_sorted_backend_GeV"],
               backend_mA_GeV=d["m_A_backend_GeV"], backend_mHp_GeV=d["m_Hp_backend_GeV"],
               backend_sba=d["sin_beta_minus_alpha_backend"], backend_v=d["backend_v_GeV"],
               backend_g_hphiphi=d["g_hphi_phi_backend_GeV"], backend_g_phiHpHm=d["g_phi_HpHm_backend_GeV"],
               backend_kappa_top=d["kappa_phi_top_backend_signed"], backend_Z7=d["Z7"],
               width_diphoton_GeV=d["width_diphoton_GeV"], width_Zgamma_GeV=d["width_Zgamma_GeV"],
               width_gg_GeV=d["width_gg_GeV"], width_bb_GeV=d["width_bb_GeV"])
    rec["B_h"] = rec["BR_bb"] + rec["BR_gg"]
    rec["B_gamma"] = rec["BR_diphoton"] + rec["BR_Zgamma"]
    rec["family_class"] = classify(rec["B_h"], rec["B_gamma"])
    if full:
        rec.update(route2_audit(m, t, l6, rec))
    _cache[key] = rec
    with open(CACHE, "a") as f:
        f.write(json.dumps(rec, default=str) + "\n")
    return rec


def route2_audit(m, t, l6, rec):
    """Independent constructor routes at higher precision + high-precision mass reconstruction."""
    out = {}
    fpD = S.build_generalized_family(tan_beta=str(t), lambda1_ref=L1, mh_GeV=MH, mphi_GeV=str(m),
                                     mA_GeV=str(int(m) + int(DELTA)), mHp_GeV=str(int(m) + int(DELTA)), v_GeV=V,
                                     rho="0", lambda6=str(l6), precision=80)
    aud = S.reconstruction_audit(fpD, 80)
    with mp.workdps(60):
        diffs = {k: abs(mp.mpf(fpD.values[k]) - mp.mpf(rec[k])) for k in
                 ("lambda2", "lambda3", "lambda4", "lambda5", "lambda7", "Z7", "g_hphiphi_GeV")}
    out["routeF50_vs_D80_max_abs"] = float(max(diffs.values()))
    out["routeF50_vs_D80_lambda7_rel"] = float(diffs["lambda7"] / abs(mp.mpf(rec["lambda7"])))
    out["linear_route_agreement_max_abs"] = aud["route_agreement_max_abs_float"]
    out["g_closure_minus_formula_GeV"] = float(mp.mpf(aud["g_residual_abs_GeV"]))
    out["Y2_identity_abs"] = float(mp.mpf(aud["Y2_identity_abs"]))
    out["closure_core_pass"] = bool(aud["closure_core_pass"])
    # precision stability of the primary constructor: 50 vs 80 digits
    fp80, _, st80 = T3.build(C.row_of("F", m, t, lam6=l6, L1=L1, Delta=DELTA), 80)
    with mp.workdps(90):
        out["prec50_vs_80_lambda7_rel"] = float(abs(mp.mpf(fp80.values["lambda7"]) - mp.mpf(rec["lambda7"])) / abs(mp.mpf(rec["lambda7"])))
        out["prec50_vs_80_Z7_rel"] = float(abs(mp.mpf(fp80.values["Z7"]) - mp.mpf(rec["Z7"])) / abs(mp.mpf(rec["Z7"])))
    out["prec80_tree_status"] = st80
    # high-precision mass reconstruction from the 50-digit quartics (independent of the double backend)
    with mp.workdps(80):
        hp = [abs(mp.mpf(rec[k]) - mp.mpf(x)) for k, x in
              (("phi_mass_reconstructed_GeV", m), ("h_mass_reconstructed_GeV", MH),
               ("A_mass_reconstructed_GeV", int(m) + int(DELTA)), ("Hp_mass_reconstructed_GeV", int(m) + int(DELTA)))]
    out["hp_mass_reconstruction_max_abs_GeV"] = float(max(hp))
    out["backend_mass_max_abs_GeV"] = max(abs(rec["backend_mphi_GeV"] - float(m)), abs(rec["backend_mh_GeV"] - float(MH)),
                                          abs(rec["backend_mA_GeV"] - float(m) - 100), abs(rec["backend_mHp_GeV"] - float(m) - 100))
    out["backend_CphiHpHm_vs_minus_vZ7_rel"] = abs(rec["backend_g_phiHpHm"] - float(rec["CphiHpHm_potential_GeV"])) / abs(float(rec["CphiHpHm_potential_GeV"]))
    out["backend_kappa_top_vs_minus_cotb_rel"] = abs(rec["backend_kappa_top"] + 1 / float(t)) * float(t)
    out["backend_g_hphiphi_vs_constructor_rel"] = abs(abs(rec["backend_g_hphiphi"]) - float(rec["g_hphiphi_GeV"])) / float(rec["g_hphiphi_GeV"])
    return out


def accounting(rec):
    """Complete native DecayTable inventory (all fermion pairs incl. off-shell t t*). The double-precision
    backend mass readback differs from the inputs by O(t*eps) at large tan(beta); the strict identity probe is run
    first, then (if it fails only on that) rerun with the measured native masses - recorded explicitly."""
    lam = [rec[f"lambda{k}"] for k in range(1, 8)]
    base = [*lam, rec["m12_sq_GeV2"], rec["tan_beta"], V, "1", "2"]
    run = lambda mh, mp_, vb: json.loads(subprocess.run([str(ACCOUNTING), *base, str(mh), str(mp_), str(vb)],
                                                        check=True, text=True, capture_output=True, timeout=60).stdout)
    probe = run(MH, rec["mphi_GeV"], rec["backend_v"])
    res = probe if probe.get("state_identity_ok") else run(repr(probe["m_h_index_GeV"]), repr(probe["m_phi_index_GeV"]), repr(probe["backend_v_GeV"]))
    names = {("hdd", 1, 1): "dd", ("hdd", 2, 2): "ss", ("hdd", 3, 3): "bb", ("huu", 1, 1): "uu", ("huu", 2, 2): "cc",
             ("huu", 3, 3): "ttstar", ("hll", 1, 1): "ee", ("hll", 2, 2): "mumu", ("hll", 3, 3): "tautau",
             ("hvv", 1, 0): "gammagamma", ("hZga", 0, 0): "Zgamma", ("hgg", 0, 0): "gluongluon"}
    inv, other = {}, 0.0
    for term in res["terms"]:
        w = term["width_GeV"] or 0.0
        n = names.get((term["family"], term["i"], term["j"]))
        if n:
            inv[n] = w
        else:
            other += w
    tot = res["total_width_GeV"]
    return dict(strict_identity_ok=bool(probe.get("state_identity_ok")), native_mphi_GeV=res["m_phi_index_GeV"],
                total_width_GeV=tot, ordered_term_sum_GeV=res["ordered_term_sum_GeV"],
                BR={k: w / tot for k, w in inv.items()}, BR_other_terms=other / tot,
                total_vs_bridge_rel=abs(tot - rec["Gamma_GeV"]) / rec["Gamma_GeV"])
