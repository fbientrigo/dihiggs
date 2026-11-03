"""Phase A: adversarial backend audit -> audit/backend_audit.json  (run: <venv>/bin/python -I scripts/backend_audit.py)"""
import cmath, json, math, subprocess, sys
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import core
from loops import *  # noqa: F401,F403

OUT = core.HERE / "audit/backend_audit.json"
BRIDGE = core.ROOT / "dihiggs/2hdmc/analytic-bridge"


def raw_bridge(rec, override=None):
    lam = [rec[f"lambda{k}"] for k in range(1, 8)]
    for k, x in (override or {}).items():
        lam[k - 1] = x
    m = rec["mphi_GeV"]
    argv = [str(BRIDGE), "analytic", core.MH, str(m), str(m + 100), str(m + 100), rec["tan_beta"], *lam,
            rec["m12_sq_GeV2"], core.V, "1"]
    return json.loads(subprocess.run(argv, text=True, capture_output=True, timeout=60).stdout or "{}")


audit = {}

# 1) lambda7 forced to zero with everything else frozen: the point must stop being the requested aligned spectrum
r = core.evaluate(200, "3000000", "-2.6364e-4")
z = raw_bridge(r, {7: "0"})
audit["T1_force_lambda7_zero"] = dict(
    reference=dict(mphi=r["backend_mphi_GeV"], sba=r["backend_sba"], ctau_mm=r["ctau_mm"]),
    lambda7_zeroed=dict(constructor_ok=z.get("constructor_ok"), ordering=z.get("physical_state_ordering"),
                        m_H=z.get("m_H_sorted_backend_GeV"), m_h=z.get("m_h_sorted_backend_GeV"),
                        sba=z.get("sin_beta_minus_alpha_backend"), failure=z.get("exact_failure_reason")),
    lambda1_required_by_identity_if_lambda7_0=float(mp.mpf(core.identity_checks(200, "3000000", "-2.6364e-4", "0")["lambda7_identity"]) * 0)
    + float((mp.mpf(core.MH) ** 2 / 246**2) + (200**2 - mp.mpf(core.MH) ** 2 / 2) / 246**2 * mp.mpf(3e6) ** 2 - 1.5 * mp.mpf("-2.6364e-4") * 3e6),
    verdict=None)
audit["T1_force_lambda7_zero"]["verdict"] = (
    "PASS (backend reacts: input spectrum/alignment not reproduced)" if (not z.get("constructor_ok")) or
    abs((z.get("m_H_sorted_backend_GeV") or 0) - 200) > 1 or abs((z.get("sin_beta_minus_alpha_backend") or 0) - 1) > 1e-6
    else "FAIL (backend insensitive to lambda7)")

# 2) lambda6 scan at fixed t: Z7 and photonic widths must respond; amplitude linear in Z7 (3-point fit, 4th predicted)
pts = [core.evaluate(200, "3000000", l6) for l6 in ("-2.6364e-4", "-1e-4", "0", "1e-4", "3e-4")]
Z = [float(p["Z7"]) for p in pts]
res = {}
for ch in ("width_diphoton_GeV", "width_Zgamma_GeV"):
    # Gamma = |a + b Z7|^2 with complex a (light-fermion absorptive parts) and real b is exactly quadratic in Z7:
    # fit the quadratic on 3 points (Z7 = +2.6e-4, ~0, -3e-4), predict the 2 others.
    import numpy as np
    fit = [0, 2, 4]
    c = np.polyfit([Z[k] for k in fit], [pts[k][ch] for k in fit], 2)
    pred = [float(np.polyval(c, zz)) for zz in Z]
    res[ch] = dict(Z7=Z, width=[p[ch] for p in pts], quadratic_prediction=pred, fit_points=fit,
                   max_rel_dev_predicted=max(abs(pred[k] - pts[k][ch]) / pts[k][ch] for k in (1, 3)),
                   discriminant_4ac_minus_b2_over_b2=float((4 * c[0] * c[2] - c[1] ** 2) / c[1] ** 2))
audit["T2_lambda6_drives_Z7_and_photon_widths"] = res
audit["T2_verdict"] = ("PASS" if all(v["max_rel_dev_predicted"] < 1e-4 for v in res.values()) and len(set(Z)) == 5 else "FAIL")
audit["T2_widths_gg_bb_independent_of_lambda6_rel_spread"] = dict(
    gg=(max(p["width_gg_GeV"] for p in pts) - min(p["width_gg_GeV"] for p in pts)) / pts[2]["width_gg_GeV"],
    bb=(max(p["width_bb_GeV"] for p in pts) - min(p["width_bb_GeV"] for p in pts)) / pts[2]["width_bb_GeV"])

# 3) large-t stability of lambda7 and the ill-conditioning of lambda1(lambda7)
rows = []
for t in ("1000", "1000000", "3000000", "10000000", "30000000"):
    p = core.evaluate(200, t, "-1e-4", full=True)
    l7f = float(p["lambda7"])  # double rounding
    ic = core.identity_checks(200, t, "-1e-4", repr(l7f))
    rows.append(dict(tan_beta=t, lambda7=p["lambda7"], identity_rel_residual_50digit=p["lambda7_identity_rel_residual"],
                     lambda1_from_double_rounded_lambda7_abs_err=ic["lambda1_from_lambda7_abs_residual"],
                     dlambda1_dlambda7=float(t) ** 3 / 2,
                     prec50_vs_80_lambda7_rel=p["prec50_vs_80_lambda7_rel"],
                     backend_mass_max_abs_GeV=p["backend_mass_max_abs_GeV"],
                     backend_mass_over_t_eps=p["backend_mass_max_abs_GeV"] / (float(t) * 2.22e-16 * 200),
                     hp_mass_reconstruction_max_abs_GeV=p["hp_mass_reconstruction_max_abs_GeV"],
                     backend_kappa_top_rel=p["backend_kappa_top_vs_minus_cotb_rel"],
                     backend_CphiHpHm_rel=p["backend_CphiHpHm_vs_minus_vZ7_rel"],
                     route_F50_D80=p["routeF50_vs_D80_max_abs"], linear_route=p["linear_route_agreement_max_abs"],
                     tree_status=p["tree_status"], ctau_mm=p["ctau_mm"]))
audit["T3_large_tanb_stability"] = rows

# 4) loop-width cross-check against an independent LO formula (scheme sensitivity)
cmp_ = []
for m, t, l6 in ((200, "3000000", "0"), (200, "3000000", "-2.6364e-4"), (200, "10000000", "-3e-4"),
                 (250, "3000000", "-3e-4"), (300, "3000000", "-1.5346e-4"), (300, "10000000", "-3e-4")):
    p = core.evaluate(m, t, l6)
    tot, ff, hh = gamma_gaga_LO(m, float(t), float(p["Z7"]))
    cmp_.append(dict(mphi=m, tan_beta=t, lambda6=l6, backend=p["width_diphoton_GeV"], LO_independent=tot,
                     LO_fermion_only=ff, LO_Hpm_only=hh, ratio_backend_over_LO=p["width_diphoton_GeV"] / tot,
                     Hpm_fraction_of_LO_amplitude=math.sqrt(hh) / (math.sqrt(hh) + math.sqrt(ff)) if hh + ff else None))
audit["T4_diphoton_vs_independent_LO"] = cmp_

json.dump(audit, open(OUT, "w"), indent=1, default=str)
print(json.dumps({k: (v if k.endswith("verdict") else "...") for k, v in audit.items()}, indent=1))
print(json.dumps(audit["T1_force_lambda7_zero"], indent=1, default=str))
for x in rows: print(x)
for x in cmp_: print(x)
print(audit["T2_widths_gg_bb_independent_of_lambda6_rel_spread"], {k: v["max_rel_dev_predicted"] for k, v in res.items()})


sm = gamma_zga_LO(125.13, 1.0, with_W=True)
cmpz = dict(SM_h125_Zgamma_LO_GeV=sm, SM_reference_GeV="~6.3e-6 (LHCHXSWG BR 1.5e-3 x Gamma_h 4.1e-3)", points=[])
for m, t, l6 in ((200, "3000000", "-2.6364e-4"), (200, "10000000", "-3e-4"), (250, "3000000", "-3e-4"),
                 (300, "3000000", "-1.5346e-4"), (300, "10000000", "-3e-4"), (200, "3000000", "0")):
    p = core.evaluate(m, t, l6)
    lo = gamma_zga_LO(m, -1 / float(t), float(p["Z7"]))
    lo_os = gamma_zga_LO(m, -1 / float(t), float(p["Z7"]), sw2=1 - MW**2 / MZ**2)
    cmpz["points"].append(dict(mphi=m, tan_beta=t, lambda6=l6, backend=p["width_Zgamma_GeV"], LO_independent=lo,
                               ratio_backend_over_LO=p["width_Zgamma_GeV"] / lo, ratio_backend_over_LO_onshell_sw2=p["width_Zgamma_GeV"] / lo_os,
                               backend_ratio_Zga_over_gaga=p["width_Zgamma_GeV"] / p["width_diphoton_GeV"]))
audit["T5_Zgamma_vs_independent_LO"] = cmpz
json.dump(audit, open(OUT, "w"), indent=1, default=str)
print("SM h->Zgamma LO", sm)
for x in cmpz["points"]: print(x)
