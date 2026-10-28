"""Validation gates, cell selection and all deliverable tables from ONE structured dataset.

  <venv>/bin/python -I scripts/finalize.py      (after search.py and production.py)

Writes candidate_inventory.csv, validated_benchmarks.csv, rejected_and_inconclusive.csv, table_inputs.csv,
table_A.csv, table_B.csv, status_matrix.csv, raw/final_points.json and tables/*.tex.
"""
import os
import csv, json, math, sys
from decimal import Decimal
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import core
from loops import gamma_gaga_LO, gamma_zga_LO, MW, MZ
from search import TARGETS, dist

H = core.HERE
OLD = core.ROOT / "reports/lambda7_tables_20261014/BENCHMARKS.csv"
LUMI = 300.0  # fb^-1
ZBR = {"jj": 0.69911, "ee": 0.033632, "mumu": 0.033662}  # PDG 2024
MASSES, FAMILIES = (200, 250, 300), ("hadronic", "mixed", "photonic")
EPS = 2.220446049250313e-16
# "gluglu" = gluon pair (never a bare "gg": legacy files used gg for gamma gamma)
DIRECT = ["bb", "gluglu", "Zgamma", "gammagamma"]
EXP = ["bb", "gluglu", "jjgamma", "eegamma", "mumugamma", "gammagamma"]
LABEL = {"bb": r"b\bar b", "gluglu": "gg", "Zgamma": r"Z\gamma", "gammagamma": r"\gamma\gamma", "jjgamma": r"jj\gamma",
         "eegamma": r"ee\gamma", "mumugamma": r"\mu\mu\gamma"}


def pairs(modes):
    return [(a, b) for i, a in enumerate(modes) for b in modes[i:]]


def yields(sigma_fb, br):
    """N_ab = L sigma BR_a^2 (a=b) or 2 L sigma BR_a BR_b (a!=b), unordered."""
    return {f"({a})({b})": LUMI * sigma_fb * (br[a] ** 2 if a == b else 2 * br[a] * br[b]) for a, b in pairs(list(br))}


def expanded_br(br4):
    return dict(bb=br4["bb"], gluglu=br4["gluglu"], jjgamma=br4["Zgamma"] * ZBR["jj"], eegamma=br4["Zgamma"] * ZBR["ee"],
                mumugamma=br4["Zgamma"] * ZBR["mumu"], gammagamma=br4["gammagamma"])


def gates(r, family=None):
    """Hard validation gates on a fully evaluated record (evaluate(full=True) + accounting). Returns {gate: (ok, detail)}."""
    m, t = r["mphi_GeV"], float(r["tan_beta"])
    l7 = Decimal(r["lambda7"])
    g = {}
    g["tree_tests_8_PASS_50digit"] = (r["tree_status"] == "PASS_ALL_TESTED", r["tree_status"])
    g["tree_tests_8_PASS_80digit"] = (r.get("prec80_tree_status") == "PASS_ALL_TESTED", r.get("prec80_tree_status"))
    g["Y2_exact_zero"] = (bool(r["Y2_is_exact_zero"]) and Decimal(r["Y2_GeV2"]) == 0, r["Y2_GeV2"])
    g["M2_eq_mh2_over_2"] = (abs(Decimal(r["M2_minus_mh2_over_2"])) < Decimal("1e-30"), r["M2_minus_mh2_over_2"])
    g["g_eq_2m2_over_v"] = (abs(Decimal(r["g_minus_2m2_over_v"])) < Decimal("1e-30"), r["g_minus_2m2_over_v"])
    lead = Decimal(r["lambda7_leading_term"])
    g["lambda7_genuinely_nonzero"] = (l7 != 0 and r["lambda7_identity_rel_residual"] < 1e-30 and Decimal("0.5") < l7 / lead < 2,
                                      f"lambda7={r['lambda7'][:14]} lead={r['lambda7_leading_term'][:12]} rel_id={r['lambda7_identity_rel_residual']:.1e}")
    g["lambda1_reconstructed_from_lambda7_50digit"] = (r["lambda1_from_lambda7_abs_residual"] < 1e-20, r["lambda1_from_lambda7_abs_residual"])
    g["constructor_routes_agree"] = (r["routeF50_vs_D80_max_abs"] < 1e-40 and r["linear_route_agreement_max_abs"] < 1e-40
                                     and r["closure_core_pass"], f"{r['routeF50_vs_D80_max_abs']:.1e}/{r['linear_route_agreement_max_abs']:.1e}")
    g["trilinear_two_routes"] = (r["g_closure_minus_formula_GeV"] < 1e-30 and r["backend_g_hphiphi_vs_constructor_rel"] < 1e-12,
                                 f"closure {r['g_closure_minus_formula_GeV']:.1e}; backend {r['backend_g_hphiphi_vs_constructor_rel']:.1e}")
    g["precision_50_vs_80"] = (r["prec50_vs_80_lambda7_rel"] < 1e-40 and r["prec50_vs_80_Z7_rel"] < 1e-40, r["prec50_vs_80_lambda7_rel"])
    g["mass_reconstruction_high_precision"] = (r["hp_mass_reconstruction_max_abs_GeV"] < 1e-30, r["hp_mass_reconstruction_max_abs_GeV"])
    g["backend_double_mass_within_roundoff"] = (r["backend_mass_max_abs_GeV"] < max(t, 1) * EPS * (m + 100),
                                                f"{r['backend_mass_max_abs_GeV']:.2e} GeV vs bound {max(t, 1) * EPS * (m + 100):.2e}")
    g["alignment_backend"] = (abs(r["backend_sba"] - 1) < 1e-12, r["backend_sba"])
    g["charged_loop_coupling_CV_eq_minus_vZ7"] = (r["backend_CphiHpHm_vs_minus_vZ7_rel"] < 1e-9, r["backend_CphiHpHm_vs_minus_vZ7_rel"])
    g["kappa_top_eq_minus_cotb"] = (r["backend_kappa_top_vs_minus_cotb_rel"] < 10 * t * EPS, r["backend_kappa_top_vs_minus_cotb_rel"])
    acc = r["accounting"]
    s = sum(acc["BR"].values()) + acc["BR_other_terms"]
    g["complete_inventory_sums_to_1"] = (abs(s - 1) < 1e-12 and acc["total_vs_bridge_rel"] < 1e-12, f"sum-1={s - 1:.1e}")
    g["diphoton_vs_independent_LO"] = (abs(r["ratio_gaga_backend_over_LO"] - 1) < (0.02 if r["B_gamma"] > 0.5 else 0.10), r["ratio_gaga_backend_over_LO"])
    g["Zgamma_vs_independent_LO_onshell_sw"] = (abs(r["ratio_Zga_backend_over_LO_os"] - 1) < (0.02 if r["B_gamma"] > 0.5 else 0.10), r["ratio_Zga_backend_over_LO_os"])
    g["ctau_80_120_mm"] = (80 <= r["ctau_mm"] <= 120, r["ctau_mm"])
    g["ctau_eq_hbarc_over_Gamma"] = (abs(r["ctau_mm"] * r["Gamma_GeV"] / core.HBARC - 1) < 1e-14, "")
    if family:
        g["composition_class"] = (r["family_class"] == family, r["family_class"])
    return g


def full_eval(m, t, l6):
    r = dict(core.evaluate(m, t, l6, full=True))
    r["accounting"] = core.accounting(r)
    zr = float(r["Z7"])
    r["ratio_gaga_backend_over_LO"] = r["width_diphoton_GeV"] / gamma_gaga_LO(m, float(t), zr)[0]
    r["ratio_Zga_backend_over_LO_os"] = r["width_Zgamma_GeV"] / gamma_zga_LO(m, -1 / float(t), zr, sw2=1 - MW**2 / MZ**2)
    r["ratio_Zga_backend_over_LO_msbar"] = r["width_Zgamma_GeV"] / gamma_zga_LO(m, -1 / float(t), zr)
    return r


def main():
    prod = {int(k): v for k, v in json.load(open(H / "raw/production.json")).items()}
    sel = json.load(open(H / "raw/search_selected.json"))
    cands = []  # (cid, source, m, t, l6, origin)
    for row in csv.DictReader(open(OLD)):
        src = f"{OLD}#id={row['id']} (consolidation point {row['source_id']})"
        if row["family"] == "A":
            cands.append((row["id"], src, int(row["mphi_GeV"]), row["tan_beta"], row["lambda6"], "recovered_October_A"))
    for k, s in sel.items():
        parts = k.split("_")
        fam, m = parts[0], int(parts[1])
        sign = parts[2] if len(parts) > 2 else "0"
        cid = f"{fam[0].upper()}{m}{'' if sign == '0' else ('n' if sign == '-' else 'p')}"
        cands.append((cid, f"raw/search_selected.json#{k}; raw/scan_ledger.csv (phase final_rounded/hadronic_solve)", m,
                      s["tan_beta"], s["lambda6"], "new_search_2026-10-19"))
    recs, seen = {}, {}
    for cid, src, m, t, l6, origin in cands:
        key = (m, Decimal(t), Decimal(l6))
        if key in seen:  # identical input point found twice (e.g. A3 re-found by the hadronic solve)
            recs[seen[key]]["also_found_as"].append(cid)
            continue
        r = full_eval(m, t, l6)
        r.update(cid=cid, source=src, origin=origin, also_found_as=[])
        r["gates"] = gates(r)
        r["dist_hist"] = {f: dist(r, TARGETS[f][m]) for f in FAMILIES if m in TARGETS[f]}
        recs[cid] = r; seen[key] = cid

    # ---------- cell selection: same class, all gates, closest to the historical mass-specific composition ----------
    final, matrix, reasons = {}, {}, {}
    for f in FAMILIES:
        for m in MASSES:
            pool = [r for r in recs.values() if r["mphi_GeV"] == m and r["family_class"] == f]
            ok = [r for r in pool if all(v[0] for v in gates(r, f).values())]
            if ok:
                best = min(ok, key=lambda r: (r["dist_hist"].get(f, 0.0), r["cid"]))
                final[(f, m)] = best; matrix[(f, m)] = "PASS"
            elif pool:
                matrix[(f, m)] = "INCONCLUSIVE"
                reasons[(f, m)] = "; ".join(f"{r['cid']}: " + ",".join(k for k, v in gates(r, f).items() if not v[0]) for r in pool)
            else:
                matrix[(f, m)] = "MISSING"
                reasons[(f, m)] = "no candidate of this class found in the searched range"

    # ---------- independent consistency assertions ----------
    for (f, m), r in final.items():
        assert r["mphi_GeV"] == m and Decimal(r["Y2_GeV2"]) == 0 and Decimal(r["lambda7"]) != 0
        assert abs(prod[m]["g_hphiphi_GeV"] - float(r["g_hphiphi_GeV"])) < 1e-9, "sigma kernel coupling != point coupling"
        assert len({rr["family_class"] for (ff, mm), rr in final.items() if ff == f}) == 1

    # ---------- per-point rejection / inconclusive ledger ----------
    rej = []
    for r in recs.values():
        fails = [k for k, v in r["gates"].items() if not v[0]]
        chosen = [f"{f}@{m}" for (f, m), x in final.items() if x is r]
        if chosen:
            continue
        why = []
        if r["ctau_mm"] < 80:
            why.append(f"prompt/not LLP (ctau={r['ctau_mm']:.3g} mm)")
        if r["family_class"] == "transitional":
            why.append(f"transitional composition (B_h={r['B_h']:.3f}, B_gamma={r['B_gamma']:.3f}) - matches no family definition")
        fam_m = (r["family_class"], int(r["mphi_GeV"]))
        if fam_m in final and final[fam_m] is not r and not fails:
            why.append(f"valid {r['family_class']} point but {final[fam_m]['cid']} is closer to the historical composition "
                       f"(d={r['dist_hist'].get(r['family_class'], float('nan')):.4f} vs {final[fam_m]['dist_hist'].get(r['family_class'], float('nan')):.4f})")
        if fails:
            why.append("failed gates: " + ",".join(fails))
        rej.append(dict(candidate_id=r["cid"], mphi_GeV=int(r["mphi_GeV"]), tan_beta=r["tan_beta"], lambda6=r["lambda6"],
                        family_class=r["family_class"], ctau_mm=f"{r['ctau_mm']:.6g}", B_h=f"{r['B_h']:.4f}", B_gamma=f"{r['B_gamma']:.4f}",
                        verdict="REJECTED" if (r["ctau_mm"] < 80 or r["family_class"] == "transitional") else "NOT_SELECTED",
                        reasons=" | ".join(why) or "not selected", source=r["source"]))
    for row in csv.DictReader(open(OLD)):
        if row["family"] == "B":
            rej.append(dict(candidate_id=row["id"], mphi_GeV=int(row["mphi_GeV"]), tan_beta=row["tan_beta"], lambda6=row["lambda6"],
                            family_class="(not evaluated)", ctau_mm="", B_h="", B_gamma="", verdict="EXCLUDED_BY_CONSTRUCTION",
                            reasons="construction B (free M^2 / lambda7=0 control): Y2!=0 or lambda7~0; excluded from final tables by rule",
                            source=f"{OLD}#id={row['id']}"))
    for (f, m), why in reasons.items():
        rej.append(dict(candidate_id=f"CELL_{f}_{m}", mphi_GeV=m, tan_beta="", lambda6="", family_class=f, ctau_mm="", B_h="", B_gamma="",
                        verdict=matrix[(f, m)], reasons=why, source=""))
    write_csv(H / "rejected_and_inconclusive.csv", rej)

    # ---------- inventory: every evaluated candidate + coarse scan rows ----------
    inv = []
    for r in recs.values():
        inv.append(dict(candidate_id=r["cid"], also_found_as=";".join(r["also_found_as"]), origin=r["origin"], mphi_GeV=int(r["mphi_GeV"]),
                        tan_beta=r["tan_beta"], lambda6=r["lambda6"], X=r["X_lambda6_tanbeta"], lambda7=r["lambda7"], Z7=r["Z7"],
                        ctau_mm=f"{r['ctau_mm']:.8g}", BR_bb=f"{r['BR_bb']:.6e}", BR_gg=f"{r['BR_gg']:.6e}", BR_Zgamma=f"{r['BR_Zgamma']:.6e}",
                        BR_gammagamma=f"{r['BR_diphoton']:.6e}", B_h=f"{r['B_h']:.6f}", B_gamma=f"{r['B_gamma']:.6f}", family_class=r["family_class"],
                        gates_failed=",".join(k for k, v in r["gates"].items() if not v[0]) or "none",
                        selected_for=";".join(f"{f}@{m}" for (f, m), x in final.items() if x is r), source=r["source"]))
    for row in csv.DictReader(open(H / "raw/scan_ledger.csv")):
        if row["phase"] == "coarse":
            inv.append(dict(candidate_id=f"scan_{row['mphi_GeV']}_{row['sign']}_{row['tan_beta']}", also_found_as="", origin="coarse_scan_2026-10-19",
                            mphi_GeV=row["mphi_GeV"], tan_beta=row["tan_beta"], lambda6=row["lambda6"], X=row["X"], lambda7=row["lambda7"], Z7=row["Z7"],
                            ctau_mm=row["ctau_mm"], BR_bb=row["BR_bb"], BR_gg=row["BR_gg"], BR_Zgamma=row["BR_Zgamma"], BR_gammagamma=row["BR_gammagamma"],
                            B_h=row["B_h"], B_gamma=row["B_gamma"], family_class=row["family_class"],
                            gates_failed="(tree tests only: " + row["tree_status"] + ")", selected_for="", source="raw/scan_ledger.csv"))
    write_csv(H / "candidate_inventory.csv", inv)

    # ---------- validated benchmarks (full precision) ----------
    vb = []
    for f in FAMILIES:
        for m in MASSES:
            if (f, m) not in final:
                continue
            r = final[(f, m)]; acc = r["accounting"]["BR"]; p = prod[m]
            vb.append(dict(family=f, mphi_GeV=m, candidate_id=r["cid"], also_found_as=";".join(r["also_found_as"]), status=matrix[(f, m)],
                           mA_GeV=m + 100, mHp_GeV=m + 100, mh_GeV=core.MH, v_GeV=core.V, sin_beta_minus_alpha=1,
                           tan_beta=r["tan_beta"], **{f"lambda{k}": r[f"lambda{k}"] for k in range(1, 8)}, lambda6_input=r["lambda6"],
                           X=r["X_lambda6_tanbeta"], Z7=r["Z7"], M2_GeV2=r["M_sq_GeV2"], m12_sq_GeV2=r["m12_sq_GeV2"], Y2_GeV2=r["Y2_GeV2"],
                           g_hphiphi_GeV=r["g_hphiphi_GeV"], C_phiHpHm_GeV=r["CphiHpHm_potential_GeV"],
                           Gamma_phi_GeV=repr(r["Gamma_GeV"]), ctau_mm=repr(r["ctau_mm"]),
                           **{f"BR_{k}": repr(v) for k, v in acc.items()}, BR_other_terms=repr(r["accounting"]["BR_other_terms"]),
                           B_h=repr(r["B_h"]), B_gamma=repr(r["B_gamma"]), dist_to_historical=repr(r["dist_hist"].get(f)),
                           sigma_LO_HEFT_fb=repr(p["sigma_LO_HEFT_fb"]), sigma_mc_err_fb=repr(p["mc_err_fb"]),
                           N_pairs_produced=repr(LUMI * p["sigma_LO_HEFT_fb"]),
                           ratio_gaga_backend_over_LO=repr(r["ratio_gaga_backend_over_LO"]),
                           ratio_Zga_backend_over_LO_onshell=repr(r["ratio_Zga_backend_over_LO_os"]),
                           ratio_Zga_backend_over_LO_msbar=repr(r["ratio_Zga_backend_over_LO_msbar"]),
                           backend_double_mass_residual_GeV=repr(r["backend_mass_max_abs_GeV"]),
                           experimental_constraints="NOT_TESTED", detector_efficiency="NOT_APPLIED (truth level)",
                           source=r["source"], backend=core.C.BACKEND))
    write_csv(H / "validated_benchmarks.csv", vb)
    json.dump({f"{f}_{m}": {k: v for k, v in r.items()} for (f, m), r in final.items()} | {"_all": {k: v for k, v in recs.items()}},
              open(H / "raw/final_points.json", "w"), indent=1, default=str)

    # ---------- tables ----------
    ti, ta, tb = [], [], []
    for f in FAMILIES:
        col = {m: final.get((f, m)) for m in MASSES}
        def cell(fn, m):
            return fn(col[m], m) if col[m] is not None else "MISSING"
        BR = lambda r, k: r["accounting"]["BR"][k]
        covA = lambda r: (BR(r, "bb") + BR(r, "gluongluon") + BR(r, "Zgamma") + BR(r, "gammagamma")) ** 2
        covB = lambda r: (BR(r, "bb") + BR(r, "gluongluon") + BR(r, "Zgamma") * sum(ZBR.values()) + BR(r, "gammagamma")) ** 2
        rows = [  # (role, quantity, formatter); inputs printed exactly as used
            ("input", "m_phi [GeV]", lambda r, m: f"{m}"),
            ("input", "tan(beta)", lambda r, m: f"{float(r['tan_beta']):.2e}"),
            ("input", "lambda6", lambda r, m: f"{float(r['lambda6']):.5e}"),
            ("fixed", "m_A = m_H+- [GeV]", lambda r, m: f"{m + 100}"), ("fixed", "m_h [GeV]", lambda r, m: core.MH),
            ("fixed", "v [GeV]", lambda r, m: core.V),
            ("fixed", "sin(beta-alpha)", lambda r, m: "1 (exact)"), ("fixed", "lambda1", lambda r, m: f"{float(r['lambda1']):.4g}"),
            ("fixed", "Y2 [GeV^2]", lambda r, m: "0 (exact)" if Decimal(r["Y2_GeV2"]) == 0 else r["Y2_GeV2"]),
            ("derived", "lambda7", lambda r, m: f"{float(r['lambda7']):.6e}"),
            ("derived", "X = lambda6 tan(beta)", lambda r, m: f"{float(r['X_lambda6_tanbeta']):.6g}"),
            ("derived", "Z7", lambda r, m: f"{float(r['Z7']):.6e}"), ("derived", "M^2 [GeV^2]", lambda r, m: f"{float(r['M_sq_GeV2']):.6f}"),
            ("derived", "g_hphiphi [GeV]", lambda r, m: f"{float(r['g_hphiphi_GeV']):.2f}"),
            ("output", "Gamma_phi [GeV]", lambda r, m: f"{r['Gamma_GeV']:.5e}"),
            ("output", "ctau_phi [mm]", lambda r, m: f"{r['ctau_mm']:.2f}"),
            ("output", "BR(bb)", lambda r, m: f"{BR(r, 'bb'):.4g}"), ("output", "BR(gluon gluon)", lambda r, m: f"{BR(r, 'gluongluon'):.4g}"),
            ("output", "BR(Zgamma)", lambda r, m: f"{BR(r, 'Zgamma'):.4g}"), ("output", "BR(gammagamma)", lambda r, m: f"{BR(r, 'gammagamma'):.4g}"),
            ("output", "BR(tautau)", lambda r, m: f"{BR(r, 'tautau'):.3g}"), ("output", "BR(cc)", lambda r, m: f"{BR(r, 'cc'):.3g}"),
            ("output", "BR(ss)", lambda r, m: f"{BR(r, 'ss'):.2g}"), ("output", "BR(mumu)", lambda r, m: f"{BR(r, 'mumu'):.2g}"),
            ("output", "BR(t t*) off-shell", lambda r, m: f"{BR(r, 'ttstar'):.3g}"),
            ("output", "BR(all other terms)", lambda r, m: f"{r['accounting']['BR_other_terms'] + sum(v for k, v in r['accounting']['BR'].items() if k in ('dd', 'uu', 'ee')):.1e}"),
            ("output", "sigma_LO_HEFT [fb]", lambda r, m: f"{prod[m]['sigma_LO_HEFT_fb']:.4f}"),
            ("output", "sigma MC error [fb]", lambda r, m: f"{prod[m]['mc_err_fb']:.4f}"),
            ("output", "phi pairs produced (300 fb^-1)", lambda r, m: f"{LUMI * prod[m]['sigma_LO_HEFT_fb']:.1f}"),
            ("output", "Table A coverage (fraction of phi-phi pairs)", lambda r, m: f"{covA(r):.4f}"),
            ("output", "Table B coverage (fraction of phi-phi pairs)", lambda r, m: f"{covB(r):.4f}"),
            ("output", "theoretical status", lambda r, m: f"{matrix[(f, m)]} (tree level; collider limits NOT_TESTED)"),
            ("output", "candidate ID / provenance", lambda r, m: f"{r['cid']}" + (f" (= {','.join(r['also_found_as'])})" if r["also_found_as"] else "")),
        ]
        NOTE = {"lambda7": "derived; do NOT re-input: dlambda1/dlambda7 = tan^3(beta)/2, use the 50-digit value in validated_benchmarks.csv",
                "Z7": "derived (Higgs basis); 50 digits in validated_benchmarks.csv", "X = lambda6 tan(beta)": "derived",
                "Table A coverage (fraction of phi-phi pairs)": "(BR_bb+BR_gluglu+BR_Zgamma+BR_gammagamma)^2 = sum(Table A)/pairs produced",
                "Table B coverage (fraction of phi-phi pairs)": "(BR_bb+BR_gluglu+BR_Zgamma*(jj+ee+mumu)+BR_gammagamma)^2 = sum(Table B)/pairs produced",
                "Y2 [GeV^2]": "branch condition M^2 = m_h^2/2", "sigma MC error [fb]": "statistical, inverse-variance mean of 2 seeds"}
        for role, q, fn in rows:
            ti.append(dict(family=f, role=role, quantity=q, **{f"{m} GeV": cell(fn, m) for m in MASSES}, note=NOTE.get(q, "")))
        for m in MASSES:
            r = col[m]
            if r is None:
                continue
            a = r["accounting"]["BR"]
            br4 = dict(bb=a["bb"], gluglu=a["gluongluon"], Zgamma=a["Zgamma"], gammagamma=a["gammagamma"])
            r["_yA"] = yields(prod[m]["sigma_LO_HEFT_fb"], br4)
            r["_yB"] = yields(prod[m]["sigma_LO_HEFT_fb"], expanded_br(br4))
        for key in yields(1, {k: 0 for k in DIRECT}):
            ta.append(dict(family=f, channel=key, **{f"{m} GeV": (repr(col[m]["_yA"][key]) if col[m] else "MISSING") for m in MASSES}))
        for key in yields(1, {k: 0 for k in EXP}):
            tb.append(dict(family=f, channel=key, **{f"{m} GeV": (repr(col[m]["_yB"][key]) if col[m] else "MISSING") for m in MASSES}))
    assert len(ta) == 30 and len(tb) == 63
    write_csv(H / "table_inputs.csv", ti); write_csv(H / "table_A.csv", ta); write_csv(H / "table_B.csv", tb)
    write_csv(H / "status_matrix.csv", [dict(family=f, **{f"{m} GeV": matrix[(f, m)] + (f" ({final[(f, m)]['cid']})" if (f, m) in final else "")
                                                         for m in MASSES}) for f in FAMILIES])
    latex(ti, ta, tb, matrix, final)
    for f in FAMILIES:
        print(f, {m: (matrix[(f, m)], final[(f, m)]["cid"] if (f, m) in final else None) for m in MASSES})
    for k, why in reasons.items():
        print("NOT FILLED", k, why[:300])


def write_csv(p, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)


def tex_num(x):
    if x in ("MISSING",):
        return r"\textit{missing}"
    v = float(x)
    if v == 0:
        return "0"
    if v >= 100:
        return f"{v:.0f}"
    if v >= 0.01:
        return f"{v:#.3g}" if v < 10 else f"{v:.1f}"
    e = int(math.floor(math.log10(v)))
    return rf"${v / 10**e:.2f}\times10^{{{e}}}$"


def tex_q(s):
    rep = {"m_phi [GeV]": r"$m_\phi$ [GeV]", "m_A = m_H+- [GeV]": r"$m_A=m_{H^\pm}$ [GeV]", "m_h [GeV]": r"$m_h$ [GeV]",
           "tan(beta)": r"$\tan\beta$", "sin(beta-alpha)": r"$\sin(\beta-\alpha)$", "lambda1": r"$\lambda_1$", "lambda6": r"$\lambda_6$",
           "lambda7": r"$\lambda_7$", "X = lambda6 tan(beta)": r"$X=\lambda_6\tan\beta$", "Z7": r"$Z_7$", "M^2 [GeV^2]": r"$M^2$ [GeV$^2$]",
           "Y2 [GeV^2]": r"$Y_2$ [GeV$^2$]", "g_hphiphi [GeV]": r"$g_{h\phi\phi}$ [GeV]", "Gamma_phi [GeV]": r"$\Gamma_\phi$ [GeV]",
           "ctau_phi [mm]": r"$c\tau_\phi$ [mm]", "BR(bb)": r"BR($b\bar b$)", "BR(gluon gluon)": r"BR($gg$, gluons)", "BR(Zgamma)": r"BR($Z\gamma$)",
           "BR(gammagamma)": r"BR($\gamma\gamma$)", "BR(tautau)": r"BR($\tau\tau$)", "BR(cc)": r"BR($c\bar c$)", "BR(ss)": r"BR($s\bar s$)",
           "BR(mumu)": r"BR($\mu\mu$)", "BR(t t*) off-shell": r"BR($t\bar t^{*}$) (off-shell)", "BR(all other terms)": r"BR(all other)",
           "sigma_LO_HEFT [fb]": r"$\sigma_{\rm LO\text{-}HEFT}$ [fb]", "phi pairs produced (300 fb^-1)": r"$\phi\phi$ pairs, 300 fb$^{-1}$", "v [GeV]": r"$v$ [GeV]",
           "sigma MC error [fb]": r"$\delta\sigma_{\rm MC}$ [fb]",
           "Table A coverage (fraction of phi-phi pairs)": r"Table A coverage of $\phi\phi$ pairs",
           "Table B coverage (fraction of phi-phi pairs)": r"Table B coverage of $\phi\phi$ pairs",
           "theoretical status": "theory status", "candidate ID / provenance": "candidate ID"}
    return rep.get(s, s)


def tex_val(s):
    s = str(s)
    if s == "MISSING":
        return r"\textit{missing}"
    if "e" in s and s.replace("e", "").replace("-", "").replace("+", "").replace(".", "").isdigit():
        mant, ex = s.split("e")
        if float(mant) == 0:
            return "0"
        return rf"${mant}\times10^{{{int(ex)}}}$"
    return s.replace(" (tree level; collider limits NOT_TESTED)", "$^{\\dagger}$").replace("_", r"\_")


ROLE_TEX = {"input": "Inputs (exact decimals)", "fixed": "Fixed by convention / branch condition",
            "derived": "Derived (do not re-input; 50-digit values in validated\\_benchmarks.csv)", "output": "Outputs"}


def sum_row(tab, ti, f):
    """Sum of the tabulated channels and its share of all produced phi-phi pairs."""
    npairs = {m: float(next(r for r in ti if r["family"] == f and r["quantity"].startswith("phi pairs"))[f"{m} GeV"]) for m in MASSES}
    tot = {m: sum(float(r[f"{m} GeV"]) for r in tab if r["family"] == f) for m in MASSES}
    return "Sum (share of all $\\phi\\phi$ pairs) & " + " & ".join(
        rf"{tex_num(repr(tot[m]))} ({100 * tot[m] / npairs[m]:.1f}\%)" for m in MASSES) + r"\\"


def chan_tex(k):
    a, b = k[1:-1].split(")(")
    return rf"$({LABEL[a]})({LABEL[b]})$"


def latex(ti, ta, tb, matrix, final):
    out = H / "tables"
    head = r"\begin{tabular}{lrrr}\toprule" + "\n" + r"Quantity & 200 GeV & 250 GeV & 300 GeV\\\midrule" + "\n"
    for f in FAMILIES:
        with open(out / f"inputs_{f}.tex", "w") as fh:
            lines, last = [], None
            for r in (r for r in ti if r["family"] == f):
                if r["role"] != last:
                    lines.append((r"\midrule" if last else "") + rf"\multicolumn{{4}}{{l}}{{\textit{{{ROLE_TEX[r['role']]}}}}}\\")
                    last = r["role"]
                q = tex_q(r["quantity"]) + (r"$^{\ddagger}$" if r["quantity"] == "lambda7" else "")
                lines.append(f"{q} & " + " & ".join(tex_val(r[f'{m} GeV']) for m in MASSES) + r"\\")
            fh.write(head + "\n".join(lines) + "\n" + r"\bottomrule\end{tabular}" + "\n")
        for name, tab in (("A", ta), ("B", tb)):
            with open(out / f"table{name}_{f}.tex", "w") as fh:
                fh.write(r"\begin{tabular}{lrrr}\toprule" + "\nChannel & 200 GeV & 250 GeV & 300 GeV\\\\\\midrule\n" +
                         "\n".join(chan_tex(r["channel"]) + " & " + " & ".join(tex_num(r[f"{m} GeV"]) for m in MASSES) + r"\\"
                                   for r in tab if r["family"] == f) + "\n" + r"\midrule" + "\n" + sum_row(tab, ti, f) + "\n" + r"\bottomrule\end{tabular}" + "\n")
    with open(out / "status_matrix.tex", "w") as fh:
        fh.write(r"\begin{tabular}{lccc}\toprule Family & 200 GeV & 250 GeV & 300 GeV\\\midrule" + "\n" +
                 "\n".join(f"{f.capitalize()} & " + " & ".join(matrix[(f, m)] + (f" ({final[(f, m)]['cid']})" if (f, m) in final else "")
                                                              for m in MASSES) + r"\\" for f in FAMILIES) + "\n" + r"\bottomrule\end{tabular}" + "\n")


if __name__ == "__main__":
    main()
