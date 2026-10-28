"""Adversarial checks of the c tau control study.  Exit code != 0 on any failure.

  <venv>/bin/python -I scripts/adversarial.py     (after scan.py; writes raw/adversarial.json)
"""
import copy, json, math, sys
from decimal import Decimal, localcontext
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ctl
from ctl import core, F

H = ctl.HERE
fails, rep = [], {}


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
    rep[name] = dict(ok=bool(cond), detail=str(detail))
    if not cond:
        fails.append(name)


scan = [json.loads(l) for l in open(H / "raw/scan_points.jsonl")]
refs = ctl.references()
tgt = [r for r in scan if r["experiment"] == "A_target"]

# 1. no cache: a fresh, empty cache file; every point recomputed by the backend in this process
fresh = H / "raw/adversarial_fresh_cache.jsonl"
fresh.unlink(missing_ok=True)
core.CACHE, core._cache = fresh, None
worst = 0.0
for r in tgt:
    n = ctl.evaluate(r["X"], r["tan_beta"], r["family"])
    for k in ("Gamma_bb_GeV", "Gamma_gluglu_GeV", "Gamma_Zgamma_GeV", "Gamma_gammagamma_GeV", "Gamma_GeV"):
        worst = max(worst, abs(n[k] / r[k] - 1))
check("fresh-process, empty-cache recomputation reproduces every target point", worst == 0.0, f"max rel diff {worst:.1e}")
check("the fresh cache holds one distinct backend call per target point",
      len({json.loads(l)["key"] for l in open(fresh)}) == len({(r["X"], r["tan_beta"]) for r in tgt}))

# 2. widths come from the backend, not from a scaling formula: each partial width scales individually and
#    the photon widths agree with the independent LO loop formula (scripts/loops.py of 2026-10-19) at every t
for fam in refs:
    rows = [r for r in scan if r["family"] == fam and r["experiment"] == "A_grid" and r["status"] == "PASS"
            and 1e3 <= float(r["tan_beta"]) <= 1e10]
    rr = [r["ratio_gaga_backend_over_LO"] for r in rows]
    check(f"{fam}: backend Gamma(gamma gamma) / independent LO is constant along the trajectory",
          max(rr) - min(rr) < 1e-5, f"{min(rr):.7f}..{max(rr):.7f} over {len(rows)} points")
    rz = [r["ratio_Zga_backend_over_LO_os"] for r in rows]
    check(f"{fam}: backend Gamma(Z gamma) / independent LO is constant along the trajectory", max(rz) - min(rz) < 1e-5,
          f"{min(rz):.7f}..{max(rz):.7f}")
    ref = next(r for r in scan if r["family"] == fam and r["experiment"] == "reference")
    far = max(rows, key=lambda r: float(r["tan_beta"]))
    t0, t1 = float(ref["tan_beta"]), float(far["tan_beta"])
    # Z7 frozen at its reference value (i.e. a code that rescaled only the fermion couplings) would move B_gamma:
    z0, z1 = float(ref["Z7"]), float(far["Z7"])
    lo_frozen = F.gamma_gaga_LO(250, t1, z0)[0] / F.gamma_gaga_LO(250, t1, z1)[0]
    if float(ref["X"]) != 0:
        check(f"{fam}: the H+- loop matters (Z7 frozen at the reference would change Gamma(gamma gamma) by x{lo_frozen:.3g})",
              abs(lo_frozen - 1) > 0.5)
    check(f"{fam}: t*Z7 + X equals kappa = 2d - (lambda1 - a) along the trajectory",
          max(abs(float(r["Z7"]) * float(r["tan_beta"]) + float(r["X"]) - kappa) for r in rows
              for kappa in [2 * (250**2 - 125.13**2 / 2) / 246**2 - (1.2 - 125.13**2 / 246**2)]) < 1e-5 * max(1, abs(float(ref["X"]))))

# 3. negative controls on a far trajectory point (mixed, c tau = 1e4 mm)
far = next(r for r in tgt if r["family"] == "mixed" and r["target_ctau_mm"] == 1e4)
full = F.full_eval(250, far["tan_beta"], far["lambda6"])
full["accounting"] = core.accounting(full)
check("far mixed point passes every validity gate", all(v[0] for k, v in F.gates(full, "mixed").items() if k not in ctl.NOT_VALIDITY))
bad = copy.deepcopy(full); bad.update(core.identity_checks(250, far["tan_beta"], far["lambda6"], repr(float(full["lambda7"]))))
check("double-rounded lambda7 at t=6.38e7 is caught (lambda1 reconstruction)",
      not F.gates(bad, "mixed")["lambda1_reconstructed_from_lambda7_50digit"][0], f"lambda1 off by {bad['lambda1_from_lambda7_abs_residual']:.2e}")
bad = copy.deepcopy(full); bad["lambda7"] = "0"; bad.update(core.identity_checks(250, far["tan_beta"], far["lambda6"], "1e-300"))
check("lambda7 = 0 is rejected", not F.gates(bad, "mixed")["lambda7_genuinely_nonzero"][0])
l6r = f"{float(far['lambda6']):.3e}"  # a lambda6 rounded to 4 digits silently changes X
r = ctl.evaluate(far["X"], far["tan_beta"], "mixed")
with localcontext() as c:
    c.prec = 60
    dx = abs(Decimal(l6r) * Decimal(far["tan_beta"]) - Decimal(far["X"]))
check("a rounded lambda6 would change X (the X gate is sensitive)", dx > Decimal("1e-40") and r["status"] == "PASS", f"dX={dx:.2e}")
# simplified checks are insufficient: a point with an exact BFB violation passes the 4 pi screen
weak = next(r for r in scan if r["family"] == "mixed" and r["experiment"] == "A_grid" and r["tree_BFB_general"] == "FAIL"
            and r["tree_quartic_screen_4pi"] == "PASS" and r["tree_full_tree_unitarity"] == "PASS")
check("a BFB-violating point passes the 4 pi screen and unitarity (so these alone would accept it) and is FAIL here",
      weak["status"] == "FAIL", f"t={weak['tan_beta']} lambda6={float(weak['lambda6']):.3g}")
# composition responds to X (BRs are not frozen by construction)
b2 = {r["note"]: r for r in scan if r["experiment"] == "B2" and r["family"] == "mixed"}
check("B_gamma responds to X at fixed t (X/X0 = 0.5 vs 1.1)", abs(b2["X/X0=0.5"]["B_gamma"] - b2["X/X0=1.1"]["B_gamma"]) > 0.1)

# 4. every scan row with status PASS kept X, Y2 = 0, lambda7 != 0 and BR normalisation
ok = [r for r in scan if r["status"] == "PASS"]
check(f"all {len(ok)} PASS rows: Y2 = 0, lambda7 != 0, |sum BR - 1| < 1e-12, ctau = hbar c / Gamma",
      all(Decimal(r["Y2_GeV2"]) == 0 and Decimal(r["lambda7"]) != 0 and abs(r["BR_inventory_sum_minus_1"]) < 1e-12
          and abs(r["ctau_mm"] * r["Gamma_GeV"] / core.HBARC - 1) < 1e-14 for r in ok))
with localcontext() as c:
    c.prec = 60
    dg = max(abs(Decimal(r["g_hphiphi_GeV"]) - 2 * Decimal(250) ** 2 / Decimal(246)) for r in ok)
check("all PASS rows have g_hphiphi = 2 m^2 / v to 1e-30 (production coupling unchanged along the trajectory)",
      dg < Decimal("1e-30"), f"max dev {dg:.1e} GeV")

json.dump(rep, open(H / "raw/adversarial.json", "w"), indent=1)
print(f"\n{len(fails)} failure(s)")
sys.exit(1 if fails else 0)
