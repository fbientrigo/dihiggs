"""Regression / adversarial tests of the validation system and of the published tables.

  <venv>/bin/python -I scripts/regression_tests.py      (after finalize.py; exits non-zero on any failure)

1. deliberately wrong lambda7 (zeroed / double-rounded)       -> gates must fail
2. deliberately wrong Y2 (free-M^2 construction, rho != 0)     -> gates must fail
3. legacy UFO normalisation (+8 i m^2/v, i.e. GHphiphi=+4 g)   -> production gate must fail
4. sigma from a different mass                                  -> kernel/point coupling check must fail
5. tables recomputed independently from validated_benchmarks.csv (numpy outer product), combinatorics, Z expansion
"""
import copy, csv, json, sys, textwrap
from decimal import Decimal
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import core
import finalize as F

H = core.HERE
fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


pts = json.load(open(H / "raw/final_points.json"))
good_key = next(k for k in pts if not k.startswith("_"))
good = pts[good_key]
fam = good_key.split("_")[0]
check(f"reference point {good_key} passes all gates", all(v[0] for v in F.gates(good, fam).values()))

# 1) wrong lambda7
bad = copy.deepcopy(good); bad["lambda7"] = "0"
bad.update(core.identity_checks(int(good["mphi_GeV"]), good["tan_beta"], good["lambda6"], "1e-300"))
check("lambda7 = 0 is rejected", not F.gates(bad, fam)["lambda7_genuinely_nonzero"][0])
r = repr(float(good["lambda7"]))  # double-rounded lambda7 -> lambda1 reconstruction at large tan(beta)
bad = copy.deepcopy(good); bad.update(core.identity_checks(int(good["mphi_GeV"]), good["tan_beta"], good["lambda6"], r))
check("double-rounded lambda7 is caught by the lambda1 reconstruction gate",
      not F.gates(bad, fam)["lambda1_reconstructed_from_lambda7_50digit"][0])

# 2) wrong Y2: genuine free-M^2 point (rho = 0.01) through the same constructor family
fpD = core.S.build_generalized_family(tan_beta=good["tan_beta"], lambda1_ref=core.L1, mh_GeV=core.MH, mphi_GeV=str(int(good["mphi_GeV"])),
                                      mA_GeV=str(int(good["mphi_GeV"]) + 100), mHp_GeV=str(int(good["mphi_GeV"]) + 100), v_GeV=core.V,
                                      rho="0.01", lambda6=good["lambda6"], precision=50)
bad = copy.deepcopy(good)
from decimal import localcontext
with localcontext() as c:
    c.prec = 60
    bad["Y2_GeV2"] = str(fpD.values["Y2_GeV2"]); bad["Y2_is_exact_zero"] = Decimal(bad["Y2_GeV2"]) == 0
    bad["M2_minus_mh2_over_2"] = str(Decimal(str(fpD.values["M_sq_GeV2"])) - Decimal(core.MH) ** 2 / 2)
    bad["g_minus_2m2_over_v"] = str(Decimal(str(fpD.values["g_hphiphi_GeV"])) - 2 * Decimal(int(good["mphi_GeV"])) ** 2 / Decimal(core.V))
g = F.gates(bad, fam)
check("Y2 != 0 (rho=0.01) is rejected", not g["Y2_exact_zero"][0] and not g["M2_eq_mh2_over_2"][0] and not g["g_eq_2m2_over_v"][0])

# 3) legacy UFO normalisation in the production gate
src = open(H / "scripts/production.py").read()
gate = src[src.index("    for r in rs:  # hard provenance gates"):src.index("    assert len({r[\"model_couplings_sha256\"]")]
m = 200; gphys = 2 * m * m / 246.0
base = dict(GHphiphi=-gphys, MH=125.13, Mphi=m, pdlabel="nn23lo1", scale="-1", ebeam=6500.0, process="g g > H > h2 h2",
            gc90_is_i_times_GHphiphi=True, run="x")


def gate_ok(rec):
    try:
        exec(textwrap.dedent(gate), dict(rs=[rec], gphys=gphys, MH=125.13, m=m, abs=abs))
        return True
    except AssertionError:
        return False


check("physical kernel passes the production gate", gate_ok(dict(base)))
check("legacy +8 i m^2/v (GHphiphi = +4 g) is rejected", not gate_ok(dict(base, GHphiphi=+4 * gphys)))
check("sign-flipped coupling (GHphiphi = +g) is rejected", not gate_ok(dict(base, GHphiphi=+gphys)))
check("m_h = 125.20 kernel is rejected", not gate_ok(dict(base, MH=125.20)))

# 4) sigma from another mass
prod = {int(k): v for k, v in json.load(open(H / "raw/production.json")).items()}
check("kernel at 250 GeV does not match a 200 GeV point coupling", abs(prod[250]["g_hphiphi_GeV"] - gphys) > 1)

# 5) independent table recomputation
vb = list(csv.DictReader(open(H / "validated_benchmarks.csv")))
TA = {(r["family"], r["channel"]): r for r in csv.DictReader(open(H / "table_A.csv"))}
TB = {(r["family"], r["channel"]): r for r in csv.DictReader(open(H / "table_B.csv"))}
TI = list(csv.DictReader(open(H / "table_inputs.csv")))
check("no bare 'gg' channel key (gluons are 'gluglu')", not any("(gg)" in k[1] for k in (*TA, *TB)))
check("lambda7 row is flagged derived / do-not-re-input",
      all(r["role"] == "derived" and "NOT re-input" in r["note"] for r in TI if r["quantity"] == "lambda7"))
check("table A has 10 rows per family and B has 21", all(sum(1 for k in TA if k[0] == f) == 10 and sum(1 for k in TB if k[0] == f) == 21
                                                          for f in F.FAMILIES))
worst = 0.0
for b in vb:
    L, s = 300.0, float(b["sigma_LO_HEFT_fb"])
    v4 = np.array([float(b["BR_bb"]), float(b["BR_gluongluon"]), float(b["BR_Zgamma"]), float(b["BR_gammagamma"])])
    zb = np.array([0.69911, 0.033632, 0.033662])
    v6 = np.concatenate([v4[:2], v4[2] * zb, v4[3:]])
    for vec, names, T in ((v4, F.DIRECT, TA), (v6, F.EXP, TB)):
        M = L * s * np.outer(vec, vec)
        U = np.triu(M + M.T - np.diag(np.diag(M)))  # unordered pairs: off-diagonal doubled, diagonal once
        for i, a in enumerate(names):
            for j in range(i, len(names)):
                got = float(T[(b["family"], f"({a})({names[j]})")][f"{b['mphi_GeV']} GeV"])
                worst = max(worst, abs(got - U[i, j]) / max(U[i, j], 1e-300))
        check(f"{b['family']}@{b['mphi_GeV']}: total over pairs = L*sigma*(sum BR)^2",
              abs(U.sum() - L * s * vec.sum() ** 2) < 1e-9 * L * s)
        cov = next(r for r in TI if r["family"] == b["family"] and r["quantity"].startswith(f"Table {'A' if T is TA else 'B'} coverage"))
        check(f"{b['family']}@{b['mphi_GeV']}: published Table {'A' if T is TA else 'B'} coverage = (sum BR)^2",
              abs(float(cov[f"{b['mphi_GeV']} GeV"]) - vec.sum() ** 2) < 5e-5)
    check(f"{b['family']}@{b['mphi_GeV']}: Y2 = 0 and lambda7 != 0 in the published row",
          Decimal(b["Y2_GeV2"]) == 0 and Decimal(b["lambda7"]) != 0)
    check(f"{b['family']}@{b['mphi_GeV']}: sigma kernel belongs to this mass",
          abs(float(b["sigma_LO_HEFT_fb"]) - prod[int(b["mphi_GeV"])]["sigma_LO_HEFT_fb"]) == 0)
    check(f"{b['family']}@{b['mphi_GeV']}: ctau = hbar c / Gamma", abs(float(b["ctau_mm"]) * float(b["Gamma_phi_GeV"]) / core.HBARC - 1) < 1e-14)
    tot = sum(float(b[k]) for k in b if k.startswith("BR_") and k not in ("BR_other_terms",)) + float(b["BR_other_terms"])
    check(f"{b['family']}@{b['mphi_GeV']}: complete BR inventory sums to 1 ({tot - 1:+.1e})", abs(tot - 1) < 1e-12)
check(f"all table entries reproduced independently (max rel dev {worst:.1e})", worst < 1e-12)
check("Z-decay expansion: jj+ee+mumu fractions as documented", abs(sum(F.ZBR.values()) - (0.69911 + 0.033632 + 0.033662)) < 1e-15)
for f in F.FAMILIES:
    cls = {pts[k]["family_class"] for k in pts if k.split("_")[0] == f}
    check(f"{f}: every filled column is classified '{f}' by the same rule ({sorted(cls)})", cls <= {f})

print(f"\n{len(fails)} failure(s)")
sys.exit(1 if fails else 0)
