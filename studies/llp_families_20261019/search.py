"""Phases B/C: composition search at fixed ctau = 100 mm on the Y2 = 0 branch.

For each (m_phi, sign of lambda6, tan beta) lambda6 is solved so that ctau = 100 mm (bisection in log|lambda6|).
Coarse log grid in tan beta -> raw/scan_ledger.csv; then tan beta is refined (bisection in log t on B_gamma)
to the historical mass-specific composition of each family, tan beta rounded to 3 significant digits and
lambda6 re-solved and rounded to 6 significant digits (exact decimal inputs).  Hadronic: lambda6 = 0, t solved.
Output: raw/scan_ledger.csv, raw/search_selected.json.   Run: <venv>/bin/python -I scripts/search.py
"""
import csv, json, math, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import core

CT = 100.0
MASSES = (200, 250, 300)
TARGETS = {  # historical compositions (bb, gluglu, Zgamma, gammagamma), transcribed from the task specification of 2026-10-19
    # (no frozen artifact in the repo carries them); used only to choose tan(beta) within a family, never as inputs
    "hadronic": {200: (0.675, 0.222, 0.0000799, 0.000509), 250: (0.578, 0.332, 0.000167, 0.000809)},
    "mixed": {200: (0.1053, 0.0350, 0.2705, 0.5732), 250: (0.0585, 0.0339, 0.3423, 0.5562), 300: (0.0344, 0.0345, 0.3828, 0.5411)},
    "photonic": {200: (0.009990, 0.002938, 0.2961, 0.6895), 250: (0.004933, 0.002663, 0.3704, 0.6212),
                 300: (0.002794, 0.002622, 0.4069, 0.5871)},
}
LEDGER = core.HERE / "raw/scan_ledger.csv"
FIELDS = ["phase", "mphi_GeV", "sign", "tan_beta", "lambda6", "X", "lambda7", "Z7", "ctau_mm", "BR_bb", "BR_gg",
          "BR_Zgamma", "BR_gammagamma", "B_h", "B_gamma", "family_class", "tree_status"]


def sig(x, n):
    return float(f"{x:.{n - 1}e}")


def fmt_t(t):
    return f"{sig(t, 3):.0f}"


def log(phase, r):
    new = not LEDGER.exists()
    with open(LEDGER, "a", newline="") as f:
        w = csv.DictWriter(f, FIELDS)
        if new:
            w.writeheader()
        w.writerow(dict(phase=phase, mphi_GeV=int(r["mphi_GeV"]), sign=("-" if float(r["lambda6"]) < 0 else "+" if float(r["lambda6"]) > 0 else "0"),
                        tan_beta=r["tan_beta"], lambda6=r["lambda6"], X=r["X_lambda6_tanbeta"], lambda7=f"{float(r['lambda7']):.10e}",
                        Z7=f"{float(r['Z7']):.10e}", ctau_mm=f"{r['ctau_mm']:.6g}", BR_bb=f"{r['BR_bb']:.6g}", BR_gg=f"{r['BR_gg']:.6g}",
                        BR_Zgamma=f"{r['BR_Zgamma']:.6g}", BR_gammagamma=f"{r['BR_diphoton']:.6g}", B_h=f"{r['B_h']:.6g}",
                        B_gamma=f"{r['B_gamma']:.6g}", family_class=r["family_class"], tree_status=r["tree_status"]))


def solve_l6(m, t, sign, digits=6):
    """|lambda6| such that ctau = CT; None if the fermionic width alone already exceeds the target (t too small)."""
    if core.evaluate(m, t, "0")["ctau_mm"] < CT:
        return None
    lo, hi = -9.0, -0.5
    if core.evaluate(m, t, f"{sign}{10 ** hi:.6g}")["ctau_mm"] > CT:
        return None
    for _ in range(45):
        mid = (lo + hi) / 2
        c = core.evaluate(m, t, f"{sign}{10 ** mid:.12g}")["ctau_mm"]
        lo, hi = (mid, hi) if c > CT else (lo, mid)
    return f"{sign}{sig(10 ** ((lo + hi) / 2), digits):.{digits - 1}e}"


def at(m, t, sign, phase):
    l6 = solve_l6(m, t, sign)
    if l6 is None:
        return None
    r = core.evaluate(m, t, l6)
    log(phase, r)
    return r


def solve_hadronic_t(m):
    """lambda6 = 0: ctau ~ t^2; fixed-point iteration on the measured ctau, t rounded to 3 digits."""
    t = 2.0e6
    for _ in range(6):
        c = core.evaluate(m, fmt_t(t), "0")["ctau_mm"]
        t = sig(t, 3) * math.sqrt(CT / c)
    r = core.evaluate(m, fmt_t(t), "0")
    log("hadronic_solve", r)
    return r


def dist(r, tgt):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip((r["BR_bb"], r["BR_gg"], r["BR_Zgamma"], r["BR_diphoton"]), tgt)))


def refine(m, sign, Bg_target, t_lo, t_hi):
    """B_gamma increases with t at fixed ctau (fermionic width ~ cot^2 beta): bisection in log t, then round t."""
    a, b = math.log(t_lo), math.log(t_hi)
    for _ in range(30):
        mid = (a + b) / 2
        r = at(m, f"{math.exp(mid):.6g}", sign, "refine")
        if r is None or r["B_gamma"] < Bg_target:
            a = mid
        else:
            b = mid
    return at(m, fmt_t(math.exp((a + b) / 2)), sign, "final_rounded")


if __name__ == "__main__":
    sel = {}
    for m in MASSES:
        had = solve_hadronic_t(m)
        sel[f"hadronic_{m}"] = dict(tan_beta=had["tan_beta"], lambda6="0", ctau_mm=had["ctau_mm"], family_class=had["family_class"],
                                    dist_to_historical=dist(had, TARGETS["hadronic"][m]) if m in TARGETS["hadronic"] else None)
        t_min = float(had["tan_beta"])
        grid = [t_min * 10 ** (k / 6) for k in range(1, 18)]  # coarse: t_min ... ~ 6e2 t_min (only t > t_min admits ctau=100 mm)
        for sign in ("-", "+"):
            for t in grid:
                at(m, f"{sig(t, 3):.0f}", sign, "coarse")
        for fam in ("mixed", "photonic"):
            tgt = TARGETS[fam][m]
            Bg = tgt[2] + tgt[3]
            for sign in ("-", "+"):
                r = refine(m, sign, Bg, t_min * 1.0001, t_min * 1000)
                sel[f"{fam}_{m}_{sign}"] = dict(tan_beta=r["tan_beta"], lambda6=r["lambda6"], ctau_mm=r["ctau_mm"],
                                                B_gamma=r["B_gamma"], B_h=r["B_h"], family_class=r["family_class"],
                                                dist_to_historical=dist(r, tgt), tree_status=r["tree_status"])
                print(fam, m, sign, sel[f"{fam}_{m}_{sign}"], flush=True)
        print("hadronic", m, sel[f"hadronic_{m}"], flush=True)
    json.dump(sel, open(core.HERE / "raw/search_selected.json", "w"), indent=1)
