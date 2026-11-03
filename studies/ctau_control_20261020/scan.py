"""Experiments A (constant-X trajectories, c tau targets) and B (local controls) at m_phi = 250 GeV.

  <venv>/bin/python -I scripts/scan.py        -> raw/scan_points.jsonl, scan_points.csv, raw/summary.json
Every row is a full backend evaluation (no formula shortcut); see ctl.evaluate.
"""
import csv, json, math, sys, traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ctl

OUT = ctl.HERE
TARGETS = (1e-3, 1e-2, 0.1, 1, 10, 100, 1e3, 1e4, 1e5, 1e6)  # 0.1 ... 1e4 requested; the outer decades probe the range
ROWS = []


def run(fam, X, t, experiment, note="", ref=None):
    try:
        r = ctl.evaluate(X, t, fam)
    except Exception as e:  # constructor/backend breakdown is a result, not a crash
        r = dict(family=fam, mphi_GeV=ctl.M, X=X, tan_beta=t if isinstance(t, str) else ctl.tstr(t), status="INCONCLUSIVE",
                 failed_gates="exception", failed_gate_details=f"{type(e).__name__}: {str(e)[:200]}")
        traceback.print_exc(limit=1)
    r.update(experiment=experiment, note=note)
    if ref is not None and "B_gamma" in r:
        r.update(ctl.dev(r, ref))
    ROWS.append(r)
    return r


def solve_t(fam, X, t0, c0, target, ref):
    """Secant in ln t on the MEASURED c tau, started from the t^2 guess; t kept at 6 significant digits."""
    lt, best = math.log(t0) + 0.5 * math.log(target / c0), None
    hist = []
    for _ in range(6):
        r = run(fam, X, math.exp(lt), "A_target_iter", f"target={target:g}", ref)
        if "ctau_mm" not in r:
            return r
        hist.append((lt, math.log(r["ctau_mm"])))
        if best is None or abs(math.log(r["ctau_mm"] / target)) < abs(math.log(best["ctau_mm"] / target)):
            best = r
        if abs(r["ctau_mm"] / target - 1) < 2e-5:
            break
        slope = 2.0 if len(hist) < 2 or hist[-1][0] == hist[-2][0] else (hist[-1][1] - hist[-2][1]) / (hist[-1][0] - hist[-2][0])
        lt += (math.log(target) - hist[-1][1]) / slope
    best = dict(best, experiment="A_target", note=f"target={target:g}", target_ctau_mm=target)
    ROWS.append(best)
    return best


def main():
    refs = ctl.references()
    summary = {}
    for fam, csvref in refs.items():
        X, t0 = csvref["X"], csvref["tan_beta"]
        ref = run(fam, X, t0, "reference", f"re-evaluation of {csvref['candidate_id']}")
        chk = {k: (ref[k], csvref[k]) for k in ("lambda7", "Z7", "lambda2", "m12_sq_GeV2")}
        summary[fam] = dict(ref_id=csvref["candidate_id"], X=X, tan_beta=t0, ctau_csv=float(csvref["ctau_mm"]), ctau_now=ref["ctau_mm"],
                            exact_string_match={k: a == b for k, (a, b) in chk.items()},
                            BR_match=max(abs(ref[k] - float(csvref[c])) for k, c in
                                         (("BR_bb", "BR_bb"), ("BR_gluglu", "BR_gluongluon"), ("BR_Zgamma", "BR_Zgamma"), ("BR_gammagamma", "BR_gammagamma"))),
                            status=ref["status"])
        print(fam, summary[fam], flush=True)
        # A: logarithmic trajectory, log10 t = 0 ... 13 in steps of 0.05 (reference included above)
        for k in range(0, 261):
            run(fam, X, 10 ** (k / 20), "A_grid", ref=ref)
        # A: lifetime targets
        for tg in TARGETS:
            r = solve_t(fam, X, float(t0), ref["ctau_mm"], tg, ref)
            print(f"  {fam} target {tg:g}: t={r.get('tan_beta')} ctau={r.get('ctau_mm', float('nan')):.6g} "
                  f"Bg={r.get('B_gamma', float('nan')):.6f} {r['status']} {r.get('failed_gates')}", flush=True)
        # B1: local slope d ln ctau / d ln t at fixed X, at the reference and at the 0.1 / 1e4 mm targets
        for base in [ref] + [r for r in ROWS if r["experiment"] == "A_target" and r["family"] == fam and r["target_ctau_mm"] in (0.1, 1e4)]:
            tb = float(base["tan_beta"])
            lo, hi = run(fam, X, tb * 0.999, "B1", f"base t={base['tan_beta']}", ref), run(fam, X, tb * 1.001, "B1", f"base t={base['tan_beta']}", ref)
            dl = math.log(float(hi["tan_beta"]) / float(lo["tan_beta"]))
            summary[fam].setdefault("B1", []).append(dict(
                base_t=base["tan_beta"], dln_ctau_dln_t=math.log(hi["ctau_mm"] / lo["ctau_mm"]) / dl,
                dB_gamma_dln_t=(hi["B_gamma"] - lo["B_gamma"]) / dl,
                **{f"dln_{k}_dln_t": math.log(hi[k] / lo[k]) / dl for k in ("Gamma_bb_GeV", "Gamma_gluglu_GeV", "Gamma_Zgamma_GeV", "Gamma_gammagamma_GeV")}))
        # B2: X at fixed tan beta (reference t); X = 0 is the lambda6 = 0 control
        if float(X) != 0:
            for fx in (-1.0, 0.0, 0.25, 0.5, 0.9, 0.999, 1.001, 1.1, 2.0, 4.0):
                run(fam, str(float(X) * fx), t0, "B2", f"X/X0={fx:g}", ref)
            lo = next(r for r in ROWS if r["experiment"] == "B2" and r["family"] == fam and r["note"] == "X/X0=0.999")
            hi = next(r for r in ROWS if r["experiment"] == "B2" and r["family"] == fam and r["note"] == "X/X0=1.001")
            dl = math.log(1.001 / 0.999)
            summary[fam]["B2_local"] = dict(dln_ctau_dln_X=math.log(hi["ctau_mm"] / lo["ctau_mm"]) / dl,
                                            dB_gamma_dln_X=(hi["B_gamma"] - lo["B_gamma"]) / dl,
                                            dln_Z7_dln_X=math.log(float(hi["Z7"]) / float(lo["Z7"])) / dl)
        print(fam, "B:", {k: summary[fam].get(k) for k in ("B1", "B2_local")}, flush=True)

    with open(OUT / "raw/scan_points.jsonl", "w") as f:
        for r in ROWS:
            f.write(json.dumps(r, default=str) + "\n")
    keys = []
    for r in ROWS:
        keys += [k for k in r if k not in keys]
    with open(OUT / "scan_points.csv", "w", newline="") as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        w.writerows(ROWS)
    json.dump(summary, open(OUT / "raw/summary.json", "w"), indent=1, default=str)
    print("rows", len(ROWS))


if __name__ == "__main__":
    main()
