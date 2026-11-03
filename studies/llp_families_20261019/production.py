"""Audit the existing physical-coupling MG5 kernels (read-only) -> production_provenance.csv + raw/production.json.
One kernel per mass: gg > H > h2 h2 with h2 = phi left stable, so the phi width/BRs do not enter and the kernel is
reusable for every Y2=0 benchmark at that mass (g = 2 m_phi^2 / v is then fixed).  No scaling between masses.
Run: python3 scripts/production.py
"""
import os
import csv, hashlib, json, re
from pathlib import Path

ROOT = Path(os.environ.get("DIHIGGS_ROOT") or Path(__file__).resolve().parents[3])
HERE = ROOT / "reports/llp_families_20261019"  # output package (lake), not code
MG = ROOT / "reports/lambda7_tables_20261014/mg5"
RUNS = {200: ("A_g0_s101", "A_g0_s107"), 250: ("A_m250_g0_s101", "A_m250_g0_s107"), 300: ("A_m300_g0_s101", "A_m300_g0_s107")}
V, MH = 246.0, 125.13


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def one(run):
    d = MG / "runs" / run
    banner = next((d / "Events").glob("*/*banner.txt"))
    b = banner.read_text()
    log = (MG / "logs" / f"{run}.log").read_text()
    xs, err = map(float, re.findall(r"Cross-section :\s+([\d.eE+-]+) \+- ([\d.eE+-]+) pb", log)[-1])
    g = lambda pat: re.search(pat, b, re.M).group(1)
    coup = (d / "Source/MODEL/couplings1.f").read_text()
    return dict(run=run, sigma_fb=xs * 1e3, err_fb=err * 1e3, banner=str(banner), banner_sha256=sha(banner),
                GHphiphi=float(g(r"^\s+3\s+(\S+)\s+# GHphiphi")), MH=float(g(r"^\s+25\s+(\S+)\s+# MH")),
                Mphi=float(g(r"^\s+9000006\s+(\S+)\s+# Mh2")), pdlabel=g(r"^\s+(\S+)\s+= pdlabel"),
                lhaid=g(r"^\s+(\S+)\s+= lhaid"), scale=g(r"^\s+(\S+)\s+= dynamical_scale_choice"),
                ebeam=float(g(r"^\s+(\S+)\s+= ebeam1")), seed=g(r"^\s+(\d+)\s+= iseed"),
                process=g(r"^generate (.+)$").strip(), nevents=int(g(r"Number of Events\s+:\s+(\d+)")),
                gc90_is_i_times_GHphiphi="GC_90 = MDL_COMPLEXI*MDL_GHPHIPHI" in coup,
                model_couplings_sha256=sha(d / "Source/MODEL/couplings1.f"),
                mg5_version=g(r"^#(\d+\.\d+\.\d+)\s*$"))


rows, out = [], {}
for m, runs in RUNS.items():
    rs = [one(r) for r in runs]
    gphys = 2 * m * m / V
    for r in rs:  # hard provenance gates: physical coupling, conventions, no legacy +8 i m^2/v vertex
        assert abs(r["GHphiphi"] + gphys) < 1e-9 * gphys, (r["run"], r["GHphiphi"])
        assert r["MH"] == MH and r["Mphi"] == m and r["pdlabel"] == "nn23lo1" and r["scale"] == "-1"
        assert r["ebeam"] == 6500.0 and r["process"] == "g g > H > h2 h2" and r["gc90_is_i_times_GHphiphi"]
    assert len({r["model_couplings_sha256"] for r in rs}) == 1 and all(r["mg5_version"] == "3.5.3" for r in rs)
    w = [1 / r["err_fb"] ** 2 for r in rs]
    mean = sum(wi * r["sigma_fb"] for wi, r in zip(w, rs)) / sum(w)
    err = sum(w) ** -0.5
    chi2 = sum(wi * (r["sigma_fb"] - mean) ** 2 for wi, r in zip(w, rs))
    out[m] = dict(sigma_LO_HEFT_fb=mean, mc_err_fb=err, seed_chi2_ndf1=chi2, g_hphiphi_GeV=gphys, runs=rs)
    rows.append(dict(mphi_GeV=m, sigma_LO_HEFT_fb=f"{mean:.6g}", mc_err_fb=f"{err:.2g}", seed_chi2=f"{chi2:.2f}",
                     g_hphiphi_GeV=f"{gphys:.6f}", GHphiphi_param_card=rs[0]["GHphiphi"],
                     feynman_rule="-i g (GC_90 = i*GHphiphi, GHphiphi = -g)", process="g g > H > h2 h2 (phi stable)",
                     ggh_vertex="HEFT, GH = -g_s^2/(12 pi^2 v) (m_t -> infinity)", sqrt_s_TeV=13, pdf="nn23lo1 (lhaid 230000)",
                     scale="dynamical_scale_choice=-1", m_h_GeV=MH, mg5="MG5_aMC 3.5.3",
                     seeds="+".join(r["seed"] for r in rs), events_per_seed=rs[0]["nevents"],
                     model_couplings_sha256=rs[0]["model_couplings_sha256"],
                     source_runs=";".join(str(MG / "runs" / r["run"]) for r in rs),
                     banner_sha256=";".join(r["banner_sha256"] for r in rs),
                     reuse_justification="same mass and g=2m^2/v for all Y2=0 points; phi stable in ME so width/BRs do not enter",
                     status="conditional LO-HEFT truth-level estimate (non-authoritative per CURRENT_CAMPAIGN.yaml)"))
with open(HERE / "production_provenance.csv", "w", newline="") as f:
    w = csv.DictWriter(f, list(rows[0])); w.writeheader(); w.writerows(rows)
json.dump(out, open(HERE / "raw/production.json", "w"), indent=1)
for m in out:
    print(m, f"{out[m]['sigma_LO_HEFT_fb']:.5f} +- {out[m]['mc_err_fb']:.5f} fb  chi2={out[m]['seed_chi2_ndf1']:.2f}")
