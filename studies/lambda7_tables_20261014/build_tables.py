"""Build BENCHMARKS/DECAYS/PRODUCTION/EVENT_YIELDS_* CSVs, provenance manifest and LaTeX
tables from raw/decays_recomputed.json (2HDMC) and mg5/runs/* (MadGraph banners).

  python3 scripts/build_tables.py        (stdlib only; asserts the physics checks)
"""
import os
import csv, glob, hashlib, itertools, json, math, re, sys, tempfile
from pathlib import Path

ROOT = Path(os.environ.get("DIHIGGS_ROOT") or Path(__file__).resolve().parents[3])
HERE = ROOT / "reports/lambda7_tables_20261014"  # output package (lake), not code
CONS = ROOT / "dihiggs/campaigns/lambda7_consolidation/20261011T1530Z"
UFO = Path("/tmp/dihiggs-honly-replay-20261007T012700Z/ufo_externalized")
HBARC = 1.973269804e-13  # GeV mm
LUMI = 300.0             # fb^-1
# PDG 2024 (S. Navas et al., Phys. Rev. D 110 (2024) 030001), Z boson summary table
BR_Z = {"jj": 0.69911, "ee": 0.033632, "mumu": 0.033662}

# report id -> (source point in raw json, family, MG5 runs giving sigma, comparison role)
POINTS = [
    ("A1", "B02", "A", ["A_g0_s101", "A_g0_s107"], "Y2=0, prompt, lambda6=0 (hadronic reference)"),
    ("A2", "B05", "A", ["A_g0_s101", "A_g0_s107"], "Y2=0, prompt, lambda6=-0.2 (Z7=0.20, H+ loop)"),
    ("A3", "B09", "A", ["A_g0_s101", "A_g0_s107"], "Y2=0, ctau~100 mm, lambda6=0"),
    ("A4", "B17", "A", ["A_g0_s101", "A_g0_s107"], "Y2=0, ctau~100 mm, lambda6=-2.64e-4 (photon-rich)"),
    ("A5", "B03", "A", ["A_m250_g0_s101", "A_m250_g0_s107"], "Y2=0, mphi=250 GeV prompt reference"),
    ("A6", "B10", "A", ["A_m250_g0_s101", "A_m250_g0_s107"], "Y2=0, mphi=250 GeV mixed LLP, ctau~100 mm"),
    ("A7", "B18", "A", ["A_m300_g0_s101", "A_m300_g0_s107"], "Y2=0, mphi=300 GeV photon-rich LLP, ctau~100 mm"),
    ("B1", "B02_L7ZERO", "B", ["B02L7Z_s101"], "lambda7=0 control of A1 (lambda1, lambda6, spectrum fixed; M2 solved)"),
    ("B2", "B13", "B", ["B13_g0p5_s101"], "spectrum of A1, M2 such that g=0.5 g0"),
    ("B3", "B14", "B", ["B14_g1p4_s101"], "spectrum of A1, M2 such that g=1.4 g0"),
    ("B4", "B17_L7ZERO", "B", ["B17L7Z_s101"], "lambda7=0 control of A4 (lambda1, lambda6, spectrum fixed; M2 solved)"),
    ("B5", "B18_L7ZERO", "B", ["B_m300_l7zero_s101", "B_m300_l7zero_s107"],
     "mphi=300 GeV, pre-existing lambda7=0 control of B18; lambda1/lambda6/spectrum fixed, M2 solved"),
]
CH = ("bb", "cc", "tautau", "mumu", "gg", "diphoton", "Zgamma", "WW", "ZZ", "hh")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def tree_sha(d):
    h = hashlib.sha256()
    for f in sorted(Path(d).rglob("*.py")):
        h.update(hashlib.sha256(f.read_bytes()).hexdigest().encode() + b"  ./" + str(f.relative_to(d)).encode() + b"\n")
    return h.hexdigest()


def mg5_run(rid, expected_g, expected_mphi):
    run = HERE / "mg5/runs" / rid
    banner = next(iter(glob.glob(str(run / "Events" / rid / "*banner.txt"))))
    txt = open(banner).read()
    xs = float(re.search(r"Integrated weight \(pb\)\s*:\s*(\S+)", txt).group(1))
    res = open(run / "SubProcesses/results.dat").readline().split()
    err = float(res[1])
    assert math.isclose(float(res[0]), xs, rel_tol=1e-4), (rid, res[0], xs)
    gh = float(re.search(r"^\s+3\s+(\S+)\s+# GHphiphi", txt, re.M).group(1))
    mh = float(re.search(r"^\s+25\s+(\S+)\s+# MH", txt, re.M).group(1))
    mphi = float(re.search(r"^\s+9000006\s+(\S+)\s+# Mh2", txt, re.M).group(1))
    card = (run / "Cards/param_card.dat").read_text()
    width = float(re.search(r"^DECAY\s+9000006\s+(\S+)", card, re.M).group(1))
    seed = int(re.search(r"^\s*(\d+)\s*=\s*iseed", txt, re.M).group(1))
    for key in ("nn23lo1 = pdlabel", "230000 = lhaid", "6500.0 = ebeam1", "6500.0 = ebeam2",
                "False = fixed_ren_scale", "False = fixed_fac_scale", "-1 = dynamical_scale_choice",
                "generate g g > H > h2 h2"):
        assert re.search(re.escape(key).replace(r"\ ", r"\s+"), txt), (rid, key)
    assert math.isclose(-gh, expected_g, rel_tol=1e-10) and math.isclose(mphi, expected_mphi, abs_tol=1e-9), rid
    source_couplings = (run / "Source/MODEL/couplings1.f").read_text()
    assert re.search(r"GC_90\s*=\s*MDL_COMPLEXI\*MDL_GHPHIPHI", source_couplings, re.I), rid
    # The project's banner checker requires an exact SLHA comment. Normalize only
    # that comment in a temporary copy; the original banner hash remains recorded.
    sys.path.insert(0, str(ROOT / "dihiggs_hep_cross/scripts"))
    from madgraph_banner_verification import banner_checks
    normalized = re.sub(r"(^\s*3\s+[\d.eE+\-]+\s+#\s*)GHphiphi.*$", r"\1GHphiphi", txt, flags=re.M)
    with tempfile.NamedTemporaryFile(mode="w", suffix="_banner.txt") as temp:
        temp.write(normalized); temp.flush()
        checks = banner_checks(Path(temp.name), {"g_hH2H2_GeV": expected_g, "mh_GeV": 125.13,
                                                  "mH2_GeV": expected_mphi, "total_width_GeV": width})
    assert all(checks.values()), (rid, checks)
    return dict(run=rid, sigma_pb=xs, err_pb=err, GHphiphi=gh, mh=mh, mphi=mphi, seed=seed,
                total_width_card_GeV=width, banner_checks=checks,
                model_vertex_source="GC_90 = i*GHphiphi (verified in compiled Source/MODEL/couplings1.f)",
                banner=str(Path(banner).relative_to(HERE)), banner_sha256=sha(banner))


def combine(runs):
    w = [1 / r["err_pb"] ** 2 for r in runs]
    m = sum(wi * r["sigma_pb"] for wi, r in zip(w, runs)) / sum(w)
    return m * 1e3, (1 / sum(w)) ** 0.5 * 1e3  # fb


raw = {r["id"]: r for r in json.load(open(HERE / "raw/decays_recomputed.json"))}
ttstar_audit = json.load(open(HERE / "raw/decay_accounting_ttstar.json"))
cons_bm = {r["benchmark_id"]: r for r in csv.DictReader(open(CONS / "BENCHMARKS.csv"))}
cons_dc = {r["benchmark_id"]: r for r in csv.DictReader(open(CONS / "DECAYS.csv"))}
ufo_hash = tree_sha(UFO)

bm_rows, dc_rows, pr_rows, ya_rows, yb_rows, manifest = [], [], [], [], [], {}
for pid, src, fam, runs, role in POINTS:
    r = raw[src]; i, v = r["inputs"], r["values"]
    mh, vv, mphi = float(i["mh_GeV"]), float(i["v_GeV"]), float(i["mphi_GeV"])
    M2, Y2, g = float(v["M_sq_GeV2"]), float(v["Y2_GeV2"]), float(v["g_hphiphi_GeV"])
    g0 = 2 * mphi**2 / vv
    # --- validation: Y2 identity, trilinear identity, mass reconstruction, alignment, tree tests
    assert abs(Y2 - (M2 - mh**2 / 2)) < 1e-6 * max(1, abs(M2)), pid
    assert math.isclose(g, (mh**2 + 2 * mphi**2 - 2 * M2) / vv, rel_tol=1e-9), pid
    assert math.isclose(g, 2 * (mphi**2 - Y2) / vv, rel_tol=1e-9, abs_tol=1e-9), pid
    for k, ref in (("h", mh), ("phi", mphi), ("A", float(i["mA_GeV"])), ("Hp", float(i["mHp_GeV"]))):
        if f"{k}_mass_reconstructed_GeV" in v:  # study-D constructor only exposes the PASS flag below
            assert math.isclose(float(v[f"{k}_mass_reconstructed_GeV"]), ref, rel_tol=1e-9), (pid, k)
    assert r["tree_status"] == "PASS_ALL_TESTED" and all(s == "PASS" for s in r["constraints"].values()), pid
    reconstructed = (("h_mass_reconstructed_GeV", mh), ("phi_mass_reconstructed_GeV", mphi),
                     ("A_mass_reconstructed_GeV", float(i["mA_GeV"])),
                     ("Hp_mass_reconstructed_GeV", float(i["mHp_GeV"])))
    if all(k in v for k, _ in reconstructed):
        for k, ref in reconstructed:
            assert math.isclose(float(v[k]), ref, rel_tol=1e-12), (pid, k, v[k], ref)
    else:
        # Study-D controls do not export reconstructed masses; retain their separate
        # constructor spectrum flag and mark the missing backend mass audit in provenance.
        assert all(s == "PASS" for key, s in r["constraints"].items() if key == "spectrum_reconstruction"), pid
    assert math.isclose(float(v["CphiHpHm_potential_GeV"]), -vv * float(v["Z7"]), rel_tol=1e-9), pid  # C_V=-v Z7
    if fam == "A":
        assert abs(Y2) < 1e-6 and math.isclose(g, g0, rel_tol=1e-12), pid
    status = cons_bm[src]["status"] if src in cons_bm else "PROVISIONAL"
    if src in cons_dc:  # reproduces the 2026-10-11 consolidation row
        assert math.isclose(float(cons_dc[src]["ctau_mm"]), r["ctau_mm"], rel_tol=1e-10), pid
    if src == "B18_L7ZERO":
        ab = next(x for x in csv.DictReader(open(CONS / "LAMBDA7_AB_COMPARISON.csv"))
                  if x["benchmark_id"] == "B18" and x["case"] == "C2")
        assert math.isclose(float(ab["ctau_mm"]), r["ctau_mm"], rel_tol=1e-10), pid
        assert math.isclose(float(ab["Gamma_tot_GeV"]), r["width_total_GeV"], rel_tol=1e-10), pid
        assert math.isclose(float(ab["BR_diphoton"]), float(r["widths_GeV"]["diphoton"]) / r["width_total_GeV"], rel_tol=1e-10), pid
    bm_rows.append(dict(
        id=pid, family=fam, role=role, source_id=src, mh_GeV=i["mh_GeV"], v_GeV=i["v_GeV"], mphi_GeV=i["mphi_GeV"],
        mA_GeV=i["mA_GeV"], mHp_GeV=i["mHp_GeV"], sin_beta_minus_alpha=1, tan_beta=i["tan_beta"],
        **{f"lambda{k}": f"{float(v[f'lambda{k}']):.10e}" for k in range(1, 8)},
        M2_GeV2=f"{M2:.10e}", m12_sq_GeV2=f"{float(v['m12_sq_GeV2']):.10e}", Y2_GeV2=f"{Y2:.10e}",
        Z7=f"{float(v['Z7']):.10e}", C_phiHpHm_GeV=f"{float(v['CphiHpHm_potential_GeV']):.10e}",
        g_hphiphi_GeV=f"{g:.10f}", g0_GeV=f"{g0:.10f}", delta_g_over_g0=f"{g / g0 - 1:.6f}",
        lambda7_T1_identity_abs_residual=f"{r['t1_lambda7_abs_residual']:.1e}",
        **{f"flag_{k}": s for k, s in r["constraints"].items()},
        flag_independent_backend_mass_audit=("FAIL_MASS_TOLERANCE" if pid == "A4" else
                                             "FAIL_MASS_TOLERANCE" if pid in ("A7", "B5") else "NOT_RUN"),
        flag_independent_backend_alignment_audit=("PASS_ANALYTIC_ROUTE" if pid == "A4" else "NOT_RUN"),
        flag_experimental_limits="NOT_EVALUATED", status=status,
        provenance=f"2HDMC analytic-bridge via lambda7_consolidation/20261011T1530Z consolidate.py (repo dihiggs@087ba24); "
                   f"source point {src}; re-evaluated 2026-10-14 by scripts/recompute_decays.py"))

    W, tot = r["widths_GeV"], r["width_total_GeV"]
    br = {c: W[c] / tot for c in CH}
    ttstar_width = (ttstar_audit["records"][src]["ttstar_width_GeV"]
                    if src in ttstar_audit["records"] else 0.0)
    br_ttstar = ttstar_width / tot
    br_unrep = 1 - sum(br.values()) - br_ttstar
    ctau = HBARC / tot
    assert math.isclose(ctau, r["ctau_mm"], rel_tol=1e-8), pid
    assert abs(sum(br.values()) + br_ttstar + br_unrep - 1) < 1e-12 and 0 <= br_unrep < 3e-2, (pid, br_unrep)  # remainder is not renormalized
    other = sum(br[c] for c in CH if c not in ("bb", "gg", "diphoton", "Zgamma")) + br_ttstar
    dc_rows.append(dict(id=pid, family=fam, Gamma_total_GeV=f"{tot:.6e}", ctau_mm=f"{ctau:.6e}",
                        ctau_2HDMC_mm=f"{r['ctau_mm']:.6e}",
                        BR_bb=f"{br['bb']:.6e}", BR_gg=f"{br['gg']:.6e}", BR_gammagamma=f"{br['diphoton']:.6e}",
                        BR_Zgamma=f"{br['Zgamma']:.6e}", BR_cc=f"{br['cc']:.6e}", BR_tautau=f"{br['tautau']:.6e}",
                        BR_mumu=f"{br['mumu']:.6e}", BR_WW=f"{br['WW']:.6e}", BR_ZZ=f"{br['ZZ']:.6e}", BR_hh=f"{br['hh']:.6e}",
                        BR_ttstar=f"{br_ttstar:.6e}", BR_other_sum=f"{other:.6e}", BR_unreported_residual=f"{br_unrep:.2e}",
                        provenance="2HDMC 1.8 patched; tt* from native indexed DecayTable accounting audit for 300 GeV; widths in GeV"))

    mg = [mg5_run(x, g, mphi) for x in runs]
    if len(mg) > 1:
        from madgraph_banner_verification import compatible
        ok, diff, err_combined = compatible(mg[0]["sigma_pb"], mg[0]["err_pb"],
                                            mg[1]["sigma_pb"], mg[1]["err_pb"])
        assert ok, (pid, diff, err_combined)
    for m in mg:
        assert math.isclose(-m["GHphiphi"], g, rel_tol=1e-9) or fam == "A" and math.isclose(-m["GHphiphi"], g0, rel_tol=1e-9), (pid, m)
        assert m["mh"] == 125.13 and m["mphi"] == mphi, (pid, m)
    sig, err = combine(mg)
    direct_owner = {"A1": "A1 (200 GeV, g0; shared kernel A1-A4)",
                    "B1": "B1 (200 GeV, g=63.648906 GeV)",
                    "B2": "B2 (200 GeV, g=0.5 g0)", "B3": "B3 (200 GeV, g=1.4 g0)",
                    "A5": "A5 (250 GeV, g0; shared kernel A5-A6)",
                    "A7": "A7 (300 GeV, g0)", "B4": "B4 (200 GeV, g=63.648443 GeV)",
                    "B5": "B5 (300 GeV, lambda7=0 control, direct integration)"}
    reused_from = {"A2": "A1", "A3": "A1", "A4": "A1", "A6": "A5"}
    pr_rows.append(dict(id=pid, family=fam, sigma_LO_fb=f"{sig:.5g}", sigma_MC_err_fb=f"{err:.2g}",
                        type="direct integration" if pid in direct_owner else "reused identical kernel",
                        process="p p (g g) > H > h2 h2, s-channel h* only",
                        UFO=f"LLscalar_v2 externalized (GC_90 = i*GHphiphi), py-tree sha256 {ufo_hash[:16]}",
                        coupling_applied=f"GHphiphi = {mg[0]['GHphiphi']:.6f} GeV  (FR = -i g, g = {g:.6f} GeV)",
                        anchor=(direct_owner[pid] if pid in direct_owner else
                                f"identical mphi,g production kernel; reused from {reused_from[pid]}"),
                        seeds="+".join(str(m["seed"]) for m in mg),
                        settings="MG5_aMC 3.5.3, LO, sqrt(s)=13 TeV, nn23lo1 (lhaid 230000), dynamical_scale_choice=-1, "
                                 f"mh=125.13, Gamma_h=4.3655e-3 GeV, phi_width_card={mg[0]['total_width_card_GeV']:.6e} GeV, "
                                 "external final-state width absent from production kernel; ggh vertex in heavy-top limit (GH=-g_s^2/(12 pi^2 v))",
                        runs=";".join(m["run"] for m in mg)))
    manifest[pid] = dict(source=src, family=fam, mg5_runs=mg, decays_raw="raw/decays_recomputed.json",
                         consolidation_status=status)

    # --- yields
    st = f"{status}_POINT; sigma LO-HEFT integrated, non-authoritative"
    mA = {"bb": br["bb"], "gg": br["gg"], "Zgamma": br["Zgamma"], "gammagamma": br["diphoton"]}
    for a, b in itertools.combinations_with_replacement(list(mA), 2):
        f = mA[a] ** 2 if a == b else 2 * mA[a] * mA[b]
        ya_rows.append(dict(id=pid, family=fam, pair=f"({a})({b})", BR_factor=f"{f:.6e}", sigma_fb=f"{sig:.5g}",
                            N_300fb=f"{LUMI * sig * f:.4g}", status=st))
    mB = {"bb": br["bb"], "gg": br["gg"], "jjgamma": br["Zgamma"] * BR_Z["jj"], "eegamma": br["Zgamma"] * BR_Z["ee"],
          "mumugamma": br["Zgamma"] * BR_Z["mumu"], "gammagamma": br["diphoton"]}
    for a, b in itertools.combinations_with_replacement(list(mB), 2):
        f = mB[a] ** 2 if a == b else 2 * mB[a] * mB[b]
        yb_rows.append(dict(id=pid, family=fam, pair=f"({a})({b})", BR_factor=f"{f:.6e}", sigma_fb=f"{sig:.5g}",
                            N_300fb=f"{LUMI * sig * f:.4g}", status=st))

assert len(ya_rows) == 10 * len(POINTS) and len(yb_rows) == 21 * len(POINTS)
for y in ya_rows + yb_rows:
    assert math.isclose(float(y["N_300fb"]), LUMI * float(y["sigma_fb"]) * float(y["BR_factor"]),
                        rel_tol=2e-3), (y["id"], y["pair"])


def write(name, rows):
    with open(HERE / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


write("BENCHMARKS.csv", bm_rows); write("DECAYS.csv", dc_rows); write("PRODUCTION.csv", pr_rows)
write("EVENT_YIELDS_DIRECT.csv", ya_rows); write("EVENT_YIELDS_EXPANDED.csv", yb_rows)
av = next(r for r in json.load(open(HERE / "raw/A4_backend_validation.json"))
          if r.get("mode") == "analytic" and r.get("precision_digits_constructed") == 50)
json.dump(dict(generated_by="scripts/build_tables.py", luminosity_fb=LUMI, BR_Z=BR_Z,
               BR_Z_source="PDG 2024, Phys. Rev. D 110 (2024) 030001 (Z -> hadrons, e+e-, mu+mu-)",
               ufo_dir=str(UFO), ufo_py_tree_sha256=ufo_hash,
               mg5_template="/tmp/dihiggs-honly-replay-cards-20261007T0130/mphi200_physical_seed101 (prepared 2026-10-07, never run)",
               mg5_engine="experiments/trilinear_normalization_20261010T191851Z/engine/mg5 (MG5_aMC 3.5.3)",
               consolidation=str(CONS), consolidation_manifest_sha256=sha(CONS / "MANIFEST.json"),
               inputs={p.name: sha(p) for p in sorted((ROOT / "_handoff/2hdm_lambda7_handoff/inputs").glob("*.tex"))},
               benchmark_validation={
                   "A4_analytic_backend": {
                       "artifact": "raw/A4_backend_validation.json",
                       "precision_digits": 50,
                       "constructor_status": av["constructor_status"],
                       "physical_state_ordering": av["physical_state_ordering"],
                       "mphi_abs_residual_GeV": av["backend_abs_residuals"]["m_H_backend_GeV"],
                       "mA_abs_residual_GeV": av["backend_abs_residuals"]["m_A_backend_GeV"],
                       "mHp_abs_residual_GeV": av["backend_abs_residuals"]["m_Hp_backend_GeV"],
                       "alignment_sin_abs_residual": av["backend_abs_residuals"]["sin_beta_minus_alpha"],
                       "alignment_cos_abs_residual": av["backend_abs_residuals"]["cos_beta_minus_alpha"],
                       "declared_mass_tolerance_GeV": "1e-8",
                       "interpretation": "alignment passes; returned mphi/mA/mHp residuals exceed tolerance; retain PROVISIONAL status"},
                   "B5_reused_2HDMC_control": {
                       "source": "dihiggs/campaigns/lambda7_consolidation/20261011T1530Z/LAMBDA7_AB_COMPARISON.csv, B18 C2",
                       "full_channel_re_evaluation": "raw/decays_recomputed.json; Gamma_total, c_tau and BR(gammagamma) matched the stored exact control",
                       "lambda7": "zero; lambda6 and Z7 retained as distinct nonzero values"},
                   "300_GeV_ttstar_accounting": {
                       "artifact": "raw/decay_accounting_ttstar.json",
                       "native_executable_sha256": ttstar_audit["executable_sha256"],
                       "source": ttstar_audit["source"],
                       "branching_fractions": {rid: ttstar_audit["records"][rid]["ttstar_width_GeV"] /
                                                ttstar_audit["records"][rid]["accounting_rerun"]["total_width_GeV"]
                                                for rid in ("B18", "B18_L7ZERO")},
                       "identity_probe": "strict input-mass check differs by 4.85e-8 GeV; accounting-only rerun uses returned native masses and is not a mass validation"}},
               benchmarks=manifest), open(HERE / "PROVENANCE.json", "w"), indent=1)


# ---------------- LaTeX tables ----------------
def sci(x, n=3):
    x = float(x)
    if x == 0:
        return "0"
    e = math.floor(math.log10(abs(x)))
    if -2 <= e <= 3:
        return f"{x:.{max(n - 1 - e, 0)}f}"
    return f"{x / 10**e:.{n - 1}f}\\times10^{{{e}}}"


T = HERE / "note/tables"; T.mkdir(parents=True, exist_ok=True)
ids = [p[0] for p in POINTS]
bmd = {r["id"]: r for r in bm_rows}; dcd = {r["id"]: r for r in dc_rows}; prd = {r["id"]: r for r in pr_rows}
cols = "l" + "r" * len(ids)
head = " & ".join(["Quantity"] + [f"\\textbf{{{x}}}" for x in ids]) + r" \\"


def row(label, fn):
    return " & ".join([label] + [f"${fn(x)}$" for x in ids]) + r" \\"


with open(T / "benchmarks.tex", "w") as f:
    f.write("\\begin{tabular}{" + cols + "}\n\\toprule\n" + head + "\n\\midrule\n")
    f.write(row(r"$\tan\beta$", lambda x: sci(bmd[x]["tan_beta"])) + "\n")
    f.write(row(r"$m_\phi$ [GeV]", lambda x: f"{float(bmd[x]['mphi_GeV']):.0f}") + "\n")
    f.write(row(r"$m_A$ [GeV]", lambda x: f"{float(bmd[x]['mA_GeV']):.0f}") + "\n")
    f.write(row(r"$m_{H^\pm}$ [GeV]", lambda x: f"{float(bmd[x]['mHp_GeV']):.0f}") + "\n")
    f.write(row(r"$\lambda_1$", lambda x: sci(bmd[x]["lambda1"])) + "\n")
    f.write(row(r"$\lambda_6$", lambda x: sci(bmd[x]["lambda6"])) + "\n")
    f.write(row(r"$\lambda_7$", lambda x: sci(bmd[x]["lambda7"]) if abs(float(bmd[x]["lambda7"])) > 1e-40 else "0")
            + "\n")
    f.write(row(r"$Z_7$", lambda x: sci(bmd[x]["Z7"])) + "\n")
    f.write(row(r"$M^2$ [GeV$^2$]", lambda x: sci(bmd[x]["M2_GeV2"], 4)) + "\n")
    f.write(row(r"$Y_2$ [GeV$^2$]", lambda x: sci(bmd[x]["Y2_GeV2"], 4) if abs(float(bmd[x]["Y2_GeV2"])) > 1e-6 else "0")
            + "\n")
    f.write(row(r"$g_{h\phi\phi}$ [GeV]", lambda x: f"{float(bmd[x]['g_hphiphi_GeV']):.1f}") + "\n")
    f.write(row(r"$(g-g_0)/g_0$", lambda x: f"{float(bmd[x]['delta_g_over_g0']):+.3f}") + "\n")
    f.write("\\midrule\n")
    f.write(row(r"$\Gamma_\phi$ [GeV]", lambda x: sci(dcd[x]["Gamma_total_GeV"])) + "\n")
    f.write(row(r"$c\tau_\phi$ [mm]", lambda x: sci(dcd[x]["ctau_mm"])) + "\n")
    for k, lab in (("bb", r"b\bar b"), ("gg", "gg"), ("gammagamma", r"\gamma\gamma"), ("Zgamma", r"Z\gamma"),
                   ("ttstar", r"t\bar t^*"), ("other_sum", r"\text{other (incl. }t\bar t^*\text{)}"),
                   ("unreported_residual", r"\text{unreported}")):
        csv_key = "BR_unreported_residual" if k == "unreported_residual" else f"BR_{k}"
        f.write(row(rf"BR$({lab})$", lambda x, csv_key=csv_key: sci(dcd[x][csv_key])) + "\n")
    f.write("\\midrule\n")
    f.write(row(r"$\sigma_{\rm LO}$ [fb]", lambda x: f"{float(prd[x]['sigma_LO_fb']):.4f}") + "\n")
    f.write(row(r"Independent backend mass audit", lambda x: "fail" if bmd[x]["flag_independent_backend_mass_audit"] == "FAIL_MASS_TOLERANCE" else "---") + "\n")
    f.write(row(r"Status", lambda x: r"\text{" + bmd[x]["status"][:4].lower() + "}") + "\n")
    f.write("\\bottomrule\n\\end{tabular}\n")


def yield_table(rows, name, labels):
    by = {}
    for r in rows:
        by.setdefault(r["pair"], {})[r["id"]] = r["N_300fb"]
    with open(T / name, "w") as f:
        f.write("\\begin{tabular}{" + cols + "}\n\\toprule\n" + head.replace("Quantity", "Channel") + "\n\\midrule\n")
        for pair, d in by.items():
            a, b = pair[1:-1].split(")(")
            f.write(" & ".join([f"$({labels[a]})({labels[b]})$"] + [f"${sci(d[x])}$" for x in ids]) + r" \\" + "\n")
        f.write("\\bottomrule\n\\end{tabular}\n")


L = {"bb": r"b\bar b", "gg": "gg", "Zgamma": r"Z\gamma", "gammagamma": r"\gamma\gamma", "jjgamma": r"jj\gamma",
     "eegamma": r"e^+e^-\gamma", "mumugamma": r"\mu^+\mu^-\gamma"}
yield_table(ya_rows, "yields_direct.tex", L)
yield_table(yb_rows, "yields_expanded.tex", L)
print("OK:", len(bm_rows), "benchmarks;", len(ya_rows), "+", len(yb_rows), "yield rows")
for pid in ids:
    print(pid, bmd[pid]["status"], "g/g0-1=", bmd[pid]["delta_g_over_g0"], "ctau=", dcd[pid]["ctau_mm"],
          "BRaa=", dcd[pid]["BR_gammagamma"], "BRZa=", dcd[pid]["BR_Zgamma"], "BRbb=", dcd[pid]["BR_bb"],
          "BRgg=", dcd[pid]["BR_gg"], "sigma=", prd[pid]["sigma_LO_fb"])
