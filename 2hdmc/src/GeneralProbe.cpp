// GeneralProbe.cpp -- NEW file, part of the physical_vs_effective_coupling_v1
// campaign. Forked from mission_runs/20260831_trilinear_lambda1_audit's
// AuditProbe.cpp (copied verbatim as a starting point, then extended below --
// AuditProbe.cpp itself is left untouched in that prior audit's own tree).
// Not part of upstream 2HDMC; added alongside the unmodified 2HDMC sources in
// this campaign's private copy to expose, via the public THDM API only, the
// quantities needed to test:
//     g_hphiphi =? 4*(m_phi^2 - m_2^2)/v
// for several candidate m_2^2 definitions, to scan lambda1 dependence, AND
// (new in this file) to test a general (non-aligned) closed-form coupling
// formula and an M2 exact-identity formula -- see the new columns appended
// at the end of HEADER below.
//
// It calls ONLY public THDM/Constraints methods. The one exception is the
// m22_2 "candidate" column, which is NOT read off a private member (2HDMC's
// `m22_2` is `private`) -- it is independently RE-DERIVED here using the same
// public algebra 2HDMC itself uses in THDM::set_param_phys (verbatim formula,
// see THDM.cpp:359), applied to the physical inputs. This file does not
// modify THDM.{h,cpp} in any way.
//
// Two modes:
//   phys   <mh> <mphi> <mA> <mHp> <sba> <lambda6> <lambda7> <m12_2>   <tan_beta> <yukawa_type>
//   lam1   <mh> <mphi> <mA> <mHp> <sba> <lambda1> <lambda6> <lambda7> <tan_beta> <yukawa_type>
//
// Output: one CSV line to stdout with a fixed schema (see HEADER below).

#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

static const char* HEADER =
    "mode,mh_in,mphi_in,mA_in,mHp_in,sba_in,lambda6_in,lambda7_in,"
    "m12_2_in,lambda1_in,tan_beta_in,yukawa_type,"
    "construction_ok,lam1_validation_available,lambda1_input_echo,"
    "lambda1_recomputed,lambda1_roundtrip_abs_error,lambda1_roundtrip_warning,"
    "lambda1_out,lambda2_out,lambda3_out,lambda4_out,lambda5_out,lambda6_out,"
    "lambda7_out,m12_2_out,tan_beta_out,"
    "m12_2_exact_used,M2_from_m12_2_exact_used,m22_2_candidate_GeV2,"
    "HiggsBasis_Z3,HiggsBasis_Z4,HiggsBasis_Z5,Z3_plus_Z4_plus_Z5,"
    "m11_2_candidate_GeV2,M22_2_Y2_candidate_GeV2,mHp_selfcheck_GeV,"
    "g_hphiphi_real_GeV,g_hphiphi_imag_GeV,g_from_Z345_alignment_formula_GeV,"
    "unitarity_ok,perturbativity_ok,stability_ok,positivity_ok,masses_ok,"
    "total_width_h2_GeV,ctau_h2_mm,v_GeV,"
    "Z1_higgs_basis,Z2_higgs_basis,Z6_higgs_basis,Z7_higgs_basis,"
    "sba_used,cba_used,g_general_formula_GeV,"
    "M2_computed_from_lam1_formula_GeV2,M2_actual_GeV2,"
    "hmass_h2_reconstructed_GeV,hmass_hA_reconstructed_GeV,hmass_hHp_reconstructed_GeV";

int main(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "--header") {
    std::cout << HEADER << "\n";
    return 0;
  }

  if (argc < 12) {
    std::cerr << "Usage:\n"
              << "  AuditProbe phys <mh> <mphi> <mA> <mHp> <sba> <lambda6> <lambda7> <m12_2> <tan_beta> <yukawa_type>\n"
              << "  AuditProbe lam1 <mh> <mphi> <mA> <mHp> <sba> <lambda1> <lambda6> <lambda7> <tan_beta> <yukawa_type>\n"
              << "  AuditProbe --header\n";
    return 2;
  }

  const std::string mode = argv[1];
  const double mh_in = std::atof(argv[2]);
  const double mphi_in = std::atof(argv[3]);
  const double mA_in = std::atof(argv[4]);
  const double mHp_in = std::atof(argv[5]);
  const double sba_in = std::atof(argv[6]);

  double lambda6_in = 0, lambda7_in = 0, m12_2_in = 0, lambda1_in = 0, tan_beta_in = 0;
  int yukawa_type = 1;

  THDM model;
  SM sm;
  model.set_SM(sm);

  bool construction_ok = false;

  if (mode == "phys") {
    lambda6_in = std::atof(argv[7]);
    lambda7_in = std::atof(argv[8]);
    m12_2_in = std::atof(argv[9]);
    tan_beta_in = std::atof(argv[10]);
    yukawa_type = std::atoi(argv[11]);
    construction_ok = model.set_param_phys(mh_in, mphi_in, mA_in, mHp_in, sba_in,
                                            lambda6_in, lambda7_in, m12_2_in, tan_beta_in);
  } else if (mode == "lam1") {
    lambda1_in = std::atof(argv[7]);
    lambda6_in = std::atof(argv[8]);
    lambda7_in = std::atof(argv[9]);
    tan_beta_in = std::atof(argv[10]);
    yukawa_type = std::atoi(argv[11]);
    construction_ok = model.set_param_phys_lam1(mh_in, mphi_in, mA_in, mHp_in, sba_in,
                                                 lambda1_in, lambda6_in, lambda7_in, tan_beta_in);
  } else {
    std::cerr << "unknown mode: " << mode << "\n";
    return 2;
  }

  // Defaults for a failed construction -- still emit a row (all-NaN payload)
  // so failed/theory-invalid points are visible in the CSV, not silently
  // dropped.
  auto nan = []() { return std::nan(""); };

  double lambda1_out = nan(), lambda2_out = nan(), lambda3_out = nan(), lambda4_out = nan(),
         lambda5_out = nan(), lambda6_out = nan(), lambda7_out = nan(), m12_2_out = nan(),
         tan_beta_out = nan();
  double m12_2_exact_used = nan(), M2_from_exact = nan(), m22_2_candidate = nan();
  double Z3 = nan(), Z4 = nan(), Z5 = nan(), Z345 = nan();
  double m11_2_candidate = nan(), M22_2_candidate = nan(), mHp_selfcheck = nan();
  double g_real = nan(), g_imag = nan(), g_from_Z345 = nan();
  bool unitarity_ok = false, perturbativity_ok = false, stability_ok = false,
       positivity_ok = false, masses_ok = false;
  double total_width_h2 = nan(), ctau_h2_mm = nan();
  // Added per validation/SOL_FINAL_REVIEW.md Additional issue #3: get_hmass()
  // routes through get_m12_2(), which round-trips to ~4.5% error at
  // tan_beta=4e7 (catastrophic cancellation, THDM.cpp:2786). These echo
  // 2HDMC's OWN reconstructed masses so every consumer of total_width_h2/
  // ctau_h2_mm can see directly whether they were computed at the requested
  // mass or a numerically-drifted one, instead of trusting mphi_in silently.
  double hmass_h2_reconstructed = nan(), hmass_hA_reconstructed = nan(),
         hmass_hHp_reconstructed = nan();

  bool lam1_val_avail = false;
  double lam1_val_input = nan(), lam1_val_recomputed = nan(), lam1_val_abs_err = nan();
  bool lam1_val_warning = false;

  // New columns (physical_vs_effective_coupling_v1 campaign, GeneralProbe
  // extension over AuditProbe): Higgs-basis Z1,Z2,Z6,Z7, the general
  // (non-aligned) coupling closed-form candidate, and the M2 exact-identity
  // closed-form candidate. See file header comment.
  double Z1_hb = nan(), Z2_hb = nan(), Z6_hb = nan(), Z7_hb = nan();
  double sba_used = nan(), cba_used = nan();
  double g_general_formula = nan();
  double M2_from_lam1_formula = nan();
  double M2_actual = nan();

  if (construction_ok) {
    model.set_yukawas_type(yukawa_type);

    model.get_param_gen(lambda1_out, lambda2_out, lambda3_out, lambda4_out, lambda5_out,
                         lambda6_out, lambda7_out, m12_2_out, tan_beta_out);
    // lambda1_out..lambda7_out are the class members lambda[1..7] verbatim
    // (THDM::get_param_gen just copies them out, THDM.cpp:1030-1042) --
    // exact, not a round-tripped/re-derived quantity. Only m12_2_out is
    // re-derived (via the private get_m12_2(), a DIFFERENT formula than the
    // one that constructed the model -- see the m12_2_exact_used note below).

    const double v = sm.get_v();
    const double v2 = v * v;
    const double beta = std::atan(tan_beta_in);
    const double sb = std::sin(beta), cb = std::cos(beta);
    const double sb2 = sb * sb, cb2 = cb * cb;
    const double s2b = std::sin(2. * beta), c2b = std::cos(2. * beta);
    const double tb = std::tan(beta), ctb = 1.0 / tb;

    // m12_2_exact_used: the ACTUAL m12_2 that fed set_param_phys internally
    // for THIS point -- not get_param_gen's re-derived m12_2_out (which used
    // a different formula, THDM.cpp:2774-2792 via get_m12_2(), and can
    // disagree with the true input at the ~1e-3 relative level for the
    // extreme tan_beta used in this project's benchmarks -- itself a real,
    // separately reported finding, see 05_M2_IDENTITY_TEST_NOTES).
    // - "phys" mode: m12_2 is a direct CLI input -> m12_2_exact_used = m12_2_in.
    // - "lam1" mode: replicate THDM::set_param_phys_lam1's own reconstruction
    //   formula verbatim (THDM.cpp:414-415), from lambda1_in and the same
    //   physical inputs -- this is the exact value THDM internally computed
    //   and fed to set_param_phys, re-derived here from public inputs only.
    if (mode == "phys") {
      m12_2_exact_used = m12_2_in;
    } else {
      const double alpha_l = -std::asin(sba_in) + beta;
      const double sa_l = std::sin(alpha_l), ca_l = std::cos(alpha_l);
      m12_2_exact_used = (mphi_in * mphi_in * ca_l * ca_l + mh_in * mh_in * sa_l * sa_l
                           - v2 * cb2 * (lambda1_in + 1.5 * lambda6_in * tb
                                         - 0.5 * lambda7_in * tb * tb * tb)) / tb;
    }

    // M2 = m12^2 / (sin(beta) cos(beta)) -- this project's own convention
    // (dihiggs_hep_cross/scripts/r9_run_model_scan.py,
    // dihiggs/dihiggs/src/M2PointEvaluator.cpp), NOT a native 2HDMC quantity.
    M2_from_exact = m12_2_exact_used / (sb * cb);

    // m22_2 candidate: THDM::set_param_phys's own public formula
    // (THDM.cpp:359), applied with the EXACT m12_2 used internally (not the
    // get_param_gen round-trip). This is NOT read from THDM's private
    // m22_2 member -- it is independently re-derived from public inputs.
    const double alpha = -std::asin(sba_in) + beta;
    const double sa = std::sin(alpha), ca = std::cos(alpha);
    const double cba = std::sqrt(1.0 - sba_in * sba_in);
    m22_2_candidate = -0.5 / sb * (mh_in * mh_in * ca * sba_in + mphi_in * mphi_in * sa * cba)
                       + m12_2_exact_used * ctb;

    // Higgs-basis Z3,Z4,Z5 (called Lambda3,4,5 in THDM::get_param_higgs,
    // THDM.cpp:1105-1138) and the Higgs-basis "Y2"-equivalent M22_2
    // (THDM.cpp:1133-1137), all re-derived verbatim from the same public
    // algebra applied to lambda1_out..lambda7_out (exact) and
    // m12_2_exact_used (exact). "See hep-ph/0504050" per THDM.cpp:1124.
    const double lambda345 = lambda3_out + lambda4_out + lambda5_out;
    const double s2b2 = s2b * s2b;
    Z3 = 0.25 * s2b2 * (lambda1_out + lambda2_out - 2 * lambda345) + lambda3_out
         - s2b * c2b * (lambda6_out - lambda7_out);
    Z4 = 0.25 * s2b2 * (lambda1_out + lambda2_out - 2 * lambda345) + lambda4_out
         - s2b * c2b * (lambda6_out - lambda7_out);
    Z5 = 0.25 * s2b2 * (lambda1_out + lambda2_out - 2 * lambda345) + lambda5_out
         - s2b * c2b * (lambda6_out - lambda7_out);
    Z345 = Z3 + Z4 + Z5;

    // Higgs-basis Z1, Z2, Z6, Z7 -- same re-derivation pattern as Z3/Z4/Z5
    // above, verbatim from THDM::get_param_higgs (THDM.cpp:1125-1131), plus
    // the sign flip on 6/7 only that THDM::get_coupling_hhh itself applies
    // (THDM.cpp:1580: "Z6=-l6; Z7=-l7;", Z1..Z5 unflipped).
    const double cb4 = cb2 * cb2, sb4 = sb2 * sb2;
    const double c3b = std::cos(3. * beta), s3b = std::sin(3. * beta);
    Z1_hb = lambda1_out * cb4 + lambda2_out * sb4 + 0.5 * lambda345 * s2b2
            + 2. * s2b * (cb2 * lambda6_out + sb2 * lambda7_out);
    Z2_hb = lambda1_out * sb4 + lambda2_out * cb4 + 0.5 * lambda345 * s2b2
            - 2. * s2b * (sb2 * lambda6_out + cb2 * lambda7_out);
    {
      const double Lambda6_hb = -0.5 * s2b * (lambda1_out * cb2 - lambda2_out * sb2 - lambda345 * c2b)
                                 + cb * c3b * lambda6_out + sb * s3b * lambda7_out;
      const double Lambda7_hb = -0.5 * s2b * (lambda1_out * sb2 - lambda2_out * cb2 + lambda345 * c2b)
                                 + sb * s3b * lambda6_out + cb * c3b * lambda7_out;
      Z6_hb = -Lambda6_hb;
      Z7_hb = -Lambda7_hb;
    }

    // General (non-aligned) closed-form coupling candidate, derived by hand
    // this mission from get_coupling_hhh's literal permutation sum
    // (physical_trilinear_derivation.md Section 1). sba_used/cba_used use
    // THDM::get_cba()'s own convention (THDM.cpp:2828-2831): cba = +sqrt(1-sba^2),
    // always the non-negative root.
    sba_used = sba_in;
    cba_used = std::sqrt(1.0 - sba_in * sba_in);
    {
      const double s = sba_used, cc = cba_used;
      // Sign fixed per validation/SOL_FINAL_REVIEW.md Additional issue #1:
      // this must be -v*(...) to equal Im[get_coupling_hhh(1,2,2)], not
      // +v*(...). The old +v*(...) sign was masked by a +/--tolerant
      // comparison in logs/run_all_probes.py (also fixed).
      g_general_formula = -v * (3. * s * cc * cc * Z1_hb
                                 + (s * s * s - 2. * s * cc * cc) * Z345
                                 + 3. * cc * (2. * s * s - cc * cc) * Z6_hb
                                 - 3. * cc * s * s * Z7_hb);
    }

    // M2 exact-identity candidate (assumptions_and_limits.md Section 2),
    // valid at exact alignment (sba_in==1) under the lambda1-input
    // construction route. lambda1_column: whichever lambda1 actually built
    // this point -- lambda1_in for "lam1" mode, the model's own
    // reconstructed lambda1_out for "phys" mode (report regardless of mode;
    // only meaningful when sba_in==1, left for the caller to gate on that).
    {
      const double lambda1_column = (mode == "lam1") ? lambda1_in : lambda1_out;
      M2_from_lam1_formula = mphi_in * mphi_in
          + (mh_in * mh_in - v2 * lambda1_column) / (tan_beta_in * tan_beta_in)
          - 1.5 * v2 * lambda6_in / tan_beta_in
          + 0.5 * v2 * lambda7_in * tan_beta_in;
    }
    M2_actual = M2_from_exact;

    m11_2_candidate = m12_2_exact_used * tb
        - 0.5 * v2 * (lambda1_out * cb2 + lambda345 * sb2
                      + 3. * lambda6_out * sb * cb + lambda7_out * sb2 * tb);
    M22_2_candidate = m11_2_candidate * sb2 + m22_2_candidate * cb2 + m12_2_exact_used * s2b;
    // THDM.cpp:1137: m_Hp = sqrt(M22_2 + 0.5*v2*Lambda3) -- self-consistency
    // check against the mHp we actually fed in (mHp_in). Large disagreement
    // here would mean an error in this re-derivation, not in 2HDMC.
    const double mhp2_check = M22_2_candidate + 0.5 * v2 * Z3;
    mHp_selfcheck = (mhp2_check >= 0) ? std::sqrt(mhp2_check) : nan();

    // h-phi-phi coupling: 2HDMC internal numbering is (1,2,3,4=h,H,A,H+) with
    // m_h < m_H always. Since mh_in < mphi_in in every point audited here,
    // our SM-like mediator (125.2 GeV) is 2HDMC index 1 ("h") and our dark
    // scalar (m_phi) is 2HDMC index 2 ("H"). get_coupling_hhh(1,2,2,c) is
    // therefore exactly the h-phi-phi vertex coefficient in 2HDMC's own
    // convention (hep-ph/0602242, per THDM.cpp's comment at get_coupling_hhh).
    std::complex<double> c;
    model.get_coupling_hhh(1, 2, 2, c);
    g_real = c.real();
    g_imag = c.imag();
    // Alignment-limit closed form derived by hand (subagent 2, Part IV/V
    // analytical findings, cross-checked against get_coupling_hhh's literal
    // permutation-sum code): at sba=1 exactly,
    //   get_coupling_hhh(1,2,2) = -i*v*(Z3+Z4+Z5).
    // Only meaningful/expected to match when sba_in == 1 exactly; reported
    // regardless so the CSV shows exactly where it stops matching away from
    // alignment.
    g_from_Z345 = -v * Z345;

    Constraints check(model);
    unitarity_ok = check.check_unitarity();
    perturbativity_ok = check.check_perturbativity();
    stability_ok = check.check_stability();
    positivity_ok = check.check_positivity();
    masses_ok = check.check_masses();
    hmass_h2_reconstructed = model.get_hmass(2);
    hmass_hA_reconstructed = model.get_hmass(3);
    hmass_hHp_reconstructed = model.get_hmass(4);

    if (mode == "lam1") {
      lam1_val_avail = model.has_param_phys_lam1_validation();
      if (lam1_val_avail) {
        model.get_param_phys_lam1_validation(lam1_val_input, lam1_val_recomputed,
                                              lam1_val_abs_err, lam1_val_warning);
      }
    }

    // Decay table for h2 (index 2 = our dark scalar) width/ctau, for
    // completeness (Part V-B asks for width/ctau columns). hbar_c in mm*GeV
    // matches physics_conventions.yaml's hbar_c_gev_mm.
    DecayTable table(model);
    total_width_h2 = table.get_gammatot_h(2);
    const double hbar_c_gev_mm = 1.973269804e-13;
    if (total_width_h2 > 0) {
      ctau_h2_mm = hbar_c_gev_mm / total_width_h2;
    }
  }

  std::printf(
      "%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,"
      "%d,%d,%.17g,%.17g,%.17g,%d,"
      "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
      "%.17g,%.17g,%.17g,"
      "%.17g,%.17g,%.17g,%.17g,"
      "%.17g,%.17g,%.17g,"
      "%.17g,%.17g,%.17g,"
      "%d,%d,%d,%d,%d,"
      "%.17g,%.17g,%.17g,"
      "%.17g,%.17g,%.17g,%.17g,"
      "%.17g,%.17g,%.17g,"
      "%.17g,%.17g,"
      "%.17g,%.17g,%.17g\n",
      mode.c_str(), mh_in, mphi_in, mA_in, mHp_in, sba_in, lambda6_in, lambda7_in, m12_2_in,
      lambda1_in, tan_beta_in, yukawa_type,
      construction_ok ? 1 : 0, lam1_val_avail ? 1 : 0, lam1_val_input, lam1_val_recomputed,
      lam1_val_abs_err, lam1_val_warning ? 1 : 0,
      lambda1_out, lambda2_out, lambda3_out, lambda4_out, lambda5_out, lambda6_out, lambda7_out,
      m12_2_out, tan_beta_out,
      m12_2_exact_used, M2_from_exact, m22_2_candidate,
      Z3, Z4, Z5, Z345,
      m11_2_candidate, M22_2_candidate, mHp_selfcheck,
      g_real, g_imag, g_from_Z345,
      unitarity_ok ? 1 : 0, perturbativity_ok ? 1 : 0, stability_ok ? 1 : 0,
      positivity_ok ? 1 : 0, masses_ok ? 1 : 0,
      total_width_h2, ctau_h2_mm, sm.get_v(),
      Z1_hb, Z2_hb, Z6_hb, Z7_hb,
      sba_used, cba_used, g_general_formula,
      M2_from_lam1_formula, M2_actual,
      hmass_h2_reconstructed, hmass_hA_reconstructed, hmass_hHp_reconstructed);

  return 0;
}
