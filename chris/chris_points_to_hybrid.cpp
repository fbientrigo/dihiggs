// chris_points_to_hybrid.cpp
// ==========================
//
// Deterministic converter/validator for Christopher's 7 benchmark points.
//
// For each point given in the physical basis with a *target* lambda1
// (instead of m12^2):
//
//   1. Compute m12_sq (this is m12 SQUARED, in GeV^2 -- never the bare m12)
//      from the target lambda1 by inverting the lambda1 formula of
//      THDM::set_param_phys (THDM.cpp:354). This is the same inversion used
//      by the local patched THDM::set_param_phys_lam1 (THDM.cpp:414) and by
//      Christopher's standalone chris/CalcLambda1.c (m12sq_from_lambda1).
//   2. Build the model with THDM::set_param_phys(mh,mH,mA,mHp,sba,l6,l7,
//      m12_sq,tan_beta) and set the Yukawa type.
//   3. Extract the generic basis (get_param_gen), the Higgs basis
//      (get_param_higgs) and the hybrid basis (get_param_hybrid).
//      NB: in this local 2HDMC, get_param_hybrid's 3rd output is cba
//      (THDM.cpp:1141, cba = sign(sinba)*sqrt(1-sinba^2)); the header
//      parameter name "sba" at THDM.h:318 is misleading. Z7 = Lambda7 of
//      the Higgs basis, no extra sign (THDM.cpp:1162).
//   4. Run Constraints: positivity, unitarity, perturbativity, stability.
//   5. Compute minimal decay outputs for h2=H with DecayTable.
//   6. Round-trip: feed the extracted hybrid parameters into a second THDM
//      via set_param_hybrid(mh,mH,cba,Z4,Z5,Z7,tanb) and compare recovered
//      generic/physical parameters against the input.
//
// Known caveats (measured, not hidden):
//   * All 7 points have lambda6 != 0 (hard Z2 violation in the generic
//     basis), while the hybrid basis of 1507.04281 assumes a softly broken
//     Z2. get_param_hybrid prints "Model has hard Z_2-violation" and the
//     round-trip can only land on the nearest soft-Z2 model.
//   * tan_beta ~ 1e7 makes the lambda1 <-> m12_sq inversion ill-conditioned
//     (condition number ~ M11/(v^2 cb^2 lambda1) ~ 1e14), so lambda1_gen
//     can deviate from the target at the few-1e-2 absolute level in
//     IEEE double. The deltas are reported per point.
//
// Usage: ./chris_points_to_hybrid <input.csv> <output.csv>
// All numeric output uses std::setprecision(17), scientific.

#include "THDM.h"
#include "SM.h"
#include "Constraints.h"
#include "DecayTable.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

double qnan() { return std::numeric_limits<double>::quiet_NaN(); }

struct InputPoint {
  int point_id;
  double mh, mH, mA, mHp;
  double sba;
  double lambda6, lambda7;
  double lambda1_target;
  double tan_beta;
  int yukawa_type;
};

std::vector<std::string> split_csv(const std::string& line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(item);
  return out;
}

// m12^2 from the target lambda1: exact inverse of the lambda[1] expression
// in THDM::set_param_phys. Angles computed exactly as set_param_phys does
// (beta = atan(tb), alpha = beta - asin(sba)) so the inversion is consistent
// with what the library will recompute.
double m12sq_from_lambda1(const InputPoint& p, double v2) {
  const double beta = std::atan(p.tan_beta);
  const double cb = std::cos(beta);
  const double cb2 = cb * cb;
  const double tb = std::tan(beta);

  const double alpha = -std::asin(p.sba) + beta;
  const double sa = std::sin(alpha);
  const double ca = std::cos(alpha);

  const double M11 = p.mH * p.mH * ca * ca + p.mh * p.mh * sa * sa;

  return (M11 -
          v2 * cb2 * (p.lambda1_target + 1.5 * p.lambda6 * tb -
                      0.5 * p.lambda7 * tb * tb * tb)) /
         tb;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input.csv> <output.csv>\n";
    return 1;
  }

  std::ifstream fin(argv[1]);
  if (!fin) {
    std::cerr << "ERROR: cannot open input " << argv[1] << "\n";
    return 1;
  }

  std::vector<InputPoint> points;
  std::string line;
  std::getline(fin, line);  // header
  int id = 1;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    std::vector<std::string> f = split_csv(line);
    if (f.size() != 10) {
      std::cerr << "ERROR: row with " << f.size() << " fields: " << line << "\n";
      return 1;
    }
    InputPoint p;
    p.point_id = id++;
    p.mh = std::stod(f[0]);
    p.mH = std::stod(f[1]);
    p.mA = std::stod(f[2]);
    p.mHp = std::stod(f[3]);
    p.sba = std::stod(f[4]);
    p.lambda6 = std::stod(f[5]);
    p.lambda7 = std::stod(f[6]);
    p.lambda1_target = std::stod(f[7]);
    p.tan_beta = std::stod(f[8]);
    p.yukawa_type = static_cast<int>(std::stod(f[9]));
    points.push_back(p);
  }
  fin.close();

  std::ofstream fout(argv[2]);
  if (!fout) {
    std::cerr << "ERROR: cannot open output " << argv[2] << "\n";
    return 1;
  }
  fout.setf(std::ios::scientific);
  fout << std::setprecision(17);

  fout << "point_id,mh_input,mH_input,mA_input,mHp_input,sin_ba_input,"
          "cba_used,lambda6_input,lambda7_input,lambda1_target,"
          "tan_beta_input,yukawa_type,m12_sq_computed,set_param_phys_ok,"
          "lambda1_gen,lambda2_gen,lambda3_gen,lambda4_gen,lambda5_gen,"
          "lambda6_gen,lambda7_gen,m12_sq_gen,tan_beta_gen,"
          "hybrid_mh,hybrid_mH,hybrid_cba,hybrid_Z4,hybrid_Z5,hybrid_Z7,"
          "hybrid_tan_beta,set_param_hybrid_ok,"
          "roundtrip_lambda1_gen,delta_lambda1_roundtrip,"
          "roundtrip_mh,roundtrip_mH,roundtrip_mA,roundtrip_mHp,"
          "roundtrip_sba,roundtrip_m12_sq,"
          "positivity_ok_phys,unitarity_ok_phys,perturbativity_ok_phys,"
          "stability_ok_phys,positivity_ok_hybrid,unitarity_ok_hybrid,"
          "perturbativity_ok_hybrid,stability_ok_hybrid,"
          "total_width_h2,br_h2_gammagamma,gamma_h2_gammagamma,gamma_h2_hh,"
          "notes\n";

  std::cout.setf(std::ios::scientific);
  std::cout << std::setprecision(17);

  int n_rows = 0;
  for (size_t i = 0; i < points.size(); ++i) {
    const InputPoint& p = points[i];
    std::ostringstream notes;

    SM sm;
    const double v2 = sm.get_v2();

    // cba consistent with set_param_phys convention (cba >= 0 branch).
    const double cba_used = std::sqrt(std::max(0.0, 1.0 - p.sba * p.sba));
    if (std::fabs(p.sba) == 1.0) {
      notes << "exact alignment sba=1 (cba=0);";
    }
    if (p.tan_beta > 1.0e6) {
      notes << "huge tan_beta=" << std::setprecision(3) << p.tan_beta
            << std::setprecision(17)
            << " -> lambda1<->m12_sq inversion ill-conditioned;";
    }
    if (p.lambda6 != 0.0 || p.lambda7 != 0.0) {
      notes << "lambda6/7 nonzero: hard Z2 violation, hybrid basis is "
               "soft-Z2 approximation;";
    }

    const double m12_sq = m12sq_from_lambda1(p, v2);  // GeV^2 (m12 squared)

    THDM model;
    model.set_SM(sm);
    const bool phys_ok = model.set_param_phys(p.mh, p.mH, p.mA, p.mHp, p.sba,
                                              p.lambda6, p.lambda7, m12_sq,
                                              p.tan_beta);
    // Yukawa type after the potential: type-dependent rho matrices need beta.
    model.set_yukawas_type(p.yukawa_type);

    double l1g = qnan(), l2g = qnan(), l3g = qnan(), l4g = qnan(),
           l5g = qnan(), l6g = qnan(), l7g = qnan(), m12g = qnan(),
           tbg = qnan();
    double hyb_mh = qnan(), hyb_mH = qnan(), hyb_cba = qnan(), Z4 = qnan(),
           Z5 = qnan(), Z7 = qnan(), hyb_tb = qnan();
    int pos_p = -1, uni_p = -1, pert_p = -1, stab_p = -1;
    double w_tot = qnan(), w_gaga = qnan(), w_hh = qnan(), br_gaga = qnan();

    if (phys_ok) {
      model.get_param_gen(l1g, l2g, l3g, l4g, l5g, l6g, l7g, m12g, tbg);
      model.get_param_hybrid(hyb_mh, hyb_mH, hyb_cba, Z4, Z5, Z7, hyb_tb);

      if (std::fabs(l1g - p.lambda1_target) > 1.0e-6) {
        notes << "lambda1_gen deviates from target by "
              << std::setprecision(3) << (l1g - p.lambda1_target)
              << std::setprecision(17) << " (double-precision cancellation);";
      }

      Constraints check(model);
      pos_p = check.check_positivity() ? 1 : 0;
      uni_p = check.check_unitarity() ? 1 : 0;
      pert_p = check.check_perturbativity() ? 1 : 0;
      stab_p = check.check_stability() ? 1 : 0;

      DecayTable table(model);
      w_tot = table.get_gammatot_h(2);
      w_gaga = table.get_gamma_hgaga(2);
      w_hh = table.get_gamma_hhh(2, 1, 1);
      br_gaga = (std::isfinite(w_tot) && w_tot > 0.0) ? w_gaga / w_tot : qnan();
    } else {
      notes << "set_param_phys returned false;";
    }

    // ---- Hybrid round-trip ----
    bool hyb_ok = false;
    double rt_l1 = qnan(), d_l1 = qnan();
    double rt_mh = qnan(), rt_mH = qnan(), rt_mA = qnan(), rt_mHp = qnan(),
           rt_sba = qnan(), rt_m12sq = qnan();
    int pos_h = -1, uni_h = -1, pert_h = -1, stab_h = -1;

    if (phys_ok) {
      THDM hyb;
      hyb.set_SM(sm);
      hyb_ok = hyb.set_param_hybrid(hyb_mh, hyb_mH, hyb_cba, Z4, Z5, Z7, hyb_tb);
      if (hyb_ok) {
        hyb.set_yukawas_type(p.yukawa_type);

        double hl2, hl3, hl4, hl5, hl6, hl7, htb;
        hyb.get_param_gen(rt_l1, hl2, hl3, hl4, hl5, hl6, hl7, rt_m12sq, htb);
        d_l1 = rt_l1 - p.lambda1_target;

        double rl6, rl7, rtb;
        hyb.get_param_phys(rt_mh, rt_mH, rt_mA, rt_mHp, rt_sba, rl6, rl7,
                           rt_m12sq, rtb);

        Constraints hcheck(hyb);
        pos_h = hcheck.check_positivity() ? 1 : 0;
        uni_h = hcheck.check_unitarity() ? 1 : 0;
        pert_h = hcheck.check_perturbativity() ? 1 : 0;
        stab_h = hcheck.check_stability() ? 1 : 0;
      } else {
        notes << "set_param_hybrid returned false, no round-trip;";
      }
    }

    fout << p.point_id << ',' << p.mh << ',' << p.mH << ',' << p.mA << ','
         << p.mHp << ',' << p.sba << ',' << cba_used << ',' << p.lambda6
         << ',' << p.lambda7 << ',' << p.lambda1_target << ',' << p.tan_beta
         << ',' << p.yukawa_type << ',' << m12_sq << ',' << (phys_ok ? 1 : 0)
         << ',' << l1g << ',' << l2g << ',' << l3g << ',' << l4g << ',' << l5g
         << ',' << l6g << ',' << l7g << ',' << m12g << ',' << tbg << ','
         << hyb_mh << ',' << hyb_mH << ',' << hyb_cba << ',' << Z4 << ','
         << Z5 << ',' << Z7 << ',' << hyb_tb << ',' << (hyb_ok ? 1 : 0) << ','
         << rt_l1 << ',' << d_l1 << ',' << rt_mh << ',' << rt_mH << ','
         << rt_mA << ',' << rt_mHp << ',' << rt_sba << ',' << rt_m12sq << ','
         << pos_p << ',' << uni_p << ',' << pert_p << ',' << stab_p << ','
         << pos_h << ',' << uni_h << ',' << pert_h << ',' << stab_h << ','
         << w_tot << ',' << br_gaga << ',' << w_gaga << ',' << w_hh << ','
         << '"' << notes.str() << '"' << '\n';
    ++n_rows;

    std::cout << "[point " << p.point_id << "] mH=" << p.mH
              << " m12_sq=" << m12_sq << " phys_ok=" << phys_ok
              << " lambda1_gen=" << l1g << " hybrid_ok=" << hyb_ok
              << " roundtrip_lambda1=" << rt_l1 << "\n";
  }

  fout.close();
  std::cout << "Wrote " << n_rows << " rows to " << argv[2] << "\n";
  return (n_rows == 7) ? 0 : 2;
}
