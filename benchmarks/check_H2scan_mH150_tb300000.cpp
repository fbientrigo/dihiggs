// One-use independent-construction check for point H2scan_mH150_tb300000.
//
// Unlike the canonical evaluator (dihiggs/app/Lambda1EvaluatorV2), which calls
// THDM::set_param_phys_lam1 (an internal m12_2-from-lambda1 inversion followed
// by a lambda1 round-trip check), this program calls THDM::set_param_phys
// directly with an externally supplied m12_2 (computed at high precision by
// benchmarks/check_H2scan_mH150_tb300000.py using the same closed-form
// inversion formula 2HDMC uses, but without double-precision cancellation).
// No 2HDMC source is modified; no Yukawa convention is changed
// (set_yukawas_type(1), identical to dihiggs::install_yukawa_type(model, 1)).
//
// Usage: check_H2scan_mH150_tb300000 mh mH mA mHp sba lambda6 lambda7 m12_2 tan_beta
#include "Constraints.h"
#include "DecayTable.h"
#include "THDM.h"

#include <cstdlib>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>

int main(int argc, char** argv) {
    if (argc != 10) {
        std::cerr << "Usage: " << argv[0]
                   << " mh mH mA mHp sba lambda6 lambda7 m12_2 tan_beta\n";
        return 1;
    }
    const double mh = std::atof(argv[1]);
    const double mH = std::atof(argv[2]);
    const double mA = std::atof(argv[3]);
    const double mHp = std::atof(argv[4]);
    const double sba = std::atof(argv[5]);
    const double lambda6 = std::atof(argv[6]);
    const double lambda7 = std::atof(argv[7]);
    const double m12_2 = std::atof(argv[8]);
    const double tan_beta = std::atof(argv[9]);

    THDM model;
    SM sm;
    model.set_SM(sm);

    const bool ok = model.set_param_phys(mh, mH, mA, mHp, sba, lambda6, lambda7, m12_2, tan_beta);
    std::cout << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    if (!ok) {
        std::cout << "construction_ok,0\n";
        return 0;
    }

    model.set_yukawas_type(1);
    const int yukawa_installed = model.get_yukawas_type();

    double l1, l2, l3, l4, l5, l6, l7, m12_2_rt, tanb_rt;
    model.get_param_gen(l1, l2, l3, l4, l5, l6, l7, m12_2_rt, tanb_rt);

    Constraints constraints(model);
    const int positivity_ok = constraints.check_positivity() ? 1 : 0;
    const int unitarity_ok = constraints.check_unitarity() ? 1 : 0;
    const int perturbativity_ok = constraints.check_perturbativity() ? 1 : 0;
    const int theory_ok = positivity_ok && unitarity_ok && perturbativity_ok;

    DecayTable decays(model);
    std::complex<double> coupling_h1_h2_h2;
    model.get_coupling_hhh(1, 2, 2, coupling_h1_h2_h2);
    const double width_bb = decays.get_gamma_hdd(2, 3, 3);
    const double width_cc = decays.get_gamma_huu(2, 2, 2);
    const double width_tautau = decays.get_gamma_hll(2, 3, 3);
    const double width_WW = decays.get_gamma_hvv(2, 3);
    const double width_ZZ = decays.get_gamma_hvv(2, 2);
    const double width_gammagamma = decays.get_gamma_hgaga(2);
    const double width_Zgamma = decays.get_gamma_hZga(2);
    const double width_gg = decays.get_gamma_hgg(2);
    const double width_hh = decays.get_gamma_hhh(2, 1, 1);
    const double total_width = decays.get_gammatot_h(2);
    const double hbar_c_gev_mm = 1.973269804e-13;
    const double ctau_mm = hbar_c_gev_mm / total_width;

    std::cout << "construction_ok,1"
              << ",yukawa_type_installed," << yukawa_installed
              << ",lambda1_reconstructed," << l1
              << ",m12_sq_input_gev2," << m12_2
              << ",m12_sq_reconstructed_gev2," << m12_2_rt
              << ",tan_beta_reconstructed," << tanb_rt
              << ",positivity_ok," << positivity_ok
              << ",unitarity_ok," << unitarity_ok
              << ",perturbativity_ok," << perturbativity_ok
              << ",theory_ok," << theory_ok
              << ",width_bb_gev," << width_bb
              << ",width_cc_gev," << width_cc
              << ",width_tautau_gev," << width_tautau
              << ",width_WW_gev," << width_WW
              << ",width_ZZ_gev," << width_ZZ
              << ",width_gammagamma_gev," << width_gammagamma
              << ",width_Zgamma_gev," << width_Zgamma
              << ",width_gg_gev," << width_gg
              << ",width_hh_gev," << width_hh
              << ",total_width_gev," << total_width
              << ",ctau_mm," << ctau_mm
              << ",br_bb," << (width_bb / total_width)
              << ",br_gammagamma," << (width_gammagamma / total_width)
              << ",br_Zgamma," << (width_Zgamma / total_width)
              << ",br_tautau," << (width_tautau / total_width)
              << ",br_gg," << (width_gg / total_width)
              << ",g_h1h2h2_real_gev," << std::real(coupling_h1_h2_h2)
              << ",g_h1h2h2_imag_gev," << std::imag(coupling_h1_h2_h2)
              << "\n";
    return 0;
}
