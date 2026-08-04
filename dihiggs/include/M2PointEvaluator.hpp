#pragma once

#include "THDM.h"
#include <complex>
#include <limits>

struct PointResult {
    double m_phi = std::numeric_limits<double>::quiet_NaN();
    double M2_input = std::numeric_limits<double>::quiet_NaN();
    double m12_sq_input = std::numeric_limits<double>::quiet_NaN();
    bool construction_ok = false;
    double yukawa_type_installed = std::numeric_limits<double>::quiet_NaN();

    double m12_sq_out = std::numeric_limits<double>::quiet_NaN();
    double lam1_out = std::numeric_limits<double>::quiet_NaN();
    double lam2_out = std::numeric_limits<double>::quiet_NaN();
    double lam3_out = std::numeric_limits<double>::quiet_NaN();
    double lam4_out = std::numeric_limits<double>::quiet_NaN();
    double lam5_out = std::numeric_limits<double>::quiet_NaN();
    double lam6_out = std::numeric_limits<double>::quiet_NaN();
    double lam7_out = std::numeric_limits<double>::quiet_NaN();

    bool positivity_ok = false;
    bool unitarity_ok = false;
    bool perturbativity_ok = false;
    bool stability_ok = false;
    bool theory_ok = false;
    bool triple_ok = false;

    double width_bb = std::numeric_limits<double>::quiet_NaN();
    double width_tautau = std::numeric_limits<double>::quiet_NaN();
    double width_WW = std::numeric_limits<double>::quiet_NaN();
    double width_ZZ = std::numeric_limits<double>::quiet_NaN();
    double width_gaga = std::numeric_limits<double>::quiet_NaN();
    double width_Zga = std::numeric_limits<double>::quiet_NaN();
    double width_gg = std::numeric_limits<double>::quiet_NaN();
    double width_hh = std::numeric_limits<double>::quiet_NaN();
    double total_width = std::numeric_limits<double>::quiet_NaN();
    double br_gaga = std::numeric_limits<double>::quiet_NaN();
    // 2HDMC returns -i times the real scalar coefficient for h1-h2-h2.
    std::complex<double> coupling_h1h2h2_native = {
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};
};

PointResult evaluate_m2_point(
    THDM& model,
    double mh,
    double mH,
    double mA,
    double mHp,
    double sin_ba,
    double tan_beta,
    double lambda6,
    double lambda7,
    double M2_input
);
