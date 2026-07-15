#pragma once

#include "THDM.h"

struct PointResult {
    double m_phi;
    double M2_input;
    double m12_sq_input;

    double m12_sq_out;
    double lam1_out;
    double lam2_out;
    double lam3_out;
    double lam4_out;
    double lam5_out;
    double lam6_out;
    double lam7_out;

    bool positivity_ok;
    bool unitarity_ok;
    bool perturbativity_ok;
    bool stability_ok;
    bool theory_ok;
    bool triple_ok;

    double width_bb;
    double width_tautau;
    double width_WW;
    double width_ZZ;
    double width_gaga;
    double width_Zga;
    double width_gg;
    double width_hh;
    double total_width;
    double br_gaga;
};

PointResult evaluate_m2_point(
    THDM& model,
    double m_phi,
    double mH,
    double mA,
    double mHp,
    double sin_ba,
    double tan_beta,
    double lambda6,
    double lambda7,
    double M2_input
);
