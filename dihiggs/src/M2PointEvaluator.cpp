#include "M2PointEvaluator.hpp"
#include "Constraints.h"
#include "DecayTable.h"
#include <cmath>

PointResult evaluate_m2_point(
    THDM& model, 
    double mh,
    double mH_input, 
    double mA, 
    double mHp,
    double sin_ba, 
    double tan_beta,
    double lambda6, 
    double lambda7, 
    double M2_input
) {
    PointResult res;
    res.m_phi = mH_input;
    res.M2_input = M2_input;
    
    // High-precision trigonometric reconstruction to avoid loss at extreme tan_beta
    double beta = atan(tan_beta);
    double sin_beta = sin(beta);
    double cos_beta = cos(beta);
    
    // Explicit conversion M2 -> m12^2
    res.m12_sq_input = M2_input * sin_beta * cos_beta;
    
    model.set_yukawas_type(1);
    res.construction_ok = model.set_param_phys(
        mh, mH_input, mA, mHp, sin_ba, lambda6, lambda7, res.m12_sq_input, tan_beta);
    if (!res.construction_ok) return res;
    
    // Extract parameters
    double lam1_g, lam2_g, lam3_g, lam4_g, lam5_g, lam6_g, lam7_g, m12_2_g, tanb_g;
    model.get_param_gen(lam1_g, lam2_g, lam3_g, lam4_g, lam5_g, lam6_g, lam7_g, m12_2_g, tanb_g);
    
    res.m12_sq_out = m12_2_g;
    res.lam1_out = lam1_g;
    res.lam2_out = lam2_g;
    res.lam3_out = lam3_g;
    res.lam4_out = lam4_g;
    res.lam5_out = lam5_g;
    res.lam6_out = lam6_g;
    res.lam7_out = lam7_g;
    
    // Evaluate constraints
    Constraints check(model);
    res.positivity_ok = check.check_positivity();
    res.unitarity_ok = check.check_unitarity();
    res.perturbativity_ok = check.check_perturbativity();
    res.stability_ok = check.check_stability();
    
    res.theory_ok = res.positivity_ok && res.unitarity_ok && res.perturbativity_ok;
    res.triple_ok = res.theory_ok; // Conceptually triple_ok == theory_ok for now
    
    if (res.theory_ok) {
        DecayTable decays(model);
        res.width_bb = decays.get_gamma_hdd(2, 3, 3);
        res.width_tautau = decays.get_gamma_hll(2, 3, 3);
        res.width_WW = decays.get_gamma_hvv(2, 3);
        res.width_ZZ = decays.get_gamma_hvv(2, 2);
        res.width_gaga = decays.get_gamma_hgaga(2);
        res.width_Zga = decays.get_gamma_hZga(2);
        res.width_gg = decays.get_gamma_hgg(2);
        res.width_hh = decays.get_gamma_hhh(2, 1, 1);
        res.total_width = decays.get_gammatot_h(2);
        if (std::isfinite(res.total_width) && res.total_width > 0.0) {
            res.br_gaga = res.width_gaga / res.total_width;
        }
    }
    
    return res;
}
