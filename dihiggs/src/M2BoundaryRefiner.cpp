#include "M2BoundaryRefiner.hpp"
#include "SM.h"
#include <gsl/gsl_errno.h>

bool evaluate_is_valid(const BatchParams& params, double m2) {
    gsl_set_error_handler_off(); // Prevent GSL from aborting the program on domain errors
    SM sm;
    THDM local_model;
    local_model.set_SM(sm);
    
    PointResult res = evaluate_m2_point(
        local_model, params.m_phi, params.mH, params.mA, params.mHp,
        params.sin_ba, params.tan_beta, params.lambda6, params.lambda7, m2
    );
    return res.triple_ok;
}

ValidInterval refine_interval_boundaries(
    const BatchParams& params,
    const ValidInterval& interval,
    double tolerance_m2,
    int max_iterations
) {
    ValidInterval refined = interval;
    if (interval.is_point) return refined;

    // 1. Refine the Lower Boundary
    double low_valid = interval.M2_low;
    double low_invalid = interval.M2_outer_low;
    
    // Only bisect if we have a proper invalid outer bound
    if (low_invalid < low_valid) {
        for (int i = 0; i < max_iterations; ++i) {
            if ((low_valid - low_invalid) <= tolerance_m2) break;
            
            double mid = (low_valid + low_invalid) / 2.0;
            if (evaluate_is_valid(params, mid)) {
                low_valid = mid;
            } else {
                low_invalid = mid;
            }
        }
        refined.M2_low = low_valid;
        refined.M2_outer_low = low_invalid;
    }

    // 2. Refine the Higher Boundary
    double high_valid = interval.M2_high;
    double high_invalid = interval.M2_outer_high;
    
    // Only bisect if we have a proper invalid outer bound
    if (high_invalid > high_valid) {
        for (int i = 0; i < max_iterations; ++i) {
            if ((high_invalid - high_valid) <= tolerance_m2) break;
            
            double mid = (high_valid + high_invalid) / 2.0;
            if (evaluate_is_valid(params, mid)) {
                high_valid = mid;
            } else {
                high_invalid = mid;
            }
        }
        refined.M2_high = high_valid;
        refined.M2_outer_high = high_invalid;
    }
    
    refined.M2_width = refined.M2_high - refined.M2_low;
    refined.M2_center = (refined.M2_high + refined.M2_low) / 2.0;

    return refined;
}
