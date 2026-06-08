#include "M2BatchEvaluator.hpp"
#include "SM.h"
#include <omp.h>
#include <gsl/gsl_errno.h>

std::vector<PointResult> evaluate_m2_batch(
    const BatchParams& params,
    const std::vector<double>& m2_candidates
) {
    gsl_set_error_handler_off();
    
    std::vector<PointResult> results(m2_candidates.size());
    
    // Set up standard model parameters that will be copied to local models
    SM sm;

    #pragma omp parallel
    {
        // Crucial: Create a strictly thread-local THDM model instance
        // to prevent OpenMP data races across matrix inversions and cache states.
        THDM local_model;
        local_model.set_SM(sm);

        #pragma omp for
        for (size_t i = 0; i < m2_candidates.size(); ++i) {
            results[i] = evaluate_m2_point(
                local_model,
                params.m_phi,
                params.mH,
                params.mA,
                params.mHp,
                params.sin_ba,
                params.tan_beta,
                params.lambda6,
                params.lambda7,
                m2_candidates[i]
            );
        }
    }

    return results;
}
