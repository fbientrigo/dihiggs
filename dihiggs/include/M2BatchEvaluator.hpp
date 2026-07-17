#pragma once

#include "M2PointEvaluator.hpp"
#include <vector>

struct BatchParams {
    double mh;
    double mH;
    double mA;
    double mHp;
    double sin_ba;
    double tan_beta;
    double lambda6;
    double lambda7;
};

std::vector<PointResult> evaluate_m2_batch(
    const BatchParams& params,
    const std::vector<double>& m2_candidates
);
