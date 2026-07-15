#pragma once

#include "M2BatchEvaluator.hpp"
#include "M2IntervalDetector.hpp"

ValidInterval refine_interval_boundaries(
    const BatchParams& params,
    const ValidInterval& interval,
    double tolerance_m2,
    int max_iterations
);
