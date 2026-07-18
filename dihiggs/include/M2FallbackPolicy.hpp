#pragma once

#include "M2ContinuationPredictor.hpp"
#include "M2BatchEvaluator.hpp"
#include "M2IntervalDetector.hpp"

struct FallbackConfig {
    int fallback_dense_count;
    double fallback_dense_pad;
    bool fallback_full_enable;
    double full_scan_min;
    double full_scan_max;
    int full_scan_count;
    bool allow_mass_step_halving;
    double min_mphi_step;
};

struct FallbackResult {
    bool success;
    std::vector<ValidInterval> intervals;
    bool step_halved;
    double new_step_mphi;
};

FallbackResult execute_fallback_policy(
    const BatchParams& params,
    double current_step_mphi,
    const PredictedBounds& failed_bounds,
    const FallbackConfig& config,
    const ValidInterval& prior_interval,
    const ValidInterval& prior_prior_interval
);
