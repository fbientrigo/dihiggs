#pragma once

#include "M2IntervalDetector.hpp"
#include <vector>

struct PredictorConfig {
    double prior_pad_abs;
    double bootstrap_pad_abs;
    double min_pad_to_bandwidth_ratio;
    double global_m2_min;
    double global_m2_max;
};

struct PredictedBounds {
    double expected_center;
    double search_low;
    double search_high;
};

PredictedBounds predict_constant(
    const ValidInterval& prior_interval,
    const PredictorConfig& config,
    bool is_bootstrap
);

PredictedBounds predict_linear(
    const ValidInterval& prior_interval,
    const ValidInterval& prior_prior_interval,
    double target_m_phi,
    const PredictorConfig& config
);

std::vector<double> generate_search_grid(const PredictedBounds& bounds, int n_points);
