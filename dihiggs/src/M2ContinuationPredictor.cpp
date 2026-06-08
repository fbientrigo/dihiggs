#include "M2ContinuationPredictor.hpp"
#include <algorithm>
#include <cmath>

PredictedBounds predict_constant(
    const ValidInterval& prior_interval,
    const PredictorConfig& config,
    bool is_bootstrap
) {
    PredictedBounds b;
    b.expected_center = prior_interval.M2_center;
    
    double pad = is_bootstrap ? config.bootstrap_pad_abs : config.prior_pad_abs;
    double dynamic_pad = config.min_pad_to_bandwidth_ratio * prior_interval.M2_width;
    pad = std::max(pad, dynamic_pad);
    
    b.search_low = prior_interval.M2_low - pad;
    b.search_high = prior_interval.M2_high + pad;
    
    b.search_low = std::max(b.search_low, config.global_m2_min);
    b.search_high = std::min(b.search_high, config.global_m2_max);
    return b;
}

PredictedBounds predict_linear(
    const ValidInterval& prior_interval,
    const ValidInterval& prior_prior_interval,
    double target_m_phi,
    const PredictorConfig& config
) {
    PredictedBounds b;
    double dmphi = prior_interval.m_phi - prior_prior_interval.m_phi;
    if (std::abs(dmphi) < 1e-5) {
        return predict_constant(prior_interval, config, false);
    }
    
    double dM2_low_dmphi = (prior_interval.M2_low - prior_prior_interval.M2_low) / dmphi;
    double dM2_high_dmphi = (prior_interval.M2_high - prior_prior_interval.M2_high) / dmphi;
    
    double step = target_m_phi - prior_interval.m_phi;
    
    double expected_low = prior_interval.M2_low + dM2_low_dmphi * step;
    double expected_high = prior_interval.M2_high + dM2_high_dmphi * step;
    
    b.expected_center = (expected_low + expected_high) / 2.0;
    
    double pad = config.prior_pad_abs;
    double dynamic_pad = config.min_pad_to_bandwidth_ratio * prior_interval.M2_width;
    pad = std::max(pad, dynamic_pad);
    
    b.search_low = expected_low - pad;
    b.search_high = expected_high + pad;
    
    b.search_low = std::max(b.search_low, config.global_m2_min);
    b.search_high = std::min(b.search_high, config.global_m2_max);
    
    return b;
}

std::vector<double> generate_search_grid(const PredictedBounds& bounds, int n_points) {
    std::vector<double> grid;
    if (n_points <= 0) return grid;
    if (n_points == 1) {
        grid.push_back(bounds.expected_center);
        return grid;
    }
    
    grid.reserve(n_points);
    double step = (bounds.search_high - bounds.search_low) / (n_points - 1);
    for (int i = 0; i < n_points; ++i) {
        grid.push_back(bounds.search_low + i * step);
    }
    return grid;
}
