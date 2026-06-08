#include "M2FallbackPolicy.hpp"
#include <iostream>
#include <cmath>

FallbackResult execute_fallback_policy(
    const BatchParams& params,
    double current_step_mphi,
    const PredictedBounds& failed_bounds,
    const FallbackConfig& config,
    const ValidInterval& prior_interval,
    const ValidInterval& prior_prior_interval
) {
    FallbackResult res;
    res.success = false;
    res.step_halved = false;
    res.new_step_mphi = current_step_mphi;
    res.intervals.clear();

    // 1. Fallback Local Dense Scan
    // Expand the failed bounds and scan with high density.
    PredictedBounds dense_bounds = failed_bounds;
    dense_bounds.search_low -= config.fallback_dense_pad;
    dense_bounds.search_high += config.fallback_dense_pad;
    
    std::cout << "[Fallback] Triggering Local Dense Scan at m_phi=" << params.m_phi << std::endl;
    auto dense_grid = generate_search_grid(dense_bounds, config.fallback_dense_count);
    auto dense_results = evaluate_m2_batch(params, dense_grid);
    auto dense_intervals = detect_intervals(dense_results);
    
    if (!dense_intervals.empty()) {
        std::cout << "[Fallback] Local Dense Scan succeeded!" << std::endl;
        res.success = true;
        res.intervals = dense_intervals;
        return res;
    }

    // 2. Dynamic Mass Step Halving
    // If local dense failed, curvature is likely too sharp. Halve the mass step.
    if (config.allow_mass_step_halving) {
        double halved = current_step_mphi / 2.0;
        if (halved >= config.min_mphi_step) {
            std::cout << "[Fallback] Triggering Dynamic Mass Step Halving: " 
                      << current_step_mphi << " -> " << halved << std::endl;
            res.step_halved = true;
            res.new_step_mphi = halved;
            return res; // Main loop handles the rewinding
        } else {
            std::cout << "[Fallback] Mass step halving reached minimum limit (" << config.min_mphi_step << ")." << std::endl;
        }
    }

    // 3. Fallback Full Mass-Slice Grid
    if (config.fallback_full_enable) {
        std::cout << "[Fallback] Triggering Full Mass-Slice Grid at m_phi=" << params.m_phi << std::endl;
        PredictedBounds full_bounds;
        full_bounds.search_low = config.full_scan_min;
        full_bounds.search_high = config.full_scan_max;
        full_bounds.expected_center = (config.full_scan_min + config.full_scan_max) / 2.0;
        
        auto full_grid = generate_search_grid(full_bounds, config.full_scan_count);
        auto full_results = evaluate_m2_batch(params, full_grid);
        auto full_intervals = detect_intervals(full_results);
        
        if (!full_intervals.empty()) {
            std::cout << "[Fallback] Full Mass-Slice Grid succeeded!" << std::endl;
            res.success = true;
            res.intervals = full_intervals;
            return res;
        }
    }

    return res;
}
