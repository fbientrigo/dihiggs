#pragma once

#include "M2PointEvaluator.hpp"
#include <cstddef>
#include <vector>

struct ValidInterval {
    double m_phi;
    double M2_low;
    double M2_high;
    double M2_center;
    double M2_width;
    double M2_outer_low;
    double M2_outer_high;
    std::size_t start_index;
    std::size_t end_index;
    bool is_point;
};

std::vector<ValidInterval> detect_intervals(const std::vector<PointResult>& results);
ValidInterval select_interval_nearest(const std::vector<ValidInterval>& intervals, double predicted_center);
