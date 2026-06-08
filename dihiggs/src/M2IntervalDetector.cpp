#include "M2IntervalDetector.hpp"

std::vector<ValidInterval> detect_intervals(const std::vector<PointResult>& results) {
    std::vector<ValidInterval> intervals;
    if (results.empty()) return intervals;

    bool in_interval = false;
    size_t start_idx = 0;

    for (size_t i = 0; i < results.size(); ++i) {
        if (results[i].triple_ok) {
            if (!in_interval) {
                start_idx = i;
                in_interval = true;
            }
        } else {
            if (in_interval) {
                ValidInterval iv;
                iv.m_phi = results[start_idx].m_phi;
                iv.M2_low = results[start_idx].M2_input;
                iv.M2_high = results[i - 1].M2_input;
                iv.M2_center = (iv.M2_low + iv.M2_high) / 2.0;
                iv.M2_width = iv.M2_high - iv.M2_low;
                iv.start_index = start_idx;
                iv.end_index = i - 1;
                iv.is_point = (start_idx == i - 1);
                
                iv.M2_outer_low = (start_idx > 0) ? results[start_idx - 1].M2_input : iv.M2_low;
                iv.M2_outer_high = results[i].M2_input; // i is the first invalid point after interval
                
                intervals.push_back(iv);
                in_interval = false;
            }
        }
    }

    if (in_interval) {
        size_t end_idx = results.size() - 1;
        ValidInterval iv;
        iv.m_phi = results[start_idx].m_phi;
        iv.M2_low = results[start_idx].M2_input;
        iv.M2_high = results[end_idx].M2_input;
        iv.M2_center = (iv.M2_low + iv.M2_high) / 2.0;
        iv.M2_width = iv.M2_high - iv.M2_low;
        iv.start_index = start_idx;
        iv.end_index = end_idx;
        iv.is_point = (start_idx == end_idx);
        
        iv.M2_outer_low = (start_idx > 0) ? results[start_idx - 1].M2_input : iv.M2_low;
        iv.M2_outer_high = iv.M2_high; // no outer high available
        
        intervals.push_back(iv);
    }

    return intervals;
}
