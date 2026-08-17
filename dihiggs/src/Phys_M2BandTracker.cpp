#include "M2BatchEvaluator.hpp"
#include "M2IntervalDetector.hpp"
#include "M2ContinuationPredictor.hpp"
#include "M2BoundaryRefiner.hpp"
#include "M2FallbackPolicy.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

using namespace std;
using namespace std::chrono;

// --- Argument parsing ---
map<string, string> parse_args(int argc, char* argv[]) {
    map<string, string> args;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg.substr(0, 2) == "--") {
            auto pos = arg.find('=');
            if (pos != string::npos) {
                args[arg.substr(2, pos - 2)] = arg.substr(pos + 1);
            } else if (i + 1 < argc && string(argv[i+1]).substr(0, 2) != "--") {
                args[arg.substr(2)] = argv[++i];
            } else {
                args[arg.substr(2)] = "true";
            }
        }
    }
    return args;
}

double get_arg(const map<string, string>& args, const string& key, double def) {
    return args.count(key) ? stod(args.at(key)) : def;
}
int get_arg_i(const map<string, string>& args, const string& key, int def) {
    return args.count(key) ? stoi(args.at(key)) : def;
}
string get_arg_s(const map<string, string>& args, const string& key, const string& def) {
    return args.count(key) ? args.at(key) : def;
}
bool get_arg_b(const map<string, string>& args, const string& key, bool def) {
    return args.count(key) ? (args.at(key) == "true" || args.at(key) == "1") : def;
}

// Write the header for points.csv
void write_points_header(ofstream& out) {
    out << "m_phi,M2_input,m12_sq_input,m12_sq_out,lam1_out,lam2_out,lam3_out,lam4_out,lam5_out,lam6_out,lam7_out,positivity_ok,unitarity_ok,perturbativity_ok,stability_ok,theory_ok,triple_ok,construction_ok,yukawa_type_installed\n";
}

// Write a batch of points to points.csv
void write_points_batch(ofstream& out, const vector<PointResult>& results) {
    for (const auto& r : results) {
        out << std::fixed << std::setprecision(6)
            << r.m_phi << "," << r.M2_input << "," << r.m12_sq_input << "," << r.m12_sq_out << ","
            << r.lam1_out << "," << r.lam2_out << "," << r.lam3_out << "," << r.lam4_out << "," << r.lam5_out << "," << r.lam6_out << "," << r.lam7_out << ","
            << r.positivity_ok << "," << r.unitarity_ok << "," << r.perturbativity_ok << "," << r.stability_ok << ","
            << r.theory_ok << "," << r.triple_ok << "," << r.construction_ok << "," << r.yukawa_type_installed << "\n";
    }
}

// Write the header for intervals.csv
void write_intervals_header(ofstream& out) {
    out << "m_phi,M2_low,M2_high,M2_center,M2_width,band_width_estimate,prediction_error_M2\n";
}

int main(int argc, char* argv[]) {
    auto args = parse_args(argc, argv);
    
    // Physical Parameters
    double mphi_min = get_arg(args, "mphi-min", 140.0);
    double mphi_max = get_arg(args, "mphi-max", 500.0);
    double step_mphi = get_arg(args, "mphi-step", 1.0);
    
    double mA = get_arg(args, "ma", 500.0);
    double mh = get_arg(args, "mh", 125.20);
    double mHp = get_arg(args, "mhp", mA);
    double sin_ba = get_arg(args, "sin-ba", 1.0);
    double tan_beta = get_arg(args, "tan-beta", 50.0);
    double lambda6 = get_arg(args, "lam6", 0.001);
    double lambda7 = get_arg(args, "lam7", 0.0);
    
    // Predictor Config
    PredictorConfig config;
    config.prior_pad_abs = get_arg(args, "prior-pad-abs", 200.0);
    config.bootstrap_pad_abs = get_arg(args, "bootstrap-pad-abs", 5000.0);
    config.min_pad_to_bandwidth_ratio = get_arg(args, "min-pad-to-bandwidth-ratio", 3.0);
    config.global_m2_min = get_arg(args, "m2-min", -1e7);
    config.global_m2_max = get_arg(args, "m2-max", 1e7);
    
    // Fallback Config
    FallbackConfig fb_config;
    fb_config.fallback_dense_count = get_arg_i(args, "fallback-dense-count", 500);
    fb_config.fallback_dense_pad = get_arg(args, "fallback-dense-pad", 500.0);
    fb_config.fallback_full_enable = get_arg_b(args, "fallback-full-enable", true);
    fb_config.full_scan_min = get_arg(args, "m2-min", 10000);
    fb_config.full_scan_max = get_arg(args, "m2-max", 300000);
    fb_config.full_scan_count = get_arg_i(args, "fallback-full-count", 1000);
    fb_config.allow_mass_step_halving = get_arg_b(args, "allow-mass-step-halving", true);
    fb_config.min_mphi_step = get_arg(args, "min-mphi-step", 0.125);
    
    // Seed Config
    double seed_m2_center = get_arg(args, "seed-M2-center", 19500);
    double seed_m2_halfwidth = get_arg(args, "seed-M2-halfwidth", 500);
    int seed_n_m2 = get_arg_i(args, "seed-n-M2", 200);
    
    // Outputs
    string out_points = get_arg_s(args, "out-points", "../scripts/out/points.csv");
    string out_intervals = get_arg_s(args, "out-intervals", "../scripts/out/intervals.csv");
    string out_summary = get_arg_s(args, "out-summary", "../scripts/out/mass_slice_summary.jsonl");
    string out_meta = get_arg_s(args, "out-meta", "../scripts/out/tracker_state_final.json");

    ofstream file_points(out_points);
    ofstream file_intervals(out_intervals);
    ofstream file_summary(out_summary);
    
    if(!file_points.is_open() || !file_intervals.is_open() || !file_summary.is_open()) {
        cerr << "Failed to open output files!" << endl;
        return 1;
    }
    
    write_points_header(file_points);
    write_intervals_header(file_intervals);
    
    ValidInterval prior_interval;
    ValidInterval prior_prior_interval;
    bool has_prior = false;
    bool has_prior_prior = false;
    
    auto t_start = high_resolution_clock::now();
    int total_evaluations = 0;
    int total_triple_ok = 0;
    
    double current_step_mphi = step_mphi;
    double m_phi = mphi_min;

    while (m_phi <= mphi_max) {
        BatchParams params = {mh, m_phi, mA, mHp, sin_ba, tan_beta, lambda6, lambda7};
        
        PredictedBounds bounds;
        if (!has_prior) {
            bounds.search_low = seed_m2_center - seed_m2_halfwidth;
            bounds.search_high = seed_m2_center + seed_m2_halfwidth;
            bounds.expected_center = seed_m2_center;
        } else if (!has_prior_prior) {
            bounds = predict_constant(prior_interval, config, true);
        } else {
            bounds = predict_linear(prior_interval, prior_prior_interval, m_phi, config);
        }
        
        int points_per_slice = has_prior ? get_arg_i(args, "local-m2-count", 200) : seed_n_m2;
        auto grid = generate_search_grid(bounds, points_per_slice);
        auto results = evaluate_m2_batch(params, grid);
        
        write_points_batch(file_points, results);
        
        total_evaluations += points_per_slice;
        int slice_triple_ok = 0;
        for (const auto& r : results) if (r.triple_ok) slice_triple_ok++;
        total_triple_ok += slice_triple_ok;
        
        auto intervals = detect_intervals(results);
        int fallbacks = 0;
        
        if (intervals.empty()) {
            fallbacks++;
            FallbackResult fb_res = execute_fallback_policy(params, current_step_mphi, bounds, fb_config, prior_interval, prior_prior_interval);
            
            if (fb_res.step_halved) {
                current_step_mphi = fb_res.new_step_mphi;
                m_phi = prior_interval.m_phi + current_step_mphi;
                continue;
            } else if (fb_res.success) {
                intervals = fb_res.intervals;
            } else {
                break;
            }
        }
        
        ValidInterval best_iv = select_interval_nearest(intervals, bounds.expected_center);
        
        // Bisection Refinement
        double tol = get_arg(args, "edge-tol", 0.1);
        int max_iter = get_arg_i(args, "max-edge-iter", 15);
        best_iv = refine_interval_boundaries(params, best_iv, tol, max_iter);

        file_intervals << std::fixed << std::setprecision(6)
            << best_iv.m_phi << "," << best_iv.M2_low << "," << best_iv.M2_high << "," << best_iv.M2_center << ","
            << best_iv.M2_width << "," << best_iv.M2_width << "," << (best_iv.M2_center - bounds.expected_center) << "\n";
            
        file_summary << "{\"m_phi\": " << m_phi << ", \"intervals\": " << intervals.size() 
                     << ", \"evals\": " << points_per_slice << "}\n";
        
        cout << "m_phi = " << m_phi << ", center = " << best_iv.M2_center << ", width = " << best_iv.M2_width << endl;
        
        prior_prior_interval = prior_interval;
        has_prior_prior = has_prior;
        prior_interval = best_iv;
        has_prior = true;
        
        if (current_step_mphi < step_mphi) {
            current_step_mphi *= 2.0;
            if (current_step_mphi > step_mphi) current_step_mphi = step_mphi;
        }

        m_phi += current_step_mphi;
    }
    
    auto t_end = high_resolution_clock::now();
    
    ofstream file_meta(out_meta);
    file_meta << "{\n"
              << "  \"total_evaluations\": " << total_evaluations << ",\n"
              << "  \"global_triple_ok\": " << total_triple_ok << ",\n"
              << "  \"time_ms\": " << duration_cast<milliseconds>(t_end - t_start).count() << "\n"
              << "}\n";
              
    return 0;
}
