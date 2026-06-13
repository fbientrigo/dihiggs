/*
 * CalcLambda1ScanFixings.cpp
 * ===========================
 *
 * OpenMP-parallel "Stage 1" fast scanner over the 2HDM Type-I dihiggs
 * parameter space, using the self-contained physics core in
 * CalcLambda1Core.{h,c} (lifted from CalcLambda1.c, no 2HDMC/GSL deps).
 *
 * Fixed (compile-time) constants: mh=125.0, yukawa_type=1, lambda7=0,
 * lambda1=0, sin(beta-alpha)=1.0, m_Hp=m_A.
 *
 * Scanned: tan_beta (log-spaced), m_A=m_Hp (linear), lambda_6 (log-spaced),
 * m_H (linear in [130, m_A], computed per (tan_beta,m_A,lambda_6) shard).
 *
 * One (tan_beta, m_A, lambda_6) "fixing" tuple = one bronze shard CSV file
 * = all m_H points for one BR/ctau plot.  Each OpenMP thread writes its own
 * shard file(s); no cross-thread I/O contention.
 *
 * Usage:
 *   CalcLambda1ScanFixings \
 *       <tan_beta_min> <tan_beta_max> <N_tan_beta> \
 *       <mA_min> <mA_max> <N_mA> \
 *       <lambda6_min> <lambda6_max> <N_lambda6> \
 *       <N_mH> \
 *       <output_dir>
 */

#include "CalcLambda1Core.h"
#include "PathTags.hpp"

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;

// ---- compile-time fixed constants (per spec) -----------------------------
constexpr double MH_FIXED = 125.0;
constexpr int YUKAWA_TYPE_FIXED = 1;
constexpr double LAMBDA7_FIXED = 0.0;
constexpr double LAMBDA1_FIXED = 0.0;
constexpr double SIN_BA_FIXED = 1.0;
constexpr double MH_LOWER_BOUND = 130.0;

namespace {

double logspace_value(double lo, double hi, int n, int i) {
    if (n <= 1) return lo;
    const double log_lo = std::log10(lo);
    const double log_hi = std::log10(hi);
    const double t = log_lo + (log_hi - log_lo) * static_cast<double>(i) / (n - 1);
    return std::pow(10.0, t);
}

double linspace_value(double lo, double hi, int n, int i) {
    if (n <= 1) return lo;
    return lo + (hi - lo) * static_cast<double>(i) / (n - 1);
}

const std::vector<std::string>& bronze_columns() {
    static const std::vector<std::string> cols = {
        "tan_beta", "m_A", "m_Hp", "lambda6", "lambda7", "lambda1_input", "sin_ba", "mh_input", "yukawa_type",
        "m_H_input",
        "m12sq_recon",
        "lambda1_recon", "lambda2_recon", "lambda3_recon", "lambda4_recon", "lambda5_recon", "lambda6_recon", "lambda7_recon",
        "alpha_recon", "beta_recon",
        "mh_check", "mH_check", "mA_check", "mHp_check", "sba_check",
        "delta_mh", "delta_mH", "delta_mA", "delta_mHp", "delta_sba",
        "masses_positive", "perturbativity_ok", "unitarity_ok", "triple_ok_fast",
        "max_abs_lambda", "max_abs_unitarity_eigenvalue",
        "chris_width_bb", "chris_width_cc", "chris_width_tautau", "chris_width_gg", "chris_width_gaga", "chris_width_Zga",
        "chris_total_width", "chris_br_bb", "chris_br_cc", "chris_br_tautau", "chris_br_gg", "chris_br_gaga", "chris_br_Zga", "chris_br_loop",
        "chris_ctau_m"
    };
    return cols;
}

// One row = one m_H point evaluated for a fixed (tan_beta, m_A, lambda_6).
void write_row(std::ostream& out, double tan_beta, double m_A, double lambda6, double m_H) {
    InputPoint in;
    in.mh = MH_FIXED;
    in.mH = m_H;
    in.mA = m_A;
    in.mHp = m_A;          // m_Hp = m_A per spec
    in.sba = SIN_BA_FIXED;
    in.lambda6 = lambda6;
    in.lambda7 = LAMBDA7_FIXED;
    in.lambda1 = LAMBDA1_FIXED;
    in.tanb = tan_beta;
    in.yukawa_type = YUKAWA_TYPE_FIXED;

    LambdaSet L = reconstruct_lambdas(in);
    MassCheck M = masses_from_lambdas(in, L);
    WidthResult W = compute_widths(in);
    ConstraintResultFast C = check_constraints_fast(in, L);

    out << tan_beta << "," << m_A << "," << in.mHp << "," << lambda6 << "," << LAMBDA7_FIXED << ","
        << LAMBDA1_FIXED << "," << SIN_BA_FIXED << "," << MH_FIXED << "," << YUKAWA_TYPE_FIXED << ","
        << m_H << ","
        << L.m12sq << ","
        << L.lambda1 << "," << L.lambda2 << "," << L.lambda3 << "," << L.lambda4 << "," << L.lambda5 << "," << L.lambda6 << "," << L.lambda7 << ","
        << L.alpha << "," << L.beta << ","
        << M.mh << "," << M.mH << "," << M.mA << "," << M.mHp << "," << M.sba << ","
        << (M.mh - in.mh) << "," << (M.mH - in.mH) << "," << (M.mA - in.mA) << "," << (M.mHp - in.mHp) << "," << (M.sba - in.sba) << ","
        << C.masses_positive << "," << C.perturbativity_ok << "," << C.unitarity_ok << "," << C.triple_ok_fast << ","
        << C.max_abs_lambda << "," << C.max_abs_unitarity_eigenvalue << ","
        << W.G_bb << "," << W.G_cc << "," << W.G_tautau << "," << W.G_gg << "," << W.G_gaga << "," << W.G_Zga << ","
        << W.G_total << "," << W.BR_bb << "," << W.BR_cc << "," << W.BR_tautau << "," << W.BR_gg << "," << W.BR_gaga << "," << W.BR_Zga << "," << W.BR_loop << ","
        << W.ctau_m
        << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 12) {
        std::cerr << "Usage: " << argv[0]
                  << " tan_beta_min tan_beta_max N_tan_beta"
                     " mA_min mA_max N_mA"
                     " lambda6_min lambda6_max N_lambda6"
                     " N_mH"
                     " output_dir\n";
        return 1;
    }

    const double tb_min = std::atof(argv[1]);
    const double tb_max = std::atof(argv[2]);
    const int N_tb = std::atoi(argv[3]);
    const double mA_min = std::atof(argv[4]);
    const double mA_max = std::atof(argv[5]);
    const int N_mA = std::atoi(argv[6]);
    const double l6_min = std::atof(argv[7]);
    const double l6_max = std::atof(argv[8]);
    const int N_l6 = std::atoi(argv[9]);
    const int N_mH = std::atoi(argv[10]);
    const std::string output_dir = argv[11];

    const fs::path shards_root = fs::path(output_dir) / "bronze" / "points_shards";

    const long long total_shards =
        static_cast<long long>(N_tb) * static_cast<long long>(N_mA) * static_cast<long long>(N_l6);

    std::cout << "Total fixing tuples (shards): " << total_shards << std::endl;
    std::cout << "N_mH per shard: " << N_mH << std::endl;

    std::atomic<long long> global_attempts(0);
    std::atomic<long long> global_triple_ok_fast(0);
    std::atomic<long long> global_shards_written(0);
    std::atomic<long long> global_shards_skipped(0);

    auto start_time = Clock::now();
    auto last_report_time = start_time;
    const double REPORT_INTERVAL_SEC = 2.5;

    #pragma omp parallel
    {
        long long local_attempts = 0;
        long long local_triple_ok_fast = 0;

        #pragma omp for schedule(dynamic)
        for (long long idx = 0; idx < total_shards; ++idx) {
            const long long i_l6 = idx % N_l6;
            const long long i_mA = (idx / N_l6) % N_mA;
            const long long i_tb = idx / (static_cast<long long>(N_l6) * N_mA);

            const double tan_beta = logspace_value(tb_min, tb_max, N_tb, static_cast<int>(i_tb));
            const double m_A = linspace_value(mA_min, mA_max, N_mA, static_cast<int>(i_mA));
            const double lambda6 = logspace_value(l6_min, l6_max, N_l6, static_cast<int>(i_l6));

            if (m_A <= MH_LOWER_BOUND) {
                #pragma omp critical
                {
                    std::cerr << "WARNING: skipping m_A=" << m_A
                              << " (<= " << MH_LOWER_BOUND << "), no valid m_H range\n";
                }
                ++global_shards_skipped;
                continue;
            }

            const std::string mA_tag = path_tags::format_float_tag(m_A, 1);
            const std::string l6_tag = path_tags::format_float_tag_sci(lambda6, 6);
            const std::string tb_tag = path_tags::format_tanbeta_tag(tan_beta);

            const std::string fixed_dir_name =
                "fixed_sinba=" + path_tags::format_float_tag(SIN_BA_FIXED, 4) +
                "_l6=" + path_tags::format_float_tag(0.0, 7) +
                "_l7=" + path_tags::format_float_tag(LAMBDA7_FIXED, 4) +
                "_l1=" + path_tags::format_float_tag(LAMBDA1_FIXED, 4) +
                "_mA=" + mA_tag;

            const fs::path shard_dir = shards_root / fixed_dir_name / ("tb_" + tb_tag);
            const fs::path shard_path = shard_dir / ("shard_l6=" + l6_tag + "_mA=" + mA_tag + "_tb=" + tb_tag + ".csv");

            std::error_code ec;
            fs::create_directories(shard_dir, ec);
            if (ec) {
                #pragma omp critical
                {
                    std::cerr << "ERROR: cannot create directory " << shard_dir << ": " << ec.message() << "\n";
                }
                continue;
            }

            std::ofstream shard_file(shard_path);
            if (!shard_file.is_open()) {
                #pragma omp critical
                {
                    std::cerr << "ERROR: cannot open shard file " << shard_path << "\n";
                }
                continue;
            }
            shard_file.setf(std::ios::scientific);
            shard_file << std::setprecision(17);

            std::vector<std::string> header = bronze_columns();
            for (size_t c = 0; c < header.size(); ++c) {
                shard_file << header[c] << (c + 1 < header.size() ? "," : "\n");
            }

            for (int i_mH = 0; i_mH < N_mH; ++i_mH) {
                const double m_H = linspace_value(MH_LOWER_BOUND, m_A, N_mH, i_mH);
                write_row(shard_file, tan_beta, m_A, lambda6, m_H);

                ++local_attempts;
                // re-derive triple_ok_fast for the global counter without
                // recomputing the full point: cheapest is to recompute via
                // check_constraints_fast on the same inputs.
                InputPoint in;
                in.mh = MH_FIXED; in.mH = m_H; in.mA = m_A; in.mHp = m_A; in.sba = SIN_BA_FIXED;
                in.lambda6 = lambda6; in.lambda7 = LAMBDA7_FIXED; in.lambda1 = LAMBDA1_FIXED;
                in.tanb = tan_beta; in.yukawa_type = YUKAWA_TYPE_FIXED;
                LambdaSet L = reconstruct_lambdas(in);
                ConstraintResultFast C = check_constraints_fast(in, L);
                if (C.triple_ok_fast) ++local_triple_ok_fast;
            }

            shard_file.close();
            ++global_shards_written;

            #pragma omp critical
            {
                global_attempts += local_attempts;
                global_triple_ok_fast += local_triple_ok_fast;
                local_attempts = 0;
                local_triple_ok_fast = 0;

                auto now = Clock::now();
                Duration diff = now - last_report_time;
                if (diff.count() > REPORT_INTERVAL_SEC) {
                    const long long shards_done = global_shards_written.load() + global_shards_skipped.load();
                    const double progress = (total_shards > 0)
                        ? (100.0 * static_cast<double>(shards_done) / static_cast<double>(total_shards))
                        : 0.0;
                    std::cerr << "\r[ " << std::fixed << std::setprecision(1) << progress << "% ] "
                              << "Shards: " << shards_done << "/" << total_shards
                              << " | Attempts: " << global_attempts.load()
                              << " | TripleOK: " << global_triple_ok_fast.load()
                              << "   " << std::flush;
                    last_report_time = now;
                }
            }
        }

        #pragma omp critical
        {
            if (local_attempts > 0) global_attempts += local_attempts;
            if (local_triple_ok_fast > 0) global_triple_ok_fast += local_triple_ok_fast;
        }
    }

    const double total_time = Duration(Clock::now() - start_time).count();

    std::cout << "\n\n--- Stage 1 Scan Finished ---" << std::endl;
    std::cout << "Total Attempts: " << global_attempts.load() << std::endl;
    std::cout << "Shards Written: " << global_shards_written.load() << std::endl;
    std::cout << "Shards Skipped: " << global_shards_skipped.load() << std::endl;
    std::cout << "Total Time: " << total_time << " s" << std::endl;
    std::cout << "TRIPLE_OK_POINTS " << global_triple_ok_fast.load() << std::endl;

    return 0;
}
