/**
 * @file PhysLam1Scan.cpp
 * @brief Parameter scan over m_phi and lambda1 using reconstructed m12^2 + set_param_phys
 *
 * Single-threaded implementation that reconstructs m12^2 from lambda1 and
 * calls THDM::set_param_phys. This keeps compatibility with 2HDMC builds where
 * set_param_phys_lam1 is not present in the linked library.
 *
 * CLI contract: 12 positional args
 *   mphi_min mphi_max N_mphi lam1_min lam1_max N_lam1 mA sin_ba tan_beta lambda6 lambda7 output.csv
 *
 * Output: CSV with 29 columns, atomic write via .tmp rename
 * Stdout markers: "Total Attempts: <int>" and "TRIPLE_OK_POINTS <int>"
 * 
 *  ./dihiggs/app/PhysLam1Scan 
 *  OMP_NUM_THREADS=12 ./dihiggs/app/PhysLam1Scan 130 290 15 0 6 100 300 1 10000 0.0001 0 test.csv
 */

#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <chrono>
#include <cstdio>  // For std::rename

using namespace std;
using namespace std::chrono;

constexpr double M_H = 125.0;  // Light Higgs mass (fixed)
constexpr double VEV = 246.0;

using Clock    = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;

static double reconstruct_m12_from_lam1(
    double mh,
    double mH,
    double sin_ba,
    double tan_beta,
    double lambda1,
    double lambda6,
    double lambda7
) {
    const double inv = 1.0 / std::sqrt(1.0 + tan_beta * tan_beta);
    const double cb = inv;
    const double sb = tan_beta * inv;
    const double cba = std::sqrt(std::max(0.0, 1.0 - sin_ba * sin_ba));
    const double ca = cb * cba + sb * sin_ba;
    const double sa = sb * cba - cb * sin_ba;

    const double term1 = (mH * mH * ca * ca + mh * mh * sa * sa) / (VEV * VEV * cb * cb);
    const double coeff = (VEV * VEV) * cb * cb / tan_beta;
    return (term1 - 1.5 * lambda6 * tan_beta + 0.5 * lambda7 * tan_beta * tan_beta * tan_beta - lambda1) * coeff;
}

void perform_param_scan_lam1(
    double mphi_min, double mphi_max, int N_mphi,
    double lam1_min, double lam1_max, int N_lam1,
    double mA_fixed, double sin_ba,   double tan_beta,
    double lambda6,  double lambda7,
    const string &output_file
) {
    // Write to temporary file first for atomic output
    string tmp_path = output_file + ".tmp";
    ofstream results_tmp(tmp_path);
    if (!results_tmp.is_open()) {
        cerr << "Error: cannot open " << tmp_path << endl;
        exit(2);
    }

    // CSV header - exact match with PhysScanWithFixings
    vector<string> columns = {
        "m_phi","mA","alpha","beta","lambda6","lambda7","m12",
        "sin_ba","tan_beta","positivity_ok","unitarity_ok","perturbativity_ok",
        "width_bb","width_tautau","width_WW","width_ZZ",
        "width_gaga","width_Zga","width_gg","width_hh",
        "total_width","br_gaga", "lam1", "computed_lam1",
        "lam2","computed_lam2","lam3","lam4","lam5"
    };
    write_csv_header(results_tmp, columns);

    // Calculate step sizes
    double step_mphi = (N_mphi > 1) ? (mphi_max - mphi_min)/(N_mphi - 1) : 0.0;
    double step_lam1 = (N_lam1 > 1) ? (lam1_max - lam1_min)/(N_lam1 - 1) : 0.0;

    // Counters for statistics
    long long total_attempts = 0;
    long long total_valid = 0;
    long long total_triple_ok = 0;

    auto start_time = Clock::now();
    auto last_report_time = start_time;

    const double REPORT_INTERVAL_SEC = 2.5;
    long long total_tasks = (long long)N_mphi * (long long)N_lam1;

    cerr << "step_lam1=" << step_lam1 << endl;
    cerr << "step_mphi=" << step_mphi << endl;
    cerr << "Total iterations to compute: " << total_tasks << endl;

    // Precalculate angles for reporting/logging
    
    // Angles for reporting/logging
    const double beta  = std::atan(tan_beta);
    const double alpha = beta - std::asin(sin_ba);

    // Single-threaded scan (no OpenMP)
    for(int i_phi = 0; i_phi < N_mphi; ++i_phi) {
        // Create model instance per outer loop
        THDM model;
        SM sm; 
        model.set_SM(sm);
        model.set_yukawas_type(1);  // Hard-coded yukawa type

        double m_phi = mphi_min + i_phi * step_mphi;

        for(int i_l1 = 0; i_l1 < N_lam1; ++i_l1) {
            ++total_attempts;

            double lambda1 = lam1_min + i_l1 * step_lam1;

            const double m12_2 = reconstruct_m12_from_lam1(
                M_H, m_phi, sin_ba, tan_beta, lambda1, lambda6, lambda7
            );
            if (!std::isfinite(m12_2)) {
                continue;
            }

            // Compatibility path: use set_param_phys with reconstructed m12^2.
            bool pset = model.set_param_phys(
                M_H, m_phi, mA_fixed, mA_fixed,
                sin_ba, lambda6, lambda7, m12_2, tan_beta
            );

            if (!pset) continue;

            // Extract generic basis parameters
            double lam1_g, lam2_g, lam3_g, lam4_g, lam5_g;
            double lam6_g, lam7_g, m12_2_g, tanb_g;
            model.get_param_gen(
                lam1_g, lam2_g, lam3_g, lam4_g, lam5_g,
                lam6_g, lam7_g, m12_2_g, tanb_g
            );

            // Check theoretical constraints
            Constraints check(model);
            bool pos  = check.check_positivity();
            bool uni  = check.check_unitarity();
            bool pert = check.check_perturbativity();

            // Count triple-OK points
            if (pos && uni && pert) {
                ++total_triple_ok;
            }

            // Compute decay widths
            DecayTable tab(model);
            double w_bb     = tab.get_gamma_hdd(2,3,3);
            double w_tautau = tab.get_gamma_hll(2,3,3);
            double w_WW     = tab.get_gamma_hvv(2,3);
            double w_ZZ     = tab.get_gamma_hvv(2,2);
            double w_gaga   = tab.get_gamma_hgaga(2);
            double w_Zga    = tab.get_gamma_hZga(2);
            double w_gg     = tab.get_gamma_hgg(2);
            double w_hh     = tab.get_gamma_hhh(2,1,1);
            double w_tot    = tab.get_gammatot_h(2);
            double br_gaga  = (w_tot > 1e-15) ? w_gaga / w_tot : 0.0;

            // Write CSV row with exact column order
            vector<double> row = {
                m_phi, mA_fixed, alpha, beta, lambda6, lambda7, m12_2_g,
                sin_ba, tan_beta,
                (pos  ? 1.0 : 0.0),
                (uni  ? 1.0 : 0.0),
                (pert ? 1.0 : 0.0),
                w_bb, w_tautau, w_WW, w_ZZ,
                w_gaga, w_Zga, w_gg, w_hh,
                w_tot, br_gaga,
                lambda1, lam1_g,  // requested vs computed lam1
                lam2_g, lam2_g,   // Note: no separate l2 calculation, use computed
                lam3_g, lam4_g, lam5_g
            };
            write_csv_row(results_tmp, row);

            ++total_valid;

            // Progress reporting to stderr
            auto now = Clock::now();
            Duration diff = now - last_report_time;
            
            if (diff.count() > REPORT_INTERVAL_SEC) {
                double total_elapsed = Duration(now - start_time).count();
                double pps = (total_elapsed > 0) ? (double)total_attempts / total_elapsed : 0.0;
                double vppm = (total_elapsed > 0) ? (double)total_valid / (total_elapsed / 60.0) : 0.0;
                double progress = (total_tasks > 0) ? (100.0 * total_attempts / (double)total_tasks) : 0.0;

                std::cerr << "\r" 
                          << "[ " << std::fixed << std::setprecision(1)
                          << progress << "% ] "
                          << "Speed: " << (long)pps << " pts/s | "
                          << "ValidCSV: " << (long)vppm << "/min | "
                          << "FoundCSV: " << total_valid << "/"
                          << total_attempts
                          << " | TripleOK: " << total_triple_ok
                          << "   " << std::flush;
                
                last_report_time = now;
            }
        } // i_l1
    } // i_phi

    // Close and flush temporary file
    results_tmp.close();
    if (!results_tmp.good()) {
        cerr << "Error: failed to close " << tmp_path << endl;
        exit(2);
    }

    // Atomic rename
    if (std::rename(tmp_path.c_str(), output_file.c_str()) != 0) {
        cerr << "Error: failed to rename " << tmp_path << " to " << output_file << endl;
        exit(2);
    }

    auto total_time = Duration(Clock::now() - start_time).count();

    std::cout << "\n\n--- Scan Finalizado ---" << endl;
    std::cout << "Total Attempts: " << total_attempts << endl;
    std::cout << "Total CSV Rows: " << total_valid << endl;
    std::cout << "Total Triple-OK Points: " << total_triple_ok << endl;
    std::cout << "Total Time: " << total_time << " s" << endl;
    if (total_time > 0) {
        std::cout << "Average Speed: "
                  << (long)(total_attempts / total_time) << " pts/s" << endl;
        std::cout << "Final Valid Rate (CSV): "
                  << (long)(total_valid / (total_time/60.0)) << " rows/min" << endl;
    }

    // TRIPLE_OK marker for bash/pipeline parsing
    std::cout << "TRIPLE_OK_POINTS " << total_triple_ok << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 13) {
        cerr << "Usage: " << argv[0]
             << " mphi_min mphi_max N_mphi lam1_min lam1_max N_lam1"
                " mA sin(b-a) tan(beta) lambda6 lambda7 output.csv\n";
        return 1;
    }

    double mphi_min = atof(argv[1]);
    double mphi_max = atof(argv[2]);
    int    N_mphi   = atoi(argv[3]);
    double lam1_min = atof(argv[4]);
    double lam1_max = atof(argv[5]);
    int    N_lam1   = atoi(argv[6]);
    double mA_fixed = atof(argv[7]);
    double sin_ba   = atof(argv[8]);
    double tan_beta = atof(argv[9]);
    double lambda6  = atof(argv[10]);
    double lambda7  = atof(argv[11]);
    string output   = argv[12];

    perform_param_scan_lam1(
        mphi_min, mphi_max, N_mphi,
        lam1_min, lam1_max, N_lam1,
        mA_fixed, sin_ba,   tan_beta,
        lambda6,  lambda7,
        output
    );
    return 0;
}
