/**
 * @file Phys_M2BoundaryScan.cpp
 * @brief Parallel physical-parameter scan over m_phi and M2.
 *
 * M2 is the externally controlled scan coordinate:
 *
 *   M2 = m12_sq / (sin(beta) * cos(beta))
 *
 * The engine derives m12_sq from M2 in long double precision, casts once to
 * double, and passes that value to THDM::set_param_phys.
 *
 * HARDENING NOTES (v2):
 *  - A fresh THDM model and SM object are created per grid point inside the
 *    inner M2 loop to eliminate any possible state carry-over between points.
 *  - Strict CLI/input validation rejects invalid inputs with nonzero exit code.
 */

#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;

constexpr double M_H_LIGHT = 125.0;
constexpr const char* CONTRACT_ID = "phys_m2_boundary_v1";

// ---------------------------------------------------------------------------
// Strict input validation: returns false and prints a clear message to stderr
// if any constraint is violated.
// ---------------------------------------------------------------------------
static bool validate_inputs(
    double mphi_min, double mphi_max, int N_mphi,
    double M2_min,   double M2_max,   int N_M2,
    double /*mA_fixed*/, double sin_ba, double tan_beta,
    double lambda6, double lambda7,
    double sin_beta_cos_beta
) {
    bool ok = true;

    auto fail = [&](const char* msg) {
        std::cerr << "[ERROR] " << msg << "\n";
        ok = false;
    };

    // Grid counts
    if (N_mphi < 1) fail("N_mphi < 1: must have at least one m_phi point.");
    if (N_M2   < 1) fail("N_M2 < 1: must have at least one M2 point.");

    // Range ordering
    if (mphi_max < mphi_min) fail("mphi_max < mphi_min.");
    if (M2_max   < M2_min)   fail("M2_max < M2_min.");

    // Physical bounds
    if (tan_beta <= 0.0)        fail("tan_beta <= 0: must be positive.");
    if (std::abs(sin_ba) > 1.0) fail("|sin(b-a)| > 1: unphysical.");

    // Finiteness checks
    auto chk_finite = [&](double v, const char* name) {
        if (!std::isfinite(v)) {
            std::string msg = std::string(name) + " is not finite.";
            std::cerr << "[ERROR] " << msg << "\n";
            ok = false;
        }
    };
    chk_finite(mphi_min, "mphi_min");
    chk_finite(mphi_max, "mphi_max");
    chk_finite(M2_min,   "M2_min");
    chk_finite(M2_max,   "M2_max");
    chk_finite(sin_ba,   "sin_ba");
    chk_finite(tan_beta, "tan_beta");
    chk_finite(lambda6,  "lambda6");
    chk_finite(lambda7,  "lambda7");

    // sin_beta * cos_beta must be nonzero and finite (M2 -> m12_sq conversion)
    if (!std::isfinite(sin_beta_cos_beta) || sin_beta_cos_beta == 0.0) {
        fail("sin_beta_cos_beta is zero or non-finite: M2->m12_sq conversion undefined.");
    }

    return ok;
}

void perform_param_scan_m2_boundary(
    double mphi_min, double mphi_max, int N_mphi,
    double M2_min, double M2_max, int N_M2,
    double mA_fixed, double sin_ba, double tan_beta,
    double lambda6, double lambda7,
    const string &output_file
) {
    ofstream results(output_file);
    if (!results.is_open()) {
        cerr << "No pude abrir: " << output_file << endl;
        return;
    }

    vector<string> columns = {
        "contract_id",
        "m_phi","mH","mA","mHp",
        "alpha","beta","sin_ba","tan_beta",
        "lambda6","lambda7",
        "M2_input","m12_sq_input",
        "sin_beta","cos_beta","sin_beta_cos_beta",
        "positivity_ok","unitarity_ok","perturbativity_ok","stability_ok",
        "lam1_out","lam2_out","lam3_out","lam4_out","lam5_out",
        "lam6_out","lam7_out","m12_sq_out","tanb_out",
        "width_bb","width_tautau","width_WW","width_ZZ",
        "width_gaga","width_Zga","width_gg","width_hh",
        "total_width","br_gaga"
    };
    write_csv_header(results, columns);

    const double step_mphi = (N_mphi > 1) ? (mphi_max - mphi_min) / (N_mphi - 1) : 0.0;
    const double step_M2 = (N_M2 > 1) ? (M2_max - M2_min) / (N_M2 - 1) : 0.0;

    const long long total_tasks = static_cast<long long>(N_mphi) * static_cast<long long>(N_M2);
    std::atomic<long long> global_attempts(0);
    std::atomic<long long> global_valid(0);
    std::atomic<long long> global_triple_ok(0);

    auto start_time = Clock::now();
    auto last_report_time = start_time;

    const double REPORT_INTERVAL_SEC = 2.5;
    const int BATCH_SIZE = 200;

    // Precompute trigonometric quantities in long double precision.
    const long double tan_beta_ld = static_cast<long double>(tan_beta);
    const long double beta_ld = std::atan(tan_beta_ld);
    const long double alpha_ld = beta_ld - std::asin(static_cast<long double>(sin_ba));
    const long double sin_beta_ld = std::sin(beta_ld);
    const long double cos_beta_ld = std::cos(beta_ld);
    const long double sin_beta_cos_beta_ld = sin_beta_ld * cos_beta_ld;

    const double beta = static_cast<double>(beta_ld);
    const double alpha = static_cast<double>(alpha_ld);
    const double sin_beta = static_cast<double>(sin_beta_ld);
    const double cos_beta = static_cast<double>(cos_beta_ld);
    const double sin_beta_cos_beta = static_cast<double>(sin_beta_cos_beta_ld);

    cout << "step_M2=" << step_M2 << endl;
    cout << "step_mphi=" << step_mphi << endl;
    cout << "Total iterations to compute: " << total_tasks << endl;

    #pragma omp parallel
    {
        long long local_attempts = 0;
        long long local_valid = 0;
        long long local_triple_ok = 0;

        std::ostringstream local_buffer;
        local_buffer.setf(std::ios::scientific);
        local_buffer << std::setprecision(17);

        #pragma omp for schedule(dynamic)
        for (int i_phi = 0; i_phi < N_mphi; ++i_phi) {
            const double m_phi = mphi_min + i_phi * step_mphi;
            const double mH = m_phi;
            const double mHp = mA_fixed;

            for (int i_M2 = 0; i_M2 < N_M2; ++i_M2) {
                ++local_attempts;

                // --- Fresh THDM + SM per grid point ---
                // A new model is instantiated here to prevent any state
                // carry-over from the previous M2 evaluation.
                THDM model;
                SM sm;
                model.set_SM(sm);
                model.set_yukawas_type(1);

                const double M2_input = M2_min + i_M2 * step_M2;
                const long double M2_input_ld = static_cast<long double>(M2_input);
                const long double m12_sq_input_ld =
                    M2_input_ld * sin_beta_ld * cos_beta_ld;
                const double m12_sq_input = static_cast<double>(m12_sq_input_ld);

                const bool pset = model.set_param_phys(
                    M_H_LIGHT, mH, mA_fixed, mHp,
                    sin_ba, lambda6, lambda7, m12_sq_input, tan_beta
                );
                if (!pset) {
                    continue;
                }
                model.set_yukawas_type(1);

                double lam1_g, lam2_g, lam3_g, lam4_g, lam5_g;
                double lam6_g, lam7_g, m12_sq_g, tanb_g;
                model.get_param_gen(
                    lam1_g, lam2_g, lam3_g, lam4_g, lam5_g,
                    lam6_g, lam7_g, m12_sq_g, tanb_g
                );

                Constraints check(model);
                const bool pos = check.check_positivity();
                const bool uni = check.check_unitarity();
                const bool pert = check.check_perturbativity();
                const bool stab = check.check_stability();

                if (pos && uni && pert) {
                    ++local_triple_ok;
                }

                DecayTable tab(model);
                const double w_bb = tab.get_gamma_hdd(2, 3, 3);
                const double w_tautau = tab.get_gamma_hll(2, 3, 3);
                const double w_WW = tab.get_gamma_hvv(2, 3);
                const double w_ZZ = tab.get_gamma_hvv(2, 2);
                const double w_gaga = tab.get_gamma_hgaga(2);
                const double w_Zga = tab.get_gamma_hZga(2);
                const double w_gg = tab.get_gamma_hgg(2);
                const double w_hh = tab.get_gamma_hhh(2, 1, 1);
                const double w_tot = tab.get_gammatot_h(2);
                const double br_gaga = (w_tot > 1e-15) ? (w_gaga / w_tot) : 0.0;

                local_buffer
                    << CONTRACT_ID << ","
                    << m_phi << "," << mH << "," << mA_fixed << "," << mHp << ","
                    << alpha << "," << beta << "," << sin_ba << "," << tan_beta << ","
                    << lambda6 << "," << lambda7 << ","
                    << M2_input << "," << m12_sq_input << ","
                    << sin_beta << "," << cos_beta << "," << sin_beta_cos_beta << ","
                    << (pos ? 1.0 : 0.0) << ","
                    << (uni ? 1.0 : 0.0) << ","
                    << (pert ? 1.0 : 0.0) << ","
                    << (stab ? 1.0 : 0.0) << ","
                    << lam1_g << "," << lam2_g << "," << lam3_g << ","
                    << lam4_g << "," << lam5_g << ","
                    << lam6_g << "," << lam7_g << "," << m12_sq_g << "," << tanb_g << ","
                    << w_bb << "," << w_tautau << "," << w_WW << "," << w_ZZ << ","
                    << w_gaga << "," << w_Zga << "," << w_gg << "," << w_hh << ","
                    << w_tot << "," << br_gaga
                    << "\n";

                ++local_valid;

                if (local_valid >= BATCH_SIZE) {
                    #pragma omp critical
                    {
                        results << local_buffer.str();
                        global_valid += local_valid;
                        global_attempts += local_attempts;
                        global_triple_ok += local_triple_ok;

                        auto now = Clock::now();
                        Duration diff = now - last_report_time;
                        if (diff.count() > REPORT_INTERVAL_SEC) {
                            const long long current_attempts = global_attempts.load();
                            const long long current_valid = global_valid.load();
                            const long long current_triple = global_triple_ok.load();
                            const double total_elapsed = Duration(now - start_time).count();

                            const double pps = (total_elapsed > 0)
                                ? static_cast<double>(current_attempts) / total_elapsed
                                : 0.0;
                            const double vppm = (total_elapsed > 0)
                                ? static_cast<double>(current_valid) / (total_elapsed / 60.0)
                                : 0.0;
                            const double progress = (total_tasks > 0)
                                ? (100.0 * current_attempts / static_cast<double>(total_tasks))
                                : 0.0;

                            std::cerr << "\r"
                                      << "[ " << std::fixed << std::setprecision(1)
                                      << progress << "% ] "
                                      << "Speed: " << static_cast<long>(pps) << " pts/s | "
                                      << "ValidCSV: " << static_cast<long>(vppm) << "/min | "
                                      << "FoundCSV: " << current_valid << "/" << current_attempts
                                      << " | TripleOK: " << current_triple
                                      << "   " << std::flush;

                            last_report_time = now;
                        }
                    }

                    local_buffer.str("");
                    local_buffer.clear();
                    local_valid = 0;
                    local_attempts = 0;
                    local_triple_ok = 0;
                }
            }
        }

        #pragma omp critical
        {
            if (local_valid > 0) {
                results << local_buffer.str();
                global_valid += local_valid;
            }
            if (local_attempts > 0) {
                global_attempts += local_attempts;
            }
            if (local_triple_ok > 0) {
                global_triple_ok += local_triple_ok;
            }
        }
    }

    const double total_time = Duration(Clock::now() - start_time).count();
    const long long final_attempts = global_attempts.load();
    const long long final_valid = global_valid.load();
    const long long final_triple_ok = global_triple_ok.load();

    std::cout << "\n\n--- Scan Finalizado ---" << endl;
    std::cout << "Total Attempts: " << final_attempts << endl;
    std::cout << "Total CSV Rows: " << final_valid << endl;
    std::cout << "Total Triple-OK Points: " << final_triple_ok << endl;
    std::cout << "Total Time: " << total_time << " s" << endl;
    if (total_time > 0) {
        std::cout << "Average Speed: "
                  << static_cast<long>(final_attempts / total_time) << " pts/s" << endl;
        std::cout << "Final Valid Rate (CSV): "
                  << static_cast<long>(final_valid / (total_time / 60.0)) << " rows/min" << endl;
    }

    std::cout << "TRIPLE_OK_POINTS " << final_triple_ok << std::endl;
    results.close();
}

int main(int argc, char* argv[]) {
    if (argc != 13) {
        cerr << "Usage: " << argv[0]
             << " mphi_min mphi_max N_mphi M2_min M2_max N_M2"
                " mA sin(b-a) tan(beta) lambda6 lambda7 output.csv\n";
        return 1;
    }

    const double mphi_min  = atof(argv[1]);
    const double mphi_max  = atof(argv[2]);
    const int    N_mphi    = atoi(argv[3]);
    const double M2_min    = atof(argv[4]);
    const double M2_max    = atof(argv[5]);
    const int    N_M2      = atoi(argv[6]);
    const double mA_fixed  = atof(argv[7]);
    const double sin_ba    = atof(argv[8]);
    const double tan_beta  = atof(argv[9]);
    const double lambda6   = atof(argv[10]);
    const double lambda7   = atof(argv[11]);
    const string output    = argv[12];

    // Precompute sin_beta_cos_beta here so validate_inputs can check it.
    const long double tan_beta_ld = static_cast<long double>(tan_beta);
    const long double beta_ld = std::atan(tan_beta_ld);
    const long double sin_beta_ld = std::sin(beta_ld);
    const long double cos_beta_ld = std::cos(beta_ld);
    const double sin_beta_cos_beta =
        static_cast<double>(sin_beta_ld * cos_beta_ld);

    if (!validate_inputs(
            mphi_min, mphi_max, N_mphi,
            M2_min, M2_max, N_M2,
            mA_fixed, sin_ba, tan_beta,
            lambda6, lambda7,
            sin_beta_cos_beta)) {
        cerr << "Input validation failed. Aborting.\n";
        return 2;
    }

    perform_param_scan_m2_boundary(
        mphi_min, mphi_max, N_mphi,
        M2_min, M2_max, N_M2,
        mA_fixed, sin_ba, tan_beta,
        lambda6, lambda7,
        output
    );
    return 0;
}
