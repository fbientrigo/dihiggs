/**
 * @file PhysScanWithFixings.cpp
 * @brief Parallel m_phi/lambda1 scan with fixed physical parameters.
 *
 * IMPORTANT:
 * - Uses THDM::set_param_phys_lam1 (same lambda1 path used by validation/autoresearch)
 *   to avoid manual m12^2 reconstruction drift.
 * - Keeps CSV contract and stdout counters expected by orchestrators.
 */

#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"
#include "ReplaySafeOutput.hpp"

#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <vector>

using namespace std;
using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;

constexpr double M_H = 125.0;

void perform_param_scan_fixings(
    double mphi_min, double mphi_max, int N_mphi,
    double lam1_min, double lam1_max, int N_lam1,
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
        "m_phi","mA","alpha","beta","lambda6","lambda7","m12"
    };
    const auto replay_columns = replay_safe_output::metadata_column_names();
    columns.insert(columns.end(), replay_columns.begin(), replay_columns.end());
    const vector<string> trailing_columns = {
        "sin_ba","tan_beta","positivity_ok","unitarity_ok","perturbativity_ok",
        "width_bb","width_tautau","width_WW","width_ZZ",
        "width_gaga","width_Zga","width_gg","width_hh",
        "total_width","br_gaga", "lam1", "computed_lam1",
        "lam2","computed_lam2","lam3","lam4","lam5"
    };
    columns.insert(columns.end(), trailing_columns.begin(), trailing_columns.end());
    write_csv_header(results, columns);

    const double step_mphi = (N_mphi > 1) ? (mphi_max - mphi_min) / (N_mphi - 1) : 0.0;
    const double step_lam1 = (N_lam1 > 1) ? (lam1_max - lam1_min) / (N_lam1 - 1) : 0.0;

    const long long total_tasks = static_cast<long long>(N_mphi) * static_cast<long long>(N_lam1);
    std::atomic<long long> global_attempts(0);
    std::atomic<long long> global_valid(0);
    std::atomic<long long> global_triple_ok(0);

    auto start_time = Clock::now();
    auto last_report_time = start_time;

    const double REPORT_INTERVAL_SEC = 2.5;
    const int BATCH_SIZE = 200;

    const double beta = std::atan(tan_beta);
    const double alpha = beta - std::asin(sin_ba);

    cout << "step_lam1=" << step_lam1 << endl;
    cout << "step_mphi=" << step_mphi << endl;
    cout << "Total iterations to compute: " << total_tasks << endl;

    #pragma omp parallel
    {
        long long local_attempts = 0;
        long long local_valid = 0;
        long long local_triple_ok = 0;

        std::ostringstream local_buffer;
        // High precision is critical: m12 truncation causes large replay drift in CalcPhys.
        replay_safe_output::configure_scientific_17(local_buffer);

        #pragma omp for schedule(dynamic)
        for (int i_phi = 0; i_phi < N_mphi; ++i_phi) {
            const double m_phi = mphi_min + i_phi * step_mphi;

            THDM model;
            SM sm;
            model.set_SM(sm);
            model.set_yukawas_type(1);

            for (int i_l1 = 0; i_l1 < N_lam1; ++i_l1) {
                ++local_attempts;
                const double lambda1 = lam1_min + i_l1 * step_lam1;

                // THDM API order is: sba, lambda1, lambda6, lambda7, tan_beta.
                const bool pset_lam1 = model.set_param_phys_lam1(
                    M_H, m_phi, mA_fixed, mA_fixed,
                    sin_ba, lambda1, lambda6, lambda7, tan_beta
                );
                if (!pset_lam1) {
                    continue;
                }

                // Canonicalize through set_param_phys using the generated m12^2 so output
                // is directly replayable in CalcPhys (ground truth path).
                double lam1_g, lam2_g, lam3_g, lam4_g, lam5_g;
                double lam6_g, lam7_g, m12_2_g, tanb_g;
                model.get_param_gen(
                    lam1_g, lam2_g, lam3_g, lam4_g, lam5_g,
                    lam6_g, lam7_g, m12_2_g, tanb_g
                );
                const double m12_2_used = m12_2_g;

                THDM canonical_model;
                canonical_model.set_SM(sm);
                const bool pset_canonical = canonical_model.set_param_phys(
                    M_H, m_phi, mA_fixed, mA_fixed,
                    sin_ba, lambda6, lambda7, m12_2_used, tan_beta
                );
                if (!pset_canonical) {
                    continue;
                }
                canonical_model.set_yukawas_type(1);

                // Re-read generic params from canonical model.
                canonical_model.get_param_gen(
                    lam1_g, lam2_g, lam3_g, lam4_g, lam5_g,
                    lam6_g, lam7_g, m12_2_g, tanb_g
                );
                const double m12_2_gen_after_set = m12_2_g;
                const replay_safe_output::Metadata replay_meta =
                    replay_safe_output::make_metadata(m12_2_used, m12_2_gen_after_set, canonical_model.get_yukawas_type());

                Constraints check(canonical_model);
                const bool pos = check.check_positivity();
                const bool uni = check.check_unitarity();
                const bool pert = check.check_perturbativity();

                if (pos && uni && pert) {
                    ++local_triple_ok;
                }

                DecayTable tab(canonical_model);
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
                    << m_phi << "," << mA_fixed << "," << alpha << "," << beta << ","
                    << lambda6 << "," << lambda7 << "," << m12_2_used << ","
                    ;
                replay_safe_output::append_metadata_csv(local_buffer, replay_meta);
                local_buffer
                    << ","
                    << sin_ba << "," << tan_beta << ","
                    << (pos ? 1.0 : 0.0) << ","
                    << (uni ? 1.0 : 0.0) << ","
                    << (pert ? 1.0 : 0.0) << ","
                    << w_bb << "," << w_tautau << "," << w_WW << "," << w_ZZ << ","
                    << w_gaga << "," << w_Zga << "," << w_gg << "," << w_hh << ","
                    << w_tot << "," << br_gaga << ","
                    << lambda1 << "," << lam1_g << ","
                    << lam2_g << "," << lam2_g << ","
                    << lam3_g << "," << lam4_g << "," << lam5_g
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
        cerr << "Uso: " << argv[0]
             << " mphi_min mphi_max N_mphi lam1_min lam1_max N_lam1"
                " mA sin(b-a) tan(beta) lambda6 lambda7 output.csv\n";
        return 1;
    }

    double mphi_min = atof(argv[1]);
    double mphi_max = atof(argv[2]);
    int N_mphi = atoi(argv[3]);
    double lam1_min = atof(argv[4]);
    double lam1_max = atof(argv[5]);
    int N_lam1 = atoi(argv[6]);
    double mA_fixed = atof(argv[7]);
    double sin_ba = atof(argv[8]);
    double tan_beta = atof(argv[9]);
    double lambda6 = atof(argv[10]);
    double lambda7 = atof(argv[11]);
    string output = argv[12];

    perform_param_scan_fixings(
        mphi_min, mphi_max, N_mphi,
        lam1_min, lam1_max, N_lam1,
        mA_fixed, sin_ba, tan_beta,
        lambda6, lambda7,
        output
    );
    return 0;
}
