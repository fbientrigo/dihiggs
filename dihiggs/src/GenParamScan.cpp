#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"  // CSV & config utilities
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <omp.h>  // OpenMP for parallelization

using namespace std;
using namespace std::chrono;

struct ParamSet {
    double lambda1;
    double lambda2;
    double lambda3;
    double lambda4;
    double lambda5;
    double m12;
    double beta;
};

// Compute branching ratio BR(H→γγ)
double computeBRGammaGamma(const ParamSet &p) {
    THDM model;
    SM sm;
    model.set_SM(sm);

    double tb = std::tan(p.beta);

    bool pset = model.set_param_gen(
        p.lambda1, p.lambda2, p.lambda3, p.lambda4, p.lambda5,
        0.0, 0.0, p.m12, tb
    );
    if (!pset) return -1.0;

    Constraints check(model);
    if (!check.check_positivity() || !check.check_unitarity() || !check.check_perturbativity()) {
        return -1.0;
    }

    DecayTable table(model);
    double w_hgaga = table.get_gamma_hgaga(2);
    double w_htot  = table.get_gammatot_h(2);

    if (std::isnan(w_hgaga) || std::isnan(w_htot) || w_htot <= 0.0) {
        return -1.0;
    }

    return w_hgaga / w_htot;
}

// Global best parameter storage
static ParamSet g_bestParams;
static double   g_bestBR = -1.0;


void perform_param_scan_withGD(const Config &cfg, const string &output_file) {
    ofstream results(output_file);
    if (!results.is_open()) {
        cerr << "Failed to open output file: " << output_file << endl;
        return;
    }

    vector<string> columns = {
        "lambda1", "lambda2", "lambda3", "lambda4", "lambda5",
        "m12_squared", "beta", "tan_beta",
        "positivity_ok", "unitarity_ok", "perturbativity_ok",
        "width_h2_bb", "width_h2_tautau", "width_h2_WW", "width_h2_ZZ", "width_hgaga",
        "total_width_h2", "branching_ratio_hgaga"
    };
    write_csv_header(results, columns);

    // Convert floating loops to integer steps
    int steps_lambda5 = (cfg.lambda5_max - cfg.lambda5_min) / cfg.step_lambda5 + 1;
    int steps_m12     = (cfg.m12_squared_max - cfg.m12_squared_min) / cfg.step_m12_squared + 1;
    int steps_beta    = (cfg.beta_max - cfg.beta_min) / cfg.step_beta + 1;

    double total_iterations = computeTotalIterations(cfg);
    double current_iteration = 0.0;
    auto start_time = high_resolution_clock::now();

    g_bestBR = -1.0;

    for (double lambda1 = cfg.lambda1_min; lambda1 <= cfg.lambda1_max; lambda1 += cfg.step_lambda1) {
        for (double lambda2 = cfg.lambda2_min; lambda2 <= cfg.lambda2_max; lambda2 += cfg.step_lambda2) {
            for (double lambda3 = cfg.lambda3_min; lambda3 <= cfg.lambda3_max; lambda3 += cfg.step_lambda3) {
                for (double lambda4 = cfg.lambda4_min; lambda4 <= cfg.lambda4_max; lambda4 += cfg.step_lambda4) {

                    #pragma omp parallel
                    {
                        double local_bestBR = -1.0;
                        ParamSet local_bestParams;
                        vector<vector<double>> thread_local_rows;  // Buffer for batch writing

                        #pragma omp for collapse(3) schedule(dynamic)
                        for (int i_lambda5 = 0; i_lambda5 < steps_lambda5; i_lambda5++) {
                            for (int i_m12 = 0; i_m12 < steps_m12; i_m12++) {
                                for (int i_beta = 0; i_beta < steps_beta; i_beta++) {

                                    // Compute the floating values
                                    double lambda5 = cfg.lambda5_min + i_lambda5 * cfg.step_lambda5;
                                    double m12     = cfg.m12_squared_min + i_m12 * cfg.step_m12_squared;
                                    double beta    = cfg.beta_min + i_beta * cfg.step_beta;

                                    THDM model;
                                    SM sm;
                                    model.set_SM(sm);

                                    bool pset = model.set_param_gen(
                                        lambda1, lambda2, lambda3, lambda4, lambda5,
                                        0.0, 0.0, m12, std::tan(beta)
                                    );
                                    if (!pset) continue;

                                    Constraints check(model);
                                    bool positivity_ok     = check.check_positivity();
                                    bool unitarity_ok      = check.check_unitarity();
                                    bool perturbativity_ok = check.check_perturbativity();

                                    DecayTable table(model);
                                    double w_h2_bb     = table.get_gamma_hdd(2, 3, 3);
                                    double w_h2_tautau = table.get_gamma_hll(2, 3, 3);
                                    double w_h2_WW     = table.get_gamma_hvv(2, 3);
                                    double w_h2_ZZ     = table.get_gamma_hvv(2, 2);
                                    double w_hgaga     = table.get_gamma_hgaga(2);
                                    double w_total_h2  = table.get_gammatot_h(2);

                                    double br_hgaga = (w_total_h2 > 1e-12) ? (w_hgaga / w_total_h2) : 0.0;

                                    vector<double> row = {
                                        lambda1, lambda2, lambda3, lambda4, lambda5,
                                        m12, beta, std::tan(beta),
                                        positivity_ok ? 1.0 : 0.0, unitarity_ok ? 1.0 : 0.0, perturbativity_ok ? 1.0 : 0.0,
                                        w_h2_bb, w_h2_tautau, w_h2_WW, w_h2_ZZ, w_hgaga,
                                        w_total_h2, br_hgaga
                                    };

                                    thread_local_rows.push_back(row);

                                    if (positivity_ok && unitarity_ok && w_total_h2 > 0.0 && !std::isnan(br_hgaga)) {
                                        if (br_hgaga > local_bestBR) {
                                            local_bestBR = br_hgaga;
                                            local_bestParams = {lambda1, lambda2, lambda3, lambda4, lambda5, m12, beta};
                                        }
                                    }

                                    #pragma omp critical
                                    {
                                        current_iteration += 1.0;
                                        auto now_time = high_resolution_clock::now();
                                        double elapsed = duration<double>(now_time - start_time).count();
                                        double progress = (current_iteration / total_iterations) * 100.0;
                                        double estimated_time_left = (elapsed / current_iteration) * (total_iterations - current_iteration);

                                        cout << fixed << setprecision(2)
                                             << "Progress: " << progress << "% | "
                                             << "Elapsed: " << elapsed / 60.0 << " min | "
                                             << "ETA: " << estimated_time_left / 60.0 << " min"
                                             << "\r" << flush;
                                    }
                                }
                            }
                        }

                        #pragma omp critical
                        {
                            for (const auto& row : thread_local_rows) {
                                write_csv_row(results, row);
                            }
                            thread_local_rows.clear();

                            if (local_bestBR > g_bestBR) {
                                g_bestBR = local_bestBR;
                                g_bestParams = local_bestParams;
                            }
                        }
                    }
                }
            }
        }
    }

    results.close();
    cout << "\nGrid search completed. Results saved to " << output_file << endl;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <config_file> <output_csv>\n";
        return 1;
    }

    string config_file = argv[1];
    string output_file = argv[2];

    try {
        Config cfg = readConfig(config_file);
        perform_param_scan_withGD(cfg, output_file);
    } catch(const exception &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
