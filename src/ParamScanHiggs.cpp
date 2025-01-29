#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"  // Your CSV & config utils
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <iomanip> // if needed for std::setprecision

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

double computeBRGammaGamma(const ParamSet &p) {
    THDM model;
    SM sm;
    model.set_SM(sm);

    double tb = std::tan(p.beta);

    bool pset = model.set_param_gen(
        p.lambda1, p.lambda2, p.lambda3, p.lambda4, p.lambda5,
        0.0, 0.0,
        p.m12,
        tb
    );
    if (!pset) return -1.0;

    Constraints check(model);
    bool positivity_ok     = check.check_positivity();
    bool unitarity_ok      = check.check_unitarity();
    bool perturbativity_ok = check.check_perturbativity();

    if (!positivity_ok || !unitarity_ok || !perturbativity_ok) {
        return -1.0;
    }

    DecayTable table(model);
    double w_hgaga = table.get_gamma_hgaga(2);
    double w_htot  = table.get_gammatot_h(2);

    if (std::isnan(w_hgaga) || std::isnan(w_htot) || w_htot <= 0.0) {
        return -1.0;
    }

    double br = w_hgaga / w_htot;
    if (std::isnan(br) || br < 0.0) {
        return -1.0;
    }
    return br;
}

ParamSet computeGradient(const ParamSet &p, double eps = 1e-4) {
    ParamSet grad = {0,0,0,0,0,0,0};

    // helper to modify a parameter in a copy
    auto brWithModified = [&](const ParamSet &origP, double ParamSet::* mem, double newVal) -> double {
        ParamSet tmp = origP;
        tmp.*mem = newVal;
        return computeBRGammaGamma(tmp);
    };

    auto finiteDiff = [&](double ParamSet::* mem) {
        double x0 = p.*mem;
        double fplus  = brWithModified(p, mem, x0 + eps);
        double fminus = brWithModified(p, mem, x0 - eps);

        if (fplus < 0.0 || fminus < 0.0) {
            return 0.0;
        }
        return (fplus - fminus) / (2.0 * eps);
    };

    grad.lambda1 = finiteDiff(&ParamSet::lambda1);
    grad.lambda2 = finiteDiff(&ParamSet::lambda2);
    grad.lambda3 = finiteDiff(&ParamSet::lambda3);
    grad.lambda4 = finiteDiff(&ParamSet::lambda4);
    grad.lambda5 = finiteDiff(&ParamSet::lambda5);
    grad.m12     = finiteDiff(&ParamSet::m12);
    grad.beta    = finiteDiff(&ParamSet::beta);

    return grad;
}

ParamSet gradientDescent(const ParamSet &start, double &bestBR, double alpha=0.005, int maxIters=70) {
    ParamSet current = start;
    bestBR           = computeBRGammaGamma(current);
    if (bestBR < 0) {
        return current;
    }

    for (int i = 0; i < maxIters; i++) {
        ParamSet g = computeGradient(current, 1e-4);

        ParamSet trial = {
            current.lambda1 + alpha * g.lambda1,
            current.lambda2 + alpha * g.lambda2,
            current.lambda3 + alpha * g.lambda3,
            current.lambda4 + alpha * g.lambda4,
            current.lambda5 + alpha * g.lambda5,
            current.m12     + alpha * g.m12,
            current.beta    + alpha * g.beta
        };

        double brNext = computeBRGammaGamma(trial);
        if (brNext < 0) {
            alpha *= 0.5;
        } else {
            current = trial;
            if (brNext > bestBR) {
                bestBR = brNext;
            }
        }
    }
    return current;
}

// globals
static ParamSet g_bestParams;
static double   g_bestBR = -1.0;

static void perform_param_scan_withGD(const Config &cfg, const string &output_file) {
    ofstream results(output_file);
    if (!results.is_open()) {
        cerr << "Failed to open output file: " << output_file << endl;
        return;
    }

    // Define columns, but we'll store only doubles in the row
    vector<string> columns = {
        "lambda1","lambda2","lambda3","lambda4","lambda5",
        "m12_squared","beta","tan_beta",
        "positivity_ok","unitarity_ok","perturbativity_ok",
        "width_h2_bb","width_h2_tautau","width_h2_WW","width_h2_ZZ","width_hgaga",
        "total_width_h2","branching_ratio_hgaga"
    };
    // This function now presumably writes out column headers (as text).
    write_csv_header(results, columns);

    double total_iterations = computeTotalIterations(cfg);
    double current_iteration = 0.0;
    auto start_time = high_resolution_clock::now();

    g_bestBR = -1.0;

    for (double lambda1 = cfg.lambda1_min; lambda1 <= cfg.lambda1_max; lambda1 += cfg.step_lambda1) {
        for (double lambda2 = cfg.lambda2_min; lambda2 <= cfg.lambda2_max; lambda2 += cfg.step_lambda2) {
            for (double lambda3 = cfg.lambda3_min; lambda3 <= cfg.lambda3_max; lambda3 += cfg.step_lambda3) {
                for (double lambda4 = cfg.lambda4_min; lambda4 <= cfg.lambda4_max; lambda4 += cfg.step_lambda4) {
                    for (double lambda5 = cfg.lambda5_min; lambda5 <= cfg.lambda5_max; lambda5 += cfg.step_lambda5) {
                        for (double m12 = cfg.m12_squared_min; m12 <= cfg.m12_squared_max; m12 += cfg.step_m12_squared) {
                            for (double beta = cfg.beta_min; beta <= cfg.beta_max; beta += cfg.step_beta) {

                                current_iteration += 1.0;
                                auto now_time = high_resolution_clock::now();
                                double elapsed = duration<double>(now_time - start_time).count();
                                double progress = (current_iteration / total_iterations);
                                print_progress(progress, elapsed, total_iterations, current_iteration);

                                // Build model
                                THDM model;
                                SM sm;
                                model.set_SM(sm);

                                bool pset = model.set_param_gen(
                                    lambda1, lambda2, lambda3, lambda4, lambda5,
                                    0.0, 0.0,
                                    m12,
                                    std::tan(beta)
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

                                // Convert bool -> double (1.0/0.0)
                                double pos_val = positivity_ok     ? 1.0 : 0.0;
                                double uni_val = unitarity_ok      ? 1.0 : 0.0;
                                double per_val = perturbativity_ok ? 1.0 : 0.0;

                                double br_hgaga = 0.0;
                                if (w_total_h2 > 1e-12) {
                                    br_hgaga = w_hgaga / w_total_h2;
                                }

                                // Build the row as doubles
                                vector<double> row;
                                row.reserve(columns.size()); // for efficiency

                                // scanned params
                                row.push_back(lambda1);
                                row.push_back(lambda2);
                                row.push_back(lambda3);
                                row.push_back(lambda4);
                                row.push_back(lambda5);
                                row.push_back(m12);
                                row.push_back(beta);
                                row.push_back(std::tan(beta));

                                // validity checks as 0/1
                                row.push_back(pos_val);
                                row.push_back(uni_val);
                                row.push_back(per_val);

                                // partial widths & total
                                row.push_back(w_h2_bb);
                                row.push_back(w_h2_tautau);
                                row.push_back(w_h2_WW);
                                row.push_back(w_h2_ZZ);
                                row.push_back(w_hgaga);

                                row.push_back(w_total_h2);
                                row.push_back(br_hgaga);

                                // Now write the row of doubles
                                write_csv_row(results, row);

                                // track best
                                // if (positivity_ok && unitarity_ok && w_total_h2 > 0.0 && !std::isnan(br_hgaga)) {
                                if (w_total_h2 > 0.0 && !std::isnan(br_hgaga)) {
                                    if (br_hgaga > g_bestBR) {
                                        g_bestBR = br_hgaga;
                                        g_bestParams = {lambda1, lambda2, lambda3, lambda4, lambda5, m12, beta};
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    results.close();
    cout << "\nGrid search completed. Results saved to " << output_file << endl;

    if (g_bestBR <= 0.0) {
        cout << "No valid points found in grid. Gradient Descent skipped.\n";
        return;
    }

    cout << "Best BR found in grid: " << g_bestBR << " at beta=" << g_bestParams.beta << "\n"
         << "Starting gradient descent...\n";

    double alpha      = 0.01;
    int    maxIters   = 50;
    double bestAfterGD = g_bestBR;

    ParamSet finalGD = gradientDescent(g_bestParams, bestAfterGD, alpha, maxIters);

    cout << "Gradient Descent done.\n"
         << "   Initial BR: " << g_bestBR << "\n"
         << "   Final   BR: " << bestAfterGD << "\n"
         << "Params after GD:\n"
         << "   lambda1=" << finalGD.lambda1 
         << ", lambda2=" << finalGD.lambda2
         << ", lambda3=" << finalGD.lambda3
         << ", lambda4=" << finalGD.lambda4
         << ", lambda5=" << finalGD.lambda5
         << ", m12="     << finalGD.m12
         << ", beta="    << finalGD.beta
         << endl;
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
