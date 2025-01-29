#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"  // Header común

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;

static void perform_param_scan(const Config &cfg, const string &output_file) {

    ofstream results(output_file);
    if (!results.is_open()) {
        cerr << "Failed to open output file: " << output_file << endl;
        return;
    }

    // Columnas CSV
    vector<string> columns = {
        "lambda1", "lambda2", "lambda3", "lambda4", "lambda5",
        "m12_squared", "beta", "tan_beta",
        "width_hgaga", "total_width_h2", "branching_ratio_hgaga"
    };
    write_csv_header(results, columns);

    // Para mostrar progreso
    double total_iterations = computeTotalIterations(cfg);
    double current_iteration = 0.0;

    auto start_time = high_resolution_clock::now();

    // Iterar en los rangos
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

                                // Chequeo de positividad
                                if (!check_positivity(lambda1, lambda2, lambda3, lambda4, lambda5)) {
                                    continue;
                                }

                                // Inicializar modelo 2HDMC
                                THDM model;
                                SM sm;
                                model.set_SM(sm);

                                bool pset = model.set_param_gen(
                                    lambda1, lambda2, lambda3, lambda4, lambda5,
                                    0.0, 0.0,  // v1, v2 en base gen? (depende tu uso)
                                    m12,       // m12^2
                                    tan(beta)  // tan_beta
                                );
                                if (!pset) {
                                    continue;
                                }

                                Constraints check(model);
                                if (!check.check_positivity() || !check.check_unitarity() || !check.check_perturbativity()) {
                                    continue;
                                }

                                DecayTable table(model);
                                double width_hgaga = table.get_gamma_hgaga(2); // H2 -> gamma gamma
                                double total_w_h2  = table.get_gammatot_h(2);

                                if (std::isnan(width_hgaga) || std::isnan(total_w_h2) || total_w_h2 <= 0.0) {
                                    continue;
                                }

                                double br_hgaga = width_hgaga / total_w_h2;

                                // Guardar resultados
                                vector<double> row = {
                                    lambda1, lambda2, lambda3, lambda4, lambda5,
                                    m12, beta, tan(beta),
                                    width_hgaga, total_w_h2, br_hgaga
                                };
                                write_csv_row(results, row);
                            }
                        }
                    }
                }
            }
        }
    }

    results.close();
    cout << endl << "Parameter scan completed. Results saved to " << output_file << endl;
}

int main(int argc, char *argv[]) {

    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <config_file> <output_csv>\n";
        return 1;
    }

    string config_file = argv[1];
    string output_file = argv[2];

    try {
        // Leer config
        Config cfg = readConfig(config_file);
        // Lanzar escaneo
        perform_param_scan(cfg, output_file);
    }
    catch(const exception &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}
