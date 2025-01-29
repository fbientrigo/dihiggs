#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"  // Incluimos nuestro header común

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;

// Función que hace el escaneo
void perform_param_scan(const string &output_file, double step) {
    ofstream results(output_file);
    if (!results.is_open()) {
        cerr << "Failed to open output file." << endl;
        return;
    }

    // Columnas CSV (puedes seguir usando tu vector<string>)
    vector<string> columns = {
        "lambda1", "lambda2", "lambda3", "lambda4", "lambda5",
        "m12_squared", "beta", "tan_beta",
        "width_hgaga", "total_width_h2", "branching_ratio_hgaga"
    };
    write_csv_header(results, columns); // Llamada a la librería

    // Rangos...
    double infty_neg = -50.0, infty_pos = 100.0;
    double beta_min = 0.01, beta_max = M_PI / 2;
    double beta_step = 0.5;
    double lambda_min = 0.01, lambda_max = 16.0 * M_PI;

    // Calculamos total_iterations si quieres
    // double total_iterations = ...

    double current_iteration = 0.0;
    auto start_time = high_resolution_clock::now();

    for (double lambda1 = lambda_min; lambda1 <= lambda_max; lambda1 += step) {
        for (double lambda2 = lambda_min; lambda2 <= lambda_max; lambda2 += step) {
            ...
            for (double beta = beta_min; beta <= beta_max; beta += beta_step) {

                current_iteration += 1.0;
                auto current_time = high_resolution_clock::now();
                double elapsed_time = duration<double>(current_time - start_time).count();
                double progress = (current_iteration / /* total_iterations */);
                print_progress(progress, elapsed_time, /* total_iterations */, current_iteration);

                // check positivity (usa la librería)
                if (!check_positivity(lambda1, lambda2, lambda3, lambda4, lambda5)) {
                    continue;
                }

                // 2HDMC code ...
                // ...
                
                // Al final, para escribir CSV:
                vector<double> row_values = {
                    lambda1, lambda2, lambda3, lambda4, lambda5,
                    m12_squared, beta, tan(beta),
                    width_hgaga, total_width_h2, branching_ratio
                };
                write_csv_row(results, row_values); // de la librería
            }
        }
    }

    results.close();
    cout << endl << "Parameter scan completed. Results saved to " << output_file << endl;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "Usage: ./ParamScan output_filename [step_size]\n";
        return -1;
    }
    string output_file = argv[1];
    double step = 1.0;
    if (argc > 2) {
        step = stod(argv[2]);
    }
    perform_param_scan(output_file, step);
    return 0;
}
