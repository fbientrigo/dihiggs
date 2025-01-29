#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>

/*
 * ParamScanGen - Parameter Scan for the General Basis of the Two-Higgs Doublet Model (2HDM)
 * 
 * This program performs a parameter scan for the general basis of the 2HDM, varying the 
 * parameters \( \lambda_1, \lambda_2, \lambda_3, \lambda_4, \lambda_5 \), \( m_{12}^2 \), 
 * and \( \beta \) according to user-defined ranges. The constraints for positivity, unitarity, 
 * and perturbativity are applied to filter valid parameter sets. The program calculates 
 * decay widths, total widths, and branching ratios for \( H_2 \to \gamma \gamma \) using 2HDMC.
 * 
 * Parameters:
 *  - Output file: CSV file to store the results of the scan.
 *  - Step size: Step size for varying the parameters (optional, default = 1.0).
 *
 * Output:
 *  The program outputs the parameter combinations and their corresponding physical observables, 
 *  including decay widths and branching ratios, into the specified CSV file.
 */

using namespace std;
using namespace std::chrono;

// Function to check positivity of the potential constraints
bool check_positivity(double lambda1, double lambda2, double lambda3, double lambda4, double lambda5) {
    if (lambda1 <= 0 || lambda2 <= 0) return false;
    if (lambda3 <= -sqrt(lambda1 * lambda2)) return false;
    if (lambda3 + lambda4 - abs(lambda5) <= -sqrt(lambda1 * lambda2)) return false;
    return true;
}

// Function to print progress bar with estimated time
void print_progress(double progress, double elapsed_time, double total_iterations, double current_iteration) {
    int bar_width = 50; // Width of the progress bar
    cout << "[";
    int pos = bar_width * progress;
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) cout << "=";
        else if (i == pos) cout << ">";
        else cout << " ";
    }
    double time_per_iteration = elapsed_time / current_iteration;
    double remaining_time = time_per_iteration * (total_iterations - current_iteration);
    cout << "] " << int(progress * 100.0) << "% | Elapsed: "
         << int(elapsed_time) << "s | Remaining: " << int(remaining_time) << "s\r";
    cout.flush();
}

// Function to write header to CSV
void write_csv_header(ofstream &results, const vector<string> &columns) {
    for (size_t i = 0; i < columns.size(); ++i) {
        results << columns[i];
        if (i < columns.size() - 1) results << ",";
    }
    results << endl;
    results.flush();
}

void write_csv_row(ofstream &results, const vector<double> &values) {
    if (!results.is_open()) {
        cerr << "Error: Output file is not open for writing." << endl;
        return;
    }
    for (size_t i = 0; i < values.size(); ++i) {
        results << values[i];
        if (i < values.size() - 1) results << ",";
    }
    results << endl;
    results.flush(); // Ensure data is written immediately
}

// Function to perform parameter scan
void perform_param_scan(const string &output_file, double step) {
    ofstream results(output_file);

    if (!results.is_open()) {
        cerr << "Failed to open output file." << endl;
        return;
    }

    // -------------------- Define CSV columns ------------------------
    vector<string> columns = {
        "lambda1", "lambda2", "lambda3", "lambda4", "lambda5", "m12_squared", "beta", "tan_beta",
        "width_hgaga", "total_width_h2", "branching_ratio_hgaga"
    };

    write_csv_header(results, columns);

    // Ranges for lambda parameters and other variables
    double infty_neg = -50.0, infty_pos = 100.0;
    double beta_min = 0.01, beta_max = M_PI / 2; // 0.1
    double beta_step = 0.5;
    double lambda_min = 0.01, lambda_max = 16.0 * M_PI;
    double lambda_345_min = -lambda_max;
    double lambda_345_max = infty_pos;

    // Calculate total iterations for progress tracking
    double total_iterations = ((lambda_max - lambda_min) / step + 1) * // lambda1
                              ((lambda_max - lambda_min) / step + 1) * // lambda2
                              ((lambda_345_max - lambda_345_min) / step + 1) * // lambda3
                              ((lambda_345_max - lambda_345_min) / step + 1) * // lambda4
                              ((lambda_345_max - lambda_345_min) / step + 1) * // lambda5
                              ((infty_pos - infty_neg) / step + 1) * //m12 es cercan a 0 TODO
                              ((beta_max - beta_min) / beta_step + 1); // beta es cercano a 0 TODO
    double current_iteration = 0;

    auto start_time = high_resolution_clock::now();

    for (double lambda1 = lambda_min; lambda1 <= lambda_max; lambda1 += step) {
        for (double lambda2 = lambda_min; lambda2 <= lambda_max; lambda2 += step) {
            for (double lambda3 = lambda_345_min; lambda3 <= lambda_345_max; lambda3 += step) {
                for (double lambda4 = lambda_345_min; lambda4 <= lambda_345_max; lambda4 += step) {
                    for (double lambda5 = lambda_345_min; lambda5 <= lambda_345_max; lambda5 += step) {
                        for (double m12_squared = infty_neg; m12_squared <= infty_pos; m12_squared += step) {
                            for (double beta = beta_min; beta <= beta_max; beta += beta_step) {
                                // // Debugging: Print the current parameters
                                // cout << "Current parameters: "
                                //      << "lambda1=" << lambda1 << ", lambda2=" << lambda2 << ", lambda3=" << lambda3
                                //      << ", lambda4=" << lambda4 << ", lambda5=" << lambda5 << ", m12_squared=" << m12_squared
                                //      << ", beta=" << beta << endl;

                                // Calculate tan(beta)
                                double tan_beta = tan(beta);

                                // Update progress
                                current_iteration++;
                                auto current_time = high_resolution_clock::now();
                                double elapsed_time = duration<double>(current_time - start_time).count();
                                print_progress(current_iteration / total_iterations, elapsed_time, total_iterations, current_iteration);

                                // Check positivity constraints
                                if (!check_positivity(lambda1, lambda2, lambda3, lambda4, lambda5)) {
                                    continue;
                                }

                                // Initialize model
                                THDM model;
                                SM sm;
                                model.set_SM(sm);

                                bool pset = model.set_param_gen(lambda1, lambda2, lambda3, lambda4, lambda5, 0.0, 0.0, m12_squared, tan_beta);

                                if (!pset) {
                                    continue;
                                }

                                Constraints check(model);
                                if (!check.check_positivity() || !check.check_unitarity() || !check.check_perturbativity()) {
                                    continue;
                                }

                                // Calculate decay widths
                                DecayTable table(model);
                                double width_hgaga = table.get_gamma_hgaga(2); // H2 -> gamma gamma
                                double total_width_h2 = table.get_gammatot_h(2); // Total width of H2

                                // Validate decay widths
                                if (std::isnan(width_hgaga) || std::isnan(total_width_h2) || total_width_h2 <= 0) {
                                    continue;
                                }

                                // Calculate branching ratio
                                double branching_ratio_hgaga = width_hgaga / total_width_h2;

                                // Write results to CSV
                                vector<double> values = {
                                    lambda1, lambda2, lambda3, lambda4, lambda5, m12_squared, beta, tan_beta,
                                    width_hgaga, total_width_h2, branching_ratio_hgaga
                                };
                          try {
                                    write_csv_row(results, values);
                                } catch (const exception &e) {
                                    cerr << "Error writing to CSV: " << e.what() << endl;
                                }
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
    if (argc < 2) {
        cout << "Usage: ./ParamScan output_filename [step_size]\n";
        return -1;
    }

    string output_file = argv[1];
    double step = 1.0;
    if (argc > 2) {
        try {
            step = stod(argv[2]);
        } catch (const exception &e) {
            cerr << "Invalid step size provided. Using default step size of 1.0." << endl;
        }
    }

    perform_param_scan(output_file, step);
    return 0;
}
