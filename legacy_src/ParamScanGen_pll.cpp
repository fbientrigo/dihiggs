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
#include <omp.h>
#include <sstream> // Para manejar búferes de escritura
#include <filesystem> // Para manejar archivos temporales

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

// Function to write header to CSV
void write_csv_header(ofstream &results, const vector<string> &columns) {
    for (size_t i = 0; i < columns.size(); ++i) {
        results << columns[i];
        if (i < columns.size() - 1) results << ",";
    }
    results << endl;
}

// Function to perform parameter scan
void perform_param_scan(const string &output_file, double step) {
    namespace fs = std::filesystem;

    // -------------------- Define CSV columns ------------------------
    vector<string> columns = {
        "lambda1", "lambda2", "lambda3", "lambda4", "lambda5", "m12_squared", "beta", "tan_beta",
        "width_hgaga", "total_width_h2", "branching_ratio_hgaga"
    };

    // Ranges for lambda parameters and other variables
    double infty_neg = -50.0, infty_pos = 100.0;
    double beta_min = 0.01, beta_max = M_PI / 2;
    double beta_step = 0.5;
    double lambda_min = 0.01, lambda_max = 16.0 * M_PI;
    double lambda_345_min = -lambda_max;
    double lambda_345_max = infty_pos;

    auto start_time = high_resolution_clock::now();

    // Precalcular los límites para cada bucle
    int steps_lambda = static_cast<int>((lambda_max - lambda_min) / step) + 1;
    int steps_lambda_345 = static_cast<int>((lambda_345_max - lambda_345_min) / step) + 1;
    int steps_m12 = static_cast<int>((infty_pos - infty_neg) / step) + 1;
    int steps_beta = static_cast<int>((beta_max - beta_min) / beta_step) + 1;

    double completed_outer_loops = 0;
    double total_outer_loops = steps_lambda * steps_lambda * steps_lambda_345 * steps_lambda_345;

    string temp_dir = "temp_results";
    fs::create_directory(temp_dir);

    for (int i1 = 0; i1 < steps_lambda; i1++) {
        for (int i2 = 0; i2 < steps_lambda; i2++) {
            for (int i3 = 0; i3 < steps_lambda_345; i3++) {
                for (int i4 = 0; i4 < steps_lambda_345; i4++) {

                    double lambda1 = lambda_min + i1 * step;
                    double lambda2 = lambda_min + i2 * step;
                    double lambda3 = lambda_345_min + i3 * step;
                    double lambda4 = lambda_345_min + i4 * step;

                    if (!check_positivity(lambda1, lambda2, lambda3, lambda4, 0.0)) {
                        continue;
                    }

                    string temp_filename = temp_dir + "/temp_" + to_string(i1) + "_" + to_string(i2) + ".csv";

                    #pragma omp parallel
                    {
                        stringstream local_buffer;

                        #pragma omp for collapse(3) schedule(dynamic)
                        for (int i5 = 0; i5 < steps_lambda_345; i5++) {
                            for (int i_m12 = 0; i_m12 < steps_m12; i_m12++) {
                                for (int i_beta = 0; i_beta < steps_beta; i_beta++) {
                                    double lambda5 = lambda_345_min + i5 * step;
                                    double m12_squared = infty_neg + i_m12 * step;
                                    double beta = beta_min + i_beta * beta_step;
                                    double tan_beta = tan(beta);

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

                                    DecayTable table(model);
                                    double width_hgaga = table.get_gamma_hgaga(2);
                                    double total_width_h2 = table.get_gammatot_h(2);

                                    if (std::isnan(width_hgaga) || std::isnan(total_width_h2) || total_width_h2 <= 0) {
                                        continue;
                                    }

                                    double branching_ratio_hgaga = width_hgaga / total_width_h2;

                                    local_buffer << lambda1 << "," << lambda2 << "," << lambda3 << ","
                                                 << lambda4 << "," << lambda5 << "," << m12_squared << ","
                                                 << beta << "," << tan_beta << "," << width_hgaga << ","
                                                 << total_width_h2 << "," << branching_ratio_hgaga << endl;
                                }
                            }
                        }

                        #pragma omp critical
                        {
                            ofstream temp_file(temp_filename, ios::app);
                            if (temp_file.is_open()) {
                                temp_file << local_buffer.str();
                                temp_file.close();
                            } else {
                                cerr << "Failed to write to temporary file: " << temp_filename << endl;
                            }
                        }
                    }

                    #pragma omp critical
                    {
                        completed_outer_loops++;
                        double progress = (completed_outer_loops / total_outer_loops) * 100.0;
                        auto current_time = high_resolution_clock::now();
                        double elapsed_time = duration<double>(current_time - start_time).count();
                        double remaining_time = (elapsed_time / completed_outer_loops) * (total_outer_loops - completed_outer_loops);
                        cout << "Progress: " << fixed << setprecision(2) << progress << "% | Elapsed: " << elapsed_time << "s | Remaining: " << remaining_time << "s\r" << flush;
                    }
                }
            }
        }
    }

    // Consolidar archivos temporales
    ofstream results(output_file);
    write_csv_header(results, columns);

    for (const auto &entry : fs::directory_iterator(temp_dir)) {
        ifstream temp_file(entry.path());
        results << temp_file.rdbuf();
        temp_file.close();
    }

    results.close();
    fs::remove_all(temp_dir);

    auto end_time = high_resolution_clock::now();
    double elapsed_time = duration<double>(end_time - start_time).count();
    cout << endl << "Parameter scan completed in " << elapsed_time << " seconds. Results saved to " << output_file << endl;
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