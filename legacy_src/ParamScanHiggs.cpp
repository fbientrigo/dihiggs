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

using namespace std;
using namespace std::chrono;

/*

Scanning of Higgs Parameters
Assumes $\beta=0$ and changes $\lambda_j$ assuming $\lambda_{6,7}=0$


*/

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
// Expandable changing the columns
void write_csv_header(ofstream &results, const vector<string> &columns) {
    for (size_t i = 0; i < columns.size(); ++i) {
        results << columns[i];
        if (i < columns.size() - 1) results << ",";
    }
    results << endl;
}

// Function to write a row to CSV
void write_csv_row(ofstream &results, const vector<double> &values) {
    for (size_t i = 0; i < values.size(); ++i) {
        results << values[i];
        if (i < values.size() - 1) results << ",";
    }
    results << endl;
}

// Function to get Higgs masses
void get_higgs_masses(THDM &model, double &mass_h, double &mass_H, double &mass_A, double &mass_Hp) {
    double mh[5], a, l6, l7, tb, m12_2;
    model.get_param_phys(mh[1], mh[2], mh[3], mh[4], a, l6, l7, m12_2, tb);
    mass_h = mh[1];
    mass_H = mh[2];
    mass_A = mh[3];
    mass_Hp = mh[4];
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
        "lambda1", "lambda2", "lambda3", "lambda4", "lambda5",
        "width_hgaga", "total_width_h2", "branching_ratio_hgaga",
        "mass_h", "mass_H", "mass_A", "mass_Hp"
    };

    write_csv_header(results, columns);

    // Ranges for lambda parameters
    double sixteenPI = 50.24;
    double infinity = 300;
    double lambda1_min = 0.01, lambda1_max = sixteenPI;
    double lambda2_min = 0.01, lambda2_max = sixteenPI;
    double lambda3_min = -sixteenPI, lambda3_max = infinity;
    double lambda4_min = -sixteenPI, lambda4_max = infinity;
    double lambda5_min = -sixteenPI, lambda5_max = infinity;

    // Calculate total iterations for progress tracking
    double total_iterations = ((lambda1_max - lambda1_min) / step + 1) *
                              ((lambda2_max - lambda2_min) / step + 1) *
                              ((lambda3_max - lambda3_min) / step + 1) *
                              ((lambda4_max - lambda4_min) / step + 1) *
                              ((lambda5_max - lambda5_min) / step + 1);
    double current_iteration = 0;

    auto start_time = high_resolution_clock::now();

    for (double lambda1 = lambda1_min; lambda1 <= lambda1_max; lambda1 += step) {
        for (double lambda2 = lambda2_min; lambda2 <= lambda2_max; lambda2 += step) {
            for (double lambda3 = lambda3_min; lambda3 <= lambda3_max; lambda3 += step) {
                for (double lambda4 = lambda4_min; lambda4 <= lambda4_max; lambda4 += step) {
                    for (double lambda5 = lambda5_min; lambda5 <= lambda5_max; lambda5 += step) {

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

                        bool pset = model.set_param_higgs(lambda1, lambda2, lambda3, lambda4, lambda5, 0.0, 0.0, 200.0); // Example mHp = 200 GeV

                        if (!pset) {
                            continue;
                        }

                        Constraints check(model);
                        if (!check.check_positivity() || !check.check_unitarity() || !check.check_perturbativity()) {
                            continue;
                        }

                        // Calculate decay widths and masses
                        DecayTable table(model);
                        double width_hgaga = table.get_gamma_hgaga(2); // H2 -> gamma gamma
                        double total_width_h2 = table.get_gammatot_h(2); // Total width of H2

                        // Calculate branching ratio
                        double branching_ratio_hgaga = 0.0;
                        if (total_width_h2 > 0) {
                            branching_ratio_hgaga = width_hgaga / total_width_h2;
                        }

                        // Get Higgs masses
                        double mass_h, mass_H, mass_A, mass_Hp;
                        get_higgs_masses(model, mass_h, mass_H, mass_A, mass_Hp);

                        // Write results to CSV
                        vector<double> values = {
                            lambda1, lambda2, lambda3, lambda4, lambda5,
                            width_hgaga, total_width_h2, branching_ratio_hgaga,
                            mass_h, mass_H, mass_A, mass_Hp
                        };

                        write_csv_row(results, values);
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
        cout << "Usage: ./ParamScan output_filename\n";
        return -1;
    }

    string output_file = argv[1];
    // Conversión segura usando std::stod
    double step = 1.0;
    try {
        step = std::stod(argv[2]); // Step size for parameter variation
        std::cout << "Step size: " << step << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: Invalid step size provided. Please provide a valid number." << std::endl;
        return 1;
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: Step size out of range." << std::endl;
        return 1;
    }

    perform_param_scan(output_file, step);
    return 0;
}
