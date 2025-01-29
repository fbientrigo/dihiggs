#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"  // Header con tus funciones de escritura CSV
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;

// -----------------------------------------------------------------------------
// Función principal de escaneo
// -----------------------------------------------------------------------------
static void perform_param_scan(const Config &cfg, const string &output_file) {

    ofstream results(output_file);
    if (!results.is_open()) {
        cerr << "Failed to open output file: " << output_file << endl;
        return;
    }

    // -------------------------------------------------------------------------
    // Definir las columnas CSV (puedes añadir/quitar según tus necesidades)
    // -------------------------------------------------------------------------
    vector<string> columns = {
        "lambda1", "lambda2", "lambda3", "lambda4", "lambda5",
        "m12_squared", "beta", "tan_beta",

        // Checks de validez
        "positivity_ok", "unitarity_ok", "perturbativity_ok",

        // Anchos parciales del estado h=2 (H)
        "width_h2_bb",       // H -> b b
        "width_h2_tautau",   // H -> tau tau
        "width_h2_WW",       // H -> W W
        "width_h2_ZZ",       // H -> Z Z
        "width_hgaga",       // H -> gamma gamma

        // Ancho total y BR específico
        "total_width_h2",
        "branching_ratio_hgaga" // BR(H->gamma gamma)
    };

    // Cabecera CSV
    write_csv_header(results, columns);

    // -------------------------------------------------------------------------
    // Calcular número total de iteraciones para mostrar progreso
    // -------------------------------------------------------------------------
    double total_iterations = computeTotalIterations(cfg);
    double current_iteration = 0.0;
    auto start_time = high_resolution_clock::now();

    // -------------------------------------------------------------------------
    // Bucle principal de escaneo de parámetros
    // -------------------------------------------------------------------------
    for (double lambda1 = cfg.lambda1_min; lambda1 <= cfg.lambda1_max; lambda1 += cfg.step_lambda1) {
        for (double lambda2 = cfg.lambda2_min; lambda2 <= cfg.lambda2_max; lambda2 += cfg.step_lambda2) {
            for (double lambda3 = cfg.lambda3_min; lambda3 <= cfg.lambda3_max; lambda3 += cfg.step_lambda3) {
                for (double lambda4 = cfg.lambda4_min; lambda4 <= cfg.lambda4_max; lambda4 += cfg.step_lambda4) {
                    for (double lambda5 = cfg.lambda5_min; lambda5 <= cfg.lambda5_max; lambda5 += cfg.step_lambda5) {
                        for (double m12 = cfg.m12_squared_min; m12 <= cfg.m12_squared_max; m12 += cfg.step_m12_squared) {
                            for (double beta = cfg.beta_min; beta <= cfg.beta_max; beta += cfg.step_beta) {

                                // Actualizar progreso
                                current_iteration += 1.0;
                                auto now_time = high_resolution_clock::now();
                                double elapsed = duration<double>(now_time - start_time).count();
                                double progress = (current_iteration / total_iterations);
                                print_progress(progress, elapsed, total_iterations, current_iteration);

                                // ---------------------------------------------------------------------
                                // Check inicial de "positividad manual" (si lo deseas)
                                // ---------------------------------------------------------------------
                                if (!check_positivity(lambda1, lambda2, lambda3, lambda4, lambda5)) {
                                    // Si quieres descartar directamente, haz 'continue;'
                                    // En caso de querer registrar igualmente, no hagas continue
                                    continue;
                                }

                                // ---------------------------------------------------------------------
                                // Construir el modelo 2HDM
                                // ---------------------------------------------------------------------
                                THDM model;
                                SM sm;
                                model.set_SM(sm);

                                bool pset = model.set_param_gen(
                                    lambda1, lambda2, lambda3, lambda4, lambda5,
                                    0.0, 0.0,   // v1, v2 en la base gen, según tu uso
                                    m12,        // m12^2
                                    std::tan(beta)  // tan_beta
                                );
                                if (!pset) {
                                    // Parámetros no válidos para 2HDM
                                    continue;
                                }

                                // ---------------------------------------------------------------------
                                // Evaluar los "checks" de validez usando la clase Constraints
                                // ---------------------------------------------------------------------
                                Constraints check(model);
                                bool positivity_ok      = check.check_positivity();
                                bool unitarity_ok       = check.check_unitarity();
                                bool perturbativity_ok  = check.check_perturbativity();

                                // Si quisieras descartar de plano cuando algo falla, usar:
                                // if (!positivity_ok || !unitarity_ok || !perturbativity_ok) {
                                //     continue;
                                // }

                                // ---------------------------------------------------------------------
                                // Calcular anchos con DecayTable
                                // ---------------------------------------------------------------------
                                DecayTable table(model);

                                // Por ejemplo: anchos parciales de h=2 (H).
                                // Si uno de los checks falló, puedes decidir poner anchuras = 0
                                // o NaN. Aquí, en el ejemplo, sólo calculamos si positivity_ok.
                                // (Ajusta según tu preferencia)
                                double w_h2_bb      = positivity_ok
                                                      ? table.get_gamma_hdd(2, 3, 3)   // H -> b b
                                                      : 0.0;
                                double w_h2_tautau  = positivity_ok
                                                      ? table.get_gamma_hll(2, 3, 3)   // H -> tau tau (l=3 -> tau)
                                                      : 0.0;
                                double w_h2_WW      = positivity_ok
                                                      ? table.get_gamma_hvv(2, 3)      // H -> W W  (V=3 -> W)
                                                      : 0.0;
                                double w_h2_ZZ      = positivity_ok
                                                      ? table.get_gamma_hvv(2, 2)      // H -> Z Z  (V=2 -> Z)
                                                      : 0.0;
                                double w_hgaga      = positivity_ok
                                                      ? table.get_gamma_hgaga(2)       // H -> gamma gamma
                                                      : 0.0;

                                // Ancho total
                                double total_w_h2 = positivity_ok
                                                    ? table.get_gammatot_h(2)
                                                    : 0.0;

                                // BR(H -> gamma gamma)
                                double br_hgaga = 0.0;
                                if (total_w_h2 > 1e-12) { 
                                    br_hgaga = w_hgaga / total_w_h2;
                                }

                                // ---------------------------------------------------------------------
                                // Preparar fila de datos para escribir en el CSV
                                // ---------------------------------------------------------------------
                                // Si tu "write_csv_row" admite directamente strings,
                                // puedes convertir manualmente cada elemento a string.
                                // (Aquí se asume una versión simplificada donde la librería
                                //  admite un vector<double> o vector<string>).
                                // ---------------------------------------------------------------------
                                vector<string> row;

                                // 1) Parámetros escaneados
                                row.push_back(to_string(lambda1));
                                row.push_back(to_string(lambda2));
                                row.push_back(to_string(lambda3));
                                row.push_back(to_string(lambda4));
                                row.push_back(to_string(lambda5));
                                row.push_back(to_string(m12));
                                row.push_back(to_string(beta));
                                row.push_back(to_string(std::tan(beta)));

                                // 2) Booleans en forma "true"/"false"
                                row.push_back(positivity_ok     ? "true" : "false");
                                row.push_back(unitarity_ok      ? "true" : "false");
                                row.push_back(perturbativity_ok ? "true" : "false");

                                // 3) Anchos parciales y total
                                row.push_back(to_string(w_h2_bb));
                                row.push_back(to_string(w_h2_tautau));
                                row.push_back(to_string(w_h2_WW));
                                row.push_back(to_string(w_h2_ZZ));
                                row.push_back(to_string(w_hgaga));

                                // 4) Ancho total y BR
                                row.push_back(to_string(total_w_h2));
                                row.push_back(to_string(br_hgaga));

                                // Finalmente, se escribe la fila
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


// -----------------------------------------------------------------------------
// main()
// -----------------------------------------------------------------------------
int main(int argc, char *argv[]) {

    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <config_file> <output_csv>\n";
        return 1;
    }

    string config_file = argv[1];
    string output_file = argv[2];

    try {
        // 1) Leer config
        Config cfg = readConfig(config_file);

        // 2) Lanzar escaneo
        perform_param_scan(cfg, output_file);

    } catch(const exception &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}
