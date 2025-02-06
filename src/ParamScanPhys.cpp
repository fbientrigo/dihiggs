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

// Estructura para almacenar temporalmente el conjunto de parámetros
struct ParamSet {
    double m_phi;
    double mA;       // m_A = m_Hp
    double alpha;
    double beta;
    double lambda6;
    double lambda7;
    double m12;
};

// Variable global para guardar la mejor BR(h->gamma gamma)
static ParamSet g_bestParams;
static double   g_bestBR = -1.0;

// ----------------------------------------------------------------------------------------------
// Función principal de escaneo
// ----------------------------------------------------------------------------------------------
void perform_param_scan_withGD(const ConfigPhys &cfg, const string &output_file) {

    // Abrir archivo de resultados
    ofstream results(output_file);
    if (!results.is_open()) {
        cerr << "Failed to open output file: " << output_file << endl;
        return;
    }

    // Definir cabeceras de columnas para el CSV
    vector<string> columns = {
        "m_phi", "mA",      // (mA = m_Hp)
        "alpha", "beta",
        "lambda6", "lambda7",
        "m12",
        "sin_ba", "tan_beta",
        "positivity_ok", "unitarity_ok", "perturbativity_ok",
        "width_h2_bb", "width_h2_tautau", "width_h2_WW", "width_h2_ZZ", "width_h2_gaga",
        "total_width_h2",
        "branching_ratio_h2_gaga"
    };
    write_csv_header(results, columns);

    // ------------------------------------------------------------------------------------------
    // Calcular el número de pasos (ejemplo, ajustando a cómo se definan en el Config)
    // ------------------------------------------------------------------------------------------
    int steps_lambda6 = (cfg.lambda6_max - cfg.lambda6_min) / cfg.step_lambda6 + 1;
    int steps_lambda7 = (cfg.lambda7_max - cfg.lambda7_min) / cfg.step_lambda7 + 1;
    int steps_m12     = (cfg.m12_squared_max - cfg.m12_squared_min) / cfg.step_m12_squared + 1;
    int steps_alpha   = (cfg.alpha_max - cfg.alpha_min) / cfg.step_alpha + 1;
    int steps_beta    = (cfg.beta_max - cfg.beta_min) / cfg.step_beta + 1;

    // Pasos para las masas
    int steps_mphi = (cfg.mphi_max - cfg.mphi_min) / cfg.step_mphi + 1;
    int steps_mA   = (cfg.mA_max   - cfg.mA_min)   / cfg.step_mA   + 1;

    if (cfg.step_lambda6 == 0.0) {
        steps_lambda6 = 1;
    }
    if (cfg.step_lambda7 == 0.0) {
        steps_lambda7 = 1;
    }
    if (cfg.step_m12_squared == 0.0) {
        steps_m12 = 1;
    }
    if (cfg.step_alpha == 0.0) {
        steps_alpha = 1;
    }
    if (cfg.step_beta == 0.0) {
        steps_beta = 1;
    }
    if (cfg.step_mphi == 0.0) {
        steps_mphi = 1;
    }
    if (cfg.step_mA == 0.0) {
        steps_mA = 1;
    }

    // Para estimar iteraciones totales (productorio)
    double total_iterations = static_cast<double>(steps_lambda6)
                            * steps_lambda7
                            * steps_m12
                            * steps_alpha
                            * steps_beta
                            * steps_mphi
                            * steps_mA;
    cout << "Total iterations: " << total_iterations << "\n";
    double current_iteration = 0.0;
    auto start_time = high_resolution_clock::now();

    // Inicializar la mejor BR global
    g_bestBR = -1.0;

    // ------------------------------------------------------------------------------------------
    // Bucles anidados externos para parámetros "poco variables" o fijos
    // (lambda6, lambda7, m12, alpha, beta).
    // ------------------------------------------------------------------------------------------
    for (int i_l6 = 0; i_l6 < steps_lambda6; i_l6++) {
        double lambda6 = cfg.lambda6_min + i_l6 * cfg.step_lambda6;

        for (int i_l7 = 0; i_l7 < steps_lambda7; i_l7++) {
            double lambda7 = cfg.lambda7_min + i_l7 * cfg.step_lambda7;

            for (int i_m12 = 0; i_m12 < steps_m12; i_m12++) {
                double m12 = cfg.m12_squared_min + i_m12 * cfg.step_m12_squared;

                for (int i_alpha = 0; i_alpha < steps_alpha; i_alpha++) {
                    double alpha = cfg.alpha_min + i_alpha * cfg.step_alpha;

                    for (int i_beta = 0; i_beta < steps_beta; i_beta++) {
                        double beta = cfg.beta_min + i_beta * cfg.step_beta;

                        double tan_beta = std::tan(beta);
                        double sin_ba   = std::sin(beta - alpha);

                        // ----------------------------------------------------------------------------------
                        // Paralelización: bucles para m_phi y m_A ( = m_{H^+} )
                        // Se desea que  m_phi <= m_A  y  ambos hasta ~400.
                        // ----------------------------------------------------------------------------------
                        #pragma omp parallel 
                        {
                            // Cada hilo guardará temporalmente sus filas para luego volcarlas a "results".
                            vector<vector<double>> thread_local_rows;
                            double local_bestBR = -1.0;
                            ParamSet local_bestParams;

                            // collapse(2) para que OpenMP itere en 2D
                            #pragma omp for collapse(2) schedule(dynamic)
                            
                            for (int i_mphi = 0; i_mphi < steps_mphi; i_mphi++) {
                                for (int i_mA = 0; i_mA < steps_mA; i_mA++) {

                                    // Determinar valores físicos de las masas
                                    double m_phi = cfg.mphi_min + i_mphi * cfg.step_mphi;
                                    double mA    = cfg.mA_min   + i_mA   * cfg.step_mA;

                                    // Imponer m_phi <= mA (mA = m_{H^+})
                                    if (mA > m_phi) {
                                        // Si se requiere estrictamente que mA >= m_phi, se ignora este punto
                                        // continue;
                                    }

                                    // Fijar m_h = 125.0
                                    double m_h  = 125.0; 
                                    double m_Hp = mA;     // por la condición m_A = m_{H^+}

                                    // ======= Construir el modelo 2HDM e intentar fijar parámetros físicos
                                    THDM model;
                                    SM sm;
                                    model.set_SM(sm);

                                    bool pset = model.set_param_phys(
                                        m_h,
                                        m_phi,
                                        mA,
                                        m_Hp,
                                        sin_ba,
                                        lambda6,
                                        lambda7,
                                        m12,
                                        tan_beta
                                    );
                                    // if (!pset) {
                                    //     // Si no se pudo setear este punto, descartamos
                                    //     #pragma omp atomic
                                    //     current_iteration += 1.0;
                                    //     //continue;
                                    // }

                                    // Chequeos de estabilidad, unitariedad, etc.
                                    Constraints check(model);
                                    bool positivity_ok     = check.check_positivity();
                                    bool unitarity_ok      = check.check_unitarity();
                                    bool perturbativity_ok = check.check_perturbativity();

                                    // Calcular anchos de decaimiento
                                    DecayTable table(model);
                                    // "2" corresponde al Higgs pesado en la librería
                                    double w_h2_bb     = table.get_gamma_hdd(2, 3, 3);
                                    double w_h2_tautau = table.get_gamma_hll(2, 3, 3);
                                    double w_h2_WW     = table.get_gamma_hvv(2, 3);
                                    double w_h2_ZZ     = table.get_gamma_hvv(2, 2);
                                    double w_h2_gaga   = table.get_gamma_hgaga(2);
                                    double w_total_h2  = table.get_gammatot_h(2);

                                    double br_h2_gaga  = (w_total_h2 > 1e-12) ? w_h2_gaga / w_total_h2 : 0.0;

                                    // Guardar la fila de resultados en el buffer local
                                    vector<double> row = {
                                        m_phi, mA,
                                        alpha, beta,
                                        lambda6, lambda7,
                                        m12,
                                        sin_ba, tan_beta,
                                        positivity_ok     ? 1.0 : 0.0,
                                        unitarity_ok      ? 1.0 : 0.0,
                                        perturbativity_ok ? 1.0 : 0.0,
                                        w_h2_bb,
                                        w_h2_tautau,
                                        w_h2_WW,
                                        w_h2_ZZ,
                                        w_h2_gaga,
                                        w_total_h2,
                                        br_h2_gaga
                                    };
                                    thread_local_rows.push_back(row);

                                    // Actualizar mejor BR local
                                    if (positivity_ok && unitarity_ok && perturbativity_ok && w_total_h2 > 0.0) {
                                        if (br_h2_gaga > local_bestBR) {
                                            local_bestBR = br_h2_gaga;
                                            local_bestParams = {m_phi, mA, alpha, beta, lambda6, lambda7, m12};
                                        }
                                    }


                                } // end for i_mA
                            } // end for i_mphi

                            // Escribir en el archivo principal de resultados
                            #pragma omp critical
                            {
                                for (auto &row : thread_local_rows) {
                                    write_csv_row(results, row);
                                }
                                thread_local_rows.clear();

                                // Actualizar mejor BR global
                                if (local_bestBR > g_bestBR) {
                                    g_bestBR = local_bestBR;
                                    g_bestParams = local_bestParams;
                                }
                            }
                                                        // Mostrar avance y ETA
                            #pragma omp critical
                            {
                                current_iteration += 1.0;
                                double elapsed = duration<double>(
                                    high_resolution_clock::now() - start_time
                                ).count();
                                double progress = (current_iteration / total_iterations) * 100.0;
                                double est_time_left = (elapsed / current_iteration)
                                                        * (total_iterations - current_iteration);

                                cout << fixed << setprecision(2)
                                        << "Progress: " << progress << "% | "
                                        << "Elapsed: " << (elapsed / 60.0) << " min | "
                                        << "ETA: " << (est_time_left / 60.0) << " min"
                                        << "\r" << flush;
                            }
                        } // end #pragma omp parallel
                    } // end for i_beta
                } // end for i_alpha
            } // end for i_m12
        } // end for i_l7
    } // end for i_l6

    // Cerrar archivo de resultados
    results.close();
    cout << "\n\nScan completed. Results saved to " << output_file << endl;
    cout << "Best BR(h->gaga) found = " << g_bestBR << endl;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <config_file> <output_csv>\n";
        return 1;
    }

    string config_file = argv[1];
    string output_file = argv[2];

    try {
        ConfigPhys cfg = readPhysConfig(config_file);
        perform_param_scan_withGD(cfg, output_file);
    } catch(const exception &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
