/**
 * @file PhysScanWithFixings.cpp
 * @brief Escaneo de parámetros m_phi y lambda1 con parámetros físicos fijos.
 *
 * Refactored to use Analytical Inversion Strategy for m12^2 based on lambda1.
 * Included Real-Time Statistical Monitoring (VPPM & PPS) with optimized granularity.
 */

#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <vector>
#include <chrono>
#include <atomic> // [MONITORING] Necesario para contadores thread-safe sin bloqueos

using namespace std;
using namespace std::chrono;

constexpr double PI = std::acos(-1.0);
constexpr double VEV = 246.0;

using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;

// Estructura para resultados
struct ParamSet {
    double m_phi, mA, alpha, beta, lambda6, lambda7, m12;
};

static double g_bestBR = -1.0;

void perform_param_scan_fixings(
    double mphi_min, double mphi_max, int N_mphi,
    double lam1_min, double lam1_max, int N_lam1,
    double mA_fixed, double sin_ba,   double tan_beta,
    double lambda6,  double lambda7,
    const string &output_file
) {
    ofstream results(output_file);
    if (!results.is_open()) {
        cerr << "No pude abrir: " << output_file << endl;
        return;
    }

    // Cabecera CSV
    vector<string> columns = {
        "m_phi","mA","alpha","beta","lambda6","lambda7","m12",
        "sin_ba","tan_beta","positivity_ok","unitarity_ok","perturbativity_ok",
        "width_bb","width_tautau","width_WW","width_ZZ",
        "width_gaga","width_Zga","width_gg","width_hh",
        "total_width","br_gaga", "lam1", "computed_lam1", "lam2","computed_lam2", "lam3", "lam4", "lam5"
    };
    write_csv_header(results, columns);

    // Calcular pasos
    double step_mphi = (N_mphi > 1) ? (mphi_max - mphi_min)/(N_mphi - 1) : 0.0;
    double step_lam1 = (N_lam1 > 1) ? (lam1_max - lam1_min)/(N_lam1 - 1) : 0.0;

    // [MONITORING] Variables estadísticas globales y temporizadores
    long long total_tasks = (long long)N_mphi * (long long)N_lam1; 
    std::atomic<long long> global_attempts(0); // Total puntos calculados (válidos + inválidos)
    std::atomic<long long> global_valid(0);    // Total puntos guardados
    
    auto start_time = Clock::now();
    auto last_report_time = start_time;

    // [HPC ADJUSTMENT] Configuracion de granularidad de reporte
    const double REPORT_INTERVAL_SEC = 0.5; // Reportar cada 0.5s para feedback rapido
    const int BATCH_SIZE = 50;              // Forzar escritura/reporte cada 50 puntos validos por hilo

    cout << "step_lam1=" << step_lam1 << endl;
    cout << "step_mphi=" << step_mphi << endl;
    cout << "Total iterations to compute: " << total_tasks << endl;

    g_bestBR = -1.0;

    const double mh = 125.0;

    // Step A: Precise Angle Initialization
    const double inv = 1.0 / std::sqrt(1.0 + tan_beta*tan_beta); 
    const double cb  = inv;                  
    const double sb  = tan_beta * inv;       
    const double cba = std::sqrt(1.0 - sin_ba*sin_ba); 

    // Physics basis angles needed for lambda inversion
    const double ca = cb * cba + sb * sin_ba; 
    const double sa = sb * cba - cb * sin_ba; 
    
    // Angles for reporting/logging
    const double beta  = std::atan(tan_beta);
    const double alpha = beta - std::asin(sin_ba);

    #pragma omp parallel
    {
        // [MONITORING] Contador local para evitar contención en el bus de memoria (False Sharing)
        long long local_attempts = 0; 
        
        vector<vector<double>> buffer;
        buffer.reserve(BATCH_SIZE); // Pre-reservar memoria

        #pragma omp for schedule(dynamic)
        for(int i_phi = 0; i_phi < N_mphi; ++i_phi) {
            
            double m_phi = mphi_min + i_phi * step_mphi;

            for(int i_l1 = 0; i_l1 < N_lam1; ++i_l1) {
                
                local_attempts++;
                
                // Actualizar contadores globales en lotes mas grandes (reduce contencion atomica)
                if (local_attempts >= 200) {
                    global_attempts += local_attempts;
                    local_attempts = 0;
                }

                double lambda1 = lam1_min + i_l1 * step_lam1;

                // Step C: Analytical Inversion (m12^2 Calculation)
                // m12^2 = [ (m_phi^2 ca^2 + mh^2 sa^2)/(v^2 cb^2) - 1.5*l6*tb + 0.5*l7*tb^3 - l1 ] * (v^2 cb^2 / tb)
                // Se despeja m12 de la ecuacion de lambda1
                double term_mass = (m_phi*m_phi * ca*ca + mh*mh * sa*sa);
                double term_l6   = 1.5 * lambda6 * tan_beta;
                double term_l7   = 0.5 * lambda7 * std::pow(tan_beta, 3);
                
                // Pre-factor sin multiplicar por el ultimo termino geometrico
                // Nota: lambda1 esta restando en la formula despejada
                double pre_factor = (term_mass / (VEV*VEV * cb*cb)) - term_l6 + term_l7 - lambda1;

                double m12_sq = pre_factor * (VEV*VEV * cb*cb / tan_beta);

                // Step 3: Validation & Model Setting
                // Calculate lambda2 con el m12_sq calculado
                double l2 = calc_lambda2(mh, m_phi, m12_sq, sa, ca, sb, cb, tan_beta, lambda6, lambda7, VEV);

                // Filter: Check Perturbativity of lambda2
                if (std::abs(l2) > 4 * PI) continue;

                // Set Model
                THDM model;
                SM sm; model.set_SM(sm);
                
                // Usamos m12_sq calculado
                bool pset = model.set_param_phys(mh, m_phi, mA_fixed, mA_fixed, sin_ba, lambda6, lambda7, m12_sq, tan_beta);
                
                if (!pset) continue;

                model.set_yukawas_type(1);

                double lam1_g, lam2_g, lam3_g, lam4_g, lam5_g, lam6_g, lam7_g, m12_2_g, tanb_g;
                model.get_param_gen(lam1_g, lam2_g, lam3_g, lam4_g, lam5_g, lam6_g, lam7_g, m12_2_g, tanb_g);

                Constraints check(model);
                bool pos = check.check_positivity();
                bool uni = check.check_unitarity();
                bool pert = check.check_perturbativity();

                DecayTable tab(model);
                double w_bb     = tab.get_gamma_hdd(2,3,3);
                double w_tautau = tab.get_gamma_hll(2,3,3);
                double w_WW     = tab.get_gamma_hvv(2,3);
                double w_ZZ     = tab.get_gamma_hvv(2,2);
                double w_gaga   = tab.get_gamma_hgaga(2);
                double w_Zga    = tab.get_gamma_hZga(2);
                double w_gg     = tab.get_gamma_hgg(2);
                double w_hh     = tab.get_gamma_hhh(2,1,1);
                double w_tot    = tab.get_gammatot_h(2);
                double br_gaga  = (w_tot > 1e-15) ? w_gaga / w_tot : 0.0;

                // Guardar fila en buffer local
                buffer.push_back(std::vector<double>{
                    m_phi, mA_fixed, alpha, beta, lambda6, lambda7, m12_sq,
                    sin_ba, tan_beta,
                    pos?1.0:0.0, uni?1.0:0.0, pert?1.0:0.0,
                    w_bb, w_tautau, w_WW, w_ZZ,
                    w_gaga, w_Zga, w_gg, w_hh,
                    w_tot, br_gaga, lambda1, lam1_g, l2, lam2_g, lam3_g, lam4_g, lam5_g
                });

                // [HPC FLUSH] Si el buffer local se llena, forzamos escritura y reporte
                // Esto soluciona el problema de granularidad en scans cortos
                if (buffer.size() >= BATCH_SIZE) {
                    #pragma omp critical
                    {
                        global_valid += buffer.size();
                        for (const auto &row : buffer) {
                            write_csv_row(results, row);
                        }
                        buffer.clear(); 

                        // Check timer
                        auto now = Clock::now();
                        std::chrono::duration<double> diff = now - last_report_time;
                        
                        if (diff.count() > REPORT_INTERVAL_SEC) {
                            // Sincronizar contador de intentos para el reporte
                            long long current_attempts = global_attempts.load() + local_attempts; // Sumar local pendiente
                            long long current_valid = global_valid.load();
                            
                            auto total_elapsed = std::chrono::duration<double>(now - start_time).count();
                            
                            double pps = (total_elapsed > 0) ? (double)current_attempts / total_elapsed : 0.0;
                            double vppm = (total_elapsed > 0) ? (double)current_valid / (total_elapsed / 60.0) : 0.0;
                            double progress = (total_tasks > 0) ? (100.0 * current_attempts / total_tasks) : 0.0;

                            std::cerr << "\r" 
                                      << "[ " << std::fixed << std::setprecision(1) << progress << "% ] "
                                      << "Speed: " << (long)pps << " pts/s | "
                                      << "Valid: " << (long)vppm << "/min | "
                                      << "Found: " << current_valid << "/" << current_attempts
                                      << "   " << std::flush;
                            
                            last_report_time = now;
                        }
                    }
                }

            } // Fin loop interno
        } // Fin loop externo

        // [LIMPIEZA FINAL] Volcar remanentes al terminar el hilo
        if (local_attempts > 0) {
            global_attempts += local_attempts;
        }

        #pragma omp critical
        {
            if (!buffer.empty()) {
                global_valid += buffer.size(); 
                for (const auto &row : buffer) {
                    write_csv_row(results, row);
                }
                buffer.clear();
            }
        }
    } // End parallel

    // [MONITORING] Reporte final
    auto total_time = std::chrono::duration<double>(Clock::now() - start_time).count();
    long long final_attempts = global_attempts.load();
    long long final_valid = global_valid.load();

    std::cout << "\n\n--- Scan Finalizado ---" << endl;
    std::cout << "Total Attempts: " << final_attempts << endl;
    std::cout << "Total Valid Points: " << final_valid << endl;
    std::cout << "Total Time: " << total_time << " s" << endl;
    if (total_time > 0) {
        std::cout << "Average Speed: " << (long)(final_attempts / total_time) << " pts/s" << endl;
        std::cout << "Final Valid Rate: " << (long)(final_valid / (total_time/60.0)) << " valid/min" << endl;
    }

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
    int    N_mphi   = atoi(argv[3]);
    double lam1_min = atof(argv[4]);
    double lam1_max = atof(argv[5]);
    int    N_lam1   = atoi(argv[6]);
    double mA_fixed = atof(argv[7]);
    double sin_ba   = atof(argv[8]);
    double tan_beta = atof(argv[9]);
    double lambda6  = atof(argv[10]);
    double lambda7  = atof(argv[11]);
    string output   = argv[12];

    perform_param_scan_fixings(
        mphi_min, mphi_max, N_mphi,
        lam1_min, lam1_max, N_lam1,
        mA_fixed, sin_ba,   tan_beta,
        lambda6,  lambda7,
        output
    );
    return 0;
}