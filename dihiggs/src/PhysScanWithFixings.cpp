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
#include <sstream>   // <--- NUEVO: para buffer de texto por hilo
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <vector>
#include <chrono>
#include <atomic> // [MONITORING] Necesario para contadores thread-safe sin bloqueos

using namespace std;
using namespace std::chrono;

constexpr double PI  = std::acos(-1.0);
constexpr double VEV = 246.0;

using Clock    = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;

// Estructura para resultados (no esencial, pero la dejamos por compatibilidad)
struct ParamSet {
    double m_phi, mA, alpha, beta, lambda6, lambda7, m12;
};

static double g_bestBR = -1.0;
//

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
#include <sstream>   // <--- NUEVO: para buffer de texto por hilo
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <vector>
#include <chrono>
#include <atomic> // [MONITORING] Necesario para contadores thread-safe sin bloqueos

using namespace std;
using namespace std::chrono;



using Clock    = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;


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
        "total_width","br_gaga", "lam1", "computed_lam1",
        "lam2","computed_lam2","lam3","lam4","lam5"
    };
    write_csv_header(results, columns);

    // Calcular pasos
    double step_mphi = (N_mphi > 1) ? (mphi_max - mphi_min)/(N_mphi - 1) : 0.0;
    double step_lam1 = (N_lam1 > 1) ? (lam1_max - lam1_min)/(N_lam1 - 1) : 0.0;

    // [MONITORING] Variables estadísticas globales y temporizadores
    long long total_tasks = (long long)N_mphi * (long long)N_lam1;
    std::atomic<long long> global_attempts(0); // Total puntos calculados (válidos + inválidos)
    std::atomic<long long> global_valid(0);    // Total puntos guardados en CSV
    std::atomic<long long> global_triple_ok(0); // [TRIPLE_OK] puntos con pos&&uni&&pert

    auto start_time       = Clock::now();
    auto last_report_time = start_time;

    const double REPORT_INTERVAL_SEC = 2.5;
    const int    BATCH_SIZE          = 200;

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
        long long local_attempts   = 0;
        long long local_valid      = 0;
        long long local_triple_ok  = 0; // [TRIPLE_OK] contador local

        std::ostringstream local_buffer;
        local_buffer.setf(std::ios::scientific);
        local_buffer << std::setprecision(6);

        #pragma omp for schedule(dynamic)
        for(int i_phi = 0; i_phi < N_mphi; ++i_phi) {
            // Set Model
            THDM model;
            SM sm; 
            model.set_SM(sm);
            double m_phi = mphi_min + i_phi * step_mphi;

            for(int i_l1 = 0; i_l1 < N_lam1; ++i_l1) {

                ++local_attempts;

                double lambda1 = lam1_min + i_l1 * step_lam1;

                double term_mass = (m_phi*m_phi * ca*ca + mh*mh * sa*sa);
                double term_l6   = 1.5 * lambda6 * tan_beta;
                double term_l7   = 0.5 * lambda7 * std::pow(tan_beta, 3);
                
                double pre_factor = (term_mass / (VEV*VEV * cb*cb)) - term_l6 + term_l7 - lambda1;

                double m12_sq = pre_factor * (VEV*VEV * cb*cb / tan_beta);

                double l2 = calc_lambda2(mh, m_phi, m12_sq,
                                         sa, ca, sb, cb,
                                         tan_beta, lambda6, lambda7, VEV);

                if (std::abs(l2) > 4 * PI) continue;

                bool pset = model.set_param_phys(
                    mh, m_phi, mA_fixed, mA_fixed,
                    sin_ba, lambda6, lambda7,
                    m12_sq, tan_beta
                );
                if (!pset) continue;

                model.set_yukawas_type(1);

                double lam1_g, lam2_g, lam3_g, lam4_g, lam5_g;
                double lam6_g, lam7_g, m12_2_g, tanb_g;
                model.get_param_gen(
                    lam1_g, lam2_g, lam3_g, lam4_g, lam5_g,
                    lam6_g, lam7_g, m12_2_g, tanb_g
                );

                Constraints check(model);
                bool pos  = check.check_positivity();
                bool uni  = check.check_unitarity();
                bool pert = check.check_perturbativity();

                // [TRIPLE_OK] contar sólo los que pasan las 3
                if (pos && uni && pert) {
                    ++local_triple_ok;
                }

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

                local_buffer
                    << m_phi   << "," << mA_fixed << "," << alpha  << "," << beta    << ","
                    << lambda6 << "," << lambda7  << "," << m12_sq << ","
                    << sin_ba  << "," << tan_beta << ","
                    << (pos  ? 1.0 : 0.0) << ","
                    << (uni  ? 1.0 : 0.0) << ","
                    << (pert ? 1.0 : 0.0) << ","
                    << w_bb   << "," << w_tautau << "," << w_WW << "," << w_ZZ << ","
                    << w_gaga << "," << w_Zga    << "," << w_gg << "," << w_hh << ","
                    << w_tot  << "," << br_gaga  << ","
                    << lambda1 << "," << lam1_g << "," << l2 << "," << lam2_g << ","
                    << lam3_g  << "," << lam4_g << "," << lam5_g
                    << "\n";

                ++local_valid;

                if (local_valid >= BATCH_SIZE) {
                    #pragma omp critical
                    {
                        results << local_buffer.str();

                        global_valid    += local_valid;
                        global_attempts += local_attempts;
                        global_triple_ok += local_triple_ok; // [TRIPLE_OK] flush al global

                        auto now  = Clock::now();
                        Duration diff = now - last_report_time;
                        
                        if (diff.count() > REPORT_INTERVAL_SEC) {
                            long long current_attempts = global_attempts.load();
                            long long current_valid    = global_valid.load();
                            long long current_triple   = global_triple_ok.load();

                            double total_elapsed = Duration(now - start_time).count();
                            
                            double pps    = (total_elapsed > 0)
                                            ? (double)current_attempts / total_elapsed
                                            : 0.0;
                            double vppm   = (total_elapsed > 0)
                                            ? (double)current_valid / (total_elapsed / 60.0)
                                            : 0.0;
                            double progress = (total_tasks > 0)
                                              ? (100.0 * current_attempts / (double)total_tasks)
                                              : 0.0;

                            std::cerr << "\r" 
                                      << "[ " << std::fixed << std::setprecision(1)
                                      << progress << "% ] "
                                      << "Speed: " << (long)pps << " pts/s | "
                                      << "ValidCSV: " << (long)vppm << "/min | "
                                      << "FoundCSV: " << current_valid << "/"
                                      << current_attempts
                                      << " | TripleOK: " << current_triple
                                      << "   " << std::flush;
                            
                            last_report_time = now;
                        }
                    }

                    local_buffer.str("");
                    local_buffer.clear();
                    local_valid      = 0;
                    local_attempts   = 0;
                    local_triple_ok  = 0; // [TRIPLE_OK] reset local
                }

            } // i_l1
        } // i_phi

        // [LIMPIEZA FINAL] Volcar remanentes al terminar el hilo
        #pragma omp critical
        {
            if (local_valid > 0) {
                results << local_buffer.str();
                global_valid += local_valid;
            }
            if (local_attempts > 0) {
                global_attempts += local_attempts;
            }
            if (local_triple_ok > 0) {
                global_triple_ok += local_triple_ok; // [TRIPLE_OK] remanente
            }
        }
    } // parallel

    auto total_time = Duration(Clock::now() - start_time).count();
    long long final_attempts  = global_attempts.load();
    long long final_valid     = global_valid.load();
    long long final_triple_ok = global_triple_ok.load(); // [TRIPLE_OK]

    std::cout << "\n\n--- Scan Finalizado ---" << endl;
    std::cout << "Total Attempts: "         << final_attempts  << endl;
    std::cout << "Total CSV Rows: "        << final_valid     << endl;
    std::cout << "Total Triple-OK Points: "<< final_triple_ok << endl;
    std::cout << "Total Time: "            << total_time      << " s" << endl;
    if (total_time > 0) {
        std::cout << "Average Speed: "
                  << (long)(final_attempts / total_time) << " pts/s" << endl;
        std::cout << "Final Valid Rate (CSV): "
                  << (long)(final_valid / (total_time/60.0)) << " rows/min" << endl;
    }

    // [TRIPLE_OK] Línea final limpia para bash / pipeline
    std::cout << "TRIPLE_OK_POINTS " << final_triple_ok << std::endl;

    results.close();
}


//
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
