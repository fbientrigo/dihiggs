/**
 * @file PhysScanWithFixings.cpp
 * @brief Escaneo de parámetros m_phi y lambda1 con parámetros físicos fijos.
 *
 * Refactored to use Analytical Inversion Strategy for m12^2 based on lambda1.
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

    double total_iters = double(N_mphi) * double(N_lam1);
    double iters_done  = 0.0;

    cout << "step_lam1=" << step_lam1 << endl;
    cout << "step_mphi=" << step_mphi << endl;

    auto start = Clock::now();
    g_bestBR = -1.0;

    const double mh = 125.0;

    // Step A: Precise Angle Initialization
    const double inv = 1.0 / std::sqrt(1.0 + tan_beta*tan_beta); // 1/sec(beta) = cos(beta) factor
    const double cb  = inv;                  // cos(beta)
    const double sb  = tan_beta * inv;       // sin(beta)
    const double cba = std::sqrt(1.0 - sin_ba*sin_ba); // cos(beta - alpha)

    // Physics basis angles needed for lambda inversion
    const double ca = cb * cba + sb * sin_ba; // cos(alpha)
    const double sa = sb * cba - cb * sin_ba; // sin(alpha)
    
    // Angles for reporting/logging
    const double beta  = std::atan(tan_beta);
    const double alpha = beta - std::asin(sin_ba);

    #pragma omp parallel for schedule(dynamic)
    for(int i_phi = 0; i_phi < N_mphi; ++i_phi) {
        vector<vector<double>> buffer;
        
        double m_phi = mphi_min + i_phi * step_mphi;

        for(int i_l1 = 0; i_l1 < N_lam1; ++i_l1) {
            double lambda1 = lam1_min + i_l1 * step_lam1;

            // Step C: Analytical Inversion (m12^2 Calculation)
            double term_mass = (m_phi*m_phi * ca*ca + mh*mh * sa*sa);
            double term_l6   = 1.5 * lambda6 * tan_beta;
            double term_l7   = 0.5 * lambda7 * std::pow(tan_beta, 3);
            double pre_factor = (term_mass / (VEV*VEV * cb*cb)) - term_l6 + term_l7 - lambda1;

            double m12_sq = pre_factor * (VEV*VEV * cb*cb / tan_beta);

            // Step 3: Validation & Model Setting
            // Calculate lambda2
            double l2 = calc_lambda2(mh, m_phi, m12_sq, sa, ca, sb, cb, tan_beta, lambda6, lambda7, VEV);

            // Filter
            if (std::abs(l2) > 4 * PI) continue;

            // Set Model
            THDM model;
            SM sm; model.set_SM(sm);
            bool pset = model.set_param_phys(mh, m_phi, mA_fixed, mA_fixed, sin_ba, lambda6, lambda7, m12_sq, tan_beta);
            
            if (!pset) continue;

            model.set_yukawas_type(1);

            // Retrieve generated parameters for logging
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

            // Guardar fila
            buffer.push_back(std::vector<double>{
                m_phi, mA_fixed, alpha, beta, lambda6, lambda7, m12_sq,
                sin_ba, tan_beta,
                pos?1.0:0.0, uni?1.0:0.0, pert?1.0:0.0,
                w_bb, w_tautau, w_WW, w_ZZ,
                w_gaga, w_Zga, w_gg, w_hh,
                w_tot, br_gaga, lambda1, lam1_g, l2, lam2_g, lam3_g, lam4_g, lam5_g
            });
        }

        #pragma omp critical
        {
            double n = static_cast<double>(buffer.size());
            for (auto &row : buffer) {
                write_csv_row(results, row);
            }
            iters_done += n;
            buffer.clear();

            auto now = high_resolution_clock::now();
            double elapsed = duration<double>(now - start).count();
            
            std::cerr << fixed << setprecision(1)
                    << "Rows: " << iters_done
                    << "  elapsed: " << elapsed << " s\r";
        }
    }

    results.close();
    cout << "\n\nEscaneo completo." << endl;
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
