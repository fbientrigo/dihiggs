/**
 * @file PhysScanWithFixings.cpp
 * @brief Escaneo de parámetros m_phi y m12 con parámetros físicos fijos.
 *
 * Este programa recorre un rango de valores de m_phi y m12, mientras
 * mantiene los demás parámetros del 2HDM (mA, sin(b−a), tanβ, λ6, λ7) fijos.
 * Para cada punto calcula condiciones de positividad, unitariedad y perturbatividad,
 * construye la tabla de decaimientos y guarda en CSV:
 *   m_phi, mA, α, β, λ6, λ7, m12,
 *   sin(b−a), tanβ,
 *   positivity_ok, unitarity_ok, perturbativity_ok,
 *   widths de bb, ττ, WW, ZZ, γγ, Zγ, gg, hh,
 *   total_width, BR(h→γγ).
 *
 * Uso:
 *   ./PhysScanWithFixings
 *     <mphi_min> <mphi_max> <N_mphi>
 *     <m12_min>  <m12_max>  <N_m12>
 *     <mA_fixed> <sin(b-a)> <tan_beta>
 *     <lambda6>   <lambda7>
 *     <output_csv>
 *
 * Ejemplo:
 *   ./PhysScanWithFixings 150 300 30 0.5 3.0 30000 300 0.9 1000 0.1 0 results.csv
 *
 * Paralelización:
 *   Usa OpenMP para distribuir el bucle 2D de (m_phi × m12) con
 *   `#pragma omp for collapse(2) schedule(dynamic)`.
 *
 * Dependencias:
 *   - Biblioteca 2HDMC (–l2HDMC)
 *   - GSL (–lgsl –lgslcblas –lm)
 *   - OpenMP
 *
 * Autor: Fabian Trigo
 * Fecha: junio 2025
 */


#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"  // CSV & config utilities
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <omp.h>

using namespace std;
using namespace std::chrono;

// Estructura para resultados
struct ParamSet {
    double m_phi, mA, alpha, beta, lambda6, lambda7, m12;
};

static ParamSet g_bestParams;
static double   g_bestBR = -1.0;

void perform_param_scan_fixings(
    double mphi_min, double mphi_max, int N_mphi,
    double m12_min,  double m12_max,  int N_m12,
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
        "total_width","br_gaga"
    };
    write_csv_header(results, columns);

    // Calcular pasos
    double step_mphi = (N_mphi > 1) ? (mphi_max - mphi_min)/(N_mphi - 1) : 0.0;
    double step_m12  = (N_m12  > 1) ? (m12_max  - m12_min )/(N_m12  - 1) : 0.0;

    double total_iters = double(N_mphi) * double(N_m12);
    double iters_done  = 0.0;
    auto start = high_resolution_clock::now();
    g_bestBR = -1.0;

    #pragma omp parallel
    {
        vector<vector<double>> buffer;
        double local_bestBR = -1.0;
        ParamSet local_bestParams;

        #pragma omp for collapse(2) schedule(dynamic)
        for(int i_phi = 0; i_phi < N_mphi; ++i_phi) {
            for(int i_m12 = 0; i_m12 < N_m12; ++i_m12) {
                double m_phi = mphi_min + i_phi * step_mphi;
                double m12   = m12_min  + i_m12 * step_m12;

                // Construir modelo
                THDM model;
                SM sm; model.set_SM(sm);
                bool ok = model.set_param_phys(
                    125.0,         // m_h fijo
                    m_phi,
                    mA_fixed,
                    mA_fixed,     // mA = mHp
                    sin_ba,
                    lambda6,
                    lambda7,
                    m12,
                    tan_beta
                );
                if (!ok) continue;

                Constraints check(model);
                bool pos = check.check_positivity();
                bool uni = check.check_unitarity();
                bool pert = check.check_perturbativity();

                DecayTable tab(model);
                double w_bb     = tab.get_gamma_hdd(2,3,3);
                double w_tautau= tab.get_gamma_hll(2,3,3);
                double w_WW    = tab.get_gamma_hvv(2,3);
                double w_ZZ    = tab.get_gamma_hvv(2,2);
                double w_gaga  = tab.get_gamma_hgaga(2);
                double w_Zga   = tab.get_gamma_hZga(2);
                double w_gg    = tab.get_gamma_hgg(2);
                double w_hh    = tab.get_gamma_hhh(2,1,1);
                double w_tot   = tab.get_gammatot_h(2);
                double br_gaga = (w_tot > 1e-12) ? w_gaga / w_tot : 0.0;

                // Guardar fila
                buffer.push_back({
                    m_phi, mA_fixed, 0.0, 0.0, lambda6, lambda7, m12,
                    sin_ba, tan_beta,
                    pos?1.0:0.0, uni?1.0:0.0, pert?1.0:0.0,
                    w_bb, w_tautau, w_WW, w_ZZ,
                    w_gaga, w_Zga, w_gg, w_hh,
                    w_tot, br_gaga
                });

                // Actualizar mejor BR local
                if (pos && w_tot>0.0 && br_gaga > local_bestBR) {
                    local_bestBR = br_gaga;
                    local_bestParams = {m_phi, mA_fixed, 0., 0., lambda6, lambda7, m12};
                }
            }
        }

        // Escribir buffer y actualizar global
        #pragma omp critical
        {
            for (auto &row : buffer) write_csv_row(results, row);
            buffer.clear();
            if (local_bestBR > g_bestBR) {
                g_bestBR = local_bestBR;
                g_bestParams = local_bestParams;
            }
            iters_done += buffer.size();
            double elapsed = duration<double>(high_resolution_clock::now()-start).count();
            double prog = (iters_done/total_iters)*100.0;
            double eta  = (elapsed/iters_done)*(total_iters-iters_done);
            cout << fixed<<setprecision(1)
                 << "Prog: "<<prog<<"%  ETA: "<<eta/60.0<<" min\r";
        }
    } // end parallel

    results.close();
    cout << "\n\nEscaneo completo. Best br_gaga="<<g_bestBR<<endl;
}

int main(int argc, char* argv[]) {
    if (argc != 12) {
        cerr << "Uso: "<< argv[0]
             << " mphi_min mphi_max N_mphi m12_min m12_max N_m12"
                " mA sin(b-a) tan(beta) lambda6 lambda7 output.csv\n";
        return 1;
    }
    double mphi_min = atof(argv[1]);
    double mphi_max = atof(argv[2]);
    int    N_mphi   = atoi(argv[3]);
    double m12_min  = atof(argv[4]);
    double m12_max  = atof(argv[5]);
    int    N_m12    = atoi(argv[6]);
    double mA_fixed = atof(argv[7]);
    double sin_ba   = atof(argv[8]);
    double tan_beta = atof(argv[9]);
    double lambda6  = atof(argv[10]);
    double lambda7  = atof(argv[11]);
    string output   = argv[12];

    perform_param_scan_fixings(
        mphi_min, mphi_max, N_mphi,
        m12_min,  m12_max,  N_m12,
        mA_fixed, sin_ba,   tan_beta,
        lambda6,  lambda7,
        output
    );
    return 0;
}
