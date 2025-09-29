/**
 * @file PhysParamScan.cpp
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
 *   ./PhysScanClassic <output_csv> <N_mphi> [N_m12] [yukawa_type]
 *
 * Paralelización:
 *   Usa OpenMP para distribuir el bucle 2D de (m_phi × m12).
 *
 * Dependencias:
 *   - Biblioteca 2HDMC ( –l2HDMC )
 *   - GSL ( –lgsl –lgslcblas –lm )
 *   - OpenMP
 *
 * Autor: Fabian Trigo
 * Fecha: septiembre 2025
 */

// activa el chequeo acelerado de lam1, lam2
#define SPEED_TEST

#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"  // CSV & config utilities
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <chrono>

constexpr double PI = std::acos(-1.0);
constexpr double VEV = 246.0;  // Valor estándar
constexpr double LMAX = 4.0 * M_PI; // Límite de perturbatividad


using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;
using namespace std;
using namespace std::chrono;

struct FixedParams {
    double mA;
    double sin_ba;
    double tan_beta;
    double lambda6;
    double lambda7;
};


// ------- m12 base ---------
std::pair<double,double>
find_m12_range(double m12_guess,
               double mh, double mphi,
               double sin_ba, double tanb,
               double lambda6, double lambda7)
{
    auto valid = [&](double x){
        double l1 = calc_lambda1(mh, mphi, x, sin_ba, tanb, lambda6, lambda7);
        double l2 = calc_lambda2(mh, mphi, x, sin_ba, tanb, lambda6, lambda7);
        return (l1>0.0 && l1<LMAX && l2>0.0 && l2<LMAX);
    };

    // --- Bracket inferior ---
    double low_ok = m12_guess;
    double low_bad = std::max(0.0, m12_guess);
    while (valid(low_bad) && low_bad>0.0) {
        low_ok  = low_bad;
        low_bad *= 0.5;              // exponencial hacia abajo
    }

    // --- Bracket superior ---
    double up_ok = m12_guess;
    double up_bad = m12_guess;
    while (valid(up_bad)) {
        up_ok = up_bad;
        up_bad *= 2.0;               // exponencial hacia arriba
    }

    // --- Bisección inferior ---
    while ((low_ok - low_bad) > TOL) {
        double mid = 0.5 * (low_ok + low_bad);
        if (valid(mid)) low_ok = mid; else low_bad = mid;
    }

    // --- Bisección superior ---
    while ((up_bad - up_ok) > TOL) {
        double mid = 0.5 * (up_ok + up_bad);
        if (valid(mid)) up_ok = mid; else up_bad = mid;
    }

    return {low_ok, up_ok};
}



// --------------------------- main scan-----------------------------
void perform_param_scan(
    int y_type, double mphi_min, double mphi_max,
    int N_mphi, int N_m12, int delta_m12_exp,
    double mA_fixed, double sin_ba, double tan_beta,
    double lambda6, double lambda7,
    const string &output_file
) {
    ofstream results(output_file);
    if (!results.is_open()) {
        cerr << "No pude abrir: " << output_file << endl;
        return;
    }

    vector<string> columns = {
        "m_phi","mA","alpha","beta","lambda6","lambda7","m12",
        "sin_ba","tan_beta","positivity_ok","unitarity_ok","perturbativity_ok",
        "width_bb","width_tautau","width_WW","width_ZZ",
        "width_gaga","width_Zga","width_gg","width_hh",
        "total_width","br_gaga",
        "lam1","computed_lam1","lam2","computed_lam2","lam3","lam4","lam5"
    };
    write_csv_header(results, columns);

    double step_mphi = (N_mphi > 1) ? (mphi_max - mphi_min)/(N_mphi - 1) : 0.0;
    double total_iters = double(N_mphi) * double(N_m12);
    int iters_done  = 0;
    int iters_jumped = 0;

    cout << "step_mphi="<<step_mphi<<endl;
    auto start = Clock::now();

    double mh = 125.0;

    double inv = 1.0/std::sqrt(1.0 + tan_beta*tan_beta);
    double cb  = inv;
    double sb  = tan_beta * inv;
    double cba = std::sqrt(1.0 - sin_ba*sin_ba);

    double ca = cb * cba + sb * sin_ba;
    double sa = sb * cba - cb * sin_ba;
    double beta  = std::atan(tan_beta);
    double alpha = beta - std::asin(sin_ba);

    #pragma omp parallel for schedule(dynamic)
    for(int i_phi = 0; i_phi < N_mphi; ++i_phi) {
        double m_phi = mphi_min + i_phi * step_mphi;
        for(int i_m12 = 0; i_m12 < N_m12; ++i_m12) {

            vector<vector<double>> buffer;

            double m12_base = m12_base(mh, m_phi, sa, ca, sb, cb, 
                lambda6, lambda7,
                tan_beta, VEV);

            double delta_m12 = std::pow(10.0, - delta_m12_exp);
            double m12 = m12_base - 0.5 * N_m12 * delta_m12 + i_m12 * delta_m12;

            double l1 = calc_lambda1(mh, m_phi, m12, sa, ca, sb, cb, tan_beta, lambda6, lambda7, VEV);
            double l2 = calc_lambda2(mh, m_phi, m12, sa, ca, sb, cb, tan_beta, lambda6, lambda7, VEV);

            //#ifdef SPEED_TEST
            if (!check_lambda(l1) || !check_lambda(l2)) {
                iters_jumped++;
                continue;
            }
            //#endif

            THDM model;
            SM sm; model.set_SM(sm);
            bool ok = model.set_param_phys(
                mh, m_phi,
                mA_fixed, mA_fixed,
                sin_ba, lambda6, lambda7, m12, tan_beta
            );
            model.set_yukawas_type(y_type);

            if (!ok) continue;

            Constraints check(model);
            bool pos  = check.check_positivity();
            bool uni  = check.check_unitarity();
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



            double lam1, lam2, lam3, lam4, lam5, lam6_g, lam7_g, m12_2_g, tanb_g;
            model.get_param_gen( lam1, lam2, lam3, lam4, lam5, lam6_g, lam7_g, m12_2_g, tanb_g);

            buffer.push_back({
                m_phi, mA_fixed, alpha, beta, lambda6, lambda7, m12,
                sin_ba, tan_beta,
                pos?1.0:0.0, uni?1.0:0.0, pert?1.0:0.0,
                w_bb, w_tautau, w_WW, w_ZZ,
                w_gaga, w_Zga, w_gg, w_hh,
                w_tot, br_gaga, lam1, l1, lam2, l2, lam3, lam4, lam5
            });

            #pragma omp critical
            {
                for (auto &row : buffer) write_csv_row(results, row);
                iters_done += buffer.size();
                buffer.clear();

                auto now = high_resolution_clock::now();
                double elapsed = duration<double>(now - start).count();
                double prog = (iters_done / total_iters) * 100.0;
                double eta = (iters_done > 0.0) ? (elapsed / iters_done) * (total_iters - iters_done) : 0.0;

                std::cerr << fixed << setprecision(1)
                          << "Iter: " << iters_done << "/" << total_iters
                          << " (" << prog << "%)  elapsed: " << elapsed
                          << " s  ETA: " << eta << " s  Omited: " << iters_jumped << "\r";
            }
        }
    }

    results.close();
    cout << "\n\nEscaneo completo.\n";

    #ifdef SPEED_TEST
    auto end = Clock::now();
    double secs = Duration(end - start).count();
    std::cerr << "== Speed test: tiempo total de escaneo = " << secs << " s ==\n";
    #endif
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Uso: " << argv[0]
                  << " <output_csv> <N_mphi> [delta_m12_exp] [N_m12] [yukawa_type]\n";
        return 1;
    }

    std::string output_csv = argv[1];

    // Obligatorio
    int N_mphi = atoi(argv[2]);
    if (N_mphi < 1) {
        std::cerr << "N_mphi debe ser positivo.\n";
        return 1;
    }

    // Defaults
    int delta_m12_exp = 13;  // default → paso ~10^-13
    int N_m12  = 10;
    int y_type = 1;

    if (argc > 3) {
        delta_m12_exp = atoi(argv[3]);
        if (delta_m12_exp < 1) {
            std::cerr << "delta_m12_exp debe ser >= 1.\n";
            return 1;
        }
    }

    if (argc > 4) {
        N_m12 = atoi(argv[4]);
        if (N_m12 < 1) {
            std::cerr << "N_m12 debe ser positivo.\n";
            return 1;
        }
    }

    if (argc > 5) {
        y_type = atoi(argv[5]);
        if (y_type < 1 || y_type > 4) {
            std::cerr << "yukawa_type debe estar entre 1 y 4.\n";
            return 1;
        }
    }

    FixedParams fp;
    fp.mA       = 300.0;
    fp.sin_ba   = 0.99;
    fp.tan_beta = 10.0;
    fp.lambda6  = 0.1;
    fp.lambda7  = 0.0;

    perform_param_scan(
        y_type,
        130.0, 300.0,
        N_mphi, N_m12, delta_m12_exp,
        fp.mA, fp.sin_ba, fp.tan_beta,
        fp.lambda6, fp.lambda7,
        output_csv
    );

    return 0;
}
