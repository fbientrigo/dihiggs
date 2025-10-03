/** 
 * @file PhysParamScan.cpp
 * @brief Escaneo de parámetros m_phi y m12 con búsqueda **lineal** desde m12_base.
 *
 * Estrategia simple y trazable:
 *  - Tomamos m12_guess = m12_base(...).
 *  - Definimos un paso físico seguro: Δ = 2π v^2 (cotβ cos^2β + tanβ sin^2β).
 *  - Si el guess es válido (0 < λ1,λ2 < 4π), avanzamos **de a poco** ±Δ·f hasta que deje de cumplir
 *    y tomamos esos extremos como [m12_min, m12_max]. Si el guess no es válido, buscamos una semilla
 *    válida en un radio ±Δ.
 *  - Luego muestreamos N_m12 valores uniformes en [m12_min, m12_max] y evaluamos el modelo.
 *
 * Ventajas: código directo, sin bisección; excelente para depuración con logs.
 *
 * Uso:
 *   ./PhysParamScan <output_csv> <N_mphi> [tol_exp] [N_m12] [yukawa_type]
 *   - tol_exp: precisión de reporte (no aplica a la búsqueda lineal; default 3 → 10^-3)
 *   - N_m12: puntos a muestrear dentro del rango encontrado (default 10)
 *   - yukawa_type: 1–4 (default 1)
 *
 * Dependencias:
 *   - 2HDMC ( -l2HDMC ), GSL ( -lgsl -lgslcblas -lm ), OpenMP (opcional)
 *
 * Autor: Fabian Trigo — versión lineal simplificada (oct-2025)
 */

// Activa chequeo rápido de λ1,λ2 antes de construir el modelo completo
#define SPEED_TEST 1
// Activa logs de diagnóstico locales (0/1)
#ifndef RLOG
#define RLOG 0
#endif

#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"  // CSV & config utilities: write_csv_header, write_csv_row
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <chrono>
#ifdef _OPENMP
  #include <omp.h>
#endif

using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;
using namespace std;
using namespace std::chrono;

// ===== Constantes físicas =====
constexpr double PI   = 3.14159265358979323846;
constexpr double VEV  = 246.0;            // GeV
constexpr double LMAX = 4.0 * PI;         // 4π — límite de perturbatividad

// ===== Parámetros fijos del escaneo =====
struct FixedParams {
    double mA;
    double sin_ba;
    double tan_beta;
    double lambda6;
    double lambda7;
};

// ===== Prototipos (debes tener estas funciones en tu proyecto) =====
// calc_lambda1/2 deben existir con esta firma (ajusta si tu firma difiere)
double calc_lambda1(double mh, double mphi, double m12,
                    double sa, double ca, double sb, double cb,
                    double tanb, double lambda6, double lambda7, double v);

double calc_lambda2(double mh, double mphi, double m12,
                    double sa, double ca, double sb, double cb,
                    double tanb, double lambda6, double lambda7, double v);

inline bool check_lambda(double l){ return (l > 0.0 && l < LMAX); }

// ===== Paso físico "seguro" Δ =====
inline double delta_safe(double /*sa*/, double /*ca*/, double sb, double cb, double tanb){
    // uso de la expresion de l1
    // Δ = 2π v^2 (cotβ cos^2β + tanβ sin^2β)
    const double cotb = 1.0 / tanb;
    const double term = cotb * cb*cb; // + tanb * sb*sb;
    return 4.0 * PI * VEV * VEV * term;  // GeV^2
}

// ===== Validación de un punto dado m12 =====
struct Lpair { double l1, l2; };

inline bool valid_point(double mh, double mphi, double m12,
                        double sa, double ca, double sb, double cb,
                        double tanb, double lambda6, double lambda7,
                        Lpair* out = nullptr)
{
    if (!std::isfinite(m12)) return false;
    const double l1 = calc_lambda1(mh, mphi, m12, sa, ca, sb, cb, tanb, lambda6, lambda7, VEV);
    const double l2 = calc_lambda2(mh, mphi, m12, sa, ca, sb, cb, tanb, lambda6, lambda7, VEV);
    if (out) { out->l1 = l1; out->l2 = l2; }
    return check_lambda(l1) && check_lambda(l2);
}

// ===== Búsqueda lineal de rango [m12_min, m12_max] desde m12_guess =====
pair<double,double>
find_m12_range_linear(double m12_guess,
                      double mh, double mphi,
                      double sa, double ca, double sb, double cb,
                      double tanb, double lambda6, double lambda7,
                      double step_frac = 0.25,
                      int    max_steps_each_dir = 2048)
{
    constexpr double M12_ABS_MAX = 1.0e8; // GeV^2 — cota dura

    // Paso base físico
    double d_safe = delta_safe(0,0,sb,cb,tanb);
    if (!std::isfinite(d_safe) || d_safe <= 0) d_safe = 1.0; // fallback

    const double step = step_frac * d_safe / max_steps_each_dir ;

    // 0) si el guess no es válido, probamos un radio local ±step
#if RLOG
    cerr << "[RANGE] mphi="<<mphi<<"  guess="<<m12_guess
         <<"  l1="<<gg.l1<<" l2="<<gg.l2
         <<"  d_safe="<<d_safe<<"  step="<<step<<"\n";
#endif


    // 1) Descenso lineal
    double m12_min = m12_guess;
    for (int k=0; k<max_steps_each_dir; ++k){
        double cand = m12_min - step;
        if (cand < 0.0 || cand > M12_ABS_MAX) break;
        if (!valid_point(mh, mphi, cand, sa, ca, sb, cb, tanb, lambda6, lambda7)) break;
        m12_min = cand;
    }

    // 2) Ascenso lineal
    double m12_max = m12_guess;

#if RLOG
    Lpair llo{}, lup{};
    valid_point(mh,mphi,m12_min,sa,ca,sb,cb,tanb,lambda6,lambda7,&llo);
    valid_point(mh,mphi,m12_max,sa,ca,sb,cb,tanb,lambda6,lambda7,&lup);
    cerr << "[RANGE] found: ["<<m12_min<<","<<m12_max<<"]  | λ(min)=("<<llo.l1<<","<<llo.l2
         <<") λ(max)=("<<lup.l1<<","<<lup.l2<<")\n";
#endif

    if (m12_max <= m12_min) return {0.0, 0.0};
    return {m12_min, m12_max};
}

// ===================== ESCANEO PRINCIPAL =====================
void perform_param_scan(
    int y_type, double mphi_min, double mphi_max,
    int N_mphi, int N_m12, double /*TOL_unused_for_linear*/,
    double mA_fixed, double sin_ba, double tan_beta,
    double lambda6, double lambda7,
    const string &output_file)
{
    ofstream results(output_file);
    if (!results.is_open()) {
        cerr << "No pude abrir: " << output_file << "\n";
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

    const double step_mphi = (N_mphi > 1) ? (mphi_max - mphi_min)/(N_mphi - 1) : 0.0;
    const double total_iters = double(N_mphi) * double(N_m12);
    long long iters_done  = 0;
    long long iters_jumped = 0;

    cout << "step_mphi="<<step_mphi<<"\n";
    const auto t0 = Clock::now();

    const double mh = 125.0;

    // Ángulos derivados
    const double inv = 1.0 / std::sqrt(1.0 + tan_beta*tan_beta);
    const double cb  = inv;
    const double sb  = tan_beta * inv;
    const double cba = std::sqrt(1.0 - sin_ba*sin_ba);

    const double ca = cb * cba + sb * sin_ba;
    const double sa = sb * cba - cb * sin_ba;
    const double beta  = std::atan(tan_beta);
    const double alpha = beta - std::asin(sin_ba);

    THDM model; SM sm; model.set_SM(sm);

    // Bucle externo (puedes reactivar OpenMP cuando quieras)
    for(int i_phi = 0; i_phi < N_mphi; ++i_phi) {
        const double m_phi = mphi_min + i_phi * step_mphi;

        // Semilla física
        const double m12_guess = m12_base(mh, m_phi, sa, ca, sb, cb, lambda6, lambda7, tan_beta, VEV);

        // Rango lineal desde la semilla
        auto [m12_min, m12_max] = find_m12_range_linear(
            m12_guess, mh, m_phi, sa, ca, sb, cb, tan_beta, lambda6, lambda7,
            /*step_frac=*/0.25, /*max_steps_each_dir=*/4096);

        if (m12_max <= m12_min) {
            iters_jumped += N_m12;
            cerr << "No valid m12 range at m_phi=" << m_phi << "\n";
            continue;
        }

        for (int i_m12 = 0; i_m12 < N_m12; ++i_m12) {
            const double m12 = m12_min + (m12_max - m12_min) * ( (N_m12==1)? 0.0 : (double)i_m12 / (N_m12 - 1) );

            // Pre-filtro λ1,λ2 (rápido)
            const double l1 = calc_lambda1(mh, m_phi, m12, sa, ca, sb, cb, tan_beta, lambda6, lambda7, VEV);
            const double l2 = calc_lambda2(mh, m_phi, m12, sa, ca, sb, cb, tan_beta, lambda6, lambda7, VEV);
            if (!check_lambda(l1) || !check_lambda(l2)) { iters_jumped++; continue; }

            // Modelo completo
            model.set_yukawas_type(y_type);
            const bool ok = model.set_param_phys(mh, m_phi, mA_fixed, mA_fixed, sin_ba, lambda6, lambda7, m12, tan_beta);
            if (!ok) { iters_jumped++; continue; }

            Constraints check(model);
            const bool pos  = check.check_positivity();
            const bool uni  = check.check_unitarity();
            const bool pert = check.check_perturbativity();

            DecayTable tab(model);
            const double w_bb     = tab.get_gamma_hdd(2,3,3);
            const double w_tautau = tab.get_gamma_hll(2,3,3);
            const double w_WW     = tab.get_gamma_hvv(2,3);
            const double w_ZZ     = tab.get_gamma_hvv(2,2);
            const double w_gaga   = tab.get_gamma_hgaga(2);
            const double w_Zga    = tab.get_gamma_hZga(2);
            const double w_gg     = tab.get_gamma_hgg(2);
            const double w_hh     = tab.get_gamma_hhh(2,1,1);
            const double w_tot    = tab.get_gammatot_h(2);
            const double br_gaga  = (w_tot > 1e-15) ? w_gaga / w_tot : 0.0;

            double lam1, lam2, lam3, lam4, lam5, lam6_g, lam7_g, m12_2_g, tanb_g;
            model.get_param_gen(lam1, lam2, lam3, lam4, lam5, lam6_g, lam7_g, m12_2_g, tanb_g);

            vector<double> row = {
                m_phi, mA_fixed, alpha, beta, lambda6, lambda7, m12,
                sin_ba, tan_beta,
                pos?1.0:0.0, uni?1.0:0.0, pert?1.0:0.0,
                w_bb, w_tautau, w_WW, w_ZZ,
                w_gaga, w_Zga, w_gg, w_hh,
                w_tot, br_gaga,
                lam1, l1, lam2, l2, lam3, lam4, lam5
            };
            write_csv_row(results, row);
            ++iters_done;

            if ((iters_done % 2000) == 0) {
                const auto now = Clock::now();
                const double elapsed = Duration(now - t0).count();
                const double prog = (iters_done / total_iters) * 100.0;
                const double eta  = (iters_done > 0) ? (elapsed / iters_done) * (total_iters - iters_done) : 0.0;
                cerr << fixed << setprecision(1)
                     << "Iter: " << iters_done << "/" << total_iters
                     << " (" << prog << "%)  elapsed: " << elapsed
                     << " s  ETA: " << eta << " s  Omitted: " << iters_jumped << "\r";
            }
        }
    }

    results.close();
    cout << "\n\nEscaneo completo.\n";
}


int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Uso: " << argv[0]
                  << " <output_csv> <N_mphi> <tan_beta> [tol_exp] [N_m12] [yukawa_type]\n "
                  << "  - tan_beta: valor obligatorio\n "
                  << "  - tol_exp: precisión de reporte (default 3 → 10^-3)\n "
                  << "  - N_m12: número de puntos en m12 (default 10)\n "
                  << "  - yukawa_type: 1–4 (default 1)\n";
        return 1;
    }

    const std::string output_csv = argv[1];

    // Obligatorios
    const int N_mphi = std::atoi(argv[2]);
    if (N_mphi < 1) { std::cerr << "N_mphi debe ser positivo.\n"; return 1; }

    const double tan_beta = std::atof(argv[3]);
    if (tan_beta <= 0.0) { std::cerr << "tan_beta debe ser positivo.\n"; return 1; }

    // Defaults
    int TOL_exp = 3;  // solo informativo en esta versión lineal
    int N_m12   = 10;
    int y_type  = 1;

    if (argc > 4) {
        TOL_exp = std::atoi(argv[4]);
        if (TOL_exp < 1) { std::cerr << "TOL_exp debe ser >= 1.\n"; return 1; }
    }
    if (argc > 5) {
        N_m12 = std::atoi(argv[5]);
        if (N_m12 < 1) { std::cerr << "N_m12 debe ser positivo.\n"; return 1; }
    }
    if (argc > 6) {
        y_type = std::atoi(argv[6]);
        if (y_type < 1 || y_type > 4) { std::cerr << "yukawa_type debe estar entre 1 y 4.\n"; return 1; }
    }

    FixedParams fp{ /*mA*/300.0, /*sin_ba*/0.999, tan_beta, /*lambda6*/0.1, /*lambda7*/0.0 };

    // Nota: TOL no se usa en la búsqueda lineal, se mantiene por compatibilidad de firma
    const double TOL = std::pow(10.0, -TOL_exp);

    perform_param_scan(
        y_type,
        /*mphi_min*/130.0, /*mphi_max*/300.0,
        N_mphi, N_m12, TOL,
        fp.mA, fp.sin_ba, fp.tan_beta,
        fp.lambda6, fp.lambda7,
        output_csv
    );

    return 0;
}
