/** 
 * @file PhysParamScan.cpp
 * @brief Escaneo de parámetros m_phi y m12 con búsqueda **lineal** desde m12_base.
 * (… descripción original …)
 */

#define SPEED_TEST 1
#ifndef RLOG
#define RLOG 0
#endif

#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <thread>
#include <map>
#include <cstring>
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

// ===== Prototipos =====
double calc_lambda1(double mh, double mphi, double m12,
                    double sa, double ca, double sb, double cb,
                    double tanb, double lambda6, double lambda7, double v);
double calc_lambda2(double mh, double mphi, double m12,
                    double sa, double ca, double sb, double cb,
                    double tanb, double lambda6, double lambda7, double v);
inline bool check_lambda(double l){ return (l > 0.0 && l < LMAX); }

// ===== Paso físico "seguro" Δ =====
inline double delta_safe(double /*sa*/, double /*ca*/, double sb, double cb, double tanb){
    const double cotb = 1.0 / tanb;
    const double term = cotb * cb*cb; // + tanb * sb*sb;
    return 4.0 * PI * VEV * VEV * term;  // GeV^2
}

// ===== Validación =====
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

// ===== Búsqueda lineal de rango [m12_min, m12_max] =====
pair<double,double>
find_m12_range_linear(double m12_guess,
                      double mh, double mphi,
                      double sa, double ca, double sb, double cb,
                      double tanb, double lambda6, double lambda7,
                      double step_frac = 0.25,
                      int    max_steps_each_dir = 2048)
{
    constexpr double M12_ABS_MAX = 1.0e8; // GeV^2 — cota dura
    double d_safe = delta_safe(0,0,sb,cb,tanb);
    if (!std::isfinite(d_safe) || d_safe <= 0) d_safe = 1.0;
    const double step = step_frac * d_safe / max_steps_each_dir ;

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
    for (int k=0; k<max_steps_each_dir; ++k){
        double cand = m12_max + step;
        if (cand < 0.0 || cand > M12_ABS_MAX) break;
        if (!valid_point(mh, mphi, cand, sa, ca, sb, cb, tanb, lambda6, lambda7)) break;
        m12_max = cand;
    }

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

    for(int i_phi = 0; i_phi < N_mphi; ++i_phi) {
        const double m_phi = mphi_min + i_phi * step_mphi;

        const double m12_guess = m12_base(mh, m_phi, sa, ca, sb, cb, lambda6, lambda7, tan_beta, VEV);

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

            const double l1 = calc_lambda1(mh, m_phi, m12, sa, ca, sb, cb, tan_beta, lambda6, lambda7, VEV);
            const double l2 = calc_lambda2(mh, m_phi, m12, sa, ca, sb, cb, tan_beta, lambda6, lambda7, VEV);
            if (!check_lambda(l1) || !check_lambda(l2)) { iters_jumped++; continue; }

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

/* ===================== CLI con flags + posicional ===================== */

struct CLI {
    // valores (con defaults del programa original)
    string  output_csv;
    int     N_mphi      = -1;
    double  tan_beta    = -1.0;
    int     TOL_exp     = 3;
    int     N_m12       = 10;
    int     y_type      = 1;
    // extras opcionales
    double  mphi_min    = 130.0;
    double  mphi_max    = 300.0;
    double  mA          = 300.0;
    double  sin_ba      = 0.999;
    double  lambda6     = 0.1;
    double  lambda7     = 0.0;
    int     threads     = 0;    // 0 => no fijar
};

// parsea --key value o --key=value
static map<string,string> parse_kv_args(int argc, char* argv[]) {
    map<string,string> kv;
    for (int i=1;i<argc;++i){
        const char* s = argv[i];
        if (std::strncmp(s,"--",2)!=0) continue;
        const char* eq = std::strchr(s,'=');
        if (eq){
            string k(s+2, eq-(s+2));
            string v(eq+1);
            kv[k]=v;
        }else{
            string k(s+2);
            if (i+1<argc && std::strncmp(argv[i+1],"--",2)!=0){
                kv[k]=argv[i+1];
                ++i;
            }else{
                kv[k]="1"; // flag booleana (no usamos aquí, pero queda soportado)
            }
        }
    }
    return kv;
}

static bool to_int(const map<string,string>& kv, const string& k, int& out){
    auto it=kv.find(k); if(it==kv.end()) return false; out=std::stoi(it->second); return true;
}
static bool to_double(const map<string,string>& kv, const string& k, double& out){
    auto it=kv.find(k); if(it==kv.end()) return false; out=std::stod(it->second); return true;
}
static bool to_string(const map<string,string>& kv, const string& k, string& out){
    auto it=kv.find(k); if(it==kv.end()) return false; out=it->second; return true;
}

static void print_usage(const char* prog){
    cerr <<
    "Uso (posicional):\n  " << prog <<
    " <output_csv> <N_mphi> <tan_beta> [tol_exp] [N_m12] [yukawa_type]\n\n"
    "Flags equivalentes (opcionales, tienen prioridad):\n"
    "  --output FILE         \n"
    "  --N_mphi INT          \n"
    "  --tan_beta FLOAT      \n"
    "  --tol_exp INT         \n"
    "  --N_m12 INT           \n"
    "  --yukawa_type {1..4}  \n"
    "Extras opcionales:\n"
    "  --mphi_min FLOAT   (def 130)\n"
    "  --mphi_max FLOAT   (def 300)\n"
    "  --mA FLOAT         (def 300)\n"
    "  --sin_ba FLOAT     (def 0.999)\n"
    "  --lambda6 FLOAT    (def 0.1)\n"
    "  --lambda7 FLOAT    (def 0.0)\n"
    "  --threads INT      (fija hilos OpenMP si está habilitado)\n";
}

/* ===================== main ===================== */
int main(int argc, char* argv[]) {
    CLI cli;
    // 1) Posicional (como antes)
    if (argc >= 4 && strncmp(argv[1],"--",2)!=0) {
        cli.output_csv = argv[1];
        cli.N_mphi     = std::atoi(argv[2]);
        cli.tan_beta   = std::atof(argv[3]);
        if (argc > 4 && strncmp(argv[4],"--",2)!=0) {
            cli.TOL_exp = std::atoi(argv[4]);
            if (cli.TOL_exp < 1) { std::cerr << "TOL_exp debe ser >= 1.\n"; return 1; }
        }
        if (argc > 5 && strncmp(argv[5],"--",2)!=0) {
            cli.N_m12 = std::atoi(argv[5]);
        }
        if (argc > 6 && strncmp(argv[6],"--",2)!=0) {
            cli.y_type = std::atoi(argv[6]);
        }
    }

    // 2) Flags (tienen prioridad sobre posicional)
    auto kv = parse_kv_args(argc, argv);
    to_string(kv, "output",      cli.output_csv);
    to_int   (kv, "N_mphi",      cli.N_mphi);
    to_double(kv, "tan_beta",    cli.tan_beta);
    to_int   (kv, "tol_exp",     cli.TOL_exp);
    to_int   (kv, "N_m12",       cli.N_m12);
    to_int   (kv, "yukawa_type", cli.y_type);

    to_double(kv, "mphi_min",    cli.mphi_min);
    to_double(kv, "mphi_max",    cli.mphi_max);
    to_double(kv, "mA",          cli.mA);
    to_double(kv, "sin_ba",      cli.sin_ba);
    to_double(kv, "lambda6",     cli.lambda6);
    to_double(kv, "lambda7",     cli.lambda7);
    to_int   (kv, "threads",     cli.threads);

    // 3) Validaciones mínimas
    if (cli.output_csv.empty() || cli.N_mphi < 1 || cli.tan_beta <= 0.0) {
        print_usage(argv[0]);
        return 1;
    }
    if (cli.N_m12 < 1) { std::cerr << "N_m12 debe ser positivo.\n"; return 1; }
    if (cli.y_type < 1 || cli.y_type > 4) { std::cerr << "yukawa_type debe estar entre 1 y 4.\n"; return 1; }

    // 4) Configurar/Reportar hilos/cores
    unsigned hw = std::thread::hardware_concurrency();
#ifdef _OPENMP
    if (cli.threads > 0) omp_set_num_threads(cli.threads);
    int used = omp_get_max_threads();
    cout << "[INFO] OpenMP habilitado. Hilos: " << used
         << " (hardware_concurrency=" << (hw? hw : 0) << ")\n";
#else
    if (cli.threads > 0) {
        cout << "[WARN] OpenMP no está habilitado en compilación; --threads será ignorado.\n";
    }
    cout << "[INFO] OpenMP deshabilitado. Ejecutando en 1 hilo. "
         << "Cores detectados: " << (hw? hw : 0) << "\n";
#endif

    // 5) Parámetros fijos (permitimos override por flags)
    FixedParams fp{ cli.mA, cli.sin_ba, cli.tan_beta, cli.lambda6, cli.lambda7 };

    const double TOL = std::pow(10.0, -cli.TOL_exp);

    // 6) Ejecutar
    perform_param_scan(
        cli.y_type,
        cli.mphi_min, cli.mphi_max,
        cli.N_mphi, cli.N_m12, TOL,
        fp.mA, fp.sin_ba, fp.tan_beta,
        fp.lambda6, fp.lambda7,
        cli.output_csv
    );

    return 0;
}
