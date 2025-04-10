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
#include <vector>
#include <memory>
#include <sstream>

using namespace std;
using namespace std::chrono;

// -----------------------------------------------------------------------------
// Estructura para almacenar el conjunto de parámetros
// -----------------------------------------------------------------------------
struct ParamSet {
    double m_phi;
    double mA;       // m_A = m_Hp
    double alpha;
    double beta;
    double lambda6;
    double lambda7;
    double m12;
};

// -----------------------------------------------------------------------------
// Estructura para almacenar los resultados de la evaluación
// -----------------------------------------------------------------------------
struct EvalResult {
    ParamSet params;
    double sin_ba;
    double tan_beta;
    bool positivity_ok;
    bool unitarity_ok;
    bool perturbativity_ok;
    double w_h2_bb;
    double w_h2_tautau;
    double w_h2_WW;
    double w_h2_ZZ;
    double w_h2_gaga;
    double w_h2_Zga;
    double w_h2_gg;
    double w_h2_hh;
    double w_total_h2;
    double br_h2_gaga;
};

// -----------------------------------------------------------------------------
// Clase abstracta para estrategias de muestreo (Sampling)
// -----------------------------------------------------------------------------
class Sampler {
public:
    virtual vector<ParamSet> generateSamples(const ConfigPhys &cfg) = 0;
    virtual ~Sampler() {}
};

// -----------------------------------------------------------------------------
// Implementación de un muestreo en rejilla (Grid Sampling)
// -----------------------------------------------------------------------------
class GridSampler : public Sampler {
public:
    vector<ParamSet> generateSamples(const ConfigPhys &cfg) override {
        vector<ParamSet> samples;
        // Calcular número de pasos para cada parámetro (con chequeo para step=0)
        int steps_lambda6 = (cfg.lambda6_max - cfg.lambda6_min) / cfg.step_lambda6 + 1;
        int steps_lambda7 = (cfg.lambda7_max - cfg.lambda7_min) / cfg.step_lambda7 + 1;
        int steps_m12     = (cfg.m12_squared_max - cfg.m12_squared_min) / cfg.step_m12_squared + 1;
        int steps_alpha   = (cfg.alpha_max - cfg.alpha_min) / cfg.step_alpha + 1;
        int steps_beta    = (cfg.beta_max - cfg.beta_min) / cfg.step_beta + 1;
        int steps_mphi    = (cfg.mphi_max - cfg.mphi_min) / cfg.step_mphi + 1;
        int steps_mA      = (cfg.mA_max - cfg.mA_min) / cfg.step_mA + 1;

        if (cfg.step_lambda6 == 0.0) steps_lambda6 = 1;
        if (cfg.step_lambda7 == 0.0) steps_lambda7 = 1;
        if (cfg.step_m12_squared == 0.0) steps_m12 = 1;
        if (cfg.step_alpha == 0.0) steps_alpha = 1;
        if (cfg.step_beta == 0.0) steps_beta = 1;
        if (cfg.step_mphi == 0.0) steps_mphi = 1;
        if (cfg.step_mA == 0.0) steps_mA = 1;

        // Generación del grid mediante bucles anidados
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
                            for (int i_mphi = 0; i_mphi < steps_mphi; i_mphi++) {
                                double m_phi = cfg.mphi_min + i_mphi * cfg.step_mphi;
                                for (int i_mA = 0; i_mA < steps_mA; i_mA++) {
                                    double mA = cfg.mA_min + i_mA * cfg.step_mA;
                                    ParamSet ps = { m_phi, mA, alpha, beta, lambda6, lambda7, m12 };
                                    samples.push_back(ps);
                                }
                            }
                        }
                    }
                }
            }
        }
        return samples;
    }
};

// -----------------------------------------------------------------------------
// Función para evaluar un conjunto de parámetros (cada punto de la búsqueda)
// -----------------------------------------------------------------------------
EvalResult evaluateParamSet(const ParamSet &ps) {
    EvalResult res;
    res.params = ps;
    res.sin_ba = sin(ps.beta - ps.alpha);
    res.tan_beta = tan(ps.beta);

    // Construir el modelo 2HDM y fijar parámetros físicos
    THDM model;
    SM sm;
    model.set_SM(sm);
    double m_h = 125.0;     // Fijamos m_h
    double m_Hp = ps.mA;    // Condición: m_A = m_Hp
    // pset se usaba para comprobar si set_param_phys se realizó correctamente
    // Lo dejamos comentado si no lo necesitamos en este momento
    bool success = model.set_param_phys(m_h, ps.m_phi, ps.mA, m_Hp, res.sin_ba, ps.lambda6, ps.lambda7, ps.m12, res.tan_beta);
    // (void)success; // Para silenciar 'unused variable' si no lo usas

    // Chequear condiciones físicas
    Constraints check(model);
    res.positivity_ok     = check.check_positivity();
    res.unitarity_ok      = check.check_unitarity();
    res.perturbativity_ok = check.check_perturbativity();

    // Calcular anchos de decaimiento
    DecayTable table(model);
    res.w_h2_bb     = table.get_gamma_hdd(2, 3, 3);
    res.w_h2_tautau = table.get_gamma_hll(2, 3, 3);
    res.w_h2_WW     = table.get_gamma_hvv(2, 3);
    res.w_h2_ZZ     = table.get_gamma_hvv(2, 2);
    res.w_h2_gaga   = table.get_gamma_hgaga(2);
    res.w_h2_Zga    = table.get_gamma_hZga(2);
    res.w_h2_gg     = table.get_gamma_hgg(2);
    res.w_h2_hh     = table.get_gamma_hhh(2, 1, 1);
    res.w_total_h2  = table.get_gammatot_h(2);
    res.br_h2_gaga  = (res.w_total_h2 > 1e-12) ? res.w_h2_gaga / res.w_total_h2 : 0.0;

    return res;
}

// -----------------------------------------------------------------------------
// Función para escribir una fila de resultados en CSV
// -----------------------------------------------------------------------------
void writeCSVRow(ofstream &ofs, const EvalResult &res) {
    ofs << res.params.m_phi << ","
        << res.params.mA << ","
        << res.params.alpha << ","
        << res.params.beta << ","
        << res.params.lambda6 << ","
        << res.params.lambda7 << ","
        << res.params.m12 << ","
        << res.sin_ba << ","
        << res.tan_beta << ","
        << (res.positivity_ok ? 1 : 0) << ","
        << (res.unitarity_ok ? 1 : 0) << ","
        << (res.perturbativity_ok ? 1 : 0) << ","
        << res.w_h2_bb << ","
        << res.w_h2_tautau << ","
        << res.w_h2_WW << ","
        << res.w_h2_ZZ << ","
        << res.w_h2_gaga << ","
        << res.w_h2_Zga << ","
        << res.w_h2_gg << ","
        << res.w_h2_hh << ","
        << res.w_total_h2 << ","
        << res.br_h2_gaga << "\n";
}

// -----------------------------------------------------------------------------
// Función principal de escaneo: usa la estrategia de muestreo pasada para generar
// los puntos y luego evalúa cada uno en paralelo (OpenMP)
// -----------------------------------------------------------------------------
void perform_param_scan(const ConfigPhys &cfg, const string &output_file, Sampler &sampler) {
    // Abrir archivo CSV para resultados
    ofstream results(output_file.c_str());
    if (!results.is_open()) {
        cerr << "Failed to open output file: " << output_file << endl;
        return;
    }
    // Escribir cabecera del CSV
    results << "m_phi,mA,alpha,beta,lambda6,lambda7,m12,sin_ba,tan_beta,positivity_ok,unitarity_ok,perturbativity_ok,"
            << "w_h2_bb,w_h2_tautau,w_h2_WW,w_h2_ZZ,w_h2_gaga,w_h2_Zga,w_h2_gg,w_h2_hh,w_total_h2,br_h2_gaga\n";

    // Generar muestras usando la estrategia de muestreo definida
    vector<ParamSet> samples = sampler.generateSamples(cfg);
    cout << "Total samples generated: " << samples.size() << "\n";

    // Variables para seguimiento del mejor resultado global
    double global_bestBR = -1.0;
    ParamSet global_bestParams;
    auto start_time = high_resolution_clock::now();

    // Evaluar cada muestra en paralelo usando OpenMP
    #pragma omp parallel
    {
        vector<string> thread_local_results;
        double local_bestBR = -1.0;
        ParamSet local_bestParams;

        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < samples.size(); i++) {
            EvalResult res = evaluateParamSet(samples[i]);
            // Preparar la fila CSV para esta muestra
            ostringstream oss;
            oss << fixed << setprecision(8)
                << res.params.m_phi << ","
                << res.params.mA << ","
                << res.params.alpha << ","
                << res.params.beta << ","
                << res.params.lambda6 << ","
                << res.params.lambda7 << ","
                << res.params.m12 << ","
                << res.sin_ba << ","
                << res.tan_beta << ","
                << (res.positivity_ok ? 1 : 0) << ","
                << (res.unitarity_ok ? 1 : 0) << ","
                << (res.perturbativity_ok ? 1 : 0) << ","
                << res.w_h2_bb << ","
                << res.w_h2_tautau << ","
                << res.w_h2_WW << ","
                << res.w_h2_ZZ << ","
                << res.w_h2_gaga << ","
                << res.w_h2_Zga << ","
                << res.w_h2_gg << ","
                << res.w_h2_hh << ","
                << res.w_total_h2 << ","
                << res.br_h2_gaga;

            thread_local_results.push_back(oss.str());

            // Actualizar el mejor resultado local (si cumple condiciones)
            if (res.positivity_ok && res.w_total_h2 > 0.0 && res.br_h2_gaga > local_bestBR) {
                local_bestBR = res.br_h2_gaga;
                local_bestParams = samples[i];
            }
        }

        // Sección crítica: escribir resultados del hilo al archivo global y actualizar el mejor global
        #pragma omp critical
        {
            for (size_t j = 0; j < thread_local_results.size(); j++) {
                results << thread_local_results[j] << "\n";
            }
            if (local_bestBR > global_bestBR) {
                global_bestBR = local_bestBR;
                global_bestParams = local_bestParams;
            }
        }
    } // Fin del paralelo

    results.close();
    auto elapsed = duration<double>(high_resolution_clock::now() - start_time).count();
    cout << "Scan completed in " << elapsed << " seconds." << endl;
    cout << "Best BR(h->gaga): " << global_bestBR << " found at parameters:" << endl;
    cout << " m_phi: " << global_bestParams.m_phi 
         << ", mA: " << global_bestParams.mA 
         << ", alpha: " << global_bestParams.alpha 
         << ", beta: " << global_bestParams.beta 
         << ", lambda6: " << global_bestParams.lambda6 
         << ", lambda7: " << global_bestParams.lambda7 
         << ", m12: " << global_bestParams.m12 << endl;
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

        // Crear un Sampler con new en vez de make_unique (C++11 compatible)
        unique_ptr<Sampler> sampler(new GridSampler());

        // Ejecutar el escaneo con la estrategia de muestreo seleccionada
        perform_param_scan(cfg, output_file, *sampler);

    } catch (const exception &e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
