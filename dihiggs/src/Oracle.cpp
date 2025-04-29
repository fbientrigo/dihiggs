// Oracle.cpp
#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"  // CSV & config utilities (si se necesita)
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace std::chrono;

// Función que calcula el punto y retorna un string JSON con los resultados.
string calculatePoint(double m_phi, double mA, double sin_ba, double tan_beta,
                        double lambda6, double lambda7, double m12) {
    // Parámetros fijos y derivados
    double m_h = 125.0;         // Higgs ligero (SM) fijo
    double m_Hp = mA;           // Se asume m_Hp = mA
    //double sin_ba = sin(beta - alpha);
    //double tan_beta = tan(beta);

    // Construir el modelo 2HDM
    THDM model;
    SM sm;
    model.set_SM(sm);

    // Intentar fijar los parámetros físicos
    bool pset = model.set_param_phys(m_h, m_phi, mA, m_Hp, sin_ba, lambda6, lambda7, m12, tan_beta);
    if (!pset) {
        return "{\"error\": \"Invalid parameter set.\"}";
    }

    // Informacion de los lambdas
    double lam1, lam2, lam3, lam4, lam5, lam6_g, lam7_g, m12_2_g, tanb_g;
    model.get_param_gen( lam1, lam2, lam3, lam4, lam5, lam6_g, lam7_g, m12_2_g, tanb_g);

    // Chequeos de estabilidad, unitariedad y perturbatividad
    Constraints check(model);
    bool positivity_ok     = check.check_positivity();
    bool unitarity_ok      = check.check_unitarity();
    bool perturbativity_ok = check.check_perturbativity();

    // Calcular los anchos usando DecayTable
    DecayTable table(model);

    // Configurar precisión y notación científica en la salida
    stringstream ss;
    ss << scientific << setprecision(15);

    // Canales fermiónicos:
    double w_h2_bb     = table.get_gamma_hdd(2, 3, 3);
    double w_h2_tautau = table.get_gamma_hll(2, 3, 3);
    // Para w_h2_uu, w_h2_du y w_h2_ln se incluyen ejemplos solo del elemento (1,1)
    double w_h2_uu     = table.get_gamma_huu(2, 1, 1);
    double w_h2_du     = table.get_gamma_hdu(2, 1, 1);
    double w_h2_ln     = table.get_gamma_hln(2, 1, 1);

    // Canales vectoriales:
    double w_h2_vv1 = table.get_gamma_hvv(2, 1);
    double w_h2_vv2 = table.get_gamma_hvv(2, 2);
    double w_h2_vv3 = table.get_gamma_hvv(2, 3);

    // Decaimientos especiales:
    double w_h2_gaga   = table.get_gamma_hgaga(2);
    double w_h2_Zga    = table.get_gamma_hZga(2);
    double w_h2_gg     = table.get_gamma_hgg(2);

    // Decaimiento Higgs a Higgs (ejemplo del canal 1,1)
    double w_h2_hh     = table.get_gamma_hhh(2, 1, 1);

    // Totales:
    double w_total_h2  = table.get_gammatot_h(2);
    double w_total_top = table.get_gammatot_top();
    double br_h2_gaga  = w_h2_gaga / w_total_h2;

    // Construir salida JSON para este punto
    ss << "{\n";
    ss << "  \"positivity_ok\": " << (positivity_ok ) << ",\n";
    ss << "  \"unitarity_ok\": " << (unitarity_ok ) << ",\n";
    ss << "  \"perturbativity_ok\": " << (perturbativity_ok ) << ",\n";
    ss << "  \"w_h2_bb\": " << w_h2_bb << ",\n";
    ss << "  \"w_h2_tautau\": " << w_h2_tautau << ",\n";
    ss << "  \"w_h2_uu\": " << w_h2_uu << ",\n";
    ss << "  \"w_h2_du\": " << w_h2_du << ",\n";
    ss << "  \"w_h2_ln\": " << w_h2_ln << ",\n";
    ss << "  \"w_h2_vv\": [" << w_h2_vv1 << ", " << w_h2_vv2 << ", " << w_h2_vv3 << "],\n";
    ss << "  \"w_h2_gaga\": " << w_h2_gaga << ",\n";
    ss << "  \"w_h2_Zga\": " << w_h2_Zga << ",\n";
    ss << "  \"w_h2_gg\": " << w_h2_gg << ",\n";
    ss << "  \"w_h2_hh\": " << w_h2_hh << ",\n";
    ss << "  \"w_total_h2\": " << w_total_h2 << ",\n";
    ss << "  \"w_total_top\": " << w_total_top << ",\n";
    ss << "  \"branching_ratio_h2_gaga\": " << br_h2_gaga << "\n";
    // ────── Volcar λ’s de la base genérica ──────
    ss << ",\n  \"lambda1\": " << lam1
       << ",\n  \"lambda2\": " << lam2
       << ",\n  \"lambda3\": " << lam3
       << ",\n  \"lambda4\": " << lam4
       << ",\n  \"lambda5\": " << lam5
       << ",\n  \"lambda6\": " << lam6_g
       << ",\n  \"lambda7\": " << lam7_g
       << ",\n  \"m12_2\": "   << m12_2_g;
    // ────────────────────────────────────────────
    ss << "}";
    return ss.str();
}

int main(int argc, char* argv[]) {
    // Modo serial: si se reciben exactamente 7 parámetros.
    if (argc == 8) {
        string result = calculatePoint(atof(argv[1]), atof(argv[2]), atof(argv[3]),
                                       atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));
        cout << result;
        return 0;
    }
    
    // Modo paralelo: Se espera la bandera "--nthreads" seguida del número de hilos
    // y luego grupos de 7 parámetros.
    int startIndex = 1;
    int nthreads = 1;
    if (string(argv[1]) == "--nthreads") {
        nthreads = atoi(argv[2]);
        startIndex = 3;
    }
    
    int nParams = argc - startIndex;
    if (nParams < 7) {
        cerr << "Error: No se encontraron suficientes parámetros para formar un grupo de 7.\n";
        return 1;
    }
    
    int nPoints = nParams / 7;
    if (nParams % 7 != 0) {
        cerr << "Error: El número de parámetros debe ser múltiplo de 7.\n";
        return 1;
    }
    
    vector<string> outputs(nPoints);
    
    #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
    for (int i = 0; i < nPoints; i++) {
        double p0 = atof(argv[startIndex + i*7 + 0]);
        double p1 = atof(argv[startIndex + i*7 + 1]);
        double p2 = atof(argv[startIndex + i*7 + 2]);
        double p3 = atof(argv[startIndex + i*7 + 3]);
        double p4 = atof(argv[startIndex + i*7 + 4]);
        double p5 = atof(argv[startIndex + i*7 + 5]);
        double p6 = atof(argv[startIndex + i*7 + 6]);
        string res = calculatePoint(p0, p1, p2, p3, p4, p5, p6);
        outputs[i] = res;
    }
    
    // Construir un JSON que contenga todos los resultados en un array.
    stringstream ss;
    ss << "[\n";
    for (int i = 0; i < nPoints; i++) {
        ss << outputs[i];
        if (i < nPoints - 1)
            ss << ",\n";
    }
    ss << "\n]\n";
    cout << ss.str();
    return 0;
}
