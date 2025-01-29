
/*
Este codigo de testing es ideal para probar el archivo de configuración y además hacer pruebas si existen probelmas de escreibir datos
permitiendo mayor cantidad de debugging
*/
#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"
#include "ParamUtils.hpp"  // Header común

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main()
{
    try
    {
        string output_file = "testing.csv";
        ofstream results(output_file);
        if (!results.is_open()) {
            cerr << "Failed to open output file: " << output_file << endl;
        }

        // Columnas CSV
        vector<string> columns = {
            "lambda1", "lambda2", "lambda3", "lambda4", "lambda5",
            "m12_squared", "beta", "tan_beta",
            "width_hgaga", "total_width_h2", "branching_ratio_hgaga"
        };
        write_csv_header(results, columns);


        // 1. Leer archivo de configuración
        Config cfg = readConfig("test_parametros.conf");

        // 2. Calcular total de iteraciones
        double total_iterations = computeTotalIterations(cfg);
        cout << "Total de iteraciones: " << total_iterations << endl;

        double current_iteration = 0.0;
        auto start_time = high_resolution_clock::now();
        // 3. Recorrer los bucles usando los valores del cfg
        for(double lambda1 = cfg.lambda1_min; lambda1 <= cfg.lambda1_max; lambda1 += cfg.step_lambda1)
        {
            for(double lambda2 = cfg.lambda2_min; lambda2 <= cfg.lambda2_max; lambda2 += cfg.step_lambda2)
            {
                for(double lambda3 = cfg.lambda3_min; lambda3 <= cfg.lambda3_max; lambda3 += cfg.step_lambda3)
                {
                    for(double lambda4 = cfg.lambda4_min; lambda4 <= cfg.lambda4_max; lambda4 += cfg.step_lambda4)
                    {
                        for(double lambda5 = cfg.lambda5_min; lambda5 <= cfg.lambda5_max; lambda5 += cfg.step_lambda5)
                        {
                            for(double m12 = cfg.m12_squared_min; 
                                        m12 <= cfg.m12_squared_max; 
                                        m12 += cfg.step_m12_squared)
                            {
                                for(double beta = cfg.beta_min; beta <= cfg.beta_max; beta += cfg.step_beta)
                                {
                                    // Aquí hacemos lo que necesitemos con los parámetros (lambda1, lambda2, etc.)
                                    // ...
                                    current_iteration += 1.0;
                                    auto now_time = high_resolution_clock::now();
                                    double elapsed = duration<double>(now_time - start_time).count();
                                    double progress = (current_iteration / total_iterations);
                                    print_progress(progress, elapsed, total_iterations, current_iteration);

                                    // Chequeo de positividad
                                    if (!check_positivity(lambda1, lambda2, lambda3, lambda4, lambda5)) {
                                        continue;
                                    }

                                    // Inicializar modelo 2HDMC
                                    THDM model;
                                    SM sm;
                                    model.set_SM(sm);

                                    bool pset = model.set_param_gen(
                                        lambda1, lambda2, lambda3, lambda4, lambda5,
                                        0.0, 0.0,  // lambda6, lambda7 = 0
                                        m12,       // m12^2
                                        tan(beta)  // tan_beta
                                    );
                                    if (!pset) {
                                        continue;
                                    }

                                    Constraints check(model);
                                    //if (!check.check_positivity() || !check.check_unitarity() || !check.check_perturbativity()) {
                                    if (!check.check_positivity() || !check.check_unitarity()) {
                                            std::cout << "DEBUG: Falla check_positivity con: "
                                            << lambda1 << "," << lambda2 << "," << lambda3 << "," 
                                            << lambda4 << "," << lambda5 << std::endl;
                                        continue;
                                    }

                                    DecayTable table(model);
                                    double width_hgaga = table.get_gamma_hgaga(2); // H2 -> gamma gamma
                                    double total_w_h2  = table.get_gammatot_h(2);

                                    // validación para posibles nan
                                    // if (std::isnan(width_hgaga) || std::isnan(total_w_h2) || total_w_h2 <= 0.0) {
                                    //     continue;
                                    // }

                                    double br_hgaga = width_hgaga / total_w_h2;

                                    // Guardar resultados
                                    vector<double> row = {
                                        lambda1, lambda2, lambda3, lambda4, lambda5,
                                        m12, beta, tan(beta),
                                        width_hgaga, total_w_h2, br_hgaga
                                    };
                                    write_csv_row(results, row);



                                    // Actualizamos el contador
                                    current_iteration++;
                                    // (Opcional) imprimir avance cada cierto número de iteraciones
                                    if (fmod(current_iteration, 100000) == 0) {
                                        cout << "Progreso: " << (current_iteration / total_iterations) * 100.0
                                             << "%\r" << flush;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        cout << "\nFinalizado. Iteraciones totales procesadas: " << current_iteration << endl;
    }
    catch(const exception &e)
    {
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}
