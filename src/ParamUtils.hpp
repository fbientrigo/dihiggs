#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <chrono>

// Si quieres usar directamente un namespace o no, es decisión de estilo.
// using namespace std;

struct Config {
    // Rangos
    double lambda1_min;
    double lambda1_max;
    double step_lambda1;

    double lambda2_min;
    double lambda2_max;
    double step_lambda2;

    double lambda3_min;
    double lambda3_max;
    double step_lambda3;

    double lambda4_min;
    double lambda4_max;
    double step_lambda4;

    double lambda5_min;
    double lambda5_max;
    double step_lambda5;

    double m12_squared_min;
    double m12_squared_max;
    double step_m12_squared;

    double beta_min;
    double beta_max;
    double step_beta;

};

// Función para leer configuración desde un archivo (opcional)
Config readConfig(const std::string &filename);

// Cálculo de número de pasos
long stepsCount(double minVal, double maxVal, double step);

// Cálculo de iteraciones totales (opcional, si tienes algo genérico)
double computeTotalIterations(const Config &cfg);

// Función para chequear positividad
bool check_positivity(double lambda1, double lambda2,
                      double lambda3, double lambda4, double lambda5);

// Funciones para CSV
void write_csv_header(std::ofstream &results, const std::vector<std::string> &columns);
void write_csv_row(std::ofstream &results, const std::vector<double> &values);

// Función para imprimir progreso
void print_progress(double progress, double elapsed_time,
                    double total_iterations, double current_iteration);
