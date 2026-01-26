#pragma once

// Activar/desactivar benchmarking I/O vs Cálculos
#ifndef BENCHMARK_IO
#define BENCHMARK_IO 1
#endif

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <atomic> // Necesario para ScanMonitor
#include <regex>
#include <sstream>

// Si quieres usar directamente un namespace o no, es decisión de estilo.
// using namespace std;

// Estructuras de configuración básicas
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

struct ConfigPhys {
    double lambda6_min;
    double lambda6_max;
    double step_lambda6;

    double lambda7_min;
    double lambda7_max;
    double step_lambda7;

    double m12_squared_min;
    double m12_squared_max;
    double step_m12_squared;

    double alpha_min;
    double alpha_max;
    double step_alpha;

    double beta_min;
    double beta_max;
    double step_beta;

    double mphi_min;
    double mphi_max;
    double step_mphi;

    double mA_min;
    double mA_max;
    double step_mA;
};

// Estructuras para JSON parsing
struct QuadraticModel {
    double a, b, c;
};

struct SearchSettings {
    int N_mphi, N_m12;
    double m12_min, m12_max;
    double mphi_min, mphi_max;
};

struct FixedParameters {
    double mA, sin_ba, tan_beta, lambda6, lambda7;
};

#include <mutex> // Para proteger acumuladores de tiempo

// =============================================================
//  Clase ScanMonitor (HPC Real-Time Monitoring)
// =============================================================
class ScanMonitor {
public:
    // Constructor: recibe el total de iteraciones estimadas para calcular %
    ScanMonitor(long long total_tasks);

    // Registra intentos (cálculos realizados). Thread-safe (lock-free).
    // Llamar frecuentemente (o en batches locales).
    void record_attempts(long long n);

    // Registra puntos válidos y verifica si es hora de imprimir en consola.
    // IMPORTANTE: Esta función debe llamarse dentro de una sección crítica o manejarse con cuidado
    // si se llama concurrentemente, ya que imprime en std::cerr.
    void update_valid_and_report(long long n_valid);

    // Registra tiempo gastado en I/O (en segundos). Thread-safe.
    void record_io_time(double seconds);

    // Registra tiempo gastado en cálculos (en segundos). Thread-safe.
    void record_calc_time(double seconds);

    // Imprime el resumen final
    void finish();

private:
    long long total_tasks;
    std::atomic<long long> global_attempts;
    std::atomic<long long> global_valid;
    
    // Acumuladores de tiempo para benchmarking I/O vs Cálculos
    double total_io_time;
    double total_calc_time;
    std::mutex time_mutex; // Mutex para proteger los acumuladores de tiempo
    
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::chrono::time_point<std::chrono::high_resolution_clock> last_report_time;
    
    // Configuración interna
    const double REPORT_INTERVAL_SEC = 0.5; 
};

// =============================================================
//  Funciones Utilitarias (I/O, Config, Física)
// =============================================================

// Función para leer configuración desde un archivo (opcional)
Config readConfig(const std::string &filename);

// Lectura de config "física"
ConfigPhys readPhysConfig(const std::string &filename);

// Lanza std::runtime_error en caso de fallo
void parse_json_config(
    const std::string &filename,
    QuadraticModel    &qm,
    SearchSettings    &ss,
    FixedParameters   &fp
);

// Cálculo de número de pasos
long stepsCount(double minVal, double maxVal, double step);

// Cálculo de iteraciones totales (opcional, si tienes algo genérico)
double computeTotalIterations(const Config &cfg);

// (Opcional) Cálculo de iteraciones totales para config física
double computeTotalIterationsPhys(const ConfigPhys &cfg);

// Función para chequear positividad
bool check_positivity(double lambda1, double lambda2,
                      double lambda3, double lambda4, double lambda5);

// Funciones para CSV
void write_csv_header(std::ofstream &results, const std::vector<std::string> &columns);
void write_csv_row(std::ofstream &results, const std::vector<double> &values);
void write_csv_row(std::ofstream &results, const std::vector<std::string> &values);

// Función para imprimir progreso (Legacy)
void print_progress(double progress, double elapsed_time,
                    double total_iterations, double current_iteration);

// Cálculos físicos auxiliares
double calc_lambda1(double mh, double mH,
                    double m12_2,
                    double sa, double ca, double sb, double cb, double tan_beta,
                    double lam6, double lam7, double VEV);

double calc_lambda2(double mh, double mH,
                    double m12_2,
                    double sa, double ca, double sb, double cb, double tan_beta,
                    double lam6, double lam7, double VEV);

bool check_lambda(double l1);

double m12_base(double mh, double m_phi, double sa, double ca, 
    double sb, double cb, 
    double lambda6, double lambda7,
    double tan_beta, double VEV);