#include "ParamUtils.hpp"
#include <iomanip> // si usas std::setw, std::setprecision, etc.
#include <sstream>
#include <unordered_map> // si lees key-value

using namespace std;
using namespace std::chrono;

// Implementación de readConfig (opcional)
Config readConfig(const string &filename)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        throw runtime_error("No se pudo abrir el archivo de configuración: " + filename);
    }
    
    unordered_map<string, double> configMap;
    string key;
    double value;
    //int count=0;
    //cout << "Leyendo archivo config\n";
    // Leemos línea a línea: "key value"
    while (file >> key >> value) {
        configMap[key] = value;
    //    count++;
    }
    //cout << "Leidas " << count << " lineas del archivo de config.\n";

    // Cerramos el archivo
    file.close();

    // Verificamos que existen todas las claves que necesitamos.
    // En caso contrario, se podría usar algún valor por defecto o lanzar excepción.
    // Por simplicidad, aquí asumimos que todas están presentes.

    Config cfg;
    cfg.lambda1_min      = configMap.at("lambda1_min");
    cfg.lambda1_max      = configMap.at("lambda1_max");
    cfg.step_lambda1     = configMap.at("step_lambda1");

    cfg.lambda2_min      = configMap.at("lambda2_min");
    cfg.lambda2_max      = configMap.at("lambda2_max");
    cfg.step_lambda2     = configMap.at("step_lambda2");

    cfg.lambda3_min      = configMap.at("lambda3_min");
    cfg.lambda3_max      = configMap.at("lambda3_max");
    cfg.step_lambda3     = configMap.at("step_lambda3");

    cfg.lambda4_min      = configMap.at("lambda4_min");
    cfg.lambda4_max      = configMap.at("lambda4_max");
    cfg.step_lambda4     = configMap.at("step_lambda4");

    cfg.lambda5_min      = configMap.at("lambda5_min");
    cfg.lambda5_max      = configMap.at("lambda5_max");
    cfg.step_lambda5     = configMap.at("step_lambda5");

    cfg.m12_squared_min  = configMap.at("m12_squared_min");
    cfg.m12_squared_max  = configMap.at("m12_squared_max");
    cfg.step_m12_squared = configMap.at("step_m12_squared");

    cfg.beta_min         = configMap.at("beta_min");
    cfg.beta_max         = configMap.at("beta_max");
    cfg.step_beta        = configMap.at("beta_step");

    return cfg;
}


// ---------------------------------------------------------------------
//  Lectura de parámetros "físicos": readPhysConfig
// ---------------------------------------------------------------------
ConfigPhys readPhysConfig(const std::string &filename)
{
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("No se pudo abrir archivo de configuración física: " + filename);
    }

    unordered_map<string, double> configMap;
    string key;
    double value;
    while (file >> key >> value) {
        configMap[key] = value;
    }
    file.close();

    ConfigPhys cfgp;
    cfgp.lambda6_min      = configMap.at("lambda6_min");
    cfgp.lambda6_max      = configMap.at("lambda6_max");
    cfgp.step_lambda6     = configMap.at("step_lambda6");

    cfgp.lambda7_min      = configMap.at("lambda7_min");
    cfgp.lambda7_max      = configMap.at("lambda7_max");
    cfgp.step_lambda7     = configMap.at("step_lambda7");

    cfgp.m12_squared_min  = configMap.at("m12_squared_min");
    cfgp.m12_squared_max  = configMap.at("m12_squared_max");
    cfgp.step_m12_squared = configMap.at("step_m12_squared");

    cfgp.alpha_min        = configMap.at("alpha_min");
    cfgp.alpha_max        = configMap.at("alpha_max");
    cfgp.step_alpha       = configMap.at("step_alpha");

    cfgp.beta_min         = configMap.at("beta_min");
    cfgp.beta_max         = configMap.at("beta_max");
    cfgp.step_beta        = configMap.at("step_beta");

    cfgp.mphi_min         = configMap.at("mphi_min");
    cfgp.mphi_max         = configMap.at("mphi_max");
    cfgp.step_mphi        = configMap.at("step_mphi");

    cfgp.mA_min           = configMap.at("mA_min");
    cfgp.mA_max           = configMap.at("mA_max");
    cfgp.step_mA          = configMap.at("step_mA");

    return cfgp;
}



// Función auxiliar para calcular el número de pasos (evita repetición de código)
inline long stepsCount(double minVal, double maxVal, double step)
{
    // floor(...) + 1 asumiendo que (maxVal - minVal) es múltiplo de step o
    // que no necesitamos truncar con exactitud. Puedes ajustar la lógica según convenga.
    return static_cast<long>(floor((maxVal - minVal) / step)) + 1;
}


double computeTotalIterations(const Config &cfg)
{
    // Calculamos iteraciones en cada dimensión
    long it_lambda1 = stepsCount(cfg.lambda1_min, cfg.lambda1_max, cfg.step_lambda1);
    long it_lambda2 = stepsCount(cfg.lambda2_min, cfg.lambda2_max, cfg.step_lambda2);
    long it_lambda3 = stepsCount(cfg.lambda3_min, cfg.lambda3_max, cfg.step_lambda3);
    long it_lambda4 = stepsCount(cfg.lambda4_min, cfg.lambda4_max, cfg.step_lambda4);
    long it_lambda5 = stepsCount(cfg.lambda5_min, cfg.lambda5_max, cfg.step_lambda5);
    
    long it_m12     = stepsCount(cfg.m12_squared_min, cfg.m12_squared_max, cfg.step_m12_squared);
    long it_beta    = stepsCount(cfg.beta_min,        cfg.beta_max,        cfg.step_beta);

    // Multiplicamos para obtener total
    double total = static_cast<double>(it_lambda1) *
                   static_cast<double>(it_lambda2) *
                   static_cast<double>(it_lambda3) *
                   static_cast<double>(it_lambda4) *
                   static_cast<double>(it_lambda5) *
                   static_cast<double>(it_m12)     *
                   static_cast<double>(it_beta);
    cout << "Total iterations: " << total << "\n";
    return total;
}

// ---------------------------------------------------------------------
//  (Opcional) Cálculo de iteraciones totales (config física)
// ---------------------------------------------------------------------
double computeTotalIterationsPhys(const ConfigPhys &cfgp)
{
    long it_l6     = stepsCount(cfgp.lambda6_min,      cfgp.lambda6_max,      cfgp.step_lambda6);
    long it_l7     = stepsCount(cfgp.lambda7_min,      cfgp.lambda7_max,      cfgp.step_lambda7);
    long it_m12    = stepsCount(cfgp.m12_squared_min,  cfgp.m12_squared_max,  cfgp.step_m12_squared);
    long it_alpha  = stepsCount(cfgp.alpha_min,        cfgp.alpha_max,        cfgp.step_alpha);
    long it_beta   = stepsCount(cfgp.beta_min,         cfgp.beta_max,         cfgp.step_beta);
    long it_mphi   = stepsCount(cfgp.mphi_min,         cfgp.mphi_max,         cfgp.step_mphi);
    long it_mA     = stepsCount(cfgp.mA_min,           cfgp.mA_max,           cfgp.step_mA);

    double total = static_cast<double>(it_l6)
                 * static_cast<double>(it_l7)
                 * static_cast<double>(it_m12)
                 * static_cast<double>(it_alpha)
                 * static_cast<double>(it_beta)
                 * static_cast<double>(it_mphi)
                 * static_cast<double>(it_mA);
    cout << "Total iterations: " << total << "\n";
    return total;
}

// Chequeo de positividad
bool check_positivity(double lambda1, double lambda2,
                      double lambda3, double lambda4, double lambda5) {
    if (lambda1 < 0 || lambda2 < 0) return false;
    if (lambda3 < -sqrt(lambda1 * lambda2)) return false;
    if (lambda3 + lambda4 - fabs(lambda5) < -sqrt(lambda1 * lambda2)) return false;
    return true;
}

// Escribir cabecera CSV
void write_csv_header(std::ofstream &results, const std::vector<std::string> &columns) {
    if (!results.is_open()) return;

    for (size_t i = 0; i < columns.size(); ++i) {
        results << columns[i];
        if (i < columns.size() - 1) results << ",";
    }
    results << endl;
}

// Escribir una fila al CSV
// Versión para vectores de doubles
void write_csv_row(std::ofstream &results, const std::vector<double> &values) {
    if (!results.is_open()) return;
    results << std::fixed << std::setprecision(15);  // or whichever precision you prefer
    for (size_t i = 0; i < values.size(); ++i) {
        results << values[i];
        if (i < values.size() - 1) results << ",";
    }
    results << std::endl;
}

// Versión para vectores de strings
void write_csv_row(std::ofstream &results, const std::vector<std::string> &values) {
    if (!results.is_open()) return;
    

    for (size_t i = 0; i < values.size(); ++i) {
        results << values[i];
        if (i < values.size() - 1) results << ",";
    }
    results << std::endl;
}


// Imprimir progreso
void print_progress(double progress, double elapsed_time,
                    double total_iterations, double current_iteration) {
    int bar_width = 50;
    cout << "[";
    int pos = int(bar_width * progress);
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) cout << "=";
        else if (i == pos) cout << ">";
        else cout << " ";
    }
    double time_per_iteration = elapsed_time / current_iteration;
    double remaining_time = time_per_iteration * (total_iterations - current_iteration);
    cout << "] " << int(progress * 100.0) << "% | Elapsed: "
         << int(elapsed_time) << "s | Remaining: " << int(remaining_time) << "s\r";
    cout.flush();
}
