#include "ParamUtils.hpp"
#include <iomanip> // si usas std::setw, std::setprecision, etc.
#include <sstream>
#include <unordered_map> // si lees key-value

using namespace std;
using namespace std::chrono;

constexpr double PI = std::acos(-1.0);
constexpr double VEV = 246.0;  // Valor estándar, ajústalo si lo obtienes de SM


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

/* 
    double mh = 125.0;

    double inv = 1.0/std::sqrt(1.0 + tan_beta*tan_beta);
    double cb  = inv;
    double sb  = tan_beta * inv;
    double cba = std::sqrt(1.0 - sin_ba*sin_ba);

    double ca = cb * cba + sb * sin_ba;
    double sa = sb * cba - cb * sin_ba;
    double beta  = std::atan(tan_beta);
    double alpha = beta - std::asin(sin_ba);

*/

double calc_lambda1(double mh, double mH,
                          double m12_2,
                          double sa, double ca, double sb, double cb, double tan_beta,
                          double lam6, double lam7, double VEV) {
    // Términos de lambda 1
    double term1 = (mH*mH * ca*ca + mh*mh * sa*sa) / (VEV*VEV * cb*cb);
    double term2 = 1.5 * lam6 * tan_beta;
    double term3 = (m12_2)/(VEV*VEV) * tan_beta/(cb*cb);
    double term4 = 0.5 * lam7 * tan_beta*tan_beta*tan_beta;
    double l1    = term1 - term2 - term3 + term4;
    return l1;
}

double calc_lambda2(double mh, double mH,
                          double m12_2,
                          double sa, double ca, double sb, double cb, double tan_beta,
                          double lam6, double lam7, double VEV) {

    double ct = 1.0/tan_beta;

    // Términos de lambda_2
    double term1 = (mh*mh * ca*ca + mH*mH * sa*sa) / (VEV*VEV * sb*sb);
    double term2 = 1.5 * lam7 * ct;
    double term3 = 0.5 * lam6 * ct*ct*ct;
    double term4 = (m12_2)/(VEV*VEV) * ct/(sb*sb);
    double l2    = term1 - term2 + term3 - term4;
    return l2;
}

bool check_lambda(double l1) {
    return (l1*l1) < (16.0 * PI * PI);
}

double m12_base(double mh, double m_phi, double sa, double ca, double sb, double cb, double tan_beta, double VEV) {
    double m12_base_f2 = (std::pow(m_phi, 2)*std::pow(sa, 2) +
                    std::pow(mh, 2)*std::pow(ca, 2))*tan_beta
                    + (VEV*VEV/2)*(lambda6*std::pow(cb,2) - 3*lambda7*std::pow(sb,2));

    double m12_base_f1 = (std::pow(m_phi, 2)*std::pow(ca, 2) +
                    std::pow(mh, 2)*std::pow(sa, 2))/tan_beta
                    + (VEV*VEV/2)*(lambda7*std::pow(sb,2) - 3*lambda6*std::pow(cb,2));
    return (m12_base_f1 + m12_base_f2) / 2.0; // Promedio para simetría

}

// Scan del modelo en archivo json

static std::string load_file(const std::string &path) {
    std::ifstream in(path);
    if (!in.is_open())
        throw std::runtime_error("No pude abrir JSON: " + path);
    std::ostringstream buf;
    buf << in.rdbuf();
    return buf.str();
}

static double extract_double(const std::string &text, const std::string &key) {
    // Busca: "key" : número (entero o decimal, posible notación científica)
    std::regex re("\"" + key + "\"\\s*:\\s*([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)");
    std::smatch m;
    if (std::regex_search(text, m, re)) {
        try {
            return std::stod(m[1].str());
        } catch (...) {
            throw std::runtime_error("Conversión a double falló para clave: " + key);
        }
    }
    throw std::runtime_error("Clave no encontrada en JSON: " + key);
}

static int extract_int(const std::string &text, const std::string &key) {
    std::regex re("\"" + key + "\"\\s*:\\s*([0-9]+)");
    std::smatch m;
    if (std::regex_search(text, m, re)) {
        return std::stoi(m[1].str());
    }
    throw std::runtime_error("Clave no encontrada en JSON: " + key);
}

void parse_json_config(
    const std::string &filename,
    QuadraticModel    &qm,
    SearchSettings    &ss,
    FixedParameters   &fp
) {
    const std::string content = load_file(filename);

    // Modelo cuadrático
    qm.a = extract_double(content, "a");
    qm.b = extract_double(content, "b");
    qm.c = extract_double(content, "c");

    // search_settings
    ss.N_mphi   = extract_int   (content, "N_mphi");
    ss.N_m12    = extract_int   (content, "N_m12");
    ss.m12_min  = extract_double(content, "m12_min");
    ss.m12_max  = extract_double(content, "m12_max");

    // fixed_parameters
    fp.mA        = extract_double(content, "mA");
    fp.sin_ba    = extract_double(content, "sin_ba");
    fp.tan_beta  = extract_double(content, "tan_beta");
    fp.lambda6   = extract_double(content, "lambda6");
    fp.lambda7   = extract_double(content, "lambda7");
}
