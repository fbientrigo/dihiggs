#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"

#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

constexpr int kDefaultYukawasType = 1;
constexpr double kWidthFloor = 1e-15;
constexpr double kDefaultMh = 125.0;

double qnan() {
  return std::numeric_limits<double>::quiet_NaN();
}

std::string trim_copy(const std::string& s) {
  std::size_t begin = 0;
  while (begin < s.size() &&
         std::isspace(static_cast<unsigned char>(s[begin])) != 0) {
    ++begin;
  }

  std::size_t end = s.size();
  while (end > begin &&
         std::isspace(static_cast<unsigned char>(s[end - 1])) != 0) {
    --end;
  }

  return s.substr(begin, end - begin);
}

bool has_key(const std::unordered_map<std::string, std::string>& row,
             const std::string& key) {
  return row.find(key) != row.end();
}

bool has_column(const std::vector<std::string>& header, const std::string& key) {
  for (std::size_t i = 0; i < header.size(); ++i) {
    if (header[i] == key) {
      return true;
    }
  }
  return false;
}

std::string require_string(
    const std::unordered_map<std::string, std::string>& row,
    const std::string& key) {
  std::unordered_map<std::string, std::string>::const_iterator it = row.find(key);
  if (it == row.end()) {
    throw std::runtime_error("Missing required column: " + key);
  }

  const std::string value = trim_copy(it->second);
  if (value.empty()) {
    throw std::runtime_error("Empty required value for column: " + key);
  }
  return value;
}

double require_double(
    const std::unordered_map<std::string, std::string>& row,
    const std::string& key) {
  const std::string raw = require_string(row, key);
  const double parsed = std::stod(raw);
  if (!std::isfinite(parsed)) {
    throw std::runtime_error("Non-finite value for column: " + key + " value=" + raw);
  }
  return parsed;
}

double get_optional_double(
    const std::unordered_map<std::string, std::string>& row,
    const std::string& key,
    double fallback) {
  if (!has_key(row, key)) {
    return fallback;
  }
  const std::string raw = trim_copy(row.at(key));
  if (raw.empty()) {
    return fallback;
  }
  const double parsed = std::stod(raw);
  if (!std::isfinite(parsed)) {
    return fallback;
  }
  return parsed;
}

double require_alias_double(
    const std::unordered_map<std::string, std::string>& row,
    const std::string& primary_key,
    const std::string& secondary_key) {
  if (has_key(row, primary_key)) {
    return require_double(row, primary_key);
  }
  return require_double(row, secondary_key);
}

std::vector<std::string> split_csv_simple(const std::string& line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;

  while (std::getline(ss, item, ',')) {
    out.push_back(trim_copy(item));
  }

  if (!line.empty() && line.back() == ',') {
    out.push_back("");
  }

  return out;
}

std::unordered_map<std::string, std::string> make_row_map(
    const std::vector<std::string>& header,
    const std::vector<std::string>& values) {
  std::unordered_map<std::string, std::string> row;
  const std::size_t n = std::min(header.size(), values.size());
  for (std::size_t i = 0; i < n; ++i) {
    row[header[i]] = values[i];
  }
  return row;
}

void validate_header_contract(const std::vector<std::string>& header) {
  const bool has_mH = has_column(header, "mH") || has_column(header, "m_phi");
  const bool has_sba = has_column(header, "sba") || has_column(header, "sin_ba");
  const bool has_l1 = has_column(header, "lambda1") || has_column(header, "lam1");

  if (!has_mH) {
    throw std::runtime_error("Header must contain mH or m_phi");
  }
  if (!has_column(header, "mA")) {
    throw std::runtime_error("Header must contain mA");
  }
  if (!has_sba) {
    throw std::runtime_error("Header must contain sba or sin_ba");
  }
  if (!has_column(header, "tan_beta")) {
    throw std::runtime_error("Header must contain tan_beta");
  }
  if (!has_l1) {
    throw std::runtime_error("Header must contain lambda1 or lam1");
  }
  if (!has_column(header, "lambda6")) {
    throw std::runtime_error("Header must contain lambda6");
  }
  if (!has_column(header, "lambda7")) {
    throw std::runtime_error("Header must contain lambda7");
  }
}

struct ColabPointInput {
  std::string point_id;
  double mh;
  double mH;
  double mA;
  double mHp;
  double sba;
  double tan_beta;
  double lambda1;
  double lambda6;
  double lambda7;
  double m12sq_analytic;
  double br_bb_colab;
  double br_gaga_colab;
  double br_Zga_colab;
};

struct ValidationRow {
  std::string point_id;

  double mh;
  double mH;
  double mA;
  double mHp;
  double sba;
  double tan_beta;
  double lambda1;
  double lambda6;
  double lambda7;
  double m12sq_analytic;

  double br_bb_colab;
  double br_gaga_colab;
  double br_Zga_colab;

  int set_param_ok;
  int validation_available;

  double lambda1_input;
  double lambda1_recomputed;
  double abs_error;
  int warning_flag;

  int positivity_ok;
  int unitarity_ok;
  int perturbativity_ok;
  int triple_ok;

  double m12_2_g;
  double lam1_g;
  double lam2_g;
  double lam3_g;
  double lam4_g;
  double lam5_g;

  double width_bb;
  double width_cc;
  double width_tautau;
  double width_WW;
  double width_ZZ;
  double width_gaga;
  double width_Zga;
  double width_gg;
  double width_hh;
  double total_width;

  double br_bb_2hdmc;
  double br_gaga_2hdmc;
  double br_Zga_2hdmc;
};

ColabPointInput parse_colab_row(
    const std::unordered_map<std::string, std::string>& row,
    std::size_t row_index) {
  ColabPointInput out{};

  out.point_id = has_key(row, "point_id")
      ? require_string(row, "point_id")
      : ("row_" + std::to_string(row_index));

  out.mh = get_optional_double(row, "mh", kDefaultMh);
  out.mH = require_alias_double(row, "mH", "m_phi");
  out.mA = require_double(row, "mA");
  out.mHp = get_optional_double(row, "mHp", out.mA);

  out.sba = require_alias_double(row, "sba", "sin_ba");
  out.tan_beta = require_double(row, "tan_beta");
  out.lambda1 = require_alias_double(row, "lambda1", "lam1");
  out.lambda6 = require_double(row, "lambda6");
  out.lambda7 = require_double(row, "lambda7");

  out.m12sq_analytic = get_optional_double(row, "m12sq_analytic", qnan());
  out.br_bb_colab = get_optional_double(row, "br_bb_colab", qnan());
  out.br_gaga_colab = get_optional_double(row, "br_gaga_colab", qnan());
  out.br_Zga_colab = get_optional_double(row, "br_Zga_colab", qnan());

  return out;
}

ValidationRow make_default_output_row(const ColabPointInput& in) {
  const double nanv = qnan();
  ValidationRow out{};
  out.point_id = in.point_id;

  out.mh = in.mh;
  out.mH = in.mH;
  out.mA = in.mA;
  out.mHp = in.mHp;
  out.sba = in.sba;
  out.tan_beta = in.tan_beta;
  out.lambda1 = in.lambda1;
  out.lambda6 = in.lambda6;
  out.lambda7 = in.lambda7;
  out.m12sq_analytic = in.m12sq_analytic;

  out.br_bb_colab = in.br_bb_colab;
  out.br_gaga_colab = in.br_gaga_colab;
  out.br_Zga_colab = in.br_Zga_colab;

  out.set_param_ok = 0;
  out.validation_available = 0;

  out.lambda1_input = in.lambda1;
  out.lambda1_recomputed = nanv;
  out.abs_error = nanv;
  out.warning_flag = 0;

  out.positivity_ok = 0;
  out.unitarity_ok = 0;
  out.perturbativity_ok = 0;
  out.triple_ok = 0;

  out.m12_2_g = nanv;
  out.lam1_g = nanv;
  out.lam2_g = nanv;
  out.lam3_g = nanv;
  out.lam4_g = nanv;
  out.lam5_g = nanv;

  out.width_bb = nanv;
  out.width_cc = nanv;
  out.width_tautau = nanv;
  out.width_WW = nanv;
  out.width_ZZ = nanv;
  out.width_gaga = nanv;
  out.width_Zga = nanv;
  out.width_gg = nanv;
  out.width_hh = nanv;
  out.total_width = nanv;

  out.br_bb_2hdmc = nanv;
  out.br_gaga_2hdmc = nanv;
  out.br_Zga_2hdmc = nanv;

  return out;
}

void write_output_header(std::ofstream& out) {
  out
      << "point_id,"
      << "mh,mH,mA,mHp,sba,tan_beta,lambda1,lambda6,lambda7,m12sq_analytic,"
      << "br_bb_colab,br_gaga_colab,br_Zga_colab,"
      << "set_param_ok,validation_available,"
      << "lambda1_input,lambda1_recomputed,abs_error,warning_flag,"
      << "positivity_ok,unitarity_ok,perturbativity_ok,triple_ok,"
      << "m12_2_g,lam1_g,lam2_g,lam3_g,lam4_g,lam5_g,"
      << "width_bb,width_cc,width_tautau,width_WW,width_ZZ,width_gaga,width_Zga,width_gg,width_hh,total_width,"
      << "br_bb_2hdmc,br_gaga_2hdmc,br_Zga_2hdmc\n";
}

std::string serialize_output_row(const ValidationRow& r) {
  std::ostringstream out;
  out.setf(std::ios::scientific);
  out << std::setprecision(17)
      << r.point_id << ','
      << r.mh << ',' << r.mH << ',' << r.mA << ',' << r.mHp << ','
      << r.sba << ',' << r.tan_beta << ',' << r.lambda1 << ','
      << r.lambda6 << ',' << r.lambda7 << ',' << r.m12sq_analytic << ','
      << r.br_bb_colab << ',' << r.br_gaga_colab << ',' << r.br_Zga_colab << ','
      << r.set_param_ok << ',' << r.validation_available << ','
      << r.lambda1_input << ',' << r.lambda1_recomputed << ',' << r.abs_error << ',' << r.warning_flag << ','
      << r.positivity_ok << ',' << r.unitarity_ok << ',' << r.perturbativity_ok << ',' << r.triple_ok << ','
      << r.m12_2_g << ',' << r.lam1_g << ',' << r.lam2_g << ',' << r.lam3_g << ',' << r.lam4_g << ',' << r.lam5_g << ','
      << r.width_bb << ',' << r.width_cc << ',' << r.width_tautau << ',' << r.width_WW << ',' << r.width_ZZ << ','
      << r.width_gaga << ',' << r.width_Zga << ',' << r.width_gg << ',' << r.width_hh << ',' << r.total_width << ','
      << r.br_bb_2hdmc << ',' << r.br_gaga_2hdmc << ',' << r.br_Zga_2hdmc
      << '\n';
  return out.str();
}

bool compute_validation_row(
    const ColabPointInput& in,
    int yukawas_type,
    ValidationRow& out_row,
    std::string& error_detail) {
  out_row = make_default_output_row(in);
  error_detail.clear();

  try {
    THDM model;
    SM sm;
    model.set_SM(sm);
    model.set_yukawas_type(yukawas_type);

    const bool ok = model.set_param_phys_lam1(
        in.mh,
        in.mH,
        in.mA,
        in.mHp,
        in.sba,
        in.lambda1,
        in.lambda6,
        in.lambda7,
        in.tan_beta);

    out_row.set_param_ok = ok ? 1 : 0;

    if (!ok) {
      error_detail = "set_param_phys_lam1 returned false";
      return false;
    }

    if (model.has_param_phys_lam1_validation()) {
      double lambda1_input = out_row.lambda1_input;
      double lambda1_recomputed = qnan();
      double abs_error = qnan();
      bool warning_flag = false;

      model.get_param_phys_lam1_validation(
          lambda1_input,
          lambda1_recomputed,
          abs_error,
          warning_flag);

      out_row.validation_available = 1;
      out_row.lambda1_input = lambda1_input;
      out_row.lambda1_recomputed = lambda1_recomputed;
      out_row.abs_error = abs_error;
      out_row.warning_flag = warning_flag ? 1 : 0;
    }

    double lam1_g = qnan();
    double lam2_g = qnan();
    double lam3_g = qnan();
    double lam4_g = qnan();
    double lam5_g = qnan();
    double lam6_g = qnan();
    double lam7_g = qnan();
    double m12_2_g = qnan();
    double tanb_g = qnan();

    model.get_param_gen(
        lam1_g, lam2_g, lam3_g, lam4_g, lam5_g,
        lam6_g, lam7_g, m12_2_g, tanb_g);

    out_row.m12_2_g = m12_2_g;
    out_row.lam1_g = lam1_g;
    out_row.lam2_g = lam2_g;
    out_row.lam3_g = lam3_g;
    out_row.lam4_g = lam4_g;
    out_row.lam5_g = lam5_g;

    Constraints check(model);
    const bool pos = check.check_positivity();
    const bool uni = check.check_unitarity();
    const bool pert = check.check_perturbativity();

    out_row.positivity_ok = pos ? 1 : 0;
    out_row.unitarity_ok = uni ? 1 : 0;
    out_row.perturbativity_ok = pert ? 1 : 0;
    out_row.triple_ok = (pos && uni && pert) ? 1 : 0;

    DecayTable table(model);

    const double w_bb = table.get_gamma_hdd(2, 3, 3);
    const double w_cc = table.get_gamma_huu(2, 2, 2);
    const double w_tautau = table.get_gamma_hll(2, 3, 3);
    const double w_WW = table.get_gamma_hvv(2, 3);
    const double w_ZZ = table.get_gamma_hvv(2, 2);
    const double w_gaga = table.get_gamma_hgaga(2);
    const double w_Zga = table.get_gamma_hZga(2);
    const double w_gg = table.get_gamma_hgg(2);
    const double w_hh = table.get_gamma_hhh(2, 1, 1);
    const double w_tot = table.get_gammatot_h(2);

    out_row.width_bb = w_bb;
    out_row.width_cc = w_cc;
    out_row.width_tautau = w_tautau;
    out_row.width_WW = w_WW;
    out_row.width_ZZ = w_ZZ;
    out_row.width_gaga = w_gaga;
    out_row.width_Zga = w_Zga;
    out_row.width_gg = w_gg;
    out_row.width_hh = w_hh;
    out_row.total_width = w_tot;

    if (std::isfinite(w_tot) && w_tot > kWidthFloor) {
      if (std::isfinite(w_bb)) {
        out_row.br_bb_2hdmc = w_bb / w_tot;
      }
      if (std::isfinite(w_gaga)) {
        out_row.br_gaga_2hdmc = w_gaga / w_tot;
      }
      if (std::isfinite(w_Zga)) {
        out_row.br_Zga_2hdmc = w_Zga / w_tot;
      }
    }

    return true;
  } catch (const std::exception& ex) {
    error_detail = ex.what();
    return false;
  } catch (...) {
    error_detail = "unknown exception during recomputation";
    return false;
  }
}

int run_colab_recompute(
    const std::string& input_csv,
    const std::string& output_csv,
    int yukawas_type) {
  std::ifstream in(input_csv.c_str());
  if (!in) {
    std::cerr << "Failed to open input CSV: " << input_csv << "\n";
    return 2;
  }

  std::ofstream out(output_csv.c_str());
  if (!out) {
    std::cerr << "Failed to open output CSV: " << output_csv << "\n";
    return 3;
  }

  std::string header_line;
  if (!std::getline(in, header_line)) {
    std::cerr << "Input CSV is empty: " << input_csv << "\n";
    return 4;
  }

  const std::vector<std::string> header = split_csv_simple(header_line);
  try {
    validate_header_contract(header);
  } catch (const std::exception& ex) {
    std::cerr << "Invalid input header: " << ex.what() << "\n";
    return 5;
  }

  write_output_header(out);

  std::size_t n_rows_total = 0;
  std::size_t n_rows_written = 0;
  std::size_t n_parse_fail = 0;
  std::size_t n_set_param_ok = 0;
  std::size_t n_set_param_fail = 0;
  std::size_t n_validation_available = 0;
  std::size_t n_warning_flag = 0;
  std::size_t n_triple_ok = 0;
  std::size_t n_compute_exception = 0;

  std::string line;
  std::size_t line_number = 1;

  while (std::getline(in, line)) {
    ++line_number;
    const std::string trimmed = trim_copy(line);
    if (trimmed.empty()) {
      continue;
    }

    ++n_rows_total;

    std::vector<std::string> values = split_csv_simple(line);
    if (values.size() < header.size()) {
      values.resize(header.size(), "");
    }

    const std::unordered_map<std::string, std::string> row_map = make_row_map(header, values);

    try {
      const ColabPointInput input = parse_colab_row(row_map, n_rows_total);
      ValidationRow out_row;
      std::string error_detail;

      const bool success = compute_validation_row(
          input,
          yukawas_type,
          out_row,
          error_detail);

      out << serialize_output_row(out_row);
      ++n_rows_written;

      if (out_row.set_param_ok == 1) {
        ++n_set_param_ok;
      } else {
        ++n_set_param_fail;
      }
      if (out_row.validation_available == 1) {
        ++n_validation_available;
      }
      if (out_row.warning_flag == 1) {
        ++n_warning_flag;
      }
      if (out_row.triple_ok == 1) {
        ++n_triple_ok;
      }
      if (!success && !error_detail.empty()) {
        ++n_compute_exception;
        std::cerr << "[row " << line_number << "] compute warning: "
                  << error_detail << "\n";
      }
    } catch (const std::exception& ex) {
      ++n_parse_fail;
      std::cerr << "[row " << line_number << "] parse error: "
                << ex.what() << "\n";
    } catch (...) {
      ++n_parse_fail;
      std::cerr << "[row " << line_number << "] unknown parse error\n";
    }
  }

  std::cout << std::setprecision(17)
            << "Input CSV: " << input_csv << "\n"
            << "Output CSV: " << output_csv << "\n"
            << "Rows seen: " << n_rows_total << "\n"
            << "Rows written: " << n_rows_written << "\n"
            << "Parse failures: " << n_parse_fail << "\n"
            << "set_param_phys_lam1 OK: " << n_set_param_ok << "\n"
            << "set_param_phys_lam1 FAIL: " << n_set_param_fail << "\n"
            << "validation_available: " << n_validation_available << "\n"
            << "warning_flag rows: " << n_warning_flag << "\n"
            << "triple_ok rows: " << n_triple_ok << "\n"
            << "rows with compute warnings/exceptions: " << n_compute_exception << "\n"
            << "TRIPLE_OK_POINTS " << n_triple_ok << "\n";

  return 0;
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 3 && argc != 4) {
    std::cerr
        << "Usage: " << argv[0]
        << " input_colab.csv output_validated.csv [yukawas_type]\n";
    return 1;
  }

  const std::string input_csv = argv[1];
  const std::string output_csv = argv[2];
  const int yukawas_type = (argc == 4) ? std::atoi(argv[3]) : kDefaultYukawasType;

  return run_colab_recompute(input_csv, output_csv, yukawas_type);
}