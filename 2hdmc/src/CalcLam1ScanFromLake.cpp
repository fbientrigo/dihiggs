#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

//error 231
#include <gsl/gsl_errno.h>

namespace {

constexpr double kMh = 125.0;
constexpr int kDefaultYukawasType = 1;
constexpr double kWidthFloor = 1e-15;

struct Sample {
  std::size_t sample_index;
  std::size_t attempt_index;
  double mh;
  double mH;
  double mA;
  double mHp;
  double sba;
  double lambda6;
  double lambda7;
  double m12_2;
  double tan_beta;
  double lambda1_input;
  double lambda1_recomputed;
  double abs_error;
  bool warning_flag;
};

struct QuarantinePointInput {
  std::string point_id;
  double m_phi;
  double mA;
  double sin_ba;
  double tan_beta;
  double lambda6;
  double lambda7;
  double lam1;
};

struct RecomputeRow {
  double m_phi;
  double mA;
  double alpha;
  double beta;
  double lambda6;
  double lambda7;
  double m12;
  double sin_ba;
  double tan_beta;
  int positivity_ok;
  int unitarity_ok;
  int perturbativity_ok;
  double width_bb;
  double width_tautau;
  double width_WW;
  double width_ZZ;
  double width_gaga;
  double width_Zga;
  double width_gg;
  double width_hh;
  double total_width;
  double br_gaga;
  double lam1;
  double computed_lam1;
  double lam2;
  double computed_lam2;
  double lam3;
  double lam4;
  double lam5;
};

double qnan() {
  return std::numeric_limits<double>::quiet_NaN();
}

std::string trim_copy(const std::string& s) {
  std::size_t begin = 0;
  while (begin < s.size() && std::isspace(static_cast<unsigned char>(s[begin])) != 0) {
    ++begin;
  }
  std::size_t end = s.size();
  while (end > begin && std::isspace(static_cast<unsigned char>(s[end - 1])) != 0) {
    --end;
  }
  return s.substr(begin, end - begin);
}

bool has_key(const std::unordered_map<std::string, std::string>& row,
             const std::string& key) {
  return row.find(key) != row.end();
}

std::string require_string(const std::unordered_map<std::string, std::string>& row,
                           const std::string& key) {
  std::unordered_map<std::string, std::string>::const_iterator it = row.find(key);
  if (it == row.end()) {
    throw std::runtime_error("Missing required column: " + key);
  }
  return trim_copy(it->second);
}

double require_double(const std::unordered_map<std::string, std::string>& row,
                      const std::string& key) {
  const std::string raw = require_string(row, key);
  const double parsed = std::stod(raw);
  if (!std::isfinite(parsed)) {
    throw std::runtime_error("Non-finite value for column: " + key + " value=" + raw);
  }
  return parsed;
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

std::unordered_map<std::string, std::string>
make_row_map(const std::vector<std::string>& header,
             const std::vector<std::string>& values) {
  std::unordered_map<std::string, std::string> row;
  const std::size_t n = std::min(header.size(), values.size());
  for (std::size_t i = 0; i < n; ++i) {
    row[header[i]] = values[i];
  }
  return row;
}

void write_random_header(std::ofstream &out) {
  out << "sample_index,attempt_index,mh,mH,mA,mHp,sba,lambda6,lambda7,m12_2,tan_beta,"
         "lambda1_input,lambda1_recomputed,abs_error,warning_flag\n";
}

void write_random_sample(std::ofstream &out, const Sample &sample) {
  out << sample.sample_index << ','
      << sample.attempt_index << ','
      << std::setprecision(17)
      << sample.mh << ','
      << sample.mH << ','
      << sample.mA << ','
      << sample.mHp << ','
      << sample.sba << ','
      << sample.lambda6 << ','
      << sample.lambda7 << ','
      << sample.m12_2 << ','
      << sample.tan_beta << ','
      << sample.lambda1_input << ','
      << sample.lambda1_recomputed << ','
      << sample.abs_error << ','
      << (sample.warning_flag ? 1 : 0) << '\n';
}

void write_recompute_header(std::ofstream& out) {
  out << "m_phi,mA,alpha,beta,lambda6,lambda7,m12,"
         "sin_ba,tan_beta,positivity_ok,unitarity_ok,perturbativity_ok,"
         "width_bb,width_tautau,width_WW,width_ZZ,"
         "width_gaga,width_Zga,width_gg,width_hh,"
         "total_width,br_gaga,lam1,computed_lam1,"
         "lam2,computed_lam2,lam3,lam4,lam5\n";
}

std::string serialize_recompute_row(const RecomputeRow& row) {
  std::ostringstream out;
  out.setf(std::ios::scientific);
  out << std::setprecision(17)
      << row.m_phi << ','
      << row.mA << ','
      << row.alpha << ','
      << row.beta << ','
      << row.lambda6 << ','
      << row.lambda7 << ','
      << row.m12 << ','
      << row.sin_ba << ','
      << row.tan_beta << ','
      << row.positivity_ok << ','
      << row.unitarity_ok << ','
      << row.perturbativity_ok << ','
      << row.width_bb << ','
      << row.width_tautau << ','
      << row.width_WW << ','
      << row.width_ZZ << ','
      << row.width_gaga << ','
      << row.width_Zga << ','
      << row.width_gg << ','
      << row.width_hh << ','
      << row.total_width << ','
      << row.br_gaga << ','
      << row.lam1 << ','
      << row.computed_lam1 << ','
      << row.lam2 << ','
      << row.computed_lam2 << ','
      << row.lam3 << ','
      << row.lam4 << ','
      << row.lam5 << '\n';
  return out.str();
}

double safe_asin_input(double value) {
  if (value < -1.0) {
    return -1.0;
  }
  if (value > 1.0) {
    return 1.0;
  }
  return value;
}

RecomputeRow make_failure_row(const QuarantinePointInput& in) {
  const double beta = std::atan(in.tan_beta);
  const double alpha = beta - std::asin(safe_asin_input(in.sin_ba));
  const double nanv = qnan();
  return RecomputeRow{
      in.m_phi,
      in.mA,
      alpha,
      beta,
      in.lambda6,
      in.lambda7,
      nanv,
      in.sin_ba,
      in.tan_beta,
      0,
      0,
      0,
      nanv,
      nanv,
      nanv,
      nanv,
      nanv,
      nanv,
      nanv,
      nanv,
      nanv,
      nanv,
      in.lam1,
      nanv,
      nanv,
      nanv,
      nanv,
      nanv,
      nanv};
}

QuarantinePointInput parse_quarantine_row(
    const std::unordered_map<std::string, std::string>& row,
    std::size_t row_index) {
  QuarantinePointInput out{};
  out.point_id = has_key(row, "point_id") ? require_string(row, "point_id")
                                           : ("row_" + std::to_string(row_index));
  out.m_phi = require_double(row, "m_phi");
  out.mA = require_double(row, "mA");
  out.sin_ba = require_double(row, "sin_ba");
  out.tan_beta = require_double(row, "tan_beta");
  out.lambda6 = require_double(row, "lambda6");
  out.lambda7 = require_double(row, "lambda7");
  out.lam1 = require_double(row, "lam1");
  return out;
}

bool compute_quarantine_point(
    const QuarantinePointInput& in,
    int yukawas_type,
    RecomputeRow& out_row,
    bool& set_param_ok,
    bool& triple_ok,
    std::string& error_detail) {
  out_row = make_failure_row(in);
  set_param_ok = false;
  triple_ok = false;
  error_detail.clear();

  try {
    THDM model;
    SM sm;
    model.set_SM(sm);
    model.set_yukawas_type(yukawas_type);

    const bool ok = model.set_param_phys_lam1(
        kMh, in.m_phi, in.mA, in.mA, in.sin_ba, in.lam1, in.lambda6, in.lambda7, in.tan_beta);
    if (!ok) {
      error_detail = "set_param_phys_lam1 returned false";
      return false;
    }
    set_param_ok = true;
//------------------------------
    // double lam1_g = qnan();
    // double lam2_g = qnan();
    // double lam3_g = qnan();
    // double lam4_g = qnan();
    // double lam5_g = qnan();
    // double lam6_g = qnan();
    // double lam7_g = qnan();
    // double m12_2_g = qnan();
    // double tanb_g = qnan();
    // model.get_param_gen(
    //     lam1_g, lam2_g, lam3_g, lam4_g, lam5_g,
    //     lam6_g, lam7_g, m12_2_g, tanb_g);

    // Constraints check(model);
    // const bool pos = check.check_positivity();
    // const bool uni = check.check_unitarity();
    // const bool pert = check.check_perturbativity();
    // triple_ok = (pos && uni && pert);

    // DecayTable tab(model);
    // const double w_bb = tab.get_gamma_hdd(2, 3, 3);
    // const double w_tautau = tab.get_gamma_hll(2, 3, 3);
    // const double w_WW = tab.get_gamma_hvv(2, 3);
    // const double w_ZZ = tab.get_gamma_hvv(2, 2);
    // const double w_gaga = tab.get_gamma_hgaga(2);
    // const double w_Zga = tab.get_gamma_hZga(2);
    // const double w_gg = tab.get_gamma_hgg(2);
    // const double w_hh = tab.get_gamma_hhh(2, 1, 1);
    // const double w_tot = tab.get_gammatot_h(2);
    // const double br_gaga = (w_tot > kWidthFloor) ? (w_gaga / w_tot) : qnan();
// -----------------------------

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

    const bool gen_ok =
        std::isfinite(lam1_g) &&
        std::isfinite(lam2_g) &&
        std::isfinite(lam3_g) &&
        std::isfinite(lam4_g) &&
        std::isfinite(lam5_g) &&
        std::isfinite(m12_2_g);

    if (!gen_ok) {
      error_detail = "non-finite values returned by get_param_gen";
      return false;
    }

    Constraints check(model);
    const bool pos = check.check_positivity();
    const bool uni = check.check_unitarity();
    const bool pert = check.check_perturbativity();
    triple_ok = (pos && uni && pert);

    DecayTable tab(model);
    const double w_bb = tab.get_gamma_hdd(2, 3, 3);
    const double w_tautau = tab.get_gamma_hll(2, 3, 3);
    const double w_WW = tab.get_gamma_hvv(2, 3);
    const double w_ZZ = tab.get_gamma_hvv(2, 2);
    const double w_gaga = tab.get_gamma_hgaga(2);
    const double w_Zga = tab.get_gamma_hZga(2);
    const double w_gg = tab.get_gamma_hgg(2);
    const double w_hh = tab.get_gamma_hhh(2, 1, 1);
    const double w_tot = tab.get_gammatot_h(2);
    double br_gaga = 0.0;
    if (std::isfinite(w_tot) && w_tot > 0.0 && std::isfinite(w_gaga) && w_gaga >= 0.0) {
        br_gaga = w_gaga / w_tot;
    }

    const bool widths_ok =
        std::isfinite(w_bb) &&
        std::isfinite(w_tautau) &&
        std::isfinite(w_WW) &&
        std::isfinite(w_ZZ) &&
        std::isfinite(w_gaga) &&
        std::isfinite(w_Zga) &&
        std::isfinite(w_gg) &&
        std::isfinite(w_hh) &&
        std::isfinite(w_tot);

    if (!widths_ok) {
      error_detail = "non-finite decay widths produced during recomputation";
      return false;
    }




// ----------- end corrected------------------
    const double beta = std::atan(in.tan_beta);
    const double alpha = beta - std::asin(safe_asin_input(in.sin_ba));

    out_row = RecomputeRow{
        in.m_phi,
        in.mA,
        alpha,
        beta,
        in.lambda6,
        in.lambda7,
        m12_2_g,
        in.sin_ba,
        in.tan_beta,
        pos ? 1 : 0,
        uni ? 1 : 0,
        pert ? 1 : 0,
        w_bb,
        w_tautau,
        w_WW,
        w_ZZ,
        w_gaga,
        w_Zga,
        w_gg,
        w_hh,
        w_tot,
        br_gaga,
        in.lam1,
        lam1_g,
        lam2_g,
        lam2_g,
        lam3_g,
        lam4_g,
        lam5_g};
    return true;
  } catch (const std::exception& ex) {
    error_detail = ex.what();
    return false;
  } catch (...) {
    error_detail = "unknown error during point recomputation";
    return false;
  }
}

int run_random_mode(const std::string& output_csv,
                    std::size_t n_samples,
                    unsigned long long seed) {
  const std::size_t max_attempts = std::max<std::size_t>(n_samples * 200, 1000);

  std::ofstream out(output_csv.c_str());
  if (!out) {
    std::cerr << "Failed to open output file: " << output_csv << "\n";
    return 2;
  }

  write_random_header(out);

  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> mh_dist(120.0, 130.0);
  std::uniform_real_distribution<double> mH_dist(200.0, 600.0);
  std::uniform_real_distribution<double> mA_dist(200.0, 700.0);
  std::uniform_real_distribution<double> mHp_dist(200.0, 700.0);
  std::uniform_real_distribution<double> sba_dist(-0.9999, 0.9999);
  std::uniform_real_distribution<double> lambda6_dist(0.0, 1.0);
  std::uniform_real_distribution<double> lambda7_dist(0.0, 1.0);
  std::uniform_real_distribution<double> m12_dist(1.0, 500.0);
  std::uniform_real_distribution<double> tan_beta_dist(100.0, 50000.0);

  std::size_t accepted = 0;
  std::size_t warnings = 0;
  std::size_t attempts = 0;
  double max_abs_error = 0.0;
  Sample worst_sample{};

  while (accepted < n_samples && attempts < max_attempts) {
    ++attempts;

    const double mh = mh_dist(rng);
    const double mH = std::max(mH_dist(rng), mh + 1e-6);
    const double mA = mA_dist(rng);
    const double mHp = mHp_dist(rng);
    const double sba = sba_dist(rng);
    const double lambda6 = lambda6_dist(rng);
    const double lambda7 = lambda7_dist(rng);
    const double m12_2 = m12_dist(rng);
    const double tan_beta = tan_beta_dist(rng);

    THDM baseline;
    if (!baseline.set_param_phys(mh, mH, mA, mHp, sba, lambda6, lambda7, m12_2, tan_beta)) {
      continue;
    }

    double lambda1 = 0.0;
    double lambda2 = 0.0;
    double lambda3 = 0.0;
    double lambda4 = 0.0;
    double lambda5 = 0.0;
    double lambda6_rt = 0.0;
    double lambda7_rt = 0.0;
    double m12_rt = 0.0;
    double tan_beta_rt = 0.0;
    baseline.get_param_gen(lambda1, lambda2, lambda3, lambda4, lambda5,
                           lambda6_rt, lambda7_rt, m12_rt, tan_beta_rt);

    THDM probe;
    if (!probe.set_param_phys_lam1(mh, mH, mA, mHp, sba, lambda1, lambda6, lambda7, tan_beta)) {
      continue;
    }

    if (!probe.has_param_phys_lam1_validation()) {
      std::cerr << "Missing validation state for accepted sample " << accepted << "\n";
      return 3;
    }

    double lambda1_input = 0.0;
    double lambda1_recomputed = 0.0;
    double abs_error = 0.0;
    bool warning_flag = false;
    probe.get_param_phys_lam1_validation(lambda1_input,
                                         lambda1_recomputed,
                                         abs_error,
                                         warning_flag);

    Sample sample{accepted, attempts - 1, mh, mH, mA, mHp, sba, lambda6,
                  lambda7, m12_2, tan_beta, lambda1_input,
                  lambda1_recomputed, abs_error, warning_flag};
    write_random_sample(out, sample);

    if (warning_flag) {
      ++warnings;
    }
    if (abs_error > max_abs_error) {
      max_abs_error = abs_error;
      worst_sample = sample;
    }

    ++accepted;
  }

  if (accepted != n_samples) {
    std::cerr << "Failed to generate requested samples. accepted=" << accepted
              << " requested=" << n_samples << " attempts=" << attempts << "\n";
    return 4;
  }

  std::cout << std::setprecision(17)
            << "Wrote " << accepted << " samples to " << output_csv << "\n"
            << "Seed: " << seed << "\n"
            << "Attempts: " << attempts << "\n"
            << "Warnings: " << warnings << "\n"
            << "Warning rate: "
            << static_cast<double>(warnings) / static_cast<double>(accepted) << "\n"
            << "Max abs error: " << max_abs_error << "\n"
            << "Worst sample index: " << worst_sample.sample_index << "\n"
            << "Worst lambda1 input: " << worst_sample.lambda1_input << "\n"
            << "Worst lambda1 recomputed: " << worst_sample.lambda1_recomputed << "\n";

  return 0;
}

int run_single_point_mode(double mh,
                          double mH,
                          double mA,
                          double mHp,
                          double sba,
                          double lambda6,
                          double lambda7,
                          double lambda1_input,
                          double tan_beta,
                          int yukawas_type) {
  THDM model;
  SM sm;
  model.set_SM(sm);
  model.set_yukawas_type(yukawas_type);

  const bool ok = model.set_param_phys_lam1(
      mh, mH, mA, mHp, sba, lambda1_input, lambda6, lambda7, tan_beta);
  std::cout << std::setprecision(17);
  std::cout << "set_param_ok=" << (ok ? 1 : 0) << "\n";
  if (!ok) {
    return 2;
  }

  Constraints check(model);
  const bool pos = check.check_positivity();
  const bool uni = check.check_unitarity();
  const bool pert = check.check_perturbativity();
  std::cout << "positivity_ok=" << (pos ? 1 : 0) << "\n";
  std::cout << "unitarity_ok=" << (uni ? 1 : 0) << "\n";
  std::cout << "perturbativity_ok=" << (pert ? 1 : 0) << "\n";
  std::cout << "TRIPLE_OK_POINTS " << ((pos && uni && pert) ? 1 : 0) << "\n";

  double lam1_g = qnan();
  double lam2_g = qnan();
  double lam3_g = qnan();
  double lam4_g = qnan();
  double lam5_g = qnan();
  double lam6_g = qnan();
  double lam7_g = qnan();
  double m12_2_g = qnan();
  double tanb_g = qnan();
  model.get_param_gen(lam1_g, lam2_g, lam3_g, lam4_g, lam5_g,
                      lam6_g, lam7_g, m12_2_g, tanb_g);
  std::cout << "m12=" << m12_2_g << "\n";
  std::cout << "computed_lam1=" << lam1_g << "\n";
  std::cout << "computed_lam2=" << lam2_g << "\n";

  DecayTable table(model);
  table.print_decays(2);
  return 0;
}

int run_quarantine_mode(const std::string& input_csv,
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
  std::vector<QuarantinePointInput> inputs;
  std::string line;
  std::size_t row_index = 0;
  while (std::getline(in, line)) {
    const std::string trimmed = trim_copy(line);
    if (trimmed.empty()) {
      continue;
    }
    const std::vector<std::string> values = split_csv_simple(trimmed);
    const std::unordered_map<std::string, std::string> row = make_row_map(header, values);
    inputs.push_back(parse_quarantine_row(row, row_index));
    ++row_index;
  }

  if (inputs.empty()) {
    std::cerr << "Input CSV contains 0 data rows: " << input_csv << "\n";
    return 5;
  }

  write_recompute_header(out);

  std::vector<std::string> serialized_rows(inputs.size());
  std::vector<std::string> error_messages(inputs.size());

  long long set_param_ok_count = 0;
  long long set_param_fail_count = 0;
  long long triple_ok_count = 0;

#pragma omp parallel for schedule(dynamic) reduction(+:set_param_ok_count,set_param_fail_count,triple_ok_count)
  for (long long i = 0; i < static_cast<long long>(inputs.size()); ++i) {
    RecomputeRow row{};
    bool set_param_ok = false;
    bool triple_ok = false;
    std::string error_detail;
    const bool success = compute_quarantine_point(
        inputs[static_cast<std::size_t>(i)], yukawas_type, row, set_param_ok, triple_ok, error_detail);
    if (set_param_ok) {
      ++set_param_ok_count;
    } else {
      ++set_param_fail_count;
      if (!error_detail.empty()) {
        error_messages[static_cast<std::size_t>(i)] =
            "row=" + inputs[static_cast<std::size_t>(i)].point_id + " " + error_detail;
      }
    }
    if (triple_ok) {
      ++triple_ok_count;
    }
    serialized_rows[static_cast<std::size_t>(i)] = serialize_recompute_row(row);
    (void)success;
  }

  for (std::size_t i = 0; i < serialized_rows.size(); ++i) {
    out << serialized_rows[i];
  }
  out.close();

  // old warnings
  // for (std::size_t i = 0; i < error_messages.size(); ++i) {
  //   if (!error_messages[i].empty()) {
  //     std::cerr << error_messages[i] << "\n";
  //   }
  // }
  // better
    std::size_t printed_errors = 0;
  const std::size_t kMaxPrintedErrors = 20;

  for (std::size_t i = 0; i < error_messages.size(); ++i) {
    if (error_messages[i].empty()) {
      continue;
    }
    if (printed_errors < kMaxPrintedErrors) {
      std::cerr << error_messages[i] << "\n";
    }
    ++printed_errors;
  }

  if (printed_errors > kMaxPrintedErrors) {
    std::cerr << "suppressed "
              << (printed_errors - kMaxPrintedErrors)
              << " additional row-level errors\n";
  }


  const long long total_rows = static_cast<long long>(inputs.size());
  std::cout << std::setprecision(17)
            << "Processed rows: " << total_rows << "\n"
            << "Total Attempts: " << total_rows << "\n"
            << "Total CSV Rows: " << total_rows << "\n"
            << "set_param_phys_lam1 OK: " << set_param_ok_count << "\n"
            << "set_param_phys_lam1 FAIL: " << set_param_fail_count << "\n"
            << "TRIPLE_OK_POINTS " << triple_ok_count << "\n";
  return 0;
}

void print_usage(const char* argv0) {
  std::cout
      << "Usage:\n"
      << "  Legacy random mode (backward compatible):\n"
      << "    " << argv0 << " output_csv [n_samples] [seed]\n\n"
      << "  Explicit random mode:\n"
      << "    " << argv0 << " random output_csv [n_samples] [seed]\n\n"
      << "  Quarantine recompute mode:\n"
      << "    " << argv0 << " quarantine input_points.csv output_scan.csv [yukawas_type]\n\n"
      << "  Single-point debug mode:\n"
      << "    " << argv0 << " single mh mH mA mHp sba lambda6 lambda7 lambda1 tan_beta [yukawas_type]\n";
}

}  // namespace

int main(int argc, char *argv[]) {
  try {
    // error 231
    gsl_set_error_handler_off();
    if (argc >= 2 && argc <= 4 &&
        std::string(argv[1]) != "random" &&
        std::string(argv[1]) != "quarantine" &&
        std::string(argv[1]) != "single") {
      const std::string output_csv = argv[1];
      const std::size_t n_samples =
          (argc >= 3) ? static_cast<std::size_t>(std::strtoull(argv[2], NULL, 10))
                      : 5000;
      const unsigned long long seed =
          (argc >= 4) ? std::strtoull(argv[3], NULL, 10)
                      : 123456ULL;
      return run_random_mode(output_csv, n_samples, seed);
    }

    if (argc >= 2 && std::string(argv[1]) == "random") {
      if (argc < 3 || argc > 5) {
        print_usage(argv[0]);
        return 1;
      }
      const std::string output_csv = argv[2];
      const std::size_t n_samples =
          (argc >= 4) ? static_cast<std::size_t>(std::strtoull(argv[3], NULL, 10))
                      : 5000;
      const unsigned long long seed =
          (argc >= 5) ? std::strtoull(argv[4], NULL, 10)
                      : 123456ULL;
      return run_random_mode(output_csv, n_samples, seed);
    }

    if (argc >= 4 && std::string(argv[1]) == "quarantine") {
      if (argc != 4 && argc != 5) {
        print_usage(argv[0]);
        return 1;
      }
      const int yukawas_type = (argc == 5) ? std::stoi(argv[4]) : kDefaultYukawasType;
      return run_quarantine_mode(argv[2], argv[3], yukawas_type);
    }

    if (argc >= 11 && std::string(argv[1]) == "single") {
      if (argc != 11 && argc != 12) {
        print_usage(argv[0]);
        return 1;
      }
      const double mh = std::stod(argv[2]);
      const double mH = std::stod(argv[3]);
      const double mA = std::stod(argv[4]);
      const double mHp = std::stod(argv[5]);
      const double sba = std::stod(argv[6]);
      const double lambda6 = std::stod(argv[7]);
      const double lambda7 = std::stod(argv[8]);
      const double lambda1_input = std::stod(argv[9]);
      const double tan_beta = std::stod(argv[10]);
      const int yukawas_type = (argc == 12) ? std::stoi(argv[11]) : kDefaultYukawasType;
      return run_single_point_mode(
          mh, mH, mA, mHp, sba, lambda6, lambda7, lambda1_input, tan_beta, yukawas_type);
    }

    print_usage(argv[0]);
    return 1;
  } catch (const std::exception& ex) {
    std::cerr << "Fatal error: " << ex.what() << "\n";
    return 10;
  } catch (...) {
    std::cerr << "Fatal unknown error\n";
    return 11;
  }
}
