#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

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

double qnan() {
  return std::numeric_limits<double>::quiet_NaN();
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
  return it->second;
}

double require_double(const std::unordered_map<std::string, std::string>& row,
                      const std::string& key) {
  return std::stod(require_string(row, key));
}

double optional_double(const std::unordered_map<std::string, std::string>& row,
                       const std::string& key,
                       double default_value) {
  std::unordered_map<std::string, std::string>::const_iterator it = row.find(key);
  if (it == row.end() || it->second.empty()) {
    return default_value;
  }
  return std::stod(it->second);
}

std::vector<std::string> split_csv_simple(const std::string& line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) {
    out.push_back(item);
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

void write_points_header(std::ofstream& out) {
  out << "point_id,mh,mH,mA,mHp,sba,tan_beta,lambda1,lambda6,lambda7,"
         "m12sq_analytic,br_bb_colab,br_gaga_colab,br_Zga_colab,"
         "set_param_ok,validation_available,lambda1_input,lambda1_recomputed,abs_error,warning_flag,"
         "positivity_ok,unitarity_ok,perturbativity_ok,triple_ok,"
         "m12_2_g,lam1_g,lam2_g,lam3_g,lam4_g,lam5_g,"
         "width_bb,width_cc,width_tautau,width_WW,width_ZZ,width_gaga,width_Zga,width_gg,width_hh,total_width,"
         "br_bb_2hdmc,br_gaga_2hdmc,br_Zga_2hdmc\n";
}

void write_points_failure_row(
    std::ofstream& out,
    const std::string& point_id,
    double mh,
    double mH,
    double mA,
    double mHp,
    double sba,
    double tan_beta,
    double lambda1,
    double lambda6,
    double lambda7,
    double m12sq_analytic,
    double br_bb_colab,
    double br_gaga_colab,
    double br_Zga_colab) {
  const double nanv = qnan();
  out << std::setprecision(17)
      << point_id << ','
      << mh << ','
      << mH << ','
      << mA << ','
      << mHp << ','
      << sba << ','
      << tan_beta << ','
      << lambda1 << ','
      << lambda6 << ','
      << lambda7 << ','
      << m12sq_analytic << ','
      << br_bb_colab << ','
      << br_gaga_colab << ','
      << br_Zga_colab << ','
      << 0 << ','   // set_param_ok
      << 0 << ','   // validation_available
      << nanv << ','
      << nanv << ','
      << nanv << ','
      << 0 << ','   // warning_flag
      << 0 << ','
      << 0 << ','
      << 0 << ','
      << 0 << ','
      << nanv << ','
      << nanv << ','
      << nanv << ','
      << nanv << ','
      << nanv << ','
      << nanv << ','
      << nanv << ','  // width_bb
      << nanv << ','  // width_cc
      << nanv << ','  // width_tautau
      << nanv << ','  // width_WW
      << nanv << ','  // width_ZZ
      << nanv << ','  // width_gaga
      << nanv << ','  // width_Zga
      << nanv << ','  // width_gg
      << nanv << ','  // width_hh
      << nanv << ','  // total_width
      << nanv << ','  // br_bb_2hdmc
      << nanv << ','  // br_gaga_2hdmc
      << nanv << '\n'; // br_Zga_2hdmc
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
            << "Warning rate: " << static_cast<double>(warnings) / static_cast<double>(accepted) << "\n"
            << "Max abs error: " << max_abs_error << "\n"
            << "Worst sample index: " << worst_sample.sample_index << "\n"
            << "Worst lambda1 input: " << worst_sample.lambda1_input << "\n"
            << "Worst lambda1 recomputed: " << worst_sample.lambda1_recomputed << "\n";

  return 0;
}

int run_single_point_mode(float mh, float mH, float mA, float mHp, float sba, float lambda6, float lambda7,
                          float lambda1_input, float tan_beta, int yukawas_type) {
      // run single point mode
    // ./CalcLam1Scan mh mH mA mHp sin(beta-alpha) lambda_6 lambda_7 lambda_1 tan_beta yukawas_type output_filename
    // ./CalcLam1Scan 125 130.7018 300 300 1.0 0.01 0.0 1.7082857205846682 10000 1 worstpoint.lhe
    // 
    // if (argc == 12) {
    //   return run_single_point_mode(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11]);
    // }

  
    THDM model;
    SM sm;
    model.set_SM(sm);
    model.set_yukawas_type(1);

    const bool ok = model.set_param_phys_lam1(mh, mH, mA, mHp, sba, lambda1_input, lambda6, lambda7, tan_beta);

    Constraints check(model);
    const bool pos  = check.check_positivity();
    const bool uni  = check.check_unitarity();
    const bool pert = check.check_perturbativity();
    const bool triple_ok = (pos && uni && pert);

    DecayTable table(model);

    table.print_decays(1);
    table.print_decays(2);
    table.print_decays(3);
    table.print_decays(4);
                            
}

int run_points_mode(const std::string& input_csv,
                    const std::string& output_csv) {
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
  write_points_header(out);

  std::size_t total_rows = 0;
  std::size_t set_param_ok_count = 0;
  std::size_t triple_ok_count = 0;

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }

    const std::vector<std::string> values = split_csv_simple(line);
    std::unordered_map<std::string, std::string> row = make_row_map(header, values);

    const std::string point_id = has_key(row, "point_id")
        ? require_string(row, "point_id")
        : ("row_" + std::to_string(total_rows));

    const double mh      = require_double(row, "mh");
    const double mH      = require_double(row, "mH");
    const double mA      = require_double(row, "mA");
    const double mHp     = require_double(row, "mHp");
    const double sba     = require_double(row, "sba");
    const double tanb    = require_double(row, "tan_beta");
    const double lambda1 = require_double(row, "lambda1");
    const double lambda6 = require_double(row, "lambda6");
    const double lambda7 = require_double(row, "lambda7");

    const double m12sq_analytic = optional_double(row, "m12sq_analytic", qnan());
    const double br_bb_colab    = optional_double(row, "br_bb_colab", qnan());
    const double br_gaga_colab  = optional_double(row, "br_gaga_colab", qnan());
    const double br_Zga_colab   = optional_double(row, "br_Zga_colab", qnan());

    THDM model;
    SM sm;
    model.set_SM(sm);
    model.set_yukawas_type(1);

    const bool ok = model.set_param_phys_lam1(
        mh, mH, mA, mHp, sba, lambda1, lambda6, lambda7, tanb);

    bool validation_available = model.has_param_phys_lam1_validation();
    double lambda1_input = qnan();
    double lambda1_recomputed = qnan();
    double abs_error = qnan();
    bool warning_flag = false;

    if (validation_available) {
      model.get_param_phys_lam1_validation(lambda1_input,
                                           lambda1_recomputed,
                                           abs_error,
                                           warning_flag);
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

    model.get_param_gen(lam1_g, lam2_g, lam3_g, lam4_g, lam5_g,
                        lam6_g, lam7_g, m12_2_g, tanb_g);

    Constraints check(model);
    const bool pos  = check.check_positivity();
    const bool uni  = check.check_unitarity();
    const bool pert = check.check_perturbativity();
    const bool triple_ok = (pos && uni && pert);

    if (triple_ok) {
      ++triple_ok_count;
    }

    DecayTable tab(model);
    // TODO: Observate after using the worst point comparison
    const double w_bb     = tab.get_gamma_hdd(2, 3, 3);
    const double w_cc     = tab.get_gamma_huu(2, 2, 2);
    const double w_tautau = tab.get_gamma_hll(2, 3, 3);
    const double w_WW     = tab.get_gamma_hvv(2, 3);
    const double w_ZZ     = tab.get_gamma_hvv(2, 2);
    const double w_gaga   = tab.get_gamma_hgaga(2);
    const double w_Zga    = tab.get_gamma_hZga(2);
    const double w_gg     = tab.get_gamma_hgg(2);
    const double w_hh     = tab.get_gamma_hhh(2, 1, 1);
    const double w_tot    = tab.get_gammatot_h(2);

    const double br_bb_2hdmc   = (w_tot > 1e-15) ? (w_bb   / w_tot) : qnan();
    const double br_gaga_2hdmc = (w_tot > 1e-15) ? (w_gaga / w_tot) : qnan();
    const double br_Zga_2hdmc  = (w_tot > 1e-15) ? (w_Zga  / w_tot) : qnan();

    out << std::setprecision(17)
        << point_id << ','
        << mh << ','
        << mH << ','
        << mA << ','
        << mHp << ','
        << sba << ','
        << tanb << ','
        << lambda1 << ','
        << lambda6 << ','
        << lambda7 << ','
        << m12sq_analytic << ','
        << br_bb_colab << ','
        << br_gaga_colab << ','
        << br_Zga_colab << ','
        << 1 << ','
        << (validation_available ? 1 : 0) << ','
        << lambda1_input << ','
        << lambda1_recomputed << ','
        << abs_error << ','
        << (warning_flag ? 1 : 0) << ','
        << (pos ? 1 : 0) << ','
        << (uni ? 1 : 0) << ','
        << (pert ? 1 : 0) << ','
        << (triple_ok ? 1 : 0) << ','
        << m12_2_g << ','
        << lam1_g << ','
        << lam2_g << ','
        << lam3_g << ','
        << lam4_g << ','
        << lam5_g << ','
        << w_bb << ','
        << w_cc << ','
        << w_tautau << ','
        << w_WW << ','
        << w_ZZ << ','
        << w_gaga << ','
        << w_Zga << ','
        << w_gg << ','
        << w_hh << ','
        << w_tot << ','
        << br_bb_2hdmc << ','
        << br_gaga_2hdmc << ','
        << br_Zga_2hdmc << '\n';
  }

  std::cout << std::setprecision(17)
            << "Processed rows: " << total_rows << "\n"
            << "set_param_phys_lam1 OK: " << set_param_ok_count << "\n"
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
      << "  Points mode:\n"
      << "    " << argv0 << " points input_points.csv output_validation.csv\n";
}

}  // namespace

int main(int argc, char *argv[]) {
  try {
    if (argc >= 2 && argc <= 4 &&
        std::string(argv[1]) != "random" &&
        std::string(argv[1]) != "points") {
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

    if (argc == 4 && std::string(argv[1]) == "points") {
      return run_points_mode(argv[2], argv[3]);
    }

    // run single point mode
    // ./CalcLam1Scan mh mH mA mHp sin(beta-alpha) lambda_6 lambda_7 lambda_1 tan_beta yukawas_type
    // ./CalcLam1Scan 125 130.7018 300 300 1.0 0.01 0.0 1.7082857205846682 10000 1 worstpoint.lhe
    // 
    if (argc == 11) {

      // cast all variables to float to avoid issues with 2hdmc's internal use of floats for parameters
      float mh = std::stof(argv[1]);
      float mH = std::stof(argv[2]);
      float mA = std::stof(argv[3]);  
      float mHp = std::stof(argv[4]);
      float sba = std::stof(argv[5]);
      float lambda6 = std::stof(argv[6]);
      float lambda7 = std::stof(argv[7]);
      float lambda1_input = std::stof(argv[8]);
      float tan_beta = std::stof(argv[9]);
      int yukawas_type = std::stoi(argv[10]);
      
      return run_single_point_mode(mh, mH, mA, mHp, sba, lambda6, lambda7, lambda1_input, tan_beta, yukawas_type);
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