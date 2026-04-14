#include "THDM.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>

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

void write_header(std::ofstream &out) {
  out << "sample_index,attempt_index,mh,mH,mA,mHp,sba,lambda6,lambda7,m12_2,tan_beta,"
         "lambda1_input,lambda1_recomputed,abs_error,warning_flag\n";
}

void write_sample(std::ofstream &out, const Sample &sample) {
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

}

int main(int argc, char *argv[]) {
  if (argc < 2 || argc > 4) {
    std::cout << "Usage: ./CalcLam1Scan output_csv [n_samples] [seed]\n";
    return 1;
  }

  const std::string output_csv = argv[1];
  const std::size_t n_samples = (argc >= 3) ? static_cast<std::size_t>(std::strtoull(argv[2], nullptr, 10)) : 5000;
  const unsigned long long seed = (argc >= 4) ? std::strtoull(argv[3], nullptr, 10) : 123456ULL;
  const std::size_t max_attempts = std::max<std::size_t>(n_samples * 200, 1000);

  std::ofstream out(output_csv.c_str());
  if (!out) {
    std::cerr << "Failed to open output file: " << output_csv << "\n";
    return 2;
  }

  write_header(out);

  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> mh_dist(120.0, 130.0);
  std::uniform_real_distribution<double> mH_dist(200.0, 600.0);
  std::uniform_real_distribution<double> mA_dist(200.0, 700.0);
  std::uniform_real_distribution<double> mHp_dist(200.0, 700.0);
  std::uniform_real_distribution<double> sba_dist(-0.9999, 0.9999);
  std::uniform_real_distribution<double> lambda6_dist(0, 1.0);
  std::uniform_real_distribution<double> lambda7_dist(0, 1.0);
  std::uniform_real_distribution<double> m12_dist(1.0, 500.0);
  std::uniform_real_distribution<double> tan_beta_dist(100, 50000.0);

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
    write_sample(out, sample);

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
