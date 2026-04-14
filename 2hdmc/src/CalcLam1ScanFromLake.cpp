#include "THDM.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct LakeRow {
  std::size_t row_index;
  double m_phi;
  double mA;
  double lambda6;
  double lambda7;
  double m12;
  double sin_ba;
  double tan_beta;
  bool positivity_ok;
  bool unitarity_ok;
  bool perturbativity_ok;
};

struct OutputRow {
  std::size_t sample_index;
  std::size_t row_index;
  double m_phi;
  double mA;
  double mHp;
  double sin_ba;
  double lambda6;
  double lambda7;
  double m12_input;
  double tan_beta;
  double lambda1_input;
  double lambda1_recomputed;
  double abs_error;
  double rel_error;
  double scaled_error;
  bool exact_equal;
  std::uint64_t ulp_diff;
  bool warning_internal;
  bool warning_abs;
  bool warning_rel;
  bool warning_scaled;
  bool baseline_ok;
  bool probe_ok;
  std::string mode;
};

using HeaderMap = std::unordered_map<std::string, std::size_t>;

std::vector<std::string> split_csv_line(const std::string &line) {
  std::vector<std::string> fields;
  std::stringstream stream(line);
  std::string field;
  while (std::getline(stream, field, ',')) {
    fields.push_back(field);
  }
  if (!line.empty() && line.back() == ',') {
    fields.push_back("");
  }
  return fields;
}

HeaderMap parse_header_map(const std::string &header_line) {
  HeaderMap header_map;
  const std::vector<std::string> fields = split_csv_line(header_line);
  for (std::size_t index = 0; index < fields.size(); ++index) {
    header_map[fields[index]] = index;
  }
  return header_map;
}

bool parse_double_field(const std::vector<std::string> &fields,
                        const HeaderMap &header_map,
                        const std::string &name,
                        double &value) {
  const auto it = header_map.find(name);
  if (it == header_map.end() || it->second >= fields.size()) {
    return false;
  }

  char *end_ptr = nullptr;
  value = std::strtod(fields[it->second].c_str(), &end_ptr);
  return end_ptr != fields[it->second].c_str() && *end_ptr == '\0' && std::isfinite(value);
}

bool parse_flag_field(const std::vector<std::string> &fields,
                      const HeaderMap &header_map,
                      const std::string &name,
                      bool &value) {
  double numeric_value = 0.0;
  if (!parse_double_field(fields, header_map, name, numeric_value)) {
    return false;
  }
  value = std::abs(numeric_value - 1.0) < 1e-12;
  return true;
}

bool row_passes_mode(const LakeRow &row, const std::string &mode) {
  if (mode == "all") {
    return true;
  }

  const bool triple_ok = row.positivity_ok && row.unitarity_ok && row.perturbativity_ok;
  if (mode == "triple_ok") {
    return triple_ok;
  }
  if (mode == "triple_ok_sba1") {
    return triple_ok && std::abs(row.sin_ba - 1.0) < 1e-12;
  }
  if (mode == "triple_ok_align" || mode == "alignment") {
    return triple_ok && row.sin_ba >= 0.995;
  }
  return false;
}

std::string normalize_mode(const std::string &mode) {
  if (mode == "alignment") {
    return "triple_ok_align";
  }
  return mode;
}

bool valid_mode(const std::string &mode) {
  return mode == "triple_ok" || mode == "triple_ok_sba1" || mode == "triple_ok_align" ||
         mode == "alignment" || mode == "all";
}

std::uint64_t monotonic_double_bits(double value) {
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(double));

  constexpr std::uint64_t sign_mask = 0x8000000000000000ULL;
  return (bits & sign_mask) ? ~bits : (bits | sign_mask);
}

std::uint64_t ulp_distance(double a, double b) {
  if (!std::isfinite(a) || !std::isfinite(b)) {
    return std::numeric_limits<std::uint64_t>::max();
  }
  const std::uint64_t ma = monotonic_double_bits(a);
  const std::uint64_t mb = monotonic_double_bits(b);
  return (ma >= mb) ? (ma - mb) : (mb - ma);
}

void write_header(std::ofstream &out) {
  out << "sample_index,row_index,m_phi,mA,mHp,sin_ba,lambda6,lambda7,m12_input,tan_beta,"
         "lambda1_input,lambda1_recomputed,abs_error,rel_error,scaled_error,exact_equal,ulp_diff,"
         "warning_internal,warning_abs,warning_rel,warning_scaled,baseline_ok,probe_ok,mode\n";
}

void write_output_row(std::ofstream &out, const OutputRow &row) {
  out << row.sample_index << ','
      << row.row_index << ','
      << std::scientific << std::setprecision(17)
      << row.m_phi << ','
      << row.mA << ','
      << row.mHp << ','
      << row.sin_ba << ','
      << row.lambda6 << ','
      << row.lambda7 << ','
      << row.m12_input << ','
      << row.tan_beta << ','
      << row.lambda1_input << ','
      << row.lambda1_recomputed << ','
      << row.abs_error << ','
      << row.rel_error << ','
      << row.scaled_error << ','
      << (row.exact_equal ? 1 : 0) << ','
      << row.ulp_diff << ','
      << (row.warning_internal ? 1 : 0) << ','
      << (row.warning_abs ? 1 : 0) << ','
      << (row.warning_rel ? 1 : 0) << ','
      << (row.warning_scaled ? 1 : 0) << ','
      << (row.baseline_ok ? 1 : 0) << ','
      << (row.probe_ok ? 1 : 0) << ','
      << row.mode << '\n';
}

void print_usage() {
  std::cerr << "Usage: ./CalcLam1ScanFromLake input_csv output_csv [n_samples] [seed] [mode]\n"
            << "Modes: triple_ok, triple_ok_sba1, triple_ok_align, alignment, all\n";
}

}

int main(int argc, char *argv[]) {
  if (argc < 3 || argc > 6) {
    print_usage();
    return 1;
  }

  const std::string input_csv = argv[1];
  const std::string output_csv = argv[2];
  const std::size_t requested_samples = (argc >= 4) ? static_cast<std::size_t>(std::strtoull(argv[3], nullptr, 10)) : 5000;
  const unsigned long long seed = (argc >= 5) ? std::strtoull(argv[4], nullptr, 10) : 123456ULL;
  const std::string raw_mode = (argc >= 6) ? argv[5] : "triple_ok";
  const std::string mode = normalize_mode(raw_mode);

  if (!valid_mode(raw_mode)) {
    std::cerr << "Unknown mode: " << raw_mode << "\n";
    print_usage();
    return 2;
  }

  std::ifstream input(input_csv.c_str());
  if (!input) {
    std::cerr << "Failed to open input CSV: " << input_csv << "\n";
    return 3;
  }

  std::string header_line;
  if (!std::getline(input, header_line)) {
    std::cerr << "Input CSV is empty: " << input_csv << "\n";
    return 4;
  }

  const HeaderMap header_map = parse_header_map(header_line);
  const std::vector<std::string> required_columns = {
      "m_phi", "mA", "lambda6", "lambda7", "m12", "sin_ba", "tan_beta",
      "positivity_ok", "unitarity_ok", "perturbativity_ok"};
  for (const std::string &column : required_columns) {
    if (header_map.find(column) == header_map.end()) {
      std::cerr << "Missing required column: " << column << "\n";
      return 5;
    }
  }

  std::vector<LakeRow> eligible_rows;
  std::string line;
  std::size_t total_rows_parsed = 0;
  std::size_t parseable_rows = 0;

  while (std::getline(input, line)) {
    if (line.empty()) {
      continue;
    }
    ++total_rows_parsed;
    const std::vector<std::string> fields = split_csv_line(line);

    LakeRow row{};
    row.row_index = total_rows_parsed - 1;
    if (!parse_double_field(fields, header_map, "m_phi", row.m_phi) ||
        !parse_double_field(fields, header_map, "mA", row.mA) ||
        !parse_double_field(fields, header_map, "lambda6", row.lambda6) ||
        !parse_double_field(fields, header_map, "lambda7", row.lambda7) ||
        !parse_double_field(fields, header_map, "m12", row.m12) ||
        !parse_double_field(fields, header_map, "sin_ba", row.sin_ba) ||
        !parse_double_field(fields, header_map, "tan_beta", row.tan_beta) ||
        !parse_flag_field(fields, header_map, "positivity_ok", row.positivity_ok) ||
        !parse_flag_field(fields, header_map, "unitarity_ok", row.unitarity_ok) ||
        !parse_flag_field(fields, header_map, "perturbativity_ok", row.perturbativity_ok)) {
      continue;
    }

    ++parseable_rows;
    if (row_passes_mode(row, mode)) {
      eligible_rows.push_back(row);
    }
  }

  if (eligible_rows.empty()) {
    std::cerr << "No eligible rows found for mode " << mode << " in " << input_csv << "\n";
    return 6;
  }

  const std::size_t eligible_row_count = eligible_rows.size();

  std::size_t actual_samples = requested_samples;
  if (actual_samples > eligible_row_count) {
    std::cerr << "Requested " << actual_samples << " samples but only " << eligible_rows.size()
              << " eligible rows are available. Truncating to eligible row count.\n";
    actual_samples = eligible_row_count;
  }

  std::mt19937_64 rng(seed);
  std::shuffle(eligible_rows.begin(), eligible_rows.end(), rng);
  eligible_rows.resize(actual_samples);

  std::ofstream output(output_csv.c_str());
  if (!output) {
    std::cerr << "Failed to open output CSV: " << output_csv << "\n";
    return 7;
  }
  write_header(output);

  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::size_t baseline_failures = 0;
  std::size_t probe_failures = 0;
  std::size_t warning_abs_count = 0;
  std::size_t warning_rel_count = 0;
  std::size_t warning_scaled_count = 0;
  std::size_t exact_equal_count = 0;
  std::size_t compared_rows = 0;
  double max_abs_error = 0.0;
  double max_rel_error = 0.0;
  double max_scaled_error = 0.0;
  std::uint64_t max_ulp_diff = 0;

  for (std::size_t sample_index = 0; sample_index < eligible_rows.size(); ++sample_index) {
    const LakeRow &row = eligible_rows[sample_index];
    OutputRow out_row{sample_index, row.row_index, row.m_phi, row.mA, row.mA, row.sin_ba,
                      row.lambda6, row.lambda7, row.m12, row.tan_beta,
                      nan, nan, nan, nan, nan,
                      false, 0,
                      false, false, false, false, false, false, mode};

    THDM baseline;
    if (!baseline.set_param_phys(125.0, row.m_phi, row.mA, row.mA,
                                 row.sin_ba, row.lambda6, row.lambda7,
                                 row.m12, row.tan_beta)) {
      ++baseline_failures;
      write_output_row(output, out_row);
      continue;
    }
    baseline.set_yukawas_type(1);
    out_row.baseline_ok = true;

    double lambda1_truth = 0.0;
    double lambda2_truth = 0.0;
    double lambda3_truth = 0.0;
    double lambda4_truth = 0.0;
    double lambda5_truth = 0.0;
    double lambda6_truth = 0.0;
    double lambda7_truth = 0.0;
    double m12_truth = 0.0;
    double tan_beta_truth = 0.0;
    baseline.get_param_gen(lambda1_truth, lambda2_truth, lambda3_truth, lambda4_truth,
                           lambda5_truth, lambda6_truth, lambda7_truth, m12_truth,
                           tan_beta_truth);

    THDM probe;
    if (!probe.set_param_phys_lam1(125.0, row.m_phi, row.mA, row.mA,
                                   row.sin_ba, lambda1_truth,
                                   row.lambda6, row.lambda7, row.tan_beta)) {
      ++probe_failures;
      out_row.lambda1_input = lambda1_truth;
      write_output_row(output, out_row);
      continue;
    }
    probe.set_yukawas_type(1);
    out_row.probe_ok = true;

    if (!probe.has_param_phys_lam1_validation()) {
      ++probe_failures;
      out_row.lambda1_input = lambda1_truth;
      out_row.probe_ok = false;
      write_output_row(output, out_row);
      continue;
    }

    probe.get_param_phys_lam1_validation(out_row.lambda1_input,
                                         out_row.lambda1_recomputed,
                                         out_row.abs_error,
                                         out_row.warning_internal);
    ++compared_rows;
    out_row.exact_equal = (out_row.lambda1_input == out_row.lambda1_recomputed);
    out_row.ulp_diff = ulp_distance(out_row.lambda1_input, out_row.lambda1_recomputed);
    if (out_row.exact_equal) {
      ++exact_equal_count;
    }
    out_row.rel_error = out_row.abs_error / std::max(std::abs(out_row.lambda1_input), 1.0);
    out_row.scaled_error = out_row.abs_error /
                           std::max(std::max(std::abs(out_row.lambda1_input), std::abs(out_row.lambda1_recomputed)), 1.0);
    out_row.warning_abs = out_row.abs_error > 1e-12;
    out_row.warning_rel = out_row.rel_error > 1e-12;
    out_row.warning_scaled = out_row.scaled_error > 1e-12;

    if (out_row.warning_abs) {
      ++warning_abs_count;
    }
    if (out_row.warning_rel) {
      ++warning_rel_count;
    }
    if (out_row.warning_scaled) {
      ++warning_scaled_count;
    }
    max_abs_error = std::max(max_abs_error, out_row.abs_error);
    max_rel_error = std::max(max_rel_error, out_row.rel_error);
    max_scaled_error = std::max(max_scaled_error, out_row.scaled_error);
    max_ulp_diff = std::max(max_ulp_diff, out_row.ulp_diff);

    write_output_row(output, out_row);
  }

  const double exact_equal_fraction =
      (compared_rows > 0) ? static_cast<double>(exact_equal_count) / static_cast<double>(compared_rows) : 0.0;

  std::cout << std::scientific << std::setprecision(17)
            << "Input rows parsed: " << total_rows_parsed << "\n"
            << "Parseable rows: " << parseable_rows << "\n"
            << "Eligible rows after filtering: " << eligible_row_count << "\n"
            << "Requested samples: " << requested_samples << "\n"
            << "Actual samples written: " << actual_samples << "\n"
            << "Baseline failures: " << baseline_failures << "\n"
            << "Probe failures: " << probe_failures << "\n"
            << "Max absolute error: " << max_abs_error << "\n"
            << "Max relative error: " << max_rel_error << "\n"
            << "Max scaled error: " << max_scaled_error << "\n"
            << "Max ULP diff: " << max_ulp_diff << "\n"
            << "Exact-equal rows: " << exact_equal_count << "\n"
            << "Exact-equal fraction: " << exact_equal_fraction << "\n"
            << "Warning rate abs: " << (actual_samples > 0 ? static_cast<double>(warning_abs_count) / static_cast<double>(actual_samples) : 0.0) << "\n"
            << "Warning rate rel: " << (actual_samples > 0 ? static_cast<double>(warning_rel_count) / static_cast<double>(actual_samples) : 0.0) << "\n"
            << "Warning rate scaled: " << (actual_samples > 0 ? static_cast<double>(warning_scaled_count) / static_cast<double>(actual_samples) : 0.0) << "\n"
            << "Selected mode: " << mode << "\n"
            << "Seed: " << seed << "\n"
            << "Output CSV: " << output_csv << "\n";

  if (max_ulp_diff == 0) {
    std::cout << "ULP interpretation: all compared lambda1 values were bitwise identical in this sampled set.\n";
  } else {
    std::cout << "ULP interpretation: nonzero ULP differences were observed in this sampled set.\n";
  }

  return 0;
}
