/**
 * @file GenScanWithFixings.cpp
 * @brief "Stage 2" validator: reads a chris/CalcLambda1ScanFixings bronze
 *        shard CSV, calibrates each row's reconstructed generic-basis
 *        lambdas against 2HDMC (set_param_gen -> get_param_phys), and
 *        writes a "validated points" CSV in the legacy 29-column schema
 *        plus calibration/stability/chris-cross-check diagnostics.
 *
 * Usage:
 *   GenScanWithFixings --bronze-csv=<shard.csv> --output-csv=<validated.csv>
 *       [--calibration-n=50] [--calibration-frac=0.10] [--rng-seed=0]
 */

#include "GenScanPointEvaluator.hpp"
#include "ReplaySafeOutput.hpp"

#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;
using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;

namespace {

// Parses "--key=value" style args. Returns true and sets `value` on match.
bool parse_flag(const std::string& arg, const std::string& key, std::string& value) {
    const std::string prefix = "--" + key + "=";
    if (arg.compare(0, prefix.size(), prefix) == 0) {
        value = arg.substr(prefix.size());
        return true;
    }
    return false;
}

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    std::istringstream ss(line);
    while (std::getline(ss, field, ',')) {
        fields.push_back(field);
    }
    return fields;
}

// Header-name-keyed reader: robust to Stage 1 column-order changes.
std::vector<gen_scan::BronzeRow> read_bronze_csv(const std::string& path) {
    std::vector<gen_scan::BronzeRow> rows;

    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "ERROR: cannot open bronze CSV: " << path << "\n";
        return rows;
    }

    std::string header_line;
    if (!std::getline(in, header_line)) {
        std::cerr << "ERROR: bronze CSV is empty: " << path << "\n";
        return rows;
    }

    const std::vector<std::string> header = split_csv_line(header_line);
    std::unordered_map<std::string, size_t> col_index;
    for (size_t i = 0; i < header.size(); ++i) {
        col_index[header[i]] = i;
    }

    auto get_double = [&](const std::vector<std::string>& fields, const std::string& name) -> double {
        auto it = col_index.find(name);
        if (it == col_index.end() || it->second >= fields.size()) return 0.0;
        return std::stod(fields[it->second]);
    };
    auto get_int = [&](const std::vector<std::string>& fields, const std::string& name) -> int {
        auto it = col_index.find(name);
        if (it == col_index.end() || it->second >= fields.size()) return 0;
        return static_cast<int>(std::stod(fields[it->second]));
    };

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = split_csv_line(line);

        gen_scan::BronzeRow row;
        row.tan_beta = get_double(fields, "tan_beta");
        row.m_A = get_double(fields, "m_A");
        row.lambda6 = get_double(fields, "lambda6");
        row.lambda7 = get_double(fields, "lambda7");
        row.lambda1_input = get_double(fields, "lambda1_input");
        row.sin_ba = get_double(fields, "sin_ba");
        row.mh_input = get_double(fields, "mh_input");
        row.yukawa_type = get_int(fields, "yukawa_type");
        row.m_H_input = get_double(fields, "m_H_input");
        row.m12sq_recon = get_double(fields, "m12sq_recon");
        row.lambda1_recon = get_double(fields, "lambda1_recon");
        row.lambda2_recon = get_double(fields, "lambda2_recon");
        row.lambda3_recon = get_double(fields, "lambda3_recon");
        row.lambda4_recon = get_double(fields, "lambda4_recon");
        row.lambda5_recon = get_double(fields, "lambda5_recon");
        row.chris_width_bb = get_double(fields, "chris_width_bb");
        row.chris_width_tautau = get_double(fields, "chris_width_tautau");
        row.chris_width_gg = get_double(fields, "chris_width_gg");
        row.chris_width_gaga = get_double(fields, "chris_width_gaga");
        row.chris_width_Zga = get_double(fields, "chris_width_Zga");
        row.chris_total_width = get_double(fields, "chris_total_width");
        row.chris_ctau_m = get_double(fields, "chris_ctau_m");

        rows.push_back(row);
    }

    return rows;
}

void write_result_row(std::ostream& out, const gen_scan::GenFixingsPointResult& r) {
    out
        << r.m_phi << "," << r.mA << "," << r.alpha << "," << r.beta << ","
        << r.lambda6 << "," << r.lambda7 << "," << r.m12 << "," << r.sin_ba << "," << r.tan_beta << ","
        << r.positivity_ok << "," << r.unitarity_ok << "," << r.perturbativity_ok << ","
        << r.width_bb << "," << r.width_tautau << "," << r.width_WW << "," << r.width_ZZ << ","
        << r.width_gaga << "," << r.width_Zga << "," << r.width_gg << "," << r.width_hh << ","
        << r.total_width << "," << r.br_gaga << ","
        << r.lam1 << "," << r.computed_lam1 << "," << r.lam2 << "," << r.computed_lam2 << ","
        << r.lam3 << "," << r.lam4 << "," << r.lam5 << ","
        << r.mA_target << "," << r.mH_target << "," << r.mh_calibrated << "," << r.mHp_calibrated << ","
        << r.sba_calibrated << "," << r.calibration_score << "," << r.calibration_n_used << ","
        << r.stability_ok << ","
        << r.chris_width_bb << "," << r.chris_width_tautau << "," << r.chris_width_gg << ","
        << r.chris_width_gaga << "," << r.chris_width_Zga << "," << r.chris_ctau_mm << ","
        << r.delta_width_bb << "," << r.ratio_width_bb << ","
        << r.delta_width_tautau << "," << r.ratio_width_tautau << ","
        << r.delta_width_gg << "," << r.ratio_width_gg << ","
        << r.delta_width_gaga << "," << r.ratio_width_gaga << ","
        << r.delta_width_Zga << "," << r.ratio_width_Zga << ","
        << r.delta_ctau_mm << "," << r.ratio_ctau_mm << ",";
    replay_safe_output::append_metadata_csv(out, r.meta);
    out << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::string bronze_csv;
    std::string output_csv;
    gen_scan::CalibrationConfig cfg;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        std::string value;
        if (parse_flag(arg, "bronze-csv", value)) {
            bronze_csv = value;
        } else if (parse_flag(arg, "output-csv", value)) {
            output_csv = value;
        } else if (parse_flag(arg, "calibration-n", value)) {
            cfg.n_samples = std::stoi(value);
        } else if (parse_flag(arg, "calibration-frac", value)) {
            cfg.variation_fraction = std::stod(value);
        } else if (parse_flag(arg, "rng-seed", value)) {
            cfg.rng_seed = static_cast<unsigned int>(std::stoul(value));
        } else {
            std::cerr << "WARNING: unrecognized argument: " << arg << "\n";
        }
    }

    if (bronze_csv.empty() || output_csv.empty()) {
        std::cerr << "Usage: " << argv[0]
                  << " --bronze-csv=<shard.csv> --output-csv=<validated.csv>"
                     " [--calibration-n=50] [--calibration-frac=0.10] [--rng-seed=0]\n";
        return 1;
    }

    const std::vector<gen_scan::BronzeRow> rows = read_bronze_csv(bronze_csv);
    const long long total_rows = static_cast<long long>(rows.size());
    std::cout << "Bronze rows read: " << total_rows << std::endl;

    const fs::path out_path(output_csv);
    if (out_path.has_parent_path()) {
        std::error_code ec;
        fs::create_directories(out_path.parent_path(), ec);
        if (ec) {
            std::cerr << "ERROR: cannot create output directory " << out_path.parent_path()
                      << ": " << ec.message() << "\n";
            return 1;
        }
    }

    std::ofstream results(output_csv);
    if (!results.is_open()) {
        std::cerr << "ERROR: cannot open output CSV: " << output_csv << "\n";
        return 1;
    }

    const std::vector<std::string> columns = gen_scan::output_csv_columns();
    for (size_t c = 0; c < columns.size(); ++c) {
        results << columns[c] << (c + 1 < columns.size() ? "," : "\n");
    }

    std::atomic<long long> global_attempts(0);
    std::atomic<long long> global_triple_ok(0);

    auto start_time = Clock::now();
    auto last_report_time = start_time;
    const double REPORT_INTERVAL_SEC = 2.5;
    const int BATCH_SIZE = 200;

    #pragma omp parallel
    {
        long long local_attempts = 0;
        long long local_triple_ok = 0;
        long long local_rows = 0;

        std::ostringstream local_buffer;
        replay_safe_output::configure_scientific_17(local_buffer);

        #pragma omp for schedule(dynamic)
        for (long long idx = 0; idx < total_rows; ++idx) {
            const unsigned int thread_seed = cfg.rng_seed + static_cast<unsigned int>(idx);
            const gen_scan::GenFixingsPointResult result =
                gen_scan::evaluate_gen_fixings_point(rows[idx], cfg, thread_seed);

            write_result_row(local_buffer, result);

            ++local_attempts;
            ++local_rows;
            if (result.positivity_ok && result.unitarity_ok && result.perturbativity_ok) {
                ++local_triple_ok;
            }

            if (local_rows >= BATCH_SIZE) {
                #pragma omp critical
                {
                    results << local_buffer.str();
                    global_attempts += local_attempts;
                    global_triple_ok += local_triple_ok;

                    auto now = Clock::now();
                    Duration diff = now - last_report_time;
                    if (diff.count() > REPORT_INTERVAL_SEC) {
                        const long long current_attempts = global_attempts.load();
                        const double progress = (total_rows > 0)
                            ? (100.0 * static_cast<double>(current_attempts) / static_cast<double>(total_rows))
                            : 0.0;
                        std::cerr << "\r[ " << std::fixed << std::setprecision(1) << progress << "% ] "
                                  << "Rows: " << current_attempts << "/" << total_rows
                                  << " | TripleOK: " << global_triple_ok.load()
                                  << "   " << std::flush;
                        last_report_time = now;
                    }
                }

                local_buffer.str("");
                local_buffer.clear();
                replay_safe_output::configure_scientific_17(local_buffer);
                local_attempts = 0;
                local_triple_ok = 0;
                local_rows = 0;
            }
        }

        #pragma omp critical
        {
            if (local_rows > 0) {
                results << local_buffer.str();
            }
            if (local_attempts > 0) {
                global_attempts += local_attempts;
            }
            if (local_triple_ok > 0) {
                global_triple_ok += local_triple_ok;
            }
        }
    }

    results.close();

    const double total_time = Duration(Clock::now() - start_time).count();
    const long long final_attempts = global_attempts.load();
    const long long final_triple_ok = global_triple_ok.load();

    std::cout << "\n\n--- Stage 2 Validation Finished ---" << std::endl;
    std::cout << "Total Attempts: " << final_attempts << std::endl;
    std::cout << "Total Time: " << total_time << " s" << std::endl;
    std::cout << "TRIPLE_OK_POINTS " << final_triple_ok << std::endl;

    return 0;
}
