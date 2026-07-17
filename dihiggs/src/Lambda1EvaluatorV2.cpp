#include "Constraints.h"
#include "DecayTable.h"
#include "ReplaySafeOutput.hpp"
#include "THDM.h"
#include "YukawaType.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr const char* kSchemaVersion = "dihiggs.lambda1.v2";
constexpr const char* kEvaluatorApi = "THDM::set_param_phys_lam1+2HDMC::DecayTable";
constexpr double kHbarCGeVmm = 1.973269804e-13;
constexpr std::size_t kNumericInputs = 9;

double nan() { return std::numeric_limits<double>::quiet_NaN(); }

struct Input {
    std::string point_id;
    std::array<std::string, kNumericInputs> raw{};
    std::array<double, kNumericInputs> value{{nan(), nan(), nan(), nan(), nan(), nan(), nan(), nan(), nan()}};
};

enum InputIndex : std::size_t {
    kMh, kMH, kMA, kMHp, kSinBa, kTanBeta, kLambda1, kLambda6, kLambda7
};

struct Result {
    Input input;
    int construction_ok = 0;
    std::string rejection_stage = "input";
    std::string rejection_reason = "invalid_input";
    std::array<double, 7> lambda_reconstructed{{nan(), nan(), nan(), nan(), nan(), nan(), nan()}};
    double lambda1_abs_residual = nan();
    int lambda1_residual_warning = 0;
    double m12_sq_reconstructed_gev2 = nan();
    double tan_beta_reconstructed = nan();
    int positivity_ok = 0;
    int unitarity_ok = 0;
    int perturbativity_ok = 0;
    int stability_ok = 0;
    int triple_ok = 0;
    int theory_ok = 0;
    std::array<double, 9> widths{{nan(), nan(), nan(), nan(), nan(), nan(), nan(), nan(), nan()}};
    std::array<double, 9> branching_ratios{{nan(), nan(), nan(), nan(), nan(), nan(), nan(), nan(), nan()}};
    double total_width_gev = nan();
    double width_unaccounted_gev = nan();
    double yukawa_type_installed = nan();
    int width_ok = 0;
    double ctau_mm = nan();
};

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(field);
    if (line.empty() || line.back() == ',') fields.emplace_back();
    return fields;
}

double parse_double(const std::string& text, const char* name) {
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(value)) {
        throw std::runtime_error(std::string("invalid_") + name);
    }
    return value;
}

std::uint64_t fnv1a64(const std::string& text) {
    std::uint64_t hash = UINT64_C(14695981039346656037);
    for (const unsigned char byte : text) {
        hash ^= byte;
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

std::string generated_id(const std::string& canonical) {
    std::ostringstream out;
    out << "lambda1_" << std::hex << std::setfill('0') << std::setw(16) << fnv1a64(canonical);
    return out.str();
}

std::string canonical_values(const Input& input) {
    std::ostringstream out;
    out << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const double value : input.value) out << value << ',';
    return out.str();
}

Input parse_input(const std::vector<std::string>& fields) {
    if (fields.size() != kNumericInputs + 1) throw std::runtime_error("wrong_field_count");
    static constexpr const char* names[kNumericInputs] = {
        "mh_gev", "mH_gev", "mA_gev", "mHp_gev", "sin_beta_minus_alpha",
        "tan_beta", "lambda1_target", "lambda6_input", "lambda7_input"
    };
    Input input;
    input.point_id = fields[0];
    for (std::size_t i = 0; i < kNumericInputs; ++i) {
        input.raw[i] = fields[i + 1];
        input.value[i] = parse_double(input.raw[i], names[i]);
    }
    if (input.point_id.empty()) input.point_id = generated_id(canonical_values(input));
    return input;
}

Input partial_input(const std::vector<std::string>& fields, const std::string& raw_line) {
    Input input;
    if (!fields.empty()) input.point_id = fields[0];
    for (std::size_t i = 0; i < kNumericInputs && i + 1 < fields.size(); ++i) input.raw[i] = fields[i + 1];
    if (input.point_id.empty()) input.point_id = generated_id(raw_line);
    return input;
}

std::string construction_reason(const Input& input) {
    if (input.value[kMh] > input.value[kMH]) return "mh_gt_mH";
    if (input.value[kTanBeta] <= 0.0) return "tan_beta_nonpositive";
    if (std::abs(input.value[kSinBa]) > 1.0) return "sin_beta_minus_alpha_out_of_range";
    if (input.value[kMh] < 0.0 || input.value[kMH] < 0.0 ||
        input.value[kMA] < 0.0 || input.value[kMHp] < 0.0) return "negative_mass";
    return "set_param_phys_lam1_returned_false";
}

Result evaluate(const Input& input) {
    Result result;
    result.input = input;
    result.rejection_stage = "construction";
    result.rejection_reason = construction_reason(input);

    THDM model;
    SM sm;
    model.set_SM(sm);
    if (!model.set_param_phys_lam1(
            input.value[kMh], input.value[kMH], input.value[kMA], input.value[kMHp],
            input.value[kSinBa], input.value[kLambda1], input.value[kLambda6],
            input.value[kLambda7], input.value[kTanBeta])) return result;

    result.construction_ok = 1;
    dihiggs::install_yukawa_type(model, 1);
    result.yukawa_type_installed = model.get_yukawas_type();
    if (model.has_param_phys_lam1_validation()) {
        double ignored_input = nan();
        bool warning = false;
        model.get_param_phys_lam1_validation(
            ignored_input, result.lambda_reconstructed[0], result.lambda1_abs_residual, warning);
        result.lambda1_residual_warning = warning ? 1 : 0;
    }

    model.get_param_gen(
        result.lambda_reconstructed[0], result.lambda_reconstructed[1],
        result.lambda_reconstructed[2], result.lambda_reconstructed[3],
        result.lambda_reconstructed[4], result.lambda_reconstructed[5],
        result.lambda_reconstructed[6], result.m12_sq_reconstructed_gev2,
        result.tan_beta_reconstructed);

    Constraints constraints(model);
    result.positivity_ok = constraints.check_positivity() ? 1 : 0;
    result.unitarity_ok = constraints.check_unitarity() ? 1 : 0;
    result.perturbativity_ok = constraints.check_perturbativity() ? 1 : 0;
    result.stability_ok = constraints.check_stability() ? 1 : 0;
    result.triple_ok = result.positivity_ok && result.unitarity_ok && result.perturbativity_ok;
    result.theory_ok = result.triple_ok;

    if (!result.positivity_ok) {
        result.rejection_stage = "positivity";
        result.rejection_reason = "check_positivity_false";
    } else if (!result.unitarity_ok) {
        result.rejection_stage = "unitarity";
        result.rejection_reason = "check_unitarity_false";
    } else if (!result.perturbativity_ok) {
        result.rejection_stage = "perturbativity";
        result.rejection_reason = "check_perturbativity_false";
    } else {
        result.rejection_stage = "accepted";
        result.rejection_reason = "none";
    }

    DecayTable decays(model);
    result.widths = {{
        decays.get_gamma_hdd(2, 3, 3), decays.get_gamma_huu(2, 2, 2),
        decays.get_gamma_hll(2, 3, 3), decays.get_gamma_hvv(2, 3),
        decays.get_gamma_hvv(2, 2), decays.get_gamma_hgaga(2),
        decays.get_gamma_hZga(2), decays.get_gamma_hgg(2),
        decays.get_gamma_hhh(2, 1, 1)
    }};
    result.total_width_gev = decays.get_gammatot_h(2);
    result.width_ok = std::isfinite(result.total_width_gev) && result.total_width_gev > 0.0;
    bool selected_widths_finite = std::isfinite(result.total_width_gev);
    double selected_width_sum = 0.0;
    for (double width : result.widths) {
        selected_widths_finite = selected_widths_finite && std::isfinite(width);
        selected_width_sum += width;
    }
    if (selected_widths_finite) result.width_unaccounted_gev = result.total_width_gev - selected_width_sum;
    if (result.width_ok) {
        for (std::size_t i = 0; i < result.widths.size(); ++i) {
            result.branching_ratios[i] = result.widths[i] / result.total_width_gev;
        }
        result.ctau_mm = kHbarCGeVmm / result.total_width_gev;
    }
    return result;
}

void write_header(std::ostream& out) {
    out << "schema_version,evaluator_commit,evaluator_dirty,evaluator_api,point_id,"
        << "mh_input_gev_raw,mh_input_gev,mH_input_gev_raw,mH_input_gev,"
        << "mA_input_gev_raw,mA_input_gev,mHp_input_gev_raw,mHp_input_gev,"
        << "sin_beta_minus_alpha_input_raw,sin_beta_minus_alpha_input,"
        << "tan_beta_input_raw,tan_beta_input,lambda1_target_raw,lambda1_target,"
        << "lambda6_input_raw,lambda6_input,lambda7_input_raw,lambda7_input,"
        << "construction_ok,yukawa_type_installed,rejection_stage,rejection_reason,"
        << "lambda1_reconstructed,lambda2_reconstructed,lambda3_reconstructed,"
        << "lambda4_reconstructed,lambda5_reconstructed,lambda6_reconstructed,lambda7_reconstructed,"
        << "lambda1_abs_residual,lambda1_residual_warning,m12_sq_reconstructed_gev2,tan_beta_reconstructed,"
        << "positivity_ok,unitarity_ok,perturbativity_ok,stability_ok,triple_ok,theory_ok,"
        << "width_bb_gev,width_cc_gev,width_tautau_gev,width_WW_gev,width_ZZ_gev,"
        << "width_gammagamma_gev,width_Zgamma_gev,width_gg_gev,width_hh_gev,total_width_gev,width_unaccounted_gev,"
        << "br_bb,br_cc,br_tautau,br_WW,br_ZZ,br_gammagamma,br_Zgamma,br_gg,br_hh,"
        << "width_ok,ctau_mm\n";
}

void write_row(std::ostream& out, const Result& r, const std::string& commit, const std::string& dirty) {
    out << kSchemaVersion << ',' << commit << ',' << dirty << ',' << kEvaluatorApi << ',' << r.input.point_id;
    for (std::size_t i = 0; i < kNumericInputs; ++i) out << ',' << r.input.raw[i] << ',' << r.input.value[i];
    out << ',' << r.construction_ok << ',' << r.yukawa_type_installed << ',' << r.rejection_stage << ',' << r.rejection_reason;
    for (const double value : r.lambda_reconstructed) out << ',' << value;
    out << ',' << r.lambda1_abs_residual << ',' << r.lambda1_residual_warning
        << ',' << r.m12_sq_reconstructed_gev2 << ',' << r.tan_beta_reconstructed
        << ',' << r.positivity_ok << ',' << r.unitarity_ok << ',' << r.perturbativity_ok
        << ',' << r.stability_ok << ',' << r.triple_ok << ',' << r.theory_ok;
    for (const double value : r.widths) out << ',' << value;
    out << ',' << r.total_width_gev << ',' << r.width_unaccounted_gev;
    for (const double value : r.branching_ratios) out << ',' << value;
    out << ',' << r.width_ok << ',' << r.ctau_mm << '\n';
}

int run(const std::string& input_path, const std::string& output_path) {
    std::ifstream input(input_path);
    if (!input) throw std::runtime_error("cannot_open_input");
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty_input");
    const std::string expected = "point_id,mh_gev,mH_gev,mA_gev,mHp_gev,sin_beta_minus_alpha,tan_beta,lambda1_target,lambda6_input,lambda7_input";
    if (line != expected) throw std::runtime_error("invalid_header");

    const std::string temporary_path = output_path + ".tmp";
    std::ofstream output(temporary_path);
    if (!output) throw std::runtime_error("cannot_open_output");
    output << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    write_header(output);
    const std::string commit = replay_safe_output::detect_git_commit();
    const std::string dirty = replay_safe_output::detect_git_dirty();

    std::size_t rows = 0;
    while (std::getline(input, line)) {
        ++rows;
        const std::vector<std::string> fields = split_csv(line);
        Result result;
        result.input = partial_input(fields, line);
        try {
            result = evaluate(parse_input(fields));
        } catch (const std::exception& error) {
            if (std::string(error.what()) == "yukawa_type_installation_mismatch") throw;
            result.rejection_reason = error.what();
        }
        write_row(output, result, commit, dirty);
    }
    output.close();
    if (!output) throw std::runtime_error("output_write_failed");
    if (std::rename(temporary_path.c_str(), output_path.c_str()) != 0) {
        throw std::runtime_error("output_rename_failed");
    }
    std::cout << "ROWS_WRITTEN " << rows << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input_v2.csv output_v2.csv\n";
        return 1;
    }
    try {
        return run(argv[1], argv[2]);
    } catch (const std::exception& error) {
        std::cerr << "Lambda1EvaluatorV2: " << error.what() << '\n';
        return 2;
    }
}
