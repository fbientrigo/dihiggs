#include "Constraints.h"
#include "DecayTable.h"
#include "ReplaySafeOutput.hpp"
#include "THDM.h"
#include "YukawaType.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr const char* kSchema = "dihiggs.point.v2";
constexpr const char* kProducer = "DihiggsPointV2Evaluator";
constexpr const char* kApi = "THDM::set_param_phys+THDM::get_param_gen+THDM::get_coupling_hhh+2HDMC::DecayTable";
constexpr const char* kStabilityAlias = "2HDMC-patched:check_stability==check_positivity";
constexpr double kHbarCGeVmm = 1.973269804e-13;

double nan() { return std::numeric_limits<double>::quiet_NaN(); }

struct Config {
    std::string campaign_id, run_id, output;
    double mh, mH_min, mH_max, mA, mHp, sin_ba, tan_beta, M2_min, M2_max, lambda6, lambda7;
    int n_mH, n_M2, yukawa_type;
};

struct Result {
    double mh, mH, mA, mHp, sin_ba, tan_beta, beta, M2, m12_sq, lambda6, lambda7;
    int yukawa_type;
    double yukawa_type_installed = nan();
    int construction_ok = 0;
    double numerical_ok = nan();
    std::string rejection_stage = "construction";
    std::string rejection_reason = "set_param_phys_returned_false";
    std::array<double, 7> lambda{{nan(), nan(), nan(), nan(), nan(), nan(), nan()}};
    double tan_beta_reconstructed = nan();
    double m12_sq_reconstructed = nan();
    double M2_reconstructed = nan();
    double positivity_ok = nan(), unitarity_ok = nan(), perturbativity_ok = nan();
    double stability_ok = nan(), triple_ok = nan(), theory_ok = nan();
    double experimental_evaluated = 0.0, experimental_ok = nan();
    double g_hH2H2_GeV = nan();
    std::array<double, 9> widths{{nan(), nan(), nan(), nan(), nan(), nan(), nan(), nan(), nan()}};
    std::array<double, 9> brs{{nan(), nan(), nan(), nan(), nan(), nan(), nan(), nan(), nan()}};
    double total_width = nan(), width_unaccounted = nan(), width_ok = nan(), ctau_mm = nan();
};

double parse_double(const std::string& value, const std::string& name) {
    std::size_t used = 0;
    const double parsed = std::stod(value, &used);
    if (used != value.size() || !std::isfinite(parsed)) throw std::runtime_error("invalid_" + name);
    return parsed;
}

int parse_int(const std::string& value, const std::string& name) {
    std::size_t used = 0;
    const long parsed = std::stol(value, &used);
    if (used != value.size() || parsed < 1 || parsed > std::numeric_limits<int>::max()) {
        throw std::runtime_error("invalid_" + name);
    }
    return static_cast<int>(parsed);
}

int parse_yukawa_type(const std::string& value) {
    std::size_t used = 0;
    const long parsed = std::stol(value, &used);
    if (used != value.size() || parsed < std::numeric_limits<int>::min() ||
        parsed > std::numeric_limits<int>::max()) {
        throw std::runtime_error("invalid_yukawa_type");
    }
    return static_cast<int>(parsed);
}

Config parse_args(int argc, char** argv) {
    std::map<std::string, std::string> args;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (key.rfind("--", 0) != 0 || i + 1 == argc) throw std::runtime_error("invalid_cli");
        if (!args.emplace(key.substr(2), argv[++i]).second) throw std::runtime_error("duplicate_" + key.substr(2));
    }
    static const std::array<const char*, 17> required{{
        "campaign-id", "run-id", "mh", "mH-min", "mH-max", "n-mH", "mA", "mHp",
        "yukawa-type", "sin-ba", "tan-beta", "M2-min", "M2-max", "n-M2",
        "lambda6", "lambda7", "output"
    }};
    for (const char* key : required) if (!args.count(key)) throw std::runtime_error(std::string("missing_") + key);
    if (args.size() != required.size()) throw std::runtime_error("unknown_argument");

    Config c{
        args.at("campaign-id"), args.at("run-id"), args.at("output"),
        parse_double(args.at("mh"), "mh"), parse_double(args.at("mH-min"), "mH-min"),
        parse_double(args.at("mH-max"), "mH-max"), parse_double(args.at("mA"), "mA"),
        parse_double(args.at("mHp"), "mHp"), parse_double(args.at("sin-ba"), "sin-ba"),
        parse_double(args.at("tan-beta"), "tan-beta"), parse_double(args.at("M2-min"), "M2-min"),
        parse_double(args.at("M2-max"), "M2-max"), parse_double(args.at("lambda6"), "lambda6"),
        parse_double(args.at("lambda7"), "lambda7"), parse_int(args.at("n-mH"), "n-mH"),
        parse_int(args.at("n-M2"), "n-M2"), parse_yukawa_type(args.at("yukawa-type"))
    };
    if (c.campaign_id.empty() || c.run_id.empty() || c.output.empty()) throw std::runtime_error("empty_identifier_or_output");
    if (c.campaign_id.find(',') != std::string::npos || c.run_id.find(',') != std::string::npos) throw std::runtime_error("comma_in_identifier");
    if (c.mH_min > c.mH_max || c.M2_min > c.M2_max) throw std::runtime_error("descending_grid");
    if (c.tan_beta <= 0.0 || std::abs(c.sin_ba) > 1.0) throw std::runtime_error("invalid_physics_input");
    if (!dihiggs::supported_yukawa_type(c.yukawa_type)) throw std::runtime_error("unsupported_yukawa_type");
    return c;
}

std::vector<double> grid(double low, double high, int count) {
    if (count == 1) {
        if (low != high) throw std::runtime_error("single_point_grid_requires_equal_bounds");
        return {low};
    }
    std::vector<double> values(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) values[static_cast<std::size_t>(i)] = low + (high - low) * i / (count - 1);
    values.back() = high;
    return values;
}

std::uint64_t fnv1a64(const std::string& text) {
    std::uint64_t hash = UINT64_C(14695981039346656037);
    for (const unsigned char byte : text) { hash ^= byte; hash *= UINT64_C(1099511628211); }
    return hash;
}

std::string point_id(const Result& r) {
    std::ostringstream canonical;
    replay_safe_output::configure_scientific_17(canonical);
    canonical << r.mh << ',' << r.mH << ',' << r.mA << ',' << r.mHp << ',' << r.sin_ba << ','
              << r.tan_beta << ',' << r.M2 << ',' << r.lambda6 << ',' << r.lambda7 << ',' << r.yukawa_type;
    std::ostringstream id;
    id << "point_" << std::hex << std::setfill('0') << std::setw(16) << fnv1a64(canonical.str());
    return id.str();
}

bool all_finite(const Result& r) {
    for (double value : r.lambda) if (!std::isfinite(value)) return false;
    return std::isfinite(r.tan_beta_reconstructed) && std::isfinite(r.m12_sq_reconstructed) &&
           std::isfinite(r.M2_reconstructed);
}

Result evaluate(const Config& c, double mH, double M2) {
    const double beta = std::atan(c.tan_beta);
    Result r{c.mh, mH, c.mA, c.mHp, c.sin_ba, c.tan_beta, beta, M2,
             M2 * std::sin(beta) * std::cos(beta), c.lambda6, c.lambda7, c.yukawa_type};

    THDM model;
    SM sm;
    model.set_SM(sm);
    r.construction_ok = model.set_param_phys(
        r.mh, r.mH, r.mA, r.mHp, r.sin_ba, r.lambda6, r.lambda7, r.m12_sq, r.tan_beta) ? 1 : 0;
    if (!r.construction_ok) {
        if (r.mh > r.mH) r.rejection_reason = "mh_gt_mH";
        return r;
    }
    dihiggs::install_yukawa_type(model, c.yukawa_type);
    r.yukawa_type_installed = model.get_yukawas_type();

    model.get_param_gen(r.lambda[0], r.lambda[1], r.lambda[2], r.lambda[3], r.lambda[4],
                        r.lambda[5], r.lambda[6], r.m12_sq_reconstructed, r.tan_beta_reconstructed);
    const double reconstructed_beta = std::atan(r.tan_beta_reconstructed);
    const double denominator = std::sin(reconstructed_beta) * std::cos(reconstructed_beta);
    r.M2_reconstructed = r.m12_sq_reconstructed / denominator;
    r.numerical_ok = all_finite(r) ? 1.0 : 0.0;
    if (!r.numerical_ok) {
        r.rejection_stage = "numerical";
        r.rejection_reason = "nonfinite_reconstruction";
        return r;
    }

    Constraints constraints(model);
    r.positivity_ok = constraints.check_positivity() ? 1.0 : 0.0;
    r.unitarity_ok = constraints.check_unitarity() ? 1.0 : 0.0;
    r.perturbativity_ok = constraints.check_perturbativity() ? 1.0 : 0.0;
    r.stability_ok = constraints.check_stability() ? 1.0 : 0.0;
    r.triple_ok = (r.positivity_ok && r.unitarity_ok && r.perturbativity_ok) ? 1.0 : 0.0;
    r.theory_ok = r.triple_ok;
    if (!r.positivity_ok) { r.rejection_stage = "positivity"; r.rejection_reason = "check_positivity_false"; }
    else if (!r.unitarity_ok) { r.rejection_stage = "unitarity"; r.rejection_reason = "check_unitarity_false"; }
    else if (!r.perturbativity_ok) { r.rejection_stage = "perturbativity"; r.rejection_reason = "check_perturbativity_false"; }
    else { r.rejection_stage = "accepted"; r.rejection_reason = "none"; }

    std::complex<double> coupling_h1_h2_h2;
    model.get_coupling_hhh(1, 2, 2, coupling_h1_h2_h2);
    r.g_hH2H2_GeV = std::abs(coupling_h1_h2_h2.imag());

    DecayTable decays(model);
    r.widths = {{decays.get_gamma_hdd(2, 3, 3), decays.get_gamma_huu(2, 2, 2),
                 decays.get_gamma_hll(2, 3, 3), decays.get_gamma_hvv(2, 3),
                 decays.get_gamma_hvv(2, 2), decays.get_gamma_hgaga(2),
                 decays.get_gamma_hZga(2), decays.get_gamma_hgg(2), decays.get_gamma_hhh(2, 1, 1)}};
    r.total_width = decays.get_gammatot_h(2);
    r.width_ok = std::isfinite(r.total_width) && r.total_width > 0.0 ? 1.0 : 0.0;
    bool selected_widths_finite = std::isfinite(r.total_width);
    double selected_width_sum = 0.0;
    for (double width : r.widths) {
        selected_widths_finite = selected_widths_finite && std::isfinite(width);
        selected_width_sum += width;
    }
    if (selected_widths_finite) r.width_unaccounted = r.total_width - selected_width_sum;
    if (r.width_ok) {
        for (std::size_t i = 0; i < r.widths.size(); ++i) r.brs[i] = r.widths[i] / r.total_width;
        r.ctau_mm = kHbarCGeVmm / r.total_width;
    }
    return r;
}

void header(std::ostream& out) {
    out << "schema_version,producer,producer_commit,producer_dirty,evaluator_api,campaign_id,run_id,point_id,"
        << "yukawa_type,yukawa_type_installed,mh_input_GeV,mH_input_GeV,mA_input_GeV,mHp_input_GeV,sin_beta_minus_alpha_input,"
        << "tan_beta_input,beta_input_rad,M2_input_GeV2,m12_sq_input_GeV2,lambda6_input,lambda7_input,"
        << "lambda1_reconstructed,lambda2_reconstructed,lambda3_reconstructed,lambda4_reconstructed,"
        << "lambda5_reconstructed,lambda6_reconstructed,lambda7_reconstructed,tan_beta_reconstructed,"
        << "m12_sq_reconstructed_GeV2,M2_reconstructed_GeV2,construction_ok,numerical_ok,rejection_stage,rejection_reason,"
        << "positivity_reported_ok,unitarity_ok,perturbativity_ok,stability_reported_ok,stability_dependency_alias,"
        << "triple_ok_legacy,theory_ok_v1,experimental_evaluated,experimental_ok,g_hH2H2_GeV,"
        << "width_bb_GeV,width_cc_GeV,width_tautau_GeV,width_WW_GeV,width_ZZ_GeV,width_gammagamma_GeV,"
        << "width_Zgamma_GeV,width_gg_GeV,width_hh_GeV,total_width_GeV,width_unaccounted_GeV,br_bb,br_cc,br_tautau,br_WW,br_ZZ,"
        << "br_gammagamma,br_Zgamma,br_gg,br_hh,width_ok,ctau_mm\n";
}

void row(std::ostream& out, const Config& c, const Result& r, const std::string& commit, const std::string& dirty) {
    out << kSchema << ',' << kProducer << ',' << commit << ',' << dirty << ',' << kApi << ','
        << c.campaign_id << ',' << c.run_id << ',' << point_id(r) << ',' << r.yukawa_type << ',' << r.yukawa_type_installed << ','
        << r.mh << ',' << r.mH << ',' << r.mA << ',' << r.mHp << ',' << r.sin_ba << ',' << r.tan_beta << ','
        << r.beta << ',' << r.M2 << ',' << r.m12_sq << ',' << r.lambda6 << ',' << r.lambda7;
    for (double value : r.lambda) out << ',' << value;
    out << ',' << r.tan_beta_reconstructed << ',' << r.m12_sq_reconstructed << ',' << r.M2_reconstructed
        << ',' << r.construction_ok << ',' << r.numerical_ok << ',' << r.rejection_stage << ',' << r.rejection_reason
        << ',' << r.positivity_ok << ',' << r.unitarity_ok << ',' << r.perturbativity_ok << ',' << r.stability_ok
        << ',' << kStabilityAlias << ',' << r.triple_ok << ',' << r.theory_ok << ',' << r.experimental_evaluated
        << ',' << r.experimental_ok << ',' << r.g_hH2H2_GeV;
    for (double value : r.widths) out << ',' << value;
    out << ',' << r.total_width << ',' << r.width_unaccounted;
    for (double value : r.brs) out << ',' << value;
    out << ',' << r.width_ok << ',' << r.ctau_mm << '\n';
}

int run(const Config& c) {
    const std::string temporary = c.output + ".tmp";
    std::ofstream out(temporary);
    if (!out) throw std::runtime_error("cannot_open_output");
    replay_safe_output::configure_scientific_17(out);
    header(out);
    const std::string commit = replay_safe_output::detect_git_commit();
    const std::string dirty = replay_safe_output::detect_git_dirty();
    std::size_t rows = 0, accepted = 0;
    for (double mH : grid(c.mH_min, c.mH_max, c.n_mH)) {
        for (double M2 : grid(c.M2_min, c.M2_max, c.n_M2)) {
            const Result result = evaluate(c, mH, M2);
            row(out, c, result, commit, dirty);
            ++rows;
            accepted += result.theory_ok == 1.0;
        }
    }
    out.close();
    if (!out) throw std::runtime_error("output_write_failed");
    if (std::rename(temporary.c_str(), c.output.c_str()) != 0) throw std::runtime_error("output_rename_failed");
    std::cout << "ROWS_WRITTEN " << rows << "\nTotal Attempts: " << rows << "\nTRIPLE_OK_POINTS " << accepted << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try { return run(parse_args(argc, argv)); }
    catch (const std::exception& error) {
        std::cerr << "DihiggsPointV2Evaluator: " << error.what() << '\n';
        return 2;
    }
}