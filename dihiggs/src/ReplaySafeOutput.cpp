#include "ReplaySafeOutput.hpp"

#include <cstdlib>
#include <iomanip>
#include <ostream>
#include <sstream>

namespace replay_safe_output {

const char* default_replay_semantics_version() {
    return "pswf_m12_used_v1";
}

const char* default_higgs_state() {
    return "h2/DECAY35";
}

const char* default_model_api_path() {
    return "THDM::set_param_phys_lam1->THDM::get_param_gen->THDM::set_param_phys";
}

const char* default_calc_engine() {
    return "2HDMC::DecayTable";
}

std::string detect_git_commit() {
    const char* env = std::getenv("DIHIGGS_GIT_COMMIT");
    return (env && *env) ? std::string(env) : std::string("unknown");
}

std::string detect_git_dirty() {
    const char* env = std::getenv("DIHIGGS_GIT_DIRTY");
    return (env && *env) ? std::string(env) : std::string("unknown");
}

Metadata make_metadata(
    double m12_2_used,
    double m12_2_gen_after_set,
    int yukawa_type
) {
    Metadata md;
    md.m12_2_used = m12_2_used;
    md.m12_2_gen_after_set = m12_2_gen_after_set;
    md.delta_m12_2_gen_minus_used = m12_2_gen_after_set - m12_2_used;
    md.replay_semantics_version = default_replay_semantics_version();
    md.yukawa_type = yukawa_type;
    md.higgs_state = default_higgs_state();
    md.model_api_path = default_model_api_path();
    md.calc_engine = default_calc_engine();
    md.git_commit = detect_git_commit();
    md.git_dirty = detect_git_dirty();
    return md;
}

void configure_scientific_17(std::ostream& os) {
    os.setf(std::ios::scientific);
    os << std::setprecision(kDoublePrecisionDigits);
}

std::string format_double_scientific_17(double value) {
    std::ostringstream oss;
    configure_scientific_17(oss);
    oss << value;
    return oss.str();
}

std::vector<std::string> metadata_column_names() {
    return {
        "m12_2_used",
        "m12_2_gen_after_set",
        "delta_m12_2_gen_minus_used",
        "replay_semantics_version",
        "yukawa_type",
        "higgs_state",
        "model_api_path",
        "calc_engine",
        "git_commit",
        "git_dirty"
    };
}

void append_metadata_csv(std::ostream& os, const Metadata& md) {
    os
        << format_double_scientific_17(md.m12_2_used) << ","
        << format_double_scientific_17(md.m12_2_gen_after_set) << ","
        << format_double_scientific_17(md.delta_m12_2_gen_minus_used) << ","
        << md.replay_semantics_version << ","
        << md.yukawa_type << ","
        << md.higgs_state << ","
        << md.model_api_path << ","
        << md.calc_engine << ","
        << md.git_commit << ","
        << md.git_dirty;
}

}  // namespace replay_safe_output
