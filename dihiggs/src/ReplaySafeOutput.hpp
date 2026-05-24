#pragma once

#include <iosfwd>
#include <string>
#include <vector>

namespace replay_safe_output {

constexpr int kDoublePrecisionDigits = 17;

struct Metadata {
    double m12_2_used = 0.0;
    double m12_2_gen_after_set = 0.0;
    double delta_m12_2_gen_minus_used = 0.0;

    std::string replay_semantics_version;
    int yukawa_type = 1;
    std::string higgs_state;
    std::string model_api_path;
    std::string calc_engine;

    std::string git_commit;
    std::string git_dirty;
};

const char* default_replay_semantics_version();
const char* default_higgs_state();
const char* default_model_api_path();
const char* default_calc_engine();

std::string detect_git_commit();
std::string detect_git_dirty();

Metadata make_metadata(
    double m12_2_used,
    double m12_2_gen_after_set,
    int yukawa_type
);

void configure_scientific_17(std::ostream& os);
std::string format_double_scientific_17(double value);

std::vector<std::string> metadata_column_names();
void append_metadata_csv(std::ostream& os, const Metadata& metadata);

}  // namespace replay_safe_output
