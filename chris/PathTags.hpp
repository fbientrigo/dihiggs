/*
 * PathTags.hpp
 * ============
 *
 * Pure, header-only helpers for building bronze-shard path segments that
 * mirror dihiggs/app/orchestrator/layout.py and io_utils.py's
 * sanitize_for_path / format_float_tag / format_tanbeta_folder conventions,
 * so Stage 1 (CalcLambda1ScanFixings) output paths look familiar to the
 * existing data-lake layout.
 */

#ifndef CHRIS_PATH_TAGS_HPP
#define CHRIS_PATH_TAGS_HPP

#include <cmath>
#include <cstdio>
#include <cctype>
#include <string>

namespace path_tags {

// Mirrors io_utils.sanitize_for_path: '.' -> 'p', '-' -> 'm', then collapse
// any run of characters outside [A-Za-z0-9_=] into a single '_', and strip
// leading/trailing '_'.
inline std::string sanitize_for_path(const std::string& s_in) {
    std::string s = s_in;
    for (char& c : s) {
        if (c == '.') c = 'p';
        else if (c == '-') c = 'm';
    }

    std::string out;
    bool prev_underscore = false;
    for (char c : s) {
        if (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '=') {
            out += c;
            prev_underscore = false;
        } else if (!prev_underscore) {
            out += '_';
            prev_underscore = true;
        }
    }

    size_t start = out.find_first_not_of('_');
    if (start == std::string::npos) return "";
    size_t end = out.find_last_not_of('_');
    return out.substr(start, end - start + 1);
}

// Mirrors io_utils.format_float_tag(x, ndp) + sanitize_for_path.
// e.g. format_float_tag(300.0, 1) -> "300p0"
inline std::string format_float_tag(double x, int ndp) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", ndp, x);
    return sanitize_for_path(std::string(buf));
}

// Scientific-notation variant for log-spaced axes (e.g. lambda6).
// e.g. format_float_tag_sci(1e-5, 6) -> "1p000000em05"
inline std::string format_float_tag_sci(double x, int sig) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*e", sig, x);
    return sanitize_for_path(std::string(buf));
}

// Mirrors layout.format_tanbeta_folder's tag (not the zero-padded folder
// name): integer tan_beta -> plain decimal string; otherwise %.6g, sanitized.
// e.g. format_tanbeta_tag(10.0) -> "10", format_tanbeta_tag(3162.27766) -> "3162p28"
inline std::string format_tanbeta_tag(double tb) {
    if (tb >= 0.0 && std::fabs(tb - std::round(tb)) < 1e-12) {
        long long tb_int = static_cast<long long>(std::llround(tb));
        return std::to_string(tb_int);
    }
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.6g", tb);
    return sanitize_for_path(std::string(buf));
}

}  // namespace path_tags

#endif  // CHRIS_PATH_TAGS_HPP
