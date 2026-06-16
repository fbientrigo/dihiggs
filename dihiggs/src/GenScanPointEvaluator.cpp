#include "GenScanPointEvaluator.hpp"

#include "THDM.h"
#include "Constraints.h"
#include "DecayTable.h"

#include <cmath>
#include <limits>
#include <random>
#include <utility>

namespace gen_scan {

namespace {

// 2HDMC-style "h2" decay channel indices used throughout dihiggs/src/*.
constexpr int H2 = 2;
constexpr int H1 = 1;
constexpr int IDX_B = 3;
constexpr int IDX_TAU = 3;
constexpr int IDX_W = 3;
constexpr int IDX_Z = 2;

// Conversion constants (mlpython/2603 conventions).
constexpr double HBAR_GEV_S = 6.582119569e-25;  // GeV*s
constexpr double C_LIGHT_MM_S = 2.99792458e11;  // mm/s

double relerr(double value, double target) {
    if (target == 0.0) {
        return value;
    }
    return (value - target) / target;
}

void delta_ratio(double a, double b, double& delta, double& ratio) {
    delta = a - b;
    ratio = (b != 0.0) ? (a / b) : std::numeric_limits<double>::quiet_NaN();
}

// One generic-basis variation candidate. lambda6 is held fixed (= row value),
// lambda7 = 0, and lambda5 is tied to lambda4 (the m_A = m_Hp constraint); the
// free, ±variation_fraction-jittered quantities are lambda1..lambda4 and m12sq.
struct Candidate {
    double lambda1 = 0.0;
    double lambda2 = 0.0;
    double lambda3 = 0.0;
    double lambda4 = 0.0;
    double lambda5 = 0.0;
    double m12sq = 0.0;
};

// Sum-of-squared relative errors of the recovered physical spectrum against the
// bronze fixing's targets (mA, mH, mh=125, mHp=mA, sba=1). Smaller = the
// variation candidate reproduces Christopher's input masses more closely.
double mass_match_score(double mh, double mH, double mA, double mHp, double sba,
                        const BronzeRow& row) {
    const double e_mA = relerr(mA, row.m_A);
    const double e_mH = relerr(mH, row.m_H_input);
    const double e_mh = relerr(mh, 125.0);
    const double e_mHp = relerr(mHp, row.m_A);
    const double e_sba = relerr(sba, 1.0);
    return e_mA * e_mA + e_mH * e_mH + e_mh * e_mh + e_mHp * e_mHp + e_sba * e_sba;
}

// Build a 2HDMC model from one generic-basis candidate, evaluate constraints +
// full DecayTable, and fill all output/diagnostic fields (incl. chris
// cross-check deltas). variation_idx labels the candidate within the cloud
// (0 = unperturbed chris reconstruction).
GenFixingsPointResult evaluate_one_candidate(
    const BronzeRow& row, const Candidate& cand, const CalibrationConfig& cfg, int variation_idx) {
    GenFixingsPointResult res;
    res.variation_idx = variation_idx;
    res.mA_target = row.m_A;
    res.mH_target = row.m_H_input;
    res.calibration_n_used = cfg.n_samples;

    res.chris_width_bb = row.chris_width_bb;
    res.chris_width_tautau = row.chris_width_tautau;
    res.chris_width_gg = row.chris_width_gg;
    res.chris_width_gaga = row.chris_width_gaga;
    res.chris_width_Zga = row.chris_width_Zga;
    res.chris_ctau_mm = row.chris_ctau_m * 1000.0;

    THDM model;
    const bool set_ok = model.set_param_gen(
        cand.lambda1, cand.lambda2, cand.lambda3, cand.lambda4, cand.lambda5,
        row.lambda6, 0.0, cand.m12sq, row.tan_beta);

    if (!set_ok) {
        // 2HDMC rejected this generic-basis point. Leave the 29-col / stability
        // fields at their defaults (all zero) so the row is still emitted and is
        // trivially filterable via calibration_score==inf.
        res.tan_beta = row.tan_beta;
        res.lambda6 = row.lambda6;
        res.computed_lam1 = row.lambda1_recon;
        res.computed_lam2 = row.lambda2_recon;
        res.calibration_score = std::numeric_limits<double>::infinity();

        res.meta = replay_safe_output::make_metadata(cand.m12sq, 0.0, row.yukawa_type);
        res.meta.model_api_path = "THDM::set_param_gen->THDM::get_param_phys (variation, set_param_gen failed)";
        res.meta.replay_semantics_version = "gen_scan_variation_v1";

        delta_ratio(0.0, res.chris_width_bb, res.delta_width_bb, res.ratio_width_bb);
        delta_ratio(0.0, res.chris_width_tautau, res.delta_width_tautau, res.ratio_width_tautau);
        delta_ratio(0.0, res.chris_width_gg, res.delta_width_gg, res.ratio_width_gg);
        delta_ratio(0.0, res.chris_width_gaga, res.delta_width_gaga, res.ratio_width_gaga);
        delta_ratio(0.0, res.chris_width_Zga, res.delta_width_Zga, res.ratio_width_Zga);
        delta_ratio(0.0, res.chris_ctau_mm, res.delta_ctau_mm, res.ratio_ctau_mm);
        return res;
    }

    model.set_yukawas_type(row.yukawa_type);

    double mh, mH, mA, mHp, sba, lambda6_out, lambda7_out, m12_2_out, tanb_out;
    model.get_param_phys(mh, mH, mA, mHp, sba, lambda6_out, lambda7_out, m12_2_out, tanb_out);

    double l1, l2, l3, l4, l5, l6, l7, m12_2_gen, tanb_gen;
    model.get_param_gen(l1, l2, l3, l4, l5, l6, l7, m12_2_gen, tanb_gen);

    res.calibration_score = mass_match_score(mh, mH, mA, mHp, sba, row);

    res.mh_calibrated = mh;
    res.mHp_calibrated = mHp;
    res.sba_calibrated = sba;

    res.m_phi = mH;
    res.mA = mA;
    res.beta = std::atan(tanb_gen);
    res.alpha = res.beta - std::asin(sba);
    res.lambda6 = l6;
    res.lambda7 = l7;
    res.m12 = m12_2_gen;
    res.sin_ba = sba;
    res.tan_beta = tanb_gen;

    res.lam1 = l1;
    res.lam2 = l2;
    res.lam3 = l3;
    res.lam4 = l4;
    res.lam5 = l5;
    res.computed_lam1 = row.lambda1_recon;
    res.computed_lam2 = row.lambda2_recon;

    Constraints check(model);
    res.positivity_ok = check.check_positivity() ? 1 : 0;
    res.unitarity_ok = check.check_unitarity() ? 1 : 0;
    res.perturbativity_ok = check.check_perturbativity() ? 1 : 0;

    if (res.positivity_ok && res.unitarity_ok && res.perturbativity_ok) {
        res.stability_ok = check.check_stability() ? 1 : 0;
    } else {
        res.stability_ok = 0;
    }

    DecayTable table(model);
    res.width_bb = table.get_gamma_hdd(H2, IDX_B, IDX_B);
    res.width_tautau = table.get_gamma_hll(H2, IDX_TAU, IDX_TAU);
    res.width_WW = table.get_gamma_hvv(H2, IDX_W);
    res.width_ZZ = table.get_gamma_hvv(H2, IDX_Z);
    res.width_gaga = table.get_gamma_hgaga(H2);
    res.width_Zga = table.get_gamma_hZga(H2);
    res.width_gg = table.get_gamma_hgg(H2);
    res.width_hh = table.get_gamma_hhh(H2, H1, H1);
    res.total_width = table.get_gammatot_h(H2);
    res.br_gaga = (res.total_width > 1e-15) ? (res.width_gaga / res.total_width) : 0.0;

    const double ctau_mm = (res.total_width > 0.0)
        ? (HBAR_GEV_S * C_LIGHT_MM_S) / res.total_width
        : std::numeric_limits<double>::infinity();

    delta_ratio(res.width_bb, res.chris_width_bb, res.delta_width_bb, res.ratio_width_bb);
    delta_ratio(res.width_tautau, res.chris_width_tautau, res.delta_width_tautau, res.ratio_width_tautau);
    delta_ratio(res.width_gg, res.chris_width_gg, res.delta_width_gg, res.ratio_width_gg);
    delta_ratio(res.width_gaga, res.chris_width_gaga, res.delta_width_gaga, res.ratio_width_gaga);
    delta_ratio(res.width_Zga, res.chris_width_Zga, res.delta_width_Zga, res.ratio_width_Zga);
    delta_ratio(ctau_mm, res.chris_ctau_mm, res.delta_ctau_mm, res.ratio_ctau_mm);

    res.meta = replay_safe_output::make_metadata(cand.m12sq, m12_2_gen, row.yukawa_type);
    res.meta.model_api_path = "THDM::set_param_gen->THDM::get_param_phys (variation)";
    res.meta.replay_semantics_version = "gen_scan_variation_v1";

    return res;
}

}  // namespace

std::vector<GenFixingsPointResult> evaluate_gen_fixings_point(
    const BronzeRow& row, const CalibrationConfig& cfg, unsigned int thread_rng_seed) {
    std::mt19937 rng(thread_rng_seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    // Baseline = the unperturbed chris reconstruction. Stage 1 enforces
    // m_A = m_Hp, so lambda5_recon already equals lambda4_recon.
    const Candidate baseline{row.lambda1_recon, row.lambda2_recon, row.lambda3_recon,
                             row.lambda4_recon, row.lambda5_recon, row.m12sq_recon};

    const int n = (cfg.n_samples > 0) ? cfg.n_samples : 1;
    std::vector<GenFixingsPointResult> cloud;
    cloud.reserve(static_cast<size_t>(n));

    for (int i = 0; i < n; ++i) {
        Candidate cand;
        if (i == 0) {
            // First variation point is always the unperturbed reconstruction.
            cand = baseline;
        } else {
            auto jitter = [&](double base) {
                const double lo = (1.0 - cfg.variation_fraction) * base;
                const double hi = (1.0 + cfg.variation_fraction) * base;
                return (lo <= hi) ? (lo + (hi - lo) * unit(rng)) : (hi + (lo - hi) * unit(rng));
            };
            // Vary lambda1..lambda4 and m12sq by ±variation_fraction; lambda5 is
            // tied to lambda4 (m_A = m_Hp), lambda6 is held fixed, lambda7 = 0.
            cand.lambda1 = jitter(baseline.lambda1);
            cand.lambda2 = jitter(baseline.lambda2);
            cand.lambda3 = jitter(baseline.lambda3);
            cand.lambda4 = jitter(baseline.lambda4);
            cand.lambda5 = cand.lambda4;
            cand.m12sq = jitter(baseline.m12sq);
        }
        cloud.push_back(evaluate_one_candidate(row, cand, cfg, i));
    }

    return cloud;
}

std::vector<std::string> output_csv_columns() {
    std::vector<std::string> cols = {
        "m_phi", "mA", "alpha", "beta", "lambda6", "lambda7", "m12", "sin_ba", "tan_beta",
        "positivity_ok", "unitarity_ok", "perturbativity_ok",
        "width_bb", "width_tautau", "width_WW", "width_ZZ", "width_gaga", "width_Zga", "width_gg", "width_hh",
        "total_width", "br_gaga",
        "lam1", "computed_lam1", "lam2", "computed_lam2", "lam3", "lam4", "lam5",
        "mA_target", "mH_target", "mh_calibrated", "mHp_calibrated", "sba_calibrated",
        "calibration_score", "calibration_n_used", "variation_idx",
        "stability_ok",
        "chris_width_bb", "chris_width_tautau", "chris_width_gg", "chris_width_gaga", "chris_width_Zga", "chris_ctau_mm",
        "delta_width_bb", "ratio_width_bb",
        "delta_width_tautau", "ratio_width_tautau",
        "delta_width_gg", "ratio_width_gg",
        "delta_width_gaga", "ratio_width_gaga",
        "delta_width_Zga", "ratio_width_Zga",
        "delta_ctau_mm", "ratio_ctau_mm",
    };
    const std::vector<std::string> meta_cols = replay_safe_output::metadata_column_names();
    cols.insert(cols.end(), meta_cols.begin(), meta_cols.end());
    return cols;
}

}  // namespace gen_scan
