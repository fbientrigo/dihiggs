#pragma once

// GenScanPointEvaluator.hpp
// =========================
//
// "Stage 2" point evaluator for GenScanWithFixings: takes one row of a
// chris/CalcLambda1ScanFixings bronze shard (fast-reconstructed generic-basis
// lambdas + chris's 6-channel widths/ctau) and produces a *cloud* of
// n_samples generic-basis variation points: lambda1..lambda4 and m12sq are
// jittered by ±variation_fraction around the bronze reconstruction, with
// lambda5 tied to lambda4 (m_A = m_Hp), lambda6 held fixed, lambda7 = 0
// (variation_idx 0 = unperturbed reconstruction). Each candidate is fed to
// 2HDMC::THDM::set_param_gen, recovered to the physical basis, and fully
// re-evaluated with Constraints (incl. stability) and DecayTable (incl. H->hh).
// Every candidate is emitted (the 29-column schema plus calibration, stability,
// variation, and chris-cross-check diagnostic columns); calibration_score
// reports how closely that candidate reproduces the bronze fixing's target
// masses (mA, mH, mh=125, mHp=mA, sba=1).

#include "ReplaySafeOutput.hpp"

#include <string>
#include <vector>

namespace gen_scan {

// Subset of bronze shard columns Stage 2 needs (see chris/CalcLambda1ScanFixings.cpp
// for the full 51-column bronze schema).
struct BronzeRow {
    double tan_beta = 0.0;
    double m_A = 0.0;
    double lambda6 = 0.0;
    double lambda7 = 0.0;
    double lambda1_input = 0.0;
    double sin_ba = 0.0;
    double mh_input = 0.0;
    int    yukawa_type = 1;
    double m_H_input = 0.0;     // = mH_target

    double m12sq_recon = 0.0;
    double lambda1_recon = 0.0;
    double lambda2_recon = 0.0;
    double lambda3_recon = 0.0;
    double lambda4_recon = 0.0;
    double lambda5_recon = 0.0;

    double chris_width_bb = 0.0;
    double chris_width_tautau = 0.0;
    double chris_width_gg = 0.0;
    double chris_width_gaga = 0.0;
    double chris_width_Zga = 0.0;
    double chris_total_width = 0.0;
    double chris_ctau_m = 0.0;  // meters
};

struct CalibrationConfig {
    int n_samples = 50;
    double variation_fraction = 0.10;
    unsigned int rng_seed = 0;
};

struct GenFixingsPointResult {
    // ---- 1-29: legacy schema (PARQUET_SCHEMA.md / PhysScanWithFixings.cpp) ----
    double m_phi = 0.0, mA = 0.0, alpha = 0.0, beta = 0.0;
    double lambda6 = 0.0, lambda7 = 0.0, m12 = 0.0;
    double sin_ba = 0.0, tan_beta = 0.0;
    int positivity_ok = 0, unitarity_ok = 0, perturbativity_ok = 0;
    double width_bb = 0.0, width_tautau = 0.0, width_WW = 0.0, width_ZZ = 0.0;
    double width_gaga = 0.0, width_Zga = 0.0, width_gg = 0.0, width_hh = 0.0;
    double total_width = 0.0, br_gaga = 0.0;
    double lam1 = 0.0, computed_lam1 = 0.0;
    double lam2 = 0.0, computed_lam2 = 0.0;
    double lam3 = 0.0, lam4 = 0.0, lam5 = 0.0;

    // ---- calibration diagnostics ----
    double mA_target = 0.0, mH_target = 0.0;
    double mh_calibrated = 0.0, mHp_calibrated = 0.0, sba_calibrated = 0.0;
    double calibration_score = 0.0;
    int calibration_n_used = 0;
    int variation_idx = 0;

    // ---- stability (only computed here) ----
    int stability_ok = 0;

    // ---- chris cross-check (6 channels chris implements) ----
    double chris_width_bb = 0.0, chris_width_tautau = 0.0, chris_width_gg = 0.0;
    double chris_width_gaga = 0.0, chris_width_Zga = 0.0, chris_ctau_mm = 0.0;
    double delta_width_bb = 0.0, ratio_width_bb = 0.0;
    double delta_width_tautau = 0.0, ratio_width_tautau = 0.0;
    double delta_width_gg = 0.0, ratio_width_gg = 0.0;
    double delta_width_gaga = 0.0, ratio_width_gaga = 0.0;
    double delta_width_Zga = 0.0, ratio_width_Zga = 0.0;
    double delta_ctau_mm = 0.0, ratio_ctau_mm = 0.0;

    replay_safe_output::Metadata meta;
};

// Generate and fully evaluate the ±variation_fraction / n_samples generic-basis
// variation cloud for one bronze row. Returns one GenFixingsPointResult per
// candidate (size == cfg.n_samples, clamped to >= 1).
std::vector<GenFixingsPointResult> evaluate_gen_fixings_point(
    const BronzeRow& row, const CalibrationConfig& cfg, unsigned int thread_rng_seed);

// Full ordered list of output CSV column names (29 legacy + diagnostics +
// replay-safety metadata). Shared between the header writer and any
// downstream consumers that need the schema.
std::vector<std::string> output_csv_columns();

}  // namespace gen_scan
