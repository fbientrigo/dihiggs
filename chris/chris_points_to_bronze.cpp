/*
 * chris_points_to_bronze.cpp
 * ===========================
 *
 * Ad-hoc bronze-shard generator for the 7 fixed cross-check points in
 * chris_points.csv (mh,mH,mA,mHp,sba,lambda6,lambda7,lambda1_target,
 * tan_beta,yukawa_type). Produces a bronze CSV with the same 51-column
 * schema as CalcLambda1ScanFixings.cpp, one row per input point, so it can
 * be fed directly into dihiggs/app/GenScanWithFixings and the result
 * compared against chris_points_to_hybrid's 2HDMC ground truth for the same
 * points.
 *
 * Usage: chris_points_to_bronze <chris_points.csv> <output_bronze.csv>
 */

#include "CalcLambda1Core.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ',')) out.push_back(item);
    return out;
}

const std::vector<std::string>& bronze_columns() {
    static const std::vector<std::string> cols = {
        "tan_beta", "m_A", "m_Hp", "lambda6", "lambda7", "lambda1_input", "sin_ba", "mh_input", "yukawa_type",
        "m_H_input",
        "m12sq_recon",
        "lambda1_recon", "lambda2_recon", "lambda3_recon", "lambda4_recon", "lambda5_recon", "lambda6_recon", "lambda7_recon",
        "alpha_recon", "beta_recon",
        "mh_check", "mH_check", "mA_check", "mHp_check", "sba_check",
        "delta_mh", "delta_mH", "delta_mA", "delta_mHp", "delta_sba",
        "masses_positive", "perturbativity_ok", "unitarity_ok", "triple_ok_fast",
        "max_abs_lambda", "max_abs_unitarity_eigenvalue",
        "chris_width_bb", "chris_width_cc", "chris_width_tautau", "chris_width_gg", "chris_width_gaga", "chris_width_Zga",
        "chris_total_width", "chris_br_bb", "chris_br_cc", "chris_br_tautau", "chris_br_gg", "chris_br_gaga", "chris_br_Zga", "chris_br_loop",
        "chris_ctau_m"
    };
    return cols;
}

void write_row(std::ostream& out, const InputPoint& in) {
    LambdaSet L = reconstruct_lambdas(in);
    MassCheck M = masses_from_lambdas(in, L);
    WidthResult W = compute_widths(in);
    ConstraintResultFast C = check_constraints_fast(in, L);

    out << in.tanb << "," << in.mA << "," << in.mHp << "," << in.lambda6 << "," << in.lambda7 << ","
        << in.lambda1 << "," << in.sba << "," << in.mh << "," << in.yukawa_type << ","
        << in.mH << ","
        << L.m12sq << ","
        << L.lambda1 << "," << L.lambda2 << "," << L.lambda3 << "," << L.lambda4 << "," << L.lambda5 << "," << L.lambda6 << "," << L.lambda7 << ","
        << L.alpha << "," << L.beta << ","
        << M.mh << "," << M.mH << "," << M.mA << "," << M.mHp << "," << M.sba << ","
        << (M.mh - in.mh) << "," << (M.mH - in.mH) << "," << (M.mA - in.mA) << "," << (M.mHp - in.mHp) << "," << (M.sba - in.sba) << ","
        << C.masses_positive << "," << C.perturbativity_ok << "," << C.unitarity_ok << "," << C.triple_ok_fast << ","
        << C.max_abs_lambda << "," << C.max_abs_unitarity_eigenvalue << ","
        << W.G_bb << "," << W.G_cc << "," << W.G_tautau << "," << W.G_gg << "," << W.G_gaga << "," << W.G_Zga << ","
        << W.G_total << "," << W.BR_bb << "," << W.BR_cc << "," << W.BR_tautau << "," << W.BR_gg << "," << W.BR_gaga << "," << W.BR_Zga << "," << W.BR_loop << ","
        << W.ctau_m
        << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <chris_points.csv> <output_bronze.csv>\n";
        return 1;
    }

    std::ifstream fin(argv[1]);
    if (!fin) {
        std::cerr << "ERROR: cannot open input " << argv[1] << "\n";
        return 1;
    }

    std::ofstream fout(argv[2]);
    if (!fout) {
        std::cerr << "ERROR: cannot open output " << argv[2] << "\n";
        return 1;
    }
    fout.setf(std::ios::scientific);
    fout << std::setprecision(17);

    const std::vector<std::string> header = bronze_columns();
    for (size_t c = 0; c < header.size(); ++c) {
        fout << header[c] << (c + 1 < header.size() ? "," : "\n");
    }

    std::string line;
    std::getline(fin, line);  // header
    int n_rows = 0;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::vector<std::string> f = split_csv(line);
        if (f.size() != 10) {
            std::cerr << "ERROR: row with " << f.size() << " fields: " << line << "\n";
            return 1;
        }

        InputPoint in;
        in.mh = std::stod(f[0]);
        in.mH = std::stod(f[1]);
        in.mA = std::stod(f[2]);
        in.mHp = std::stod(f[3]);
        in.sba = std::stod(f[4]);
        in.lambda6 = std::stod(f[5]);
        in.lambda7 = std::stod(f[6]);
        in.lambda1 = std::stod(f[7]);
        in.tanb = std::stod(f[8]);
        in.yukawa_type = static_cast<int>(std::stod(f[9]));

        write_row(fout, in);
        ++n_rows;
    }

    fout.close();
    std::cout << "Wrote " << n_rows << " rows to " << argv[2] << "\n";
    return 0;
}
