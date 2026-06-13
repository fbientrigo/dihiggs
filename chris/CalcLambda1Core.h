/*
 * CalcLambda1Core.h
 * =================
 *
 * Public, dependency-free interface to the per-point physics core lifted
 * from chris/CalcLambda1.c (CalcLambda1Standalone): lambda1->m12^2 inversion,
 * generic-basis lambda reconstruction, mass self-consistency check, and the
 * 6-channel decay-width calculator.
 *
 * No 2HDMC, GSL, or _Complex types appear in this header so that it can be
 * included from C++ via `extern "C"`.  The 801x801 stability grid
 * (calc_stability_min_V4) is intentionally NOT exposed here -- Stage 1
 * (CalcLambda1ScanFixings) skips it for speed; Stage 2 (GenScanWithFixings)
 * computes the authoritative stability check via 2HDMC::Constraints.
 */

#ifndef CALC_LAMBDA1_CORE_H
#define CALC_LAMBDA1_CORE_H

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------
 * Input point: physical-basis masses + lambda_1 + tan(beta).
 * ------------------------------------------------------------ */
typedef struct {
    double mh, mH, mA, mHp, sba;
    double lambda6, lambda7, lambda1;
    double tanb;
    int yukawa_type;
} InputPoint;

/* ------------------------------------------------------------
 * Reconstructed generic-basis couplings + angles.
 * ------------------------------------------------------------ */
typedef struct {
    double lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7;
    double m12sq;
    double alpha;
    double beta;
} LambdaSet;

/* ------------------------------------------------------------
 * Trig helper set for (tan(beta), sin(beta-alpha)).
 * ------------------------------------------------------------ */
typedef struct {
    double beta;
    double alpha;
    double sb;
    double cb;
    double sa;
    double ca;
    double sba;
    double cba;
} AngleSet;

/* ------------------------------------------------------------
 * Masses recovered from the reconstructed LambdaSet -- a
 * self-consistency check against the InputPoint masses.
 * ------------------------------------------------------------ */
typedef struct {
    double mh;
    double mH;
    double mA;
    double mHp;
    double sba;
} MassCheck;

/* ------------------------------------------------------------
 * 6-channel decay widths/BRs (bb, cc, tautau, gg, gaga, Zga) and ctau.
 * ------------------------------------------------------------ */
typedef struct {
    double G_bb;
    double G_cc;
    double G_tautau;
    double G_gg;
    double G_gaga;
    double G_Zga;
    double G_total;

    double BR_bb;
    double BR_cc;
    double BR_tautau;
    double BR_gg;
    double BR_gaga;
    double BR_Zga;
    double BR_loop;
    double ctau_m;
} WidthResult;

/* ------------------------------------------------------------
 * Fast (no stability grid) constraint check: positivity,
 * perturbativity (max|lambda_i|<4*pi), unitarity (Jacobi
 * eigenvalues < 16*pi).
 * ------------------------------------------------------------ */
typedef struct {
    int masses_positive;
    int perturbativity_ok;
    int unitarity_ok;
    int triple_ok_fast;

    double max_abs_lambda;
    double max_abs_unitarity_eigenvalue;
} ConstraintResultFast;

AngleSet angles_from_tanb_sba(double tanb, double sba_in);
double m12sq_from_lambda1(InputPoint in);
LambdaSet reconstruct_lambdas(InputPoint in);
MassCheck masses_from_lambdas(InputPoint in, LambdaSet L);
WidthResult compute_widths(InputPoint in);
double calc_unitarity_max_abs(LambdaSet L);
ConstraintResultFast check_constraints_fast(InputPoint in, LambdaSet L);

#ifdef __cplusplus
}
#endif

#endif /* CALC_LAMBDA1_CORE_H */
