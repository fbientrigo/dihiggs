/*
 * CalcLambda1Core.c
 * =================
 *
 * Implementation of CalcLambda1Core.h.  This is a verbatim extraction of the
 * per-point physics core from chris/CalcLambda1.c (CalcLambda1Standalone):
 * angle conventions, the lambda1<->m12^2 inversion, generic-basis lambda
 * reconstruction, the inverse mass self-consistency check, the 6-channel
 * decay-width calculator (bb, cc, tautau, gg, gaga, Zga), and the Jacobi-based
 * tree-unitarity eigenvalue check.  chris/CalcLambda1.c itself is left
 * untouched; this file is compiled separately (gcc, C99, -fopenmp not
 * required) and linked into chris/CalcLambda1ScanFixings.cpp.
 *
 * The 801x801 stability-grid minimisation (calc_stability_min_V4) is
 * deliberately NOT included here -- see CalcLambda1Core.h.
 */

#include "CalcLambda1Core.h"

#include <math.h>
#include <complex.h>
#include <float.h>
#include <string.h>

/* ------------------------------------------------------------
 * Constants (verbatim from CalcLambda1.c)
 * ------------------------------------------------------------ */

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

#define PI M_PI
#define HBARC_GEV_M 1.973269804e-16

#define GF_2HDMC 1.16637e-5
#define V_SM (1.0/sqrt(sqrt(2.0)*GF_2HDMC))
#define MZ 91.1876
#define MW 80.379
#define ALPHA_EM (1.0/137.035999084)
#define ALPHA_S_MZ 0.1181

#define MB_MB 4.18
#define MB_POLE 4.78
#define MC_MC 1.27
#define MC_POLE 1.67
#define MTAU 1.77686
#define MT_POLE 172.76
#define MT_MSBAR_REF 163.0

#define GAMMA_B_2HDMC 0.5486721046913843
#define GAMMA_C_2HDMC 0.6092741869074962

#define PERTURBATIVITY_LIMIT (4.0*PI)
#define UNITARITY_LIMIT      (16.0*PI)

/* ------------------------------------------------------------
 * Math utilities
 * ------------------------------------------------------------ */

static double sqr(double x) { return x*x; }

static double cube(double x) { return x*x*x; }

static double clamp_double(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static double cabs2(double complex z) {
    double a = creal(z);
    double b = cimag(z);
    return a*a + b*b;
}

/* ------------------------------------------------------------
 * Angles from (tan(beta), sin(beta-alpha)).
 *
 * Stable for tan(beta) ~ 1e7: avoids computing sin(beta)/cos(beta) from
 * beta=atan(tanb) directly, instead using cb=1/sqrt(1+tanb^2).
 * 2HDMC convention: alpha = -asin(sba) + beta, cba = sqrt(1-sba^2).
 * ------------------------------------------------------------ */

AngleSet angles_from_tanb_sba(double tanb, double sba_in) {
    AngleSet A;

    A.sba = clamp_double(sba_in, -1.0, 1.0);
    A.cba = sqrt(fmax(0.0, 1.0 - A.sba*A.sba));

    A.beta = atan(tanb);

    A.cb = 1.0 / sqrt(1.0 + tanb*tanb);
    A.sb = tanb * A.cb;

    A.alpha = A.beta - asin(A.sba);

    A.sa = A.sb*A.cba - A.cb*A.sba;
    A.ca = A.cb*A.cba + A.sb*A.sba;

    return A;
}

/* ------------------------------------------------------------
 * m12^2 from lambda_1: exact inverse of THDM::set_param_phys's
 * lambda_1 = (M11 - m12^2 tanb)/(v^2 cb^2) - 1.5 lambda6 tanb + 0.5 lambda7 tanb^3
 * ------------------------------------------------------------ */

double m12sq_from_lambda1(InputPoint in) {
    AngleSet A = angles_from_tanb_sba(in.tanb, in.sba);

    double tb = in.tanb;
    double v2 = sqr(V_SM);

    double M11 = sqr(in.mH)*sqr(A.ca) + sqr(in.mh)*sqr(A.sa);

    return (
        M11
        - v2*sqr(A.cb)*(in.lambda1 + 1.5*in.lambda6*tb - 0.5*in.lambda7*cube(tb))
    ) / tb;
}

/* ------------------------------------------------------------
 * Reconstruct (lambda1..lambda7, m12sq, alpha, beta) from
 * (mh,mH,mA,mHp,sba,l6,l7,lambda1,tanb) using the same algebraic relations
 * as THDM::set_param_phys.
 * ------------------------------------------------------------ */

LambdaSet reconstruct_lambdas(InputPoint in) {
    LambdaSet L;

    AngleSet A = angles_from_tanb_sba(in.tanb, in.sba);

    double sb = A.sb;
    double cb = A.cb;
    double tb = in.tanb;
    double ctb = 1.0 / tb;

    double sa = A.sa;
    double ca = A.ca;

    double v2 = sqr(V_SM);
    double m12sq = m12sq_from_lambda1(in);

    L.lambda1 = in.lambda1;

    L.lambda2 = (sqr(in.mH)*sqr(sa) + sqr(in.mh)*sqr(ca) - m12sq*ctb)/(v2*sqr(sb))
                + 0.5*in.lambda6*cube(ctb)
                - 1.5*in.lambda7*ctb;

    L.lambda3 = ((sqr(in.mH)-sqr(in.mh))*ca*sa + 2.0*sqr(in.mHp)*sb*cb - m12sq)/(v2*sb*cb)
                - 0.5*in.lambda6*ctb
                - 0.5*in.lambda7*tb;

    L.lambda4 = ((sqr(in.mA)-2.0*sqr(in.mHp))*cb*sb + m12sq)/(v2*sb*cb)
                - 0.5*in.lambda6*ctb
                - 0.5*in.lambda7*tb;

    L.lambda5 = (m12sq - sqr(in.mA)*sb*cb)/(v2*sb*cb)
                - 0.5*in.lambda6*ctb
                - 0.5*in.lambda7*tb;

    L.lambda6 = in.lambda6;
    L.lambda7 = in.lambda7;
    L.m12sq = m12sq;
    L.alpha = A.alpha;
    L.beta = A.beta;

    return L;
}

/* ------------------------------------------------------------
 * Inverse mass self-consistency check: recompute (mh,mH,mA,mHp,sba) from
 * the reconstructed LambdaSet and compare against the InputPoint masses.
 * ------------------------------------------------------------ */

MassCheck masses_from_lambdas(InputPoint in, LambdaSet L) {
    MassCheck M;

    AngleSet A = angles_from_tanb_sba(in.tanb, in.sba);

    double sb = A.sb;
    double cb = A.cb;
    double tb = in.tanb;
    double ctb = 1.0/tb;
    double v2 = sqr(V_SM);
    double sbcb = sb*cb;

    double mA2 = L.m12sq/(sbcb)
               - 0.5*v2*(2.0*L.lambda5 + L.lambda6*ctb + L.lambda7*tb);

    double mHp2 = mA2 - 0.5*v2*(L.lambda4 - L.lambda5);

    double M11 = v2*sqr(cb)*(L.lambda1 + 1.5*L.lambda6*tb - 0.5*L.lambda7*cube(tb))
               + L.m12sq*tb;

    double M22 = v2*sqr(sb)*(L.lambda2 - 0.5*L.lambda6*cube(ctb) + 1.5*L.lambda7*ctb)
               + L.m12sq*ctb;

    double M12 = (L.lambda3 + 0.5*L.lambda6*ctb + 0.5*L.lambda7*tb)*v2*sbcb
               - 2.0*mHp2*sbcb
               + L.m12sq;

    double tr = M11 + M22;
    double det_term = sqrt(fmax(0.0, sqr(M11-M22) + 4.0*sqr(M12)));

    double mh2 = 0.5*(tr - det_term);
    double mH2 = 0.5*(tr + det_term);

    M.mh = sqrt(fmax(0.0, mh2));
    M.mH = sqrt(fmax(0.0, mH2));
    M.mA = sqrt(fmax(0.0, mA2));
    M.mHp = sqrt(fmax(0.0, mHp2));

    double alpha_diag = 0.5*atan2(2.0*M12, M11-M22);
    M.sba = sin(A.beta - alpha_diag);

    return M;
}

/* ------------------------------------------------------------
 * Running alpha_s and quark masses (LO, as used in the notebooks)
 * ------------------------------------------------------------ */

static int get_Nactivef(double mu) {
    if (mu >= MT_POLE) return 6;
    if (mu >= MB_POLE) return 5;
    if (mu >= MC_POLE) return 4;
    return 3;
}

static double alpha_s_lo(double mu, int nf) {
    if (mu < 1.0) mu = 1.0;
    double beta0 = 11.0 - 2.0*((double)nf)/3.0;
    return ALPHA_S_MZ / (1.0 + beta0*ALPHA_S_MZ/(2.0*PI)*log(mu/MZ));
}

static double mb_running(double mu) {
    double mu_eff = (mu > MB_MB ? mu : MB_MB);
    return MB_MB * pow(alpha_s_lo(mu_eff,5)/alpha_s_lo(MB_MB,5), GAMMA_B_2HDMC);
}

static double mc_running(double mu) {
    double mu_eff = (mu > 2.0 ? mu : 2.0);
    return MC_MC * pow(alpha_s_lo(mu_eff,4)/alpha_s_lo(MC_MC,4), GAMMA_C_2HDMC);
}

/* ------------------------------------------------------------
 * Complex loop functions
 * ------------------------------------------------------------ */

static double complex ftau(double complex t) {
    if (creal(t) <= 1.0 && fabs(cimag(t)) < 1e-15) {
        double tr = creal(t);
        if (tr < 0.0) tr = 0.0;
        double a = asin(sqrt(tr));
        return a*a + 0.0*I;
    }

    double complex root = csqrt(t - 1.0);
    double complex x = clog((csqrt(t) + root)/(csqrt(t) - root)) - I*PI;
    return -0.25*x*x;
}

static double complex gtau(double complex t) {
    if (creal(t) <= 1.0 && fabs(cimag(t)) < 1e-15) {
        double tr = creal(t);
        if (tr <= 0.0) tr = DBL_MIN;
        return asin(sqrt(tr)) * sqrt(1.0/tr - 1.0) + 0.0*I;
    }

    double complex x = clog((csqrt(t) + csqrt(t - 1.0))/(csqrt(t) - csqrt(t - 1.0))) - I*PI;
    return 0.5*x*csqrt(1.0 - 1.0/t);
}

static double complex F_sf(double t) {
    double ti = 1.0/t;
    return ti*(1.0 + (1.0 - ti)*ftau(t + 0.0*I));
}

static double complex F_0(double t) {
    double ti = 1.0/t;
    return ti*(-1.0 + ti*ftau(t + 0.0*I));
}

static double complex I_2(double tau, double lam) {
    double denom = tau - lam;
    if (fabs(denom) < 1e-14) denom = (denom >= 0.0 ? 1e-14 : -1e-14);
    return -tau*lam/(2.0*denom) * (ftau(1.0/tau + 0.0*I) - ftau(1.0/lam + 0.0*I));
}

static double complex I_1(double tau, double lam) {
    double denom = tau - lam;
    if (fabs(denom) < 1e-14) denom = (denom >= 0.0 ? 1e-14 : -1e-14);

    double complex term1 = tau*lam/(2.0*denom);
    double complex term2 = sqr(tau)*sqr(lam)/(2.0*sqr(denom)) *
                           (ftau(1.0/tau + 0.0*I) - ftau(1.0/lam + 0.0*I));
    double complex term3 = sqr(tau)*lam/sqr(denom) *
                           (gtau(1.0/tau + 0.0*I) - gtau(1.0/lam + 0.0*I));
    return term1 + term2 + term3;
}

static double complex FF_s(double tau, double lam) {
    return I_1(tau, lam) - I_2(tau, lam);
}

static double complex FHp(double tau, double lam) {
    return I_1(tau, lam);
}

/* ------------------------------------------------------------
 * Decay widths (Type-I-like: epsV=0, epsA=1/tan(beta))
 * ------------------------------------------------------------ */

static double K_hff_scalar(double mphi) {
    int Nf = get_Nactivef(mphi);
    double a = alpha_s_lo(mphi, Nf)/PI;
    double K = 1.0 + 5.67*a + (35.94 - 1.36*Nf)*a*a;
    return (K > 0.0 ? K : 0.0);
}

static double width_ff_scalar(double mphi, double mf_pole, double mf_run,
                              int Nc, double epsV, double epsA,
                              int apply_qcd) {
    if (mphi <= 2.0*mf_pole) return 0.0;

    double kappa_f = epsV - epsA;
    double beta_f = sqrt(1.0 - 4.0*sqr(mf_pole)/sqr(mphi));

    double G = ((double)Nc) * sqr(fabs(kappa_f)) * sqr(mf_run) * mphi
               /(8.0*PI*sqr(V_SM)) * cube(beta_f);

    if (apply_qcd && Nc == 3) G *= K_hff_scalar(mphi);

    return G;
}

static double width_bb(double mphi, double epsV, double epsA) {
    return width_ff_scalar(mphi, MB_POLE, mb_running(mphi), 3, epsV, epsA, 1);
}

static double width_cc(double mphi, double epsV, double epsA) {
    return width_ff_scalar(mphi, MC_POLE, mc_running(mphi), 3, epsV, epsA, 1);
}

static double width_tautau(double mphi, double epsV, double epsA) {
    return width_ff_scalar(mphi, MTAU, MTAU, 1, epsV, epsA, 0);
}

static double Kgg_nnlo_like(double mphi) {
    int Nf = get_Nactivef(mphi);
    double a_s = alpha_s_lo(mphi, Nf);
    return 1.0
           + a_s/PI*(95.0/4.0 - 7.0*((double)Nf)/6.0)
           + sqr(a_s/PI)*(156.808 - 5.708*log(sqr(MT_MSBAR_REF)/sqr(mphi)));
}

static double width_gg(double mphi, double epsV, double epsA) {
    double kappa_f = epsV - epsA;

    double complex S = 0.0 + 0.0*I;
    double masses[3] = {MT_POLE, MB_POLE, MC_POLE};

    for (int i=0;i<3;i++) {
        double tau = sqr(mphi)/(4.0*sqr(masses[i]));
        S += kappa_f * F_sf(tau);
    }

    double a_s = alpha_s_lo(mphi, get_Nactivef(mphi));
    double pref = cube(mphi)*sqr(a_s)/(32.0*cube(PI)*sqr(V_SM));
    return pref*Kgg_nnlo_like(mphi)*cabs2(S);
}

static double lambda_eff_Hp_loop(double lambda6, double tanb, double lambda1) {
    return lambda6 + lambda1/tanb;
}

static double width_gaga(double mphi, double mHp, double lambda6,
                         double epsV, double epsA, double tanb, double lambda1) {
    double kappa_f = epsV - epsA;

    double complex S = 0.0 + 0.0*I;

    double mf[4] = {MT_POLE, MB_POLE, MC_POLE, MTAU};
    double Qf[4] = {2.0/3.0, -1.0/3.0, 2.0/3.0, -1.0};
    int Nc[4] = {3,3,3,1};

    for (int i=0;i<4;i++) {
        double tau = sqr(mphi)/(4.0*sqr(mf[i]));
        double complex term = 2.0*((double)Nc[i])*sqr(Qf[i])*kappa_f*F_sf(tau);

        if (i == 0) {
            term *= (1.0 - alpha_s_lo(mphi, get_Nactivef(mphi))/PI);
        }

        S += term;
    }

    double tau_hp = sqr(mphi)/(4.0*sqr(mHp));
    double lambda_eff = lambda_eff_Hp_loop(lambda6, tanb, lambda1);
    double kappa_hp = lambda_eff*sqr(V_SM)/(2.0*sqr(mHp));
    S += kappa_hp*F_0(tau_hp);

    return sqr(ALPHA_EM)*cube(mphi)/(256.0*cube(PI)*sqr(V_SM))*cabs2(S);
}

static double width_Zga(double mphi, double mHp, double lambda6,
                        double epsV, double epsA, double tanb, double lambda1) {
    if (mphi <= MZ) return 0.0;

    double sW2 = 1.0 - sqr(MW/MZ);
    double cW = sqrt(1.0 - sW2);
    double kappa_f = epsV - epsA;

    double complex S = 0.0 + 0.0*I;

    double mf[4] = {MB_POLE, MC_POLE, MT_POLE, MTAU};
    double Qf[4] = {-1.0/3.0, 2.0/3.0, 2.0/3.0, -1.0};
    double I3f[4] = {-0.5, 0.5, 0.5, -0.5};
    int Nc[4] = {3,3,3,1};

    for (int i=0;i<4;i++) {
        double tau = 4.0*sqr(mf[i])/sqr(mphi);
        double lam = 4.0*sqr(mf[i])/sqr(MZ);

        S += 2.0*((double)Nc[i])*Qf[i]*(I3f[i] - 2.0*sW2*Qf[i])/cW
             * kappa_f * FF_s(tau, lam);
    }

    double tau_hp = 4.0*sqr(mHp)/sqr(mphi);
    double lam_hp = 4.0*sqr(mHp)/sqr(MZ);

    double lambda_eff = lambda_eff_Hp_loop(lambda6, tanb, lambda1);
    double kappa_hp = lambda_eff*sqr(V_SM)/(2.0*sqr(mHp));

    S -= (2.0*cW - 1.0/cW)*kappa_hp*FHp(tau_hp, lam_hp);

    double GF = 1.0/(sqrt(2.0)*sqr(V_SM));
    double phase = cube(1.0 - sqr(MZ)/sqr(mphi));
    double pref = ALPHA_EM*sqr(GF)*sqr(MW)*cube(mphi)*phase/(64.0*pow(PI,4.0));

    return pref*cabs2(S);
}

WidthResult compute_widths(InputPoint in) {
    WidthResult W;

    double epsV = 0.0;
    double epsA = 1.0/in.tanb;

    W.G_bb = width_bb(in.mH, epsV, epsA);
    W.G_cc = width_cc(in.mH, epsV, epsA);
    W.G_tautau = width_tautau(in.mH, epsV, epsA);
    W.G_gg = width_gg(in.mH, epsV, epsA);
    W.G_gaga = width_gaga(in.mH, in.mHp, in.lambda6, epsV, epsA, in.tanb, in.lambda1);
    W.G_Zga = width_Zga(in.mH, in.mHp, in.lambda6, epsV, epsA, in.tanb, in.lambda1);

    W.G_total = W.G_bb + W.G_cc + W.G_tautau + W.G_gg + W.G_gaga + W.G_Zga;

    if (W.G_total > 0.0) {
        W.BR_bb = W.G_bb/W.G_total;
        W.BR_cc = W.G_cc/W.G_total;
        W.BR_tautau = W.G_tautau/W.G_total;
        W.BR_gg = W.G_gg/W.G_total;
        W.BR_gaga = W.G_gaga/W.G_total;
        W.BR_Zga = W.G_Zga/W.G_total;
        W.BR_loop = W.BR_gaga + W.BR_Zga;
        W.ctau_m = HBARC_GEV_M/W.G_total;
    } else {
        W.BR_bb = W.BR_cc = W.BR_tautau = W.BR_gg = W.BR_gaga = W.BR_Zga = W.BR_loop = 0.0;
        W.ctau_m = INFINITY;
    }

    return W;
}

/* ------------------------------------------------------------
 * Jacobi eigenvalues for small symmetric matrices (used for the
 * tree-unitarity eigenvalue check).
 * ------------------------------------------------------------ */

static void jacobi_eigenvalues(int n, double A[4][4], double eval[4]) {
    for (int iter=0; iter<200; iter++) {
        int p=0, q=1;
        double max_off = 0.0;

        for (int i=0;i<n;i++) {
            for (int j=i+1;j<n;j++) {
                double val = fabs(A[i][j]);
                if (val > max_off) {
                    max_off = val;
                    p = i;
                    q = j;
                }
            }
        }

        if (max_off < 1e-13) break;

        double app = A[p][p];
        double aqq = A[q][q];
        double apq = A[p][q];

        double tau = (aqq - app)/(2.0*apq);
        double t = (tau >= 0.0 ? 1.0 : -1.0)/(fabs(tau) + sqrt(1.0 + tau*tau));
        double c = 1.0/sqrt(1.0 + t*t);
        double s = t*c;

        A[p][p] = app - t*apq;
        A[q][q] = aqq + t*apq;
        A[p][q] = A[q][p] = 0.0;

        for (int k=0;k<n;k++) {
            if (k == p || k == q) continue;
            double akp = A[k][p];
            double akq = A[k][q];
            A[k][p] = A[p][k] = c*akp - s*akq;
            A[k][q] = A[q][k] = s*akp + c*akq;
        }
    }

    for (int i=0;i<n;i++) eval[i] = A[i][i];
}

/* ------------------------------------------------------------
 * Tree-level unitarity: max |eigenvalue| over the S21 (3x3), S01 (4x4)
 * and S00 (4x4) scattering matrices, plus the |lambda3-lambda4| seed.
 * ------------------------------------------------------------ */

double calc_unitarity_max_abs(LambdaSet L) {
    double l1=L.lambda1, l2=L.lambda2, l3=L.lambda3, l4=L.lambda4, l5=L.lambda5, l6=L.lambda6, l7=L.lambda7;
    double egmax = fabs(l3 - l4);
    double s2 = sqrt(2.0);

    double A3[4][4] = {{0}};
    double A4[4][4] = {{0}};
    double A0[4][4] = {{0}};
    double eval[4] = {0,0,0,0};

    /* S21 3x3 */
    A3[0][0]=l1;      A3[0][1]=l5;      A3[0][2]=s2*l6;
    A3[1][0]=l5;      A3[1][1]=l2;      A3[1][2]=s2*l7;
    A3[2][0]=s2*l6;   A3[2][1]=s2*l7;   A3[2][2]=l3+l4;

    jacobi_eigenvalues(3, A3, eval);
    for (int i=0;i<3;i++) if (fabs(eval[i]) > egmax) egmax = fabs(eval[i]);

    /* S01 4x4 */
    memset(A4, 0, sizeof(A4));
    A4[0][0]=l1; A4[0][1]=l4; A4[0][2]=l6; A4[0][3]=l6;
    A4[1][0]=l4; A4[1][1]=l2; A4[1][2]=l7; A4[1][3]=l7;
    A4[2][0]=l6; A4[2][1]=l7; A4[2][2]=l3; A4[2][3]=l5;
    A4[3][0]=l6; A4[3][1]=l7; A4[3][2]=l5; A4[3][3]=l3;

    jacobi_eigenvalues(4, A4, eval);
    for (int i=0;i<4;i++) if (fabs(eval[i]) > egmax) egmax = fabs(eval[i]);

    /* S00 4x4 */
    A0[0][0]=3.0*l1;        A0[0][1]=2.0*l3+l4;     A0[0][2]=3.0*l6;        A0[0][3]=3.0*l6;
    A0[1][0]=2.0*l3+l4;     A0[1][1]=3.0*l2;        A0[1][2]=3.0*l7;        A0[1][3]=3.0*l7;
    A0[2][0]=3.0*l6;        A0[2][1]=3.0*l7;        A0[2][2]=l3+2.0*l4;    A0[2][3]=3.0*l5;
    A0[3][0]=3.0*l6;        A0[3][1]=3.0*l7;        A0[3][2]=3.0*l5;        A0[3][3]=l3+2.0*l4;

    jacobi_eigenvalues(4, A0, eval);
    for (int i=0;i<4;i++) if (fabs(eval[i]) > egmax) egmax = fabs(eval[i]);

    return egmax;
}

/* ------------------------------------------------------------
 * Fast constraint check (positivity, perturbativity, unitarity) --
 * stability (801x801 grid) is intentionally omitted; see header.
 * ------------------------------------------------------------ */

ConstraintResultFast check_constraints_fast(InputPoint in, LambdaSet L) {
    ConstraintResultFast C;

    C.masses_positive = (in.mh > 0.0 && in.mH > 0.0 && in.mA > 0.0 && in.mHp > 0.0);

    double lambdas[7] = {L.lambda1,L.lambda2,L.lambda3,L.lambda4,L.lambda5,L.lambda6,L.lambda7};
    C.max_abs_lambda = 0.0;
    int all_finite = 1;

    for (int i=0;i<7;i++) {
        if (!isfinite(lambdas[i])) all_finite = 0;
        if (fabs(lambdas[i]) > C.max_abs_lambda) C.max_abs_lambda = fabs(lambdas[i]);
    }

    C.perturbativity_ok = (all_finite && C.max_abs_lambda < PERTURBATIVITY_LIMIT);

    C.max_abs_unitarity_eigenvalue = calc_unitarity_max_abs(L);
    C.unitarity_ok = (isfinite(C.max_abs_unitarity_eigenvalue) &&
                      C.max_abs_unitarity_eigenvalue < UNITARITY_LIMIT);

    C.triple_ok_fast = (C.masses_positive && C.perturbativity_ok && C.unitarity_ok);

    return C;
}
