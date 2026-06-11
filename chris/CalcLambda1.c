/*
 * CalcLambda1Standalone.c
 * =======================
 *
 * Archivo único y autocontenido: no necesita THDM.h, DecayTable.h, Constraints.h,
 * lib2HDMC, GSL ni ningún archivo externo.
 *
 * Uso:
 *   ./CalcLambda1Standalone mh mH mA mHp sin(beta-alpha) lambda_6 lambda_7 lambda_1 tan_beta yukawas_type output_filename
 *
 * Ejemplo:
 *   ./CalcLambda1Standalone 125 150.78947368421052 300 300 1 2.0176567952506178e-04 0 1 5.4523643594368286e6 1 salida_punto.txt
 *
 * Convenciones:
 *   - El input es igual a CalcPhys, pero reemplazando m12^2 por lambda_1.
 *   - El programa calcula m12^2 internamente desde lambda_1 usando la relación
 *     del elemento M11 de la matriz CP-par.
 *   - Se reconstruyen lambda_1,...,lambda_7 usando las mismas fórmulas
 *     algebraicas de la base física CP-conservante.
 *   - Los anchos se calculan con las funciones analíticas usadas en los notebooks:
 *       H -> b bbar, c cbar, tau tau, gg, gamma gamma, Z gamma.
 *   - El "triple check" autocontenido es:
 *       masses_positive,
 *       tree-level unitarity por autovalores de las matrices escalares,
 *       perturbativity por max |lambda_i| < 4*pi,
 *       stability por minimización numérica robusta del potencial cuártico.
 *
 * Advertencia:
 *   Este programa es autocontenido y no llama a 2HDMC. Por tanto, los constraints
 *   son una implementación analítica/numerica autónoma, no una llamada literal a
 *   Constraints::print_all. Está diseñado para diagnosticar y evitar dependencia
 *   externa, no para reemplazar una validación oficial con 2HDMC cuando esta se
 *   requiera en una publicación.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <string.h>

/* ------------------------------------------------------------
 * Constantes generales
 * ------------------------------------------------------------ */

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

#define PI M_PI
#define HBARC_GEV_M 1.973269804e-16

/* Constantes SM usadas en el notebook */
#define GF_2HDMC 1.16637e-5
#define V_SM (1.0/sqrt(sqrt(2.0)*GF_2HDMC))
#define MZ 91.1876
#define MW 80.379
#define ALPHA_EM (1.0/137.035999084)
#define ALPHA_S_MZ 0.1181

/* Masas usadas en el notebook */
#define MB_MB 4.18
#define MB_POLE 4.78
#define MC_MC 1.27
#define MC_POLE 1.67
#define MTAU 1.77686
#define MT_POLE 172.76
#define MT_MSBAR_REF 163.0

/* Exponentes ajustados usados en el notebook */
#define GAMMA_B_2HDMC 0.5486721046913843
#define GAMMA_C_2HDMC 0.6092741869074962

/* Límites de diagnóstico */
#define PERTURBATIVITY_LIMIT (4.0*PI)
#define UNITARITY_LIMIT      (16.0*PI)

/* Malla para estabilidad. Aumentar mejora precisión y aumenta tiempo. */
#define STABILITY_N_R1 801
#define STABILITY_N_RHO 801
#define STABILITY_TOL (-1.0e-10)

/* ------------------------------------------------------------
 * Estructuras
 * ------------------------------------------------------------ */

typedef struct {
    double mh, mH, mA, mHp, sba;
    double lambda6, lambda7, lambda1;
    double tanb;
    int yukawa_type;
    const char *output_filename;
} InputPoint;

typedef struct {
    double lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7;
    double m12sq;
    double alpha;
    double beta;
} LambdaSet;

typedef struct {
    int masses_positive;
    int perturbativity_ok;
    int unitarity_ok;
    int stability_ok;
    int triple_ok;

    double max_abs_lambda;
    double max_abs_unitarity_eigenvalue;
    double min_V4;
    double min_r1;
    double min_rho;
    double min_x;
} ConstraintResult;

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
 * Utilidades matemáticas
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

static double safe_sqrt(double x) {
    if (x <= 0.0) return 0.0;
    return sqrt(x);
}

/* ------------------------------------------------------------
 * Parámetros geométricos
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

static double beta_from_tanb(double tanb) {
    return atan(tanb);
}

static double alpha_from_sba(double beta, double sba) {
    double s = clamp_double(sba, -1.0, 1.0);
    return beta - asin(s);
}

/*
 * Trigonometría estable para tan(beta) muy grande.
 *
 * Evita calcular sin(beta) y cos(beta) desde beta=atan(tanb) cuando tanb ~ 10^7
 * o mayor. En ese régimen cos(beta) es muy pequeño y los errores relativos
 * pueden afectar la reconstrucción de masas.
 *
 * 2HDMC usa:
 *   alpha = -asin(sba) + beta
 *   cba   = sqrt(1-sba^2)
 *
 * Por identidades trigonométricas:
 *   sin(alpha) = sin(beta)cba - cos(beta)sba
 *   cos(alpha) = cos(beta)cba + sin(beta)sba
 */
static AngleSet angles_from_tanb_sba(double tanb, double sba_in) {
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
 * m12^2 desde lambda1
 *
 * 2HDMC usa, para la base física:
 *
 * lambda1 = (mH^2 ca^2 + mh^2 sa^2 - m12^2 tanb)/(v^2 cb^2)
 *           - 3/2 lambda6 tanb + 1/2 lambda7 tanb^3
 *
 * Despejando:
 *
 * m12^2 = [M11 - lambda1 v^2 cb^2
 *              - 3/2 lambda6 v^2 cb^2 tanb
 *              + 1/2 lambda7 v^2 cb^2 tanb^3] / tanb
 *
 * equivalente a la forma con sb,cb:
 * m12^2 = [M11 - lambda1 v^2 cb^2
 *              - 3/2 lambda6 v^2 sb cb
 *              + 1/2 lambda7 v^2 sb^3/cb] / tanb
 * ------------------------------------------------------------ */

static double m12sq_from_lambda1(InputPoint in) {
    AngleSet A = angles_from_tanb_sba(in.tanb, in.sba);

    double tb = in.tanb;
    double v2 = sqr(V_SM);

    double M11 = sqr(in.mH)*sqr(A.ca) + sqr(in.mh)*sqr(A.sa);

    /*
     * Fórmula exacta inversa de THDM::set_param_phys:
     *
     * lambda1 =
     *   (M11 - m12^2 tanb)/(v^2 cb^2)
     *   - 1.5 lambda6 tanb
     *   + 0.5 lambda7 tanb^3
     *
     * Despeje:
     *
     * m12^2 =
     *   [ M11 - v^2 cb^2
     *     (lambda1 + 1.5 lambda6 tanb - 0.5 lambda7 tanb^3)
     *   ] / tanb
     */
    return (
        M11
        - v2*sqr(A.cb)*(in.lambda1 + 1.5*in.lambda6*tb - 0.5*in.lambda7*cube(tb))
    ) / tb;
}

/* ------------------------------------------------------------
 * Reconstrucción de lambdas usando las fórmulas de base física.
 * Estas son las relaciones algebraicas usadas para pasar de
 * (mh,mH,mA,mHp,sba,l6,l7,m12^2,tanb) a lambda_i.
 * ------------------------------------------------------------ */

static LambdaSet reconstruct_lambdas(InputPoint in) {
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

    /*
     * Fórmulas exactas de THDM::set_param_phys, con la única diferencia de que
     * lambda_1 es input y m12^2 se calculó primero desde lambda_1.
     *
     * Esto corrige la versión anterior: ahora lambda_2,...,lambda_5 son
     * exactamente los que 2HDMC construiría en base física si pudiera recibir
     * lambda_1 directamente sin la cancelación numérica del caso tan(beta) >> 1.
     */
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
 * Chequeo inverso de masas desde la base genérica reconstruida.
 *
 * No llama a 2HDMC: reconstruye las masas desde las mismas relaciones
 * algebraicas de la base genérica. Sirve para confirmar que las masas de
 * entrada se preservan antes de probar el punto con CalcGen.
 * ------------------------------------------------------------ */

typedef struct {
    double mh;
    double mH;
    double mA;
    double mHp;
    double sba;
} MassCheck;

static MassCheck masses_from_lambdas(InputPoint in, LambdaSet L) {
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

    /*
     * sba de chequeo. Para eigenvector del estado ligero:
     * tan(2 alpha) = 2 M12/(M11-M22), pero la convención global puede cambiar
     * por signos. Por eso este valor se imprime solo como diagnóstico.
     */
    double alpha_diag = 0.5*atan2(2.0*M12, M11-M22);
    M.sba = sin(A.beta - alpha_diag);

    return M;
}

/* ------------------------------------------------------------
 * Running simple de alpha_s y masas running como en el notebook
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
 * Funciones de lazo complejas
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
 * Anchos
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

    /* mf, Qf, Nc */
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

    /* mf, Qf, I3f, Nc */
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

static WidthResult compute_widths(InputPoint in) {
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
 * Jacobi eigenvalues para matrices simétricas reales pequeñas.
 * Se usa para reproducir la estructura de unitaridad sin GSL.
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
 * Unitaridad árbol: matrices S20, S21, S01, S00.
 * Basado en la construcción usada en 2HDMC, pero autocontenida.
 * ------------------------------------------------------------ */

static double calc_unitarity_max_abs(LambdaSet L) {
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
 * Estabilidad: minimización numérica del potencial cuártico.
 *
 * V4(r1,r2,rho,x) con r1+r2=1, rho in [0,1], x=cos(theta).
 * Para cada (r1,rho), se minimiza analíticamente en x.
 * ------------------------------------------------------------ */

static double V4_min_for_direction(LambdaSet L, double r1, double rho, double *x_best) {
    double r2 = 1.0 - r1;
    double rr = r1*r2;
    double sqrt_rr = safe_sqrt(rr);

    double C = 0.5*L.lambda1*sqr(r1)
             + 0.5*L.lambda2*sqr(r2)
             + L.lambda3*rr
             + L.lambda4*sqr(rho)*rr;

    double A = L.lambda5*sqr(rho)*rr;
    double B = 2.0*rho*sqrt_rr*(L.lambda6*r1 + L.lambda7*r2);

    double best = C - A - B + 2.0*A; /* x=-1 */
    double xb = -1.0;

    double Vp = C - A + B + 2.0*A; /* x=+1 */
    if (Vp < best) {
        best = Vp;
        xb = 1.0;
    }

    if (A > 0.0) {
        double xs = -B/(4.0*A);
        xs = clamp_double(xs, -1.0, 1.0);
        double Vs = C - A + B*xs + 2.0*A*xs*xs;
        if (Vs < best) {
            best = Vs;
            xb = xs;
        }
    }

    if (x_best) *x_best = xb;
    return best;
}

static double calc_stability_min_V4(LambdaSet L, double *r1_best, double *rho_best, double *x_best) {
    double minV = DBL_MAX;
    double br1 = 0.0, brho = 0.0, bx = 0.0;

    for (int i=0;i<STABILITY_N_R1;i++) {
        double r1 = ((double)i)/((double)(STABILITY_N_R1-1));

        for (int j=0;j<STABILITY_N_RHO;j++) {
            double rho = ((double)j)/((double)(STABILITY_N_RHO-1));
            double xb = 0.0;
            double V = V4_min_for_direction(L, r1, rho, &xb);

            if (V < minV) {
                minV = V;
                br1 = r1;
                brho = rho;
                bx = xb;
            }
        }
    }

    if (r1_best) *r1_best = br1;
    if (rho_best) *rho_best = brho;
    if (x_best) *x_best = bx;

    return minV;
}

static ConstraintResult check_constraints(InputPoint in, LambdaSet L) {
    ConstraintResult C;

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

    C.min_V4 = calc_stability_min_V4(L, &C.min_r1, &C.min_rho, &C.min_x);
    C.stability_ok = (isfinite(C.min_V4) && C.min_V4 > STABILITY_TOL);

    C.triple_ok = (C.masses_positive && C.perturbativity_ok && C.unitarity_ok && C.stability_ok);

    return C;
}

/* ------------------------------------------------------------
 * Impresión
 * ------------------------------------------------------------ */

static void print_report(FILE *out, InputPoint in, LambdaSet L, MassCheck Mcheck, ConstraintResult C, WidthResult W) {
    fprintf(out, "============================================================\n");
    fprintf(out, "CalcLambda1Standalone: 2HDM point from lambda_1 input\n");
    fprintf(out, "============================================================\n\n");

    fprintf(out, "Usage convention:\n");
    fprintf(out, "  mh mH mA mHp sin(beta-alpha) lambda_6 lambda_7 lambda_1 tan_beta yukawas_type output_filename\n\n");

    fprintf(out, "INPUT PARAMETERS\n");
    fprintf(out, "------------------------------------------------------------\n");
    fprintf(out, "mh                  = %.18e\n", in.mh);
    fprintf(out, "mH                  = %.18e\n", in.mH);
    fprintf(out, "mA                  = %.18e\n", in.mA);
    fprintf(out, "mHp                 = %.18e\n", in.mHp);
    fprintf(out, "sin(beta-alpha)     = %.18e\n", in.sba);
    fprintf(out, "lambda_6            = %.18e\n", in.lambda6);
    fprintf(out, "lambda_7            = %.18e\n", in.lambda7);
    fprintf(out, "lambda_1 input      = %.18e\n", in.lambda1);
    fprintf(out, "tan(beta)           = %.18e\n", in.tanb);
    fprintf(out, "yukawa_type         = %d\n", in.yukawa_type);
    fprintf(out, "output_filename     = %s\n\n", in.output_filename);

    if (in.yukawa_type != 1) {
        fprintf(out, "WARNING:\n");
        fprintf(out, "  The standalone decay implementation uses the notebook convention\n");
        fprintf(out, "  kappa_f = epsV - epsA with epsV=0 and epsA=1/tan(beta).\n");
        fprintf(out, "  The yukawa_type input is printed for compatibility, but only Type-I-like\n");
        fprintf(out, "  almost-inert decays are implemented here.\n\n");
    }

    if (fabs(in.lambda6) > 0.0 || fabs(in.lambda7) > 0.0) {
        fprintf(out, "Z2 WARNING:\n");
        fprintf(out, "  yukawa_type may respect Z2, but lambda6 or lambda7 is nonzero.\n");
        fprintf(out, "  This is expected for the approximate-Z2 / almost-inert setup.\n\n");
    }

    fprintf(out, "RECONSTRUCTED PARAMETERS\n");
    fprintf(out, "------------------------------------------------------------\n");
    fprintf(out, "alpha               = %.18e\n", L.alpha);
    fprintf(out, "beta                = %.18e\n", L.beta);
    fprintf(out, "m12^2               = %.18e\n", L.m12sq);
    fprintf(out, "lambda_1            = %.18e\n", L.lambda1);
    fprintf(out, "lambda_2            = %.18e\n", L.lambda2);
    fprintf(out, "lambda_3            = %.18e\n", L.lambda3);
    fprintf(out, "lambda_4            = %.18e\n", L.lambda4);
    fprintf(out, "lambda_5            = %.18e\n", L.lambda5);
    fprintf(out, "lambda_6            = %.18e\n", L.lambda6);
    fprintf(out, "lambda_7            = %.18e\n\n", L.lambda7);

    fprintf(out, "MASS CHECK FROM RECONSTRUCTED GENERIC BASIS\n");
    fprintf(out, "------------------------------------------------------------\n");
    fprintf(out, "mh_from_lambdas      = %.18e\n", Mcheck.mh);
    fprintf(out, "mH_from_lambdas      = %.18e\n", Mcheck.mH);
    fprintf(out, "mA_from_lambdas      = %.18e\n", Mcheck.mA);
    fprintf(out, "mHp_from_lambdas     = %.18e\n", Mcheck.mHp);
    fprintf(out, "sba_diag_check       = %.18e\n", Mcheck.sba);
    fprintf(out, "delta_mh             = %.18e\n", Mcheck.mh - in.mh);
    fprintf(out, "delta_mH             = %.18e\n", Mcheck.mH - in.mH);
    fprintf(out, "delta_mA             = %.18e\n", Mcheck.mA - in.mA);
    fprintf(out, "delta_mHp            = %.18e\n\n", Mcheck.mHp - in.mHp);

    fprintf(out, "TRIPLE CHECK AUTOCONTENIDO\n");
    fprintf(out, "------------------------------------------------------------\n");
    fprintf(out, "Masses positive      = %d\n", C.masses_positive ? 1 : 0);
    fprintf(out, "Tree-level unitarity = %d\n", C.unitarity_ok ? 1 : 0);
    fprintf(out, "Perturbativity       = %d\n", C.perturbativity_ok ? 1 : 0);
    fprintf(out, "Stability            = %d\n", C.stability_ok ? 1 : 0);
    fprintf(out, "TRIPLE_OK            = %d\n", C.triple_ok ? 1 : 0);
    fprintf(out, "max |lambda_i|       = %.18e\n", C.max_abs_lambda);
    fprintf(out, "max |unitarity eig|  = %.18e\n", C.max_abs_unitarity_eigenvalue);
    fprintf(out, "min V4               = %.18e\n", C.min_V4);
    fprintf(out, "min V4 r1            = %.18e\n", C.min_r1);
    fprintf(out, "min V4 rho           = %.18e\n", C.min_rho);
    fprintf(out, "min V4 cos(theta)    = %.18e\n\n", C.min_x);

    fprintf(out, "DECAY WIDTHS USED IN THE NOTEBOOK\n");
    fprintf(out, "------------------------------------------------------------\n");
    fprintf(out, "Gamma_bb             = %.18e GeV\n", W.G_bb);
    fprintf(out, "Gamma_cc             = %.18e GeV\n", W.G_cc);
    fprintf(out, "Gamma_tautau         = %.18e GeV\n", W.G_tautau);
    fprintf(out, "Gamma_gg             = %.18e GeV\n", W.G_gg);
    fprintf(out, "Gamma_gammagamma     = %.18e GeV\n", W.G_gaga);
    fprintf(out, "Gamma_Zgamma         = %.18e GeV\n", W.G_Zga);
    fprintf(out, "Gamma_total_selected = %.18e GeV\n\n", W.G_total);

    fprintf(out, "BRANCHING RATIOS\n");
    fprintf(out, "------------------------------------------------------------\n");
    fprintf(out, "BR_bb                = %.18e\n", W.BR_bb);
    fprintf(out, "BR_cc                = %.18e\n", W.BR_cc);
    fprintf(out, "BR_tautau            = %.18e\n", W.BR_tautau);
    fprintf(out, "BR_gg                = %.18e\n", W.BR_gg);
    fprintf(out, "BR_gammagamma        = %.18e\n", W.BR_gaga);
    fprintf(out, "BR_Zgamma            = %.18e\n", W.BR_Zga);
    fprintf(out, "BR_loop              = %.18e\n\n", W.BR_loop);

    fprintf(out, "LIFETIME\n");
    fprintf(out, "------------------------------------------------------------\n");
    fprintf(out, "ctau                 = %.18e m\n\n", W.ctau_m);

    fprintf(out, "MACHINE_READABLE_SUMMARY\n");
    fprintf(out, "------------------------------------------------------------\n");
    fprintf(out, "mh %.18e\n", in.mh);
    fprintf(out, "mH %.18e\n", in.mH);
    fprintf(out, "mA %.18e\n", in.mA);
    fprintf(out, "mHp %.18e\n", in.mHp);
    fprintf(out, "sba %.18e\n", in.sba);
    fprintf(out, "lambda1 %.18e\n", L.lambda1);
    fprintf(out, "lambda2 %.18e\n", L.lambda2);
    fprintf(out, "lambda3 %.18e\n", L.lambda3);
    fprintf(out, "lambda4 %.18e\n", L.lambda4);
    fprintf(out, "lambda5 %.18e\n", L.lambda5);
    fprintf(out, "lambda6 %.18e\n", L.lambda6);
    fprintf(out, "lambda7 %.18e\n", L.lambda7);
    fprintf(out, "m12sq %.18e\n", L.m12sq);
    fprintf(out, "mh_from_lambdas %.18e\n", Mcheck.mh);
    fprintf(out, "mH_from_lambdas %.18e\n", Mcheck.mH);
    fprintf(out, "mA_from_lambdas %.18e\n", Mcheck.mA);
    fprintf(out, "mHp_from_lambdas %.18e\n", Mcheck.mHp);
    fprintf(out, "unitarity %d\n", C.unitarity_ok ? 1 : 0);
    fprintf(out, "perturbativity %d\n", C.perturbativity_ok ? 1 : 0);
    fprintf(out, "stability %d\n", C.stability_ok ? 1 : 0);
    fprintf(out, "triple_ok %d\n", C.triple_ok ? 1 : 0);
    fprintf(out, "Gamma_total %.18e\n", W.G_total);
    fprintf(out, "ctau_m %.18e\n", W.ctau_m);
    fprintf(out, "BR_bb %.18e\n", W.BR_bb);
    fprintf(out, "BR_gg %.18e\n", W.BR_gg);
    fprintf(out, "BR_gammagamma %.18e\n", W.BR_gaga);
    fprintf(out, "BR_Zgamma %.18e\n", W.BR_Zga);
    fprintf(out, "BR_loop %.18e\n", W.BR_loop);
}

/* ------------------------------------------------------------
 * Main
 * ------------------------------------------------------------ */

int main(int argc, char **argv) {
    if (argc < 12) {
        printf("Usage: ./CalcLambda1Standalone mh mH mA mHp sin(beta-alpha) lambda_6 lambda_7 lambda_1 tan_beta yukawas_type output_filename\n");
        return 1;
    }

    InputPoint in;
    in.mh = atof(argv[1]);
    in.mH = atof(argv[2]);
    in.mA = atof(argv[3]);
    in.mHp = atof(argv[4]);
    in.sba = atof(argv[5]);
    in.lambda6 = atof(argv[6]);
    in.lambda7 = atof(argv[7]);
    in.lambda1 = atof(argv[8]);
    in.tanb = atof(argv[9]);
    in.yukawa_type = (int)atof(argv[10]);
    in.output_filename = argv[11];

    if (in.tanb <= 0.0 || fabs(in.sba) > 1.0) {
        fprintf(stderr, "ERROR: invalid tan(beta) or sin(beta-alpha).\n");
        return 2;
    }

    LambdaSet L = reconstruct_lambdas(in);
    MassCheck Mcheck = masses_from_lambdas(in, L);
    ConstraintResult C = check_constraints(in, L);
    WidthResult W = compute_widths(in);

    print_report(stdout, in, L, Mcheck, C, W);

    FILE *fout = fopen(in.output_filename, "w");
    if (!fout) {
        fprintf(stderr, "ERROR: cannot open output file %s\n", in.output_filename);
        return 3;
    }

    print_report(fout, in, L, Mcheck, C, W);
    fclose(fout);

    return C.triple_ok ? 0 : 10;
}
