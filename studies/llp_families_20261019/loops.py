"""Independent LO loop widths for phi -> gamma gamma and phi -> Z gamma (cross-check of the 2HDMC values)."""
import math

# ---- independent LO phi->gamma gamma (task formula; G_F = 1/(sqrt2 v^2), alpha(0)) ----
def f_loop(tau):
    if tau <= 1:
        return math.asin(math.sqrt(tau)) ** 2
    s = math.sqrt(1 - 1 / tau)
    return -0.25 * (math.log((1 + s) / (1 - s)) - 1j * math.pi) ** 2


def A_half(tau):
    return 2 * (tau + (tau - 1) * f_loop(tau)) / tau**2


def A_zero(tau):
    return -(tau - f_loop(tau)) / tau**2


def gamma_gaga_LO(m, tanb, Z7, mt=172.5, mb=4.18, mc=1.27, mtau=1.77686, mhp=None, v=246.0, alpha=1 / 137.035999):
    mhp = mhp or m + 100
    ferm = sum(Nc * Q * Q * A_half(m * m / (4 * mf * mf)) for Nc, Q, mf in
               ((3, 2 / 3, mt), (3, -1 / 3, mb), (3, 2 / 3, mc), (1, -1, mtau)))
    a_f = -ferm / tanb
    a_h = -(v * v * Z7 / (2 * mhp * mhp)) * A_zero(m * m / (4 * mhp * mhp))
    GF = 1 / (math.sqrt(2) * v * v)
    pref = GF * alpha**2 * m**3 / (128 * math.sqrt(2) * math.pi**3)
    return pref * abs(a_f + a_h) ** 2, pref * abs(a_f) ** 2, pref * abs(a_h) ** 2



# 5) independent LO phi -> Z gamma.  Fermion and W parts: Djouadi, Anatomy I (hep-ph/0503172), validated on the SM h->Z gamma.
#    H+- part: normalised by the decoupling (low-energy) limit relative to the gamma-gamma formula above: a heavy H+- must
#    give Zgamma/gammagamma = g_Z(H+H-)/e = (2c_W^2-1)/(2 s_W c_W), and a heavy fermion (I3-2Q s_W^2)/(2 s_W c_W Q).  With
#    I1 -> 1/6, I1-I2 -> -1/3, A0 -> 1/3, A1/2 -> 4/3 this fixes the H+- coefficient to +(2c_W^2-1)/c_W * v^2 Z7/(2 m_H+^2).
#    (A first transcription with -(2c_W^2-1) - as printed in Anatomy II eqs. 2.23/2.33 - is off by 1/c_W and a sign; the
#    2HDMC source documents the same c_W factor relative to Anatomy II and follows HDECAY.)  s_W^2 is varied between
#    0.2312 (MSbar) and 1-M_W^2/M_Z^2 (on-shell) as a scheme probe.
SW2, MZ, MW = 0.2312, 91.1876, 80.379


def f_inv(x):  # f evaluated at 1/x  (x = 4 m^2 / M^2)
    return f_loop(1 / x)


def g_inv(x):
    t = 1 / x
    if t <= 1:
        return math.sqrt(1 / t - 1) * math.asin(math.sqrt(t))
    s = math.sqrt(1 - 1 / t)
    return s / 2 * (math.log((1 + s) / (1 - s)) - 1j * math.pi)


def I1(t, l):
    return t * l / (2 * (t - l)) + t * t * l * l / (2 * (t - l) ** 2) * (f_inv(t) - f_inv(l)) + t * t * l / (t - l) ** 2 * (g_inv(t) - g_inv(l))


def I2(t, l):
    return -t * l / (2 * (t - l)) * (f_inv(t) - f_inv(l))


def gamma_zga_LO(m, kf, Z7=0.0, mhp=None, with_W=False, sw2=SW2, v=246.0, alpha=1 / 137.035999, mt=172.5, mb=4.18, mc=1.27, mtau=1.77686):
    SW2_ = sw2
    cw = math.sqrt(1 - SW2_)
    amp = 0
    for Nc, Q, I3, mf in ((3, 2 / 3, 0.5, mt), (3, -1 / 3, -0.5, mb), (3, 2 / 3, 0.5, mc), (1, -1, -0.5, mtau)):
        tau, lam = 4 * mf * mf / (m * m), 4 * mf * mf / (MZ * MZ)
        amp += Nc * Q * (2 * I3 - 4 * Q * SW2_) / cw * kf * (I1(tau, lam) - I2(tau, lam))
    if with_W:
        tau, lam = 4 * MW * MW / (m * m), 4 * MW * MW / (MZ * MZ)
        amp += cw * (4 * (3 - SW2_ / cw**2) * I2(tau, lam) + ((1 + 2 / tau) * SW2_ / cw**2 - (5 + 2 / tau)) * I1(tau, lam))
    if Z7:
        mhp = mhp or m + 100
        tau, lam = 4 * mhp * mhp / (m * m), 4 * mhp * mhp / (MZ * MZ)
        amp += (2 * cw * cw - 1) / cw * v * v * Z7 / (2 * mhp * mhp) * I1(tau, lam)
    GF = 1 / (math.sqrt(2) * v * v)
    return GF**2 * MW**2 * alpha * m**3 / (64 * math.pi**4) * (1 - MZ**2 / m**2) ** 3 * abs(amp) ** 2


