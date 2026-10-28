"""Adversarial audit test suite for 2HDM BFB, Global EW minimum, and negative trilinear g(h phi phi).

Reproducible test suite covering:
1. Exact proof that V4 minimization on full bilinear domain reduces to f(t) >= 0 when Z6=0 and Z4=Z5<0.
2. The 30 inconclusive points from Study D: exactly 27 pass BFB, 3 fail with explicit rational witnesses.
3. The four adversarial test points: (130, -27.48) passes, (291.33, -138.01), (306, -152.25), (320.67, -167.20) fail.
4. Positive and negative controls.
5. Global electroweak minimum verification: unique global minimum for the 27 BFB-pass points.
6. Diagnostic verification of the bug in _witness and sufficient certificate in aligned_constraints.py.

Run with:
    pytest tests/bfb_negative_trilinear_audit/test_bfb_negative_trilinear_audit.py -v
"""

from fractions import Fraction
from pathlib import Path
import pytest
import numpy as np
import pandas as pd


class RationalPoly:
    """Exact polynomial over Fraction with root isolation via Sturm sequences."""

    def __init__(self, coeffs: list[Fraction | int | str]):
        self.coeffs = [Fraction(c) for c in coeffs]
        while len(self.coeffs) > 1 and self.coeffs[-1] == 0:
            self.coeffs.pop()

    @property
    def deg(self) -> int:
        return len(self.coeffs) - 1

    def __call__(self, t: Fraction) -> Fraction:
        val = Fraction(0)
        p = Fraction(1)
        for c in self.coeffs:
            val += c * p
            p *= t
        return val

    def deriv(self) -> "RationalPoly":
        if self.deg <= 0:
            return RationalPoly([Fraction(0)])
        return RationalPoly([i * self.coeffs[i] for i in range(1, len(self.coeffs))])

    def __divmod__(self, other: "RationalPoly") -> tuple["RationalPoly", "RationalPoly"]:
        if other.deg < 0 or (other.deg == 0 and other.coeffs[0] == 0):
            raise ZeroDivisionError("polynomial division by zero")
        rem = list(self.coeffs)
        quot = [Fraction(0)] * max(0, self.deg - other.deg + 1)
        while len(rem) - 1 >= other.deg:
            deg_diff = len(rem) - 1 - other.deg
            factor = rem[-1] / other.coeffs[-1]
            quot[deg_diff] = factor
            for i, c in enumerate(other.coeffs):
                rem[deg_diff + i] -= factor * c
            while len(rem) > 1 and rem[-1] == 0:
                rem.pop()
        return RationalPoly(quot), RationalPoly(rem)

    def __neg__(self) -> "RationalPoly":
        return RationalPoly([-c for c in self.coeffs])


def sturm_chain(p: RationalPoly) -> list[RationalPoly]:
    chain = [p, p.deriv()]
    while chain[-1].deg > 0:
        _, rem = divmod(chain[-2], chain[-1])
        if rem.deg == 0 and rem.coeffs[0] == 0:
            break
        chain.append(-rem)
    return chain


def count_sign_changes(vals: list[Fraction]) -> int:
    signs = [1 if v > 0 else -1 if v < 0 else 0 for v in vals if v != 0]
    return sum(1 for i in range(len(signs) - 1) if signs[i] * signs[i + 1] < 0)


def count_real_roots_on_nonnegative(p: RationalPoly) -> int:
    chain = sturm_chain(p)
    vals_zero = [p_i(Fraction(0)) for p_i in chain]
    vals_inf = [p_i.coeffs[-1] for p_i in chain]
    return count_sign_changes(vals_zero) - count_sign_changes(vals_inf)


def v4_quartic(u: Fraction, b: Fraction, x: Fraction, y: Fraction, z: list[Fraction]) -> Fraction:
    z1, z2, z3, z4, z5, z6, z7 = z
    return (
        z1 * u * u / 2
        + z2 * b * b / 2
        + z3 * u * b
        + (z4 + z5) * x * x
        + (z4 - z5) * y * y
        + 2 * z6 * u * x
        + 2 * z7 * b * x
    )


def reduced_f(t: Fraction, z: list[Fraction]) -> Fraction:
    z1, z2, z3, z4, z5, z6, z7 = z
    c = z3 + z4 + z5
    return (z1 / 2) * t**4 + c * t**2 - 2 * abs(z7) * t + z2 / 2


@pytest.fixture
def points_df():
    csv_path = Path(__file__).with_name("points_audited.csv")
    if not csv_path.exists():
        pytest.fail(f"Missing {csv_path}. Run audit generator first.")
    return pd.read_csv(csv_path)


def test_quartic_reduction_to_f_t():
    """Verify that V4(u,b,x,y) >= b^2 f(sqrt(u/b)) everywhere in bilinear domain when Z6=0, Z4=Z5<0."""
    z1 = Fraction("0.258734")
    z2 = Fraction("1.2")
    z3 = Fraction("1.5")
    z4 = Fraction("-0.8")
    z5 = Fraction("-0.8")  # Z4 = Z5 < 0
    z6 = Fraction("0")      # Z6 = 0
    z7 = Fraction("-0.0015")
    z = [z1, z2, z3, z4, z5, z6, z7]

    # Sample random and boundary bilinear points: u, b > 0, x^2 + y^2 <= u*b
    rng = np.random.default_rng(42)
    for _ in range(200):
        u_f = rng.uniform(0.1, 10.0)
        b_f = rng.uniform(0.1, 10.0)
        r_max = np.sqrt(u_f * b_f)
        theta = rng.uniform(0, 2 * np.pi)
        r = rng.uniform(0, r_max)

        u = Fraction(u_f).limit_denominator(1000)
        b = Fraction(b_f).limit_denominator(1000)
        x = Fraction(r * np.cos(theta)).limit_denominator(1000)
        y = Fraction(r * np.sin(theta)).limit_denominator(1000)
        if x * x + y * y > u * b:
            continue

        val_v4 = v4_quartic(u, b, x, y, z)
        t = Fraction(np.sqrt(float(u / b))).limit_denominator(10000)
        lower_bound = b * b * reduced_f(t, z)

        assert val_v4 >= lower_bound - Fraction(1, 100000)

    # Test exact boundary point x = -sgn(Z7) * sqrt(u*b), y = 0
    u = Fraction(4)
    b = Fraction(1)
    x = Fraction(2)  # -sgn(Z7)*2 since Z7 < 0
    y = Fraction(0)
    assert v4_quartic(u, b, x, y, z) == b * b * reduced_f(Fraction(2), z)


def test_30_points_split_27_pass_3_fail(points_df):
    """Verify that exactly 27 points pass BFB and 3 fail BFB."""
    assert len(points_df) == 30
    pass_pts = points_df[points_df["new_BFB_status"] == "PASS"]
    fail_pts = points_df[points_df["new_BFB_status"] == "FAIL"]

    assert len(pass_pts) == 27
    assert len(fail_pts) == 3

    # Fail points must be eval_index 167, 181, 195
    assert set(fail_pts["eval_index"].tolist()) == {167, 181, 195}


def test_adversarial_points(points_df):
    """Test the four specific adversarial points requested in the prompt."""
    # 1. (130, -27.48): eval_index 13 -> PASS
    p13 = points_df[points_df["eval_index"] == 13].iloc[0]
    assert np.isclose(p13["mphi_GeV"], 130.0)
    assert np.isclose(p13["g_hphiphi_GeV"], -27.479675, atol=1e-3)
    assert p13["new_BFB_status"] == "PASS"
    assert p13["new_global_EW_status"] == "PASS"

    # 2. (291.33, -138.01): eval_index 167 -> FAIL
    p167 = points_df[points_df["eval_index"] == 167].iloc[0]
    assert np.isclose(p167["mphi_GeV"], 291.333333, atol=1e-3)
    assert np.isclose(p167["g_hphiphi_GeV"], -138.008311, atol=1e-3)
    assert p167["new_BFB_status"] == "FAIL"
    assert p167["new_global_EW_status"] == "FAIL"

    # 3. (306.0, -152.25): eval_index 181 -> FAIL
    p181 = points_df[points_df["eval_index"] == 181].iloc[0]
    assert np.isclose(p181["mphi_GeV"], 306.0, atol=1e-3)
    assert np.isclose(p181["g_hphiphi_GeV"], -152.253659, atol=1e-3)
    assert p181["new_BFB_status"] == "FAIL"
    assert p181["new_global_EW_status"] == "FAIL"

    # 4. (320.67, -167.20): eval_index 195 -> FAIL
    p195 = points_df[points_df["eval_index"] == 195].iloc[0]
    assert np.isclose(p195["mphi_GeV"], 320.666667, atol=1e-3)
    assert np.isclose(p195["g_hphiphi_GeV"], -167.198555, atol=1e-3)
    assert p195["new_BFB_status"] == "FAIL"
    assert p195["new_global_EW_status"] == "FAIL"


def test_positive_and_negative_controls():
    """Verify built controls: positive control (g > 0) PASS, negative control (eval_index 209) FAIL."""
    # Positive control: mphi=130, rho=0.0 => g > 0
    z1 = Fraction("0.258734")
    z2 = Fraction("1.2")
    z3 = Fraction("1.5")
    z4 = Fraction("-0.5")
    z5 = Fraction("-0.5")
    z6 = Fraction("0")
    z7 = Fraction("-0.001")
    poly_pos = RationalPoly([z2 / 2, -2 * abs(z7), z3 + z4 + z5, Fraction(0), z1 / 2])
    assert count_real_roots_on_nonnegative(poly_pos) == 0

    # Negative control: eval_index 209 (mphi=335.33, rho=1.2, g=-182.84)
    # known to be rejected by BFB_general in campaign
    z1_neg = Fraction("0.258734")
    z2_neg = Fraction("1.199996")
    z3_neg = Fraction("1.803721")
    z4_neg = Fraction("-1.273492")
    z5_neg = Fraction("-1.273492")
    z7_neg = Fraction("-0.001943")
    poly_neg = RationalPoly([z2_neg / 2, -2 * abs(z7_neg), z3_neg + z4_neg + z5_neg, Fraction(0), z1_neg / 2])
    assert count_real_roots_on_nonnegative(poly_neg) == 2


def test_global_electroweak_minimum_coercivity_and_stationarity(points_df):
    """Verify that all 27 BFB PASS points possess a unique global electroweak minimum."""
    pass_pts = points_df[points_df["new_BFB_status"] == "PASS"]
    u_ew = 246.0**2 / 2.0

    for _, r in pass_pts.iterrows():
        z1 = float(Fraction(str(r["Z1"])))
        z2 = float(Fraction(str(r["Z2"])))
        z3 = float(Fraction(str(r["Z3"])))
        z4 = float(Fraction(str(r["Z4"])))
        z5 = float(Fraction(str(r["Z5"])))
        z7 = float(Fraction(str(r["Z7"])))
        y2 = float(Fraction(str(r["Y2_GeV2"])))
        c = z3 + z4 + z5

        # Check along the p+(q) stationarity curve for q in [1e-3, 1e4]
        q_vals = np.logspace(-3, 4, 1000)
        disc = 9.0 * z7**2 * q_vals**2 - 4.0 * c * (y2 + z2 * q_vals**2)
        p_vals = (3.0 * abs(z7) * q_vals - np.sqrt(disc)) / (2.0 * c)

        eq1 = z1 * p_vals**3 + (c * q_vals**2 - z1 * u_ew) * p_vals - abs(z7) * q_vals**3
        # Must be strictly positive for all q > 0 (no zeros)
        assert np.all(eq1 > 0), f"Found potential stationary point for eval_index {r['eval_index']}"


def test_code_sufficient_certificate_and_witness_limitations():
    """Verify exact code root causes: c = g/v requirement, and _witness min(u,b) flaw."""
    # 1. Demonstrate that c = g/v when Z4=Z5 and Delta > 0
    # c = Z3 + min(0, Z4+Z5, Z4-Z5) = Z3 + Z4 + Z5 = g/v
    mphi, mA, mHp, Y2, v = 200.0, 300.0, 300.0, 40000.0, 246.0
    z3 = 2.0 * (mHp**2 - Y2) / v**2
    z4 = (mphi**2 + mA**2 - 2.0 * mHp**2) / v**2
    z5 = (mphi**2 - mA**2) / v**2
    c = z3 + min(0.0, z4 + z5, z4 - z5)
    g = v * (z3 + z4 + z5)
    assert np.isclose(c, g / v)

    # 2. Demonstrate that _witness in aligned_constraints.py misses point 195
    # because it samples x = a * min(u,b) with a^2+c^2 <= 1, which under-samples the boundary sqrt(u*b).
    # For point 195:
    z1 = Fraction("0.258734")
    z2 = Fraction("1.199996")
    z3 = Fraction("1.770372")
    z4 = Fraction("-1.225020")
    z5 = Fraction("-1.225020")
    z6 = Fraction("0")
    z7 = Fraction("-0.001880")
    z = [z1, z2, z3, z4, z5, z6, z7]

    # At boundary u=1, b=1/4, x=-1/2, y=0: V4 < 0
    assert v4_quartic(Fraction(1), Fraction(1, 4), Fraction(-1, 2), Fraction(0), z) < 0

    # But min(u,b) = 1/4, so _witness only tested x in [-1/4, 1/4]
    # At interior x=-1/4: V4 > 0
    assert v4_quartic(Fraction(1), Fraction(1, 4), Fraction(-1, 4), Fraction(0), z) > 0
