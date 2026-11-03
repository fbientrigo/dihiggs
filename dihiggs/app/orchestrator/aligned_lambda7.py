"""Stable arbitrary-precision T1 constructor for the aligned lambda7 family.

This module constructs a tree-level family point.  It does not certify BFB,
unitarity, global-minimum, or phenomenological constraints.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, localcontext
from typing import Any


@dataclass(frozen=True)
class FamilyPoint:
    """String-preserving result; ``values`` contains only decimal strings."""

    inputs: dict[str, str]
    trajectory: str
    values: dict[str, str]
    constraints: dict[str, dict[str, Any]]
    provenance: dict[str, Any]

    def as_dict(self) -> dict[str, Any]:
        return {
            "inputs": dict(self.inputs),
            "trajectory": self.trajectory,
            "values": dict(self.values),
            "constraints": self.constraints,
            "provenance": self.provenance,
        }


def _decimal(name: str, value: Any, *, positive: bool = False) -> Decimal:
    try:
        out = value if isinstance(value, Decimal) else Decimal(str(value))
    except (InvalidOperation, ValueError, TypeError) as exc:
        raise ValueError(f"{name} must be a finite decimal number") from exc
    if not out.is_finite() or (positive and out <= 0):
        qualifier = "finite and > 0" if positive else "finite"
        raise ValueError(f"{name} must be {qualifier}")
    return out


def _lexeme(value: Decimal) -> str:
    return str(value)


class _Jet2:
    """Small second-order forward-mode jet for an independent potential audit."""

    def __init__(self, value: Decimal, grad: list[Decimal], hess: list[list[Decimal]]):
        self.value, self.grad, self.hess = value, grad, hess

    @classmethod
    def variable(cls, value: Decimal, index: int, size: int):
        zero = Decimal(0)
        grad = [zero for _ in range(size)]
        grad[index] = Decimal(1)
        return cls(value, grad, [[zero for _ in range(size)] for _ in range(size)])

    @classmethod
    def constant(cls, value: Decimal, size: int):
        zero = Decimal(0)
        return cls(value, [zero for _ in range(size)],
                   [[zero for _ in range(size)] for _ in range(size)])

    def __add__(self, other):
        if not isinstance(other, _Jet2):
            other = _Jet2.constant(Decimal(other), len(self.grad))
        return _Jet2(self.value + other.value,
                     [a + b for a, b in zip(self.grad, other.grad)],
                     [[self.hess[i][j] + other.hess[i][j] for j in range(len(self.grad))]
                      for i in range(len(self.grad))])

    __radd__ = __add__

    def __neg__(self):
        return self * Decimal(-1)

    def __sub__(self, other):
        return self + (-other if isinstance(other, _Jet2) else -Decimal(other))

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        if not isinstance(other, _Jet2):
            other = _Jet2.constant(Decimal(other), len(self.grad))
        n = len(self.grad)
        return _Jet2(self.value * other.value,
                     [self.grad[i] * other.value + other.grad[i] * self.value for i in range(n)],
                     [[self.hess[i][j] * other.value + other.hess[i][j] * self.value
                       + self.grad[i] * other.grad[j] + other.grad[i] * self.grad[j]
                       for j in range(n)] for i in range(n)])

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, _Jet2):
            raise TypeError("_Jet2 only supports division by scalars")
        return self * (Decimal(1) / Decimal(other))


def _potential_closure(inputs: dict[str, str], values: dict[str, Decimal], precision: int):
    """Differentiate the original six-field potential and rotate its Hessians."""
    with localcontext() as ctx:
        ctx.prec = precision
        D = Decimal
        t, L = D(inputs["tan_beta"]), D(inputs["lambda1_ref"])
        mh, mp, ma, mhp, v = (D(inputs[key]) for key in
            ("mh_GeV", "mphi_GeV", "mA_GeV", "mHp_GeV", "v_GeV"))
        r = D(1) / t
        norm = (D(1) + r * r).sqrt()
        c, s = r / norm, D(1) / norm
        w1, w2 = v * c, v * s
        q11, q22, q12 = (values[key] for key in
            ("m11_sq_GeV2", "m22_sq_GeV2", "m12_sq_GeV2"))
        l1, l2, l3, l4, l5, l6, l7 = (values[f"lambda{i}"] for i in range(1, 8))
        size = 6
        f = [_Jet2.variable(x, i, size) for i, x in enumerate((w1, w2, D(0), D(0), D(0), D(0)))]
        w1j, w2j, e1, e2, x, y = f
        A = (w1j * w1j + e1 * e1 + x * x) * D("0.5")
        B = (w2j * w2j + e2 * e2 + y * y) * D("0.5")
        C = (w1j * w2j + e1 * e2 + x * y) * D("0.5")
        Q = (w1j * e2 - e1 * w2j) * D("0.5")
        V = (q11 * A + q22 * B - D(2) * q12 * C + l1 * A * A / D(2)
             + l2 * B * B / D(2) + l3 * A * B + l4 * (C * C + Q * Q)
             + l5 * (C * C - Q * Q) + D(2) * (l6 * A + l7 * B) * C)
        scale = max(v * v, mh * mh, mp * mp, ma * ma, mhp * mhp)
        tadpole_scale = v * scale
        grad_res = max(abs(item) for item in V.grad) / tadpole_scale
        even = [[V.hess[i][j] for j in (0, 1)] for i in (0, 1)]
        rot = ((c, s), (s, -c))
        rotated = [[sum(rot[i][a] * even[a][b] * rot[j][b]
                        for a in range(2) for b in range(2)) for j in range(2)] for i in range(2)]
        def block_res(start: int, mass: Decimal):
            target = ((mass * mass * s * s, -mass * mass * s * c),
                      (-mass * mass * s * c, mass * mass * c * c))
            residual = max(abs(V.hess[start+i][start+j] - target[i][j])
                           for i in range(2) for j in range(2)) / scale
            return residual
        even_target = ((mh * mh, D(0)), (D(0), mp * mp))
        even_res = max(abs(rotated[i][j] - even_target[i][j]) for i in range(2) for j in range(2)) / scale
        odd_res = block_res(2, ma)
        charged_res = block_res(4, mhp)
        def symmetric_2x2_eigenvalues(matrix):
            aa, ab, ba, dd = matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1]
            disc = ((aa-dd)*(aa-dd) + D(4)*ab*ba).sqrt()
            high = (aa+dd+disc)/D(2)
            low = ((aa*dd-ab*ba)/high) if high != 0 else (aa+dd-disc)/D(2)
            return (min(low,high), max(low,high))
        even_eigs = symmetric_2x2_eigenvalues(rotated)
        odd_eigs = symmetric_2x2_eigenvalues([[V.hess[i][j] for j in (2,3)] for i in (2,3)])
        charged_eigs = symmetric_2x2_eigenvalues([[V.hess[i][j] for j in (4,5)] for i in (4,5)])
        physical_eigs = [*even_eigs, odd_eigs[1], charged_eigs[1]]
        # Rotate the quartic polynomial independently by evaluating it in
        # Higgs-basis bilinears U=H1†H1, B=H2†H2, ReC and ImC. Symmetric
        # differences isolate all seven coefficients.
        one, two, four = D(1), D(2), D(4)
        zero = D(0)
        def quartic(U: Decimal, BB: Decimal, R: Decimal, I: Decimal) -> Decimal:
            A = c*c*U + s*s*BB - two*s*c*R
            B0 = s*s*U + c*c*BB + two*s*c*R
            Cr = s*c*(U-BB) + (c*c-s*s)*R
            return (l1*A*A/two + l2*B0*B0/two + l3*A*B0
                    + l4*(Cr*Cr+I*I) + l5*(Cr*Cr-I*I)
                    + two*(l6*A+l7*B0)*Cr)
        qU, qB = quartic(one, zero, zero, zero), quartic(zero, one, zero, zero)
        z1_ind, z2_ind = two*qU, two*qB
        z3_ind = quartic(one, one, zero, zero) - qU - qB
        z6_ind = (quartic(one, zero, one, zero)-quartic(one, zero,-one,zero))/four
        z7_ind = (quartic(zero, one, one, zero)-quartic(zero, one,-one,zero))/four
        z45_sum = (quartic(one,zero,one,zero)+quartic(one,zero,-one,zero))/two-qU
        z45_diff = quartic(zero,zero,zero,one)
        z4_ind, z5_ind = (z45_sum+z45_diff)/two, (z45_sum-z45_diff)/two
        z_ind = (z1_ind,z2_ind,z3_ind,z4_ind,z5_ind,z6_ind,z7_ind)
        z_targets = tuple(values[f"Z{i}"] for i in range(1,8))
        z_abs = max(abs(z_ind[i]-z_targets[i]) for i in range(7))
        z_norm = z_abs / max(one, *(abs(x) for x in z_ind), *(abs(x) for x in z_targets))
        z6_norm = abs(z6_ind) / max(one, *(abs(x) for x in z_ind))
        cubic_hphi = v * (z3_ind + z4_ind + z5_ind)
        cubic_phi_hp = -v * z7_ind
        cubic_hphi_target = two * mp * mp / v
        cubic_phi_hp_target = -v * values["Z7"]
        cubic_abs = max(abs(cubic_hphi-cubic_hphi_target), abs(cubic_phi_hp-cubic_phi_hp_target))
        cubic_norm = cubic_abs / max(one,abs(cubic_hphi_target),abs(cubic_phi_hp_target))
        hessian_norm = max(grad_res, even_res, odd_res, charged_res)
        tolerance = D(10) ** D(-(precision - 12))
        alignment_abs_tolerance = tolerance
        return {
            "method": "second-order automatic differentiation of generic six-field appendix potential; independent basis rotation",
            "precision_digits": precision,
            "tadpole_max_abs_normalized": str(+grad_res),
            "even_hessian_max_abs_normalized": str(+even_res),
            "CP_odd_hessian_max_abs_normalized": str(+odd_res),
            "charged_hessian_max_abs_normalized": str(+charged_res),
            "even_hessian_offdiagonal_abs_GeV2": str(+abs(rotated[0][1])),
            "independent_Z1_to_Z7": [str(+x) for x in z_ind],
            "Z1_to_Z7_max_abs_residual": str(+z_abs),
            "Z1_to_Z7_max_abs_normalized_residual": str(+z_norm),
            "Z6_independently_reconstructed": str(+z6_ind),
            "Z6_abs_normalized": str(+z6_norm),
            "independent_cubic_h_phi_phi_GeV": str(+cubic_hphi),
            "independent_cubic_phi_Hp_Hm_GeV": str(+cubic_phi_hp),
            "cubic_max_abs_residual_GeV": str(+cubic_abs),
            "cubic_max_abs_normalized_residual": str(+cubic_norm),
            "CP_even_degenerate": bool(mh == mp),
            "physical_hessian_eigenvalues_GeV2": [str(+x) for x in physical_eigs],
            "smallest_physical_hessian_eigenvalue_GeV2": str(+min(physical_eigs)),
            "CP_even_hessian_eigenvalues_GeV2": [str(+x) for x in even_eigs],
            "CP_odd_non_goldstone_eigenvalue_GeV2": str(+odd_eigs[1]),
            "charged_non_goldstone_eigenvalue_GeV2": str(+charged_eigs[1]),
            "alignment_abs_zero_tolerance": str(alignment_abs_tolerance),
            "normalized_hessian_residual": str(+hessian_norm),
            "declared_normalization": {"Hessian": "max(v^2,mh^2,mphi^2,mA^2,mHp^2) in GeV^2",
                                       "tadpoles": "v*max(v^2,mh^2,mphi^2,mA^2,mHp^2) in GeV^3"},
            "tolerance": str(tolerance),
            "pass": max(hessian_norm,z_norm,z6_norm,cubic_norm) <= tolerance and abs(z6_ind) <= alignment_abs_tolerance,
            "alignment_pass": abs(z6_ind) <= alignment_abs_tolerance and abs(rotated[0][1]) <= scale*alignment_abs_tolerance,
        }


def build_aligned_lambda7_family(
    *,
    tan_beta: Any,
    lambda1_ref: Any,
    mh_GeV: Any,
    mphi_GeV: Any,
    mA_GeV: Any,
    mHp_GeV: Any,
    v_GeV: Any,
    lambda6: Any | None = None,
    X_lambda6_tanbeta: Any | None = None,
    precision: int = 80,
) -> FamilyPoint:
    """Build the aligned family, fixing exactly one of lambda6 or X=lambda6*t.

    ``precision`` is the Decimal context precision in significant digits.
    Inputs should be decimal strings when their source has more than binary64
    precision.  The separate ``Xstar`` trajectory is named by the caller and
    can be constructed with ``X_lambda6_tanbeta=(b-L-k/t**2)/(1+1/t**2)``.
    """
    if (lambda6 is None) == (X_lambda6_tanbeta is None):
        raise ValueError("supply exactly one of lambda6 and X_lambda6_tanbeta")
    if precision < 18:
        raise ValueError("precision must be at least 18 decimal digits")
    t = _decimal("tan_beta", tan_beta, positive=True)
    L = _decimal("lambda1_ref", lambda1_ref)
    mh = _decimal("mh_GeV", mh_GeV, positive=True)
    mp = _decimal("mphi_GeV", mphi_GeV, positive=True)
    ma = _decimal("mA_GeV", mA_GeV, positive=True)
    mhp = _decimal("mHp_GeV", mHp_GeV, positive=True)
    v = _decimal("v_GeV", v_GeV, positive=True)
    with localcontext() as ctx:
        ctx.prec = precision
        one, two, three, four = Decimal(1), Decimal(2), Decimal(3), Decimal(4)
        r = one / t
        a = mh * mh / (v * v)
        d = (mp * mp - mh * mh / two) / (v * v)
        k = L - a
        b = two * mp * mp / (v * v)
        if lambda6 is not None:
            l6 = _decimal("lambda6", lambda6)
            X = l6 * t
            trajectory = "fixed_lambda6"
        else:
            X = _decimal("X_lambda6_tanbeta", X_lambda6_tanbeta)
            l6 = X / t
            trajectory = "fixed_X"

        l7 = -two * d * r + three * l6 * r**2 + two * k * r**3
        l2 = a + four * d * r**2 - four * l6 * r**3 - three * k * r**4
        common = two * l6 * r + k * r**2
        l3 = two * mhp * mhp / (v * v) - common
        l4 = (mp * mp + ma * ma - two * mhp * mhp) / (v * v) - common
        l5 = (mp * mp - ma * ma) / (v * v) - common
        z1 = a
        z2 = L - four * l6 * r + two * (b - L) * r**2 - four * l6 * r**3 - three * k * r**4
        z3 = two * mhp * mhp / (v * v)
        z4 = (mp * mp + ma * ma - two * mhp * mhp) / (v * v)
        z5 = (mp * mp - ma * ma) / (v * v)
        z6 = Decimal(0)
        z7 = -l6 * (one + r**2) + (b - L) * r - k * r**3

        # Reciprocal-based trigonometry avoids atan(t) and cancellation at large t.
        norm = (one + r * r).sqrt()
        c_beta, s_beta = r / norm, one / norm
        M_sq = mh * mh / two
        m12_sq = M_sq * s_beta * c_beta
        m11_sq = -mh * mh * c_beta * c_beta / two
        m22_sq = -mh * mh * s_beta * s_beta / two
        Y1, Y2, Y3 = -mh * mh / two, Decimal(0), Decimal(0)
        ghpp = two * mp * mp / v
        cph = -v * z7
        kap_f, kap_v = -r, Decimal(0)
        values = {
            "lambda1": L, "lambda2": l2, "lambda3": l3, "lambda4": l4,
            "lambda5": l5, "lambda6": l6, "lambda7": l7,
            "X_lambda6_tanbeta": X,
            "Z1": z1, "Z2": z2, "Z3": z3, "Z4": z4, "Z5": z5,
            "Z6": z6, "Z7": z7,
            "M_sq_GeV2": M_sq, "m12_sq_GeV2": m12_sq,
            "m11_sq_GeV2": m11_sq, "m22_sq_GeV2": m22_sq,
            "Y1_GeV2": Y1, "Y2_GeV2": Y2, "Y3_GeV2": Y3,
            "c_beta": c_beta, "s_beta": s_beta,
            "h_mass_reconstructed_GeV": mh, "phi_mass_reconstructed_GeV": mp,
            "A_mass_reconstructed_GeV": ma, "Hp_mass_reconstructed_GeV": mhp,
            "g_hphiphi_GeV": ghpp, "CphiHpHm_potential_GeV": cph,
            "kappa_phi_fermions": kap_f, "kappa_phi_V": kap_v,
        }
        lexical = {key: _lexeme(val) for key, val in values.items()}
        inputs = {
            "tan_beta": _lexeme(t), "lambda1_ref": _lexeme(L),
            "mh_GeV": _lexeme(mh), "mphi_GeV": _lexeme(mp),
            "mA_GeV": _lexeme(ma), "mHp_GeV": _lexeme(mhp), "v_GeV": _lexeme(v),
            "lambda6_input": _lexeme(l6) if trajectory == "fixed_lambda6" else "",
            "X_lambda6_tanbeta_input": _lexeme(X) if trajectory == "fixed_X" else "",
        }
        closure = _potential_closure(inputs, values, precision)
        lexical.update({
            "tadpole_max_abs_normalized": closure["tadpole_max_abs_normalized"],
            "even_hessian_max_abs_normalized": closure["even_hessian_max_abs_normalized"],
            "CP_odd_hessian_max_abs_normalized": closure["CP_odd_hessian_max_abs_normalized"],
            "charged_hessian_max_abs_normalized": closure["charged_hessian_max_abs_normalized"],
        })
        closure_status = "PASS" if closure["pass"] else "INCONCLUSIVE"
        statuses = {
            name: {"status": "NOT_TESTED", "method": "not implemented in T2 lightweight constructor",
                  "observable": name, "tolerance": None, "precision_digits": precision, "evidence": None}
            for name in ("full_tree_unitarity", "BFB_general", "global_EW_minimum", "widths",
                         "branching_fractions", "ctau_mm", "displaced_photon_selection")
        }
        statuses.update({
            "spectrum_reconstruction": {"status": closure_status, "method": closure["method"],
                "observable": "CP-even, CP-odd, and charged Hessian blocks", "tolerance": closure["tolerance"],
                "precision_digits": precision, "evidence": closure},
            "alignment": {"status": "PASS" if closure["alignment_pass"] else "INCONCLUSIVE", "method": "independent Higgs-basis quartic transform and rotated Hessian",
                "observable": "Z6 and CP-even Hessian off-diagonal", "tolerance": closure["tolerance"], "precision_digits": precision,
                "evidence": {"Z6_independently_reconstructed": closure["Z6_independently_reconstructed"],
                             "Z6_abs_normalized": closure["Z6_abs_normalized"],
                             "even_hessian_offdiagonal_abs_GeV2": closure["even_hessian_offdiagonal_abs_GeV2"]}},
            "stationarity": {"status": closure_status, "method": closure["method"],
                "observable": "all six-field potential gradients at neutral vacuum", "tolerance": closure["tolerance"],
                "precision_digits": precision, "evidence": closure},
            "local_minimum_tree": {"status": ("PASS" if closure["pass"] and min(Decimal(x) for x in closure["physical_hessian_eigenvalues_GeV2"]) > 2*max(v*v,mh*mh,mp*mp,ma*ma,mhp*mhp)*Decimal(closure["normalized_hessian_residual"])
                                                    else "INCONCLUSIVE"),
                "method": "independently differentiated physical Hessian excluding Goldstones",
                "observable": "physical CP-even, CP-odd, charged eigenvalues", "tolerance": "> 0 GeV; normalized Hessian residual <= declared tolerance",
                "precision_digits": precision, "evidence": closure},
            "quartic_magnitude_screen": {"status": "NOT_TESTED", "method": "not requested by constructor",
                "observable": "all lambda_j and Z_j", "tolerance": "configurable; default 4*pi",
                "precision_digits": precision, "evidence": None},
        })
        return FamilyPoint(inputs, trajectory, lexical, statuses, {
            "formula_source": "CODEX_CONTRACT_T1.md", "precision_digits": precision,
            "arithmetic": "decimal.Decimal; reciprocal-based trigonometry",
            "independent_potential_closure": closure,
            "decimal_inputs_preserved": True,
            "trajectory_note": "fixed_lambda6 and fixed_X are separate paths; Xstar must be explicitly evaluated and passed as X",
            "coupling_conventions": {
                "L_int": "-g_hphiphi_GeV*h*phi**2/2; vertex -i*g",
                "V_int": "CphiHpHm_potential_GeV*phi*Hplus*Hminus; vertex -i*C",
            },
        })


def xstar_lambda6_tanbeta(*, tan_beta: Any, lambda1_ref: Any, mh_GeV: Any,
                          mphi_GeV: Any, v_GeV: Any, precision: int = 80) -> str:
    """Return Xstar solving Z7=0; it is an evaluated-X path, not fixed lambda6."""
    with localcontext() as ctx:
        ctx.prec = precision
        t = _decimal("tan_beta", tan_beta, positive=True)
        L = _decimal("lambda1_ref", lambda1_ref)
        mh = _decimal("mh_GeV", mh_GeV, positive=True)
        mp = _decimal("mphi_GeV", mphi_GeV, positive=True)
        v = _decimal("v_GeV", v_GeV, positive=True)
        r = Decimal(1) / t
        a = mh * mh / (v * v)
        k = L - a
        b = Decimal(2) * mp * mp / (v * v)
        xstar = (b - L - k * r * r) / (Decimal(1) + r * r)
        return str(+xstar)
