"""Conservative, independently checkable tree-level tests for aligned points.

The BFB/global-minimum PASS is deliberately a sufficient certificate for the
aligned Higgs-basis family (Z6=0), not a general numerical optimizer claim.
All polynomial inputs are interpreted as exact decimal rationals.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any

import numpy as np
import mpmath as mp


def _q(value: Any) -> Fraction:
    return Fraction(str(value))


def scattering_blocks(lambdas: dict[str, Any]) -> dict[str, np.ndarray]:
    """Return GI Eq. (7) eigenvalue matrices (the matrices are 8*pi*S)."""
    l = [_q(lambdas[f"lambda{i}"]) for i in range(1, 8)]
    a, b, c, d, e, f, g = map(float, l)
    root2 = np.sqrt(2.0)
    return {
        "Y2_sigma1": np.array([[a, e, root2*f], [e, b, root2*g],
                                [root2*f, root2*g, c+d]], dtype=float),
        "Y2_sigma0": np.array([[c-d]], dtype=float),
        "Y0_sigma1": np.array([[a, d, f, f], [d, b, g, g],
                               [f, g, c, e], [f, g, e, c]], dtype=float),
        "Y0_sigma0": np.array([[3*a, 2*c+d, 3*f, 3*f],
                               [2*c+d, 3*b, 3*g, 3*g],
                               [3*f, 3*g, c+2*d, 3*e],
                               [3*f, 3*g, 3*e, c+2*d]], dtype=float),
    }


def unitarity(lambdas: dict[str, Any], *, strict_bound: Any = "8*pi", precision: int = 80) -> dict[str, Any]:
    """Calculate all 12 independent eigenvalues and apply |Lambda| < 8*pi."""
    if precision < 30:
        raise ValueError("unitarity precision must be at least 30 decimal digits")
    def spectrum_at(dps):
        with mp.workdps(dps):
            a,b,c,d,e,f,g = [mp.mpf(str(lambdas[f"lambda{i}"])) for i in range(1,8)]
            root2=mp.sqrt(2)
            mats={"Y2_sigma1":[[a,e,root2*f],[e,b,root2*g],[root2*f,root2*g,c+d]],
                "Y2_sigma0":[[c-d]],
                "Y0_sigma1":[[a,d,f,f],[d,b,g,g],[f,g,c,e],[f,g,e,c]],
                "Y0_sigma0":[[3*a,2*c+d,3*f,3*f],[2*c+d,3*b,3*g,3*g],
                              [3*f,3*g,c+2*d,3*e],[3*f,3*g,3*e,c+2*d]]}
            out={}; max_res=mp.mpf(0); max_orth=mp.mpf(0); max_norm=mp.mpf(0)
            for name, rows in mats.items():
                mat=mp.matrix(rows); evals, vecs=mp.eigsy(mat)
                vals=[evals[i] for i in range(evals.rows)]
                out[name]=[mp.nstr(x,dps) for x in vals]
                n=mat.rows
                max_norm=max(max_norm,max(sum(abs(mat[i,j]) for j in range(n)) for i in range(n)))
                for i in range(n):
                    for j in range(n):
                        residual=abs(sum(mat[i,k]*vecs[k,j] for k in range(n))-vals[j]*vecs[i,j])
                        orth=abs(sum(vecs[k,i]*vecs[k,j] for k in range(n))-(1 if i==j else 0))
                        max_res=max(max_res,residual);max_orth=max(max_orth,orth)
            return out,max_res,max_orth,max_norm
    spectra,resid,orth,norm=spectrum_at(precision)
    high_spectra,_,_,_=spectrum_at(precision+20)
    with mp.workdps(precision+20):
        conv=max(abs(mp.mpf(high_spectra[n][i])-mp.mpf(spectra[n][i]))
                 for n in spectra for i in range(len(spectra[n])))
        measured=max(resid,orth*norm,conv)
        bound=8*mp.pi if strict_bound=="8*pi" else mp.mpf(str(strict_bound))
        flat=[mp.mpf(x) for spec in spectra.values() for x in spec]
        max_value=max(abs(x) for x in flat); max_residual=abs(bound-max_value)
        status=("INCONCLUSIVE" if max_residual<=measured else "PASS" if max_value<bound else "FAIL")
        bound_s,max_s=mp.nstr(bound,precision),mp.nstr(max_value,precision)
        a0_s,residual_s,err_s=(mp.nstr(max_value/(16*mp.pi),precision),
            mp.nstr(max_residual,precision),mp.nstr(measured,precision))
    return {
        "status": status,
        "method": "independent real-symmetric diagonalization of GI hep-ph/0508020 Eq. (7a-d)",
        "observable": "all eigenvalues Lambda of matrices 8*pi*S",
        "tolerance": {"criterion": "strict |Re(a0)| < 1/2", "a0": "Lambda/(16*pi)",
                      "equivalent": "strict |Lambda| < 8*pi", "bound": bound_s},
        "precision": f"mpmath symmetric eigensolver at {precision} and {precision+20} decimal digits; exact decimal inputs",
        "eigenvalues": spectra,
        "eigenvalue_count": len(flat),
        "maximum_absolute_eigenvalue": max_s,
        "maximum_absolute_a0": a0_s,
        "distance_from_strict_bound": residual_s,
        "boundary_error_budget": err_s,
        "measured_eigensolver_error_components": {"max_eigenpair_residual": mp.nstr(resid,precision),
            "max_eigenvector_orthogonality_error": mp.nstr(orth,precision),
            "matrix_infinity_norm": mp.nstr(norm,precision),
            "max_spectrum_change_at_plus_20_digits": mp.nstr(conv,precision)},
        "evidence": "complete spectrum retained for all four GI blocks",
    }


def _quartic(u: Fraction, b: Fraction, x: Fraction, y: Fraction,
             z: list[Fraction]) -> Fraction:
    z1,z2,z3,z4,z5,z6,z7 = z
    return (z1*u*u/2 + z2*b*b/2 + z3*u*b + (z4+z5)*x*x
            + (z4-z5)*y*y + 2*z6*u*x + 2*z7*b*x)


def _analytic_family_Z(inputs: dict[str, Any]) -> list[Fraction] | None:
    """Evaluate the stable Higgs-basis coefficients from exact decimal inputs."""
    required=("tan_beta","lambda1_ref","mh_GeV","mphi_GeV","mA_GeV","mHp_GeV","v_GeV")
    if not all(key in inputs for key in required) or not (inputs.get("lambda6_input") or inputs.get("X_lambda6_tanbeta_input")):
        return None
    t,L=_q(inputs["tan_beta"]),_q(inputs["lambda1_ref"])
    mh,mp,ma,mhp,v=(_q(inputs[key]) for key in ("mh_GeV","mphi_GeV","mA_GeV","mHp_GeV","v_GeV"))
    l6=_q(inputs["lambda6_input"]) if inputs.get("lambda6_input") else _q(inputs["X_lambda6_tanbeta_input"])/t
    r=1/t; a=mh*mh/(v*v); k=L-a; b=2*mp*mp/(v*v)
    return [a,L-4*l6*r+2*(b-L)*r**2-4*l6*r**3-3*k*r**4,
            2*mhp*mhp/(v*v),(mp*mp+ma*ma-2*mhp*mhp)/(v*v),
            (mp*mp-ma*ma)/(v*v),Fraction(0),
            -l6*(1+r*r)+(b-L)*r-k*r**3]


def _witness(z: list[Fraction]) -> dict[str, str] | None:
    # Rational points include interior, real/complex edges, and charged points.
    vals = [Fraction(0), Fraction(-2), Fraction(-1), Fraction(-1, 2),
            Fraction(1, 2), Fraction(1), Fraction(2)]
    for u in (Fraction(0), Fraction(1, 4), Fraction(1), Fraction(4)):
        for b in (Fraction(0), Fraction(1, 4), Fraction(1), Fraction(4)):
            if u == 0 or b == 0:
                candidates = [(Fraction(0), Fraction(0))]
            else:
                candidates = [(a * min(u,b), c * min(u,b)) for a in vals for c in vals
                              if a*a + c*c <= 1]
            for x,y in candidates:
                val = _quartic(u,b,x,y,z)
                if val < 0:
                    return {"u": str(u), "b": str(b), "x": str(x), "y": str(y),
                            "quartic_value": str(val), "domain_check": "u,b>=0 and x^2+y^2<=u*b exactly"}
    return None


def bfb_and_global_minimum(values: dict[str, Any], inputs: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Exact rational negative witness or sufficient aligned-family certificate."""
    family_z=_analytic_family_Z(inputs)
    z = family_z if family_z is not None else [_q(values[f"Z{i}"]) for i in range(1,8)]
    witness = _witness(z)
    if witness is not None:
        bfb = {"status": "FAIL", "method": "exact rational bilinear-domain witness",
               "observable": "quartic potential on u,b>=0, x^2+y^2<=u*b",
               "tolerance": "strict negative rational value", "precision": "exact Fraction arithmetic",
               "evidence": witness}
        global_min = {"status": "FAIL", "method": "exact quartic instability witness implies no global minimum",
                      "observable": "full tree potential", "tolerance": "certified potential ordering",
                      "precision": "exact rational witness does not establish the vacuum ordering",
                      "evidence": "potential is unbounded below along the reported bilinear ray; no point can be a global minimum"}
        return bfb, global_min

    z1,z2,z3,z4,z5,z6,z7=z
    c = z3 + min(Fraction(0), z4+z5, z4-z5)
    # V - V_EW = .5 Z1 (u-v^2/2)^2 + .5 Z2 b^2 + c*u*b + 2 Z7*b*x,
    # and |x|<=sqrt(u*b). Completing the square in sqrt(u/b) proves the
    # following exact sufficient condition, including interior and boundary.
    margin = z2/Fraction(2) - (z7*z7/c if c > 0 else Fraction(0))
    cert_ok = z1 > 0 and c > 0 and z6 == 0 and margin > 0
    cert = {"certificate": "Z1>0, c>0, Z2/2-Z7^2/c>0, Z6=0",
            "Z1": str(z1), "c": str(c), "quadratic_margin": str(margin),
            "Z6": str(z6), "all_comparisons": "exact decimal rational arithmetic",
            "domain": "full gauge-invariant bilinear domain, including interior and boundary",
            "scope": ("exact symbolic aligned-family formulas evaluated from input rationals; no backend-double transfer claim"
                      if family_z is not None else "exact supplied polynomial; constructor provenance not established"),
            "serialized_Z_minus_exact_family_Z": ([str(_q(values[f"Z{i}"])-z[i-1]) for i in range(1,8)]
                                                   if family_z is not None else None)}
    if cert_ok:
        bfb_status = "PASS"
        bfb_method = "exact sufficient global quartic bound on bilinear domain"
        closure = (inputs.get("_constructor_closure") or {})
        if all(f"Y{i}_GeV2" in values for i in (1,2,3)):
            qY1, qY2, qY3 = (_q(values[f"Y{i}_GeV2"]) for i in (1,2,3))
            expected_y1=(-_q(inputs["mh_GeV"])**2/2) if family_z is not None else None
            stationarity_residual=(qY1+z1*_q(inputs["v_GeV"])**2/2) if family_z is not None else None
            stationarity_scale=max(Fraction(1),_q(inputs["v_GeV"])**2) if family_z is not None else Fraction(1)
            closure_tol=Fraction(str(closure.get("tolerance","0"))) if closure.get("tolerance") else Fraction(0)
            stationary_input_match=(family_z is not None and qY1==expected_y1 and qY2==0 and qY3==0
                                    and abs(stationarity_residual)/stationarity_scale<=closure_tol)
            global_gate = (family_z is not None and closure.get("pass") is True and stationary_input_match
                           and z6==0 and qY1<0 and qY2==0 and qY3==0)
        else:
            qY1, qY2, qY3 = Fraction(0), Fraction(1), Fraction(1)
            global_gate = False
        global_status = "PASS" if global_gate else "INCONCLUSIVE"
        global_method = "exact aligned-potential decomposition; constructor closure checked independently"
        uew = -qY1/z1 if z1 else Fraction(0)
        v_requested = _q(inputs["v_GeV"]) if inputs.get("v_GeV") else None
        vev_shift = None
        if v_requested is not None and uew > 0:
            with mp.workdps(50):
                v_actual = mp.sqrt(2*mp.mpf(uew.numerator)/uew.denominator)
                vev_shift = {"requested_v_GeV": str(v_requested), "stationary_v_GeV": mp.nstr(v_actual,50),
                             "absolute_shift_GeV": mp.nstr(abs(v_actual-mp.mpf(str(v_requested))),50),
                             "reason": "exact family formulas give uEW=v_requested^2/2"}
        global_evidence = ({"proof": "V-VEW >= .5*Z1*(u-uEW)^2 + .5*Z2*b^2 + c*u*b + 2*Z7*b*x, with c*ub bounding all real, complex and charged bilinears",
                           "uEW": str(uew), "Y1": str(qY1), "Y2": str(qY2), "Y3": str(qY3),
                           "requested_vs_stationary_vev": vev_shift,
                           "constructor_closure": "PASS", "stationarity_residual_exact_GeV2": str(stationarity_residual) if family_z is not None else None,
                           "gauge_orbit": "all bilinears with b=0 and u=uEW are degenerate EW representatives"}
                          if global_gate else "quartic is certified BFB, but full quadratic terms or independent constructor closure do not establish a stationary EW global minimum")
    else:
        bfb_status = global_status = "INCONCLUSIVE"
        bfb_method = "rational witness search and exact sufficient certificate"
        global_method = "no complete stationary-point enumeration/certified ordering"
        global_evidence = "no negative rational witness found; positive sufficient certificate unavailable"
    return (
        {"status": bfb_status, "method": bfb_method,
         "observable": "quartic potential on u,b>=0, x^2+y^2<=u*b",
         "tolerance": "exact positivity; no floating tolerance", "precision": "exact Fraction arithmetic",
         "evidence": cert if cert_ok else {"certificate": cert, "negative_witness": None}},
        {"status": global_status, "method": global_method, "observable": "full tree potential",
         "tolerance": "certified V >= V_EW", "precision": "exact Fraction arithmetic",
         "evidence": global_evidence if cert_ok else global_evidence},
    )


def assess_family_point(point: dict[str, Any], *, quartic_limit: Any = "4*pi") -> dict[str, Any]:
    values, inputs = point["values"], dict(point["inputs"])
    inputs["_constructor_closure"] = point.get("provenance", {}).get("independent_potential_closure", {})
    lambdas = {f"lambda{i}": values[f"lambda{i}"] for i in range(1,8)}
    zvals = {f"Z{i}": values[f"Z{i}"] for i in range(1,8)}
    bfb, global_min = bfb_and_global_minimum(values, inputs)
    unit = unitarity(lambdas)
    limit = float(4*np.pi) if quartic_limit == "4*pi" else float(quartic_limit)
    biggest = max(abs(float(_q(v))) for v in list(lambdas.values())+list(zvals.values()))
    qscreen = {"status": "PASS" if biggest < limit else "FAIL",
               "method": "configurable absolute quartic magnitude screen",
               "observable": "all lambda1..7 and Z1..7", "tolerance": f"strict max(abs(q)) < {limit}",
               "precision": "exact decimal rational parse; displayed comparison in binary64",
               "maximum_absolute_quartic": biggest,
               "evidence": {"all_lambdas": lambdas, "all_Z": zvals}}
    statuses = dict(point.get("constraints", {}))
    for key in ("spectrum_reconstruction", "stationarity", "alignment", "local_minimum_tree"):
        status = statuses.get(key, {}).get("status", "NOT_TESTED")
        statuses[key] = {"status": status, "method": "independent T1 potential closure in aligned_lambda7.py",
                         "observable": key, "tolerance": statuses.get(key, {}).get("tolerance", "see T1 closure"),
                         "precision": "Decimal arbitrary precision; six-field potential differentiated independently",
                         "evidence": statuses.get(key, {})}
    statuses.update({"quartic_screen_4pi": qscreen, "BFB_general": bfb,
                     "global_EW_minimum": global_min, "full_tree_unitarity": unit})
    for key in ("widths", "branching_fractions", "ctau_mm", "displaced_photon_selection"):
        statuses[key] = {"status": "NOT_TESTED", "method": "not evaluated in this independent constraint module",
                         "observable": key, "tolerance": "not applicable", "precision": "not applicable",
                         "evidence": "requires T4 calculator integration and applicable documented selection"}
    return {"status": "SCREENED_ONLY", "constraints": statuses}
