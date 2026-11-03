"""Exact-rational certificates for the general real 2HDM bilinear potential.

The routines in this module deliberately return INCONCLUSIVE when exact
feasibility or completeness is not established.  Floating point optimizers
are never used as proof evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Any, Mapping

import mpmath as mp
import sympy as sp


def _q(value: Any) -> sp.Rational:
    """Interpret decimal strings and integers as exact rationals."""
    if isinstance(value, sp.Expr) and not value.has(sp.Float):
        return value  # algebraic stationary-point coordinates stay exact
    if isinstance(value, sp.Rational):
        return value
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    return sp.Rational(str(value))


def _s(value: Any) -> str:
    return str(value)


@dataclass(frozen=True)
class RealPotential:
    """Coefficients of V(u,b,x,y) in the Higgs basis; all values are exact."""

    Y: tuple[sp.Rational, sp.Rational, sp.Rational]
    Z: tuple[sp.Rational, sp.Rational, sp.Rational, sp.Rational,
             sp.Rational, sp.Rational, sp.Rational]
    v2: sp.Rational | None = None

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any]) -> "RealPotential":
        return cls(tuple(_q(values[f"Y{i}_GeV2"]) for i in range(1, 4)),
                   tuple(_q(values[f"Z{i}"]) for i in range(1, 8)),
                   _q(values["v2_GeV2"]) if "v2_GeV2" in values else None)

    def as_dict(self) -> dict[str, str]:
        out = {f"Y{i}_GeV2": _s(self.Y[i - 1]) for i in range(1, 4)}
        out.update({f"Z{i}": _s(self.Z[i - 1]) for i in range(1, 8)})
        if self.v2 is not None:
            out["v2_GeV2"] = _s(self.v2)
        return out


def potential_value(p: RealPotential, u: Any, b: Any, x: Any, y: Any) -> sp.Expr:
    u, b, x, y = map(_q, (u, b, x, y))
    y1, y2, y3 = p.Y
    z1, z2, z3, z4, z5, z6, z7 = p.Z
    return sp.expand(y1*u + y2*b + 2*y3*x + z1*u**2/2 + z2*b**2/2
                     + z3*u*b + (z4+z5)*x**2 + (z4-z5)*y**2
                     + 2*z6*u*x + 2*z7*b*x)


def _quartic(p: RealPotential, u: Any, b: Any, x: Any, y: Any) -> sp.Expr:
    return sp.expand(potential_value(
        RealPotential((sp.S.Zero,)*3, p.Z), u, b, x, y))


def _bilinear_domain(u: sp.Expr, b: sp.Expr, x: sp.Expr, y: sp.Expr) -> bool:
    signs = (sp.ask(sp.Q.nonnegative(u)), sp.ask(sp.Q.nonnegative(b)),
             sp.ask(sp.Q.nonnegative(u*b-x*x-y*y)))
    return all(s is True for s in signs)


def _slemma_matrix(p: RealPotential, tau: sp.Rational) -> sp.Matrix:
    z1,z2,z3,z4,z5,z6,z7 = p.Z
    c = (z1+z2+2*z3)/8
    eta = sp.Matrix([(z6+z7)/4, 0, (z1-z2)/8])
    e = sp.Matrix([[ (z4+z5)/4, 0, (z6-z7)/4],
                   [0, (z4-z5)/4, 0],
                   [(z6-z7)/4, 0, (z1+z2-2*z3)/8]])
    return sp.BlockMatrix([[sp.Matrix([[c-tau]]), eta.T],
                           [eta, e+tau*sp.eye(3)]]).as_explicit()


def _psd_exact(matrix: sp.Matrix) -> bool:
    """Exact PSD test via all principal minors (valid for symmetric matrices)."""
    n = matrix.rows
    if matrix.cols != n or matrix != matrix.T:
        return False
    for size in range(1, n+1):
        for indices in __import__("itertools").combinations(range(n), size):
            minor = sp.factor(matrix.extract(indices, indices).det())
            if sp.ask(sp.Q.nonnegative(minor)) is not True:
                return False
    return True


def _negative_witness(p: RealPotential) -> dict[str, str] | None:
    # Exact rational samples cover neutral-real, neutral-complex, and charged
    # directions.  Absence of a witness is never treated as a proof.
    samples = (sp.Rational(0), sp.Rational(-2), sp.Rational(-1),
               sp.Rational(-1, 2), sp.Rational(1, 2), sp.Rational(1), sp.Rational(2))
    scales = (sp.Rational(1, 4), sp.Rational(1), sp.Rational(4))
    for u in scales:
        value=_quartic(p,u,0,0,0)
        if value<0:
            return {"u":_s(u),"b":"0","x":"0","y":"0","quartic_value":_s(value),
                    "domain_check":"exact u>=0, b=0, x=y=0"}
    for b in scales:
        value=_quartic(p,0,b,0,0)
        if value<0:
            return {"u":"0","b":_s(b),"x":"0","y":"0","quartic_value":_s(value),
                    "domain_check":"exact u=0, b>=0, x=y=0"}
    for u in scales:
        for b in scales:
            for ix in samples:
                for iy in samples:
                    x, y = ix*min(u,b), iy*min(u,b)
                    if x*x+y*y > u*b:
                        continue
                    value = _quartic(p, u, b, x, y)
                    if value < 0:
                        return {"u": _s(u), "b": _s(b), "x": _s(x), "y": _s(y),
                                "quartic_value": _s(value),
                                "domain_check": "exact u>=0, b>=0, u*b-x^2-y^2>=0"}
    return None


def certify_bfb(p: RealPotential) -> dict[str, Any]:
    """Return an exact negative witness or an exact S-lemma certificate.

    A PSD certificate proves nonnegative quartic stability.  PD proves strict
    positivity on the normalized ball.  Failure to find either is INCONCLUSIVE.
    """
    witness = _negative_witness(p)
    if witness:
        return {"status": "FAIL", "method": "exact rational bilinear-domain witness",
                "observable": "quartic potential", "tolerance": "strictly negative exact value",
                "precision": "exact rational arithmetic", "evidence": witness}
    # High-precision dual optimization proposes tau; exact rational principal
    # minors remain the only acceptance test. A scale-aware rational grid is
    # a bounded fallback for difficult singular/marginal cases.
    zscale=max((abs(v) for v in p.Z),default=sp.S.Zero)
    taus=[]
    if zscale:
        z1,z2,z3,z4,z5,z6,z7=p.Z
        c=(z1+z2+2*z3)/8
        eta=((z6+z7)/4,sp.S.Zero,(z1-z2)/8)
        e=sp.Matrix([[(z4+z5)/4,0,(z6-z7)/4],
                     [0,(z4-z5)/4,0],
                     [(z6-z7)/4,0,(z1+z2-2*z3)/8]])
        with mp.workdps(100):
            em=mp.matrix([[mp.mpf(str(e[i,j].p))/mp.mpf(str(e[i,j].q)) for j in range(3)] for i in range(3)])
            ev,_=mp.eigsy(em)
            low=max(mp.mpf(0),-min(ev[i] for i in range(3)))
            scale=mp.mpf(str(zscale.p))/mp.mpf(str(zscale.q))
            eps=max(scale*mp.mpf("1e-90"),mp.mpf("1e-190"))
            vv=mp.matrix([mp.mpf(str(v.p))/mp.mpf(str(v.q)) for v in eta])
            def derivative(tau):
                bm=em+tau*mp.eye(3)
                w=mp.lu_solve(bm,vv)
                return -1+sum(value*value for value in w)
            left=low+eps
            try:
                if derivative(left)<=0:
                    proposal=left
                else:
                    right=left+scale
                    while derivative(right)>0 and right<low+scale*mp.mpf("1e80"):
                        right=low+2*(right-low)
                    if derivative(right)>0:
                        right=low+scale*mp.mpf("1e80")
                    for _ in range(400):
                        mid=(left+right)/2
                        if derivative(mid)>0:
                            left=mid
                        else:
                            right=mid
                    proposal=(left+right)/2
                taus.append(sp.Rational(mp.nstr(proposal,95)))
            except (ZeroDivisionError, ValueError):
                pass
        taus.extend(zscale*sp.Rational(i,20) for i in range(0,41))
        taus.extend(zscale*sp.Rational(i) for i in (4,8,16,32,64))
    else:
        taus=[sp.S.Zero]
    taus=list(dict.fromkeys(taus))
    semidefinite = None
    for tau in taus:
        matrix = _slemma_matrix(p, tau)
        if _psd_exact(matrix):
            strict = all(sp.ask(sp.Q.positive(sp.factor(matrix.extract(ix, ix).det()))) is True
                         for ix in ((0,), (1,), (2,), (3,), (0,1), (0,2), (0,3),
                                    (1,2), (1,3), (2,3), (0,1,2), (0,1,3),
                                    (0,2,3), (1,2,3), (0,1,2,3)))
            result = {"status": "PASS", "method": "exact rational S-lemma certificate",
                    "observable": "quartic potential on full normalized bilinear ball",
                    "tolerance": "nonnegative" if not strict else "strictly positive",
                    "precision": "exact rational principal-minor comparisons",
                    "evidence": {"tau": _s(tau), "matrix": [[_s(v) for v in row] for row in matrix.tolist()],
                                 "certificate_scope": "full gauge-invariant domain, interior and boundary",
                                 "strict_certificate": strict}}
            if strict:
                return result
            if semidefinite is None:
                semidefinite = result
    if semidefinite is not None:
        return semidefinite
    return {"status": "INCONCLUSIVE", "method": "exact witness search and rational S-lemma certificate",
            "observable": "quartic potential on full bilinear domain",
            "tolerance": "nonnegative", "precision": "exact rational arithmetic",
            "evidence": "No negative rational witness or S-lemma certificate found; no optimizer PASS is used."}


def _atan_interval(q: sp.Rational, terms: int) -> tuple[sp.Rational, sp.Rational]:
    total = sum(((-1)**k*q**(2*k+1)/ (2*k+1) for k in range(terms)), sp.S.Zero)
    next_term = q**(2*terms+1)/(2*terms+1)
    return (total, total+next_term) if terms % 2 == 0 else (total-next_term, total)


def _pi_interval(terms: int) -> tuple[sp.Rational, sp.Rational]:
    # Machin's identity with rigorous alternating-series bounds.
    lo5,hi5=_atan_interval(sp.Rational(1,5),terms)
    lo239,hi239=_atan_interval(sp.Rational(1,239),terms)
    return 16*lo5-4*hi239,16*hi5-4*lo239


def certify_tree_unitarity(lambdas: Mapping[str, Any], *, max_terms: int = 120) -> dict[str, Any]:
    """Certify |Lambda|<8*pi from exact GI blocks and rigorous rational pi bounds.

    The sqrt(2) off-diagonal factors cancel from each characteristic
    determinant, so principal-minor signs are exact rationals.
    """
    if max_terms < 4:
        raise ValueError("max_terms must be at least 4")
    l=[_q(lambdas[f"lambda{i}"]) for i in range(1,8)]
    a,b,c,d,e,f,g=l
    triples={
        "Y2_sigma1":(a,b,c+d,e,f,g),
        "Y0_sigma1":(a,b,c+e,d,f,g),
        "Y0_sigma0":(3*a,3*b,c+2*d+3*e,2*c+d,3*f,3*g),
    }
    scalars={"Y2_sigma0":c-d,"Y0_sigma1_antisymmetric":c-e,
             "Y0_sigma0_antisymmetric":c+2*d-3*e}

    def minors(sign: int, diagonal: tuple[sp.Rational,sp.Rational,sp.Rational],
               off: tuple[sp.Rational,sp.Rational,sp.Rational], threshold: sp.Rational):
        aa,cc,dd=diagonal
        bb,ff,gg=off
        d1,d2,d3=threshold+sign*aa,threshold+sign*cc,threshold+sign*dd
        all_minors=(d1,d2,d3,
                    sp.expand(d1*d2-bb*bb),
                    sp.expand(d1*d3-2*ff*ff),
                    sp.expand(d2*d3-2*gg*gg))
        leading=(d1,sp.expand(d1*d2-bb*bb))
        det=sp.expand(d1*d2*d3+sign*4*bb*ff*gg
                      -2*d1*gg*gg-2*d2*ff*ff-d3*bb*bb)
        return leading+(det,),all_minors+(det,)

    terms=4
    while terms<=max_terms:
        pi_lo,pi_hi=_pi_interval(terms)
        tlo,thi=8*pi_lo,8*pi_hi
        pass_ok=True
        fail_evidence=[]
        for name,(aa,cc,dd,bb,ff,gg) in triples.items():
            plus_lo=minors(1,(aa,cc,dd),(bb,ff,gg),tlo)
            minus_lo=minors(-1,(aa,cc,dd),(bb,ff,gg),tlo)
            pass_ok &= all(v>0 for v in (*plus_lo[0],*minus_lo[0]))
            for sign,label in ((1,"+"),(-1,"-")):
                bad=minors(sign,(aa,cc,dd),(bb,ff,gg),thi)[1]
                if any(v<0 for v in bad):
                    fail_evidence.append({"block":name,"matrix":"T I "+label+" M",
                                         "negative_principal_minor":_s(next(v for v in bad if v<0))})
        for name,value in scalars.items():
            pass_ok &= abs(value)<tlo
            if abs(value)>thi:
                fail_evidence.append({"block":name,"scalar_eigenvalue":_s(value),"strict_bound_upper":_s(thi)})
        if pass_ok:
            return {"status":"PASS","method":"exact GI 3x3 principal minors and Machin pi interval",
                    "observable":"complete tree-level scattering spectrum",
                    "tolerance":"strict |Lambda| < 8*pi (a0=Lambda/(16*pi))",
                    "precision":f"exact rational coefficients; {terms} alternating terms per arctangent",
                    "evidence":{"pi_interval":[_s(pi_lo),_s(pi_hi)],"8pi_lower_bound":_s(tlo),
                               "triples":list(triples),"scalars":{k:_s(v) for k,v in scalars.items()},
                               "all_six_shifted_blocks_positive_definite":True}}
        if fail_evidence:
            return {"status":"FAIL","method":"exact GI principal-minor instability against rigorous pi upper bound",
                    "observable":"complete tree-level scattering spectrum",
                    "tolerance":"strict |Lambda| < 8*pi (a0=Lambda/(16*pi))",
                    "precision":f"exact rational coefficients; {terms} alternating terms per arctangent",
                    "evidence":{"pi_interval":[_s(pi_lo),_s(pi_hi)],"8pi_upper_bound":_s(thi),
                               "violations":fail_evidence}}
        terms*=2
    return {"status":"INCONCLUSIVE","method":"exact GI principal-minor bounds with finite Machin refinement",
            "observable":"complete tree-level scattering spectrum",
            "tolerance":"strict |Lambda| < 8*pi (a0=Lambda/(16*pi))",
            "precision":f"exact rational coefficients; refinement capped at {max_terms} terms",
            "evidence":"Principal-minor signs do not separate all eigenvalues from the strict bound within the configured rational pi interval."}


def _exact_sign(expr: sp.Expr) -> int | None:
    expr = sp.simplify(expr)
    if expr == 0 or expr.is_zero is True:
        return 0
    sign = sp.sign(expr)
    if sign == 1:
        return 1
    if sign == -1:
        return -1
    return None


def _candidate(u: sp.Expr, b: sp.Expr, x: sp.Expr, y: sp.Expr,
               energy: sp.Expr, stratum: str) -> dict[str, Any]:
    # Keep the already exact expressions intact. Simplifying large CRootOf
    # expressions here can dominate certification time without adding proof.
    return {"u": _s(u), "b": _s(b), "x": _s(x), "y": _s(y),
            "V_GeV4": _s(energy), "stratum": stratum}


def enumerate_stationary_candidates(p: RealPotential) -> dict[str, Any]:
    """Enumerate exactly solvable strata and report whether coverage is complete.

    The aligned real-neutral branch is solved by exact real-root isolation.
    A nonsingular interior affine solution is also checked exactly. General
    nonaligned neutral KKT strata and rank-deficient affine cones remain
    explicitly INCONCLUSIVE.
    """
    y1,y2,y3 = p.Y
    z1,z2,z3,z4,z5,z6,z7 = p.Z
    u,b,x,y,r = sp.symbols("u b x y r", real=True)
    candidates = [_candidate(0,0,0,0,0,"origin")]
    unresolved = []
    aligned = (y2 == 0 and y3 == 0 and z6 == 0)
    if aligned:
        c = z3+z4+z5
        P = y1
        Q = z1/2 + 2*z6*r + c*r**2 + 2*z7*r**3 + z2*r**4/2
        # The r=0 stationary point is the u-axis branch when Q(0)>0.
        roots_poly = sp.Poly(sp.expand(2*Q*sp.diff(P,r)-P*sp.diff(Q,r)), r, domain=sp.QQ)
        roots = []
        if not roots_poly.is_zero:
            try:
                roots = sp.real_roots(roots_poly.as_expr())
            except (NotImplementedError, sp.PolynomialError):
                unresolved.append("Exact real-root isolation failed for aligned neutral stationary polynomial.")
        for root in roots:
            qv = sp.simplify(Q.subs(r, root))
            pv = sp.simplify(P)
            qsign = _exact_sign(qv)
            usign = _exact_sign(-pv/(2*qv)) if qsign not in (0, None) else None
            if qsign == 0:
                if pv != 0:
                    continue  # dV/du=P+2*u*Q cannot vanish
                qprime=sp.simplify(sp.diff(Q,r).subs(r,root))
                if _exact_sign(qprime)==0:
                    uv=sp.S.One;bv=sp.simplify(root**2);xv=root
                    candidates.append(_candidate(uv,bv,xv,0,sp.S.Zero,
                                                 "neutral-real; singular stationary continuum witness"))
                    unresolved.append(f"Aligned neutral root r={root} has a stationary continuum in u.")
                continue
            if usign == 1:
                uv = sp.simplify(-pv/(2*qv))
                bv = sp.simplify(root**2*uv)
                xv = sp.simplify(root*uv)
                energy = sp.simplify(-pv*pv/(4*qv))
                candidates.append(_candidate(uv,bv,xv,0,energy,"neutral-real"))
            elif usign is None:
                unresolved.append(f"Could not certify sign of u at aligned neutral root r={root}.")
        # If Q is identically zero, a continuum may be stationary or singular.
        if roots_poly.is_zero:
            unresolved.append("Aligned neutral stationary polynomial vanishes identically; stationary continuum is unresolved.")
    # General neutral-real stationary points, including nonaligned inputs.
    # b=r^2 u and x=r u reduce the constrained potential to u P(r)+u^2 Q(r).
    P = y1+2*y3*r+y2*r**2
    Q = z1/2+2*z6*r+(z3+z4+z5)*r**2+2*z7*r**3+z2*r**4/2
    critical = sp.Poly(sp.expand(2*Q*sp.diff(P,r)-P*sp.diff(Q,r)),r,domain=sp.QQ)
    try:
        neutral_roots = [] if critical.is_zero else sp.real_roots(critical.as_expr())
    except (NotImplementedError, sp.PolynomialError):
        neutral_roots = []
        unresolved.append("Exact real-root isolation failed for the general neutral stationary polynomial.")
    if critical.is_zero:
        unresolved.append("General neutral stationary polynomial vanishes identically; singular continuum needs exact feasibility analysis.")
    elif aligned:
        # Already enumerated by the aligned branch above.
        pass
    else:
        for root in neutral_roots:
            pv,qv = sp.simplify(P.subs(r,root)),sp.simplify(Q.subs(r,root))
            ps,qs = _exact_sign(pv),_exact_sign(qv)
            if ps == 0 and qs == 0:
                pp,qp = sp.simplify(sp.diff(P,r).subs(r,root)),sp.simplify(sp.diff(Q,r).subs(r,root))
                pps,qps = _exact_sign(pp),_exact_sign(qp)
                if qps not in (0,None):
                    uv=sp.simplify(-pp/qp)
                    if _exact_sign(uv)==1:
                        bv=sp.simplify(root**2*uv);xv=sp.simplify(root*uv)
                        candidates.append(_candidate(uv,bv,xv,0,sp.S.Zero,"neutral-real singular"))
                    else:
                        unresolved.append(f"Singular neutral root r={root} lacks certified positive u.")
                elif pps == 0 and qps == 0:
                    uv=sp.S.One;bv=sp.simplify(root**2);xv=root
                    candidates.append(_candidate(uv,bv,xv,0,potential_value(p,uv,bv,xv,0),
                                                 "neutral-real; singular stationary continuum witness"))
                    unresolved.append(f"Neutral root r={root} has a stationary continuum in u.")
                else:
                    # With P=Q=0 and Q'=0 but P' nonzero there is no positive
                    # u satisfying the radial derivative equation.
                    continue
            elif qs == 0:
                # Radial stationarity requires P+2uQ=0; with Q=0 and P!=0
                # this root is exactly excluded.
                pass
            else:
                uv=sp.simplify(-pv/(2*qv))
                us=_exact_sign(uv)
                if us == 1:
                    bv=sp.simplify(root**2*uv);xv=sp.simplify(root*uv)
                    candidates.append(_candidate(uv,bv,xv,0,sp.simplify(-pv*pv/(4*qv)),"neutral-real"))
                elif us is None:
                    unresolved.append(f"Could not certify sign of u at neutral root r={root}.")

    # Both one-doublet axes are checked separately. At an axis, derivatives
    # multiplying the empty doublet's norm are quadratic in field variations;
    # only the occupied-direction and cross-bilinear tadpoles must vanish.
    if z1 != 0:
        u_axis=sp.simplify(-y1/z1)
        if _exact_sign(u_axis)==1 and sp.simplify(y3+z6*u_axis)==0:
            candidates.append(_candidate(u_axis,0,0,0,potential_value(p,u_axis,0,0,0),"b=0 axis"))
    elif y1 == 0:
        if z6 != 0:
            u_axis=sp.simplify(-y3/z6)
            if _exact_sign(u_axis)==1:
                candidates.append(_candidate(u_axis,0,0,0,potential_value(p,u_axis,0,0,0),"b=0 axis"))
        elif y3 == 0:
            candidates.append(_candidate(1,0,0,0,potential_value(p,1,0,0,0),"b=0 stationary continuum witness"))
            unresolved.append("The b=0 axis contains a stationary continuum; connected-set ordering is unresolved.")
    if z2 != 0:
        b_axis=sp.simplify(-y2/z2)
        if _exact_sign(b_axis)==1 and sp.simplify(y3+z7*b_axis)==0:
            candidates.append(_candidate(0,b_axis,0,0,potential_value(p,0,b_axis,0,0),"u=0 axis"))
    elif y2 == 0:
        if z7 != 0:
            b_axis=sp.simplify(-y3/z7)
            if _exact_sign(b_axis)==1:
                candidates.append(_candidate(0,b_axis,0,0,potential_value(p,0,b_axis,0,0),"u=0 axis"))
        elif y3 == 0:
            candidates.append(_candidate(0,1,0,0,potential_value(p,0,1,0,0),"u=0 stationary continuum witness"))
            unresolved.append("The u=0 axis contains a stationary continuum; connected-set ordering is unresolved.")

    # Interior stationarity is an affine system in (u,b,x); y is free only
    # when Z4=Z5. A unique solution is retained if it lies strictly inside.
    grad = sp.Matrix([[z1,z3,2*z6], [z3,z2,2*z7],
                      [2*z6,2*z7,2*(z4+z5)]])
    rhs = sp.Matrix([-y1,-y2,-2*y3])
    determinant = sp.factor(grad.det())
    if determinant != 0:
        sol = grad.inv()*rhs
        uv,bv,xv = map(sp.simplify, sol)
        delta = sp.simplify(uv*bv-xv*xv)
        su,sb,sd = _exact_sign(uv),_exact_sign(bv),_exact_sign(delta)
        if su == 1 and sb == 1 and sd == 1:
            if z4 == z5:
                # Any y in the open interval is stationary; energy does not
                # depend on y. The point y=0 still proves interior feasibility.
                candidates.append(_candidate(uv,bv,xv,0,potential_value(p,uv,bv,xv,0),"charge/interior; y-flat"))
            else:
                candidates.append(_candidate(uv,bv,xv,0,potential_value(p,uv,bv,xv,0),"charge/interior"))
        elif None in (su,sb,sd):
            unresolved.append("Could not certify exact bilinear-domain signs for the interior stationary solution.")
    else:
        try:
            family=sp.linsolve((grad,rhs),(u,b,x))
            if family is sp.EmptySet or family == sp.EmptySet:
                pass  # inconsistent linear system: this stratum is excluded exactly
            else:
                unresolved.append("Interior stationarity has a singular affine family; feasibility/completeness not certified.")
                for row in family:
                    free=sorted(set().union(*(entry.free_symbols for entry in row)),key=str)
                    witness=[sp.simplify(entry.subs({symbol:0 for symbol in free})) for entry in row]
                    wu,wb,wx=witness
                    delta=sp.simplify(wu*wb-wx*wx)
                    if (_exact_sign(wu)==1 and _exact_sign(wb)==1 and _exact_sign(delta)==1):
                        candidates.append(_candidate(wu,wb,wx,0,potential_value(p,wu,wb,wx,0),
                                                     "interior; singular-family exact feasible witness"))
                        break
        except (ValueError, NotImplementedError):
            unresolved.append("Interior stationarity has a singular affine family; feasibility/completeness not certified.")
            unresolved.append("Could not extract an exact feasible witness from the singular interior affine family.")

    # Complex-neutral/charge-boundary branch y != 0 has fixed multiplier
    # mu=-(Z4-Z5). Solve the resulting exact affine stationarity equations.
    mu = -(z4-z5)
    cm = grad.copy()
    cm[0,1] -= mu
    cm[1,0] -= mu
    cm[2,2] += 2*mu
    crhs = sp.Matrix([-y1,-y2,-2*y3])
    if sp.factor(cm.det()) != 0:
        cu,cb,cx = map(sp.simplify, cm.inv()*crhs)
        ysq = sp.simplify(cu*cb-cx*cx)
        if (_exact_sign(cu)==1 and _exact_sign(cb)==1 and _exact_sign(ysq)==1):
            # The potential is independent of the sign of y; retain a real
            # square-root representative and its exact energy.
            yy=sp.sqrt(ysq)
            candidates.append(_candidate(cu,cb,cx,yy,potential_value(p,cu,cb,cx,yy),"neutral-complex; fixed multiplier"))
        elif None in (_exact_sign(cu),_exact_sign(cb),_exact_sign(ysq)):
            unresolved.append("Could not certify exact feasibility of the fixed-multiplier complex stationary solution.")
    else:
        unresolved.append("Fixed-multiplier complex stationary equations are singular; affine feasibility not certified.")

    # Deduplicate identical exact candidate records while keeping source order.
    unique = []
    seen = set()
    for item in candidates:
        key = tuple(item[k] for k in ("u","b","x","y","V_GeV4"))
        if key not in seen:
            seen.add(key)
            unique.append(item)
    return {"status": "PASS" if not unresolved else "INCONCLUSIVE",
            "method": "exact rational/algebraic stationary-stratum equations",
            "observable": "competing stationary candidates of the full real potential",
            "tolerance": "exact feasibility and potential ordering",
            "precision": "exact QQ coefficients and SymPy real algebraic roots",
            "complete": not unresolved, "candidates": unique,
            "neutral_real_critical_polynomial": _s(sp.factor(critical.as_expr())),
            "neutral_real_critical_degree": critical.degree() if not critical.is_zero else None,
            "unresolved_reasons": unresolved,
            "ordering": "Candidate energies are exact expressions; completeness is required before global PASS."}


def certify_runaway(p: RealPotential) -> dict[str, Any]:
    """Check exact polynomial rays, including the flat-quartic mixed-scale ray."""
    n=sp.symbols("n", positive=True)
    curves=((n**3,n,n**2,0,"u=n^3,b=n,x=n^2,y=0"),
            (n,n**3,n**2,0,"u=n,b=n^3,x=n^2,y=0"),
            (n,n,0,0,"u=n,b=n,x=y=0"),
            (n,0,0,0,"u=n,b=x=y=0"),
            (0,n,0,0,"u=x=y=0,b=n"))
    for u0,b0,x0,y0,label in curves:
        if not _bilinear_domain(u0,b0,x0,y0):
            continue
        poly=sp.Poly(sp.expand(potential_value(p,u0,b0,x0,y0)),n,domain=sp.QQ)
        if poly.is_zero:
            continue
        coefficient=poly.LC()
        if coefficient < 0:
            return {"status":"FAIL","method":"exact polynomial runaway ray",
                    "observable":"boundedness of the full tree potential",
                    "tolerance":"V(n) tends to minus infinity for n tending to plus infinity",
                    "precision":"exact rational coefficients",
                    "evidence":{"curve":label,"V_of_n":_s(poly.as_expr()),
                               "leading_power":poly.degree(),"leading_coefficient":_s(coefficient),
                               "bilinear_domain":"verified exactly for n>=0"}}
    return {"status":"INCONCLUSIVE","method":"finite exact runaway-ray witnesses",
            "observable":"boundedness of the full tree potential",
            "tolerance":"global lower bound required",
            "precision":"exact rational coefficients",
            "evidence":"No negative-leading polynomial ray found; absence of these witnesses is not a boundedness proof."}


def certify_marginal_full_stability(p: RealPotential) -> dict[str, Any]:
    """Certify one exact marginal family using the corrected uniform-C test.

    This bounded implementation proves the identity C*J4-J2^2 == 0 when the
    only negative quadratic direction and quartic are proportional to b and
    b^2. Other marginal cases remain INCONCLUSIVE.
    """
    y1,y2,y3=p.Y
    z1,z2,z3,z4,z5,z6,z7=p.Z
    if y1==0 and y3==0 and y2<0 and z2>0 and all(v==0 for v in (z1,z3,z4,z5,z6,z7)):
        C=sp.factor(y2*y2/(z2/2))
        return {"status":"PASS","method":"exact uniform marginal bound identity",
                "observable":"full-potential boundedness at a quartic-flat boundary",
                "tolerance":"C>0 and C*J4-J2^2>=0 wherever J2<0,J4>0",
                "precision":"exact rational identity",
                "evidence":{"uniform_C":_s(C),"J2":"Y2*b","J4":"Z2*b^2/2",
                            "identity":_s(sp.expand(C*(z2*sp.Symbol("b")**2/2)-(y2*sp.Symbol("b"))**2)),
                            "quartic_zero_directions":"J2=0 on all of them",
                            "proof_scope":"exact proportional b-only marginal form"}}
    return {"status":"INCONCLUSIVE","method":"corrected uniform-C marginal criterion",
            "observable":"full-potential boundedness at quartic-flat directions",
            "tolerance":"uniform C>0 bound and nonnegative J2 on J4=0 required",
            "precision":"exact rational coefficient inspection",
            "evidence":"No implemented exact uniform-C identity matches this marginal potential."}


def _target_stationarity(p: RealPotential, target: tuple[Any,Any,Any,Any]) -> tuple[bool,str]:
    u,b,x,y=map(_q,target)
    if not _bilinear_domain(u,b,x,y):
        return False,"EW target is outside the exact bilinear domain."
    y1,y2,y3=p.Y;z1,z2,z3,z4,z5,z6,z7=p.Z
    du=sp.factor(y1+z1*u+z3*b+2*z6*x)
    db=sp.factor(y2+z2*b+z3*u+2*z7*x)
    dx=sp.factor(2*y3+2*(z4+z5)*x+2*z6*u+2*z7*b)
    dy=sp.factor(2*(z4-z5)*y)
    g=sp.factor(u*b-x*x-y*y)
    if u==0 and b==0:
        return True,"The origin is a stationary point in the original doublet fields."
    if b==0 and u>0:
        ok=(du==0 and dx==0 and x==0 and y==0)
        return ok,"b=0 axis field tadpoles: occupied norm and both cross components."
    if u==0 and b>0:
        ok=(db==0 and dx==0 and x==0 and y==0)
        return ok,"u=0 axis field tadpoles: occupied norm and both cross components."
    if g>0:
        ok=(du==0 and db==0 and dx==0 and dy==0)
        return ok,"Interior field stationarity requires all four bilinear derivatives to vanish."
    # On the neutral boundary, require one exact Lagrange multiplier for g=0.
    components=((du,b),(db,u),(dx,-2*x),(dy,-2*y))
    mus=[]
    for gradient,normal in components:
        if normal!=0:
            mus.append(sp.factor(gradient/normal))
        elif gradient!=0:
            return False,"Neutral-boundary field tadpoles are not proportional to the bilinear-domain normal."
    if not mus:
        ok=(du==db==dx==dy==0)
    else:
        ok=all(sp.simplify(mu-mus[0])==0 for mu in mus[1:])
    return ok,"Neutral-boundary KKT stationarity with an exact common multiplier."


def _aligned_global_certificate(p: RealPotential,target:tuple[Any,Any,Any,Any]) -> dict[str,Any] | None:
    y1,y2,y3=p.Y;z1,z2,z3,z4,z5,z6,z7=p.Z
    c=z3+z4+z5;d=z4+z5;q=z4-z5
    margin=sp.factor(z2/2-z7*z7/c) if c!=0 else None
    uew=sp.factor(-y1/z1) if z1!=0 else None
    conditions=(y2==0 and y3==0 and z6==0 and z1>0 and c>0 and z2>0
                and d<=0 and d<=q and margin is not None and margin>=0
                and uew is not None and uew>0)
    if not conditions or tuple(map(_q,target))!=(uew,sp.S.Zero,sp.S.Zero,sp.S.Zero):
        return None
    return {"status":"PASS","method":"exact aligned full-domain potential decomposition",
            "observable":"EW target global minimum including neutral, complex and charge-breaking bilinears",
            "tolerance":"exact V-V_EW nonnegative on full bilinear domain; equality allowed",
            "precision":"exact rational arithmetic",
            "evidence":{"certificate":"Z1>0, C>0, Z2>0, Z4+Z5<=min(0,Z4-Z5), Z2/2-Z7^2/C>=0",
                       "C":_s(c),"quadratic_margin":_s(margin),"u_EW":_s(uew),
                       "decomposition":"V-V_EW >= Z1/2*(u-u_EW)^2 + Z2/2*b^2 + C*u*b + 2*Z7*b*x",
                       "bilinear_domain":"all u,b>=0 and x^2+y^2<=u*b"}}


def certify_global_minimum(p: RealPotential, target: tuple[Any,Any,Any,Any]) -> dict[str, Any]:
    """Compare an EW target with exact competing-point witnesses/certificates."""
    stationary,stationary_reason=_target_stationarity(p,target)
    if not stationary:
        return {"status":"FAIL","method":"exact field-stationarity rejection",
                "observable":"EW target global minimum","tolerance":"target must be a stationary field configuration",
                "precision":"exact rational/algebraic comparison",
                "evidence":{"target":list(map(_s,target)),"reason":stationary_reason,
                           "stationary_enumeration":{"status":"NOT_TESTED","reason":"The supplied target is not a stationary field configuration."}}}
    aligned=_aligned_global_certificate(p,target)
    if aligned is not None:
        aligned["evidence"]["stationary_enumeration"]={
            "status":"NOT_TESTED","reason":"An exact full-domain aligned potential decomposition certifies ordering without enumerating auxiliary stationary candidates."}
        return aligned
    bfb=certify_bfb(p)
    if bfb["status"]=="FAIL":
        witness=bfb["evidence"]
        n=sp.symbols("n",positive=True)
        ray=tuple(n*sp.Rational(witness[key]) for key in ("u","b","x","y"))
        poly=sp.Poly(sp.expand(potential_value(p,*ray)),n,domain=sp.QQ)
        return {"status":"FAIL","method":"exact negative quartic ray implies unbounded full potential",
                "observable":"global minimum of the full tree potential",
                "tolerance":"full potential tends to minus infinity on an exact feasible ray",
                "precision":"exact rational polynomial in the scaling parameter",
                "evidence":{"quartic_witness":witness,"ray":"(u,b,x,y)=n*(witness)",
                            "V_of_n":_s(poly.as_expr()),"leading_power":poly.degree(),
                            "leading_coefficient":_s(poly.LC()),
                            "stationary_enumeration":{"status":"NOT_TESTED","reason":"A negative homogeneous quartic witness proves unboundedness."}}}
    runaway=certify_runaway(p)
    if runaway["status"]=="FAIL":
        return {"status":"FAIL","method":"exact runaway ray",
                "observable":"global minimum of the full tree potential",
                "tolerance":"potential must be bounded below",
                "precision":"exact rational polynomial along a verified domain ray",
                "evidence":{"runaway_witness":runaway["evidence"],
                           "stationary_enumeration":{"status":"NOT_TESTED","reason":"An exact runaway ray proves no global minimum exists."}}}
    target_energy=potential_value(p,*target)
    enumeration=enumerate_stationary_candidates(p)
    marginal=certify_marginal_full_stability(p) if not (bfb["status"]=="PASS" and bfb["evidence"].get("strict_certificate")) else None
    lower=[]
    unresolved_order=[]
    for candidate in enumeration["candidates"]:
        energy=sp.sympify(candidate["V_GeV4"])
        comparison=_exact_sign(energy-target_energy)
        if comparison == -1:
            lower.append({**candidate,"difference_from_target_GeV4":_s(sp.simplify(energy-target_energy))})
        elif comparison is None:
            unresolved_order.append(candidate)
    if lower:
        return {"status":"FAIL","method":"exact feasible competing-point witness",
                "observable":"EW target is not a global minimum of the full tree potential",
                "tolerance":"strict lower potential value",
                "precision":"exact rational/algebraic comparison",
                "evidence":{"target":list(map(_s,target)),"target_V_GeV4":_s(target_energy),
                           "lower_points":lower,"stationary_enumeration":enumeration}}
    stable=(bfb["status"]=="PASS" and bfb["evidence"].get("strict_certificate") is True)
    stable=stable or (marginal is not None and marginal["status"]=="PASS")
    status="PASS" if enumeration["complete"] and stable and not unresolved_order else "INCONCLUSIVE"
    return {"status":status,"method":"complete exact stationary enumeration and potential ordering",
            "observable":"global minimum of the full tree potential",
            "tolerance":"certified V >= V_target over the full bilinear domain",
            "precision":"exact rational/algebraic comparisons",
            "evidence":{"target":list(map(_s,target)),"target_V_GeV4":_s(target_energy),
                       "stationary_enumeration":enumeration,
                       "boundedness":{"strict_quartic":bfb,"corrected_marginal":marginal},
                       "runaway_check":runaway,
                       "unresolved_energy_comparisons":unresolved_order,
                       "reason":("All strata were completely enumerated, boundedness was certified, and exact energy orderings passed." if status=="PASS"
                                 else "No lower witness was found, but stationary completeness, boundedness, or exact energy ordering remains unresolved.")}}


def rounded_input_potential(inputs: Mapping[str, Any], *, v2_sm: Any | None = None) -> dict[str, Any]:
    """Reconstruct the exact real model represented by binary64 input lexemes.

    Floats are first rounded exactly as Python/C++ binary64 values and then
    converted with Fraction.from_float. The result certifies this rationalized
    input model, not every intermediate operation inside the backend.
    """
    def binary64(name: str) -> sp.Rational:
        value = float(inputs[name])
        frac = Fraction.from_float(value)
        return sp.Rational(frac.numerator, frac.denominator)

    l = [binary64(f"lambda{i}") for i in range(1,8)]
    m12 = binary64("m12_sq_GeV2")
    t = binary64("tan_beta")
    gf_value=None
    if v2_sm is None:
        v_backend=float(inputs["v_GeV"])
        # Match the actual SM interface call GF=1/(sqrt(2)*v*v), then the
        # backend getter get_v2()=1/(sqrt(2)*GF), each rounded to binary64.
        gf_value=float(1.0/(math.sqrt(2.0)*v_backend*v_backend))
        v2_backend=float(1.0/(math.sqrt(2.0)*gf_value))
        frac=Fraction.from_float(v2_backend)
        v2=sp.Rational(frac.numerator,frac.denominator)
    else:
        if isinstance(v2_sm, float):
            frac=Fraction.from_float(v2_sm)
            v2=sp.Rational(frac.numerator,frac.denominator)
        else:
            v2 = _q(v2_sm)
    l1,l2,l3,l4,l5,l6,l7 = l
    lam345 = l3+l4+l5
    c2, s2, sc = 1/(1+t*t), t*t/(1+t*t), t/(1+t*t)
    z1 = l1*c2*c2+l2*s2*s2+2*lam345*sc*sc+4*l6*c2*sc+4*l7*s2*sc
    z2 = l1*s2*s2+l2*c2*c2+2*lam345*sc*sc-4*l6*sc*s2-4*l7*sc*c2
    z3 = l3+(l1+l2-2*lam345)*sc*sc+2*(l7-l6)*sc*(c2-s2)
    d = 2*(l1+l2-2*l3)*sc*sc+(l4+l5)*(c2-s2)**2+4*(l7-l6)*sc*(c2-s2)
    e = l4-l5
    z4,z5 = (d+e)/2,(d-e)/2
    z6 = -l1*c2*sc+l2*s2*sc+lam345*sc*(c2-s2)+l6*(c2*c2-3*sc*sc)+l7*(3*sc*sc-s2*s2)
    z7 = -l1*s2*sc+l2*c2*sc-lam345*sc*(c2-s2)+l6*(3*sc*sc-s2*s2)+l7*(c2*c2-3*sc*sc)
    m11 = m12*t-v2/2*(l1*c2+lam345*s2+3*l6*sc+l7*s2*t)
    m22 = m12/t-v2/2*(l2*s2+lam345*c2+l6*c2/t+3*l7*sc)
    y1 = c2*m11+s2*m22-2*sc*m12
    y2 = s2*m11+c2*m22+2*sc*m12
    y3 = sc*(m22-m11)-(c2-s2)*m12
    p = RealPotential((y1,y2,y3),(z1,z2,z3,z4,z5,z6,z7),v2)
    hessian = sp.Matrix([[z1*v2,z6*v2],
                         [z6*v2,y2+(z3+z4+z5)*v2/2]])
    determinant = sp.factor(hessian.det())
    masses = {"hessian_11": z1*v2,
              "hessian_determinant": determinant,
              "A_squared": y2+(z3+z4-z5)*v2/2,
              "Hplus_squared": y2+z3*v2/2}
    stationary = {"Y1_plus_Z1_v2_over_2": sp.factor(y1+z1*v2/2),
                  "Y3_plus_Z6_v2_over_2": sp.factor(y3+z6*v2/2)}
    positive = all(_exact_sign(v) == 1 for v in masses.values())
    stationarity = all(v == 0 for v in stationary.values())
    return {"status": "PASS" if positive and stationarity else "FAIL",
            "method": "exact reconstruction after binary64 input rounding",
            "observable": "rounded-input Higgs-basis potential, stationarity, and physical Hessian",
            "tolerance": "strict exact positivity; exact stationarity identities",
            "precision": "exact rationals from Fraction.from_float",
            "scope": "rationalized binary64 parameters and supplied SM v^2; not a proof of internal backend operation rounding",
            "inputs_binary64_rational": {**{f"lambda{i}": _s(v) for i,v in enumerate(l,1)},
                                         "m12_sq_GeV2": _s(m12),"tan_beta": _s(t),"v2_GeV2": _s(v2),
                                         "SM_GF_binary64":repr(gf_value) if gf_value is not None else None},
            "potential": p.as_dict(), "stationarity_residuals": {k:_s(v) for k,v in stationary.items()},
            "physical_hessian_and_masses": {k:_s(v) for k,v in masses.items()},
            "physical_positive": positive, "stationary": stationarity}
