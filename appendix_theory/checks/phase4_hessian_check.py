#!/usr/bin/env python3
"""CAS audit for issue #81 Phase 4.

This script does not define the theory. It expands the already-frozen
DH05/BFLRS11 CP-conserving potential and checks the hand derivation.
Charged conjugate fields are treated as algebraically independent variables.
"""
import sympy as sp

v1, v2 = sp.symbols('v1 v2', nonzero=True, real=True)
m11, m22, m12 = sp.symbols('m11sq m22sq m12sq', real=True)
l1,l2,l3,l4,l5,l6,l7 = sp.symbols(
    'lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 lambda7', real=True)
r1,r2,e1,e2 = sp.symbols('rho1 rho2 eta1 eta2', real=True)
p1m,p1p,p2m,p2p = sp.symbols('phi1m phi1p phi2m phi2p')
i = sp.I

x11 = p1m*p1p + ((v1+r1)**2 + e1**2)/2
x22 = p2m*p2p + ((v2+r2)**2 + e2**2)/2
x12 = p1m*p2p + (v1+r1-i*e1)*(v2+r2+i*e2)/2
x21 = p2m*p1p + (v2+r2-i*e2)*(v1+r1+i*e1)/2

V = (m11*x11 + m22*x22 - m12*(x12+x21)
     + l1*x11**2/2 + l2*x22**2/2 + l3*x11*x22 + l4*x12*x21
     + l5*(x12**2+x21**2)/2
     + l6*x11*(x12+x21) + l7*x22*(x12+x21))

zero = {r1:0,r2:0,e1:0,e2:0,p1m:0,p1p:0,p2m:0,p2p:0}
Mpm = sp.Matrix([[sp.diff(V, mi, pj).subs(zero)
                  for pj in (p1p,p2p)] for mi in (p1m,p2m)])
MA = sp.hessian(V, (e1,e2)).subs(zero)
Mrho = sp.hessian(V, (r1,r2)).subs(zero)

l345 = l3+l4+l5
m11_tad = m12*v2/v1 - (l1*v1**2+l345*v2**2+3*l6*v1*v2+l7*v2**3/v1)/2
m22_tad = m12*v1/v2 - (l2*v2**2+l345*v1**2+l6*v1**3/v2+3*l7*v1*v2)/2
tad = {m11:m11_tad, m22:m22_tad}
Mpm = sp.simplify(Mpm.subs(tad))
MA = sp.simplify(MA.subs(tad))
Mrho = sp.simplify(Mrho.subs(tad))

Dpm = m12 - ((l4+l5)*v1*v2 + l6*v1**2 + l7*v2**2)/2
DA = m12 - l5*v1*v2 - (l6*v1**2+l7*v2**2)/2
G = sp.Matrix([[v2/v1,-1],[-1,v1/v2]])
vev = sp.Matrix([v1,v2])

assert sp.simplify(Mpm-Dpm*G) == sp.zeros(2)
assert sp.simplify(MA-DA*G) == sp.zeros(2)
assert sp.simplify(Mpm*vev) == sp.zeros(2,1)
assert sp.simplify(MA*vev) == sp.zeros(2,1)

print('charged factor:', sp.factor(Dpm))
print('CP-odd factor:', sp.factor(DA))
print('CP-even Hessian after tadpoles:')
sp.print_latex(Mrho)
print('PASS')
