#!/usr/bin/env python3
"""Algebra audit for issue #81 Phase 5.

This script does not define the Yukawa convention. It checks only algebra that
follows from the already-frozen Phase-4 rotation and the project Type-I choice
that Phi2 is Yukawa-active.
"""
import sympy as sp

sa, ca = sp.symbols('s_alpha c_alpha', real=True)
sb, cb = sp.symbols('s_beta c_beta', nonzero=True, real=True)
sba, cba = sp.symbols('s_ba c_ba', real=True)

# Phase-4 CP-even rotation.
R = sp.Matrix([[-sa, ca], [ca, sa]])
# Orthogonality uses sa^2+ca^2=1. Check the explicit inverse formula by hand:
rho1_from_inverse = -sa*sp.Symbol('h') + ca*sp.Symbol('phi')
rho2_from_inverse = ca*sp.Symbol('h') + sa*sp.Symbol('phi')

# Trigonometric expansion for alpha = beta-(beta-alpha).
sa_from_ba = sb*cba - cb*sba
ca_from_ba = cb*cba + sb*sba

kappa_h = sp.expand(ca_from_ba/sb)
kappa_phi = sp.expand(sa_from_ba/sb)

assert sp.simplify(kappa_h - (sba + (cb/sb)*cba)) == 0
assert sp.simplify(kappa_phi - (cba - (cb/sb)*sba)) == 0

# Exact project alignment: sin(beta-alpha)=+1, cos(beta-alpha)=0.
assert sp.simplify(kappa_h.subs({sba:1, cba:0}) - 1) == 0
assert sp.simplify(kappa_phi.subs({sba:1, cba:0}) + cb/sb) == 0

print('rho1 =', rho1_from_inverse)
print('rho2 =', rho2_from_inverse)
print('kappa_h =', kappa_h)
print('kappa_phi =', kappa_phi)
print('alignment kappa_h =', kappa_h.subs({sba:1, cba:0}))
print('alignment kappa_phi =', kappa_phi.subs({sba:1, cba:0}))
print('PASS')
