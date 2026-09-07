"""Independent algebra audit for appendix theory Phase 5.

This script does not define the physics convention. It only checks the
matrix inversion and exact-alignment substitutions used in the hand derivation.
"""
import sympy as sp

alpha, beta = sp.symbols("alpha beta", real=True)

R = sp.Matrix([
    [-sp.sin(alpha), sp.cos(alpha)],
    [ sp.cos(alpha), sp.sin(alpha)],
])

assert sp.simplify(R.T * R - sp.eye(2)) == sp.zeros(2)
assert sp.simplify(R.T - R) == sp.zeros(2)
assert sp.simplify(R.inv() - R) == sp.zeros(2)

kappa_h = sp.cos(alpha) / sp.sin(beta)
kappa_phi = sp.sin(alpha) / sp.sin(beta)

alignment = {alpha: beta - sp.pi/2}

assert sp.simplify(kappa_h.subs(alignment) - 1) == 0
assert sp.simplify(kappa_phi.subs(alignment) + sp.cot(beta)) == 0

print("PHASE5_YUKAWA_CHECK=PASS")
print("R_alpha inverse = R_alpha")
print("kappa_h(alignment) =", sp.simplify(kappa_h.subs(alignment)))
print("kappa_phi(alignment) =", sp.simplify(kappa_phi.subs(alignment)))
