import sympy as sp

beta, alpha = sp.symbols('beta alpha', real=True)
sb, cb = sp.sin(beta), sp.cos(beta)
sa, ca = sp.sin(alpha), sp.cos(alpha)

# Phase-4 inverse rotation
# rho1 = -sa*h + ca*phi
# rho2 =  ca*h + sa*phi
kh = sp.simplify(-cb*sa + sb*ca)
kphi = sp.simplify(cb*ca + sb*sa)

assert sp.simplify(kh - sp.sin(beta-alpha)) == 0
assert sp.simplify(kphi - sp.cos(beta-alpha)) == 0

# Exact-alignment branch alpha = beta - pi/2
align = {alpha: beta-sp.pi/2}
assert sp.simplify(kh.subs(align)) == 1
assert sp.simplify(kphi.subs(align)) == 0

# Orthogonality / quartic norm check
h, phi = sp.symbols('h phi', real=True)
rho1 = -sa*h + ca*phi
rho2 = ca*h + sa*phi
assert sp.simplify(rho1**2 + rho2**2 - h**2 - phi**2) == 0

# Gauge-mass expansion: sum_i (vi+rho_i)^2
v1, v2, r1, r2 = sp.symbols('v1 v2 r1 r2', real=True)
expr = sp.expand((v1+r1)**2 + (v2+r2)**2)
linear = sp.expand(expr - (v1**2+v2**2) - (r1**2+r2**2))
assert sp.simplify(linear - 2*(v1*r1+v2*r2)) == 0

print('PHASE6_GAUGE_CHECK=PASS')
print('kappa_V_h =', kh)
print('kappa_V_phi =', kphi)
print('alignment:', sp.simplify(kh.subs(align)), sp.simplify(kphi.subs(align)))
print('quartic_norm = h^2 + phi^2')
