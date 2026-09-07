import sympy as sp

# Symbols
v,m,cw = sp.symbols('v m cw', positive=True, nonzero=True)
CV,CL,F0,I1 = sp.symbols('CV CL F0 I1')
I = sp.I

# Phase-9 relation between potential and literal Lagrangian coefficients
assert sp.simplify(CL + CV) == CL + CV  # declaration only
CL_sub = -CV

# gamma gamma: code consumes the Feynman-rule object -i CV
FR = -I*CV
S_gaga = sp.simplify(FR * v/(2*m**2) * F0)
Ahat_gaga = sp.simplify(S_gaga/(-I))
assert sp.simplify(Ahat_gaga - CV*v/(2*m**2)*F0) == 0
assert sp.simplify(Ahat_gaga.subs(CV,-CL) + CL*v/(2*m**2)*F0) == 0

# Z gamma: active 2HDMC line
KZ = 2*cw - 1/cw
S_Zga = sp.simplify(-KZ * FR * v/(2*m**2) * I1)
Ahat_Zga = sp.simplify(S_Zga/(-I))
assert sp.simplify(Ahat_Zga + KZ*CV*v/(2*m**2)*I1) == 0
assert sp.simplify(Ahat_Zga.subs(CV,-CL) - KZ*CL*v/(2*m**2)*I1) == 0

# Large-tan(beta) exact Z7 from Phase 7, lambda7=0
t,X,l1,l2,l345,l6,u = sp.symbols('t X l1 l2 l345 l6 u', nonzero=True)
Z7 = -(l6*t**4 + (l1-l345)*t**3 - 3*l6*t**2 + (l345-l2)*t)/(1+t**2)**2

# Fixed-lambda6 and fixed-X expansions
fixed_l6 = sp.series(Z7.subs(t,1/u), u, 0, 4)
fixed_X = sp.series(Z7.subs(l6,X/t).subs(t,1/u), u, 0, 5)

# Exact-alignment relation from Z6=0, lambda7=0
Z6 = (-l1*t + l2*t**3 - l345*t**3 + l345*t - 3*l6*t**2 + l6)/(1+t**2)**2
l345_align = sp.solve(sp.Eq(sp.factor(sp.together(Z6*(1+t**2)**2)),0), l345)[0]
Z7_align = sp.factor(sp.simplify(Z7.subs(l345,l345_align)))
expected_align = -(t*(l1-l2) + l6*(t**2-1))/(1+t**2)
assert sp.simplify(Z7_align-expected_align) == 0

fixedX_align = sp.series(Z7_align.subs(l6,X/t).subs(t,1/u), u, 0, 5)

print('PHASE10_LOOP_CONVENTION_CHECK=PASS')
print('Ahat_gammagamma_Hp =', Ahat_gaga)
print('Ahat_Zgamma_Hp_2HDMC =', Ahat_Zga)
print('fixed lambda6:', fixed_l6)
print('fixed X:', fixed_X)
print('alignment exact Z7:', Z7_align)
print('alignment fixed X:', fixedX_align)
