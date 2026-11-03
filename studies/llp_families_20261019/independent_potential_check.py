"""Independent check (no project code): rebuild every published point from the generic-basis potential with sympy/mpmath
at 60 digits -> masses, h-phi mixing, M^2 - m_h^2/2 and g_hphiphi = d^3V/dh dphi^2.  Usage: python -I scripts/independent_potential_check.py .
Uses only validated_benchmarks.csv."""
import csv, sys
import sympy as sp
D = 60
H = sys.argv[1]
f = sp.symbols('r1 r2 e1 e2 a1 b1 a2 b2', real=True)
r1, r2, e1, e2, a1, b1, a2, b2 = f
L = sp.symbols('l1:8', real=True); m11, m22, m12, v1, v2 = sp.symbols('m11 m22 m12 v1 v2', real=True)
P1 = sp.Matrix([(a1 + sp.I*b1)/sp.sqrt(2), (v1 + r1 + sp.I*e1)/sp.sqrt(2)])
P2 = sp.Matrix([(a2 + sp.I*b2)/sp.sqrt(2), (v2 + r2 + sp.I*e2)/sp.sqrt(2)])
d = lambda A, B: (A.H*B)[0]
x11, x22, x12 = d(P1, P1), d(P2, P2), d(P1, P2); x21 = sp.conjugate(x12)
V = (m11*x11 + m22*x22 - m12*(x12 + x21) + L[0]/2*x11**2 + L[1]/2*x22**2 + L[2]*x11*x22 + L[3]*x12*x21
     + L[4]/2*(x12**2 + x21**2) + (L[5]*x11 + L[6]*x22)*(x12 + x21))
V = sp.expand(V)
zero = {k: 0 for k in f}
tad = [sp.diff(V, r1).subs(zero), sp.diff(V, r2).subs(zero)]
sol = sp.solve(tad, [m11, m22], dict=True)[0]
Vs = V.subs(sol)
hess = sp.hessian(Vs, f).subs(zero)
out = []
for r in csv.DictReader(open(H + '/validated_benchmarks.csv')):
    t = sp.Float(r['tan_beta'], D); v = sp.Float('246', D)
    c = 1/sp.sqrt(1 + t**2); s = t*c
    num = {v1: v*c, v2: v*s, m12: sp.Float(r['m12_sq_GeV2'], D)}
    num.update({L[i]: sp.Float(r[f'lambda{i+1}'], D) for i in range(7)})
    M = hess.subs(num).evalf(D)
    import mpmath; mpmath.mp.dps = D
    E, _ = mpmath.eigsy(mpmath.matrix([[mpmath.mpf(str(sp.N(M[i, j], D))) for j in range(8)] for i in range(8)]))
    ev = sorted(sp.Float(str(z), D) for z in E)
    # neutral CP-even block in (h, phi) Higgs-basis rotation
    Mr = sp.Matrix([[M[0, 0], M[0, 1]], [M[1, 0], M[1, 1]]])
    R = sp.Matrix([[c, s], [-s, c]])
    Mh = R*Mr*R.T
    # cubic h-phi-phi: d^3 V / dh dphi^2 with r1 = c h - s phi, r2 = s h + c phi
    h, p = sp.symbols('h p')
    V3 = Vs.subs({e1: 0, e2: 0, a1: 0, b1: 0, a2: 0, b2: 0}).subs({r1: c*h - s*p, r2: s*h + c*p}).subs(num)
    g = sp.diff(V3, h, p, p).subs({h: 0, p: 0})
    M2 = num[m12]/(s*c)
    mphi = sp.Float(r['mphi_GeV'], D)
    out.append(dict(id=r['candidate_id'], m=r['mphi_GeV'],
        mh=sp.sqrt(Mh[0, 0]), mphi=sp.sqrt(Mh[1, 1]), mix=Mh[0, 1],
        nonzero_masses=[sp.sqrt(z) for z in ev if abs(z) > 1e-20],
        M2_minus_mh2_over_2=M2 - sp.Float('125.13', D)**2/2,
        g=g, g_minus_2m2_over_v=sp.Abs(g) - 2*mphi**2/v, sign_g=sp.sign(g)))
for o in out:
    print(o['id'], o['m'], 'mh=%s mphi=%s mix=%s' % (sp.N(o['mh'], 15), sp.N(o['mphi'], 15), sp.N(o['mix'], 3)))
    print('   masses:', [sp.N(z, 12) for z in o['nonzero_masses']])
    print('   M2-mh2/2=%s  g=%s  |g|-2m2/v=%s' % (sp.N(o['M2_minus_mh2_over_2'], 3), sp.N(o['g'], 12), sp.N(o['g_minus_2m2_over_v'], 3)))
