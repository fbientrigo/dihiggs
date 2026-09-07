import sympy as sp

c,s,v = sp.symbols('c s v', real=True, nonzero=True)
A,B,C,D = sp.symbols('A B C D', commutative=True)
m11,m22,m12 = sp.symbols('m11 m22 m12', real=True)
l1,l2,l3,l4,l5,l6,l7 = sp.symbols('l1 l2 l3 l4 l5 l6 l7', real=True)

# Inverse Higgs-basis rotation:
# Phi1 = c H1 - s H2, Phi2 = s H1 + c H2.
X1 = c**2*A + s**2*B - c*s*(C+D)
X2 = s**2*A + c**2*B + c*s*(C+D)
X12 = c*s*(A-B) + c**2*C - s**2*D
X21 = c*s*(A-B) - s**2*C + c**2*D

V2 = m11*X1 + m22*X2 - m12*(X12+X21)
V4 = (sp.Rational(1,2)*l1*X1**2 + sp.Rational(1,2)*l2*X2**2
      + l3*X1*X2 + l4*X12*X21
      + sp.Rational(1,2)*l5*(X12**2+X21**2)
      + l6*X1*(X12+X21) + l7*X2*(X12+X21))
V2 = sp.expand(V2)
V4 = sp.expand(V4)

coef = lambda e, mon: sp.expand(e).coeff(mon)
Y1 = coef(V2,A)
Y2 = coef(V2,B)
Y3 = coef(V2,C)
Z1 = 2*coef(V4,A**2)
Z2 = 2*coef(V4,B**2)
Z3 = coef(V4,A*B)
Z4 = coef(V4,C*D)
Z5 = 2*coef(V4,C**2)
Z6 = coef(V4,A*C)
Z7 = coef(V4,B*C)

# Direct Higgs-basis tadpoles and Hessian.
r1,r2,aa,gg = sp.symbols('r1 r2 aa gg', real=True)
I = sp.I
h10=(v+r1+I*gg)/sp.sqrt(2)
h20=(r2+I*aa)/sp.sqrt(2)
Ah=sp.expand(sp.conjugate(h10)*h10)
Bh=sp.expand(sp.conjugate(h20)*h20)
Ch=sp.expand(sp.conjugate(h10)*h20)
Dh=sp.conjugate(Ch)
y1,y2,y3,z1,z2,z3,z4,z5,z6,z7 = sp.symbols('y1 y2 y3 z1 z2 z3 z4 z5 z6 z7', real=True)
VH=(y1*Ah+y2*Bh+y3*(Ch+Dh)+sp.Rational(1,2)*z1*Ah**2
    +sp.Rational(1,2)*z2*Bh**2+z3*Ah*Bh+z4*Ch*Dh
    +sp.Rational(1,2)*z5*(Ch**2+Dh**2)+z6*Ah*(Ch+Dh)+z7*Bh*(Ch+Dh))
origin={r1:0,r2:0,aa:0,gg:0}
T1=sp.simplify(sp.diff(VH,r1).subs(origin))
T2=sp.simplify(sp.diff(VH,r2).subs(origin))
station={y1:-z1*v**2/2,y3:-z6*v**2/2}
M11=sp.simplify(sp.diff(VH,r1,2).subs(origin).subs(station))
M12=sp.simplify(sp.diff(VH,r1,r2).subs(origin).subs(station))
M22=sp.simplify(sp.diff(VH,r2,2).subs(origin).subs(station))
MA=sp.simplify(sp.diff(VH,aa,2).subs(origin).subs(station))

assert sp.simplify(T1-v*(y1+z1*v**2/2)) == 0
assert sp.simplify(T2-v*(y3+z6*v**2/2)) == 0
assert sp.simplify(M11-z1*v**2) == 0
assert sp.simplify(M12-z6*v**2) == 0
assert sp.simplify(M22-(y2+(z3+z4+z5)*v**2/2)) == 0
assert sp.simplify(MA-(y2+(z3+z4-z5)*v**2/2)) == 0

# Generic tadpoles imply Higgs-basis stationarity.
l345=l3+l4+l5
m11_t=m12*s/c-v**2*sp.Rational(1,2)*(l1*c**2+l345*s**2+3*l6*s*c+l7*s**3/c)
m22_t=m12*c/s-v**2*sp.Rational(1,2)*(l2*s**2+l345*c**2+l6*c**3/s+3*l7*s*c)
assert sp.simplify(Y1.subs({m11:m11_t,m22:m22_t}) + Z1*v**2/2) == 0
assert sp.simplify(Y3.subs({m11:m11_t,m22:m22_t}) + Z6*v**2/2) == 0

# Match compact DH/2HDMC forms.
s2b=2*s*c
c2b=c**2-s**2
c3b=4*c**3-3*c
s3b=3*s-4*s**3
Z6_ref=-sp.Rational(1,2)*s2b*(l1*c**2-l2*s**2-l345*c2b)+c*c3b*l6+s*s3b*l7
Z7_ref=-sp.Rational(1,2)*s2b*(l1*s**2-l2*c**2+l345*c2b)+s*s3b*l6+c*c3b*l7
assert sp.simplify(sp.expand(Z6-Z6_ref).subs(s**2,1-c**2)) == 0
assert sp.simplify(sp.expand(Z7-Z7_ref).subs(s**2,1-c**2)) == 0

# Physical mass relation.
S,Cba,mh,mphi=sp.symbols('S Cba mh mphi', real=True)
R=sp.Matrix([[S,Cba],[Cba,-S]])
Mb=sp.simplify(R.T*sp.diag(mh**2,mphi**2)*R)
assert sp.simplify(Mb[0,1]-(mh**2-mphi**2)*S*Cba) == 0

# Large tan(beta), lambda7=0.
t=sp.symbols('t', positive=True)
ct=1/sp.sqrt(1+t**2); st=t/sp.sqrt(1+t**2)
Z7_t=sp.factor(Z7.subs({c:ct,s:st,l7:0}))
assert sp.simplify(sp.limit(Z7_t,t,sp.oo)+l6) == 0

print('PHASE7_HIGGS_BASIS_CHECK=PASS')
print('Y3 =', sp.factor(Y3))
print('Higgs-basis tadpoles:', T1, ',', T2)
print('CP-even offdiag =', M12)
print('Z6*v^2 physical relation = (mh^2-mphi^2) sba cba')
print('large-tanbeta lambda7=0: Z7 -> -lambda6')
