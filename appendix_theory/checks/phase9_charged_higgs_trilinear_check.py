import sympy as sp
v=sp.symbols('v', nonzero=True, real=True)
Z3,Z7=sp.symbols('Z3 Z7', real=True)
rv,rp,G0,A,Hp,Hm=sp.symbols('rv rp G0 A Hp Hm', real=True)
Y1,Y2,Y3,Z1,Z2,Z4,Z5,Z6=sp.symbols('Y1 Y2 Y3 Z1 Z2 Z4 Z5 Z6', real=True)
I=sp.I
A1=((v+rv)**2+G0**2)/2
A2=Hm*Hp+(rp**2+A**2)/2
C=(v+rv-I*G0)*(rp+I*A)/2
Cd=(v+rv+I*G0)*(rp-I*A)/2
V=(Y1*A1+Y2*A2+Y3*(C+Cd)+sp.Rational(1,2)*Z1*A1**2+
   sp.Rational(1,2)*Z2*A2**2+Z3*A1*A2+Z4*C*Cd+
   sp.Rational(1,2)*Z5*(C**2+Cd**2)+Z6*A1*(C+Cd)+Z7*A2*(C+Cd))
zero={rv:0,rp:0,G0:0,A:0,Hp:0,Hm:0}
assert sp.simplify(sp.diff(V,rv,Hp,Hm).subs(zero)-v*Z3)==0
assert sp.simplify(sp.diff(V,rp,Hp,Hm).subs(zero)-v*Z7)==0
s,c,h,phi=sp.symbols('s c h phi', real=True)
Vtri=sp.expand(v*(Z3*(s*h+c*phi)+Z7*(c*h-s*phi))*Hp*Hm)
assert sp.simplify(sp.diff(Vtri,phi,Hp,Hm)-v*(Z3*c-Z7*s))==0
assert sp.simplify(sp.diff(Vtri,h,Hp,Hm)-v*(Z3*s+Z7*c))==0
print('PHASE9_CHARGED_HIGGS_TRILINEAR_CHECK=PASS')
print('d3V_rhov_Hp_Hm =', v*Z3)
print('d3V_rhoperp_Hp_Hm =', v*Z7)
print('C_V(phi H+ H-) =', v*(Z3*c-Z7*s))
print('alignment C_V =', -v*Z7)
print('alignment C_L =', v*Z7)
print('alignment Feynman rule = + i v Z7')
