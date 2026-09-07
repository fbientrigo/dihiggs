import math

# Phase 11 exact fixed-X relation and Phase 13 physical h-phi-phi spot-check.
# This is an audit script: it does not choose the formulas it checks.

mh = 125.13
mphi = 150.0
t = 300000.0
lambda6 = 1.0e-10
X = lambda6*t
M2 = 2.24999999995003345e4
GF = 1.16637e-5
v = 1.0 / math.sqrt(math.sqrt(2.0)*GF)

Q = (mphi*mphi-M2)*t*t
Z7_XQ = (X/2.0-Q/(v*v))/t + (X/2.0+Q/(v*v))/(t**3)
Z7_historical = -X/t

# Exact physical h-phi-phi formula in exact alignment.
g_analytic = (mh*mh + 2.0*mphi*mphi - 2.0*M2)/v
g_2hdmc = 63.5914252007596588

# Symbolic identities represented numerically on the versioned benchmark.
assert abs(Z7_XQ + 2.472539016945184e-6) < 1e-18
assert abs(g_analytic-g_2hdmc) < 1e-8

print('PHASE11_13_FREEZE_CHECK=PASS')
print('Q [GeV^2] =', Q)
print('Q/v^2 =', Q/(v*v))
print('X =', X)
print('Z7 exact X,Q =', Z7_XQ)
print('historical -X/t =', Z7_historical)
print('ratio exact/historical =', Z7_XQ/Z7_historical)
print('g_hphiphi analytic [GeV] =', g_analytic)
print('g_hphiphi 2HDMC [GeV] =', g_2hdmc)
