# 03 — Vacuum and minimization

Scope: Phase 3 of issue #81.  
Input convention: `02_SCALAR_POTENTIAL.md`.  
Analytic anchor for post-derivation comparison: DH05 `hbasis.tex`, Eqs. `minconditionsa`, `minconditionsb` specialized to a real CP-conserving vacuum.

The equations below are derived first from the potential. Source equations and 2HDMC are used only afterward as independent checks.

## What is established

[DERIVED] For the real neutral vacuum

\[
\Phi_1\to\frac1{\sqrt2}\binom{0}{v_1},\qquad
\Phi_2\to\frac1{\sqrt2}\binom{0}{v_2},
\]

define

\[
\lambda_{345}\equiv\lambda_3+\lambda_4+\lambda_5.
\]

The vacuum potential is

\[
\boxed{
\begin{aligned}
V_0(v_1,v_2)={}&
\frac12m_{11}^2v_1^2
+\frac12m_{22}^2v_2^2
-m_{12}^2v_1v_2\\
&+\frac18\lambda_1v_1^4
+\frac18\lambda_2v_2^4
+\frac14\lambda_{345}v_1^2v_2^2\\
&+\frac12\lambda_6v_1^3v_2
+\frac12\lambda_7v_1v_2^3.
\end{aligned}}
\]

The two independent neutral stationarity equations are

\[
\boxed{
0=\frac{\partial V_0}{\partial v_1}
=m_{11}^2v_1-m_{12}^2v_2
+\frac12\lambda_1v_1^3
+\frac12\lambda_{345}v_1v_2^2
+\frac32\lambda_6v_1^2v_2
+\frac12\lambda_7v_2^3
}
\]

and

\[
\boxed{
0=\frac{\partial V_0}{\partial v_2}
=m_{22}^2v_2-m_{12}^2v_1
+\frac12\lambda_2v_2^3
+\frac12\lambda_{345}v_1^2v_2
+\frac12\lambda_6v_1^3
+\frac32\lambda_7v_1v_2^2.
}
\]

For nonzero `v1,v2`,

\[
\boxed{
\begin{aligned}
m_{11}^2={}&m_{12}^2\frac{v_2}{v_1}
-\frac12\left[
\lambda_1v_1^2+\lambda_{345}v_2^2
+3\lambda_6v_1v_2
+\lambda_7\frac{v_2^3}{v_1}
\right],\\[1mm]
m_{22}^2={}&m_{12}^2\frac{v_1}{v_2}
-\frac12\left[
\lambda_2v_2^2+\lambda_{345}v_1^2
+\lambda_6\frac{v_1^3}{v_2}
+3\lambda_7v_1v_2
\right].
\end{aligned}}
\]

With `s=sin(beta)`, `c=cos(beta)`,

\[
\boxed{
\begin{aligned}
m_{11}^2={}&m_{12}^2\tan\beta
-\frac{v^2}{2}
\left[
\lambda_1c^2+\lambda_{345}s^2+3\lambda_6sc
+\lambda_7\frac{s^3}{c}
\right],\\[1mm]
m_{22}^2={}&m_{12}^2\cot\beta
-\frac{v^2}{2}
\left[
\lambda_2s^2+\lambda_{345}c^2
+\lambda_6\frac{c^3}{s}
+3\lambda_7sc
\right].
\end{aligned}}
\]

These reproduce the CP-conserving `xi=0` limit of DH05's general source equations.

## Derivation

At the vacuum,

\[
\Phi_1^\dagger\Phi_1=\frac{v_1^2}{2},\qquad
\Phi_2^\dagger\Phi_2=\frac{v_2^2}{2},\qquad
\Phi_1^\dagger\Phi_2=\Phi_2^\dagger\Phi_1=\frac{v_1v_2}{2}.
\]

Substitute term by term.

### Quadratic terms

\[
m_{11}^2\Phi_1^\dagger\Phi_1\to\frac12m_{11}^2v_1^2,
\qquad
m_{22}^2\Phi_2^\dagger\Phi_2\to\frac12m_{22}^2v_2^2.
\]

For the off-diagonal term, because `m12^2` and the VEVs are real,

\[
-[m_{12}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}]
=-2m_{12}^2\frac{v_1v_2}{2}
=-m_{12}^2v_1v_2.
\]

This Hermitian-conjugate factor of two is essential.

### Quartics `lambda1,lambda2`

\[
\frac{\lambda_1}{2}(\Phi_1^\dagger\Phi_1)^2
\to\frac18\lambda_1v_1^4,
\qquad
\frac{\lambda_2}{2}(\Phi_2^\dagger\Phi_2)^2
\to\frac18\lambda_2v_2^4.
\]

### Quartics `lambda3,lambda4,lambda5`

\[
\lambda_3(\Phi_1^\dagger\Phi_1)(\Phi_2^\dagger\Phi_2)
\to\frac14\lambda_3v_1^2v_2^2,
\]

\[
\lambda_4(\Phi_1^\dagger\Phi_2)(\Phi_2^\dagger\Phi_1)
\to\frac14\lambda_4v_1^2v_2^2,
\]

and

\[
\left[\frac{\lambda_5}{2}(\Phi_1^\dagger\Phi_2)^2+\mathrm{h.c.}\right]
\to2\left[\frac{\lambda_5}{2}\frac{v_1^2v_2^2}{4}\right]
=\frac14\lambda_5v_1^2v_2^2.
\]

Hence

\[
\frac14(\lambda_3+\lambda_4+\lambda_5)v_1^2v_2^2
=\frac14\lambda_{345}v_1^2v_2^2.
\]

### Hard-breaking quartics `lambda6,lambda7`

\[
\begin{aligned}
&\left[\lambda_6(\Phi_1^\dagger\Phi_1)(\Phi_1^\dagger\Phi_2)+\mathrm{h.c.}\right]\\
&\quad\to
2\lambda_6\left(\frac{v_1^2}{2}\right)
\left(\frac{v_1v_2}{2}\right)
=\frac12\lambda_6v_1^3v_2,
\end{aligned}
\]

and

\[
\left[\lambda_7(\Phi_2^\dagger\Phi_2)(\Phi_1^\dagger\Phi_2)+\mathrm{h.c.}\right]
\to\frac12\lambda_7v_1v_2^3.
\]

Adding these contributions gives the displayed `V0`.

Now differentiate. For `v1`,

\[
\frac{\partial}{\partial v_1}\frac12m_{11}^2v_1^2=m_{11}^2v_1,
\qquad
\frac{\partial}{\partial v_1}(-m_{12}^2v_1v_2)=-m_{12}^2v_2,
\]

\[
\frac{\partial}{\partial v_1}\frac18\lambda_1v_1^4=\frac12\lambda_1v_1^3,
\qquad
\frac{\partial}{\partial v_1}\frac14\lambda_{345}v_1^2v_2^2
=\frac12\lambda_{345}v_1v_2^2,
\]

\[
\frac{\partial}{\partial v_1}\frac12\lambda_6v_1^3v_2
=\frac32\lambda_6v_1^2v_2,
\qquad
\frac{\partial}{\partial v_1}\frac12\lambda_7v_1v_2^3
=\frac12\lambda_7v_2^3.
\]

The `v2` derivative is analogous. Thus the factors `3 lambda6` and `3 lambda7` in the solved equations arise directly from differentiating the cubic VEV powers.

Because the neutral fields enter as `v_i+rho_i`,

\[
\left.\frac{\partial V}{\partial\rho_i}\right|_{\rho=\eta=\phi^\pm=0}
=
\frac{\partial V_0}{\partial v_i},
\]

so these are exactly the tree-level CP-even tadpole conditions in this basis.

### Introduce `M^2` only now

[PROJECT-DEFINITION][TRANSLATED]

\[
\boxed{M^2\equiv\frac{m_{12}^2}{\sin\beta\cos\beta}},
\qquad
\boxed{m_{12}^2=M^2\sin\beta\cos\beta}.
\]

Substituting gives

\[
\boxed{
\begin{aligned}
m_{11}^2={}&M^2s^2
-\frac{v^2}{2}
\left[\lambda_1c^2+\lambda_{345}s^2+3\lambda_6sc
+\lambda_7\frac{s^3}{c}\right],\\[1mm]
m_{22}^2={}&M^2c^2
-\frac{v^2}{2}
\left[\lambda_2s^2+\lambda_{345}c^2
+\lambda_6\frac{c^3}{s}+3\lambda_7sc\right].
\end{aligned}}
\]

Therefore

\[
\boxed{m_{22}^2\neq m_{12}^2\neq M^2}
\]

in general. `m22^2` is a diagonal quadratic coefficient; `m12^2` is the off-diagonal quadratic coefficient; `M^2` is a derived rescaling of `m12^2` by the VEV orientation.

### CAS audit

The symbolic audit used exactly

```python
V0 = (
    m11*v1**2/2 + m22*v2**2/2 - m12*v1*v2
    + l1*v1**4/8 + l2*v2**4/8
    + (l3+l4+l5)*v1**2*v2**2/4
    + l6*v1**3*v2/2 + l7*v1*v2**3/2
)
T1 = diff(V0, v1)
T2 = diff(V0, v2)
solve([T1,T2], [m11,m22])
```

It reproduces the hand-derived expressions. The CAS did not construct the potential or choose conventions.

## Convention map

### DH05/BFLRS11 form

The independently derived equations are the real `xi=0` specialization of DH05 Eqs. `minconditionsa`, `minconditionsb`:

\[
\begin{aligned}
m_{11}^2={}&m_{12}^2\tan\beta
-\frac{v^2}{2}
\left[\lambda_1c^2+\lambda_{345}s^2+3\lambda_6sc
+\lambda_7s^2\tan\beta\right],\\
m_{22}^2={}&m_{12}^2\cot\beta
-\frac{v^2}{2}
\left[\lambda_2s^2+\lambda_{345}c^2
+\lambda_6c^2\cot\beta+3\lambda_7sc\right].
\end{aligned}
\]

Using `s^2 tan(beta)=s^3/c` and `c^2 cot(beta)=c^3/s` gives exact equality.

### GHOO18 quadratic-symbol translation

From the term-by-term potential map,

\[
m_{11,\rm DH}^2=-\frac12m_{11,G}^2,\quad
m_{22,\rm DH}^2=-\frac12m_{22,G}^2,\quad
m_{12,\rm DH}^2=\frac12m_{12,G}^2.
\]

Therefore the same stationarity equations in GHOO18 quadratic symbols are

\[
\boxed{
\begin{aligned}
m_{11,G}^2={}&-m_{12,G}^2\tan\beta
+v^2\left[\lambda_1c^2+\lambda_{345}s^2
+3\lambda_6sc+\lambda_7\frac{s^3}{c}\right],\\
m_{22,G}^2={}&-m_{12,G}^2\cot\beta
+v^2\left[\lambda_2s^2+\lambda_{345}c^2
+\lambda_6\frac{c^3}{s}+3\lambda_7sc\right].
\end{aligned}}
\]

This is `[DERIVED][TRANSLATED]`, not a separately quoted GHOO18 result.

| Symbol | Definition | Basis | Dimension | Status |
|---|---|---|---:|---|
| `v1,v2` | neutral VEVs | selected generic/Type-I basis | 1 | vacuum coordinates |
| `v` | `sqrt(v1^2+v2^2)` | VEV norm | 1 | derived |
| `tan beta` | `v2/v1` | basis-dependent, project basis fixed | 0 | derived coordinate |
| `m11^2` | coefficient of `Phi1†Phi1` | generic | 2 | eliminable by stationarity |
| `m22^2` | coefficient of `Phi2†Phi2` | generic | 2 | eliminable by stationarity |
| `m12^2` | coefficient in `-[m12^2 Phi1†Phi2+h.c.]` | generic | 2 | off-diagonal mass parameter |
| `M^2` | `m12^2/(s_beta c_beta)` | derived | 2 | derived coordinate |
| `lambda1...lambda7` | quartic coefficients | generic | 0 | potential parameters |
| `lambda345` | `lambda3+lambda4+lambda5` | real CP basis shorthand | 0 | derived shorthand |

## What was checked against the source

1. DH05 `hbasis.tex` source lines 2236–2257. Setting `xi=0` and all potentially complex coefficients real gives exactly the two independently derived equations.
2. BFLRS11 explicitly states its `notation 1` follows Davidson–Haber.
3. Active 2HDMC `THDM::set_param_gen` computes

```cpp
m22_2 = m12_2*ctb
       - 0.5*v2*(lambda[2]*sb2
       + (lambda[3]+lambda[4]+lambda[5])*cb2
       + lambda[6]*cb2*ctb
       + 3.*lambda[7]*sb*cb);
```

where the C++ variable `v2` denotes the electroweak `v^2`. This is exactly the derived `m22^2` condition and independently checks the active implementation signs and `lambda6/lambda7` factors.
4. Dimensional check: each tadpole has mass dimension three; after division by a VEV, each solved `m_ii^2` expression has dimension two.
5. Limiting check: `lambda6=lambda7=0` reduces to the standard softly broken `Z2` stationarity equations.

## What remains uncertain

Stationarity establishes an extremum, not that it is the desired local/global minimum. Positivity of the Hessian in physical directions belongs to the mass-matrix/stability stages. No physical-state identity has yet been assumed. The Higgs-basis relation `Y3=-Z6 v^2/2` has not been re-derived here; it remains reserved for the explicit Higgs-basis rotation.

## Next smallest validation

Expand the potential to quadratic order around the stationary point and construct separately the charged, CP-odd, and CP-even Hessians. The first check must be that the charged and CP-odd matrices possess the VEV-direction zero mode required by electroweak symmetry breaking.
