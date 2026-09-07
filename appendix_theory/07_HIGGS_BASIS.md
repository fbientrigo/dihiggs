# Phase 7 — Explicit Higgs-basis construction and alignment criterion

This phase starts from the generic-basis fields and potential already fixed in Phases 1–3 and performs the Higgs-flavor rotation explicitly. No Higgs-basis coefficient is imported as an input.

The frozen real CP-conserving rotation is

\[
\boxed{
H_1=c_\beta\Phi_1+s_\beta\Phi_2,
\qquad
H_2=-s_\beta\Phi_1+c_\beta\Phi_2
}
\]

with inverse

\[
\boxed{
\Phi_1=c_\beta H_1-s_\beta H_2,
\qquad
\Phi_2=s_\beta H_1+c_\beta H_2
}.
\]

The sign of `H2` is now a project convention. It is not left available for later adjustment of trilinear signs.

---

## What is established

[DERIVED] The vacuum is carried entirely by `H1`:

\[
\boxed{
\langle H_1\rangle=\frac1{\sqrt2}\binom{0}{v},
\qquad
\langle H_2\rangle=0
}.
\]

With the Phase-4 physical fields,

\[
\boxed{
H_1=\begin{pmatrix}
G^+\\[2mm]
\dfrac{v+\rho_v+iG^0}{\sqrt2}
\end{pmatrix},
\qquad
H_2=\begin{pmatrix}
H^+\\[2mm]
\dfrac{\rho_\perp+iA}{\sqrt2}
\end{pmatrix}}
\]

where

\[
\rho_v=c_\beta\rho_1+s_\beta\rho_2,
\qquad
\rho_\perp=-s_\beta\rho_1+c_\beta\rho_2.
\]

On the exact-alignment branch already frozen in Phase 4,

\[
h=\rho_v,
\qquad
\phi=-\rho_\perp,
\]

so

\[
\boxed{H_2^0=\frac{-\phi+iA}{\sqrt2}}
\]

in the project convention. This sign is a downstream input to every odd-`H2` scalar coupling.

[DERIVED] We write the Higgs-basis potential in the GHOO/invariant `Y_i,Z_i` convention

\[
\begin{aligned}
V={}&Y_1 H_1^\dagger H_1+Y_2 H_2^\dagger H_2
+\left[Y_3H_1^\dagger H_2+\mathrm{h.c.}\right]\\
&+\frac12Z_1(H_1^\dagger H_1)^2
+\frac12Z_2(H_2^\dagger H_2)^2
+Z_3(H_1^\dagger H_1)(H_2^\dagger H_2)
+Z_4(H_1^\dagger H_2)(H_2^\dagger H_1)\\
&+\left\{\frac12Z_5(H_1^\dagger H_2)^2
+\left[Z_6(H_1^\dagger H_1)+Z_7(H_2^\dagger H_2)\right]H_1^\dagger H_2
+\mathrm{h.c.}\right\}.
\end{aligned}
\]

All coefficients are real in the CP-conserving convention used here.

[DERIVED] The quadratic map is

\[
\boxed{Y_1=m_{11}^2c^2+m_{22}^2s^2-2m_{12}^2sc}
\]

\[
\boxed{Y_2=m_{11}^2s^2+m_{22}^2c^2+2m_{12}^2sc}
\]

\[
\boxed{Y_3=(m_{22}^2-m_{11}^2)sc-m_{12}^2(c^2-s^2)}.
\]

[DERIVED] Defining

\[
\lambda_{345}=\lambda_3+\lambda_4+\lambda_5,
\quad
s_{2\beta}=2sc,
\quad
c_{2\beta}=c^2-s^2,
\]

the quartic coefficients are

\[
\boxed{
Z_1=\lambda_1c^4+\lambda_2s^4
+\frac12\lambda_{345}s_{2\beta}^2
+2s_{2\beta}(c^2\lambda_6+s^2\lambda_7)}
\]

\[
\boxed{
Z_2=\lambda_1s^4+\lambda_2c^4
+\frac12\lambda_{345}s_{2\beta}^2
-2s_{2\beta}(s^2\lambda_6+c^2\lambda_7)}
\]

and, for `i=3,4,5`,

\[
\boxed{
Z_i=\lambda_i
+\frac14s_{2\beta}^2(\lambda_1+\lambda_2-2\lambda_{345})
-s_{2\beta}c_{2\beta}(\lambda_6-\lambda_7)}.
\]

The coefficients that carry the residual `H2` sign are

\[
\boxed{
\begin{aligned}
Z_6={}&-\frac12s_{2\beta}
\left(\lambda_1c^2-\lambda_2s^2-\lambda_{345}c_{2\beta}\right)
+c\cos3\beta\,\lambda_6+s\sin3\beta\,\lambda_7,
\end{aligned}}
\]

\[
\boxed{
\begin{aligned}
Z_7={}&-\frac12s_{2\beta}
\left(\lambda_1s^2-\lambda_2c^2+\lambda_{345}c_{2\beta}\right)
+s\sin3\beta\,\lambda_6+c\cos3\beta\,\lambda_7.
\end{aligned}}
\]

Equivalent polynomial forms, useful for a line-by-line algebra audit, are

\[
\boxed{
\begin{aligned}
Z_6={}&c^4\lambda_6-c^3s(\lambda_1-\lambda_{345})
-3c^2s^2\lambda_6+3c^2s^2\lambda_7\\
&+cs^3(\lambda_2-\lambda_{345})-s^4\lambda_7,
\end{aligned}}
\]

\[
\boxed{
\begin{aligned}
Z_7={}&c^4\lambda_7+c^3s(\lambda_2-\lambda_{345})
+3c^2s^2(\lambda_6-\lambda_7)\\
&+cs^3(\lambda_{345}-\lambda_1)-s^4\lambda_6.
\end{aligned}}
\]

[DERIVED] Higgs-basis stationarity is

\[
\boxed{Y_1=-\frac12Z_1v^2},
\qquad
\boxed{Y_3=-\frac12Z_6v^2}.
\]

[DERIVED] The physical charged and CP-odd masses are

\[
\boxed{m_{H^\pm}^2=Y_2+\frac12Z_3v^2}
\]

and

\[
\boxed{m_A^2=Y_2+\frac12(Z_3+Z_4-Z_5)v^2}.
\]

The CP-even Hessian in the ordered basis `(rho_v,rho_perp)` is

\[
\boxed{
\mathcal M_{\rm even,H}^2=
\begin{pmatrix}
Z_1v^2 & Z_6v^2\\
Z_6v^2 & Y_2+\dfrac12(Z_3+Z_4+Z_5)v^2
\end{pmatrix}}
\]

or equivalently

\[
\boxed{
\mathcal M_{\rm even,H}^2=
\begin{pmatrix}
Z_1v^2 & Z_6v^2\\
Z_6v^2 & m_A^2+Z_5v^2
\end{pmatrix}}.
\]

Therefore, with alignment defined physically as the VEV direction `rho_v` being a CP-even mass eigenvector,

\[
\boxed{\text{exact alignment}\iff Z_6=0}.
\]

For the project mass-state convention

\[
h=s_{\beta-\alpha}\rho_v+c_{\beta-\alpha}\rho_\perp,
\qquad
\phi=c_{\beta-\alpha}\rho_v-s_{\beta-\alpha}\rho_\perp,
\]

the same matrix gives the exact identities

\[
\boxed{Z_1v^2=m_h^2s_{\beta-\alpha}^2+m_\phi^2c_{\beta-\alpha}^2}
\]

\[
\boxed{Z_6v^2=(m_h^2-m_\phi^2)s_{\beta-\alpha}c_{\beta-\alpha}}.
\]

[DERIVED] Since the exact generic-to-Higgs map is now known, the large-`tan beta`, `lambda7=0` limit is also fixed:

\[
\boxed{Z_7\longrightarrow-\lambda_6\qquad (\tan\beta\to\infty,\ \lambda_7=0)}.
\]

This is an algebraic statement about the frozen `H2` sign. It is not yet a statement about the sign of the physical `phi H+H-` Feynman rule.

---

## Derivation

### 1. Why this rotation is the Higgs basis

Start from

\[
\langle\Phi_1\rangle=\frac1{\sqrt2}\binom0{v_1},
\qquad
\langle\Phi_2\rangle=\frac1{\sqrt2}\binom0{v_2},
\]

with

\[
v_1=vc,
\qquad
v_2=vs.
\]

Then

\[
\langle H_1\rangle
=c\langle\Phi_1\rangle+s\langle\Phi_2\rangle
=\frac1{\sqrt2}\binom0{v(c^2+s^2)}
=\frac1{\sqrt2}\binom0v,
\]

while

\[
\langle H_2\rangle
=-s\langle\Phi_1\rangle+c\langle\Phi_2\rangle
=\frac1{\sqrt2}\binom0{-svc+cvs}=0.
\]

Thus the rotation is selected by the vacuum vector itself.

### 2. Component fields and the frozen `H2` sign

Insert

\[
\Phi_i=\begin{pmatrix}\phi_i^+\\(v_i+\rho_i+i\eta_i)/\sqrt2\end{pmatrix}.
\]

The upper components give

\[
H_1^+=c\phi_1^++s\phi_2^+=G^+,
\]

\[
H_2^+=-s\phi_1^++c\phi_2^+=H^+.
\]

The imaginary neutral parts give

\[
\operatorname{Im}H_1^0=\frac{c\eta_1+s\eta_2}{\sqrt2}=\frac{G^0}{\sqrt2},
\]

\[
\operatorname{Im}H_2^0=\frac{-s\eta_1+c\eta_2}{\sqrt2}=\frac{A}{\sqrt2}.
\]

The real parts are exactly `rho_v` and `rho_perp`. No additional sign is introduced.

### 3. Bilinears needed to transform the potential

Define

\[
A_H\equiv H_1^\dagger H_1,
\quad
B_H\equiv H_2^\dagger H_2,
\quad
C_H\equiv H_1^\dagger H_2,
\quad
D_H\equiv H_2^\dagger H_1.
\]

Using the inverse rotation,

\[
\Phi_1=cH_1-sH_2,
\qquad
\Phi_2=sH_1+cH_2,
\]

one obtains

\[
\boxed{
\Phi_1^\dagger\Phi_1
=c^2A_H+s^2B_H-cs(C_H+D_H)}
\]

\[
\boxed{
\Phi_2^\dagger\Phi_2
=s^2A_H+c^2B_H+cs(C_H+D_H)}
\]

\[
\boxed{
\Phi_1^\dagger\Phi_2
=cs(A_H-B_H)+c^2C_H-s^2D_H}
\]

\[
\boxed{
\Phi_2^\dagger\Phi_1
=cs(A_H-B_H)-s^2C_H+c^2D_H}.
\]

These four identities are sufficient to reconstruct every coefficient in the transformed potential.

### 4. Quadratic coefficients operator by operator

The generic quadratic potential is

\[
V_2=m_{11}^2\Phi_1^\dagger\Phi_1
+m_{22}^2\Phi_2^\dagger\Phi_2
-m_{12}^2(\Phi_1^\dagger\Phi_2+\Phi_2^\dagger\Phi_1).
\]

Collecting the coefficient of `A_H` gives

\[
Y_1=m_{11}^2c^2+m_{22}^2s^2-2m_{12}^2sc.
\]

Collecting `B_H` gives

\[
Y_2=m_{11}^2s^2+m_{22}^2c^2+2m_{12}^2sc.
\]

Collecting `C_H` (and identically `D_H` because parameters are real) gives

\[
Y_3=(m_{22}^2-m_{11}^2)sc-m_{12}^2(c^2-s^2).
\]

This last sign is useful for translating DH05: DH writes the Higgs-basis mixed quadratic operator as `-[M12^2 H1^dag H2+h.c.]`, hence

\[
\boxed{Y_3=-M_{12,\rm DH}^2}
\]

for the real `chi=xi=0` convention.

### 5. Quartic coefficients: extraction rule

Substitute the four bilinears above into

\[
\begin{aligned}
V_4={}&\frac12\lambda_1(\Phi_1^\dagger\Phi_1)^2
+\frac12\lambda_2(\Phi_2^\dagger\Phi_2)^2
+\lambda_3(\Phi_1^\dagger\Phi_1)(\Phi_2^\dagger\Phi_2)\\
&+\lambda_4(\Phi_1^\dagger\Phi_2)(\Phi_2^\dagger\Phi_1)
+\frac12\lambda_5\left[(\Phi_1^\dagger\Phi_2)^2+(\Phi_2^\dagger\Phi_1)^2\right]\\
&+\lambda_6(\Phi_1^\dagger\Phi_1)(\Phi_1^\dagger\Phi_2+\Phi_2^\dagger\Phi_1)
+\lambda_7(\Phi_2^\dagger\Phi_2)(\Phi_1^\dagger\Phi_2+\Phi_2^\dagger\Phi_1).
\end{aligned}
\]

Then identify

\[
\frac12Z_1A_H^2,
\quad
\frac12Z_2B_H^2,
\quad
Z_3A_HB_H,
\quad
Z_4C_HD_H,
\]

\[
\frac12Z_5(C_H^2+D_H^2),
\quad
Z_6A_H(C_H+D_H),
\quad
Z_7B_H(C_H+D_H).
\]

For example, the coefficient of `A_H C_H` is `Z6`, not `2 Z6`, because the hermitian conjugate multiplies `A_H D_H` rather than duplicating `A_H C_H`. This is one of the normalization points checked symbolically.

The full collection yields the boxed map in the previous section.

### 6. Higgs-basis tadpoles derived directly

For the neutral fields write

\[
H_1^0=\frac{v+\rho_v+iG^0}{\sqrt2},
\qquad
H_2^0=\frac{\rho_\perp+iA}{\sqrt2}.
\]

At linear order in `rho_v`,

\[
H_1^\dagger H_1=\frac{v^2}{2}+v\rho_v+\cdots.
\]

Therefore

\[
\left.\frac{\partial V}{\partial\rho_v}\right|_0
=v\left(Y_1+\frac12Z_1v^2\right).
\]

For `rho_perp`,

\[
H_1^\dagger H_2=\frac{v\rho_\perp}{2}+\cdots,
\]

so the real mixed quadratic term contributes `Y3 v rho_perp`, while

\[
[Z_6(H_1^\dagger H_1)(H_1^\dagger H_2)+\mathrm{h.c.}]
\supset\frac12Z_6v^3\rho_\perp.
\]

Thus

\[
\left.\frac{\partial V}{\partial\rho_\perp}\right|_0
=v\left(Y_3+\frac12Z_6v^2\right).
\]

Setting both tadpoles to zero gives

\[
Y_1=-\frac12Z_1v^2,
\qquad
Y_3=-\frac12Z_6v^2.
\]

As a second independent algebra check, substituting the Phase-3 generic-basis tadpole solutions for `m11^2,m22^2` into the derived formulas for `Y1,Y3` reproduces these two identities exactly.

### 7. Charged and CP-odd masses in the Higgs basis

Because `H2` has zero VEV, its charged field obtains at quadratic order

\[
V\supset\left(Y_2+\frac12Z_3v^2\right)H^+H^-.
\]

Hence

\[
m_{H^\pm}^2=Y_2+\frac12Z_3v^2.
\]

For the CP-odd field `A`, direct differentiation gives

\[
m_A^2=Y_2+\frac12(Z_3+Z_4-Z_5)v^2.
\]

Subtracting,

\[
m_{H^\pm}^2-m_A^2=\frac12(Z_5-Z_4)v^2,
\]

which is the Higgs-basis form of the Phase-4 generic result.

### 8. CP-even Hessian and `Z6`

The second derivatives in `(rho_v,rho_perp)` after using the Higgs-basis tadpoles are

\[
\frac{\partial^2V}{\partial\rho_v^2}=Z_1v^2,
\]

\[
\boxed{
\frac{\partial^2V}{\partial\rho_v\partial\rho_\perp}=Z_6v^2},
\]

\[
\frac{\partial^2V}{\partial\rho_\perp^2}
=Y_2+\frac12(Z_3+Z_4+Z_5)v^2.
\]

This is the missing first-principles step behind the alignment criterion.

The vector corresponding to the VEV direction in this basis is simply

\[
\binom10.
\]

It is an eigenvector of the CP-even Hessian if and only if the lower component of

\[
\mathcal M_{\rm even,H}^2\binom10
\]

vanishes, i.e.

\[
Z_6v^2=0.
\]

Since `v != 0`,

\[
\boxed{\rho_v\text{ is a mass eigenstate}\iff Z_6=0}.
\]

Combined with Phase 6, which independently showed that `rho_v` has the SM tree-level `VV` coupling, this proves

\[
\boxed{\text{exact tree-level alignment}\iff Z_6=0}
\]

within the CP-conserving neutral sector and frozen Higgs-basis sign convention.

At an exactly degenerate CP-even spectrum the mass eigenvectors are not unique, but the matrix degeneracy itself still requires the off-diagonal entry to vanish. What is lost at degeneracy is a unique mixing-angle label, not the algebraic condition `Z6=0`.

### 9. Relation to physical masses and `beta-alpha`

From Phase 4,

\[
\binom{h}{\phi}
=
\begin{pmatrix}
s_{\beta-\alpha}&c_{\beta-\alpha}\\
c_{\beta-\alpha}&-s_{\beta-\alpha}
\end{pmatrix}
\binom{\rho_v}{\rho_\perp}.
\]

Therefore

\[
\mathcal M_{\rm even,H}^2
=R^T\begin{pmatrix}m_h^2&0\\0&m_\phi^2\end{pmatrix}R.
\]

The `(1,1)` entry yields

\[
Z_1v^2=m_h^2s_{\beta-\alpha}^2+m_\phi^2c_{\beta-\alpha}^2,
\]

and the off-diagonal entry yields

\[
\boxed{Z_6v^2=(m_h^2-m_\phi^2)s_{\beta-\alpha}c_{\beta-\alpha}}.
\]

### 10. First controlled large-`tan beta` consequence

For `lambda7=0` and `t=tan beta`, the exact expression is

\[
Z_7=
-\frac{\lambda_6t^4+(\lambda_1-\lambda_{345})t^3-3\lambda_6t^2+(\lambda_{345}-\lambda_2)t}{(1+t^2)^2}.
\]

Hence

\[
\boxed{Z_7=-\lambda_6+\frac{\lambda_{345}-\lambda_1}{t}+\mathcal O(t^{-2})}.
\]

In particular,

\[
\boxed{Z_7\to-\lambda_6}.
\]

This validates the basis-translation part of the project hypothesis. The scalar trilinear still requires a separate operator expansion because `phi=-rho_perp` and because potential coefficient, Lagrangian coefficient and Feynman rule have different signs/factors.

---

## Convention map

| Object | Project convention | DH05 | GHOO18 | 2HDMC | Status |
|---|---|---|---|---|---|
| `H1` | `c Phi1+s Phi2` | same at `xi=chi=0` | same real Higgs basis | `get_param_higgs` basis | `[DERIVED][SOURCE]` |
| `H2` | `-s Phi1+c Phi2` | same at `xi=chi=0` | same real Higgs basis | same transformation in `get_param_higgs` | `[DERIVED][SOURCE][IMPLEMENTATION-CHECKED]` |
| mixed quadratic | `+[Y3 H1dag H2+h.c.]` | `-[M12_H^2 H1dag H2+h.c.]` | `+[Y3 ...]` | internal Higgs-basis quadratic object | `[TRANSLATED]` |
| relation | `Y3=-M12_H^2` | follows from potential signs | uses `Y3` | consistent with stationarity | `[DERIVED][TRANSLATED]` |
| quartics | `Z1...Z7` | `Lambda1...Lambda7` at `chi=0` | `Z1...Z7` | `get_param_higgs` returns names `Lambda1...Lambda7` | `[TRANSLATED][IMPLEMENTATION-CHECKED]` |
| `Z6` | coefficient of `A_H C_H+h.c.` | `Lambda6` | `Z6` | returned `Lambda6` matches analytic project `Z6` | `[DERIVED][SOURCE][IMPLEMENTATION-CHECKED]` |
| `Z7` | coefficient of `B_H C_H+h.c.` | `Lambda7` | `Z7` | returned `Lambda7` matches analytic project `Z7` | `[DERIVED][SOURCE][IMPLEMENTATION-CHECKED]` |
| residual sign | frozen `H2=-sPhi1+cPhi2` | `H2 -> -H2` allowed within Higgs-basis family | convention must be declared | trilinear routine adds an extra local sign layer | `[PROJECT-DEFINITION]` |
| exact alignment | `rho_v` mass eigenstate | `Lambda6=0` in real basis | `Z6=0` | hybrid input uses `Z6=(mh^2-mH^2)sba cba/v^2` | `[DERIVED][SOURCE][IMPLEMENTATION-CHECKED]` |

Under a residual field redefinition

\[
H_2\to-H_2,
\]

one has

\[
\boxed{Y_3\to-Y_3,
\quad Z_6\to-Z_6,
\quad Z_7\to-Z_7,
\quad Z_5\to Z_5}.
\]

This is why absolute signs of odd-`H2` couplings are meaningless until the `H2` convention is frozen. It is now frozen for the project.

---

## What was checked against the source

The following comparisons were made only after the analytic map above was obtained.

1. **DH05 Higgs-basis fields.** Eq. `higgsbasis` at `xi=chi=0` gives exactly the project `H1,H2`; Eq. `abbasis` gives the same component fields.
2. **DH05 quadratic map.** Eqs. `maa`–`mab` match `Y1,Y2` and `Y3=-M12_H^2` after accounting for DH's minus sign in the mixed quadratic operator.
3. **DH05 quartic map.** Eqs. `Lam1def`–`Lam7def` specialized to `xi=chi=0` agree term by term with the derived `Z1...Z7`, with `Z_i=Lambda_i`.
4. **DH05 stationarity.** `M11^2=-Lambda1 v^2/2`, `M12_H^2=+Lambda6 v^2/2`; with `Y3=-M12_H^2`, this is exactly the project stationarity equation.
5. **DH05 CP-even Higgs-basis matrix.** Its matrix matches the direct Hessian.
6. **2HDMC `get_param_higgs`.** The implementation computes `Lambda1...Lambda7` using the same compact trigonometric formulas as the project `Z1...Z7`. Thus returned `Lambda6,Lambda7` are identified with project `Z6,Z7` in the frozen `H2` convention.
7. **2HDMC hybrid input.** `set_param_hybrid_sba` computes `Z6=(mh^2-mH^2)cba*sba/v^2`, matching the physical-mass identity.

### Resolution of the Phase-6 implementation warning

Phase 6 observed that `get_coupling_hhh` performs

```cpp
get_param_higgs(..., l6, l7, ...);
double Z6 = -l6, Z7 = -l7;
```

Phase 7 establishes that `get_param_higgs` itself already returns the project/DH `Z6,Z7` signs. Therefore the minus signs inside `get_coupling_hhh` are **not** part of the generic-to-Higgs-basis transformation. They belong to a later Feynman-rule/convention layer local to that trilinear implementation and must be audited in the trilinear phase.

### Symbolic audit

`appendix_theory/checks/phase7_higgs_basis_check.py` reconstructs the bilinears and all `Y_i,Z_i`, derives the Higgs-basis tadpoles and CP-even Hessian, verifies the DH/2HDMC compact forms, checks the physical `Z6` identity, and verifies `Z7 -> -lambda6` at large `tan beta` with `lambda7=0`. It returns `PHASE7_HIGGS_BASIS_CHECK=PASS`.

---

## What remains uncertain

1. The Higgs-basis construction itself is closed: `Y_i`, `Z_i`, stationarity, charged/CP-odd masses and the CP-even mass matrix are independently derived.
2. `Y3=-Z6v^2/2` is no longer source-only; it is independently verified.
3. Exact alignment `iff Z6=0` is independently verified when alignment is defined as the VEV direction being a CP-even mass eigenstate.
4. The exact generic-basis `Z7` map is independently verified, and its `lambda7=0`, large-`tan beta` limit is `Z7 -> -lambda6`.
5. **Still open:** the sign convention inside 2HDMC `get_coupling_hhh`. Phase 7 shows that its local `Z6=-l6, Z7=-l7` is not the Higgs-basis rotation.
6. The physical `phi H+H-` trilinear has not been derived here. In particular, `Z7 -> -lambda6` alone does not determine its sign because exact alignment gives `phi=-rho_perp`.
7. `X=lambda6 tan beta` remains deliberately absent from the foundational derivation.

---

## Next smallest validation

Starting from the now-verified Higgs-basis potential and

\[
H_2=\begin{pmatrix}H^+\\(\rho_\perp+iA)/\sqrt2\end{pmatrix},
\]

extract only terms linear in `rho_v` or `rho_perp` and bilinear in `H+ H-`. Keep three objects explicitly separate:

1. coefficient in the potential `V`;
2. coefficient in `L_int=-V_int`;
3. Feynman rule including the overall `i`.

Then substitute

\[
rho_v=s_{\beta-\alpha}h+c_{\beta-\alpha}\phi,
\qquad
rho_\perp=c_{\beta-\alpha}h-s_{\beta-\alpha}\phi,
\]

and only afterward take exact alignment. This is the shortest path to a convention-safe `phi H+H-` result.
