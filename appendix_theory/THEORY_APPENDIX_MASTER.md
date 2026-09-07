# 2HDM theory appendix — canonical derivation and convention audit

Issue: `fbientrigo/dihiggs#81`  
Status: **canonical single-document audit record**  
Coverage: **Phases 0–9**

This document is the primary theory record. It is written so that an independent reviewer can reconstruct every sign, normalization and basis map without trusting project code or a copied formula. Phase-specific files and symbolic scripts are auxiliary checks only.

The governing order is

\[
\boxed{
\text{fields}\to\text{potential}\to\text{vacuum}\to\text{stationarity}
\to\text{mass matrices}\to\text{physical states}\to\text{Yukawa/gauge couplings}
\to\text{Higgs basis}\to\text{scalar trilinears}\to\text{loop amplitudes}\to X
}.
\]

No later phenomenological coordinate may be used to choose an earlier sign.

---

# A. Frozen conventions and audit rules

## A.1 Epistemic labels

- `[SOURCE]`: explicitly stated in an inspected source.
- `[DERIVED]`: obtained algebraically from definitions already frozen.
- `[TRANSLATED]`: obtained from an explicit convention map.
- `[IMPLEMENTATION-CHECKED]`: compared with active code only after the analytic result was frozen.
- `[PROJECT-DEFINITION]`: a convention chosen by this project.
- `[OPEN-QUESTION]`: not yet established strongly enough for downstream use.

## A.2 Objects that must not be conflated

\[
m_{22}^2\neq m_{12}^2\neq M^2
\]

generically, and

\[
\lambda_6,\lambda_7\neq Z_6,Z_7.
\]

For scalar trilinears three objects will always be separated:

1. coefficient in `V`;
2. coefficient in `L_int=-V_int`;
3. Feynman rule including the overall `i` and combinatorics.

The project coordinate

\[
X\equiv\lambda_6\tan\beta
\]

is deliberately deferred until the analytic scalar coupling is established.

## A.3 Source anchors

1. Davidson–Haber, hep-ph/0504050 (DH05): main convention anchor.
2. Branco et al., arXiv:1106.0034 (BFLRS11): independent review cross-check and Type-I Yukawa structure.
3. Grzadkowski–Haber–Ogreid–Osland, arXiv:1808.01472 (GHOO18): Higgs-basis/alignment/trilinear translation.
4. Active vendored 2HDMC: implementation check only.

DH05 uses `Q=T3+Y/2` with scalar-doublet `Y=1`; this is equivalent to the modern convention `Q=T3+Y` with `Y=1/2`.

---

# Phase 0 — Source and convention inventory

## What is established

The project generic-basis potential uses the DH05/BFLRS11 sign convention. GHOO18 uses a different quadratic normalization.

DH/DH-like quadratic terms:

\[
V_2=m_{11}^2\Phi_1^\dagger\Phi_1+m_{22}^2\Phi_2^\dagger\Phi_2
-[m_{12}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}].
\]

GHOO18 writes

\[
V_{2,G}=-\frac12\{m_{11,G}^2\Phi_1^\dagger\Phi_1+m_{22,G}^2\Phi_2^\dagger\Phi_2+[m_{12,G}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}]\}.
\]

## Derivation

Operator matching gives

\[
\boxed{m_{11,\rm DH}^2=-\frac12m_{11,G}^2},\qquad
\boxed{m_{22,\rm DH}^2=-\frac12m_{22,G}^2},\qquad
\boxed{m_{12,\rm DH}^2=+\frac12m_{12,G}^2}.
\]

The DH CP-even field convention is frozen as

\[
h_{\rm DH}=-s_\alpha\rho_1+c_\alpha\rho_2,
\qquad
H_{\rm DH}=c_\alpha\rho_1+s_\alpha\rho_2.
\]

## Convention map

The early simple BFLRS11 display uses the global negatives of these two fields. Such a global field flip is harmless for masses but flips every interaction with an odd number of that field.

## What was checked against the source

DH05 generic potential, VEVs, CP-even rotation and Higgs basis; BFLRS11 generic potential and Type-I table; GHOO18 generic/Higgs-basis potentials.

## What remains uncertain

No project `h/phi` map, Yukawa modifier, gauge modifier or scalar trilinear is inferred at Phase 0.

## Next smallest validation

Define the fields and vacuum before introducing physical states.

---

# Phase 1 — Fields, VEVs and beta

## What is established

\[
\boxed{
\Phi_i=\begin{pmatrix}
\phi_i^+\\[1mm]
(v_i+\rho_i+i\eta_i)/\sqrt2
\end{pmatrix}},\qquad i=1,2.
\]

The real neutral vacuum is

\[
\langle\Phi_1\rangle=\frac1{\sqrt2}\binom0{v_1},\qquad
\langle\Phi_2\rangle=\frac1{\sqrt2}\binom0{v_2}.
\]

Define

\[
\boxed{v^2=v_1^2+v_2^2},\qquad
\boxed{\tan\beta=\frac{v_2}{v_1}},
\]

so

\[
v_1=vc_\beta,\qquad v_2=vs_\beta.
\]

## Derivation

The factor `1/sqrt(2)` canonically normalizes `rho_i,eta_i`. With `Y=1/2`, the upper component has electric charge +1 and the lower component charge 0. Choosing VEVs only in neutral components preserves electromagnetism.

## Convention map

Project notation: `rho_i` are CP-even fluctuations; `eta_i` are CP-odd fluctuations. Positive `tan beta` fixes the coordinate patch used by 2HDMC/project scans.

## What was checked against the source

DH05 `potmin`, `tanbdef`; BFLRS11 component expansion.

## What remains uncertain

No physical scalar state has been assumed.

## Next smallest validation

Construct the complete renormalizable scalar potential from gauge singlets.

---

# Phase 2 — Generic CP-conserving scalar potential

## What is established

The project potential is

\[
\boxed{
\begin{aligned}
V={}&m_{11}^2\Phi_1^\dagger\Phi_1+m_{22}^2\Phi_2^\dagger\Phi_2
-[m_{12}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}]\\
&+\frac12\lambda_1(\Phi_1^\dagger\Phi_1)^2
+\frac12\lambda_2(\Phi_2^\dagger\Phi_2)^2\\
&+\lambda_3(\Phi_1^\dagger\Phi_1)(\Phi_2^\dagger\Phi_2)
+\lambda_4(\Phi_1^\dagger\Phi_2)(\Phi_2^\dagger\Phi_1)\\
&+\left\{\frac12\lambda_5(\Phi_1^\dagger\Phi_2)^2
+[\lambda_6(\Phi_1^\dagger\Phi_1)+\lambda_7(\Phi_2^\dagger\Phi_2)]\Phi_1^\dagger\Phi_2
+\mathrm{h.c.}\right\}.
\end{aligned}}
\]

All coefficients are real in the present CP-conserving branch.

## Derivation

The independent bilinears are `Phi_i^dag Phi_j`. Dimension-two Hermitian combinations give the quadratic part. Products of bilinears give all renormalizable quartics. Hermiticity determines the explicit h.c. structure and the `1/2` factors.

## Convention map

The potential is term-by-term DH05/BFLRS11. The campaign condition `lambda7=0` is not part of the model definition.

## What was checked against the source

DH05 Eq. `pot`, BFLRS11 Eq. `2_VH1`. Active 2HDMC is used only later as a posterior check.

## What remains uncertain

The potential does not yet identify a stationary vacuum or masses.

## Next smallest validation

Insert the real neutral VEVs and differentiate.

---

# Phase 3 — Vacuum potential and minimization

Define

\[
\lambda_{345}\equiv\lambda_3+\lambda_4+\lambda_5.
\]

## What is established

At the neutral vacuum,

\[
\Phi_1^\dagger\Phi_1=\frac{v_1^2}{2},\quad
\Phi_2^\dagger\Phi_2=\frac{v_2^2}{2},\quad
\Phi_1^\dagger\Phi_2=\frac{v_1v_2}{2}.
\]

Therefore

\[
\boxed{
\begin{aligned}
V_0={}&\frac12m_{11}^2v_1^2+\frac12m_{22}^2v_2^2-m_{12}^2v_1v_2
+\frac18\lambda_1v_1^4+\frac18\lambda_2v_2^4\\
&+\frac14\lambda_{345}v_1^2v_2^2
+\frac12\lambda_6v_1^3v_2+\frac12\lambda_7v_1v_2^3.
\end{aligned}}
\]

The tadpoles are

\[
\boxed{
0=m_{11}^2v_1-m_{12}^2v_2+\frac12\lambda_1v_1^3
+\frac12\lambda_{345}v_1v_2^2+\frac32\lambda_6v_1^2v_2+\frac12\lambda_7v_2^3}
\]

and

\[
\boxed{
0=m_{22}^2v_2-m_{12}^2v_1+\frac12\lambda_2v_2^3
+\frac12\lambda_{345}v_1^2v_2+\frac12\lambda_6v_1^3+\frac32\lambda_7v_1v_2^2}.
\]

## Derivation

The factors `3 lambda6` and `3 lambda7` arise directly from

\[
\partial_{v_1}(v_1^3v_2)=3v_1^2v_2,
\qquad
\partial_{v_2}(v_1v_2^3)=3v_1v_2^2.
\]

Solving for the diagonal quadratic coefficients gives

\[
\boxed{
\begin{aligned}
m_{11}^2={}&m_{12}^2\frac{v_2}{v_1}
-\frac12\left[\lambda_1v_1^2+\lambda_{345}v_2^2+3\lambda_6v_1v_2+\lambda_7\frac{v_2^3}{v_1}\right],\\
m_{22}^2={}&m_{12}^2\frac{v_1}{v_2}
-\frac12\left[\lambda_2v_2^2+\lambda_{345}v_1^2+\lambda_6\frac{v_1^3}{v_2}+3\lambda_7v_1v_2\right].
\end{aligned}}
\]

Only now define

\[
\boxed{M^2\equiv\frac{m_{12}^2}{s_\beta c_\beta}}.
\]

Then

\[
\boxed{m_{22}^2\neq m_{12}^2\neq M^2}
\]

generically.

## Convention map

`M2` in project data means this derived `M^2`, not the generic coefficient `m22_2`.

## What was checked against the source

DH05 stationarity agrees after CP-conserving specialization. A symbolic derivative audit reproduces both equations. Active 2HDMC `set_param_gen` reproduces the derived `m22^2` expression term by term.

## What remains uncertain

Stationarity is not proof of a global minimum.

## Next smallest validation

Compute charged, CP-odd and CP-even Hessians directly from the same potential.

---

# Phase 4 — Mass matrices, Goldstones and physical scalar states

## What is established

After using the Phase-3 tadpoles, both charged and CP-odd matrices factorize as

\[
\boxed{\mathcal M^2=D
\begin{pmatrix}
\tan\beta&-1\\-1&\cot\beta
\end{pmatrix}}.
\]

Thus `(v1,v2)` is a null eigenvector and `(-v2,v1)` the physical orthogonal direction.

Hence

\[
\boxed{G^+=c\phi_1^++s\phi_2^+},\qquad
\boxed{H^+=-s\phi_1^++c\phi_2^+},
\]

\[
\boxed{G^0=c\eta_1+s\eta_2},\qquad
\boxed{A=-s\eta_1+c\eta_2}.
\]

The masses are

\[
\boxed{m_{H^\pm}^2=M^2-\frac{v^2}{2}(\lambda_4+\lambda_5+\lambda_6\cot\beta+\lambda_7\tan\beta)}
\]

and

\[
\boxed{m_A^2=M^2-\frac{v^2}{2}(2\lambda_5+\lambda_6\cot\beta+\lambda_7\tan\beta)}.
\]

Therefore

\[
\boxed{m_{H^\pm}^2-m_A^2=\frac{v^2}{2}(\lambda_5-\lambda_4)}.
\]

## Derivation

The CP-even Hessian before tadpole elimination is

\[
\mathcal M_\rho^2=
\begin{pmatrix}M_{11}^2&M_{12}^2\\M_{12}^2&M_{22}^2\end{pmatrix}
\]

with

\[
\begin{aligned}
M_{11}^2&=m_{11}^2+\frac32\lambda_1v_1^2+\frac12\lambda_{345}v_2^2+3\lambda_6v_1v_2,\\
M_{22}^2&=m_{22}^2+\frac32\lambda_2v_2^2+\frac12\lambda_{345}v_1^2+3\lambda_7v_1v_2,\\
M_{12}^2&=-m_{12}^2+\lambda_{345}v_1v_2+\frac32\lambda_6v_1^2+\frac32\lambda_7v_2^2.
\end{aligned}
\]

An equivalent post-tadpole representation useful for 2HDMC is

\[
\boxed{
\begin{aligned}
M_{11}^2&=m_A^2s^2+v^2(\lambda_1c^2+\lambda_5s^2+2\lambda_6sc),\\
M_{22}^2&=m_A^2c^2+v^2(\lambda_2s^2+\lambda_5c^2+2\lambda_7sc),\\
M_{12}^2&=-m_A^2sc+v^2[(\lambda_3+\lambda_4)sc+\lambda_6c^2+\lambda_7s^2].
\end{aligned}}
\]

The DH CP-even rotation is

\[
\boxed{
\binom{h}{\phi}=
\begin{pmatrix}-s_\alpha&c_\alpha\\c_\alpha&s_\alpha\end{pmatrix}
\binom{\rho_1}{\rho_2}}
\]

with project state names `h=h_DH`, `phi=H_DH` on the branch connected to `s_{beta-alpha}=1`.

Define

\[
\rho_v=c\rho_1+s\rho_2,
\qquad
\rho_\perp=-s\rho_1+c\rho_2.
\]

Then

\[
\boxed{h=s_{\beta-\alpha}\rho_v+c_{\beta-\alpha}\rho_\perp},
\]

\[
\boxed{\phi=c_{\beta-\alpha}\rho_v-s_{\beta-\alpha}\rho_\perp}.
\]

At exact alignment,

\[
\boxed{h=\rho_v},\qquad
\boxed{\phi=-\rho_\perp=s\rho_1-c\rho_2}.
\]

## Convention map

The sign of `phi` is inherited from the DH field convention. It is not chosen to force a desired Yukawa or scalar-coupling sign.

## What was checked against the source

DH05 Goldstones/charged state/CP-even states; later general BFLRS11 Hessian; GHOO beta rotations; active 2HDMC masses and CP-even matrix. A symbolic Hessian audit checks all sectors.

## What remains uncertain

Positive masses imply a local quadratic minimum in those directions, not global-vacuum uniqueness.

## Next smallest validation

Derive Type-I fermion couplings from the Yukawa Lagrangian with the frozen scalar signs.

---

# Phase 5 — Type-I Yukawa sector

## What is established

In Type I, all charged fermions couple to `Phi2`:

\[
-\mathcal L_Y=
\bar Q_LY_d\Phi_2d_R+
\bar Q_LY_u\widetilde\Phi_2u_R+
\bar L_LY_\ell\Phi_2\ell_R+\mathrm{h.c.}
\]

After fermion mass diagonalization,

\[
\boxed{m_f=\frac{y_fv_2}{\sqrt2}=\frac{y_fvs_\beta}{\sqrt2}}.
\]

## Derivation

The CP-even neutral interaction is

\[
\mathcal L_Y^{\rm even}=-\sum_f\frac{m_f}{v_2}\rho_2\bar f f.
\]

Because the Phase-4 rotation is symmetric and orthogonal,

\[
\boxed{\rho_2=c_\alpha h+s_\alpha\phi}.
\]

Thus

\[
\boxed{
\mathcal L_Y^{\rm even}=-\sum_f\frac{m_f}{v}
\left(\frac{c_\alpha}{s_\beta}h+\frac{s_\alpha}{s_\beta}\phi\right)\bar f f}
\]

and

\[
\boxed{\kappa_f^h=\frac{c_\alpha}{s_\beta}},\qquad
\boxed{\kappa_f^\phi=\frac{s_\alpha}{s_\beta}}.
\]

At `alpha=beta-pi/2`,

\[
\boxed{\kappa_f^h=1},\qquad
\boxed{\kappa_f^\phi=-\cot\beta}.
\]

## Convention map

The minus sign follows from the already-frozen `phi=s rho1-c rho2`, not from a coupling table.

## What was checked against the source

BFLRS11 confirms Type-I uses `Phi2` for all charged fermions and gives the same functional dependence after sign translation. 2HDMC Type-I uses `set_yukawas_type(1)` and scales the orthogonal-basis Yukawa matrices by `cot beta`.

## What remains uncertain

Nothing blocks the neutral Type-I modifier result.

## Next smallest validation

Derive `hVV/phiVV` directly from scalar kinetic terms.

---

# Phase 6 — Gauge couplings and physical alignment

## What is established

Start from

\[
\boxed{\mathcal L_{\rm kin}=\sum_i(D_\mu\Phi_i)^\dagger(D^\mu\Phi_i)}.
\]

The gauge masses are

\[
\boxed{m_W^2=\frac{g^2v^2}{4}},\qquad
\boxed{m_Z^2=\frac{(g^2+g'^2)v^2}{4}}.
\]

## Derivation

For the charged gauge part,

\[
\mathcal L\supset\frac{g^2}{4}\sum_i(v_i+\rho_i)^2W^+W^-.
\]

The term linear in CP-even fields is

\[
\frac{g^2}{2}(v_1\rho_1+v_2\rho_2)W^+W^-.
\]

For the neutral gauge part, the photon cancels on the neutral VEV because `Q=0`; the `Z` term has the same scalar combination. Therefore

\[
\boxed{v_1\rho_1+v_2\rho_2=v\rho_v}
\]

is the unique tree-level neutral CP-even direction coupled linearly to `WW/ZZ`.

Since

\[
\rho_v=s_{\beta-\alpha}h+c_{\beta-\alpha}\phi,
\]

we obtain

\[
\boxed{\kappa_V^h=s_{\beta-\alpha}},\qquad
\boxed{\kappa_V^\phi=c_{\beta-\alpha}}.
\]

At exact alignment,

\[
\boxed{\kappa_V^h=1},\qquad
\boxed{\kappa_V^\phi=0}.
\]

Also

\[
|v_i+\rho_i+i\eta_i|^2=(v_i+\rho_i)^2+\eta_i^2
\]

contains no term linear in `eta_i`, hence tree-level `AWW=AZZ=0`.

## Convention map

Alignment now has a physical meaning: the CP-even mass eigenstate parallel to the VEV has exactly the SM tree-level `VV` coupling.

## What was checked against the source

DH05/BFLRS11/GHOO gauge-coupling statements; active 2HDMC `q_{k1}=(sba,cba,0,i)` and `get_coupling_vvh` reproduce the derived projection.

## What remains uncertain

Tree-level `phiVV=0` does not imply loop-induced `phi->gamma gamma` or `phi->Z gamma` vanish.

## Next smallest validation

Rotate the entire potential into the Higgs basis and identify `Y_i,Z_i` without importing their formulas.

---

# Phase 7 — Explicit Higgs-basis construction

The project freezes

\[
\boxed{H_1=c_\beta\Phi_1+s_\beta\Phi_2},\qquad
\boxed{H_2=-s_\beta\Phi_1+c_\beta\Phi_2}.
\]

The inverse is

\[
\Phi_1=cH_1-sH_2,\qquad
\Phi_2=sH_1+cH_2.
\]

## What is established

The vacuum rotates to

\[
\boxed{\langle H_1\rangle=\frac1{\sqrt2}\binom0v},\qquad
\boxed{\langle H_2\rangle=0}.
\]

The components are

\[
\boxed{
H_1=\begin{pmatrix}G^+\\(v+\rho_v+iG^0)/\sqrt2\end{pmatrix}},
\qquad
\boxed{
H_2=\begin{pmatrix}H^+\\(\rho_\perp+iA)/\sqrt2\end{pmatrix}}.
\]

Because exact alignment has `phi=-rho_perp`,

\[
\boxed{H_2^0=\frac{-\phi+iA}{\sqrt2}}
\]

in the project convention. This sign is now frozen and cannot be altered later to repair a trilinear sign.

Write the Higgs-basis potential as

\[
\boxed{
\begin{aligned}
V={}&Y_1H_1^\dagger H_1+Y_2H_2^\dagger H_2+[Y_3H_1^\dagger H_2+\mathrm{h.c.}]\\
&+\frac12Z_1(H_1^\dagger H_1)^2+\frac12Z_2(H_2^\dagger H_2)^2
+Z_3(H_1^\dagger H_1)(H_2^\dagger H_2)+Z_4(H_1^\dagger H_2)(H_2^\dagger H_1)\\
&+\left\{\frac12Z_5(H_1^\dagger H_2)^2
+[Z_6(H_1^\dagger H_1)+Z_7(H_2^\dagger H_2)]H_1^\dagger H_2+\mathrm{h.c.}\right\}.
\end{aligned}}
\]

## Derivation

Define Higgs-basis bilinears

\[
A=H_1^\dagger H_1,\quad B=H_2^\dagger H_2,\quad
C=H_1^\dagger H_2,\quad D=H_2^\dagger H_1.
\]

The inverse rotation gives

\[
\boxed{\Phi_1^\dagger\Phi_1=c^2A+s^2B-cs(C+D)},
\]

\[
\boxed{\Phi_2^\dagger\Phi_2=s^2A+c^2B+cs(C+D)},
\]

\[
\boxed{\Phi_1^\dagger\Phi_2=cs(A-B)+c^2C-s^2D},
\]

\[
\boxed{\Phi_2^\dagger\Phi_1=cs(A-B)-s^2C+c^2D}.
\]

Substituting these four identities into the generic potential and collecting independent Higgs-basis operators yields the quadratic map

\[
\boxed{Y_1=m_{11}^2c^2+m_{22}^2s^2-2m_{12}^2sc},
\]

\[
\boxed{Y_2=m_{11}^2s^2+m_{22}^2c^2+2m_{12}^2sc},
\]

\[
\boxed{Y_3=(m_{22}^2-m_{11}^2)sc-m_{12}^2(c^2-s^2)}.
\]

With `lambda345=lambda3+lambda4+lambda5`, `s2b=2sc`, `c2b=c^2-s^2`,

\[
\boxed{Z_1=\lambda_1c^4+\lambda_2s^4+\frac12\lambda_{345}s_{2\beta}^2
+2s_{2\beta}(c^2\lambda_6+s^2\lambda_7)},
\]

\[
\boxed{Z_2=\lambda_1s^4+\lambda_2c^4+\frac12\lambda_{345}s_{2\beta}^2
-2s_{2\beta}(s^2\lambda_6+c^2\lambda_7)},
\]

and for `i=3,4,5`,

\[
\boxed{Z_i=\lambda_i+\frac14s_{2\beta}^2(\lambda_1+\lambda_2-2\lambda_{345})
-s_{2\beta}c_{2\beta}(\lambda_6-\lambda_7)}.
\]

The sign-sensitive coefficients are

\[
\boxed{
Z_6=-\frac12s_{2\beta}(\lambda_1c^2-\lambda_2s^2-\lambda_{345}c_{2\beta})
+c\cos3\beta\,\lambda_6+s\sin3\beta\,\lambda_7},
\]

\[
\boxed{
Z_7=-\frac12s_{2\beta}(\lambda_1s^2-\lambda_2c^2+\lambda_{345}c_{2\beta})
+s\sin3\beta\,\lambda_6+c\cos3\beta\,\lambda_7}.
\]

A useful polynomial form is

\[
\boxed{
\begin{aligned}
Z_7={}&c^4\lambda_7+c^3s(\lambda_2-\lambda_{345})
+3c^2s^2(\lambda_6-\lambda_7)\\
&+cs^3(\lambda_{345}-\lambda_1)-s^4\lambda_6.
\end{aligned}}
\]

### Higgs-basis minimization

Set

\[
H_1^0=\frac{v+\rho_v+iG^0}{\sqrt2},\qquad
H_2^0=\frac{\rho_\perp+iA}{\sqrt2}.
\]

The linear terms give

\[
\left.\frac{\partial V}{\partial\rho_v}\right|_0
=v\left(Y_1+\frac12Z_1v^2\right),
\]

\[
\left.\frac{\partial V}{\partial\rho_\perp}\right|_0
=v\left(Y_3+\frac12Z_6v^2\right).
\]

Therefore

\[
\boxed{Y_1=-\frac12Z_1v^2},\qquad
\boxed{Y_3=-\frac12Z_6v^2}.
\]

This is independently reproduced by substituting the Phase-3 generic tadpoles into the derived `Y1,Y3,Z1,Z6` formulas.

### Masses and alignment in Higgs basis

The charged and CP-odd masses are

\[
\boxed{m_{H^\pm}^2=Y_2+\frac12Z_3v^2},
\]

\[
\boxed{m_A^2=Y_2+\frac12(Z_3+Z_4-Z_5)v^2}.
\]

The CP-even Hessian in `(rho_v,rho_perp)` is

\[
\boxed{
\mathcal M_{H,\rm even}^2=
\begin{pmatrix}
Z_1v^2&Z_6v^2\\
Z_6v^2&m_A^2+Z_5v^2
\end{pmatrix}}.
\]

The VEV direction is `(1,0)`. It is an eigenvector iff the off-diagonal entry vanishes. Since `v != 0`,

\[
\boxed{\text{exact alignment}\iff Z_6=0}.
\]

Using the Phase-4 physical rotation,

\[
\boxed{Z_6v^2=(m_h^2-m_\phi^2)s_{\beta-\alpha}c_{\beta-\alpha}},
\]

and

\[
\boxed{Z_1v^2=m_h^2s_{\beta-\alpha}^2+m_\phi^2c_{\beta-\alpha}^2}.
\]

### Controlled large-`tan beta` limit

For `lambda7=0` and `t=tan beta`, the exact `Z7` becomes

\[
Z_7=-\frac{\lambda_6t^4+(\lambda_1-\lambda_{345})t^3-3\lambda_6t^2+(\lambda_{345}-\lambda_2)t}{(1+t^2)^2}.
\]

Hence

\[
\boxed{Z_7=-\lambda_6+\frac{\lambda_{345}-\lambda_1}{t}+\mathcal O(t^{-2})},
\]

so

\[
\boxed{Z_7\to-\lambda_6}.
\]

This validates the basis-translation hypothesis only. It does **not** yet determine the physical `phi H+H-` sign because `phi=-rho_perp` in exact alignment.

## Convention map

DH05 Higgs-basis potential uses

\[
-[M_{12,H}^2H_1^\dagger H_2+\mathrm{h.c.}],
\]

while the project/GHOO `Y_i` convention uses `+[Y3 H1dag H2+h.c.]`. Thus

\[
\boxed{Y_3=-M_{12,H}^2}.
\]

At real `chi=0`, DH `Lambda_i` equal project `Z_i`. Under the residual redefinition `H2->-H2`,

\[
Y_3,Z_6,Z_7\to-(Y_3,Z_6,Z_7),\qquad Z_5\to Z_5.
\]

The project has frozen `H2=-s Phi1+c Phi2`, so this residual sign is no longer available for later adjustment.

## What was checked against the source

- DH05 `higgsbasis` and `abbasis` match the derived field rotation/components.
- DH05 `maa`–`mab` match `Y1,Y2,Y3=-M12_H^2`.
- DH05 `Lam1def`–`Lam7def` match the derived `Z1...Z7` at `xi=chi=0`.
- DH05 Higgs-basis stationarity matches `Y1=-Z1v^2/2`, `Y3=-Z6v^2/2` after the explicit quadratic-sign translation.
- DH05 CP-even Higgs-basis matrix matches the direct Hessian.
- 2HDMC `get_param_higgs` reproduces the same `Z1...Z7`, including `Z6,Z7` signs.
- 2HDMC `set_param_hybrid_sba` uses `Z6=(mh^2-mH^2)sba*cba/v^2`, matching the derived physical identity.
- `phase7_higgs_basis_check.py` reconstructs the basis map symbolically and returns `PHASE7_HIGGS_BASIS_CHECK=PASS`.

## What remains uncertain

The Higgs-basis map, stationarity, `Z6` alignment criterion, exact `Z7` map and `Z7->-lambda6` limit are closed.

One implementation issue remains isolated: `2HDMC::get_coupling_hhh` calls `get_param_higgs` and then defines local `Z6=-l6`, `Z7=-l7`. Phase 7 proves those extra minuses are **not** part of the basis transformation. Their meaning belongs to the trilinear/Feynman-rule convention layer and must be audited there.

The physical `phi H+H-` coupling remains **NOT VERIFIED**.

## Next smallest validation

Extract from the verified Higgs-basis potential only the terms

\[
\rho_vH^+H^-,\qquad \rho_\perp H^+H^-.
\]

Then convert to `h,phi`, keeping

\[
\rho_v=s_{\beta-\alpha}h+c_{\beta-\alpha}\phi,
\qquad
\rho_\perp=c_{\beta-\alpha}h-s_{\beta-\alpha}\phi,
\]

and only afterward take alignment. Keep potential coefficient, interaction-Lagrangian coefficient and Feynman rule separate.

---

# B. Validation status after Phase 9

| Claim | Status |
|---|---|
| Generic potential and signs | **VERIFIED** |
| Vacuum potential and tadpoles | **VERIFIED** |
| `m22^2`, `m12^2`, `M^2` distinction | **VERIFIED** |
| Scalar Hessians, masses and state signs | **VERIFIED** |
| Type-I `kappa_f^phi=-cot beta` at alignment | **VERIFIED** |
| `kappa_V^h=sba`, `kappa_V^phi=cba` | **VERIFIED** |
| Higgs-basis field/sign convention | **VERIFIED** |
| complete `(m_ij^2,lambda_i,beta)->(Y_i,Z_i)` map | **VERIFIED** |
| `Y3=-Z6v^2/2` | **VERIFIED** |
| exact alignment iff `Z6=0` | **VERIFIED** |
| exact `Z7(lambda_i,beta)` map | **VERIFIED** |
| `lambda7=0`, large `tan beta`: `Z7->-lambda6` | **VERIFIED** |
| exact physical `phi H+H-` coupling | **VERIFIED** |
| object-specific `X=lambda6 tan beta` translation | **VERIFIED** |

# C. Phase gates

\[
\boxed{\text{PHASE 0 PASS}}\quad
\boxed{\text{PHASE 1 PASS}}\quad
\boxed{\text{PHASE 2 PASS}}\quad
\boxed{\text{PHASE 3 PASS}}\quad
\boxed{\text{PHASE 4 PASS}}\quad
\boxed{\text{PHASE 5 PASS}}\quad
\boxed{\text{PHASE 6 PASS}}\quad
\boxed{\text{PHASE 7 PASS}}\quad
\boxed{\text{PHASE 9 PASS}}.
\]

---

# Phase 9 — Charged-Higgs trilinear extracted directly from the verified Higgs-basis potential

## What is established

The frozen Higgs basis is
\[
H_1=\begin{pmatrix}G^+\\(v+\rho_v+iG^0)/\sqrt2\end{pmatrix},\qquad
H_2=\begin{pmatrix}H^+\\(\rho_\perp+iA)/\sqrt2\end{pmatrix},
\]
with
\[
\rho_v=s_{\beta-\alpha}h+c_{\beta-\alpha}\phi,\qquad
\rho_\perp=c_{\beta-\alpha}h-s_{\beta-\alpha}\phi.
\]
At exact alignment, `sba=1,cba=0`, so `h=rho_v` and `phi=-rho_perp`.

The only Higgs-basis operators capable of producing a neutral CP-even field times `H+H-` at cubic order are
\[
Z_3(H_1^\dagger H_1)(H_2^\dagger H_2)
\]
and
\[
\{Z_7(H_2^\dagger H_2)H_1^\dagger H_2+\mathrm{h.c.}\}.
\]
Direct expansion gives
\[
\boxed{V\supset v\,(Z_3\rho_v+Z_7\rho_\perp)H^+H^-}.
\]
Therefore, before alignment,
\[
\boxed{C_V^{hH^+H^-}=v(Z_3s_{\beta-\alpha}+Z_7c_{\beta-\alpha})},
\]
\[
\boxed{C_V^{\phi H^+H^-}=v(Z_3c_{\beta-\alpha}-Z_7s_{\beta-\alpha})}.
\]
At exact alignment,
\[
\boxed{C_V^{hH^+H^-}=vZ_3},\qquad
\boxed{C_V^{\phi H^+H^-}=-vZ_7}.
\]

Here `C_V` means the coefficient of the displayed monomial in the scalar potential.

Because `L_int=-V_int`, the literal coefficient multiplying `phi H+H-` in the interaction Lagrangian is
\[
\boxed{C_{\mathcal L}^{\phi H^+H^-}=+vZ_7}
\]
in exact alignment.

If instead one defines a dimension-one coupling by
\[
\boxed{\mathcal L_{\rm int}\supset-g_{\phi H^+H^-}\,\phi H^+H^-},
\]
then
\[
\boxed{g_{\phi H^+H^-}=C_V^{\phi H^+H^-}=-vZ_7}.
\]

Since `phi`, `H+`, and `H-` are distinct external fields, there is no extra identical-particle factorial. The Feynman rule is
\[
\boxed{\phi H^+H^-:\quad -iC_V^{\phi H^+H^-}=+ivZ_7}
\]
in exact alignment.

## Derivation

Define
\[
A_1\equiv H_1^\dagger H_1,\qquad A_2\equiv H_2^\dagger H_2,\qquad C\equiv H_1^\dagger H_2.
\]
Keeping only pieces that can contain one neutral CP-even field and `H+H-`,
\[
A_1=\cdots+v\rho_v+\cdots,
\]
\[
A_2=H^-H^++\cdots,
\]
\[
C=\cdots+\frac v2(\rho_\perp+iA)+\cdots,
\qquad
C^\dagger=\cdots+\frac v2(\rho_\perp-iA)+\cdots.
\]
Thus
\[
Z_3A_1A_2\supset vZ_3\rho_vH^+H^-.
\]
For real `Z7`,
\[
Z_7A_2C+\mathrm{h.c.}=Z_7A_2(C+C^\dagger),
\]
and
\[
C+C^\dagger\supset v\rho_\perp.
\]
Hence
\[
Z_7A_2(C+C^\dagger)\supset vZ_7\rho_\perp H^+H^-.
\]
No other Higgs-basis term contributes at this field order: `Z2` has no VEV insertion; `Z4,Z5,Z6` require Goldstone/other charged fields to generate charged-scalar factors; quadratic terms do not generate a cubic scalar interaction.

Substitute the already-derived physical rotation:
\[
\rho_v=s h+c\phi,\qquad \rho_\perp=c h-s\phi,
\]
where `s=s_(beta-alpha)`, `c=c_(beta-alpha)`. Then
\[
V\supset v[(Z_3s+Z_7c)h+(Z_3c-Z_7s)\phi]H^+H^-.
\]
This proves the general result and its alignment limit without using a coupling table.

## Large-tan(beta) limit and lambda6

Phase 7 established, for `lambda7=0`,
\[
Z_7=-\lambda_6+\frac{\lambda_{345}-\lambda_1}{\tan\beta}+\mathcal O(\tan^{-2}\beta).
\]
Therefore the different trilinear objects behave as
\[
\boxed{C_V^{\phi H^+H^-}= -vZ_7
=+v\lambda_6-\frac{v(\lambda_{345}-\lambda_1)}{\tan\beta}+\cdots},
\]
\[
\boxed{C_{\mathcal L}^{\phi H^+H^-}=+vZ_7
=-v\lambda_6+\frac{v(\lambda_{345}-\lambda_1)}{\tan\beta}+\cdots},
\]
and, under `L_int=-g phi H+H-`,
\[
\boxed{g_{\phi H^+H^-}\simeq +v\lambda_6}.
\]
The Feynman rule tends to
\[
\boxed{-iv\lambda_6}.
\]

With the later project coordinate `X=lambda6 tan(beta)`, so `lambda6=X cot(beta)`,
\[
\boxed{C_V^{\phi H^+H^-}\simeq +vX\cot\beta},
\]
\[
\boxed{C_{\mathcal L}^{\phi H^+H^-}\simeq -vX\cot\beta},
\]
\[
\boxed{g_{\phi H^+H^-}\simeq +vX\cot\beta}\quad\text{if }\mathcal L_{\rm int}=-g\phi H^+H^-,
\]
and
\[
\boxed{\text{Feynman rule}\simeq-i vX\cot\beta}.
\]

## Convention map

GHOO18 explicitly states that its `q_i` is the coefficient of `H_i H+H-` in the potential. Its cubic-coupling appendix further states that potential coefficients become Feynman rules by multiplying by `-i` plus combinatorial factors for identical fields. For `H_i H+H-`, the three fields are distinct, so no extra factorial appears. This agrees with the object separation above.

Active 2HDMC implements
`c=-i v Re(q_{i1} Z3_local + q_{i2} Z7_local)` in `get_coupling_hhh`, after setting `Z7_local=-Lambda7_returned`. Phase 7 proved `Lambda7_returned=Z7_project`. Meanwhile `get_qki` uses second components `(-cba,+sba)` for `(h,H)`. Combining both sign layers gives
\[
-i v(Z_3s+Z_7c)
\]
for `h` and
\[
-i v(Z_3c-Z_7s)
\]
for `H=phi`, exactly the Feynman rules derived from the project potential. Thus the Phase-6 `Z7_local=-l7` mystery is resolved: it compensates the sign convention used in the second `qki` component and is not a different Higgs-basis `Z7`.

## Project-convention conflict discovered

Earlier project material used the shorthand
\[
g_{\phi H^+H^-}=vZ_7\simeq-vX\cot\beta.
\]
With the now-frozen physical state `phi=-rho_perp`, this equality is **not** the coefficient of `phi H+H-` in the potential and is **not** the `g` defined by `L_int=-g phi H+H-`.

It is, however, exactly the literal coefficient in `L_int`:
\[
C_{\mathcal L}^{\phi H^+H^-}=vZ_7.
\]
Therefore the old numerical expression can be retained only if its symbol is explicitly defined as the coefficient multiplying the monomial directly in `L_int`, not as the potential coefficient or as a `-g` convention. The master document must preserve this distinction.

## What was checked against the source

1. Direct symbolic differentiation of the full Higgs-basis potential gives
`d^3 V/(d rho_v dH+ dH-)=v Z3` and
`d^3 V/(d rho_perp dH+ dH-)=v Z7`.
2. Physical-state substitution gives `v(Z3 cba-Z7 sba)` for `phi`.
3. GHOO18 confirms that `q_i` denotes the coefficient in the potential and that the Feynman rule is obtained with `-i` (plus factorials only for identical particles).
4. Active 2HDMC reproduces exactly `-i` times the derived potential coefficient once its `qki` and local `Z7=-l7` conventions are combined.

## What remains uncertain

The analytic trilinear and the 2HDMC Feynman-rule translation are closed. What remains is a **project naming decision**: older plots/notes that call `vZ7` the coupling `g_phiH+H-` must be relabeled or explicitly defined as a Lagrangian monomial coefficient. No physics formula should be changed without tracking which convention those downstream loop-amplitude formulas assumed.

## Next smallest validation

Audit the charged-scalar contribution to `phi -> gamma gamma` and `phi -> Z gamma` starting from a declared interaction-Lagrangian convention. Determine whether the loop formula consumes `C_V`, `C_L`, `g` defined by `L=-g phiH+H-`, or the raw Feynman rule. Only then update any prior `g=vZ7` plotting convention.
