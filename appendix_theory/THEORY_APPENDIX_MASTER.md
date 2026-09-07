# 2HDM theory appendix — canonical derivation and convention audit

Issue: `fbientrigo/dihiggs#81`  
Status of this document: **canonical single-document audit record**.  
Current coverage: **Phases 0–5**. Future phases must be appended here as the primary narrative; phase-specific files may remain as auxiliary audit records.

The purpose of this document is not to collect formulas. It is to make every sign, normalization, field definition, basis transformation and physical conclusion traceable from the beginning.

---

# A. Audit rules and notation frozen before any physics conclusion

## A.1 Epistemic labels

- `[SOURCE]`: stated explicitly in an inspected primary source.
- `[DERIVED]`: obtained algebraically from already declared definitions.
- `[TRANSLATED]`: obtained by an explicit map between source conventions.
- `[IMPLEMENTATION-CHECKED]`: independently compared with the active code after the analytic result was frozen.
- `[PROJECT-DEFINITION]`: a convention or coordinate chosen by this project.
- `[OPEN-QUESTION]`: not yet promoted to an established result.

A source quotation is not treated as an independent derivation. Code is not allowed to choose the analytic formula it is later used to check.

## A.2 Objects that must never be conflated

The following are distinct unless an equation explicitly relates them:

\[
m_{22}^2,\qquad m_{12}^2,\qquad M^2,\qquad Y_2,\qquad Y_3.
\]

Similarly,

\[
\lambda_6,\lambda_7
\]

are generic-basis quartics and are not the Higgs-basis quantities

\[
Z_6,Z_7.
\]

The project coordinate

\[
X\equiv\lambda_6\tan\beta
\]

is a later project-defined coordinate. It is intentionally absent from the foundational derivation through Phase 5.

## A.3 Interaction-sign discipline

For every scalar interaction we distinguish three objects:

1. coefficient in the potential `V`;
2. coefficient in `L_int=-V_int`;
3. Feynman rule obtained from `i L_int`, including combinatorial factors.

For a Yukawa interaction we explicitly define

\[
\mathcal L_{\phi ff}
=-\frac{m_f}{v}\kappa_f^\phi\,\phi\bar f f.
\]

Thus a signed modifier `kappa` is meaningful only together with the declared scalar-field sign.

## A.4 Primary sources inspected

### S1 — Davidson & Haber, hep-ph/0504050

Primary TeX: `hbasis.tex`.  
Role: main analytic anchor because it contains the generic potential, vacuum definitions, minimization and Higgs-basis conventions in one internally consistent notation.

Archive SHA256:
`b6f133cc7875c71e3952386f56721f261aa51ccc9cfa1c8ee8aeeb17fbb8db01`.

### S2 — Branco et al., arXiv:1106.0034

Primary TeX: `PhysRep_large.tex`.  
Role: independent review-level cross-check, especially scalar Hessians and Type-I Yukawa structure. It explicitly states that its potential notation follows Davidson–Haber, but some displayed scalar-state signs in an early pedagogical subsection must be translated carefully.

Archive SHA256:
`6e3693b8fc46a8f26dc546a4389ab3f4f5ead704fbcf0ee5e61bd8524ed19481`.

### S3 — Grzadkowski, Haber, Ogreid, Osland, arXiv:1808.01472

Primary TeX: `paper_heavyhiggs_jhep_revised3.tex`.  
Role: Higgs-basis/alignment and cubic/Yukawa cross-check, with a different quadratic-symbol convention that must be translated explicitly.

Archive SHA256:
`03c4b3f6a32c30beb4f62645ea941c0735c38fd09b5d99f215e8ce229d3d874b`.

---

# Phase 0 — Source and convention inventory

## What is established

[SOURCE][PROJECT-DEFINITION] Davidson–Haber is used as the main generic-basis sign convention.

Its scalar potential uses

\[
+m_{11}^2\Phi_1^\dagger\Phi_1
+m_{22}^2\Phi_2^\dagger\Phi_2
-[m_{12}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}]
\]

with the quartic normalization written explicitly in Phase 2.

[SOURCE] Branco uses the same operator normalization for the generic potential.

[TRANSLATED] GHOO18 instead writes its quadratic terms as

\[
-\frac12\{m_{11,G}^2\Phi_1^\dagger\Phi_1
+m_{22,G}^2\Phi_2^\dagger\Phi_2
+[m_{12,G}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}]\},
\]

so operator matching gives

\[
\boxed{m_{11,\rm DH}^2=-\frac12m_{11,G}^2},
\qquad
\boxed{m_{22,\rm DH}^2=-\frac12m_{22,G}^2},
\qquad
\boxed{m_{12,\rm DH}^2=+\frac12m_{12,G}^2}.
\]

No formula from GHOO can therefore be copied symbol-for-symbol into the DH convention without this translation.

## Derivation

The quadratic map follows by matching coefficients of the independent operators

\[
\Phi_1^\dagger\Phi_1,
\quad
\Phi_2^\dagger\Phi_2,
\quad
\Phi_1^\dagger\Phi_2+\mathrm{h.c.}
\]

between the two potentials.

## Convention map

DH05 uses hypercharge `Y=1` with

\[
Q=T_3+\frac{Y}{2},
\]

which is equivalent to the modern convention `Y=1/2` with `Q=T_3+Y`.

DH05 CP-even source states are frozen as

\[
h_{\rm DH}=-\rho_1\sin\alpha+\rho_2\cos\alpha,
\qquad
H_{\rm DH}=\rho_1\cos\alpha+\rho_2\sin\alpha.
\]

The early simple display in Branco is the global-sign reversal of these two fields. This does not change masses but does change every interaction containing an odd number of those fields unless all signs are transformed consistently.

## What was checked against the source

- DH05 generic potential, VEV definition, `tan beta`, Higgs basis and CP-even rotation.
- Branco generic potential and Type-I model table.
- GHOO generic potential and Higgs-basis potential.

## What remains uncertain

At this stage no project `h/phi` state map, no Yukawa modifier, no gauge modifier and no trilinear coupling is promoted.

## Next smallest validation

Define the doublets and reconstruct the complete renormalizable potential before introducing physical states.

---

# Phase 1 — Fields, VEVs and the beta coordinate

## What is established

[DERIVED][SOURCE] The two scalar doublets are written

\[
\boxed{
\Phi_i=
\begin{pmatrix}
\phi_i^+\\[1mm]
\dfrac{v_i+\rho_i+i\eta_i}{\sqrt2}
\end{pmatrix}},
\qquad i=1,2.
\]

We select a neutral CP-conserving vacuum with real non-negative VEVs,

\[
\langle\Phi_1\rangle=\frac1{\sqrt2}\binom0{v_1},
\qquad
\langle\Phi_2\rangle=\frac1{\sqrt2}\binom0{v_2}.
\]

Define

\[
\boxed{v^2=v_1^2+v_2^2},
\qquad
\boxed{\tan\beta=\frac{v_2}{v_1}},
\]

and therefore

\[
v_1=v c_\beta,
\qquad
v_2=v s_\beta.
\]

No physical `H+`, `A`, `h` or `phi` is assumed yet.

## Derivation

The `1/sqrt(2)` normalization makes the real neutral fluctuations canonically normalized. The charged component has electric charge +1 and the lower component is neutral under the chosen hypercharge convention.

## Convention map

The project uses `rho_i` for CP-even real neutral fluctuations and `eta_i` for CP-odd real neutral fluctuations.

## What was checked against the source

DH05 `potmin` and `tanbdef`, Branco component expansion, GHOO VEV definitions.

## What remains uncertain

No statement about which linear combination is physical or SM-like.

## Next smallest validation

Construct the most general CP-conserving renormalizable scalar potential from gauge invariants.

---

# Phase 2 — Generic CP-conserving scalar potential

## What is established

The independent gauge-invariant bilinears are

\[
B_{ij}=\Phi_i^\dagger\Phi_j.
\]

Hermiticity and renormalizability lead to the DH/Branco potential

\[
\boxed{
\begin{aligned}
V={}&m_{11}^2\Phi_1^\dagger\Phi_1
+m_{22}^2\Phi_2^\dagger\Phi_2
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

For the CP-conserving specialization used in the current project, all displayed parameters are real.

The campaign restriction `lambda7=0` is **not** imposed here; it is a later parameter-space choice.

## Derivation

Dimension-two gauge singlets are `B11`, `B22`, `B12` and `B21=B12^dagger`. Dimension-four terms are products of these bilinears. Hermiticity fixes which coefficients are real and which terms require an explicit hermitian conjugate. CP conservation permits a basis where the remaining complex coefficients are real.

## Convention map

The exact DH/Branco operator normalization is frozen as the project generic-basis potential. GHOO quadratic symbols obey the Phase-0 factor/sign map.

## What was checked against the source

DH05 Eq. `pot` and Branco Eq. `2_VH1` agree term by term. The active 2HDMC generic-parameter path later reproduces stationarity equations derived from this potential.

## What remains uncertain

A potential definition alone does not identify the vacuum or physical masses.

## Next smallest validation

Evaluate `V` on the neutral vacuum and differentiate it rather than importing minimization formulas.

---

# Phase 3 — Neutral vacuum and minimization

## What is established

Define

\[
\lambda_{345}\equiv\lambda_3+\lambda_4+\lambda_5.
\]

Direct substitution of the neutral real VEVs gives

\[
\boxed{
\begin{aligned}
V_0={}&\frac12m_{11}^2v_1^2
+\frac12m_{22}^2v_2^2
-m_{12}^2v_1v_2\\
&+\frac18\lambda_1v_1^4
+\frac18\lambda_2v_2^4
+\frac14\lambda_{345}v_1^2v_2^2\\
&+\frac12\lambda_6v_1^3v_2
+\frac12\lambda_7v_1v_2^3.
\end{aligned}}
\]

The two neutral stationarity equations are

\[
\boxed{
0=m_{11}^2v_1-m_{12}^2v_2
+\frac12\lambda_1v_1^3
+\frac12\lambda_{345}v_1v_2^2
+\frac32\lambda_6v_1^2v_2
+\frac12\lambda_7v_2^3
}
\]

and

\[
\boxed{
0=m_{22}^2v_2-m_{12}^2v_1
+\frac12\lambda_2v_2^3
+\frac12\lambda_{345}v_1^2v_2
+\frac12\lambda_6v_1^3
+\frac32\lambda_7v_1v_2^2
}.
\]

The factors `3 lambda6` and `3 lambda7` are produced directly by differentiating the cubic powers of the VEVs; they are not convention guesses.

Solving for the diagonal quadratic coefficients,

\[
\boxed{
\begin{aligned}
m_{11}^2={}&m_{12}^2\frac{v_2}{v_1}
-\frac12\left[\lambda_1v_1^2+\lambda_{345}v_2^2
+3\lambda_6v_1v_2+\lambda_7\frac{v_2^3}{v_1}\right],\\
m_{22}^2={}&m_{12}^2\frac{v_1}{v_2}
-\frac12\left[\lambda_2v_2^2+\lambda_{345}v_1^2
+\lambda_6\frac{v_1^3}{v_2}+3\lambda_7v_1v_2\right].
\end{aligned}}
\]

## Derivation

Each factor in `V0` follows from

\[
\Phi_1^\dagger\Phi_1=v_1^2/2,
\quad
\Phi_2^\dagger\Phi_2=v_2^2/2,
\quad
\Phi_1^\dagger\Phi_2=v_1v_2/2.
\]

For example, the `lambda6` term plus h.c. is

\[
2\lambda_6\left(\frac{v_1^2}{2}\right)\left(\frac{v_1v_2}{2}\right)
=\frac12\lambda_6v_1^3v_2,
\]

whose derivative with respect to `v1` is

\[
\frac32\lambda_6v_1^2v_2.
\]

The neutral-field tadpoles satisfy

\[
\left.\frac{\partial V}{\partial\rho_i}\right|_0
=\frac{\partial V_0}{\partial v_i},
\]

because the neutral real fields enter through `v_i+rho_i`.

### Derived soft coordinate M^2

Only after minimization define

\[
\boxed{M^2\equiv\frac{m_{12}^2}{s_\beta c_\beta}}.
\]

Then

\[
\boxed{
\begin{aligned}
m_{11}^2={}&M^2s^2
-\frac{v^2}{2}\left[\lambda_1c^2+\lambda_{345}s^2+3\lambda_6sc+\lambda_7\frac{s^3}{c}\right],\\
m_{22}^2={}&M^2c^2
-\frac{v^2}{2}\left[\lambda_2s^2+\lambda_{345}c^2+\lambda_6\frac{c^3}{s}+3\lambda_7sc\right].
\end{aligned}}
\]

Therefore

\[
\boxed{m_{22}^2\neq m_{12}^2\neq M^2}
\]

generically.

## Convention map

GHOO stationarity equations can be translated with the Phase-0 quadratic-symbol map. `M^2` is a project/literature shorthand derived from the off-diagonal soft term and VEV orientation; it is not a diagonal mass parameter.

## What was checked against the source

- DH05 minimization equations agree after setting the phase to zero and parameters real.
- A symbolic differentiation audit reproduces both tadpoles.
- Active 2HDMC `set_param_gen` reproduces the derived `m22^2` equation term by term.

## What remains uncertain

Stationarity does not prove a local or global minimum.

## Next smallest validation

Compute the charged, CP-odd and CP-even Hessians from the same potential and prove the Goldstone zero modes.

---

# Phase 4 — Hessians, Goldstones, masses and scalar-state signs

## What is established

The physical charged and CP-odd directions are not assumed. They emerge because the post-tadpole Hessians have a common rank-one geometry.

### Charged sector before tadpoles

\[
\mathcal M_\pm^2=
\begin{pmatrix}X_\pm&Y_\pm\\Y_\pm&Z_\pm\end{pmatrix}
\]

with

\[
\begin{aligned}
X_\pm&=m_{11}^2+\frac12\lambda_1v_1^2+\frac12\lambda_3v_2^2+\lambda_6v_1v_2,\\
Z_\pm&=m_{22}^2+\frac12\lambda_2v_2^2+\frac12\lambda_3v_1^2+\lambda_7v_1v_2,\\
Y_\pm&=-m_{12}^2+\frac12(\lambda_4+\lambda_5)v_1v_2
+\frac12\lambda_6v_1^2+\frac12\lambda_7v_2^2.
\end{aligned}
\]

After substituting the Phase-3 tadpoles,

\[
\boxed{
\mathcal M_\pm^2
=D_\pm
\begin{pmatrix}
v_2/v_1&-1\\-1&v_1/v_2
\end{pmatrix}}
\]

where

\[
D_\pm=m_{12}^2-\frac12[(\lambda_4+\lambda_5)v_1v_2+\lambda_6v_1^2+\lambda_7v_2^2].
\]

### CP-odd sector before tadpoles

\[
\begin{aligned}
(\mathcal M_A^2)_{11}&=m_{11}^2+\frac12\lambda_1v_1^2
+\frac12(\lambda_3+\lambda_4-\lambda_5)v_2^2+\lambda_6v_1v_2,\\
(\mathcal M_A^2)_{22}&=m_{22}^2+\frac12\lambda_2v_2^2
+\frac12(\lambda_3+\lambda_4-\lambda_5)v_1^2+\lambda_7v_1v_2,\\
(\mathcal M_A^2)_{12}&=-m_{12}^2+\lambda_5v_1v_2
+\frac12\lambda_6v_1^2+\frac12\lambda_7v_2^2.
\end{aligned}
\]

After tadpoles,

\[
\boxed{
\mathcal M_A^2
=D_A
\begin{pmatrix}
v_2/v_1&-1\\-1&v_1/v_2
\end{pmatrix}}
\]

with

\[
D_A=m_{12}^2-\lambda_5v_1v_2-\frac12(\lambda_6v_1^2+\lambda_7v_2^2).
\]

### Goldstone proof

For either sector,

\[
\begin{pmatrix}
v_2/v_1&-1\\-1&v_1/v_2
\end{pmatrix}
\binom{v_1}{v_2}=0.
\]

Hence the vacuum direction is the zero mode, while the orthogonal vector `(-v2,v1)` is physical. Therefore

\[
\boxed{G^+=c_\beta\phi_1^++s_\beta\phi_2^+},
\qquad
\boxed{H^+=-s_\beta\phi_1^++c_\beta\phi_2^+},
\]

\[
\boxed{G^0=c_\beta\eta_1+s_\beta\eta_2},
\qquad
\boxed{A=-s_\beta\eta_1+c_\beta\eta_2}.
\]

Thus `beta` diagonalizes these sectors because it describes the vacuum orientation.

### Physical masses

The nonzero eigenvalues are

\[
\boxed{
m_A^2=M^2-\frac{v^2}{2}
\left(2\lambda_5+\lambda_6\cot\beta+\lambda_7\tan\beta\right)}
\]

and

\[
\boxed{
m_{H^\pm}^2=M^2-\frac{v^2}{2}
\left(\lambda_4+\lambda_5+\lambda_6\cot\beta+\lambda_7\tan\beta\right)}.
\]

Therefore

\[
\boxed{m_{H^\pm}^2-m_A^2=\frac{v^2}{2}(\lambda_5-\lambda_4)}.
\]

The `lambda6` and `lambda7` dependence cancels exactly in this splitting.

## Derivation — CP-even sector

Before tadpole elimination,

\[
\boxed{
\mathcal M_\rho^2=
\begin{pmatrix}\mathcal M_{11}^2&\mathcal M_{12}^2\\\mathcal M_{12}^2&\mathcal M_{22}^2\end{pmatrix}}
\]

with

\[
\begin{aligned}
\mathcal M_{11}^2&=m_{11}^2+\frac32\lambda_1v_1^2+\frac12\lambda_{345}v_2^2+3\lambda_6v_1v_2,\\
\mathcal M_{22}^2&=m_{22}^2+\frac32\lambda_2v_2^2+\frac12\lambda_{345}v_1^2+3\lambda_7v_1v_2,\\
\mathcal M_{12}^2&=-m_{12}^2+\lambda_{345}v_1v_2+\frac32\lambda_6v_1^2+\frac32\lambda_7v_2^2.
\end{aligned}
\]

After tadpoles,

\[
\boxed{
\begin{aligned}
\mathcal M_{11}^2&=M^2s^2+v^2\left[\lambda_1c^2+\frac32\lambda_6sc-\frac12\lambda_7\frac{s^3}{c}\right],\\
\mathcal M_{22}^2&=M^2c^2+v^2\left[\lambda_2s^2-\frac12\lambda_6\frac{c^3}{s}+\frac32\lambda_7sc\right],\\
\mathcal M_{12}^2&=-M^2sc+v^2\left[\lambda_{345}sc+\frac32\lambda_6c^2+\frac32\lambda_7s^2\right].
\end{aligned}}
\]

An implementation-friendly equivalent form is

\[
\boxed{
\begin{aligned}
\mathcal M_{11}^2&=m_A^2s^2+v^2(\lambda_1c^2+\lambda_5s^2+2\lambda_6sc),\\
\mathcal M_{22}^2&=m_A^2c^2+v^2(\lambda_2s^2+\lambda_5c^2+2\lambda_7sc),\\
\mathcal M_{12}^2&=-m_A^2sc+v^2[(\lambda_3+\lambda_4)sc+\lambda_6c^2+\lambda_7s^2].
\end{aligned}}
\]

### Alpha rotation

Freeze the DH CP-even sign convention

\[
\boxed{
\begin{pmatrix}h_{\rm DH}\\H_{\rm DH}\end{pmatrix}
=
\begin{pmatrix}-s_\alpha&c_\alpha\\c_\alpha&s_\alpha\end{pmatrix}
\begin{pmatrix}\rho_1\\\rho_2\end{pmatrix}}.
\]

The off-diagonal element after rotation vanishes when

\[
\boxed{\tan2\alpha=\frac{2\mathcal M_{12}^2}{\mathcal M_{11}^2-\mathcal M_{22}^2}},
\]

with the quadrant determined by the eigenvector branch rather than the tangent alone.

The eigenvalues are

\[
\boxed{
m_{h,H}^2=\frac12\left[
\mathcal M_{11}^2+\mathcal M_{22}^2
\mp\sqrt{(\mathcal M_{11}^2-\mathcal M_{22}^2)^2+4(\mathcal M_{12}^2)^2}
\right]}.
\]

### Vacuum-aligned CP-even direction and project state signs

Define

\[
\rho_v=c_\beta\rho_1+s_\beta\rho_2,
\qquad
\rho_\perp=-s_\beta\rho_1+c_\beta\rho_2.
\]

Then

\[
\boxed{h_{\rm DH}=s_{\beta-\alpha}\rho_v+c_{\beta-\alpha}\rho_\perp},
\]

\[
\boxed{H_{\rm DH}=c_{\beta-\alpha}\rho_v-s_{\beta-\alpha}\rho_\perp}.
\]

The project branch uses

\[
\sin(\beta-\alpha)=1.
\]

Therefore

\[
h_{\rm DH}=\rho_v,
\qquad
H_{\rm DH}=-\rho_\perp.
\]

Freeze the project names

\[
\boxed{h\equiv h_{\rm DH}},
\qquad
\boxed{\phi\equiv H_{\rm DH}}.
\]

Hence, in exact alignment,

\[
\boxed{h=c_\beta\rho_1+s_\beta\rho_2},
\qquad
\boxed{\phi=s_\beta\rho_1-c_\beta\rho_2}.
\]

This sign of `phi` is inherited from the DH field convention. It is not chosen to force a later Yukawa or scalar-trilinear sign.

## Convention map

2HDMC's active physical-input branch uses `alpha=beta-asin(sba)` on its non-negative `cba` branch, consistent with `alpha=beta-pi/2` at exact alignment.

An early pedagogical restricted subsection of Branco contains charged/CP-odd displayed mass-term factors that do not match the Hessian of its printed potential. Its later general Hessian agrees with the independent derivation and is used as the source cross-check. The discrepancy is recorded rather than silently repaired.

## What was checked against the source

DH Goldstone/charged rotations, DH CP-even states, GHOO beta rotations, Branco general scalar Hessian and active 2HDMC mass matrices.

## What remains uncertain

- Gauge interpretation of `rho_v` still requires the kinetic-term derivation.
- Exact alignment versus `Z6=0` still requires the Higgs-basis mass matrix.
- Positive physical squared masses are a local quadratic condition, not a global-minimum proof.

## Next smallest validation

Derive Type-I Yukawa modifiers from the Yukawa Lagrangian using the state signs just frozen.

---

# Phase 5 — Type-I Yukawa couplings

## What is established

Type I means all charged fermions couple to `Phi2` and not `Phi1`. In flavor-matrix notation define the sign convention

\[
\boxed{
-\mathcal L_Y^{\rm I}
=\overline Q_LY_d\Phi_2d_R
+\overline Q_LY_u\widetilde\Phi_2u_R
+\overline L_LY_\ell\Phi_2\ell_R
+\mathrm{h.c.}}
\]

with

\[
\widetilde\Phi_2=i\sigma_2\Phi_2^*.
\]

Using

\[
\Phi_2^0=\frac{v_2+\rho_2+i\eta_2}{\sqrt2},
\qquad
\widetilde\Phi_2^0=\frac{v_2+\rho_2-i\eta_2}{\sqrt2},
\]

the CP-even fluctuation `rho2` enters every charged-fermion sector with the same sign.

After diagonalizing the fermion mass matrices,

\[
\boxed{m_f=\frac{y_fv_2}{\sqrt2}}
\]

and therefore

\[
\boxed{
\mathcal L_Y^{\rm CP-even}
=-\sum_f\frac{m_f}{v_2}\rho_2\bar f f
=-\sum_f\frac{m_f}{v s_\beta}\rho_2\bar f f}.
\]

## Derivation

Phase 4 fixed

\[
\begin{pmatrix}h\\\phi\end{pmatrix}
=R_\alpha\begin{pmatrix}\rho_1\\\rho_2\end{pmatrix},
\qquad
R_\alpha=
\begin{pmatrix}-s_\alpha&c_\alpha\\c_\alpha&s_\alpha\end{pmatrix}.
\]

Since

\[
R_\alpha^TR_\alpha=I,
\qquad
R_\alpha^T=R_\alpha,
\]

the inverse transformation is

\[
\boxed{\rho_2=c_\alpha h+s_\alpha\phi}.
\]

Substitute into the interaction:

\[
\boxed{
\mathcal L_Y^{\rm CP-even}
=-\sum_f\frac{m_f}{v}
\left(
\frac{c_\alpha}{s_\beta}h
+\frac{s_\alpha}{s_\beta}\phi
\right)\bar f f}.
\]

Define

\[
\mathcal L_Y^{\rm CP-even}
=-\sum_f\frac{m_f}{v}
(\kappa_f^h h+\kappa_f^\phi\phi)\bar f f.
\]

Then

\[
\boxed{\kappa_f^h=\frac{c_\alpha}{s_\beta}},
\qquad
\boxed{\kappa_f^\phi=\frac{s_\alpha}{s_\beta}},
\qquad f=u,d,\ell.
\]

### Exact alignment

On the frozen branch

\[
\alpha=\beta-\frac\pi2,
\]

so

\[
c_\alpha=s_\beta,
\qquad
s_\alpha=-c_\beta.
\]

Therefore

\[
\boxed{\kappa_f^h=1},
\qquad
\boxed{\kappa_f^\phi=-\cot\beta},
\qquad f=u,d,\ell.
\]

The minus sign follows from the scalar state convention. Equivalently, exact alignment gave

\[
\phi=s_\beta\rho_1-c_\beta\rho_2,
\]

so the `Phi2` projection of `phi` is `-c_beta`, while the mass-generating VEV fraction is `s_beta`.

### Meaning of the signed coupling

With the declared convention,

\[
\mathcal L_{\phi ff}
=-\frac{m_f}{v}\kappa_f^\phi\phi\bar f f
=+\frac{m_f}{v}\cot\beta\,\phi\bar f f
\]

at exact alignment.

A global field redefinition `phi -> -phi` would flip every interaction containing an odd number of `phi` fields. Therefore a signed coupling should never be quoted without the field convention.

## Convention map

- Branco Type-I model table: `u_R,d_R,e_R` all couple to `Phi2`.
- Branco Yukawa table gives `xi_h=c_alpha/s_beta` and `xi_H=s_alpha/s_beta` for all charged fermions; this matches the derivation after adopting the DH-sign state convention.
- GHOO Type-I appendix sets its `eta_1^{u,d,l}=0` and its CP-even neutral couplings are controlled by the `Phi2` projection `R_{j2}/s_beta`.
- The absolute sign of a GHOO heavy state must be translated with its scalar field sign; only the dependence and projection structure are used as a cross-check here.

## What was checked against the source

The independent result was compared only after derivation with:

1. Branco `tab:3_models` Type-I assignment;
2. Branco `Eq:Yukawa` overall `-m_f/v` normalization;
3. Branco `tab:3_couplings` CP-even modifiers;
4. GHOO Appendix `Yuk_Type_I` and its `R_{j2}/s_beta` neutral structure;
5. active project use of `2HDMC::set_yukawas_type(1)` as an implementation consistency check.

A symbolic audit verifies that the chosen rotation is orthogonal/self-inverse and that `alpha=beta-pi/2` gives exactly `1` and `-cot(beta)`.

## What remains uncertain

- `kappa_V^phi=0` is still **not** inferred; it must come from gauge kinetic terms.
- `Z6=0` and its relation to exact alignment remain for the Higgs-basis phase.
- The project does not yet promote any `phi H+H-` trilinear formula.
- Pseudoscalar and charged-Higgs Yukawa vertices are derivable from the same Lagrangian but are not needed for the current high-risk claim.

## Next smallest validation

Start from

\[
\sum_{i=1}^2(D_\mu\Phi_i)^\dagger(D^\mu\Phi_i)
\]

and derive the terms linear in `rho_i` and quadratic in `W/Z`. This must show independently which CP-even direction couples to vector-boson pairs and determine `kappa_V^h` and `kappa_V^phi`.

---

# B. Validation status after Phase 5

| Claim | Status after Phase 5 |
|---|---|
| Generic scalar potential and signs | **VERIFIED** |
| Neutral tadpoles/minimization | **VERIFIED** |
| `m22^2`, `m12^2`, `M^2` distinction | **VERIFIED** |
| Charged/CP-odd Goldstone rotations | **VERIFIED** |
| `m_A^2`, `m_H+^2`, CP-even Hessian | **VERIFIED** |
| project scalar-state convention `h,phi` | **VERIFIED** |
| Type-I `kappa_f^h=c_alpha/s_beta` | **VERIFIED** |
| Type-I exact-alignment `kappa_f^phi=-cot beta` | **VERIFIED** |
| exact-alignment `kappa_V^phi=0` | **NOT YET VERIFIED** |
| `Y3=-Z6 v^2/2` | **SOURCE-VERIFIED; independent Higgs-basis derivation pending** |
| exact alignment iff `Z6=0` | **PARTIALLY VERIFIED; derivation pending** |
| exact `phi H+H-` trilinear | **NOT VERIFIED** |
| generic-to-Higgs-basis `Z7` map | **NOT VERIFIED independently** |
| large-`tan beta` trilinear approximations | **NOT VERIFIED** |
| `X=lambda6 tan beta` interpretation | **DEFERRED by design** |

---

# C. Phase gates

\[
\boxed{\text{PHASE 0 PASS}}
\quad
\boxed{\text{PHASE 1 PASS}}
\quad
\boxed{\text{PHASE 2 PASS}}
\quad
\boxed{\text{PHASE 3 PASS}}
\quad
\boxed{\text{PHASE 4 PASS}}
\quad
\boxed{\text{PHASE 5 PASS}}.
\]

The next phase is intentionally narrow: gauge couplings from kinetic terms. No Higgs-basis or trilinear result should be promoted before that derivation is complete.
