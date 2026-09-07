# 2HDM theory appendix — canonical derivation and convention audit

Issue: `fbientrigo/dihiggs#81`  
Status: **canonical single-document audit record**.  
Current coverage: **Phases 0–6**.

This is the document that must remain sufficient for an independent reviewer to reconstruct the theory chain without trusting project code or a formula copied from the literature. Phase-specific Markdown files and symbolic scripts are auxiliary audits; this file is the primary narrative.

The governing order is

\[
\boxed{
\text{fields}
\to\text{potential}
\to\text{vacuum}
\to\text{stationarity}
\to\text{mass matrices}
\to\text{physical states}
\to\text{Yukawa/gauge couplings}
\to\text{Higgs basis}
\to\text{scalar trilinears}
\to\text{loop amplitudes}
\to X
}.
\]

No later project coordinate or phenomenological expectation is allowed to determine an earlier sign or normalization.

---

# A. Audit rules, sources and conventions frozen before deriving physics

## A.1 Epistemic labels

Every nontrivial statement is conceptually assigned one or more of:

- `[SOURCE]`: explicitly stated in an inspected source.
- `[DERIVED]`: obtained algebraically from definitions already fixed in this document.
- `[TRANSLATED]`: obtained through an explicit convention map.
- `[IMPLEMENTATION-CHECKED]`: compared with project/2HDMC code only after the analytic result was frozen.
- `[PROJECT-DEFINITION]`: a convention or coordinate chosen by this project.
- `[OPEN-QUESTION]`: not yet established strongly enough to use downstream.

Source agreement alone is not an independent derivation when the object can be reconstructed from first principles. Code is never allowed to choose the analytic formula it is later used to check.

## A.2 Objects that must not be conflated

The following objects are distinct unless an equation explicitly connects them:

\[
m_{22}^2,\qquad m_{12}^2,\qquad M^2,\qquad Y_2,\qquad Y_3.
\]

Likewise,

\[
\lambda_6,\lambda_7
\]

are generic-basis quartics and are **not** the Higgs-basis quantities

\[
Z_6,Z_7.
\]

For scalar trilinears, three objects will later be kept separate:

1. coefficient in the potential `V`;
2. coefficient in `L_int=-V_int`;
3. Feynman rule including `-i` and identical-particle factors.

The project coordinate

\[
X\equiv\lambda_6\tan\beta
\]

is a later `[PROJECT-DEFINITION]`. It is intentionally absent from the foundational derivation through Phase 6.

## A.3 Primary sources used as convention anchors

The primary analytic anchor is Davidson–Haber (DH05, hep-ph/0504050), because it provides the generic potential, VEVs, minimization, scalar rotations and Higgs-basis notation in one internally connected convention.

Branco et al. (BFLRS11, arXiv:1106.0034) is used as an independent broad review and for Type-I Yukawa structure. It explicitly warns that signs/factors of symbols differ across 2HDM papers.

Grzadkowski–Haber–Ogreid–Osland (GHOO18, arXiv:1808.01472) is used mainly for Higgs-basis/alignment and later trilinear translation. Its generic quadratic symbols use a different normalization from DH05 and must not be copied symbol-for-symbol.

The active vendored 2HDMC implementation is an implementation check, never the primary analytic source.

## A.4 Generic-basis potential convention

We freeze the DH05/BFLRS11 generic-basis convention:

\[
\begin{aligned}
V={}&m_{11}^2\Phi_1^\dagger\Phi_1
+m_{22}^2\Phi_2^\dagger\Phi_2
-\left[m_{12}^2\Phi_1^\dagger\Phi_2+\text{h.c.}\right]\\
&+\frac12\lambda_1(\Phi_1^\dagger\Phi_1)^2
+\frac12\lambda_2(\Phi_2^\dagger\Phi_2)^2\\
&+\lambda_3(\Phi_1^\dagger\Phi_1)(\Phi_2^\dagger\Phi_2)
+\lambda_4(\Phi_1^\dagger\Phi_2)(\Phi_2^\dagger\Phi_1)\\
&+\left\{
\frac12\lambda_5(\Phi_1^\dagger\Phi_2)^2
+\left[\lambda_6(\Phi_1^\dagger\Phi_1)
+\lambda_7(\Phi_2^\dagger\Phi_2)\right]\Phi_1^\dagger\Phi_2
+\text{h.c.}
\right\}.
\end{aligned}
\]

For the CP-conserving branch used here all displayed coefficients are taken real. `lambda7=0` is a campaign specialization and is **not** part of the model definition.

GHOO18 instead writes its generic quadratic terms with an overall `-1/2`. Operator matching gives

\[
\boxed{m_{11,\rm DH}^2=-\frac12m_{11,\rm G}^2},
\qquad
\boxed{m_{22,\rm DH}^2=-\frac12m_{22,\rm G}^2},
\qquad
\boxed{m_{12,\rm DH}^2=+\frac12m_{12,\rm G}^2}.
\]

This map is `[DERIVED][TRANSLATED]` by matching coefficients of the same operators.

---

# Phase 1 — Fields, electroweak quantum numbers and vacuum coordinates

## What is established

We begin with two complex `SU(2)_L` scalar doublets with identical hypercharge. In modern notation

\[
Q=T^3+Y,
\qquad
Y(\Phi_i)=\frac12.
\]

DH05 uses the equivalent convention `Q=T3+Y/2` and calls the doublet hypercharge `Y=1`.

The component expansion is fixed as

\[
\boxed{
\Phi_i=
\begin{pmatrix}
\phi_i^+\\[1mm]
(v_i+\rho_i+i\eta_i)/\sqrt2
\end{pmatrix}}
\]

with a real, neutral, CP-conserving vacuum

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

and

\[
c_\beta\equiv\frac{v_1}{v},
\qquad
s_\beta\equiv\frac{v_2}{v}.
\]

No `h`, `phi`, `A` or `H^\pm` state is assumed at this stage.

## Derivation

The factor `1/sqrt(2)` makes `rho_i` and `eta_i` canonically normalized real fields. The upper component has charge `+1` and the lower component charge `0` because, for a `Y=1/2` doublet,

\[
Q_{\rm upper}=+\frac12+\frac12=1,
\qquad
Q_{\rm lower}=-\frac12+\frac12=0.
\]

The vacuum is chosen in the neutral lower components so electromagnetism remains unbroken.

## Convention map

`v1,v2 >=0` and positive `tan beta` define the coordinate patch used by the project/2HDMC. A relative CP phase is absent because the present branch is CP conserving.

## What was checked against the source

DH05 gives the same VEV normalization and `tan beta=v2/v1`; BFLRS11 uses the same component expansion. No physical-state convention was imported yet.

## What remains uncertain

Nothing needed for the next phase. The question of whether a stationary point is the global electroweak vacuum is deliberately postponed.

## Next smallest validation

Enumerate the complete renormalizable gauge-singlet scalar operator basis and reconstruct the potential before invoking any physical masses.

---

# Phase 2 — Scalar potential reconstructed from the operator basis

## What is established

The complete renormalizable scalar potential in the frozen generic basis is the expression in Section A.4.

## Derivation

The gauge-invariant bilinears are

\[
B_{ij}\equiv\Phi_i^\dagger\Phi_j,
\qquad B_{ji}=B_{ij}^\dagger.
\]

Dimension-two Hermitian invariants give `B11`, `B22` and the complex off-diagonal `B12` plus h.c. Dimension-four terms are all products of two such bilinears. Hermiticity gives the displayed `lambda1...lambda7` structure. CP conservation allows a real basis for `m12^2,lambda5,lambda6,lambda7` in the branch studied here.

The signs and factors are therefore not arbitrary conventions after they are frozen:

- `+m11^2 B11`;
- `+m22^2 B22`;
- `-[m12^2 B12+h.c.]`;
- `1/2` multiplying `lambda1`, `lambda2` and the displayed `lambda5` monomial;
- `lambda6,lambda7` appear before adding h.c.

## Convention map

DH05 and BFLRS11 match term-for-term. GHOO18 requires the quadratic translation written in A.4; its quartic normalization is otherwise compatible with the displayed operator form.

## What was checked against the source

After reconstructing the operator content independently, it was checked against DH05 Eq. `pot` and BFLRS11 Eq. `2_VH1`. The active 2HDMC generic-input path later reproduced the independently derived stationarity relation, providing an implementation-level consistency check.

## What remains uncertain

Nothing blocks minimization. Boundedness/global-minimum questions are separate from the algebraic stationarity conditions.

## Next smallest validation

Substitute the neutral VEVs term by term and differentiate the resulting vacuum potential.

---

# Phase 3 — Vacuum potential and minimization

Define

\[
\lambda_{345}\equiv\lambda_3+\lambda_4+\lambda_5.
\]

## What is established

Direct substitution of the neutral real vacuum gives

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

The stationarity equations are

\[
\boxed{
0=m_{11}^2v_1-m_{12}^2v_2
+\frac12\lambda_1v_1^3
+\frac12\lambda_{345}v_1v_2^2
+\frac32\lambda_6v_1^2v_2
+\frac12\lambda_7v_2^3}
\]

and

\[
\boxed{
0=m_{22}^2v_2-m_{12}^2v_1
+\frac12\lambda_2v_2^3
+\frac12\lambda_{345}v_1^2v_2
+\frac12\lambda_6v_1^3
+\frac32\lambda_7v_1v_2^2}.
\]

Solving for the diagonal quadratic coefficients,

\[
\boxed{
\begin{aligned}
m_{11}^2={}&m_{12}^2\frac{v_2}{v_1}
-\frac12\left[
\lambda_1v_1^2+\lambda_{345}v_2^2
+3\lambda_6v_1v_2
+\lambda_7\frac{v_2^3}{v_1}
\right],\\
m_{22}^2={}&m_{12}^2\frac{v_1}{v_2}
-\frac12\left[
\lambda_2v_2^2+\lambda_{345}v_1^2
+\lambda_6\frac{v_1^3}{v_2}
+3\lambda_7v_1v_2
\right].
\end{aligned}}
\]

## Derivation

At the vacuum,

\[
\Phi_1^\dagger\Phi_1=\frac{v_1^2}{2},
\quad
\Phi_2^\dagger\Phi_2=\frac{v_2^2}{2},
\quad
\Phi_1^\dagger\Phi_2=\frac{v_1v_2}{2}.
\]

The Hermitian conjugate doubles real `m12`, `lambda5`, `lambda6`, `lambda7` contributions where appropriate. In particular,

\[
V_0\supset\frac12\lambda_6v_1^3v_2
+\frac12\lambda_7v_1v_2^3.
\]

Differentiation directly explains the easily missed factors:

\[
\frac{\partial}{\partial v_1}(v_1^3v_2)=3v_1^2v_2,
\qquad
\frac{\partial}{\partial v_2}(v_1v_2^3)=3v_1v_2^2.
\]

These are the origin of the `3 lambda6` and `3 lambda7` terms; they are not source conventions.

The neutral CP-even tadpoles obey

\[
\left.\frac{\partial V}{\partial\rho_i}\right|_0
=\frac{\partial V_0}{\partial v_i},
\]

because `rho_i` appears through `v_i+rho_i` around the selected vacuum.

Only after minimization do we introduce the useful shorthand

\[
\boxed{M^2\equiv\frac{m_{12}^2}{s_\beta c_\beta}}.
\]

Thus

\[
\boxed{m_{22}^2\neq m_{12}^2\neq M^2}
\]

generically, despite all three having mass dimension two.

## Convention map

In beta notation,

\[
\begin{aligned}
m_{11}^2={}&M^2s^2-\frac{v^2}{2}
\left[\lambda_1c^2+\lambda_{345}s^2+3\lambda_6sc+\lambda_7\frac{s^3}{c}\right],\\
m_{22}^2={}&M^2c^2-\frac{v^2}{2}
\left[\lambda_2s^2+\lambda_{345}c^2+\lambda_6\frac{c^3}{s}+3\lambda_7sc\right].
\end{aligned}
\]

`M2` in project data means this derived `M^2`; it is not the generic-basis coefficient `m22_2`.

## What was checked against the source

The derived stationarity equations match DH05 after specializing its general phase-dependent expressions to the real CP-conserving branch. A symbolic derivative audit reproduces them. Active 2HDMC `set_param_gen` reproduces the derived `m22^2` relation term by term.

## What remains uncertain

Vanishing tadpoles prove stationarity, not global minimality. The Hessian must be examined next; global-vacuum questions remain distinct even after positive local masses are obtained.

## Next smallest validation

Compute charged, CP-odd and CP-even Hessians directly from the same potential before defining physical scalar states.

---

# Phase 4 — Hessians, Goldstones, masses and scalar-state signs

## What is established

After imposing the tadpoles, both the charged and CP-odd Hessians factorize as

\[
\boxed{
\mathcal M^2=D
\begin{pmatrix}
v_2/v_1&-1\\
-1&v_1/v_2
\end{pmatrix}}
\]

with a sector-dependent scalar `D`.

The vacuum vector is therefore an exact zero mode:

\[
\mathcal M^2\binom{v_1}{v_2}=0,
\]

while the orthogonal vector `(-v2,v1)` is physical. Consequently the beta rotation is derived from the vacuum geometry:

\[
\boxed{G^+=c_\beta\phi_1^+ +s_\beta\phi_2^+},
\qquad
\boxed{H^+=-s_\beta\phi_1^+ +c_\beta\phi_2^+},
\]

\[
\boxed{G^0=c_\beta\eta_1+s_\beta\eta_2},
\qquad
\boxed{A=-s_\beta\eta_1+c_\beta\eta_2}.
\]

The physical masses are

\[
\boxed{m_A^2=M^2-\frac{v^2}{2}
\left(2\lambda_5+\lambda_6\cot\beta+\lambda_7\tan\beta\right)},
\]

\[
\boxed{m_{H^\pm}^2=M^2-\frac{v^2}{2}
\left(\lambda_4+\lambda_5+\lambda_6\cot\beta+\lambda_7\tan\beta\right)},
\]

so

\[
\boxed{m_{H^\pm}^2-m_A^2=\frac{v^2}{2}(\lambda_5-\lambda_4)}.
\]

## Derivation: CP-even Hessian

Before tadpole elimination,

\[
\mathcal M_\rho^2=
\begin{pmatrix}\mathcal M_{11}^2&\mathcal M_{12}^2\\\mathcal M_{12}^2&\mathcal M_{22}^2\end{pmatrix}
\]

with

\[
\begin{aligned}
\mathcal M_{11}^2={}&m_{11}^2+\frac32\lambda_1v_1^2+\frac12\lambda_{345}v_2^2+3\lambda_6v_1v_2,\\
\mathcal M_{22}^2={}&m_{22}^2+\frac32\lambda_2v_2^2+\frac12\lambda_{345}v_1^2+3\lambda_7v_1v_2,\\
\mathcal M_{12}^2={}&-m_{12}^2+\lambda_{345}v_1v_2+\frac32\lambda_6v_1^2+\frac32\lambda_7v_2^2.
\end{aligned}
\]

After tadpoles it is especially useful to express the same matrix in terms of `m_A^2`:

\[
\boxed{
\begin{aligned}
\mathcal M_{11}^2={}&m_A^2s^2+v^2(\lambda_1c^2+\lambda_5s^2+2\lambda_6sc),\\
\mathcal M_{22}^2={}&m_A^2c^2+v^2(\lambda_2s^2+\lambda_5c^2+2\lambda_7sc),\\
\mathcal M_{12}^2={}&-m_A^2sc+v^2[(\lambda_3+\lambda_4)sc+\lambda_6c^2+\lambda_7s^2].
\end{aligned}}
\]

This is algebraically equivalent to the `M^2` representation and matches active 2HDMC after the independent derivation.

## CP-even diagonalization and frozen sign convention

Adopt the DH05 convention

\[
\boxed{
\begin{pmatrix}h\\\phi\end{pmatrix}
=\begin{pmatrix}-s_\alpha&c_\alpha\\c_\alpha&s_\alpha\end{pmatrix}
\begin{pmatrix}\rho_1\\\rho_2\end{pmatrix}}
\]

on the project branch, with `h=h_DH` and `phi=H_DH`.

Zeroing the off-diagonal element of the rotated Hessian gives

\[
\boxed{\tan2\alpha=\frac{2\mathcal M_{12}^2}{\mathcal M_{11}^2-\mathcal M_{22}^2}},
\]

with the quadrant fixed by the eigenvector branch, not the tangent alone.

The eigenvalues are

\[
\boxed{
m_{h,H}^2=\frac12\left[
\mathcal M_{11}^2+\mathcal M_{22}^2
\mp\sqrt{(\mathcal M_{11}^2-\mathcal M_{22}^2)^2+4(\mathcal M_{12}^2)^2}
\right]}.
\]

Define the CP-even vacuum-aligned and orthogonal directions

\[
\boxed{\rho_v=c_\beta\rho_1+s_\beta\rho_2},
\qquad
\boxed{\rho_\perp=-s_\beta\rho_1+c_\beta\rho_2}.
\]

Then

\[
\boxed{h=s_{\beta-\alpha}\rho_v+c_{\beta-\alpha}\rho_\perp},
\]

\[
\boxed{\phi=c_{\beta-\alpha}\rho_v-s_{\beta-\alpha}\rho_\perp}.
\]

On the project exact-alignment branch

\[
s_{\beta-\alpha}=1,
\]

so

\[
\boxed{h=\rho_v=c_\beta\rho_1+s_\beta\rho_2},
\]

\[
\boxed{\phi=-\rho_\perp=s_\beta\rho_1-c_\beta\rho_2}.
\]

This sign of `phi` is frozen **before** deriving Yukawa or scalar trilinear couplings. It is not chosen to obtain a desired `-cot beta` or charged-Higgs-loop sign.

## Convention map

BFLRS11's early simple display uses global negatives of the DH CP-even fields. That is a field-sign convention and cannot be mixed selectively with its coupling tables. The project consistently retains the DH signs above.

## What was checked against the source

DH05 gives the same charged/Goldstone beta rotations and CP-even alpha convention. GHOO18 uses the same beta geometry. The later general scalar-sector Hessian in BFLRS11 agrees with the independently derived entries. Active 2HDMC reproduces `m_A^2`, the charged–odd splitting and the CP-even matrix.

A source-internal caution is retained: an early restricted BFLRS11 pedagogical subsection displays charged/odd “mass terms” with factors that do not match the Hessian of its own printed potential; those early formulas are not used as normalization evidence.

## What remains uncertain

Positive physical squared masses are local quadratic conditions, not proof of the global electroweak minimum. At an exact CP-even degeneracy the mixing angle is not uniquely defined. The gauge interpretation of `rho_v` is deliberately left for Phase 6 rather than inferred here.

## Next smallest validation

Derive Type-I Yukawa modifiers from the gauge-invariant Yukawa Lagrangian using the already-fixed scalar signs.

---

# Phase 5 — Type-I Yukawa sector

## What is established

The Type-I assignment places all charged-fermion Yukawa couplings on `Phi2`. Write the sign convention explicitly:

\[
\boxed{
-\mathcal L_Y^{\rm I}
=\bar Q_LY_d\Phi_2d_R
+\bar Q_LY_u\widetilde\Phi_2u_R
+\bar L_LY_\ell\Phi_2\ell_R
+\text{h.c.}}
\]

where `tilde Phi2=i sigma2 Phi2*` for the up-type operator.

After fermion mass diagonalization,

\[
\boxed{m_f=\frac{y_fv_2}{\sqrt2}=\frac{y_fvs_\beta}{\sqrt2}},
\qquad f=u,d,\ell.
\]

The neutral CP-even interaction therefore begins as

\[
\boxed{
\mathcal L_Y^{\rm CP-even}
=-\sum_f\frac{m_f}{v_2}\rho_2\bar f f}.
\]

## Derivation

The Phase-4 rotation matrix is orthogonal and symmetric, so its inverse is itself:

\[
\rho_1=-s_\alpha h+c_\alpha\phi,
\qquad
\boxed{\rho_2=c_\alpha h+s_\alpha\phi}.
\]

Substitution gives

\[
\boxed{
\mathcal L_Y^{\rm CP-even}
=-\sum_f\frac{m_f}{v}
\left[
\frac{c_\alpha}{s_\beta}h
+\frac{s_\alpha}{s_\beta}\phi
\right]\bar f f}.
\]

Define modifiers through

\[
\mathcal L_Y^{\rm CP-even}
\equiv-\sum_f\frac{m_f}{v}
(\kappa_f^h h+\kappa_f^\phi\phi)\bar f f.
\]

Then

\[
\boxed{\kappa_f^h=\frac{c_\alpha}{s_\beta}},
\qquad
\boxed{\kappa_f^\phi=\frac{s_\alpha}{s_\beta}}.
\]

At exact alignment,

\[
\alpha=\beta-\frac\pi2,
\qquad
c_\alpha=s_\beta,
\qquad
s_\alpha=-c_\beta,
\]

so

\[
\boxed{\kappa_f^h=1},
\qquad
\boxed{\kappa_f^\phi=-\cot\beta},
\qquad f=u,d,\ell.
\]

The minus sign is therefore a consequence of the already-frozen field convention

\[
\phi=s_\beta\rho_1-c_\beta\rho_2,
\]

not an imported coupling-table sign.

If one instead defines

\[
\mathcal L_{\phi ff}=-g_{\phi ff}\phi\bar f f,
\]

then

\[
\boxed{g_{\phi ff}=-\frac{m_f}{v}\cot\beta}
\]

and the corresponding scalar Feynman rule is

\[
\boxed{+i\frac{m_f}{v}\cot\beta}.
\]

This explicit separation prevents later confusion between a dimensionless modifier, a Lagrangian coefficient and a vertex factor.

## Convention map

BFLRS11 Type I gives `c_alpha/s_beta` for the DH-sign light-state convention after translation and `s_alpha/s_beta` for the companion CP-even state. GHOO18 uses a Higgs-basis presentation and must be translated at the level of its explicit `-L_Y`, not by guessing the meaning of a table sign.

## What was checked against the source

Only after the derivation, BFLRS11 was used to confirm that all `u_R,d_R,e_R` couple to `Phi2` in Type I and that the functional dependence agrees. Active 2HDMC Type-I paths use `set_yukawas_type(1)` and scale its orthogonal-basis Yukawa matrices by `cot beta`, consistent with the derived structure after accounting for field-sign conventions.

## What remains uncertain

Nothing blocks the neutral Type-I modifier result. A global field redefinition `phi -> -phi` would flip every odd-`phi` coupling simultaneously; the project has already frozen one consistent sign convention.

## Next smallest validation

Derive the `WW/ZZ` modifiers from the scalar kinetic term instead of using a coupling table.

---

# Phase 6 — Gauge couplings and the physical meaning of alignment

## What is established

Start from

\[
\boxed{\mathcal L_{\rm kin}=\sum_{i=1}^2(D_\mu\Phi_i)^\dagger(D^\mu\Phi_i)}.
\]

Use modern hypercharge notation

\[
Q=T^3+Y,
\qquad
Y(\Phi_i)=\frac12,
\]

and

\[
\boxed{D_\mu=\partial_\mu+i g\frac{\sigma^a}{2}W_\mu^a+i g'\frac12B_\mu}.
\]

The gauge masses are

\[
\boxed{m_W^2=\frac{g^2v^2}{4}},
\qquad
\boxed{m_Z^2=\frac{(g^2+g'^2)v^2}{4}}.
\]

Most importantly, every tree-level neutral CP-even `SVV` coupling is controlled by one scalar combination:

\[
\boxed{v_1\rho_1+v_2\rho_2=v\rho_v}.
\]

Therefore

\[
\boxed{\kappa_V^h=\sin(\beta-\alpha)},
\qquad
\boxed{\kappa_V^\phi=\cos(\beta-\alpha)}.
\]

At exact alignment,

\[
\boxed{\kappa_V^h=1},
\qquad
\boxed{\kappa_V^\phi=0}.
\]

This is the first point where the Phase-4 geometrical statement `h=rho_v`, `phi=-rho_perp` obtains its physical gauge interpretation.

## Derivation: charged gauge part

Set temporarily

\[
\Phi_i\rightarrow\frac1{\sqrt2}\binom0{x_i},
\qquad x_i=v_i+\rho_i.
\]

With

\[
W^\pm=\frac{W^1\mp iW^2}{\sqrt2},
\]

one finds

\[
(T^1W^1+T^2W^2)\frac1{\sqrt2}\binom0{x_i}
=\binom{x_iW^+/2}{0}.
\]

Hence

\[
(D_\mu\Phi_i)^\dagger(D^\mu\Phi_i)
\supset\frac{g^2}{4}x_i^2W_\mu^+W^{-\mu}.
\]

Summing the doublets and expanding,

\[
\mathcal L\supset
\frac{g^2}{4}(v_1^2+v_2^2)W^+W^-
+\frac{g^2}{2}(v_1\rho_1+v_2\rho_2)W^+W^-+\cdots.
\]

Thus

\[
\boxed{\mathcal L_{WW,\rm linear}=\frac{2m_W^2}{v}\rho_vW_\mu^+W^{-\mu}}.
\]

## Derivation: neutral gauge part and photon cancellation

For the neutral lower component,

\[
T^3=-\frac12,
\qquad Y=+\frac12.
\]

Define

\[
A=s_WW^3+c_WB,
\qquad
Z=c_WW^3-s_WB,
\]

with

\[
s_W=\frac{g'}{\sqrt{g^2+g'^2}},
\qquad
c_W=\frac{g}{\sqrt{g^2+g'^2}},
\qquad
g_Z=\sqrt{g^2+g'^2}.
\]

The lower-component gauge factor is

\[
-i\frac g2W^3+i\frac{g'}2B
=-i\frac{g_Z}{2}Z.
\]

The photon cancels exactly because the neutral VEV has `Q=0`. Therefore

\[
(D_\mu\Phi_i)^\dagger(D^\mu\Phi_i)
\supset\frac{g_Z^2}{8}x_i^2Z_\mu Z^\mu.
\]

After summing and expanding,

\[
\boxed{\mathcal L_{ZZ,\rm linear}=\frac{m_Z^2}{v}\rho_vZ_\mu Z^\mu}.
\]

The factor-of-two difference between the displayed `WW` and `ZZ` Lagrangian coefficients is due to the identical real `Z` fields. The physical three-point rules are

\[
SW_\mu^+W_\nu^-:\quad i\kappa_V^S\frac{2m_W^2}{v}g_{\mu\nu},
\]

\[
SZ_\mu Z_\nu:\quad i\kappa_V^S\frac{2m_Z^2}{v}g_{\mu\nu}.
\]

## Projection onto mass eigenstates

Using

\[
\rho_1=-s_\alpha h+c_\alpha\phi,
\qquad
\rho_2=c_\alpha h+s_\alpha\phi,
\]

we obtain

\[
\begin{aligned}
\rho_v
&=c_\beta\rho_1+s_\beta\rho_2\\
&=(-c_\beta s_\alpha+s_\beta c_\alpha)h
 +(c_\beta c_\alpha+s_\beta s_\alpha)\phi\\
&=\boxed{s_{\beta-\alpha}h+c_{\beta-\alpha}\phi}.
\end{aligned}
\]

Therefore, defining

\[
\mathcal L_{SVV}\equiv\kappa_V^S
\left[
\frac{2m_W^2}{v}SW^+W^-+\frac{m_Z^2}{v}SZZ
\right],
\]

immediately gives the modifiers above.

The conceptual result is basis-geometric:

\[
\boxed{
\text{tree-level }SVV\text{ strength}
=\text{projection of }S\text{ onto the VEV direction}.}
\]

It depends only on the canonical kinetic terms and VEV geometry, not on the scalar potential.

## CP-odd and quartic checks

With imaginary fields retained,

\[
|v_i+\rho_i+i\eta_i|^2=(v_i+\rho_i)^2+\eta_i^2.
\]

No term linear in `eta_i` exists, so

\[
\boxed{AWW=AZZ=0}
\]

at tree level for the linear `AVV` vertices.

Orthogonality of the CP-even rotation gives

\[
\boxed{\rho_1^2+\rho_2^2=h^2+\phi^2}.
\]

Thus the diagonal CP-even quartic gauge interactions are angle-independent and the mixed `VVhphi` term cancels. This is an internal consistency check independent of the trilinear projection.

## Convention map

The project/DH state convention gives

\[
\kappa_V^h=s_{\beta-\alpha},
\qquad
\kappa_V^\phi=c_{\beta-\alpha}.
\]

BFLRS11 states the same physical modifiers but its early displayed scalar fields carry global sign differences relative to DH; those displays are not mixed with the project convention. GHOO18 defines alignment as a mass eigenstate parallel to the VEV and gives it the SM `VV` coupling.

Active 2HDMC uses

\[
q_{k1}=(s_{\beta-\alpha},c_{\beta-\alpha},0,i)
\]

for `(h,H,A,H+)`, and its `get_coupling_vvh` multiplies `Re(q_k1)` by the SM `ZZ` or `WW` vertex. This is a posterior implementation check of the derived geometry.

## What was checked against the source

DH05 gives `h_DH VV ∝ sin(beta-alpha)` and `H_DH VV ∝ cos(beta-alpha)`. BFLRS11 and GHOO18 agree after convention translation. 2HDMC reproduces `h:sba`, `H:cba`, `A:0` in its gauge-vertex implementation.

The symbolic audit `checks/phase6_gauge_check.py` verifies

\[
-c_\beta s_\alpha+s_\beta c_\alpha=\sin(\beta-\alpha),
\]

\[
c_\beta c_\alpha+s_\beta s_\alpha=\cos(\beta-\alpha),
\]

the exact-alignment values `(1,0)`, and the angle-independent norm `rho1^2+rho2^2=h^2+phi^2`.

## What remains uncertain

Phase 6 does **not** prove

\[
\text{exact alignment}\Longleftrightarrow Z_6=0.
\]

That statement belongs to the Higgs-basis potential/mass matrix and must be derived in Phase 7.

A new implementation caution is also frozen: 2HDMC's `get_param_higgs` exposes variables called `Lambda6,Lambda7`, while its scalar-trilinear routine later defines local `Z6=-l6`, `Z7=-l7`. No project statement about the sign of `Z6`, `Z7` or `phi H^+H^-` may use those names without an explicit Phase-7/9 convention audit.

Loop-induced `phi -> gamma gamma` and `phi -> Z gamma` remain possible; `kappa_V^phi=0` only states that the tree-level linear `phi WW` and `phi ZZ` vertices vanish in exact alignment.

## Next smallest validation

Construct the Higgs basis explicitly:

\[
\boxed{H_1=c_\beta\Phi_1+s_\beta\Phi_2},
\qquad
\boxed{H_2=-s_\beta\Phi_1+c_\beta\Phi_2}.
\]

Then independently derive:

1. `H1` carries the entire VEV and `H2` none;
2. all quadratic coefficients `Y1,Y2,Y3`;
3. all quartics `Z1...Z7`;
4. the stationarity equations in this basis;
5. the CP-even Higgs-basis mass matrix and the precise relation between its off-diagonal element and alignment.

Only after those steps may `Z7` and the charged-Higgs trilinear be interpreted.

---

# B. Validation status after Phase 6

| Claim | Status |
|---|---|
| field/VEV normalization | **VERIFIED** |
| complete generic CP-conserving potential | **VERIFIED** |
| vacuum potential and tadpoles | **VERIFIED** |
| `m22^2`, `m12^2`, `M^2` distinction | **VERIFIED** |
| charged, CP-odd and CP-even Hessians | **VERIFIED** |
| Goldstone and physical beta rotations | **VERIFIED** |
| DH/project CP-even state signs | **VERIFIED** |
| Type-I `kappa_f^h=c_alpha/s_beta` | **VERIFIED** |
| Type-I exact-alignment `kappa_f^phi=-cot beta` | **VERIFIED** |
| `kappa_V^h=sin(beta-alpha)` | **VERIFIED** |
| `kappa_V^phi=cos(beta-alpha)` | **VERIFIED** |
| exact alignment: `kappa_V^h=1`, `kappa_V^phi=0` | **VERIFIED** |
| tree-level `AVV=0` | **VERIFIED** |
| Higgs-basis `Y3=-Z6 v^2/2` | **PARTIALLY VERIFIED: source only; derivation pending** |
| exact alignment `iff Z6=0` | **PARTIALLY VERIFIED; Higgs-basis derivation pending** |
| exact `phi H+H-` trilinear | **NOT VERIFIED** |
| exact generic-to-Higgs-basis `Z7` map | **PARTIALLY VERIFIED: source only** |
| large-`tan beta`, `lambda7=0` limit of `Z7` | **NOT VERIFIED** |
| `g_(phi H+H-) ~ -v lambda6` | **NOT VERIFIED** |
| re-expression in `X=lambda6 tan beta` | **NOT VERIFIED** |

The distinction between the verified gauge alignment statements and the pending `Z6` statement is intentional. Phase 6 establishes alignment operationally through VEV projection and gauge interactions; Phase 7 must establish its Higgs-basis parameter criterion.

## Phase gates

| Gate | Status | Evidence |
|---|---|---|
| `PHASE_0_PASS` | **PASS** | source/convention inventory frozen |
| `PHASE_1_PASS` | **PASS** | fields, charges, VEVs and beta established |
| `PHASE_2_PASS` | **PASS** | complete generic potential reconstructed |
| `PHASE_3_PASS` | **PASS** | `V0`, tadpoles and soft-coordinate distinction derived |
| `PHASE_4_PASS` | **PASS** | all scalar Hessians, Goldstones, masses and state signs derived |
| `PHASE_5_PASS` | **PASS** | Type-I modifiers derived from `-L_Y` |
| `PHASE_6_PASS` | **PASS** | `WW/ZZ` masses and scalar modifiers derived from kinetic terms; exact-alignment gauge meaning established |
| `PHASE_7_PASS` | NOT RUN | Higgs-basis potential must be independently reconstructed |
| `PHASE_9_PASS` | NOT RUN | exact charged-Higgs trilinear not yet derived |

---

# C. Open issues that remain intentionally open

1. **Global vacuum:** positive physical masses test local quadratic curvature; they do not alone prove the selected neutral stationary point is the global minimum.
2. **Exact CP-even degeneracy:** at exact degeneracy a unique mixing angle is not defined; state-identification statements must state their non-degeneracy assumption.
3. **BFLRS11 source-internal sign/factor cautions:** project signs come from its own derivation, not selective source snippets.
4. **Residual Higgs-basis sign:** `H2 -> -H2` changes `Y3,Z6,Z7` and odd-`H2` couplings together. Phase 7 must freeze the project Higgs-basis sign by the explicit rotation above.
5. **2HDMC Higgs-basis sign layer:** `get_param_higgs`/`get_coupling_hhh` naming and sign translations require an explicit audit before trilinear comparison.
6. **GHOO18 approximate-alignment prose typo:** a prose sentence suggesting large `|Z6|` is not used; the explicit Higgs-basis mass matrix will decide the condition.

---

# D. Maintenance rule for all later phases

A phase is not considered closed merely because an auxiliary file exists. Before its gate becomes `PASS`, this canonical document must be enriched with:

1. the assumptions imported from earlier phases;
2. the derivation in sufficient detail to reproduce signs/factors;
3. the explicit convention map;
4. source checks performed after the derivation;
5. implementation checks, if relevant, performed after the analytic result is frozen;
6. unresolved questions and the smallest next validation.

This prevents the appendix from becoming a chain of disconnected formulas and keeps the complete reasoning auditable by an independent reviewer.
