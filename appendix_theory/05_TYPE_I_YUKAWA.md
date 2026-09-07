# Phase 5 — Type-I Yukawa sector from the fixed field convention

This phase derives the neutral CP-even Yukawa modifiers from the gauge-invariant Yukawa Lagrangian. No coupling table is used as input.

The field/sign convention imported from Phase 4 is

\[
\Phi_i=\begin{pmatrix}\phi_i^+\\(v_i+\rho_i+i\eta_i)/\sqrt2\end{pmatrix},
\qquad
v_1=v c_\beta,\quad v_2=v s_\beta,
\]

and

\[
\begin{pmatrix}h\\\phi\end{pmatrix}
=
\begin{pmatrix}-s_\alpha&c_\alpha\\c_\alpha&s_\alpha\end{pmatrix}
\begin{pmatrix}\rho_1\\\rho_2\end{pmatrix}.
\]

On the project alignment branch,

\[
\sin(\beta-\alpha)=1,
\qquad
\alpha=\beta-\frac\pi2,
\]

so that

\[
h=c_\beta\rho_1+s_\beta\rho_2,
\qquad
\phi=s_\beta\rho_1-c_\beta\rho_2.
\]

---

## What is established

[DERIVED] In Type I, all charged-fermion masses arise from the neutral VEV of `Phi2`. After diagonalizing the fermion mass matrices,

\[
\boxed{m_f=\frac{y_f v_2}{\sqrt2}=\frac{y_f v s_\beta}{\sqrt2}}
\]

for each mass eigenstate `f=u,d,l`.

[DERIVED] The CP-even neutral interaction is therefore

\[
\boxed{\mathcal L_Y^{\rm CP-even}
=-\sum_f\frac{m_f}{v_2}\rho_2\,\bar f f}.
\]

[DERIVED] The Phase-4 CP-even rotation is orthogonal and symmetric, so its inverse is itself. Hence

\[
\boxed{\rho_2=c_\alpha h+s_\alpha\phi}.
\]

Substituting gives

\[
\boxed{
\mathcal L_Y^{\rm CP-even}
=-\sum_f\frac{m_f}{v}
\left(
\frac{c_\alpha}{s_\beta}h+
\frac{s_\alpha}{s_\beta}\phi
\right)\bar f f
}.
\]

Defining the modifiers by

\[
\mathcal L_Y^{\rm CP-even}
\equiv
-\sum_f\frac{m_f}{v}
\left(\kappa_f^h h+\kappa_f^\phi\phi\right)\bar f f,
\]

we obtain

\[
\boxed{\kappa_f^h=\frac{c_\alpha}{s_\beta}},
\qquad
\boxed{\kappa_f^\phi=\frac{s_\alpha}{s_\beta}},
\qquad f=u,d,\ell.
\]

At exact alignment,

\[
c_\alpha=c_{\beta-\pi/2}=s_\beta,
\qquad
s_\alpha=s_{\beta-\pi/2}=-c_\beta,
\]

therefore

\[
\boxed{\kappa_f^h=1},
\qquad
\boxed{\kappa_f^\phi=-\cot\beta},
\qquad f=u,d,\ell.
\]

The minus sign is not imported from a table. It follows from the Phase-4 state convention `phi=s_beta rho1-c_beta rho2`, equivalently from `s_alpha=-c_beta` on the selected alignment branch.

---

## Derivation

### 1. Gauge-invariant Yukawa Lagrangian before electroweak symmetry breaking

For one Higgs doublet with hypercharge appropriate to the down-type and charged-lepton Yukawa operators, the up-type operator uses the conjugate doublet

\[
\widetilde\Phi_2\equiv i\sigma_2\Phi_2^*.
\]

The Type-I restriction is that `Phi1` does not appear in any charged-fermion Yukawa operator. In flavor-matrix notation,

\[
\boxed{
-\mathcal L_Y^{\rm I}
=
\overline Q_L Y_d\Phi_2 d_R
+
\overline Q_L Y_u\widetilde\Phi_2 u_R
+
\overline L_L Y_\ell\Phi_2\ell_R
+\mathrm{h.c.}
}.
\]

This equation fixes the global sign convention used below. Statements such as “the Yukawa coupling is plus/minus ...” are meaningless unless this Lagrangian sign and the scalar field sign are both stated.

### 2. Neutral component and the origin of the fermion masses

Using

\[
\Phi_2^0=\frac{v_2+\rho_2+i\eta_2}{\sqrt2},
\qquad
\widetilde\Phi_2^0=\frac{v_2+\rho_2-i\eta_2}{\sqrt2},
\]

we see that the real CP-even fluctuation `rho2` enters with the same sign in the up-, down- and charged-lepton sectors. The sign difference in the imaginary part is relevant for pseudoscalar couplings, but not for the CP-even result derived here.

Before flavor diagonalization, the mass matrices are

\[
M_d=\frac{v_2}{\sqrt2}Y_d,
\qquad
M_u=\frac{v_2}{\sqrt2}Y_u,
\qquad
M_\ell=\frac{v_2}{\sqrt2}Y_\ell.
\]

Bi-unitary transformations diagonalize each matrix. Rephasing the mass eigenfields makes the physical masses positive. For a mass eigenstate `f`,

\[
y_f=\frac{\sqrt2 m_f}{v_2}.
\]

The neutral CP-even term then becomes

\[
\mathcal L_Y^{\rm CP-even}
=-\sum_f\frac{y_f}{\sqrt2}\rho_2\bar f f
=-\sum_f\frac{m_f}{v_2}\rho_2\bar f f.
\]

Using `v2=v s_beta`,

\[
\boxed{
\mathcal L_Y^{\rm CP-even}
=-\sum_f\frac{m_f}{v s_\beta}\rho_2\bar f f
}.
\]

No mixing angle `alpha` has entered yet.

### 3. Invert the already-derived CP-even scalar rotation

Phase 4 fixed

\[
R_\alpha=
\begin{pmatrix}-s_\alpha&c_\alpha\\c_\alpha&s_\alpha\end{pmatrix},
\qquad
\begin{pmatrix}h\\\phi\end{pmatrix}
=R_\alpha
\begin{pmatrix}\rho_1\\\rho_2\end{pmatrix}.
\]

Direct multiplication gives

\[
R_\alpha^T R_\alpha=I,
\qquad
R_\alpha^T=R_\alpha.
\]

Therefore

\[
\begin{pmatrix}\rho_1\\\rho_2\end{pmatrix}
=R_\alpha
\begin{pmatrix}h\\\phi\end{pmatrix}
\]

and explicitly

\[
\rho_1=-s_\alpha h+c_\alpha\phi,
\qquad
\boxed{\rho_2=c_\alpha h+s_\alpha\phi}.
\]

This inversion is the sign-sensitive step. If `phi` were globally redefined by `phi -> -phi`, every odd-`phi` interaction would reverse sign. We do not make such a redefinition here.

### 4. Substitute rho2 into the Yukawa interaction

Insert the inverse rotation into

\[
-\frac{m_f}{v s_\beta}\rho_2\bar f f.
\]

Then

\[
\mathcal L_Y^{\rm CP-even}
=-\frac{m_f}{v s_\beta}
(c_\alpha h+s_\alpha\phi)\bar f f.
\]

Comparing with the declared modifier convention

\[
\mathcal L_Y^{\rm CP-even}
=-\frac{m_f}{v}
(\kappa_f^h h+\kappa_f^\phi\phi)\bar f f,
\]

gives

\[
\boxed{\kappa_f^h=c_\alpha/s_\beta},
\qquad
\boxed{\kappa_f^\phi=s_\alpha/s_\beta}.
\]

Because all three charged-fermion sectors use the same `Phi2`, these expressions are universal in Type I.

### 5. Exact alignment without assuming the answer

The selected project branch is

\[
\sin(\beta-\alpha)=1.
\]

For the continuous branch used by 2HDMC and fixed in Phase 4,

\[
\beta-\alpha=\frac\pi2,
\qquad
\alpha=\beta-\frac\pi2.
\]

Hence

\[
c_\alpha=\sin\beta,
\qquad
s_\alpha=-\cos\beta.
\]

Therefore

\[
\kappa_f^h=\frac{s_\beta}{s_\beta}=1,
\]

and

\[
\boxed{
\kappa_f^\phi
=\frac{-c_\beta}{s_\beta}
=-\cot\beta
}.
\]

The same result can be seen directly from the exact-alignment field identity

\[
\phi=s_\beta\rho_1-c_\beta\rho_2.
\]

Since Type-I fermions couple only to `rho2`, the coefficient of `phi` inside `rho2` is `-c_beta`; division by the mass-generating VEV fraction `s_beta` gives `-c_beta/s_beta`.

### 6. What this sign means physically

In the convention declared above,

\[
\mathcal L_{\phi ff}
=-\frac{m_f}{v}\kappa_f^\phi\phi\bar f f.
\]

Thus at exact alignment

\[
\mathcal L_{\phi ff}
=+\frac{m_f}{v}\cot\beta\,\phi\bar f f.
\]

With the standard convention that a term

\[
\mathcal L_{\rm int}=-g_{\phi ff}\,\phi\bar f f
\]

gives the vertex `-i g_{phi ff}`, we have

\[
g_{\phi ff}=\frac{m_f}{v}\kappa_f^\phi
=-\frac{m_f}{v}\cot\beta,
\]

and therefore the Feynman rule is

\[
\boxed{+i\frac{m_f}{v}\cot\beta}.
\]

To avoid ambiguity downstream, the project should quote `kappa_f^phi` and/or the explicit `L_int`, not an isolated signed “coupling” symbol.

---

## Convention map

| Object | Project convention | Branco et al. | GHOO18 | Status |
|---|---|---|---|---|
| Type-I basis restriction | only `Phi2` couples to `u_R,d_R,e_R` | Table `tab:3_models`: all three use `Phi2` | `eta_1^{u,d,l}=0` in Appendix `Yuk_Type_I` | `[SOURCE][PROJECT-DEFINITION]` |
| Fermion mass | `m_f=y_f v2/sqrt2` | implicit in Type-I table and `m_f/v` normalization | encoded through `kappa/rho` matrices | `[DERIVED]` |
| CP-even modifier of `h` | `c_alpha/s_beta` | Table `tab:3_couplings`: same | general `R_{j2}/s_beta` structure after state translation | `[DERIVED][SOURCE-CHECKED]` |
| CP-even modifier of `phi` | `s_alpha/s_beta` | Type-I `xi_H=s_alpha/s_beta` when source `H` is identified with the DH-sign state | same `R_{j2}/s_beta` structure, but heavy-state sign map must be made explicit | `[DERIVED][TRANSLATED]` |
| exact alignment | `kappa_h=1`, `kappa_phi=-cot beta` | follows from its table if the DH-sign state convention is used | magnitude `cot beta`; absolute heavy-state sign depends on state convention | `[DERIVED][SOURCE-CHECKED]` |

### Important Branco sign caution

The early scalar-field display in BFLRS11 uses CP-even fields that are global-sign reversals of the DH05 display, whereas its later Yukawa table gives the standard `c_alpha/s_beta` and `s_alpha/s_beta` modifiers. Therefore the table is used only as a post-derivation check after the project state convention is fixed. It is not used to define the project scalar signs.

### Important GHOO sign caution

GHOO18 uses its own neutral-state rotation matrix `R` and writes Yukawa couplings as coefficients of its `H_alpha` fields. A global sign change of a heavy scalar changes every odd-heavy-scalar interaction. Therefore its `cot beta` dependence is a strong source check, while the absolute sign must be translated together with the scalar-state convention. No GHOO sign is used to override the Phase-4 project definition.

---

## What was checked against the source

1. BFLRS11 Table `tab:3_models` explicitly assigns `Phi2` to `u_R`, `d_R`, and `e_R` in Type I.
2. BFLRS11 Eq. `Eq:Yukawa` defines its neutral Yukawa normalization with an overall `-m_f/v` in the Lagrangian.
3. BFLRS11 Table `tab:3_couplings` gives, for Type I, `xi_h^u=xi_h^d=xi_h^l=c_alpha/s_beta` and `xi_H^u=xi_H^d=xi_H^l=s_alpha/s_beta`.
4. GHOO18 Appendix `Yuk_Type_I` states `eta_1^{u,0}=eta_1^{d,0}=eta_1^{l,0}=0`, i.e. the Type-I basis restriction.
5. GHOO18 gives neutral Type-I couplings proportional to `R_{j2}/s_beta` for the CP-even component, confirming that the fermion coupling is controlled by the `Phi2` projection.
6. The project implementation explicitly installs `set_yukawas_type(1)` in the active evaluators; 2HDMC documents `set_yukawas_type` as using the hep-ph/0504050 convention. This is an implementation consistency check, not part of the analytic derivation.

### Algebra audit

The Phase-4 rotation matrix satisfies `R_alpha^T R_alpha=I` and `R_alpha^T=R_alpha`. Symbolic substitution of `alpha=beta-pi/2` gives exactly

\[
\frac{c_\alpha}{s_\beta}=1,
\qquad
\frac{s_\alpha}{s_\beta}=-\cot\beta.
\]

---

## What remains uncertain

1. The gauge coupling `kappa_V^phi=0` is not inferred from the Yukawa result; it remains for Phase 6 and must be derived from the kinetic terms.
2. The exact relation between alignment and Higgs-basis `Z6=0` remains for the Higgs-basis derivation.
3. The absolute sign of a GHOO heavy-state Yukawa coefficient is not a project convention until the GHOO heavy state is explicitly mapped to the project `phi` including any global field sign.
4. Pseudoscalar and charged-Higgs Yukawa vertices can be derived from the same Type-I Lagrangian, but they are not needed to validate `kappa_f^phi` and are intentionally not used as a shortcut here.

---

## Next smallest validation

Proceed to the gauge sector from

\[
\sum_{i=1}^2(D_\mu\Phi_i)^\dagger(D^\mu\Phi_i).
\]

Expand only the terms linear in `rho_i` and quadratic in `W/Z`. Show that the gauge-boson mass term and the CP-even scalar coupling both select the vacuum direction `c_beta rho1+s_beta rho2`. Then translate that direction using the Phase-4 rotation. This will independently determine `kappa_V^h` and `kappa_V^phi` without importing a coupling table.
