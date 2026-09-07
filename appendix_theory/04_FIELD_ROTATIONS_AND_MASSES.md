# 04 — Field rotations and scalar masses

Scope: Phase 4 of issue #81.  
Input convention: `02_SCALAR_POTENTIAL.md` and the tadpole relations of `03_VACUUM_AND_MINIMIZATION.md`.  
Primary analytic anchor for post-derivation comparison: Davidson–Haber (DH05), especially Eqs. `gpm`, `goldn`, `scalareigenstates`, `massmhh`, `chhiggsmass`.  
Independent cross-source check: Branco et al. (BFLRS11), general scalar-sector discussion around Eqs. `2_eq:mpseu`, `2_eq:mneut`; GHOO18 charged/neutral rotation in section 2.

The mass matrices below are obtained from the same potential already frozen in Phase 2. The source formulas are checked only after the Hessians have been constructed.

Throughout this file,

\[
\lambda_{345}\equiv\lambda_3+\lambda_4+\lambda_5,
\qquad
s_\beta\equiv s,\quad c_\beta\equiv c,
\qquad
v_1=vc,\quad v_2=vs.
\]

The soft coordinate is introduced only as the derived shorthand

\[
M^2\equiv\frac{m_{12}^2}{s c}.
\]

---

## What is established

[DERIVED] The charged and CP-odd Hessians become rank-one matrices after the neutral tadpole equations are imposed. Both possess the null eigenvector proportional to `(v1,v2)`, so the rotation angle `beta` is forced by the vacuum direction rather than postulated.

[DERIVED] The physical charged state and CP-odd state are

\[
\boxed{\begin{aligned}
G^+ &= c\,\phi_1^+ + s\,\phi_2^+,\\
H^+ &= -s\,\phi_1^+ + c\,\phi_2^+,
\end{aligned}}
\qquad
\boxed{\begin{aligned}
G^0 &= c\,\eta_1+s\,\eta_2,\\
A   &= -s\,\eta_1+c\,\eta_2.
\end{aligned}}
\]

with

\[
\boxed{m_{H^\pm}^2=M^2-\frac{v^2}{2}
\left(\lambda_4+\lambda_5+\lambda_6\cot\beta+\lambda_7\tan\beta\right)}
\]

and

\[
\boxed{m_A^2=M^2-\frac{v^2}{2}
\left(2\lambda_5+\lambda_6\cot\beta+\lambda_7\tan\beta\right)}.
\]

Consequently,

\[
\boxed{m_{H^\pm}^2-m_A^2=\frac{v^2}{2}(\lambda_5-\lambda_4)}.
\]

[DERIVED] The CP-even real fields `(rho1,rho2)` form a symmetric `2 x 2` Hessian. After tadpole elimination it is diagonalized by an angle `alpha`. In the DH05 sign convention,

\[
\boxed{\begin{pmatrix}h_{\rm DH}\\H_{\rm DH}\end{pmatrix}
=\begin{pmatrix}-\sin\alpha & \cos\alpha\\\cos\alpha & \sin\alpha\end{pmatrix}
\begin{pmatrix}\rho_1\\\rho_2\end{pmatrix}}.
\]

For the project branch continuously connected to `sin(beta-alpha)=1`, we freeze the state-sign convention

\[
\boxed{h\equiv h_{\rm DH},\qquad \phi\equiv H_{\rm DH}}.
\]

This map is based on the mixing direction, not numerical mass ordering. The statement that `h` has the SM gauge coupling will be independently derived from the kinetic terms in Phase 6.

---

## Derivation

### 1. Quadratic terms and what is meant by a mass matrix

For the real neutral fluctuations,

\[
V^{(2)}_{\rm neutral}=\frac12\,\eta^T\mathcal M_A^2\eta
+\frac12\,\rho^T\mathcal M_\rho^2\rho,
\]

with

\[
\eta=(\eta_1,\eta_2)^T,\qquad \rho=(\rho_1,\rho_2)^T.
\]

Thus

\[
(\mathcal M_A^2)_{ij}=\left.\frac{\partial^2 V}{\partial\eta_i\partial\eta_j}\right|_0,
\qquad
(\mathcal M_\rho^2)_{ij}=\left.\frac{\partial^2 V}{\partial\rho_i\partial\rho_j}\right|_0.
\]

For charged complex fields,

\[
V^{(2)}_\pm=(\phi_1^-,\phi_2^-)\mathcal M_\pm^2
\begin{pmatrix}\phi_1^+\\\phi_2^+\end{pmatrix},
\]

so

\[
(\mathcal M_\pm^2)_{ij}=\left.\frac{\partial^2 V}{\partial\phi_i^-\partial\phi_j^+}\right|_0.
\]

These matrices are coefficients of `V`. The mass terms in the Lagrangian carry the opposite sign because `L_int = -V_int`.

### 2. Charged sector before minimization

Direct expansion of the frozen potential gives

\[
\mathcal M_\pm^2=\begin{pmatrix}X_\pm & Y_\pm\\Y_\pm & Z_\pm\end{pmatrix},
\]

with

\[
\begin{aligned}
X_\pm={}&m_{11}^2+\frac12\lambda_1v_1^2+\frac12\lambda_3v_2^2+\lambda_6v_1v_2,\\
Z_\pm={}&m_{22}^2+\frac12\lambda_2v_2^2+\frac12\lambda_3v_1^2+\lambda_7v_1v_2,\\
Y_\pm={}&-m_{12}^2+\frac12(\lambda_4+\lambda_5)v_1v_2+\frac12\lambda_6v_1^2+\frac12\lambda_7v_2^2.
\end{aligned}
\]

No zero mode is manifest yet because `m11^2,m22^2` have not been eliminated with the vacuum equations.

Insert the Phase-3 tadpoles and define

\[
D_\pm\equiv m_{12}^2-\frac12\left[(\lambda_4+\lambda_5)v_1v_2+\lambda_6v_1^2+\lambda_7v_2^2\right].
\]

Then

\[
\boxed{\mathcal M_\pm^2=D_\pm
\begin{pmatrix}\dfrac{v_2}{v_1} & -1\\[2mm]-1 & \dfrac{v_1}{v_2}\end{pmatrix}}.
\]

Act on the vacuum direction explicitly:

\[
\mathcal M_\pm^2\begin{pmatrix}v_1\\v_2\end{pmatrix}
=D_\pm\begin{pmatrix}v_2-v_2\\-v_1+v_1\end{pmatrix}=0.
\]

Hence

\[
G^+=\frac{v_1\phi_1^++v_2\phi_2^+}{v}=c\phi_1^++s\phi_2^+.
\]

The normalized orthogonal direction is fixed up to an overall sign. We choose the DH05/GHOO18 sign,

\[
H^+=-s\phi_1^++c\phi_2^+.
\]

The nonzero eigenvalue is

\[
\begin{aligned}
m_{H^\pm}^2
&=D_\pm\left(\frac{v_2}{v_1}+\frac{v_1}{v_2}\right)\\
&=D_\pm\frac{v^2}{v_1v_2},
\end{aligned}
\]

therefore

\[
\boxed{m_{H^\pm}^2=M^2-\frac{v^2}{2}(\lambda_4+\lambda_5+\lambda_6\cot\beta+\lambda_7\tan\beta)}.
\]

### 3. CP-odd sector before minimization

The Hessian in `(eta1,eta2)` is

\[
\mathcal M_A^2=\begin{pmatrix}X_A&Y_A\\Y_A&Z_A\end{pmatrix},
\]

where

\[
\begin{aligned}
X_A={}&m_{11}^2+\frac12\lambda_1v_1^2+\frac12(\lambda_3+\lambda_4-\lambda_5)v_2^2+\lambda_6v_1v_2,\\
Z_A={}&m_{22}^2+\frac12\lambda_2v_2^2+\frac12(\lambda_3+\lambda_4-\lambda_5)v_1^2+\lambda_7v_1v_2,\\
Y_A={}&-m_{12}^2+\lambda_5v_1v_2+\frac12\lambda_6v_1^2+\frac12\lambda_7v_2^2.
\end{aligned}
\]

After the same tadpole substitution, define

\[
D_A\equiv m_{12}^2-\lambda_5v_1v_2-\frac12(\lambda_6v_1^2+\lambda_7v_2^2).
\]

Then

\[
\boxed{\mathcal M_A^2=D_A
\begin{pmatrix}\dfrac{v_2}{v_1} & -1\\[2mm]-1 & \dfrac{v_1}{v_2}\end{pmatrix}}.
\]

Therefore

\[
G^0=c\eta_1+s\eta_2,
\qquad
\boxed{A=-s\eta_1+c\eta_2}.
\]

The physical eigenvalue is

\[
\boxed{m_A^2=M^2-\frac{v^2}{2}(2\lambda_5+\lambda_6\cot\beta+\lambda_7\tan\beta)}.
\]

Subtracting the charged result gives

\[
\boxed{m_{H^\pm}^2=m_A^2+\frac{v^2}{2}(\lambda_5-\lambda_4)}.
\]

The `lambda6` and `lambda7` dependence cancels in this splitting.

### 4. Why beta diagonalizes both sectors

Both post-tadpole matrices have the form

\[
K\begin{pmatrix}\tan\beta&-1\\-1&\cot\beta\end{pmatrix}.
\]

Define

\[
R_\beta=\begin{pmatrix}c&s\\-s&c\end{pmatrix}.
\]

Then

\[
\boxed{R_\beta\mathcal M^2R_\beta^T=\operatorname{diag}(0,m_{\rm phys}^2)}
\]

for both charged and CP-odd sectors. Thus `beta` follows because the Goldstone direction is parallel to the vacuum vector.

### 5. CP-even sector before tadpole elimination

The direct Hessian is

\[
\boxed{\mathcal M_\rho^2=\begin{pmatrix}\mathcal M_{11}^2&\mathcal M_{12}^2\\\mathcal M_{12}^2&\mathcal M_{22}^2\end{pmatrix}}
\]

with

\[
\begin{aligned}
\mathcal M_{11}^2={}&m_{11}^2+\frac32\lambda_1v_1^2+\frac12\lambda_{345}v_2^2+3\lambda_6v_1v_2,\\
\mathcal M_{22}^2={}&m_{22}^2+\frac32\lambda_2v_2^2+\frac12\lambda_{345}v_1^2+3\lambda_7v_1v_2,\\
\mathcal M_{12}^2={}&-m_{12}^2+\lambda_{345}v_1v_2+\frac32\lambda_6v_1^2+\frac32\lambda_7v_2^2.
\end{aligned}
\]

The later general scalar-sector section of BFLRS11 gives these same three entries before tadpole simplification.

After inserting the tadpoles and rewriting with `M^2`,

\[
\boxed{\begin{aligned}
\mathcal M_{11}^2={}&M^2s^2+v^2\left[\lambda_1c^2+\frac32\lambda_6sc-\frac12\lambda_7\frac{s^3}{c}\right],\\
\mathcal M_{22}^2={}&M^2c^2+v^2\left[\lambda_2s^2-\frac12\lambda_6\frac{c^3}{s}+\frac32\lambda_7sc\right],\\
\mathcal M_{12}^2={}&-M^2sc+v^2\left[\lambda_{345}sc+\frac32\lambda_6c^2+\frac32\lambda_7s^2\right].
\end{aligned}}
\]

An equivalent representation useful for matching 2HDMC follows by eliminating `M^2` in favor of `m_A^2`:

\[
\boxed{\begin{aligned}
\mathcal M_{11}^2={}&m_A^2s^2+v^2(\lambda_1c^2+\lambda_5s^2+2\lambda_6sc),\\
\mathcal M_{22}^2={}&m_A^2c^2+v^2(\lambda_2s^2+\lambda_5c^2+2\lambda_7sc),\\
\mathcal M_{12}^2={}&-m_A^2sc+v^2[(\lambda_3+\lambda_4)sc+\lambda_6c^2+\lambda_7s^2].
\end{aligned}}
\]

The equivalence is algebraic. For example, inserting

\[
m_A^2=M^2-\frac{v^2}{2}(2\lambda_5+\lambda_6 c/s+\lambda_7 s/c)
\]

into the first line cancels `lambda5 s^2` and converts `2 lambda6 sc` into `3 lambda6 sc/2`, reproducing the first form.

### 6. CP-even diagonalization and alpha

Adopt the DH05 field-sign convention

\[
\begin{pmatrix}h_{\rm DH}\\H_{\rm DH}\end{pmatrix}
=R_\alpha\begin{pmatrix}\rho_1\\\rho_2\end{pmatrix},
\qquad
R_\alpha=\begin{pmatrix}-s_\alpha&c_\alpha\\c_\alpha&s_\alpha\end{pmatrix}.
\]

The off-diagonal entry of `R_alpha M_rho^2 R_alpha^T` is

\[
\mathcal M_{12}^2\cos2\alpha+\frac12(\mathcal M_{22}^2-\mathcal M_{11}^2)\sin2\alpha.
\]

Thus

\[
\boxed{\tan2\alpha=\frac{2\mathcal M_{12}^2}{\mathcal M_{11}^2-\mathcal M_{22}^2}}
\]

with the quadrant of `alpha` fixed by the eigenvector branch, not by this tangent alone.

The eigenvalues are

\[
\boxed{m_{h,H}^2=\frac12\left[\mathcal M_{11}^2+\mathcal M_{22}^2
\mp\sqrt{(\mathcal M_{11}^2-\mathcal M_{22}^2)^2+4(\mathcal M_{12}^2)^2}\right]}.
\]

DH05 conventionally calls the lower eigenvalue `h` and the upper `H`. This ordering convention is not what defines the project state identity below.

### 7. Vacuum-aligned and orthogonal CP-even directions

Define

\[
\rho_v=c\rho_1+s\rho_2,
\qquad
\rho_\perp=-s\rho_1+c\rho_2.
\]

Substituting this beta rotation into the DH states gives

\[
\boxed{h_{\rm DH}=\sin(\beta-\alpha)\rho_v+\cos(\beta-\alpha)\rho_\perp},
\]

\[
\boxed{H_{\rm DH}=\cos(\beta-\alpha)\rho_v-\sin(\beta-\alpha)\rho_\perp}.
\]

On the project campaign branch

\[
\sin(\beta-\alpha)=1,
\]

this becomes

\[
\boxed{h_{\rm DH}=\rho_v,\qquad H_{\rm DH}=-\rho_\perp}.
\]

The project therefore freezes

\[
\boxed{h\equiv h_{\rm DH},\qquad \phi\equiv H_{\rm DH}},
\]

so in exact alignment

\[
\boxed{h=c\rho_1+s\rho_2},
\qquad
\boxed{\phi=s\rho_1-c\rho_2=-\rho_\perp}.
\]

This sign is inherited from the DH05 CP-even convention selected in Phase 0. It is not selected to obtain a desired Yukawa or scalar-trilinear sign.

---

## Convention map

| Object | Project / DH-reference convention | BFLRS11 | GHOO18 | Status |
|---|---|---|---|---|
| Charged Goldstone | `G+=c beta phi1+ + s beta phi2+` | same beta direction | same explicit rotation | `[DERIVED][SOURCE]` |
| Physical charged state | `H+=-s beta phi1+ + c beta phi2+` | same geometry | same explicit sign | `[DERIVED][SOURCE]` |
| Neutral Goldstone | `G0=c beta eta1+s beta eta2` | same geometry | `G0=c beta chi1+s beta chi2` | `[DERIVED][SOURCE]` |
| CP-odd state | `A=-s beta eta1+c beta eta2` | pseudoscalar direction | called `eta3` before CP specialization | `[DERIVED][TRANSLATED]` |
| CP-even source states | `h_DH=-s alpha rho1+c alpha rho2`; `H_DH=c alpha rho1+s alpha rho2` | early display uses their global negatives | general `H_i=R_ij eta_j` | `[SOURCE][TRANSLATED]` |
| Project CP-even names | `h=h_DH`, `phi=H_DH` on branch connected to `s_(beta-alpha)=1` | do not import displayed signs directly | precise `H_i` map deferred | `[PROJECT-DEFINITION][TRANSLATED]` |
| Exact-alignment directions | `h=rho_v`, `phi=-rho_perp` | compatible after sign translation | vacuum-aligned state explicit in AL | `[DERIVED]`; gauge interpretation pending Phase 6 |
| `mH+^2` | expression above | later general section agrees | Higgs-basis `Y2,Z3` form | `[DERIVED][SOURCE-CHECKED]` |
| `mA^2` | expression above | later general section agrees | neutral Higgs-basis matrix | `[DERIVED][SOURCE-CHECKED]` |

### Internal BFLRS11 caution

The early pedagogical restricted subsection around source lines 317–330 writes charged/pseudoscalar “mass terms” whose displayed factors do not match the Hessian of the potential printed immediately above them. The later general scalar-sector derivation around lines 8640–8710 defines the matrices explicitly as second derivatives of `V` and agrees with the independent derivation here.

Those early restricted formulas are therefore not used as normalization evidence. The discrepancy is recorded rather than silently repaired.

---

## What was checked against the source

1. DH05 Eq. `gpm`: `G^+=c_beta Phi1^+ + s_beta Phi2^+` at `xi=0`.
2. DH05 Eq. `goldn`: the neutral Goldstone follows the same vacuum direction.
3. DH05 text below `chhiggsmass`: `H^+=-s_beta Phi1^+ + c_beta Phi2^+`.
4. DH05 Eq. `scalareigenstates`: exactly the `R_alpha` convention used here.
5. DH05 Eqs. `hbasis`, `Hbasis`: same beta-alpha decomposition into vacuum-aligned and orthogonal CP-even directions.
6. BFLRS11 later general Hessian: same pre-tadpole CP-even entries, including `3 lambda6`, `3 lambda7`, and the `3/2` off-diagonal coefficients.
7. GHOO18 section 2: explicit `(c_beta,s_beta;-s_beta,c_beta)` rotation for charged and CP-odd fields.

### Independent implementation check against active 2HDMC

After the analytic formulas were frozen, `THDM.cpp` was inspected.

The CP-odd mass path uses

\[
m_A^2=\frac{m_{12}^2}{sc}-\frac{v^2}{2}(2\lambda_5+\lambda_6\cot\beta+\lambda_7\tan\beta),
\]

and the charged mass relation

\[
m_{H^\pm}^2=m_A^2+\frac12v^2(\lambda_5-\lambda_4),
\]

both identical to the derivation.

For the CP-even sector 2HDMC constructs

\[
\begin{aligned}
M_{11}&=m_A^2s^2+v^2(\lambda_1c^2+2\lambda_6sc+\lambda_5s^2),\\
M_{12}&=-m_A^2sc+v^2[(\lambda_3+\lambda_4)sc+\lambda_6c^2+\lambda_7s^2],\\
M_{22}&=m_A^2c^2+v^2(\lambda_2s^2+2\lambda_7sc+\lambda_5c^2),
\end{aligned}
\]

which matches the second analytic representation term by term.

The active physical-input branch uses

\[
\alpha=\beta-\arcsin[\sin(\beta-\alpha)]
\]

with the non-negative `cos(beta-alpha)` branch. At exact alignment this gives `alpha=beta-pi/2`, consistent with `h=rho_v` and `phi=-rho_perp`.

### CAS audit

`appendix_theory/checks/phase4_hessian_check.py` reconstructs the frozen potential, computes all three Hessians, substitutes the Phase-3 tadpoles, and checks the charged and CP-odd factorized matrices and their zero eigenvectors. The script is an algebra audit only; it is not a theoretical source.

---

## What remains uncertain

1. The gauge-coupling interpretation of the vacuum-aligned state is not yet used as proof. Phase 6 must derive it from the kinetic terms.
2. C4 (`kappa_V^phi=0`) therefore remains unverified here.
3. C3 (`kappa_f^phi=-cot beta`) remains unverified until the Type-I Yukawa Lagrangian is expanded using the exact `phi=s rho1-c rho2` convention fixed here.
4. Exact alignment versus `Z6=0` remains only partially verified until the Higgs-basis mass matrix is independently derived.
5. Degenerate CP-even masses require special treatment because `alpha` is then not uniquely fixed by diagonalization.
6. Positive scalar squared masses establish a local quadratic minimum in the physical scalar directions, not by themselves the global electroweak minimum.

---

## Next smallest validation

Proceed to the Type-I Yukawa sector with no coupling table imported. Start from the Yukawa Lagrangian in the fixed basis where only `Phi2` couples to charged fermions, expand

\[
\Phi_2^0=\frac{1}{\sqrt2}(v_2+\rho_2+i\eta_2),
\]

and substitute the already-derived inverse CP-even rotation. Derive `kappa_f^h` and `kappa_f^phi` first at arbitrary `alpha,beta`, then at `sin(beta-alpha)=1`.

Correct if the signs follow from the frozen field convention rather than from a literature coupling table.
