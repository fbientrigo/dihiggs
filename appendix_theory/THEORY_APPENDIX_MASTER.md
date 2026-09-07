# 2HDM theory appendix — master derivation notebook

Issue: `#81 — [paper] Add theory appendix: 2HDM potential, conventions, and relevant couplings`  
Status of this master: **Phases 0–5 incorporated**  
Role: single continuously updated derivation document.  
Rule: no downstream formula is admitted unless its upstream field, sign, basis, and normalization conventions have already been fixed here.

---

# 0. Audit protocol and source hierarchy

The appendix is built in the order

```text
fields
-> scalar potential
-> vacuum
-> minimization
-> quadratic Hessians
-> physical states
-> Yukawa couplings
-> gauge couplings
-> Higgs basis
-> generic-to-Higgs-basis map
-> scalar trilinears
-> loop amplitudes
-> project coordinate X
```

The order is part of the scientific control. In particular, the project coordinate

\[
X\equiv\lambda_6\tan\beta
\]

is **not** allowed to define the theory upstream. It will only be interpreted after the relevant scalar coupling is derived in the original generic basis and translated to the Higgs basis.

## 0.1 Evidence labels

Important statements are classified as one or more of

- `[SOURCE]`: explicitly supported by a cited source;
- `[DERIVED]`: obtained algebraically from previously frozen definitions;
- `[TRANSLATED]`: obtained by an explicit convention map;
- `[IMPLEMENTATION-CHECKED]`: compared only after analytic derivation with the active implementation;
- `[NUMERICALLY-CHECKED]`: checked on explicit numerical points;
- `[PROJECT-DEFINITION]`: a convention/coordinate selected by the project;
- `[NUMERICAL-OBSERVATION]`: empirical behavior of a scan or benchmark;
- `[HYPOTHESIS]`: plausible but not established;
- `[OPEN-QUESTION]`: unresolved.

A source quotation alone is not enough to promote a claim designated for independent derivation.

## 0.2 Primary source packs

Three TeX source packs are the initial analytic corpus:

1. **DH05** — Davidson & Haber, `hep-ph/0504050`, source file `hbasis.tex`.
2. **BFLRS11** — Branco et al., `arXiv:1106.0034`, source file `PhysRep_large.tex`.
3. **GHOO18** — Grzadkowski, Ogreid & Osland, `arXiv:1808.01472`, source file `paper_heavyhiggs_jhep_revised3.tex`.

Source priority for equations is

```text
source TeX
> rendered equation in the same paper
> prose in the same paper
> review discussion
> implementation
> general recollection
```

Implementation is used as a cross-check, never as the origin of a theoretical formula.

---

# Phase 0 — Source and convention inventory

## What is established

The project generic-basis potential convention is chosen to coincide with DH05/BFLRS11:

\[
\begin{aligned}
V={}&m_{11}^2\Phi_1^\dagger\Phi_1
+m_{22}^2\Phi_2^\dagger\Phi_2
-\left[m_{12}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}\right]
+\cdots
\end{aligned}
\]

where the quartic normalization is detailed in Phase 2.

GHOO18 uses instead

\[
V_G\supset-\frac12\left\{
 m_{11,G}^2\Phi_1^\dagger\Phi_1
+m_{22,G}^2\Phi_2^\dagger\Phi_2
+[m_{12,G}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}]
\right\}.
\]

Matching coefficients operator by operator gives

\[
\boxed{
 m_{11,\rm DH}^2=-\frac12m_{11,G}^2,
\qquad
 m_{22,\rm DH}^2=-\frac12m_{22,G}^2,
\qquad
 m_{12,\rm DH}^2=+\frac12m_{12,G}^2.
}
\]

The quartic normalization is the same once the same fields and basis are identified.

## Type-I convention map

The **project** follows the BFLRS11/GHOO18 assignment:

\[
\boxed{\Phi_2\text{ couples to all charged fermions in Type I}.}
\]

BFLRS11 explicitly states this convention. GHOO18 writes

\[
\eta_1^{u,0}=\eta_1^{d,0}=\eta_1^{l,0}=0,
\]

which again means only `Phi2` couples.

DH05 gives a basis-independent definition of Type I and also exhibits a special basis with the opposite doublet assignment. The two special-basis descriptions are related by

\[
\Phi_1\leftrightarrow\Phi_2,
\qquad
\tan\beta\leftrightarrow\cot\beta.
\]

Therefore the phrase “Type I” does **not** by itself identify which symbol `Phi1` or `Phi2` is the Yukawa doublet. That mapping must always be stated.

## CP-even sign map

The DH05 CP-even convention selected as the project reference is

\[
\boxed{
h_{\rm DH}=-s_\alpha\rho_1+c_\alpha\rho_2,
\qquad
H_{\rm DH}=c_\alpha\rho_1+s_\alpha\rho_2.
}
\]

The early BFLRS11 pedagogical display uses the global negatives of these two states. Its nearby coupling table/prose does not consistently carry that field-sign flip, so BFLRS11 is not used to fix project scalar signs.

Phase 4 will derive the project state identity from the mass matrix. The resulting map is

\[
\boxed{h=h_{\rm DH},\qquad \phi=H_{\rm DH}.}
\]

On the matched GHOO18 alignment branch the CP-even non-SM state satisfies

\[
\boxed{H_2^{\rm GHOO}=-\phi_{\rm project}.}
\]

This single field redefinition must be carried simultaneously into Yukawa and scalar-trilinear signs.

## Higgs-basis off-diagonal sign map

DH05 writes

\[
V\supset-[M_{12}^2H_1^\dagger H_2+\mathrm{h.c.}],
\]

while GHOO18 writes

\[
V\supset+[Y_3\mathcal H_1^\dagger\mathcal H_2+\mathrm{h.c.}].
\]

At real Higgs-basis phase `chi=0`, DH05 supplies

\[
Y_3=-M_{12}^2,
\qquad Z_6=\Lambda_6,
\qquad Z_7=\Lambda_7.
\]

Thus source statements

\[
M_{12}^2=+\frac12\Lambda_6v^2
\]

and

\[
Y_3=-\frac12Z_6v^2
\]

are the same physical stationarity equation written with different coefficient names.

This relation is source-verified but remains scheduled for an independent Higgs-basis derivation in Phase 7.

---

# Phase 1 — Scalar doublets, vacuum coordinates, and normalization

## What is established

The model contains two complex electroweak scalar doublets with identical gauge quantum numbers. In the DH05 hypercharge normalization,

\[
\Phi_i\sim(\mathbf 2,Y=1),
\qquad
Q=T_3+\frac{Y}{2}.
\]

Equivalently, with the modern convention `Q=T3+Y`, both doublets have `Y=1/2`.

For the real CP-conserving neutral vacuum,

\[
\boxed{
\Phi_i=
\begin{pmatrix}
\phi_i^+\\[2mm]
\dfrac{v_i+\rho_i+i\eta_i}{\sqrt2}
\end{pmatrix}
}
\]

with

\[
\boxed{
\langle\Phi_i\rangle=
\frac1{\sqrt2}
\begin{pmatrix}0\\v_i\end{pmatrix}.
}
\]

Define

\[
\boxed{v^2=v_1^2+v_2^2},
\qquad
\boxed{\tan\beta=\frac{v_2}{v_1}},
\]

so

\[
v_1=v\cos\beta,
\qquad
v_2=v\sin\beta.
\]

No physical `h`, `phi`, `A`, or `H+` state is assumed at this point.

## Why the lower component receives the VEV

For `Y=1`,

\[
Q_{\rm upper}=+\frac12+\frac12=+1,
\qquad
Q_{\rm lower}=-\frac12+\frac12=0.
\]

An electromagnetic-preserving vacuum therefore lies in the neutral lower component.

## Why the factor `1/sqrt(2)` is fixed

For

\[
\Phi_i^0=\frac{v_i+\rho_i+i\eta_i}{\sqrt2},
\]

the kinetic term gives

\[
|\partial_\mu\Phi_i^0|^2
=\frac12(\partial_\mu\rho_i)^2
+\frac12(\partial_\mu\eta_i)^2,
\]

which is the canonical normalization of real scalar fields.

`tan(beta)` is basis-dependent in a generic 2HDM. It becomes a project parameter only because the generic/Yukawa basis has been explicitly fixed.

---

# Phase 2 — Complete CP-conserving scalar potential

## Operator construction

The gauge-singlet bilinears

\[
B_{ij}=\Phi_i^\dagger\Phi_j
\]

have mass dimension two. Renormalizability permits terms containing one bilinear or a product of two bilinears.

A complete Hermitian operator basis is

\[
B_{11},\ B_{22},\ B_{12}+B_{21},
\]

and

\[
B_{11}^2,
\ B_{22}^2,
\ B_{11}B_{22},
\ B_{12}B_{21},
\ B_{12}^2+\mathrm{h.c.},
\ B_{11}B_{12}+\mathrm{h.c.},
\ B_{22}B_{12}+\mathrm{h.c.}.
\]

Using the DH05/BFLRS11 coefficient convention gives

\[
\boxed{
\begin{aligned}
V={}&m_{11}^2\Phi_1^\dagger\Phi_1
+m_{22}^2\Phi_2^\dagger\Phi_2
-\left[m_{12}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}\right]\\
&+\frac{\lambda_1}{2}(\Phi_1^\dagger\Phi_1)^2
+\frac{\lambda_2}{2}(\Phi_2^\dagger\Phi_2)^2
+\lambda_3(\Phi_1^\dagger\Phi_1)(\Phi_2^\dagger\Phi_2)\\
&+\lambda_4(\Phi_1^\dagger\Phi_2)(\Phi_2^\dagger\Phi_1)\\
&+\left\{
\frac{\lambda_5}{2}(\Phi_1^\dagger\Phi_2)^2
+\left[\lambda_6(\Phi_1^\dagger\Phi_1)
+\lambda_7(\Phi_2^\dagger\Phi_2)\right]\Phi_1^\dagger\Phi_2
+\mathrm{h.c.}
\right\}.
\end{aligned}}
\]

Hermiticity requires

\[
m_{11}^2,m_{22}^2,\lambda_{1,2,3,4}\in\mathbb R,
\]

while `m12^2, lambda5, lambda6, lambda7` may be complex in the general theory.

The project specializes to a real CP-conserving basis:

\[
\boxed{m_{12}^2,\lambda_5,\lambda_6,\lambda_7\in\mathbb R}.
\]

The campaign choice `lambda7=0` is **not** imposed in the definition of the theory.

---

# Phase 3 — Vacuum potential and minimization

Define

\[
\lambda_{345}=\lambda_3+\lambda_4+\lambda_5.
\]

At the neutral real vacuum,

\[
\Phi_1^\dagger\Phi_1=\frac{v_1^2}{2},
\qquad
\Phi_2^\dagger\Phi_2=\frac{v_2^2}{2},
\qquad
\Phi_1^\dagger\Phi_2=\frac{v_1v_2}{2}.
\]

Term by term,

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

The Hermitian-conjugate factors are essential. For example,

\[
-[m_{12}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}]
\to-m_{12}^2v_1v_2,
\]

and

\[
[\lambda_6(\Phi_1^\dagger\Phi_1)(\Phi_1^\dagger\Phi_2)+\mathrm{h.c.}]
\to\frac12\lambda_6v_1^3v_2.
\]

## Tadpole equations

Differentiate explicitly:

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

The asymmetric factors `3 lambda6` and `3 lambda7` arise directly from

\[
\frac{\partial}{\partial v_1}(v_1^3v_2)=3v_1^2v_2,
\qquad
\frac{\partial}{\partial v_2}(v_1v_2^3)=3v_1v_2^2.
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
\right],\\[1mm]
m_{22}^2={}&m_{12}^2\frac{v_1}{v_2}
-\frac12\left[
\lambda_2v_2^2+\lambda_{345}v_1^2
+\lambda_6\frac{v_1^3}{v_2}
+3\lambda_7v_1v_2
\right].
\end{aligned}}
\]

With `s=sin beta`, `c=cos beta`,

\[
\boxed{
\begin{aligned}
m_{11}^2={}&m_{12}^2\tan\beta
-\frac{v^2}{2}\left[
\lambda_1c^2+\lambda_{345}s^2+3\lambda_6sc
+\lambda_7\frac{s^3}{c}
\right],\\[1mm]
m_{22}^2={}&m_{12}^2\cot\beta
-\frac{v^2}{2}\left[
\lambda_2s^2+\lambda_{345}c^2
+\lambda_6\frac{c^3}{s}+3\lambda_7sc
\right].
\end{aligned}}
\]

Only now define the project soft coordinate

\[
\boxed{M^2\equiv\frac{m_{12}^2}{s_\beta c_\beta}},
\qquad
m_{12}^2=M^2s_\beta c_\beta.
\]

Thus

\[
\boxed{m_{22}^2\neq m_{12}^2\neq M^2}
\]

generically. `m22^2` is a diagonal potential coefficient, `m12^2` is the off-diagonal coefficient, and `M^2` is a derived rescaling of the latter.

The derived stationarity equations were subsequently checked against DH05 and the active 2HDMC `set_param_gen` path.

---

# Phase 4 — Quadratic expansion, masses, Goldstones, and physical states

## Charged sector

Define

\[
V_\pm^{(2)}=(\phi_1^-,\phi_2^-)\mathcal M_\pm^2
\begin{pmatrix}\phi_1^+\\\phi_2^+\end{pmatrix}.
\]

Before tadpole elimination,

\[
\mathcal M_\pm^2=
\begin{pmatrix}X_\pm&Y_\pm\\Y_\pm&Z_\pm\end{pmatrix}
\]

with

\[
\begin{aligned}
X_\pm={}&m_{11}^2+\frac12\lambda_1v_1^2+\frac12\lambda_3v_2^2+\lambda_6v_1v_2,\\
Z_\pm={}&m_{22}^2+\frac12\lambda_2v_2^2+\frac12\lambda_3v_1^2+\lambda_7v_1v_2,\\
Y_\pm={}&-m_{12}^2+\frac12(\lambda_4+\lambda_5)v_1v_2
+\frac12\lambda_6v_1^2+\frac12\lambda_7v_2^2.
\end{aligned}
\]

After inserting the Phase-3 stationarity equations,

\[
D_\pm=m_{12}^2-\frac12\left[(\lambda_4+\lambda_5)v_1v_2+\lambda_6v_1^2+\lambda_7v_2^2\right]
\]

and

\[
\boxed{
\mathcal M_\pm^2=D_\pm
\begin{pmatrix}
v_2/v_1&-1\\
-1&v_1/v_2
\end{pmatrix}.
}
\]

Therefore

\[
\mathcal M_\pm^2
\begin{pmatrix}v_1\\v_2\end{pmatrix}=0.
\]

The zero and orthogonal directions are

\[
\boxed{G^+=c_\beta\phi_1^+ + s_\beta\phi_2^+},
\qquad
\boxed{H^+=-s_\beta\phi_1^+ + c_\beta\phi_2^+}.
\]

The physical mass is

\[
\boxed{
m_{H^\pm}^2
=M^2-\frac{v^2}{2}
\left(\lambda_4+\lambda_5+\lambda_6\cot\beta+\lambda_7\tan\beta\right).
}
\]

## CP-odd sector

Similarly,

\[
V_A^{(2)}=\frac12(\eta_1,\eta_2)\mathcal M_A^2
\begin{pmatrix}\eta_1\\\eta_2\end{pmatrix}.
\]

After tadpole elimination,

\[
D_A=m_{12}^2-\lambda_5v_1v_2
-\frac12(\lambda_6v_1^2+\lambda_7v_2^2)
\]

and

\[
\boxed{
\mathcal M_A^2=D_A
\begin{pmatrix}
v_2/v_1&-1\\
-1&v_1/v_2
\end{pmatrix}.
}
\]

Thus

\[
\boxed{G^0=c_\beta\eta_1+s_\beta\eta_2},
\qquad
\boxed{A=-s_\beta\eta_1+c_\beta\eta_2},
\]

with

\[
\boxed{
m_A^2=M^2-\frac{v^2}{2}
\left(2\lambda_5+\lambda_6\cot\beta+\lambda_7\tan\beta\right).
}
\]

The charged/odd splitting is

\[
\boxed{m_{H^\pm}^2-m_A^2=\frac{v^2}{2}(\lambda_5-\lambda_4)}.
\]

The `lambda6,lambda7` terms cancel exactly in this difference.

## CP-even sector

For

\[
\rho=(\rho_1,\rho_2)^T,
\qquad
V_\rho^{(2)}=\frac12\rho^T\mathcal M_\rho^2\rho,
\]

the pre-tadpole Hessian entries are

\[
\begin{aligned}
\mathcal M_{11}^2={}&m_{11}^2+\frac32\lambda_1v_1^2
+\frac12\lambda_{345}v_2^2+3\lambda_6v_1v_2,\\
\mathcal M_{22}^2={}&m_{22}^2+\frac32\lambda_2v_2^2
+\frac12\lambda_{345}v_1^2+3\lambda_7v_1v_2,\\
\mathcal M_{12}^2={}&-m_{12}^2+\lambda_{345}v_1v_2
+\frac32\lambda_6v_1^2+\frac32\lambda_7v_2^2.
\end{aligned}
\]

After tadpole elimination,

\[
\boxed{
\begin{aligned}
\mathcal M_{11}^2={}&M^2s^2+v^2\left[\lambda_1c^2+\frac32\lambda_6sc-\frac12\lambda_7\frac{s^3}{c}\right],\\
\mathcal M_{22}^2={}&M^2c^2+v^2\left[\lambda_2s^2-\frac12\lambda_6\frac{c^3}{s}+\frac32\lambda_7sc\right],\\
\mathcal M_{12}^2={}&-M^2sc+v^2\left[\lambda_{345}sc+\frac32\lambda_6c^2+\frac32\lambda_7s^2\right].
\end{aligned}}
\]

An equivalent representation is

\[
\boxed{
\begin{aligned}
\mathcal M_{11}^2={}&m_A^2s^2+v^2(\lambda_1c^2+\lambda_5s^2+2\lambda_6sc),\\
\mathcal M_{22}^2={}&m_A^2c^2+v^2(\lambda_2s^2+\lambda_5c^2+2\lambda_7sc),\\
\mathcal M_{12}^2={}&-m_A^2sc+v^2[(\lambda_3+\lambda_4)sc+\lambda_6c^2+\lambda_7s^2].
\end{aligned}}
\]

This second form matches the active 2HDMC construction term by term.

## CP-even diagonalization

Freeze the DH05 sign convention

\[
\boxed{
\begin{pmatrix}h\\\phi\end{pmatrix}
=
\begin{pmatrix}
-s_\alpha&c_\alpha\\
c_\alpha&s_\alpha
\end{pmatrix}
\begin{pmatrix}\rho_1\\\rho_2\end{pmatrix}.
}
\]

Diagonalization yields

\[
\boxed{
\tan2\alpha=\frac{2\mathcal M_{12}^2}{\mathcal M_{11}^2-\mathcal M_{22}^2}
}
\]

with the quadrant fixed by the chosen eigenvector branch rather than by the tangent alone.

Define the CP-even vacuum and orthogonal directions

\[
\rho_v=c_\beta\rho_1+s_\beta\rho_2,
\qquad
\rho_\perp=-s_\beta\rho_1+c_\beta\rho_2.
\]

Then

\[
\boxed{
h=s_{\beta-\alpha}\rho_v+c_{\beta-\alpha}\rho_\perp,
}
\]

\[
\boxed{
\phi=c_{\beta-\alpha}\rho_v-s_{\beta-\alpha}\rho_\perp.
}
\]

On the project exact-alignment branch

\[
\boxed{s_{\beta-\alpha}=+1},
\]

we have

\[
\boxed{h=\rho_v=c_\beta\rho_1+s_\beta\rho_2},
\]

\[
\boxed{\phi=-\rho_\perp=s_\beta\rho_1-c_\beta\rho_2}.
\]

This fixes the real-scalar sign that all later Yukawa and trilinear translations must respect.

---

# Phase 5 — Type-I Yukawa sector

## Type-I Lagrangian in the project basis

The project basis is

\[
\boxed{\Phi_2\text{ Yukawa-active},\qquad \Phi_1\text{ Yukawa-inert for charged fermions}.}
\]

Suppressing flavour indices,

\[
\boxed{
-\mathcal L_Y=
 y_d\bar Q_L\Phi_2d_R
+y_u\bar Q_L\widetilde\Phi_2u_R
+y_\ell\bar L_L\Phi_2\ell_R
+\mathrm{h.c.}
}
\]

with

\[
\widetilde\Phi_2=i\sigma_2\Phi_2^*.
\]

Using

\[
\Phi_2=
\begin{pmatrix}
\phi_2^+\\
(v_2+\rho_2+i\eta_2)/\sqrt2
\end{pmatrix},
\qquad
\widetilde\Phi_2=
\begin{pmatrix}
(v_2+\rho_2-i\eta_2)/\sqrt2\\
-\phi_2^-
\end{pmatrix},
\]

the CP-even neutral factor is identical for up, down, and lepton sectors.

## Fermion masses

The VEV term gives

\[
\mathcal L_{\rm mass}=-\frac{y_fv_2}{\sqrt2}\bar ff,
\]

hence

\[
\boxed{m_f=\frac{y_fv_2}{\sqrt2}}
\]

and therefore

\[
\boxed{y_f=\frac{\sqrt2m_f}{vs_\beta}}.
\]

The CP-even interaction before rotating to physical scalars is

\[
\boxed{
\mathcal L_Y^{\rm CP-even}
=-\frac{m_f}{vs_\beta}\rho_2\bar ff.
}
\]

## Inverse CP-even rotation

Because the Phase-4 matrix is orthogonal and symmetric,

\[
\boxed{\rho_1=-s_\alpha h+c_\alpha\phi},
\qquad
\boxed{\rho_2=c_\alpha h+s_\alpha\phi}.
\]

Thus

\[
\begin{aligned}
\mathcal L_Y^{\rm CP-even}
&=-\frac{m_f}{vs_\beta}
(c_\alpha h+s_\alpha\phi)\bar ff\\
&=-\frac{m_f}{v}
\left[\frac{c_\alpha}{s_\beta}h
+\frac{s_\alpha}{s_\beta}\phi\right]\bar ff.
\end{aligned}
\]

Define

\[
\mathcal L_Y^{\rm CP-even}
\equiv-\frac{m_f}{v}
(\kappa_f^h h+\kappa_f^\phi\phi)\bar ff.
\]

Then

\[
\boxed{\kappa_f^h=\frac{c_\alpha}{s_\beta}},
\qquad
\boxed{\kappa_f^\phi=\frac{s_\alpha}{s_\beta}},
\qquad f=u,d,\ell.
\]

## Exact expressions in `beta-alpha`

Using

\[
\cos\alpha
=\cos\beta\cos(\beta-\alpha)
+\sin\beta\sin(\beta-\alpha),
\]

we obtain

\[
\boxed{
\kappa_f^h
=s_{\beta-\alpha}+c_{\beta-\alpha}\cot\beta.
}
\]

Using

\[
\sin\alpha
=\sin\beta\cos(\beta-\alpha)
-\cos\beta\sin(\beta-\alpha),
\]

we obtain

\[
\boxed{
\kappa_f^\phi
=c_{\beta-\alpha}-s_{\beta-\alpha}\cot\beta.
}
\]

These are exact tree-level identities in the project Type-I basis.

## Exact alignment

For

\[
s_{\beta-\alpha}=+1,
\qquad
c_{\beta-\alpha}=0,
\]

we obtain

\[
\boxed{\kappa_f^h=1},
\qquad
\boxed{\kappa_f^\phi=-\cot\beta}.
\]

The same sign is visible directly from the aligned fields. Since

\[
\rho_2=s_\beta h-c_\beta\phi,
\]

then

\[
\mathcal L_Y^{\rm CP-even}
=-\frac{m_f}{v}
\left(h-\cot\beta\,\phi\right)\bar ff.
\]

Therefore C3 is **VERIFIED** without using a literature coupling table as input.

## CP-odd sign check

From Phase 4,

\[
\eta_2=s_\beta G^0+c_\beta A.
\]

The imaginary part of `Phi2` and `tilde Phi2` enters with opposite signs. In the BFLRS11 convention

\[
\mathcal L_Y\supset+i\frac{m_f}{v}\xi_A^fA\bar f\gamma_5f,
\]

this gives

\[
\boxed{
\xi_A^u=+\cot\beta,
\qquad
\xi_A^d=\xi_A^\ell=-\cot\beta.
}
\]

This is a secondary internal consistency check, not needed to establish C3.

## Conditional fermionic-width scaling

For a CP-even scalar `S` with

\[
\mathcal L\supset-\frac{m_f}{v}\kappa_f^SS\bar ff,
\]

the tree-level width is

\[
\Gamma(S\to f\bar f)
=N_c\frac{m_S}{8\pi}\frac{m_f^2}{v^2}|\kappa_f^S|^2
\left(1-\frac{4m_f^2}{m_S^2}\right)^{3/2}.
\]

Hence, at fixed scalar mass, fermion-mass prescription, and radiative corrections,

\[
\boxed{\Gamma(\phi\to f\bar f)\propto\cot^2\beta}
\]

in exact alignment.

This does **not** imply

\[
\Gamma_{\rm total}\propto\cot^2\beta,
\qquad
\mathrm{BR}\propto\cot^2\beta,
\qquad
c\tau\propto\tan^2\beta
\]

globally, because other open channels have independent parameter dependence.

---

# Cross-source convention table after Phase 5

| Object | Project | DH05 | BFLRS11 | GHOO18 |
|---|---|---|---|---|
| generic potential | DH normalization | reference | same notation-1 normalization | quadratic symbols rescaled/sign-flipped |
| `tan beta` | `v2/v1` | same ratio | same | same |
| Type-I Yukawa doublet | `Phi2` | Type-I is basis-invariant; displayed special basis may use opposite assignment | `Phi2` | `Phi2` via `eta1^f=0` |
| CP-even `h` | `h_DH` | `-s_alpha rho1+c_alpha rho2` | early displayed `h_B=-h_DH` | matched vacuum-aligned state requires separate `H_i` map |
| CP-even `phi` | `H_DH` | `c_alpha rho1+s_alpha rho2` | early displayed `H_B=-H_DH` | matched AL `H2=-phi_project` |
| exact-alignment `kappa_f^phi` | `-cot beta` | compatible after Type-I basis map | table magnitude/form agrees but early field sign display is inconsistent | source `H2` sign flips to project `phi` |
| `Y3` convention | future Phase 7: `+Y3 H1†H2` | `Y3=-M12^2` | negative off-diagonal mass convention | `+Y3` directly |
| `Z7` | Higgs-basis quantity | `Lambda7 -> Z7` after phase choice | barred quartic layer | direct `Z7` |

---

# Current validation ledger after Phase 5

| Claim | Status | Result now |
|---|---|---|
| C1 project-compatible scalar potential | **VERIFIED** | DH05/BFLRS11 convention reconstructed and implementation-checked |
| C2 physical scalar rotations/state map | **VERIFIED** | `h=h_DH`, `phi=H_DH`; exact AL `phi=s rho1-c rho2` |
| C3 Type-I exact-alignment `kappa_f^phi` | **VERIFIED** | `kappa_f^phi=-cot beta` for `u,d,l` |
| C4 exact-alignment `kappa_V^phi` | **NOT VERIFIED** | Phase 6 kinetic-term derivation required |
| C5 `Y3=-Z6 v^2/2` | **PARTIALLY VERIFIED** | source-verified, independent Phase-7 derivation pending |
| C6 exact alignment and `Z6` | **PARTIALLY VERIFIED** | generic mass geometry derived; Higgs-basis off-diagonal identification pending |
| C7 exact `phi H+H-` coupling | **NOT VERIFIED** | Phase 9 |
| C8 exact generic-to-Higgs `Z7` | **PARTIALLY VERIFIED** | source relation known; independent rotation pending |
| C9 large-`tan beta`, `lambda7=0` limit of `Z7` | **NOT VERIFIED** | blocked on C8 |
| C10 large-`tan beta` trilinear approximation | **NOT VERIFIED** | blocked on C7–C9 |
| C11 re-expression with `X=lambda6 tan beta` | **NOT VERIFIED** | intentionally deferred |

Phase gates:

```text
PHASE_0_PASS = PASS
PHASE_1_PASS = PASS
PHASE_2_PASS = PASS
PHASE_3_PASS = PASS
PHASE_4_PASS = PASS
PHASE_5_PASS = PASS
```

---

# Open issues relevant to the derivation

1. **Gauge coupling C4:** derive from the kinetic terms; do not infer solely from vacuum orthogonality.
2. **Exact alignment and `Z6`:** the finite-mass Higgs-basis statement must be derived from the explicitly rotated mass matrix. Decoupling with nonzero `Z6` must be phrased separately as an asymptotic alignment route.
3. **Residual Higgs-basis sign:** Phase 7 will freeze
   \[
   H_1=c_\beta\Phi_1+s_\beta\Phi_2,
   \qquad
   H_2=-s_\beta\Phi_1+c_\beta\Phi_2.
   \]
   With this choice, the relation between `phi` and `sqrt(2) Re H2^0` is no longer adjustable downstream.
4. **Global vacuum:** positive physical scalar masses are local quadratic conditions, not a proof that the neutral stationary point is the global minimum.
5. **Exact CP-even degeneracy:** at exact degeneracy the mixing angle is not unique; state identity must then be phrased via subspace/projector information.
6. **Loop amplitudes:** the current three source packs do not contain a complete charged-scalar-loop normalization for both `phi -> gamma gamma` and `phi -> Z gamma`; a dedicated primary source is required before the loop phase.

---

# Next phase

Phase 6 must start from

\[
\mathcal L_{\rm kin}
=\sum_{i=1}^2(D_\mu\Phi_i)^\dagger(D^\mu\Phi_i)
\]

and independently derive the terms linear in `rho_i` and quadratic in `W` or `Z`.

The target is to derive, rather than quote,

\[
\kappa_V^h,
\qquad
\kappa_V^\phi,
\]

first at arbitrary `alpha,beta`, and only then apply `s_{beta-alpha}=1`.

No Higgs-basis trilinear, loop amplitude, or `X` interpretation is allowed before that chain is complete.
