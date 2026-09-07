# Phase 6 — Gauge couplings from the scalar kinetic terms

This phase derives the neutral CP-even couplings to electroweak vector-boson pairs directly from the canonical kinetic terms. No `hVV/HVV` coupling table is used as input.

The field convention inherited from Phases 1–5 is

\[
\Phi_i=\begin{pmatrix}\phi_i^+\\(v_i+\rho_i+i\eta_i)/\sqrt2\end{pmatrix},
\qquad v_1=vc_\beta,\quad v_2=vs_\beta,
\qquad v^2=v_1^2+v_2^2,
\]

and the CP-even mass states are

\[
\begin{pmatrix}h\\\phi\end{pmatrix}
=\begin{pmatrix}-s_\alpha&c_\alpha\\c_\alpha&s_\alpha\end{pmatrix}
\begin{pmatrix}\rho_1\\\rho_2\end{pmatrix}.
\]

Hence

\[
\rho_1=-s_\alpha h+c_\alpha\phi,
\qquad
\rho_2=c_\alpha h+s_\alpha\phi.
\]

---

## What is established

[DERIVED] Electroweak gauge-boson masses depend only on the total vacuum norm,

\[
\boxed{m_W^2=\frac{g^2v^2}{4}},
\qquad
\boxed{m_Z^2=\frac{(g^2+g'^2)v^2}{4}}.
\]

[DERIVED] The terms linear in a neutral CP-even fluctuation and quadratic in gauge fields depend only on

\[
\boxed{v_1\rho_1+v_2\rho_2=v\rho_v},
\qquad
\rho_v\equiv c_\beta\rho_1+s_\beta\rho_2.
\]

Therefore the unique tree-level CP-even direction with a linear `VV` coupling is the direction parallel to the VEV. The orthogonal direction

\[
\rho_\perp=-s_\beta\rho_1+c_\beta\rho_2
\]

has no tree-level linear `WW` or `ZZ` coupling.

[DERIVED] In the mass basis,

\[
\boxed{\rho_v=s_{\beta-\alpha}h+c_{\beta-\alpha}\phi}.
\]

Defining

\[
\mathcal L_{SVV}=\kappa_V^S\left[
\frac{2m_W^2}{v}SW_\mu^+W^{-\mu}
+\frac{m_Z^2}{v}SZ_\mu Z^\mu
\right],
\]

we obtain

\[
\boxed{\kappa_V^h=\sin(\beta-\alpha)},
\qquad
\boxed{\kappa_V^\phi=\cos(\beta-\alpha)}.
\]

On the project exact-alignment branch,

\[
\sin(\beta-\alpha)=1,
\qquad
\cos(\beta-\alpha)=0,
\]

so

\[
\boxed{\kappa_V^h=1},
\qquad
\boxed{\kappa_V^\phi=0}.
\]

[DERIVED] No term linear in `eta_i` occurs in the gauge-boson mass operator, so the CP-odd state `A` has no tree-level `AWW` or `AZZ` coupling.

[DERIVED] Orthogonality also implies

\[
\rho_1^2+\rho_2^2=h^2+\phi^2,
\]

so the diagonal CP-even quartic gauge couplings are angle-independent and the mixed `VVhphi` quartic cancels.

---

## Derivation

### 1. Covariant derivative and hypercharge convention

Use modern hypercharge notation

\[
Q=T^3+Y,
\qquad Y(\Phi_i)=\frac12,
\]

with

\[
\boxed{D_\mu=\partial_\mu+i g\frac{\sigma^a}{2}W_\mu^a+i g'\frac12B_\mu}.
\]

Davidson–Haber call the doublet hypercharge `Y=1` while defining `Q=T3+Y/2`; this is the same physical normalization because the covariant derivative again contains `g'/2` on the doublet.

Define

\[
W_\mu^\pm=\frac{W_\mu^1\mp iW_\mu^2}{\sqrt2},
\]

and

\[
A_\mu=s_WW_\mu^3+c_WB_\mu,
\qquad
Z_\mu=c_WW_\mu^3-s_WB_\mu,
\]

where

\[
s_W=\frac{g'}{\sqrt{g^2+g'^2}},
\qquad
c_W=\frac{g}{\sqrt{g^2+g'^2}},
\qquad
g_Z\equiv\sqrt{g^2+g'^2}=\frac{g}{c_W}.
\]

For the neutral CP-even extraction set

\[
\Phi_i\rightarrow\frac1{\sqrt2}\begin{pmatrix}0\\x_i\end{pmatrix},
\qquad x_i=v_i+\rho_i.
\]

### 2. Charged gauge sector

Using `T^a=sigma^a/2`,

\[
(T^1W^1+T^2W^2)
\frac1{\sqrt2}\begin{pmatrix}0\\x_i\end{pmatrix}
=\begin{pmatrix}x_iW^+/2\\0\end{pmatrix}.
\]

Therefore

\[
(D_\mu\Phi_i)^\dagger(D^\mu\Phi_i)
\supset\frac{g^2}{4}x_i^2W_\mu^+W^{-\mu}.
\]

Summing both doublets,

\[
\mathcal L_{\rm kin}\supset
\frac{g^2}{4}\left[(v_1+\rho_1)^2+(v_2+\rho_2)^2\right]W^+W^-.
\]

Expanding,

\[
\mathcal L_{\rm kin}\supset
\frac{g^2}{4}(v_1^2+v_2^2)W^+W^-
+\frac{g^2}{2}(v_1\rho_1+v_2\rho_2)W^+W^-+\cdots.
\]

Hence

\[
\boxed{m_W^2=\frac{g^2v^2}{4}},
\]

and

\[
\boxed{\mathcal L_{WW,\rm linear}=\frac{2m_W^2}{v}\rho_vW_\mu^+W^{-\mu}}.
\]

### 3. Neutral gauge sector and photon cancellation

The lower neutral component has

\[
T^3=-\frac12,
\qquad Y=+\frac12.
\]

Thus its neutral gauge factor is

\[
-i\frac g2W^3+i\frac{g'}2B
=\frac i2(-gW^3+g'B)
=-i\frac{g_Z}{2}Z.
\]

The photon cancels exactly because the neutral vacuum component has `Q=0`.

Therefore

\[
(D_\mu\Phi_i)^\dagger(D^\mu\Phi_i)
\supset\frac{g_Z^2}{8}x_i^2Z_\mu Z^\mu.
\]

After summing the doublets,

\[
\mathcal L_{\rm kin}\supset
\frac{g_Z^2}{8}\left[(v_1+\rho_1)^2+(v_2+\rho_2)^2\right]Z^2.
\]

Matching the vacuum term to `(1/2)m_Z^2 Z_mu Z^mu` gives

\[
\boxed{m_Z^2=\frac{g_Z^2v^2}{4}},
\]

and the linear term is

\[
\boxed{\mathcal L_{ZZ,\rm linear}=\frac{m_Z^2}{v}\rho_vZ_\mu Z^\mu}.
\]

The displayed factor-of-two difference between the `WW` and `ZZ` Lagrangian coefficients is the standard identical-field normalization for the real `Z` pair. The corresponding Feynman rules are

\[
SW_\mu^+W_\nu^-:\quad i\kappa_V^S\frac{2m_W^2}{v}g_{\mu\nu},
\]

\[
SZ_\mu Z_\nu:\quad i\kappa_V^S\frac{2m_Z^2}{v}g_{\mu\nu}.
\]

### 4. Vacuum geometry

Combining both sectors,

\[
\mathcal L_{VV,\rm linear}
=(v_1\rho_1+v_2\rho_2)
\left[\frac{g^2}{2}W^+W^-+\frac{g_Z^2}{4}ZZ\right].
\]

But

\[
v_1\rho_1+v_2\rho_2=v(c_\beta\rho_1+s_\beta\rho_2)=v\rho_v.
\]

Thus

\[
\boxed{\text{tree-level }SVV\text{ coupling}\propto
\text{projection of }S\text{ onto the VEV direction}.}
\]

This result uses only gauge invariance, canonical kinetic terms and VEV geometry; it does not use the scalar potential.

### 5. Projection into `h,phi`

Insert

\[
\rho_1=-s_\alpha h+c_\alpha\phi,
\qquad
\rho_2=c_\alpha h+s_\alpha\phi.
\]

Then

\[
\begin{aligned}
\rho_v
&=c_\beta\rho_1+s_\beta\rho_2\\
&=(-c_\beta s_\alpha+s_\beta c_\alpha)h
 +(c_\beta c_\alpha+s_\beta s_\alpha)\phi\\
&=\sin(\beta-\alpha)h+\cos(\beta-\alpha)\phi.
\end{aligned}
\]

Therefore

\[
\boxed{\kappa_V^h=s_{\beta-\alpha}},
\qquad
\boxed{\kappa_V^\phi=c_{\beta-\alpha}}.
\]

### 6. Exact alignment acquires its physical meaning

Phase 4 established purely from the mass matrix that, on the project branch,

\[
s_{\beta-\alpha}=1
\Longrightarrow
h=\rho_v,
\qquad
\phi=-\rho_\perp.
\]

Phase 6 independently shows

\[
\boxed{h=\rho_v\Rightarrow hVV\text{ is exactly SM-like at tree level}},
\]

\[
\boxed{\phi\perp\rho_v\Rightarrow\phi VV=0\text{ at tree level}}.
\]

This does **not** yet prove `Z6=0`; that is a separate Higgs-basis statement reserved for Phase 7.

### 7. CP-odd field

Keeping the imaginary neutral fields,

\[
|v_i+\rho_i+i\eta_i|^2=(v_i+\rho_i)^2+\eta_i^2.
\]

There is no term linear in `eta_i`. Hence no linear combination of `eta_1,eta_2`, including the physical pseudoscalar `A`, has a tree-level `AVV` vertex of the mass-generated form.

### 8. Quartic gauge consistency check

The CP-even quadratic part of the kinetic expansion contains

\[
\frac{g^2}{4}(\rho_1^2+\rho_2^2)W^+W^-
+\frac{g_Z^2}{8}(\rho_1^2+\rho_2^2)ZZ.
\]

Because the CP-even rotation is orthogonal,

\[
\boxed{\rho_1^2+\rho_2^2=h^2+\phi^2}.
\]

This makes the diagonal `VVhh` and `VVphiphi` quartics angle-independent and cancels the mixed `VVhphi` term.

---

## Convention map

| Object | Project/DH convention | Literature translation | 2HDMC | Status |
|---|---|---|---|---|
| Hypercharge | modern `Y=1/2`, `Q=T3+Y`; DH equivalent `Y=1`, `Q=T3+Y/2` | same EW representation | SM `g,g'` | `[TRANSLATED]` |
| vacuum CP-even direction | `rho_v=c_beta rho1+s_beta rho2` | alignment direction parallel to VEV | first `q_{k1}` component | `[DERIVED][SOURCE]` |
| `kappa_V^h` | `sin(beta-alpha)` | DH/BFLRS standard modifier | `q_11=sba` | `[DERIVED][SOURCE][IMPLEMENTATION-CHECKED]` |
| `kappa_V^phi` | `cos(beta-alpha)` | DH heavy-state modifier | `q_21=cba` | `[DERIVED][SOURCE][IMPLEMENTATION-CHECKED]` |
| exact alignment | `h=rho_v`, `phi=-rho_perp` | aligned state SM-like | `sba=1,cba=0` | `[DERIVED][TRANSLATED]` |
| `AVV` | zero at tree level | standard CP-odd result | `q_31=0` | `[DERIVED][IMPLEMENTATION-CHECKED]` |

BFLRS11's early simple field display uses global negatives of the DH `h,H` fields while nearby phenomenology prose uses the standard positive modifiers. The project does not choose between those source snippets: its signs come from the kinetic derivation in the already-frozen DH field convention.

---

## What was checked against the source

The source comparison was performed only after the analytic derivation was frozen.

1. DH05 gives `h_DH VV` proportional to `sin(beta-alpha)` and `H_DH VV` proportional to `cos(beta-alpha)`.
2. BFLRS11 gives the same physical modifiers and no `AVV` coupling, after accounting for its displayed field-sign convention.
3. GHOO18 defines exact alignment geometrically as a mass eigenstate aligned with the VEV and gives the aligned state the SM `WW/ZZ` vertices while orthogonal neutral states have zero `VV` in exact alignment.
4. Active 2HDMC defines

\[
q_{k1}=(s_{\beta-\alpha},c_{\beta-\alpha},0,i)
\]

for `(h,H,A,H+)`, and `THDM::get_coupling_vvh` returns

\[
c(ZZh_k)=i\frac{gM_Z}{c_W}\operatorname{Re}q_{k1},
\qquad
c(WWh_k)=igM_W\operatorname{Re}q_{k1}.
\]

Since `gM_W=2m_W^2/v` and `gM_Z/c_W=2m_Z^2/v`, this matches the independently derived modifier structure exactly.

The symbolic audit `appendix_theory/checks/phase6_gauge_check.py` verifies the trigonometric projection, exact-alignment limit and orthogonal quartic norm.

---

## What remains uncertain

1. `kappa_V^phi=0` is now **VERIFIED** for exact alignment in the declared project convention.
2. Exact alignment `iff Z6=0` is **not yet promoted**. Phase 7 must rotate the potential explicitly to Higgs basis and identify the CP-even off-diagonal mass entry.
3. The residual Higgs-basis sign `H2 -> -H2` still changes `Y3,Z6,Z7` and all odd-`H2` scalar couplings together; it must be frozen before interpreting trilinear signs.
4. A new implementation caution is recorded for Phase 7/9: 2HDMC's `get_param_higgs` returns variables named `Lambda6,Lambda7`, whereas `get_coupling_hhh` subsequently defines local `Z6=-l6` and `Z7=-l7`. This sign layer must be audited before mapping 2HDMC trilinears to GHOO/DH notation.
5. `phi -> gamma gamma` and `phi -> Z gamma` can still occur through loops; `kappa_V^phi=0` only removes the tree-level linear `phi WW/ZZ` vertices.

---

## Next smallest validation

Proceed to the Higgs-basis construction using the fixed rotation

\[
\boxed{H_1=c_\beta\Phi_1+s_\beta\Phi_2},
\qquad
\boxed{H_2=-s_\beta\Phi_1+c_\beta\Phi_2}.
\]

Phase 7 must derive, not quote:

1. `H1` carries the full VEV and `H2` has zero VEV;
2. the complete quadratic and quartic Higgs-basis potential;
3. the map `(m_ij^2,lambda_i,beta)->(Y_i,Z_i)`;
4. the Higgs-basis stationarity relations, including the sign of `Y3=-Z6 v^2/2` in the selected convention;
5. the CP-even Higgs-basis mass matrix and the precise non-degenerate condition connecting its off-diagonal entry to alignment.

Only after those steps should `Z7` or any `phi H+H-` trilinear be interpreted.
