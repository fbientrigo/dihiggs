# Phase 6 — Gauge couplings from the scalar kinetic terms

This phase derives the couplings of the CP-even neutral scalars to electroweak gauge-boson pairs directly from the gauge-covariant kinetic terms. No `hVV/HVV` coupling table is used as input.

The scalar-field convention inherited from the previous phases is

\[
\Phi_i=
\begin{pmatrix}
\phi_i^+\\[1mm]
(v_i+\rho_i+i\eta_i)/\sqrt2
\end{pmatrix},
\qquad
v_1=v c_\beta,\quad v_2=v s_\beta,
\qquad
v^2=v_1^2+v_2^2.
\]

The CP-even mass eigenstates are defined by the Davidson--Haber sign convention frozen in Phase 4,

\[
\begin{pmatrix}h\\\phi\end{pmatrix}
=
\begin{pmatrix}
-s_\alpha&c_\alpha\\
c_\alpha&s_\alpha
\end{pmatrix}
\begin{pmatrix}\rho_1\\\rho_2\end{pmatrix}.
\]

Equivalently,

\[
\rho_1=-s_\alpha h+c_\alpha\phi,
\qquad
\rho_2=c_\alpha h+s_\alpha\phi.
\]

---

## What is established

[DERIVED] The electroweak gauge-boson masses arise from the single combination

\[
\boxed{v^2=v_1^2+v_2^2},
\]

with

\[
\boxed{m_W^2=\frac{g^2v^2}{4}},
\qquad
\boxed{m_Z^2=\frac{(g^2+g'^2)v^2}{4}}.
\]

[DERIVED] The terms linear in a CP-even fluctuation and quadratic in gauge fields depend only on

\[
\boxed{v_1\rho_1+v_2\rho_2=v\rho_v},
\qquad
\rho_v\equiv c_\beta\rho_1+s_\beta\rho_2.
\]

Thus the direction parallel to the VEV is the unique CP-even direction with a tree-level linear `VV` coupling.

[DERIVED] The orthogonal direction

\[
\rho_\perp=-s_\beta\rho_1+c_\beta\rho_2
\]

has no tree-level coupling linear in the scalar and quadratic in `W` or `Z`.

[DERIVED] In the Phase-4 mass basis,

\[
\boxed{\rho_v=\sin(\beta-\alpha)h+\cos(\beta-\alpha)\phi}.
\]

Therefore, with the modifier convention

\[
\mathcal L_{SVV}\equiv
\kappa_V^S\left[
\frac{2m_W^2}{v}S W_\mu^+W^{-\mu}
+\frac{m_Z^2}{v}S Z_\mu Z^\mu
\right],
\]

we obtain

\[
\boxed{\kappa_V^h=\sin(\beta-\alpha)},
\qquad
\boxed{\kappa_V^\phi=\cos(\beta-\alpha)}.
\]

On the project exact-alignment branch

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

Hence the project state `h` has exactly the SM tree-level `WW` and `ZZ` couplings, while the project state `phi` has no tree-level `phi WW` or `phi ZZ` vertex in exact alignment.

[DERIVED] No term linear in `eta_i` appears in the gauge-boson mass operator, so the CP-odd state `A` has no tree-level `AVV` vertex of the form `AWW` or `AZZ`.

[DERIVED] The CP-even quartic gauge interaction is angle-independent because

\[
\rho_1^2+\rho_2^2=h^2+\phi^2.
\]

This provides an independent structural cross-check of the orthogonality of the CP-even rotation.

---

## Derivation

### 1. Covariant-derivative convention

We use the modern hypercharge convention

\[
Q=T^3+Y,
\qquad Y(\Phi_i)=\frac12,
\]

so

\[
\boxed{D_\mu=\partial_\mu+i g\,\frac{\sigma^a}{2}W_\mu^a+i g'\frac12 B_\mu}.
\]

Davidson--Haber instead call the doublet hypercharge `Y=1` while using `Q=T3+Y/2`; their covariant derivative therefore contains the same physical factor `g'Y/2=g'/2`. No physics or sign changes under this notation translation.

Define

\[
W_\mu^\pm=\frac{W_\mu^1\mp iW_\mu^2}{\sqrt2},
\]

and

\[
\begin{pmatrix}A_\mu\\Z_\mu\end{pmatrix}
=
\begin{pmatrix}
s_W&c_W\\
c_W&-s_W
\end{pmatrix}
\begin{pmatrix}W_\mu^3\\B_\mu\end{pmatrix},
\]

with

\[
s_W=\frac{g'}{\sqrt{g^2+g'^2}},
\qquad
c_W=\frac{g}{\sqrt{g^2+g'^2}},
\qquad
g_Z\equiv\sqrt{g^2+g'^2}=\frac{g}{c_W}.
\]

Only the neutral lower component is needed to obtain the gauge-boson masses and the neutral CP-even `VV` trilinears. Set temporarily

\[
\Phi_i\longrightarrow
\frac1{\sqrt2}
\begin{pmatrix}0\\x_i\end{pmatrix},
\qquad
x_i\equiv v_i+\rho_i,
\]

with `eta_i=phi_i^+=0` for this extraction.

### 2. Charged gauge part explicitly

Using `T^a=sigma^a/2`,

\[
\left(T^1W^1+T^2W^2\right)
\frac1{\sqrt2}\begin{pmatrix}0\\x_i\end{pmatrix}
=
\begin{pmatrix}x_iW^+/2\\0\end{pmatrix}.
\]

Therefore

\[
D_\mu\Phi_i\supset
\begin{pmatrix}i g x_iW_\mu^+/2\\0\end{pmatrix},
\]

and

\[
(D_\mu\Phi_i)^\dagger(D^\mu\Phi_i)
\supset
\boxed{\frac{g^2}{4}x_i^2W_\mu^+W^{-\mu}}.
\]

Summing the two doublets,

\[
\mathcal L_{\rm kin}\supset
\frac{g^2}{4}
\left[(v_1+\rho_1)^2+(v_2+\rho_2)^2\right]
W_\mu^+W^{-\mu}.
\]

Expanding only to first order in the CP-even fluctuations,

\[
\mathcal L_{\rm kin}\supset
\frac{g^2}{4}(v_1^2+v_2^2)W^+W^-
+\frac{g^2}{2}(v_1\rho_1+v_2\rho_2)W^+W^-+\cdots.
\]

The first term identifies

\[
\boxed{m_W^2=\frac{g^2v^2}{4}}.
\]

The second becomes

\[
\boxed{
\mathcal L_{WW,\,\rm linear}
=\frac{2m_W^2}{v}\rho_v W_\mu^+W^{-\mu}
}.
\]

### 3. Neutral gauge part explicitly and photon cancellation

The lower neutral component has

\[
T^3=-\frac12,
\qquad
Y=+\frac12.
\]

Hence its neutral gauge factor is

\[
-i\frac{g}{2}W_\mu^3+i\frac{g'}{2}B_\mu
=\frac{i}{2}(-gW_\mu^3+g'B_\mu).
\]

Using the `A,Z` definitions,

\[
-gW^3+g'B=-g_Z Z.
\]

The photon cancels exactly because the neutral VEV has `Q=0`. Therefore

\[
D_\mu\Phi_i\supset
\frac1{\sqrt2}
\begin{pmatrix}0\\-i g_Zx_i Z_\mu/2\end{pmatrix},
\]

and

\[
(D_\mu\Phi_i)^\dagger(D^\mu\Phi_i)
\supset
\boxed{\frac{g_Z^2}{8}x_i^2Z_\mu Z^\mu}.
\]

Summing the two doublets,

\[
\mathcal L_{\rm kin}\supset
\frac{g_Z^2}{8}
\left[(v_1+\rho_1)^2+(v_2+\rho_2)^2\right]Z_\mu Z^\mu.
\]

The vacuum term must be written as `(1/2)m_Z^2 Z_mu Z^mu`, giving

\[
\boxed{m_Z^2=\frac{g_Z^2v^2}{4}}.
\]

The linear CP-even term is

\[
\boxed{
\mathcal L_{ZZ,\,\rm linear}
=\frac{m_Z^2}{v}\rho_v Z_\mu Z^\mu
}.
\]

The apparent factor-of-two difference between the displayed `WW` and `ZZ` Lagrangian coefficients is only due to `Z` being a real identical field pair. The corresponding SM Feynman rules are

\[
\boxed{hW^+_\mu W^-_\nu:\quad i\frac{2m_W^2}{v}g_{\mu\nu}},
\]

\[
\boxed{hZ_\mu Z_\nu:\quad i\frac{2m_Z^2}{v}g_{\mu\nu}}.
\]

Thus the normalized modifier is the same for `WW` and `ZZ`.

### 4. Why only the vacuum direction couples linearly

Combining the two sectors,

\[
\mathcal L_{VV,\,\rm linear}
=
\left(v_1\rho_1+v_2\rho_2\right)
\left[
\frac{g^2}{2}W^+W^-+
\frac{g_Z^2}{4}ZZ
\right].
\]

Since

\[
v_1\rho_1+v_2\rho_2
=v(c_\beta\rho_1+s_\beta\rho_2)
=v\rho_v,
\]

we obtain the geometrical result

\[
\boxed{\text{tree-level }SVV\text{ coupling}\propto
\text{projection of }S\text{ onto }\rho_v}.
\]

This statement uses only gauge invariance, canonical kinetic terms and the VEV geometry. It does not use the scalar potential.

### 5. Projection onto the mass eigenstates

From Phase 4,

\[
\rho_1=-s_\alpha h+c_\alpha\phi,
\qquad
\rho_2=c_\alpha h+s_\alpha\phi.
\]

Therefore

\[
\begin{aligned}
\rho_v
&=c_\beta\rho_1+s_\beta\rho_2\\
&=(-c_\beta s_\alpha+s_\beta c_\alpha)h
 +(c_\beta c_\alpha+s_\beta s_\alpha)\phi\\
&=\sin(\beta-\alpha)h+\cos(\beta-\alpha)\phi.
\end{aligned}
\]

Substituting into the gauge interaction gives

\[
\boxed{
\mathcal L_{VV,\,\rm linear}
=
\left[s_{\beta-\alpha}h+c_{\beta-\alpha}\phi\right]
\left[
\frac{2m_W^2}{v}W^+W^-+
\frac{m_Z^2}{v}ZZ
\right]
}.
\]

Thus

\[
\boxed{\kappa_V^h=s_{\beta-\alpha}},
\qquad
\boxed{\kappa_V^\phi=c_{\beta-\alpha}}.
\]

### 6. Exact alignment is now physically identified

Phase 4 established that on the project branch

\[
s_{\beta-\alpha}=1
\quad\Longrightarrow\quad
h=\rho_v,
\qquad
\phi=-\rho_\perp.
\]

The kinetic derivation now supplies the physical meaning:

\[
\boxed{h=\rho_v\Rightarrow hVV\text{ exactly SM-like at tree level}},
\]

\[
\boxed{\phi\perp\rho_v\Rightarrow \phi VV=0\text{ at tree level}}.
\]

This closes the previous logical gap: before Phase 6, `h=rho_v` and `phi=-rho_perp` were purely mixing statements; they are now linked independently to the gauge couplings.

This does **not** yet prove `Z6=0`. That is a separate Higgs-basis statement reserved for Phase 7.

### 7. CP-odd state and absence of `AVV`

If the neutral component is kept complex,

\[
\Phi_i^0=\frac{v_i+\rho_i+i\eta_i}{\sqrt2},
\]

the gauge-boson mass operator depends on

\[
|v_i+\rho_i+i\eta_i|^2
=(v_i+\rho_i)^2+\eta_i^2.
\]

There is no term linear in `eta_i`. Hence no linear combination of `eta_1,eta_2`, including the physical pseudoscalar `A`, has a tree-level `AVV` vertex of this mass-generated form.

### 8. Quartic CP-even gauge couplings as a consistency check

The quadratic CP-even part of the same operator is

\[
\mathcal L\supset
\frac{g^2}{4}(\rho_1^2+\rho_2^2)W^+W^-
+\frac{g_Z^2}{8}(\rho_1^2+\rho_2^2)ZZ.
\]

Because the CP-even rotation is orthogonal,

\[
\boxed{\rho_1^2+\rho_2^2=h^2+\phi^2}.
\]

Thus the diagonal `VVhh` and `VVphiphi` quartics do not depend on `alpha` or `beta`, while the mixed `VVhphi` term cancels. This is exactly the structural behavior expected from a basis rotation of canonical kinetic terms.

---

## Convention map

| Object | Project / DH convention | BFLRS11 | GHOO18 | 2HDMC | Status |
|---|---|---|---|---|---|
| Hypercharge | modern `Y=1/2`, `Q=T3+Y`; equivalent DH `Y=1`, `Q=T3+Y/2` | same physical convention | same EW representation | SM gauge couplings `g,g'` | `[TRANSLATED]` |
| vacuum direction | `rho_v=c_beta rho1+s_beta rho2` | implicit in `hVV/HVV` modifiers | alignment state parallel to VEV | `q_{k1}` first component | `[DERIVED][SOURCE]` |
| `kappa_V^h` | `sin(beta-alpha)` | prose gives same | aligned `H1` is SM-like | `q_11=sba` | `[DERIVED][SOURCE][IMPLEMENTATION-CHECKED]` |
| `kappa_V^phi` | `cos(beta-alpha)` | heavy-state prose gives `cos(alpha-beta)=cos(beta-alpha)` | non-aligned heavy CP-even state has zero `VV` in exact AL | `q_21=cba` | `[DERIVED][SOURCE][IMPLEMENTATION-CHECKED]` |
| exact alignment | `h=rho_v`, `phi=-rho_perp` | coupling-level statement consistent after sign translation | `e1=v,e2=e3=0` | `sba=1,cba=0` | `[DERIVED][TRANSLATED]` |
| `AVV` | zero at tree level | stated zero | CP-odd aligned heavy state has no `VV` | `q_31=0` | `[DERIVED][IMPLEMENTATION-CHECKED]` |

### BFLRS11 sign caution retained

The early simple scalar-state display in BFLRS11 uses global negatives of the DH05 `h,H` fields, while its nearby gauge-coupling prose uses the standard positive modifiers. Phase 6 does not resolve that source-internal presentation by choosing one sentence over another. The project coupling signs come from the kinetic-term derivation in the already-frozen DH field convention.

---

## What was checked against the source

The source comparison was performed only after the analytic result above was frozen.

1. **DH05** states that scalar kinetic terms are canonically gauge-covariant and later summarizes `h_DH WW, h_DH ZZ` as proportional to `sin(beta-alpha)` and `H_DH WW, H_DH ZZ` as proportional to `cos(beta-alpha)`.
2. **BFLRS11** states that the light CP-even state has SM `WW/ZZ` coupling times `sin(beta-alpha)` and the heavy state times `cos(alpha-beta)=cos(beta-alpha)`, and that the pseudoscalar has no `VV` coupling.
3. **GHOO18** defines exact alignment geometrically as a mass eigenstate aligned with the VEV. Its exact-alignment gauge table gives the aligned state the SM `H1ZZ` and `H1WW` vertices and zero `H2VV/H3VV` at exact alignment.
4. **2HDMC implementation check:** `THDM::get_qki` contains
   
   \[
   q_{k1}=(s_{\beta-\alpha},\ c_{\beta-\alpha},\ 0,\ i)
   \]
   
   for `(h,H,A,H+)`, while `THDM::get_coupling_vvh` returns
   
   \[
   c(ZZh_k)=i\,\frac{gM_Z}{c_W}\,\mathrm{Re}(q_{k1}),
   \qquad
   c(WWh_k)=i\,gM_W\,\mathrm{Re}(q_{k1}).
   \]
   
   Since `g M_W=2m_W^2/v` and `g M_Z/c_W=2m_Z^2/v`, this is exactly the derived modifier structure. Thus 2HDMC gives `h: sba`, `H: cba`, `A:0`.

A separate symbolic audit verifies the trigonometric projection, the exact-alignment limit and the angle-independence of the CP-even quartic norm.

---

## What remains uncertain

1. `kappa_V^phi=0` is now **VERIFIED** for exact alignment in the declared project convention.
2. The equivalence between exact alignment and `Z6=0` is **not yet promoted**. Phase 7 must rotate the potential explicitly into the Higgs basis and show that the CP-even off-diagonal mass entry is `Z6 v^2` in the selected Higgs-basis sign convention.
3. The residual sign `H2 -> -H2` in Higgs basis still matters for `Y3`, `Z6`, `Z7` and odd-`H2` scalar couplings. It must be frozen explicitly before any trilinear sign is interpreted.
4. A newly relevant implementation caution is visible in 2HDMC: `get_param_higgs` reports quantities named `Lambda6,Lambda7`, whereas `get_coupling_hhh` subsequently defines local `Z6=-l6`, `Z7=-l7`. This sign layer must be audited in Phase 7/9; no direct identification of printed 2HDMC `Lambda6,Lambda7` with GHOO `Z6,Z7` is allowed yet.
5. Loop-induced `phi -> gamma gamma` or `phi -> Z gamma` amplitudes are not contradicted by `kappa_V^phi=0`; the statement concerns tree-level linear `phi WW/ZZ` couplings only.

---

## Next smallest validation

Proceed to the Higgs-basis construction from the already-frozen generic fields:

\[
\boxed{H_1=c_\beta\Phi_1+s_\beta\Phi_2},
\qquad
\boxed{H_2=-s_\beta\Phi_1+c_\beta\Phi_2}.
\]

The next phase must independently derive, not quote:

1. `H1` carries the full VEV and `H2` has zero VEV;
2. the complete quadratic and quartic Higgs-basis potential;
3. the map `(m_ij^2,lambda_i,beta) -> (Y_i,Z_i)`;
4. the Higgs-basis stationarity relations, including whether `Y3=-Z6 v^2/2` in the selected sign convention;
5. the CP-even Higgs-basis mass matrix and the exact condition connecting its off-diagonal entry to alignment.

Only after those steps should `Z7` or any `phi H+H-` trilinear be interpreted.
