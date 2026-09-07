# 02 — Scalar fields and generic CP-conserving potential

Scope: Phases 1–2 of issue #81.  
Analytic anchor: Davidson–Haber (DH05), `hbasis.tex`, Eq. `pot` and Eq. `potmin`.  
Cross-source checks: Branco et al. (BFLRS11), `PhysRep_large.tex`, Eq. `2_VH1`; GHOO18, `paper_heavyhiggs_jhep_revised3.tex`, Eq. `Eq:pot` and Eq. `vevs`.

The project convention is chosen to coincide with the DH05/BFLRS11 generic-basis convention. GHOO18 is retained as an explicit translation layer because its quadratic mass-parameter normalization differs.

---

# Phase 1 — Start from the two scalar doublets

## What is established

[SOURCE][TRANSLATED] The scalar sector contains two complex scalar doublets with identical electroweak quantum numbers. In the Davidson–Haber hypercharge normalization,

\[
\Phi_i \sim (\mathbf 2,Y=1),\qquad i=1,2,
\]

with

\[
Q=T_3+\frac{Y}{2}.
\]

Equivalently, in the common convention `Q=T_3+Y`, the same doublets carry `Y=+1/2`. These are two hypercharge normalizations, not different models.

[PROJECT-DEFINITION][TRANSLATED] For the CP-conserving project basis,

\[
\boxed{
\Phi_i=
\begin{pmatrix}
\phi_i^+\\[2mm]
\dfrac{v_i+\rho_i+i\eta_i}{\sqrt2}
\end{pmatrix}
},\qquad i=1,2,
\]

with real `rho_i, eta_i` and real non-negative VEVs `v_i`. Hence

\[
\langle\Phi_i\rangle=
\frac1{\sqrt2}
\begin{pmatrix}0\\v_i\end{pmatrix},
\qquad
\boxed{v^2=v_1^2+v_2^2},
\qquad
\boxed{\tan\beta=\frac{v_2}{v_1}},
\]

so

\[
v_1=v\cos\beta,\qquad v_2=v\sin\beta,
\]

with `0 <= beta <= pi/2` when `v_1,v_2 >= 0`.

No physical `h`, `phi`, `A`, or `H^+` state is defined in this phase.

## Derivation

For `Y=1` and `Q=T_3+Y/2`,

\[
Q_{\rm upper}=+\frac12+\frac12=+1,
\qquad
Q_{\rm lower}=-\frac12+\frac12=0.
\]

An electromagnetic-preserving vacuum therefore lies in the neutral lower components. DH05 explicitly uses an electroweak gauge rotation to place the VEVs there.

The neutral field decomposition

\[
\Phi_i^0=\frac{1}{\sqrt2}(v_i+\rho_i+i\eta_i)
\]

ensures canonical normalization because

\[
|\partial_\mu\Phi_i^0|^2
=\frac12(\partial_\mu\rho_i)^2
+\frac12(\partial_\mu\eta_i)^2.
\]

The existence of two equal-hypercharge doublets is model content. The neutral vacuum orientation is gauge choice. Taking both VEVs real and non-negative is a basis/phase choice appropriate to the CP-conserving vacuum. `tan(beta)` is basis-dependent in a fully generic 2HDM and becomes a meaningful project coordinate only after the generic/Yukawa basis is fixed.

## Convention map

| Object | DH05 | BFLRS11 | GHOO18 | Project choice |
|---|---|---|---|---|
| doublet | `Phi_i`, `Y=1` | `Phi_i` | `Phi_j` | `Phi_i` |
| CP-even neutral fluctuation | `sqrt(2) Re Phi_i^0-v_i` | `rho_i` | `eta_i` | `rho_i` |
| CP-odd neutral fluctuation | imaginary neutral component | `eta_i` | `chi_i` | `eta_i` |
| VEV phase | `v_2 e^{i xi}` generally | real in simple CP-conserving chapter | `e^{i xi_j}` generally | `xi=0`, real VEVs |
| `tan beta` | `v_2/v_1` | `v_2/v_1` | `v_2/v_1` | `v_2/v_1` |

## What was checked against the source

- DH05 `hbasis.tex` lines 2189–2190: two complex `Y=1`, `SU(2)_L` scalar doublets.
- DH05 Eq. `potmin`, source lines 2225–2235: neutral VEVs and `v^2=v_1^2+v_2^2`.
- DH05 Eq. `tanbdef`, source lines 2252–2259: `s_beta=v_2/v`, `c_beta=v_1/v`, `tan beta=v_2/v_1`.
- BFLRS11 source lines 547–557: real non-negative VEV convention and explicit field decomposition.
- GHOO18 Eq. `vevs`, source lines 264–274: general phased doublet decomposition.

## What remains uncertain

Nothing in the field normalization blocks Phase 2. The physical state map remains intentionally undefined until quadratic terms are derived.

## Next smallest validation

Construct every independent gauge-invariant operator of mass dimension <=4 from `Phi_i^dagger Phi_j`, impose Hermiticity, and then specialize the general potential to a real CP-conserving basis.

---

# Phase 2 — Generic CP-conserving scalar potential

## What is established

[SOURCE][TRANSLATED] In the selected DH05/BFLRS11 convention,

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
&+\left\{\frac{\lambda_5}{2}(\Phi_1^\dagger\Phi_2)^2
+\left[\lambda_6(\Phi_1^\dagger\Phi_1)
+\lambda_7(\Phi_2^\dagger\Phi_2)\right]\Phi_1^\dagger\Phi_2
+\mathrm{h.c.}\right\}.
\end{aligned}}
\]

For the CP-conserving project basis,

\[
\boxed{m_{12}^2,\lambda_5,\lambda_6,\lambda_7\in\mathbb R},
\]

in addition to `m11^2,m22^2,lambda1...lambda4` being real by Hermiticity. The campaign choice `lambda7=0` is not imposed here.

## Derivation

Each scalar has mass dimension one, so the gauge-singlet bilinears

\[
B_{ij}=\Phi_i^\dagger\Phi_j
\]

have dimension two. Renormalizability permits one bilinear or a product of two bilinears.

The quadratic operator basis is

\[
B_{11},\qquad B_{22},\qquad B_{12}+B_{21}.
\]

The minus sign multiplying `m12^2` is not fixed by gauge symmetry; it is the selected DH05/BFLRS11/2HDMC convention.

A complete convenient quartic basis is

\[
B_{11}^2,\ B_{22}^2,\ B_{11}B_{22},\ B_{12}B_{21},\
B_{12}^2+\mathrm{h.c.},\ B_{11}B_{12}+\mathrm{h.c.},\ B_{22}B_{12}+\mathrm{h.c.}
\]

with coefficients normalized as in the boxed potential.

Hermiticity requires

\[
m_{11}^2,m_{22}^2,\lambda_{1,2,3,4}\in\mathbb R,
\]

while `m12^2,lambda5,lambda6,lambda7` may be complex in the general model because their operators are accompanied by Hermitian conjugates. The project specializes to a basis where these parameters and the VEVs are real, excluding a spontaneously CP-violating vacuum from the working branch.

| Term | Operator dim. | Coefficient dim. | Hermiticity/CP |
|---|---:|---:|---|
| `Phi1†Phi1`, `Phi2†Phi2` | 2 | 2 | self-adjoint |
| `Phi1†Phi2` | 2 | 2 | add h.c.; real coefficient in project CP basis |
| `lambda1,lambda2` terms | 4 | 0 | self-adjoint |
| `lambda3` term | 4 | 0 | self-adjoint |
| `lambda4` term | 4 | 0 | self-adjoint |
| `lambda5` term | 4 | 0 | add h.c.; real in project CP basis |
| `lambda6` term | 4 | 0 | add h.c.; real in project CP basis |
| `lambda7` term | 4 | 0 | add h.c.; real in project CP basis |

## Convention map

GHOO18 writes

\[
V_G\supset-\frac12\left\{
m_{11,G}^2\Phi_1^\dagger\Phi_1
+m_{22,G}^2\Phi_2^\dagger\Phi_2
+[m_{12,G}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}]
\right\}.
\]

Matching the coefficients of the same operators gives

\[
\boxed{
m_{11,\rm DH}^2=-\frac12m_{11,G}^2,
\qquad
m_{22,\rm DH}^2=-\frac12m_{22,G}^2,
\qquad
m_{12,\rm DH}^2=+\frac12m_{12,G}^2.
}
\]

The displayed quartic normalization is the same,

\[
\lambda_{i,\rm DH}=\lambda_{i,G},\qquad i=1,\ldots,7,
\]

once the same field basis is identified. This is `[DERIVED][TRANSLATED]` coefficient matching, not a physical sign statement.

## What was checked against the source

- DH05 Eq. `pot`, source lines 2193–2210: complete generic potential and reality properties.
- BFLRS11 Eq. `2_VH1`, source lines 6252–6295: identical convention and explicit statement that it follows Davidson–Haber definitions.
- GHOO18 Eq. `Eq:pot`, source lines 244–258: same quartic normalization but different quadratic normalization.
- The active vendored 2HDMC `THDM::set_param_gen` accepts `(lambda1,...,lambda7,m12_2,tan_beta)`; the actual stationarity formula implemented for `m22_2` is checked independently in `03_VACUUM_AND_MINIMIZATION.md`.

## What remains uncertain

No sign/factor ambiguity remains in the selected generic CP-conserving potential. Whether every downstream 2HDMC helper preserves the convention is an implementation audit question, not a blocker for the theoretical definition.

## Next smallest validation

Insert the real neutral VEVs into the selected potential, evaluate each invariant separately, and differentiate `V(v1,v2)` without using a pre-tabulated tadpole formula.

`PHASE_2_PASS = PASS` for the analytic potential convention.
