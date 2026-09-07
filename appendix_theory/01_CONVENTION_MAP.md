# 01 — Convention map

Current state: source inventory from Phase 0 plus convention resolutions through Phase 4.  
Rule: a symbol is translated only after its operator, field sign, and basis are explicit.

## A. Generic-basis scalar potential

The project generic-basis convention is Davidson–Haber (DH05) / Branco notation 1:

\[
\begin{aligned}
V={}&m_{11}^2\Phi_1^\dagger\Phi_1+m_{22}^2\Phi_2^\dagger\Phi_2
-[m_{12}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}]\\
&+\frac{\lambda_1}{2}(\Phi_1^\dagger\Phi_1)^2
+\frac{\lambda_2}{2}(\Phi_2^\dagger\Phi_2)^2
+\lambda_3(\Phi_1^\dagger\Phi_1)(\Phi_2^\dagger\Phi_2)\\
&+\lambda_4(\Phi_1^\dagger\Phi_2)(\Phi_2^\dagger\Phi_1)
+\left\{\frac{\lambda_5}{2}(\Phi_1^\dagger\Phi_2)^2
+[\lambda_6(\Phi_1^\dagger\Phi_1)+\lambda_7(\Phi_2^\dagger\Phi_2)]\Phi_1^\dagger\Phi_2
+\mathrm{h.c.}\right\}.
\end{aligned}
\]

GHOO18 instead writes the quadratic sector as

\[
-\frac12\{m_{11,G}^2\Phi_1^\dagger\Phi_1+m_{22,G}^2\Phi_2^\dagger\Phi_2+[m_{12,G}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}]\}.
\]

Term-by-term matching gives

\[
\boxed{
m_{11,\rm DH}^2=-\frac12m_{11,G}^2,
\qquad
m_{22,\rm DH}^2=-\frac12m_{22,G}^2,
\qquad
m_{12,\rm DH}^2=+\frac12m_{12,G}^2.}
\]

The quartic normalization is the same once the same generic basis is identified.

## B. Fields, VEVs, and letters

| Item | DH05 | BFLRS11 | GHOO18 | Project |
|---|---|---|---|---|
| Doublets | `Phi1,Phi2`, `Y=1` | same | same representation | `Phi1,Phi2` |
| CP-even neutral fluctuation | `sqrt(2) Re Phi_i^0-v_i` | `rho_i` | `eta_i` | `rho_i` |
| CP-odd neutral fluctuation | imaginary neutral component | `eta_i` | `chi_i` | current derivation files use `eta_i`; final writing should prefer an unambiguous project letter such as `a_i` |
| VEVs | `v1`, `v2 e^{i xi}` generally | real in simple CP-conserving treatment | phased `v_i` generally | real non-negative `v1,v2` |
| `tan beta` | `v2/v1` | `v2/v1` | `v2/v1` | `v2/v1` in the fixed project/Yukawa basis |

Hard collision: `eta_i` is CP-odd in BFLRS11 and CP-even in GHOO18. Any cross-source transcription must rename fields before comparing formulas.

## C. Type-I Yukawa basis — corrected map

This is a real convention difference and must not be suppressed.

- **DH05:** its explicitly chosen Type-I preferred basis is the one with the second Yukawa matrices zero, so in that displayed choice **Phi1 couples to fermions**. DH05 also notes that the exchanged basis `Phi1 <-> Phi2` is physically equivalent.
- **BFLRS11:** Type I is written with **Phi2 coupling to all right-handed fermions**.
- **GHOO18:** `eta_1^{u,d,l,0}=0`, hence **Phi2 couples**.
- **Project:** uses the BFLRS11/GHOO18 assignment, **Phi2 couples**, with `tan beta=v2/v1`.

Therefore a DH05 Type-I statement cannot be imported at fixed symbol `tan beta` without the doublet exchange; in the exchanged convention the corresponding ratio is mapped by `tan beta <-> cot beta`.

This corrects the earlier Phase-0 shorthand that treated all three sources as using the same Type-I assignment.

## D. Charged and CP-odd rotations — resolved in Phase 4

The independently derived Hessians force

\[
R_\beta=\begin{pmatrix}c_\beta&s_\beta\\-s_\beta&c_\beta\end{pmatrix}.
\]

Project convention:

\[
\boxed{
G^+=c_\beta\phi_1^++s_\beta\phi_2^+,
\qquad
H^+=-s_\beta\phi_1^++c_\beta\phi_2^+}
\]

and

\[
\boxed{
G^0=c_\beta\eta_1+s_\beta\eta_2,
\qquad
A=-s_\beta\eta_1+c_\beta\eta_2.}
\]

These directions agree with the real-phase DH05 and GHOO18 rotations. Their origin in the project is the Hessian zero mode, not a copied convention.

## E. CP-even rotation and project state signs — resolved in Phase 4

DH05 defines

\[
h_{\rm DH}=-s_\alpha\rho_1+c_\alpha\rho_2,
\qquad
H_{\rm DH}=c_\alpha\rho_1+s_\alpha\rho_2.
\]

The project freezes

\[
\boxed{h\equiv h_{\rm DH},\qquad \phi\equiv H_{\rm DH}}
\]

on the branch continuously connected to

\[
\sin(\beta-\alpha)=+1.
\]

This is a field/mixing convention, not a mass-ordering rule. The project does not inherit DH05's conventional `m_H>=m_h` naming.

At exact alignment,

\[
\boxed{h=c_\beta\rho_1+s_\beta\rho_2},
\qquad
\boxed{\phi=s_\beta\rho_1-c_\beta\rho_2}.
\]

Thus `h` is the CP-even VEV direction and `phi` the orthogonal project state. The statement that these overlaps equal normalized `VV` couplings remains reserved for the kinetic-term derivation in Phase 6.

### BFLRS11 sign map

The early pedagogical display uses

\[
h_B=s_\alpha\rho_1-c_\alpha\rho_2=-h,
\qquad
H_B=-c_\alpha\rho_1-s_\alpha\rho_2=-\phi.
\]

Nearby BFLRS11 prose/Yukawa tables use the standard opposite sign convention. Therefore those displayed field definitions are not used as project sign authority.

### GHOO18 sign map

In the CP-conserving exact-alignment branch with `alpha3=0`, GHOO18 has

\[
H_2^{\rm GHOO}=-s_\beta\rho_1+c_\beta\rho_2.
\]

Hence

\[
\boxed{H_2^{\rm GHOO}=-\phi_{\rm project}.}
\]

This is the single field redefinition behind the correlated sign flips later seen between GHOO18 and DH-project conventions in both the non-SM Type-I Yukawa coupling and the `phi H+ H-` potential coefficient. Those signs must always be translated together.

## F. Higgs basis

In the real CP-conserving limit all three sources support the geometric rotation

\[
\boxed{
H_1=c_\beta\Phi_1+s_\beta\Phi_2,
\qquad
H_2=-s_\beta\Phi_1+c_\beta\Phi_2.}
\]

DH05/S2 write the mixed quadratic coefficient with an explicit minus,

\[
V\supset-[M_{12}^2H_1^\dagger H_2+\mathrm{h.c.}],
\]

whereas GHOO18 writes

\[
V\supset+[Y_3H_1^\dagger H_2+\mathrm{h.c.}].
\]

In the real `chi=0` convention,

\[
\boxed{Y_3=-M_{12}^2,\qquad Z_i=\Lambda_i.}
\]

The source stationarity statements

\[
M_{12}^2=+\frac12\Lambda_6v^2
\]

and

\[
Y_3=-\frac12Z_6v^2
\]

are therefore the same convention-translated equation. Independent derivation remains Phase 7.

At exact alignment, the CP-even real part of the above `H2` is

\[
\rho_{H_2}=-s_\beta\rho_1+c_\beta\rho_2=-\phi_{\rm project}.
\]

So when Phase 7 freezes this exact Higgs-basis rotation, the residual real sign is already constrained by the Phase-4 project state convention.

## G. Scalar trilinear bookkeeping

- GHOO18 explicitly quotes cubic quantities as **coefficients in the potential**.
- BFLRS11 writes the trilinear potential `V3` explicitly.
- DH05 scalar `g_...` normalization is not imported without an independent convention audit.

The project rule is always

\[
\text{coefficient in }V
\quad\to\quad
\text{coefficient in }\mathcal L_{\rm int}=-V_{\rm int}
\quad\to\quad
\text{Feynman rule}.
\]

No bare symbol `g` is accepted without saying which of these three objects it denotes.

## H. Symbols that must never be silently identified

| Pair | Required distinction |
|---|---|
| `lambda6` vs `Z6` | generic-basis quartic vs Higgs-basis quartic |
| `lambda6` vs `Z7` | different generic/Higgs-basis objects |
| `lambda7` vs `Z7` | generic-basis quartic vs rotated Higgs-basis combination |
| `Lambda7` vs `lambda7` | DH05 Higgs-basis vs generic-basis quartic |
| `M^2` vs `m12^2` | `M^2=m12^2/(s_beta c_beta)` in the project DH convention |
| `M^2` vs `m22^2` | derived soft coordinate vs diagonal quadratic coefficient |
| S1 `M_ij^2` vs physical `M_i^2` | Higgs-basis potential coefficients vs squared masses |
| four-index `Z_abcd` vs one-index `Z_i` | covariant quartic tensor vs Higgs-basis quartics |

## I. Frozen project choices through Phase 4

1. Generic scalar potential: DH05/BFLRS11 notation 1.
2. Real CP-conserving vacuum: `v1,v2>=0`, `tan beta=v2/v1`.
3. Type I: `Phi2` couples to all fermion species.
4. Charged/odd physical directions: `H+ = -s beta phi1+ + c beta phi2+`, `A=-s beta eta1+c beta eta2`.
5. CP-even signs: `h=h_DH`, `phi=H_DH`.
6. Exact project alignment branch: `sin(beta-alpha)=+1`, so `h=rho_v`, `phi=-rho_perp`.
7. No state is identified merely by mass ordering.
8. Higgs-basis notation later: `Y_i,Z_i` with `+[Y3 H1†H2+h.c.]`.
9. Trilinears: coefficient in `V` first, then `L_int`, then vertex.
10. `X` remains absent until the foundational coupling chain is derived.

`PHASE_4_CONVENTION_MAP = LOCKED FOR PHASES 5–6`.
