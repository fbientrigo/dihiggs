# 01 — Convention map

Mission phase: `PHASE 0`  
Rule: equations are not translated into project notation when the source object or sign convention is still ambiguous.

## Paper-by-paper convention table

| Item | DH05 — hep-ph/0504050 | BFLRS11 — 1106.0034 | GHOO18 — 1808.01472 | Phase-0 assessment |
|---|---|---|---|---|
| Doublets | `Phi_1,Phi_2`, complex `SU(2)_L`, `Y=1` | same convention in scalar-sector review | `Phi_1,Phi_2`; general phases `e^{i xi_j}` | `[SOURCE]` compatible after phase choice |
| Neutral expansion | generic VEV first; CP-conserving real basis later | `Phi_a=(phi_a^+,(v_a+rho_a+i eta_a)/sqrt2)^T` | `Phi_j=e^{i xi_j}(varphi_j^+,(v_j+eta_j+i chi_j)/sqrt2)^T` | `[SOURCE]` normalization compatible; names differ |
| `v` | `v^2=v_1^2+v_2^2=(246 GeV)^2` | `v=sqrt(v_1^2+v_2^2)` | same intended definition | `[SOURCE]` compatible |
| `tan beta` | `v_2/v_1`; basis-dependent in general | `v_2/v_1`; phenomenology basis fixed | `v_2/v_1` | `[SOURCE]` same ratio, different discussion of basis meaning |
| Generic quadratic potential | `+m11^2 Phi1†Phi1 +m22^2 Phi2†Phi2 -[m12^2 Phi1†Phi2+h.c.]` | exactly same “notation 1” | `-1/2{m11^2 Phi1†Phi1+m22^2 Phi2†Phi2+[m12^2 Phi1†Phi2+h.c.]}` | `[SOURCE][TRANSLATED]` GHOO18 quadratic symbols are not numerically identical to DH05/BFLRS11 symbols |
| Quartic normalization | `lambda1/2`, `lambda2/2`, `lambda5/2`; `lambda3,4,6,7` as written | same | same quartic normalization | `[SOURCE]` compatible quartic normalization |
| CP assumption | general complex potential first; sec. 4 chooses real CP-conserving basis | multiple sections; simple phenomenology chapter assumes real VEVs; scalar review treats general case | paper treats generic 2HDM67 and alignment; CP-conserving Type-I subcases are identified separately | `[SOURCE]` project CP conservation is a specialization, not a property of all three source treatments |
| CP-even rotation | `h^0=-rho_1 sin alpha+rho_2 cos alpha`; `H^0=rho_1 cos alpha+rho_2 sin alpha` | `h=rho_1 sin alpha-rho_2 cos alpha`; `H=-rho_1 cos alpha-rho_2 sin alpha` | neutral states `H_i=R_ij eta_j`; ordered by mass; in AL `alpha_1=beta, alpha_2=0` | `[SOURCE][TRANSLATED]` BFLRS11's displayed `h,H` are both global-sign reversals of DH05's displayed states; GHOO18 requires a separate state map |
| Mass ordering | notation `h^0,H^0`; conventional light/heavy interpretation but physical mapping uses couplings | explicitly calls `h` lighter and `H` heavier in the simple chapter | `M_1<=M_2<=M_3` by definition | `[SOURCE]` project state identity must not be assigned from ordering alone |
| Goldstone/charged rotation | `G^+=c_beta Phi_1^+ + e^{-i xi}s_beta Phi_2^+`; orthogonal `H^+` | beta diagonalizes charged and CP-odd sectors | explicit orthogonal matrix `(v1/v,v2/v;-v2/v,v1/v)` | `[SOURCE]` compatible at `xi=0` |
| Type-I Yukawa basis | discussed invariantly; Type-I gives a preferred basis | all `u_R,d_R,e_R` couple to `Phi_2` | `eta_1^{u,d,l}=0`, hence `Phi_2` is the Yukawa doublet | `[SOURCE]` same Type-I basis convention |
| Gauge-coupling convention | `h^0 VV ~ sin(beta-alpha)`, `H^0 VV ~ cos(beta-alpha)` in Eq. `littletable` | prose states `h VV ~ sin(beta-alpha)`, `H VV ~ cos(alpha-beta)` | `e_i/v` is the `H_i VV` modifier; exact AL `e_1=v,e_2=e_3=0` | `[SOURCE]` physical coupling pattern compatible, but BFLRS11 field-sign display requires care |
| Higgs-basis fields | `H_1,H_2`, Eq. `higgsbasis`, with residual phase `chi` | `H_1,H_2`, Eq. `2_eq:HBT`; residual `H_2` rephasing | `mathcal H_1, mathcal H_2`; residual `mathcal H_2 -> e^{i chi} mathcal H_2` | `[SOURCE]` same geometric construction up to phase convention |
| Higgs-basis quadratic notation | `M_11^2,M_22^2,-M_12^2`; invariants `Y_1,Y_2,Y_3` with `Y_3=-M_12^2 e^{-2i chi}` | barred `m_ij^2`; same negative off-diagonal convention | `Y_1,Y_2,+Y_3` directly in potential | `[SOURCE][TRANSLATED]` do not identify GHOO18 `Y_3` with a positive `M_12^2` coefficient |
| Higgs-basis quartics | `Lambda_i`; invariant/pseudo-invariant `Z_i`, with `Z_i=Lambda_i` in the real `chi=0` Higgs basis | barred `lambda_i` | `Z_i` directly | `[SOURCE]` `lambda_i`, `Lambda_i`, barred `lambda_i`, and `Z_i` are distinct notation layers until translated |
| Higgs-basis stationarity | `M_11^2=-Lambda_1 v^2/2`, `M_12^2=+Lambda_6 v^2/2`; equivalently invariant `Y_3=-Z_6v^2/2` | `bar m_11^2=-bar lambda_1v^2/2`, `bar m_12^2=+bar lambda_6v^2/2` | `Y_1=-Z_1v^2/2`, `Y_3=-Z_6v^2/2` | `[SOURCE]` source agreement once coefficient definitions are respected; independent derivation still required later |
| `lambda_6,lambda_7` | generic-basis couplings in Eq. `pot`; exact rotation to `Lambda_6,Lambda_7` supplied | same generic-basis meaning in notation 1 | generic `lambda_6,lambda_7`; Higgs-basis `Z_6,Z_7` separately | `[SOURCE]` never silently identify generic `lambda_6/7` with Higgs-basis `Z_6/7` |
| Scalar trilinear object | `g_...` notation used for self-couplings but normalization must be audited before reuse | Appendix writes trilinear **potential** `V_3` explicitly | explicitly states cubic tables are **coefficients of the potential**; Feynman rule needs `-i` and combinatorics | `[SOURCE]` GHOO18 convention is unambiguous; project will separately report `V`, `L_int=-V_int`, and vertex rule |

## Explicit quadratic-parameter translation: DH05/BFLRS11 ↔ GHOO18

Let the symbols on the left denote DH05/BFLRS11 Eq. `pot`/`2_VH1`, and symbols with superscript `(G)` denote GHOO18 Eq. `Eq:pot`.

Matching the same operators gives

\[
 m_{11,\mathrm{DH}}^2=-\frac12 m_{11,(G)}^2,\qquad
 m_{22,\mathrm{DH}}^2=-\frac12 m_{22,(G)}^2,\qquad
 m_{12,\mathrm{DH}}^2=+\frac12 m_{12,(G)}^2.
\]

Epistemic class: `[DERIVED][TRANSLATED]` by term-by-term coefficient matching of the two SOURCE potentials. This is a notation map only; no physical statement follows from the sign of a symbol before the convention is specified.

## Higgs-basis off-diagonal coefficient map

DH05 writes

\[
V\supset-[M_{12}^2 H_1^\dagger H_2+\mathrm{h.c.}],
\]

and Eq. `hbasisinv` gives, in the real `chi=0` Higgs basis,

\[
Y_3=-M_{12}^2,\qquad Z_6=\Lambda_6,\qquad Z_7=\Lambda_7.
\]

GHOO18 writes directly

\[
V\supset+[Y_3\mathcal H_1^\dagger\mathcal H_2+\mathrm{h.c.}].
\]

Therefore the two sources use the same physical Higgs-basis operator with different names for its coefficient. Epistemic class: `[SOURCE][TRANSLATED]`.

## CP-even sign map between DH05 and the simple BFLRS11 phenomenology chapter

Using `rho_i=sqrt(2) Re Phi_i^0-v_i`, DH05 Eq. `scalareigenstates` defines

\[
h_{DH}=-s_\alpha\rho_1+c_\alpha\rho_2,\qquad
H_{DH}=c_\alpha\rho_1+s_\alpha\rho_2.
\]

BFLRS11 displays

\[
h_B=s_\alpha\rho_1-c_\alpha\rho_2,\qquad
H_B=-c_\alpha\rho_1-s_\alpha\rho_2,
\]

so

\[
h_B=-h_{DH},\qquad H_B=-H_{DH}.
\]

Epistemic class: `[DERIVED][TRANSLATED]`. A global sign of a real scalar field is convention freedom, but coupling signs must be transformed consistently. For that reason BFLRS11's quoted modifier signs are not imported directly into project notation in Phase 0.

## Project notation: allowed statements at Phase 0

| Project symbol | Status now | Phase-0 meaning |
|---|---|---|
| `Phi_1,Phi_2` | `[PROJECT-DEFINITION][TRANSLATED]` | Will use the DH05/BFLRS11 generic-basis normalization |
| `v_1,v_2,v` | `[PROJECT-DEFINITION][TRANSLATED]` | Real non-negative VEVs in the CP-conserving basis; `v^2=v_1^2+v_2^2` |
| `tan beta` | `[PROJECT-DEFINITION][TRANSLATED]` | `v_2/v_1` in the fixed Type-I/generic basis |
| `h` | `[PROJECT-DEFINITION]` | SM-like CP-even state near 125.13 GeV; exact source-state map deferred until mass/mixing derivation |
| `phi` | `[PROJECT-DEFINITION]` | additional CP-even state; source-state map deferred until mass/mixing derivation |
| `A,H^±` | `[PROJECT-DEFINITION]` | physical CP-odd and charged states after diagonalization; not assumed before deriving the mass matrices |
| `lambda_6,lambda_7` | `[PROJECT-DEFINITION][TRANSLATED]` | generic-basis quartics in the DH05/BFLRS11 operator normalization |
| `Z_6,Z_7` | `[SOURCE]` | Higgs-basis quantities; no equality to generic `lambda_6,lambda_7` is asserted |
| `X` | intentionally absent | Not introduced in Phase 0 |

## Phase-0 convention gate

`PHASE_0_PASS = YES`.

The next phase may define the two scalar doublets from first principles. It may not yet quote the project `phi f f`, `phi VV`, `phi H^+H^-`, `h phi phi`, or loop-amplitude formulas.
