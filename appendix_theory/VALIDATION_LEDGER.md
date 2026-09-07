# Validation ledger — theory appendix

Issue: `#81`  
Rule: `VERIFIED` requires the independent validation appropriate to the claim. A source quotation alone is insufficient when the quantity is reconstructible from first principles.

## F1 — Scalar doublets and VEV normalization

Claim: `Phi_i=(phi_i^+,(v_i+rho_i+i eta_i)/sqrt(2))^T`, with `v^2=v1^2+v2^2` and `tan beta=v2/v1`.  
Status: **VERIFIED**  
Epistemic class: `SOURCE`, `DERIVED`, `TRANSLATED`, `PROJECT-DEFINITION`  
Evidence: electroweak charges/canonical normalization derived in Phase 1; DH05/BFLRS11 source comparison agrees.  
Notes: physical scalar states were not assumed at this stage.

## F2 — Vacuum potential and stationarity

Claim: the `V0(v1,v2)` and both tadpole equations in Phase 3 follow from direct substitution/differentiation of the frozen generic potential.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Evidence: hand derivation, symbolic audit, DH05 comparison and active 2HDMC `m22^2` implementation all agree.  
Notes: stationarity is not global-vacuum proof.

## F3 — Distinction `m22^2`, `m12^2`, `M^2`

Claim: `M^2=m12^2/(s_beta c_beta)` is a derived soft coordinate and is generically distinct from the generic-basis diagonal coefficient `m22^2` and the off-diagonal coefficient `m12^2`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `PROJECT-DEFINITION`, `TRANSLATED`  
Project translation: code/data `M2` means this `M^2`, not `m22_2`.

## F4 — Scalar Hessians and masses

Claim: the charged, CP-odd and CP-even mass matrices in Phase 4 are the second derivatives of the frozen potential at the Phase-3 stationary point.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Evidence: direct Hessians; tadpoles inserted only afterward; symbolic audit; source comparison; active 2HDMC reproduces `m_A^2`, charged–odd splitting and CP-even matrix.  
Notes: positive eigenvalues are local quadratic conditions, not global-minimum proof.

## F5 — Goldstone and physical scalar rotations

Claim: the vacuum direction is the exact zero eigenvector of charged and CP-odd Hessians, giving
`G+=c_beta phi1+ + s_beta phi2+`, `H+=-s_beta phi1+ + c_beta phi2+`, `G0=c_beta eta1+s_beta eta2`, `A=-s_beta eta1+c_beta eta2`; CP-even states use the frozen DH alpha convention.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`, `PROJECT-DEFINITION`  
Evidence: explicit eigenvectors and diagonalization in Phase 4.

## F6 — Gauge masses and VEV-projection theorem for neutral CP-even scalars

Claim: canonical scalar kinetic terms give
`mW^2=g^2 v^2/4`, `mZ^2=(g^2+g'^2)v^2/4`, and every tree-level neutral CP-even `SVV` coupling is proportional to the projection of `S` onto `rho_v=c_beta rho1+s_beta rho2`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Assumptions: canonical electroweak kinetic terms; neutral CP-conserving vacuum; Phase-4 scalar signs.  
Independent derivation: charged and neutral gauge pieces expanded explicitly in Phase 6; photon cancellation follows from `Q=0` of the VEV.  
Source check: DH05/BFLRS11/GHOO18 agree after convention translation.  
Implementation check: 2HDMC `get_qki` and `get_coupling_vvh` implement the same VEV projection.  
Structural checks: no linear `eta_i` term; `rho1^2+rho2^2=h^2+phi^2` makes CP-even diagonal quartic gauge couplings angle-independent.

## C1 — Exact project-compatible generic scalar potential

Claim: the complete CP-conserving project generic-basis potential is the DH05/BFLRS11 operator normalization frozen in Phase 2.  
Status: **VERIFIED**  
Epistemic class: `SOURCE`, `DERIVED`, `TRANSLATED`, `IMPLEMENTATION-CHECKED`  
Notes: campaign `lambda7=0` is not part of the potential definition; GHOO quadratic symbols require the recorded factor/sign map.

## C2 — Physical scalar state convention

Claim: on the project DH-sign branch, `h=h_DH`, `phi=H_DH`; at exact alignment `h=rho_v=c_beta rho1+s_beta rho2` and `phi=-rho_perp=s_beta rho1-c_beta rho2`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`, `PROJECT-DEFINITION`, `IMPLEMENTATION-CHECKED`  
Assumptions: CP conservation and non-degenerate CP-even state identification.  
Notes: state identity is fixed by mixing direction/convention, not solely mass ordering.

## C3 — Type-I exact-alignment fermion modifier

Claim: with
`L_(phi ff)=-(m_f/v) kappa_f^phi phi fbar f`, the project convention gives `kappa_f^phi=-cot(beta)` for `u,d,l`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`, `IMPLEMENTATION-CHECKED`  
Independent derivation: Type I gives `m_f=y_f v2/sqrt2`; `rho2=c_alpha h+s_alpha phi`; hence `kappa_f^phi=s_alpha/s_beta`, which becomes `-cot beta` at `alpha=beta-pi/2`.  
Notes: the sign follows from the previously frozen `phi` field convention.

## C4 — Neutral CP-even gauge modifiers

Claim: `kappa_V^h=sin(beta-alpha)` and `kappa_V^phi=cos(beta-alpha)`, so exact alignment gives `kappa_V^h=1` and `kappa_V^phi=0`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`, `IMPLEMENTATION-CHECKED`  
Independent derivation: Phase 6 expands `sum_i (D_mu Phi_i)^dagger D^mu Phi_i` and proves the linear gauge term depends only on `v1 rho1+v2 rho2=v rho_v`; then `rho_v=s_(beta-alpha)h+c_(beta-alpha)phi`.  
Implementation check: 2HDMC uses `q_{k1}=(sba,cba,0,i)` and `get_coupling_vvh` multiplies `Re(q_{k1})` by the SM `WW/ZZ` vertex.  
Limiting check: exact alignment gives `(1,0)`.  
Notes: this is a tree-level statement; loop-induced `phi->gamma gamma/Z gamma` remain possible.

## C4b — CP-odd tree-level `AVV`

Claim: the CP-odd state has no tree-level linear `AWW` or `AZZ` vertex.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Independent derivation: `|v_i+rho_i+i eta_i|^2=(v_i+rho_i)^2+eta_i^2` contains no term linear in `eta_i`.  
Implementation check: 2HDMC has `q_31=0` in the `VVh_k` projection.

## C5 — Higgs-basis stationarity `Y3=-Z6 v^2/2`

Claim: in the selected Higgs-basis potential convention, `Y3=-Z6 v^2/2`.  
Status: **PARTIALLY VERIFIED**  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Source evidence: DH05 invariant stationarity and GHOO18 agree.  
Missing validation: derive it from the explicit Phase-7 field rotation and frozen `H2` sign.

## C6 — Exact alignment and `Z6`

Claim: exact alignment is equivalent to `Z6=0` under the project assumptions.  
Status: **PARTIALLY VERIFIED**  
Epistemic class: `SOURCE`, `DERIVED`, `OPEN-QUESTION`  
Established so far: Phase 4 identifies the aligned mass direction; Phase 6 independently proves its SM gauge coupling.  
Missing validation: derive the Higgs-basis CP-even mass matrix and identify its off-diagonal entry as `Z6 v^2`, treating degeneracy explicitly.

## C7 — Exact `phi H+H-` trilinear

Claim: exact Higgs-basis expression and exact-alignment reduction in a declared potential/Lagrangian/Feynman-rule convention.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Blocked on: Phase-7 Higgs-basis sign/map and independent cubic expansion.

## C8 — Exact generic-basis expression for `Z7`

Claim: exact analytic map `(lambda_i,beta)->Z7` in the frozen real Higgs-basis convention.  
Status: **PARTIALLY VERIFIED**  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Source evidence: DH05.  
Missing validation: explicit basis rotation/CAS audit.

## C9 — Large-`tan beta`, `lambda7=0` limit of `Z7`

Claim: determine whether `Z7 ~ -lambda6` or another relation holds.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Blocked on C8.

## C10 — Large-`tan beta` approximation for `g_(phi H+H-)`

Claim: determine whether the declared trilinear object satisfies approximately `-v lambda6`.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Blocked on C7–C9.

## C11 — Re-expression with `X=lambda6 tan(beta)`

Claim: determine whether the derived trilinear may be rewritten approximately as `-v X cot(beta)`.  
Status: **NOT VERIFIED**  
Epistemic class: `PROJECT-DEFINITION`, `OPEN-QUESTION`  
Blocked on C10.  
Notes: `X` remains absent from the foundational derivation through Phase 6.

## Implementation caution carried into Phase 7/9

Active 2HDMC `get_param_higgs` returns objects named `Lambda6,Lambda7`, while `get_coupling_hhh` subsequently defines local `Z6=-l6`, `Z7=-l7`. This is not being interpreted yet. Phase 7/9 must establish the exact sign translation before any 2HDMC scalar-trilinear comparison is used as evidence.

## Phase gates

| Gate | Status | Evidence |
|---|---|---|
| `PHASE_0_PASS` | **PASS** | source/convention inventory frozen |
| `PHASE_1_PASS` | **PASS** | field quantum numbers, VEV normalization and beta derived |
| `PHASE_2_PASS` | **PASS** | complete generic potential reconstructed |
| `PHASE_3_PASS` | **PASS** | vacuum potential, tadpoles and `M^2` distinction derived/audited |
| `PHASE_4_PASS` | **PASS** | scalar Hessians, Goldstones, masses and state signs derived/audited |
| `PHASE_5_PASS` | **PASS** | Type-I modifiers derived from the Yukawa Lagrangian |
| `PHASE_6_PASS` | **PASS** | gauge masses/modifiers derived from kinetic terms; `kappa_V^phi=0` and `AVV=0` verified |
| `PHASE_7_PASS` | NOT RUN | explicit Higgs-basis potential derivation required |
| `PHASE_9_PASS` | NOT RUN | exact charged-Higgs trilinear derivation required |
