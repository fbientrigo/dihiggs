# Validation ledger — theory appendix

Issue: `#81`  
Rule: `VERIFIED` requires the independent validation appropriate to the claim. A source quotation alone is insufficient when the quantity is reconstructible from first principles.

## F1 — Scalar doublets and VEV normalization

Claim: `Phi_i=(phi_i^+,(v_i+rho_i+i eta_i)/sqrt(2))^T`, with `v^2=v1^2+v2^2` and `tan beta=v2/v1`.  
Status: **VERIFIED**  
Epistemic class: `SOURCE`, `DERIVED`, `TRANSLATED`, `PROJECT-DEFINITION`  
Evidence: electroweak charges/canonical normalization derived in Phase 1; DH05/BFLRS11 source comparison agrees.

## F2 — Vacuum potential and stationarity

Claim: `V0(v1,v2)` and both generic-basis tadpole equations follow from direct substitution/differentiation of the frozen potential.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Evidence: hand derivation, symbolic audit, DH05 comparison and active 2HDMC `m22^2` implementation all agree.  
Notes: stationarity is not global-vacuum proof.

## F3 — Distinction `m22^2`, `m12^2`, `M^2`

Claim: `M^2=m12^2/(s_beta c_beta)` is a derived soft coordinate and is generically distinct from `m22^2` and `m12^2`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `PROJECT-DEFINITION`, `TRANSLATED`  
Project translation: code/data `M2` means this `M^2`, not `m22_2`.

## F4 — Scalar Hessians and masses

Claim: charged, CP-odd and CP-even mass matrices are the second derivatives of the frozen generic potential evaluated at the Phase-3 stationary point.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Evidence: direct Hessians; tadpoles inserted only afterward; symbolic audit; source comparison; active 2HDMC reproduces `m_A^2`, charged–odd splitting and CP-even matrix.

## F5 — Goldstone and physical scalar rotations

Claim: the vacuum direction is the exact zero eigenvector of charged and CP-odd Hessians, giving `G+,H+,G0,A`; CP-even states use the frozen DH alpha convention.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`, `PROJECT-DEFINITION`

## F6 — Gauge masses and VEV-projection theorem

Claim: canonical scalar kinetic terms give `mW^2=g^2v^2/4`, `mZ^2=(g^2+g'^2)v^2/4`, and every tree-level neutral CP-even `SVV` coupling is proportional to the projection of `S` onto `rho_v=c_beta rho1+s_beta rho2`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`

## F7 — Higgs-basis field rotation and operator map

Claim: with the frozen project sign
`H1=c_beta Phi1+s_beta Phi2`, `H2=-s_beta Phi1+c_beta Phi2`, one has `<H1^0>=v/sqrt2`, `<H2^0>=0`, and the complete generic potential transforms into the declared `Y_i,Z_i` Higgs-basis potential.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`, `IMPLEMENTATION-CHECKED`, `PROJECT-DEFINITION`  
Independent derivation: Phase 7 reconstructs the four generic bilinears in terms of `H1,H2` bilinears and extracts all quadratic/quartic coefficients operator by operator.  
Symbolic audit: `phase7_higgs_basis_check.py` automatically reconstructs `Y1,Y2,Y3,Z1...Z7`.  
Source check: DH05 `higgsbasis`, `maa`–`mab`, `Lam1def`–`Lam7def`.  
Implementation check: active 2HDMC `get_param_higgs` matches the project `Z1...Z7` formulas, including the signs of `Z6,Z7`.

## C1 — Exact project-compatible generic scalar potential

Claim: the complete CP-conserving generic-basis potential is the DH05/BFLRS11 operator normalization frozen in Phase 2.  
Status: **VERIFIED**  
Epistemic class: `SOURCE`, `DERIVED`, `TRANSLATED`, `IMPLEMENTATION-CHECKED`

## C2 — Physical scalar state convention

Claim: on the project DH-sign branch, `h=h_DH`, `phi=H_DH`; at exact alignment `h=rho_v` and `phi=-rho_perp`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`, `PROJECT-DEFINITION`, `IMPLEMENTATION-CHECKED`

## C3 — Type-I exact-alignment fermion modifier

Claim: with `L_(phi ff)=-(m_f/v) kappa_f^phi phi fbar f`, the project convention gives `kappa_f^phi=-cot(beta)` for `u,d,l`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`, `IMPLEMENTATION-CHECKED`

## C4 — Neutral CP-even gauge modifiers

Claim: `kappa_V^h=sin(beta-alpha)` and `kappa_V^phi=cos(beta-alpha)`, so exact alignment gives `(1,0)`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`, `IMPLEMENTATION-CHECKED`

## C4b — CP-odd tree-level `AVV`

Claim: the CP-odd state has no tree-level linear `AWW` or `AZZ` vertex.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`

## C5 — Higgs-basis stationarity `Y3=-Z6 v^2/2`

Claim: in the selected Higgs-basis convention with `+[Y3 H1dag H2+h.c.]`,
`Y1=-Z1 v^2/2` and `Y3=-Z6 v^2/2`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`, `IMPLEMENTATION-CHECKED`  
Independent derivation A: direct Higgs-basis tadpoles in Phase 7.  
Independent derivation B: substitute Phase-3 generic tadpole solutions into the derived `Y1,Y3,Z1,Z6`; both identities vanish exactly.  
Translation: DH05 uses `-[M12_H^2 H1dag H2+h.c.]`, so `Y3=-M12_H^2`; DH stationarity `M12_H^2=+Lambda6 v^2/2` becomes the project equation.

## C6 — Exact alignment and `Z6`

Claim: exact tree-level alignment is equivalent to `Z6=0` when alignment means the VEV direction `rho_v` is a CP-even mass eigenstate.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Independent derivation: the Phase-7 CP-even Higgs-basis Hessian is
`[[Z1 v^2, Z6 v^2],[Z6 v^2, m_A^2+Z5 v^2]]`; the VEV vector `(1,0)` is an eigenvector iff `Z6=0`.  
Physical identity: `Z6 v^2=(m_h^2-m_phi^2)sba*cba`.  
Degeneracy note: exact degeneracy removes uniqueness of the mixing-angle label, not the off-diagonal condition.

## C7 — Exact `phi H+H-` trilinear

Claim: exact Higgs-basis expression and exact-alignment reduction in a declared potential/Lagrangian/Feynman-rule convention.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Prerequisites now resolved: Higgs-basis sign and `Z7` definition.  
Remaining validation: explicit cubic operator extraction, physical-state substitution and implementation sign translation.

## C8 — Exact generic-basis expression for `Z7`

Claim: exact analytic map `(lambda_i,beta)->Z7` in the frozen real Higgs-basis convention.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Independent derivation: operator expansion in Phase 7; symbolic coefficient extraction.  
Source check: DH05 `Lam7def`.  
Implementation check: 2HDMC `get_param_higgs` returns the same expression and sign.

## C9 — Large-`tan beta`, `lambda7=0` limit of `Z7`

Claim: for the frozen project `H2` sign,
`Z7 -> -lambda6` as `tan beta -> infinity` with `lambda7=0`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`  
Exact finite-`t` form:
`Z7=-[lambda6 t^4+(lambda1-lambda345)t^3-3 lambda6 t^2+(lambda345-lambda2)t]/(1+t^2)^2`.  
Expansion: `Z7=-lambda6+(lambda345-lambda1)/t+O(t^-2)`.  
Notes: this does not yet determine the physical `phi H+H-` sign because exact alignment has `phi=-rho_perp`.

## C10 — Large-`tan beta` approximation for the physical `phi H+H-` trilinear

Claim: determine whether the declared trilinear object satisfies approximately `-v lambda6`.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Blocked only on C7/trilinear convention audit; C8–C9 are now closed.

## C11 — Re-expression with `X=lambda6 tan(beta)`

Claim: determine whether the derived trilinear may be rewritten approximately as `-v X cot(beta)`.  
Status: **NOT VERIFIED**  
Epistemic class: `PROJECT-DEFINITION`, `OPEN-QUESTION`  
Blocked on C7/C10.  
Notes: `X` remains absent from the foundational derivation through Phase 7.

## Implementation caution carried into the trilinear phase

Phase 7 establishes that `THDM::get_param_higgs` returns `Lambda6,Lambda7` with the **same signs** as the project/DH `Z6,Z7` under the frozen `H2=-s Phi1+c Phi2` convention. However, `THDM::get_coupling_hhh` subsequently defines local `Z6=-l6`, `Z7=-l7`. Therefore that extra minus sign is not part of the basis transformation. It is a local trilinear/Feynman-rule convention layer that must be mapped operator by operator before implementation-level sign comparisons are used.

## Phase gates

| Gate | Status | Evidence |
|---|---|---|
| `PHASE_0_PASS` | **PASS** | source/convention inventory frozen |
| `PHASE_1_PASS` | **PASS** | field quantum numbers, VEV normalization and beta derived |
| `PHASE_2_PASS` | **PASS** | complete generic potential reconstructed |
| `PHASE_3_PASS` | **PASS** | vacuum potential, tadpoles and `M^2` distinction derived/audited |
| `PHASE_4_PASS` | **PASS** | scalar Hessians, Goldstones, masses and state signs derived/audited |
| `PHASE_5_PASS` | **PASS** | Type-I modifiers derived from the Yukawa Lagrangian |
| `PHASE_6_PASS` | **PASS** | gauge masses/modifiers derived from kinetic terms |
| `PHASE_7_PASS` | **PASS** | full Higgs-basis rotation, `Y_i,Z_i`, stationarity, mass matrix and alignment criterion derived/audited |
| `PHASE_9_PASS` | NOT RUN | exact charged-Higgs trilinear derivation required |
