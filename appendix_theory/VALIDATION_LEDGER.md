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

Claim: with the frozen project sign `H1=c_beta Phi1+s_beta Phi2`, `H2=-s_beta Phi1+c_beta Phi2`, one has `<H1^0>=v/sqrt2`, `<H2^0>=0`, and the complete generic potential transforms into the declared `Y_i,Z_i` Higgs-basis potential.  
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

Claim: in the selected Higgs-basis convention with `+[Y3 H1dag H2+h.c.]`, `Y1=-Z1 v^2/2` and `Y3=-Z6 v^2/2`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`, `IMPLEMENTATION-CHECKED`  
Independent derivation A: direct Higgs-basis tadpoles in Phase 7.  
Independent derivation B: substitute Phase-3 generic tadpole solutions into the derived `Y1,Y3,Z1,Z6`; both identities vanish exactly.  
Translation: DH05 uses `-[M12_H^2 H1dag H2+h.c.]`, so `Y3=-M12_H^2`; DH stationarity `M12_H^2=+Lambda6 v^2/2` becomes the project equation.

## C6 — Exact alignment and `Z6`

Claim: exact tree-level alignment is equivalent to `Z6=0` when alignment means the VEV direction `rho_v` is a CP-even mass eigenstate.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Independent derivation: the Phase-7 CP-even Higgs-basis Hessian is `[[Z1 v^2, Z6 v^2],[Z6 v^2, m_A^2+Z5 v^2]]`; the VEV vector `(1,0)` is an eigenvector iff `Z6=0`.  
Physical identity: `Z6 v^2=(m_h^2-m_phi^2)sba*cba`.  
Degeneracy note: exact degeneracy removes uniqueness of the mixing-angle label, not the off-diagonal condition.

## C7 — Exact `phi H+H-` trilinear

Claim: the exact potential coefficient is `C_V=v(Z3 cba-Z7 sba)` and at alignment `C_V=-vZ7`; the literal coefficient in `L_int=-V_int` is `+vZ7`; with `L_int=-g phi H+H-`, `g=-vZ7`; the Feynman rule is `-i C_V=+i vZ7` at alignment.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Independent derivation: direct cubic extraction from the Phase-7 Higgs-basis potential and symbolic third derivatives.  
Source check: GHOO18 defines `q_i` as the coefficient in the potential and prescribes multiplication by `-i`, with combinatorics only for identical fields.  
Implementation check: active 2HDMC reproduces `-i v(Z3 cba-Z7 sba)` after combining `get_qki` with its local `Z7=-l7` convention.

## C8 — Exact generic-basis expression for `Z7`

Claim: exact analytic map `(lambda_i,beta)->Z7` in the frozen real Higgs-basis convention.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Independent derivation: operator expansion in Phase 7; symbolic coefficient extraction.  
Source check: DH05 `Lam7def`.  
Implementation check: 2HDMC `get_param_higgs` returns the same expression and sign.

## C9 — Large-`tan beta`, `lambda7=0` limit of `Z7`

Claim: for the frozen project `H2` sign, `Z7 -> -lambda6` as `tan beta -> infinity` with `lambda7=0`.  
Status: **VERIFIED AT FIXED `lambda6`**  
Epistemic class: `DERIVED`  
Exact finite-`t` form: `Z7=-[lambda6 t^4+(lambda1-lambda345)t^3-3 lambda6 t^2+(lambda345-lambda2)t]/(1+t^2)^2`.  
Expansion: `Z7=-lambda6+(lambda345-lambda1)/t+O(t^-2)`.  
Notes: Phase 10 shows that this asymptotic statement cannot be converted without qualification into a fixed-`X=lambda6 t` statement.

## C10 — Large-`tan beta` approximation for the physical `phi H+H-` trilinear

Claim: object-dependent large-`tan beta` limit at fixed `lambda6`, `lambda7=0`.  
Status: **VERIFIED WITH OBJECT LABELS AT FIXED `lambda6`**  
Epistemic class: `DERIVED`  
Results: `C_V=-vZ7 -> +v lambda6`; literal `C_L=+vZ7 -> -v lambda6`; if `L_int=-g phiH+H-`, then `g=C_V -> +v lambda6`; Feynman rule `-> -i v lambda6`.  
Notes: the old unqualified statement `g~-v lambda6` is convention-ambiguous and must not be used without declaring which object is meant.

## C11 — Re-expression with `X=lambda6 tan(beta)`

Claim: determine the leading large-`tan beta` trilinear when `X=lambda6 tan(beta)` itself is held fixed.  
Status: **CONDITIONAL; PREVIOUS UNQUALIFIED FORM DOWNGRADED**  
Epistemic class: `PROJECT-DEFINITION`, `DERIVED`, `OPEN-QUESTION`  
Exact observation: the fixed-`lambda6` limit `Z7->-lambda6` and the fixed-`X` limit are different asymptotic procedures. With `lambda7=0`, exact alignment and fixed `X`, `Z7=-(X+lambda1-lambda2) cot(beta)+O(cot^3 beta)`.  
Therefore `Z7~-X cot(beta)` additionally requires `|lambda1-lambda2| << |X|` or an equivalent numerical cancellation.  
Consequences: at fixed `X`, `C_V=-vZ7 ~ v(X+lambda1-lambda2)cot(beta)` and the literal `C_L=+vZ7` carries the opposite sign.  
Notes: use exact `Z7` point-by-point unless the extra hierarchy is explicitly verified.

## C12 — Charged-scalar loop convention in `phi -> gamma gamma`

Claim: the active project/2HDMC path consumes the Feynman-rule object returned by `get_coupling_hhh`, and in a standard reduced-amplitude convention the charged-scalar term is `C_V v/(2 mHp^2) A0`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Evidence: project evaluators call `DecayTable::get_gamma_hgaga(2)`; `DecayTable::hgaga` calls `get_coupling_hhh(h,4,4,...)` and multiplies it by `v/(2mHp^2) F_0`; Phase 9 established `get_coupling_hhh=-i C_V`; 2HDMC documents `F_0` with the Djouadi sign.  
Exact aligned Type-I result: `Ahat_gammagamma = -cot(beta) sum_f Nc Qf^2 A_1/2 - [v^2 Z7/(2mHp^2)] A_0`.  
Historical translation: if `g_old=vZ7=C_L`, the scalar term is `-g_old v/(2mHp^2) A_0`.  
Notes: existing widths produced through `DecayTable` are not invalidated by the historical naming ambiguity.

## C13 — Charged-scalar loop convention in `phi -> Z gamma`

Claim: the active 2HDMC object/sign mapping for the charged-Higgs term is fixed.  
Status: **VERIFIED FOR OBJECT/SIGN; EXTERNAL NORMALIZATION CROSS-CHECK OPEN**  
Epistemic class: `DERIVED`, `IMPLEMENTATION-CHECKED`, `OPEN-QUESTION`  
Active code result after factoring the common phase: `Ahat_Hp^(Zgamma)=-K_Z C_V v/(2mHp^2) I_1`, with `K_Z=2cW-1/cW`. At alignment, `C_V=-vZ7`, hence the term is `+K_Z v^2 Z7/(2mHp^2) I_1`.  
Caution: the source code itself states that its chosen normalization reproduces HDECAY but is not consistent with Anatomy II Eqs. 2.23/2.33. The project must not modify this implementation until that normalization difference is independently reproduced.

## Implementation caution carried into the trilinear phase

Phase 9 resolves the implementation layer: `get_qki` uses second components `(-cba,+sba)` for `(h,H)` while `get_coupling_hhh` uses local `Z7=-Lambda7_returned`. The two sign layers combine to reproduce exactly the project Feynman rules `-i v(Z3 sba+Z7 cba)` and `-i v(Z3 cba-Z7 sba)`. The local minus is therefore an implementation convention, not a different Higgs-basis `Z7`.

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
| `PHASE_9_PASS` | **PASS** | charged-Higgs trilinear extracted directly; potential/Lagrangian/Feynman-rule objects and 2HDMC sign layer resolved |
| `PHASE_10_GAMMAGAMMA_PASS` | **PASS** | charged-scalar gamma-gamma loop object, sign and normalization mapped to active 2HDMC |
| `PHASE_10_ZGAMMA_MAPPING_PASS` | **PASS** | active Z-gamma object/sign mapping derived; external normalization discrepancy remains open |
