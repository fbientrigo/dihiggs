# Validation ledger — theory appendix

Issue: `#81`  
Ledger rule: `VERIFIED` requires source evidence plus the independent validation required by the mission. Source quotation alone is not sufficient when the claim is designated for derivation.

## F1 — Scalar-doublet and VEV normalization

Claim: `Phi_i=(phi_i^+,(v_i+rho_i+i eta_i)/sqrt(2))^T`, with real neutral VEVs, `v^2=v1^2+v2^2`, `tan beta=v2/v1`.  
Status: **VERIFIED**  
Epistemic class: `SOURCE`, `DERIVED`, `TRANSLATED`, `PROJECT-DEFINITION`  
Assumptions: electromagnetic-preserving CP-conserving neutral vacuum; fixed generic/Type-I basis.  
Primary source: DH05 `potmin`, `tanbdef`; BFLRS11 field expansion; GHOO18 `vevs`.  
Independent derivation: component charges and canonical `1/sqrt(2)` normalization in `02_SCALAR_POTENTIAL.md`.  
Dimensional check: fields and VEVs dimension 1.  
Notes: no physical scalar state was assumed at this stage.

## F2 — Vacuum potential and stationarity

Claim: `V0(v1,v2)` and both neutral tadpole equations in `03_VACUUM_AND_MINIMIZATION.md` follow directly from the frozen potential.  
Status: **VERIFIED**  
Epistemic class: `SOURCE`, `DERIVED`, `IMPLEMENTATION-CHECKED`  
Assumptions: real CP-conserving DH05/BFLRS11 convention.  
Primary source: DH05 `minconditionsa`, `minconditionsb`, specialized only after the independent derivation.  
Independent derivation: every vacuum invariant and both derivatives evaluated term by term.  
Limiting check: `lambda6=lambda7=0` reduces to the softly broken `Z2` form.  
Implementation check: active 2HDMC `set_param_gen` reproduces the derived `m22^2` relation.  
Notes: stationarity alone is not global-vacuum proof.

## F3 — `m22^2`, `m12^2`, and `M^2`

Claim: `M^2=m12^2/(sin beta cos beta)` is a derived coordinate and is generically distinct from `m12^2` and the diagonal coefficient `m22^2`.  
Status: **VERIFIED**  
Epistemic class: `PROJECT-DEFINITION`, `DERIVED`, `TRANSLATED`  
Assumptions: `sin beta cos beta != 0`.  
Independent derivation: substitution into the tadpole equations.  
Dimensional check: all dimension 2; dimensional equality does not imply object identity.  
Project translation: code/data `M2` is this `M^2`, not `m22_2`.

## F4 — Charged, CP-odd and CP-even Hessians

Claim: The three scalar mass matrices in `04_FIELD_ROTATIONS_AND_MASSES.md` are the Hessians of the frozen potential evaluated at the Phase-3 stationary point.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Assumptions: CP-conserving neutral stationary point; real `v1,v2`; nonzero `v1,v2` for the `tan beta` patch.  
Primary source after derivation: DH05 `gpm`, `goldn`, `scalareigenstates`, `chhiggsmass`; BFLRS11 general scalar-sector Hessian; GHOO18 beta rotations.  
Independent derivation: all second derivatives taken directly from `V`; tadpoles inserted only afterward.  
CAS audit: `appendix_theory/checks/phase4_hessian_check.py` reconstructs the potential and verifies both factorized Goldstone-sector matrices and zero modes.  
Implementation check: active 2HDMC reproduces `m_A^2`, `m_H+^2-m_A^2`, and the CP-even matrix in its `m_A^2` representation.  
Dimensional check: every matrix entry has dimension 2.  
Limiting checks: `lambda6=lambda7=0` reduces to the standard softly broken `Z2` matrices.  
Notes: positivity of physical eigenvalues is a local quadratic condition, not global-minimum proof.

## F5 — Goldstone and physical charged/CP-odd rotations

Claim: The vacuum direction is the exact zero eigenvector of both charged and CP-odd Hessians, giving `G+=c_beta phi1+ + s_beta phi2+`, `H+=-s_beta phi1+ + c_beta phi2+`, `G0=c_beta eta1+s_beta eta2`, `A=-s_beta eta1+c_beta eta2`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`  
Independent derivation: direct multiplication of the post-tadpole matrices by `(v1,v2)^T`; orthogonal eigenvector `(-v2,v1)^T`.  
Primary source: DH05 `gpm`, `goldn` and charged-Higgs text; GHOO18 section-2 beta rotation.  
Notes: the beta rotation emerges from the vacuum geometry; it is not postulated from a coupling table.

## F6 — Type-I Yukawa mass generation and neutral CP-even interaction

Claim: In the frozen project Type-I basis only `Phi2` couples to charged fermions, giving `m_f=y_f v2/sqrt(2)` and `L_Y^CP-even=-(m_f/(v s_beta)) rho2 fbar f`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`, `PROJECT-DEFINITION`, `IMPLEMENTATION-CHECKED`  
Assumptions: Type I; project basis `Phi2` Yukawa-active; fermion mass basis with positive real masses; tree level.  
Primary source after derivation: BFLRS11 Type-I assignment and `Eq:Yukawa`; GHOO18 weak-basis Yukawa Lagrangian and `Yuk_Type_I`.  
Independent derivation: expand `Phi2^0=(v2+rho2+i eta2)/sqrt(2)` and `tilde Phi2^0=(v2+rho2-i eta2)/sqrt(2)`, separate mass and interaction terms, then add h.c.  
Project translation: DH05's alternative displayed Type-I special basis is related by `Phi1<->Phi2`; it is not used as the project basis.  
Implementation check: active 2HDMC `set_yukawas_type(1)` scales `rho_D`, `rho_U`, and `rho_L` uniformly by `cot(beta)`.  
Notes: this foundational statement precedes the physical-scalar rotation used in C3.

## C1 — Exact project-compatible scalar potential

Claim: The project uses the complete CP-conserving DH05/BFLRS11 generic-basis potential and active 2HDMC is consistent with it.  
Status: **VERIFIED**  
Epistemic class: `SOURCE`, `DERIVED`, `TRANSLATED`, `IMPLEMENTATION-CHECKED`  
Assumptions: real generic basis; no campaign restriction `lambda7=0` in the definition.  
Primary source: DH05 `pot`; BFLRS11 `2_VH1`.  
Independent derivation: complete gauge-singlet operator basis reconstructed in Phase 2.  
Notes: GHOO18 quadratic symbols require the explicit factor/sign translation already recorded.

## C2 — Definition and rotation of physical scalar states

Claim: Physical charged and CP-odd states arise from the beta rotation, and the CP-even states arise from diagonalizing the derived symmetric `(rho1,rho2)` Hessian with the DH05 alpha convention. On the project branch connected continuously to `sin(beta-alpha)=1`, `h=h_DH` and `phi=H_DH`.  
Status: **VERIFIED**  
Epistemic class: `SOURCE`, `DERIVED`, `TRANSLATED`, `PROJECT-DEFINITION`, `IMPLEMENTATION-CHECKED`  
Assumptions: CP conservation; non-degenerate CP-even eigenvalues for unique mixing direction; project uses the DH05 field signs fixed in Phase 0.  
Primary source: DH05 `scalareigenstates`, `hbasis`, `Hbasis`; GHOO18 beta rotations.  
Independent derivation: charged/odd eigenvectors derived from Hessians; CP-even `alpha` follows by zeroing the off-diagonal element of `R_alpha M_rho^2 R_alpha^T`.  
Project translation: `h=h_DH`, `phi=H_DH`; at exact alignment `h=c_beta rho1+s_beta rho2`, `phi=s_beta rho1-c_beta rho2=-rho_perp`.  
Implementation check: active 2HDMC physical-input branch uses `alpha=beta-asin(sba)` on its non-negative-`cba` branch.  
Notes: state identity is fixed by mixing direction and convention, not mass ordering alone. The SM gauge-coupling interpretation is independently deferred to C4/Phase 6.

## C3 — Type-I exact-alignment result for `kappa_f^phi`

Claim: In the adopted project convention, `kappa_f^phi=-cot(beta)` for up quarks, down quarks and charged leptons.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `TRANSLATED`, `PROJECT-DEFINITION`, `IMPLEMENTATION-CHECKED`  
Assumptions: Type I with `Phi2` Yukawa-active; tree level; Phase-4 scalar signs; exact project alignment `s_(beta-alpha)=+1`.  
Primary source after derivation: BFLRS11 `Eq:Yukawa` and Type-I table; GHOO18 `Yuk_Type_I`.  
Independent derivation: `rho2=c_alpha h+s_alpha phi` and `L_Y=-(m_f/(v s_beta))rho2 fbar f` imply `kappa_f^phi=s_alpha/s_beta=c_(beta-alpha)-s_(beta-alpha)cot(beta)`. Exact alignment gives `-cot(beta)`.  
Cross-source translation: the matched GHOO18 AL field satisfies `H2_GHOO=-phi_project`; its displayed `+cot(beta)` magnitude therefore maps to the project `-cot(beta)` coefficient.  
Implementation check: 2HDMC Type-I `rho_F` magnitude is common `cot(beta)` and uses the DH05 scalar convention; combined with the frozen project scalar sign this reproduces `-cot(beta)`.  
Limiting check: `tan beta -> infinity` suppresses the tree-level non-SM fermion coupling; `kappa_f^h -> 1` exactly in alignment.  
Notes: the conditional width scaling is `Gamma(phi->f fbar) proportional cot^2(beta)` at fixed masses and radiative prescription; no global total-width/lifetime scaling is inferred.

## C4 — Exact-alignment result for `kappa_V^phi`

Claim: `kappa_V^phi=0` in exact alignment.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Assumptions: exact alignment; Phase-4 state identity.  
Primary source: DH05 `littletable`; GHOO18 alignment section.  
Independent derivation: required from `sum_i (D_mu Phi_i)^dagger D^mu Phi_i` in Phase 6.  
Notes: the vacuum-orthogonal mixing geometry anticipates the result but is not being used as a substitute for the requested kinetic-term derivation.

## C5 — Higgs-basis minimization relation for `Y3`

Claim: In the selected Higgs-basis convention, `Y3=-Z6 v^2/2`.  
Status: **PARTIALLY VERIFIED**  
Epistemic class: `SOURCE`  
Assumptions: Higgs-basis potential written with `+[Y3 H1†H2+h.c.]`.  
Primary source: DH05 invariant stationarity; GHOO18 `yz`.  
Independent derivation: still required from the explicit Higgs-basis rotation in Phase 7.  
Notes: generic-basis minimization does not substitute for this convention audit.

## C6 — Exact alignment and `Z6`

Claim: finite-mass exact alignment is controlled by the vanishing Higgs-basis CP-even off-diagonal term and hence by `Z6=0` once the chosen Higgs basis is explicitly derived.  
Status: **PARTIALLY VERIFIED**  
Epistemic class: `SOURCE`, `DERIVED`, `OPEN-QUESTION`  
Assumptions: CP-conserving neutral sector; non-degenerate finite-mass state identification.  
Primary source: DH05 Higgs-basis CP-even matrix; GHOO18 alignment appendix.  
Independent derivation completed so far: generic-basis CP-even mass matrix and alignment directions are derived.  
Remaining validation: derive the Higgs-basis matrix from the explicit field rotation and identify its off-diagonal entry as `Z6 v^2`.  
Notes: decoupling with nonzero `Z6` is an asymptotic route to approximate alignment and must be phrased separately from the finite-mass exact statement.

## C7 — Exact `phi H+ H-` coupling in the Higgs basis

Claim: exact Higgs-basis expression and exact-alignment reduction.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Assumptions: CP conservation; Phase-4 scalar signs retained.  
Primary source: GHOO18 cubic-potential section.  
Independent derivation: required after Higgs-basis construction.  
Notes: potential coefficient, Lagrangian coefficient and Feynman rule remain separate objects.

## C8 — Exact generic-basis expression for `Z7`

Claim: exact analytic map `(lambda_i,beta)->Z7`.  
Status: **PARTIALLY VERIFIED**  
Epistemic class: `SOURCE`  
Primary source: DH05 `Lam7def`, `hbasisinv`.  
Independent derivation: explicit rotation/CAS audit still required.  
Notes: generic `lambda7` is not Higgs-basis `Z7`.

## C9 — Large-`tan beta`, `lambda7=0` limit of `Z7`

Claim: determine whether `Z7 ~ -lambda6` or another relation holds.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Independent derivation: blocked on C8.

## C10 — Large-`tan beta` approximation for `g_(phi H+H-)`

Claim: determine whether `g_(phi H+H-) ~ -v lambda6`.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Assumptions: exact alignment, `lambda7=0`, large `tan beta`, explicit trilinear-object convention.  
Independent derivation: blocked on C7–C9.

## C11 — Re-expression using `X=lambda6 tan(beta)`

Claim: determine whether the derived coupling may be written approximately as `-v X cot(beta)`.  
Status: **NOT VERIFIED**  
Epistemic class: `PROJECT-DEFINITION`, `OPEN-QUESTION`  
Independent derivation: blocked on C10.  
Notes: `X` remains absent from the foundational mass/mixing/Yukawa derivation.

## Phase gates

| Gate | Status | Evidence |
|---|---|---|
| `PHASE_0_PASS` | **PASS** | Source and convention inventory complete enough to begin derivations. |
| `PHASE_1_PASS` | **PASS** | Field quantum numbers, component normalization, VEV convention and beta coordinate derived. |
| `PHASE_2_PASS` | **PASS** | Complete CP-conserving generic potential reconstructed and convention-checked. |
| `PHASE_3_PASS` | **PASS** | Vacuum potential, tadpoles, solved quadratic coefficients and `M^2` translation derived and audited. |
| `PHASE_4_PASS` | **PASS** | Charged, CP-odd and CP-even Hessians derived; Goldstone zero modes proven; beta/alpha rotations obtained; project `h/phi/A/H±` sign convention fixed; source, CAS and 2HDMC checks agree. |
| `PHASE_5_PASS` | **PASS** | Type-I Yukawa Lagrangian expanded in the project `Phi2` basis; masses and neutral CP-even interactions derived; exact `kappa_f^h` and `kappa_f^phi` obtained; C3 closed; source and 2HDMC checks agree. |
| `PHASE_9_PASS` | NOT RUN | Requires exact charged-Higgs trilinear derivation. |
