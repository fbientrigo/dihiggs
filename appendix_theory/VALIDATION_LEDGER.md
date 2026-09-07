# Validation ledger — theory appendix

Issue: `#81`  
Ledger rule: `VERIFIED` requires source evidence plus the independent validation required by the mission. Source quotation alone is not sufficient when the claim is designated for derivation.

## F1 — Scalar-doublet and VEV normalization

Claim: The project starts from two equal-hypercharge electroweak doublets with `Phi_i=(phi_i^+,(v_i+rho_i+i eta_i)/sqrt(2))^T`, real non-negative VEVs, `v^2=v1^2+v2^2`, and `tan beta=v2/v1`.  
Status: **VERIFIED**  
Epistemic class: `SOURCE`, `DERIVED`, `TRANSLATED`, `PROJECT-DEFINITION`  
Assumptions: electromagnetic-preserving, CP-conserving neutral vacuum; fixed generic/Type-I basis.  
Primary source: DH05 `hbasis.tex`, Eq. `potmin` and `tanbdef`; BFLRS11 explicit field expansion; GHOO18 Eq. `vevs`.  
Source convention: DH05 uses `Y=1` with `Q=T3+Y/2`, equivalent to modern `Y=1/2` with `Q=T3+Y`.  
Independent derivation: component charges and the `1/sqrt(2)` canonical normalization were checked explicitly in `02_SCALAR_POTENTIAL.md`.  
Project translation: neutral real fluctuations are called `rho_i`, CP-odd fluctuations `eta_i`.  
Dimensional check: scalar fields and VEVs have mass dimension 1.  
Limiting checks: not applicable.  
Implementation check: not required for field normalization.  
Numerical check: not required.  
Notes: no physical `h`, `phi`, `A`, `H+` state has been assumed.

## F2 — Vacuum potential and neutral stationarity equations

Claim: Direct substitution of the real neutral VEVs gives the `V0(v1,v2)` and tadpole equations recorded in `03_VACUUM_AND_MINIMIZATION.md`.  
Status: **VERIFIED**  
Epistemic class: `DERIVED`, `SOURCE`, `IMPLEMENTATION-CHECKED`  
Assumptions: selected DH05/BFLRS11 CP-conserving real convention; `v1,v2 != 0` only when solving for `m11^2,m22^2`.  
Primary source: DH05 Eqs. `minconditionsa`, `minconditionsb`, used only after the independent derivation by setting `xi=0` and all coefficients real.  
Source convention: DH05/BFLRS11 generic basis.  
Independent derivation: every invariant was evaluated at the vacuum, Hermitian-conjugate factors were kept explicitly, and `dV0/dv1`, `dV0/dv2` were computed term by term.  
Project translation: `lambda345=lambda3+lambda4+lambda5`.  
Dimensional check: tadpoles have dimension 3; solved quadratic parameters have dimension 2.  
Limiting checks: `lambda6=lambda7=0` reduces to the softly broken `Z2` stationarity form.  
Implementation check: active 2HDMC `THDM::set_param_gen` implements exactly the derived `m22^2` equation, including the `lambda6 c_beta^2 cot(beta)` and `3 lambda7 s_beta c_beta` terms.  
Numerical check: symbolic CAS audit reproduces the hand algebra; no model-point numeric check is needed for this identity.  
Notes: stationarity is not yet a proof of local/global vacuum stability.

## F3 — Distinction among `m22^2`, `m12^2`, and `M^2`

Claim: `M^2 = m12^2/(sin beta cos beta)` is a derived soft-scale coordinate and is generically distinct from both `m12^2` and the diagonal potential coefficient `m22^2`.  
Status: **VERIFIED**  
Epistemic class: `PROJECT-DEFINITION`, `DERIVED`, `TRANSLATED`  
Assumptions: `sin beta cos beta != 0`.  
Primary source: generic-potential definition of `m12^2`; project notation contract for the `M^2` coordinate.  
Source convention: DH05/BFLRS11 sign convention for the off-diagonal bilinear.  
Independent derivation: substitution of `m12^2=M^2 s_beta c_beta` into the independently derived stationarity equations.  
Project translation: `M2` in code/data means this `M^2`; it is not `m22_2`.  
Dimensional check: all three quantities have mass dimension 2, which does not make them equal.  
Limiting checks: definition becomes singular at `v1=0` or `v2=0`, outside the present `tan beta` coordinate patch.  
Implementation check: active canonical project producers already distinguish `M2` from `m12_sq`; no implementation statement is used to derive the definition.  
Numerical check: not required.  
Notes: this resolves the conceptual conflation, not the later physical role of `m22^2` in other literature bases.

## C1 — Exact project-compatible scalar potential

Claim: The complete CP-conserving scalar potential used by the project has the Davidson–Haber/Branco operator normalization and the active 2HDMC generic-basis path is consistent with this convention.  
Status: **VERIFIED**  
Epistemic class: `SOURCE`, `DERIVED`, `TRANSLATED`, `IMPLEMENTATION-CHECKED`  
Assumptions: real CP-conserving generic basis; no campaign restriction `lambda7=0` imposed in the definition.  
Primary source: DH05 `hbasis.tex`, Eq. `pot`; BFLRS11 `PhysRep_large.tex`, Eq. `2_VH1`.  
Source convention: `+m11^2`, `+m22^2`, `-[m12^2 Phi1†Phi2+h.c.]`; `lambda1,2,5` carry `1/2`; `lambda3,4,6,7` have the displayed unit normalization before h.c.  
Independent derivation: the complete dimension-2/dimension-4 gauge-singlet operator basis was enumerated from `Bij=Phi_i†Phi_j`; Hermiticity and CP specialization were applied explicitly.  
Project translation: this convention is frozen as the generic-basis project convention. GHOO18 quadratic symbols are mapped separately and are not copied symbol-for-symbol.  
Dimensional check: quadratic coefficients dimension 2; quartics dimension 0.  
Limiting checks: setting `lambda6=lambda7=0` gives the standard softly broken `Z2` potential.  
Implementation check: active vendored 2HDMC `THDM::set_param_gen` accepts the same generic parameter set and its independently derived `m22_2` stationarity formula exactly matches this potential.  
Numerical check: not required for the operator identity; later model-point checks remain useful for downstream translations.  
Notes: GHOO18 Eq. `Eq:pot` differs only in its quadratic-symbol normalization; the explicit map is in `01_CONVENTION_MAP.md` and `02_SCALAR_POTENTIAL.md`.

## C2 — Definition and rotation of physical scalar states

Claim: The project physical states can be obtained from the generic doublets with a fixed `beta` rotation in charged/CP-odd sectors and an `alpha` rotation in the CP-even sector.  
Status: **PARTIALLY VERIFIED**  
Epistemic class: `SOURCE`  
Assumptions: neutral CP-conserving stationary point.  
Primary source: DH05 Eqs. `gpm`, `goldn`, `scalareigenstates`; BFLRS11 scalar-sector discussion; GHOO18 Eqs. `Eq:R-def`–`Rmatrix`.  
Source convention: DH05 CP-even signs are selected as the reference convention for later derivation.  
Independent derivation: not yet performed from the quadratic potential.  
Project translation: `h` and `phi` are not yet mapped to source `h/H/H_i`; mapping will be by derived mixing/gauge projection, not mass ordering alone.  
Dimensional check: fields dimension 1.  
Limiting checks: pending.  
Implementation check: pending.  
Numerical check: pending.  
Notes: BFLRS11's displayed simple `h,H` states are global-sign reversals of the DH05 displayed states.

## C3 — Type-I exact-alignment result for kappa_f^phi

Claim: In the adopted project convention, `kappa_f^phi=-cot(beta)` for up quarks, down quarks, and charged leptons in exact alignment.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Assumptions: Type I; exact alignment; project state/sign convention.  
Primary source: BFLRS11 Eq. `Eq:Yukawa` and Table `tab:3_couplings`; GHOO18 Appendix `Yuk_Type_I`.  
Source convention: source signs differ with scalar field-sign choices.  
Independent derivation: required from the Type-I Yukawa Lagrangian.  
Project translation: not yet allowed.  
Dimensional check: modifier dimensionless.  
Limiting checks: pending.  
Implementation check: pending.  
Numerical check: pending.  
Notes: high-risk claim remains unverified.

## C4 — Exact-alignment result for kappa_V^phi

Claim: The non-SM project CP-even state has `kappa_V^phi=0` in exact alignment.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Assumptions: exact alignment; state identity fixed by the derived CP-even rotation.  
Primary source: DH05 Eq. `littletable`; GHOO18 alignment section.  
Source convention: source evidence exists, but project state translation has not been derived.  
Independent derivation: required from `sum_i (D_mu Phi_i)^dagger D^mu Phi_i`.  
Project translation: pending.  
Dimensional check: modifier dimensionless.  
Limiting checks: pending.  
Implementation check: pending.  
Numerical check: pending.  
Notes: high-risk claim remains unverified.

## C5 — Higgs-basis minimization relation for Y3

Claim: In the selected Higgs-basis convention, `Y3=-Z6 v^2/2`.  
Status: **PARTIALLY VERIFIED**  
Epistemic class: `SOURCE`  
Assumptions: Higgs-basis potential written with `+[Y3 H1†H2+h.c.]`.  
Primary source: DH05 text preceding Eq. `potmininv2`; GHOO18 Eq. `yz`.  
Source convention: DH05 also writes the coefficient as `-M12^2`, with `Y3=-M12^2` at `chi=0`.  
Independent derivation: still required directly from the explicit Higgs-basis rotation; generic-basis minimization does not substitute for that audit.  
Project translation: pending until Phase 7.  
Dimensional check: both sides dimension 2.  
Limiting checks: `Z6=0 -> Y3=0` is source-consistent but not yet used as proof of alignment.  
Implementation check: pending.  
Numerical check: pending.  
Notes: source agreement remains insufficient for `VERIFIED` under the mission rule.

## C6 — Relation between exact alignment and Z6

Claim: Exact alignment is equivalent to `Z6=0` under the project assumptions.  
Status: **PARTIALLY VERIFIED**  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Assumptions: CP-conserving neutral sector; degeneracies handled explicitly.  
Primary source: GHOO18 Appendix `app:al`; DH05 CP-even Higgs-basis mass matrix Eq. `massmhh`.  
Source convention: Higgs-basis `Z6/Lambda6`.  
Independent derivation: pending diagonalization.  
Project translation: pending.  
Dimensional check: `Z6` dimensionless.  
Limiting checks: pending.  
Implementation check: pending.  
Numerical check: pending.  
Notes: no promotion before mass-matrix derivation.

## C7 — Exact phi H+ H- coupling in the Higgs basis

Claim: Exact Higgs-basis expression and exact-alignment reduction of the project `phi H+H-` trilinear.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Assumptions: CP conserving; state sign fixed.  
Primary source: GHOO18 cubic-potential section.  
Source convention: `q_i` is a coefficient in `V`, not directly a Lagrangian coefficient or Feynman rule.  
Independent derivation: required by expanding the Higgs-basis potential.  
Project translation: pending.  
Dimensional check: trilinear coefficient dimension 1.  
Limiting checks: pending.  
Implementation check: pending.  
Numerical check: pending.  
Notes: no statement `g proportional to Z7 v` is yet promoted.

## C8 — Exact generic-basis expression for Z7

Claim: Exact analytic map `(lambda_i,beta) -> Z7` in the chosen real convention.  
Status: **PARTIALLY VERIFIED**  
Epistemic class: `SOURCE`  
Assumptions: CP-conserving real basis and fixed Higgs-basis sign convention.  
Primary source: DH05 Eq. `Lam7def` together with Eq. `hbasisinv`.  
Source convention: DH05 uses `Lambda7` for the Higgs-basis quartic before identifying it with `Z7` in the real `chi=0` Higgs basis.  
Independent derivation: required by explicit field rotation/CAS-audited expansion.  
Project translation: pending.  
Dimensional check: dimensionless.  
Limiting checks: pending.  
Implementation check: pending.  
Numerical check: pending.  
Notes: generic `lambda7` is not Higgs-basis `Z7`.

## C9 — Large-tan(beta), lambda7=0 limit of Z7

Claim: Determine whether `Z7 ~ -lambda6` or another relation holds.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Assumptions: `lambda7=0`, `tan beta >>1`, project generic-basis convention.  
Primary source: none used as a shortcut.  
Independent derivation: required after C8.  
Project translation: pending.  
Dimensional check: dimensionless.  
Limiting checks: required.  
Implementation check: pending.  
Numerical check: pending.  
Notes: high-risk claim remains unverified.

## C10 — Large-tan(beta) approximation for g_(phi H+ H-)

Claim: Determine whether `g_(phi H+H-) ~ -v lambda6` in the explicitly declared coupling convention.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Assumptions: exact alignment, `lambda7=0`, `tan beta >>1`, state sign and trilinear object fixed.  
Primary source: none sufficient without C7–C9 derivation.  
Independent derivation: required.  
Project translation: pending.  
Dimensional check: expected dimension 1.  
Limiting checks: pending.  
Implementation check: pending.  
Numerical check: pending.  
Notes: high-risk claim remains unverified.

## C11 — Re-expression using X=lambda6 tan(beta)

Claim: Determine whether the derived trilinear may be written approximately as `-v X cot(beta)`.  
Status: **NOT VERIFIED**  
Epistemic class: `PROJECT-DEFINITION`, `OPEN-QUESTION`  
Assumptions: only meaningful after C10.  
Primary source: none; `X` is not literature notation.  
Independent derivation: blocked on C10.  
Project translation: intentionally deferred.  
Dimensional check: `X` dimensionless.  
Limiting checks: pending.  
Implementation check: pending.  
Numerical check: pending.  
Notes: `X` remains absent from the foundational derivation.

## Phase gates

| Gate | Status | Evidence |
|---|---|---|
| `PHASE_0_PASS` | **PASS** | Source TeX inventory fixes the convention inventory. |
| `PHASE_1_PASS` | **PASS** | Field quantum numbers, component normalization, neutral VEV convention and `beta` definition established without physical-state assumptions. |
| `PHASE_2_PASS` | **PASS** | Exact CP-conserving DH05/BFLRS11 potential reconstructed from the operator basis; GHOO18 translation explicit; active 2HDMC generic path independently consistent. |
| `PHASE_3_PASS` | **PASS** | `V0`, both neutral tadpoles, solved `m11^2,m22^2`, and `M^2` translation derived by hand, source-checked, CAS-checked, and implementation-checked for the active `m22_2` path. |
| `PHASE_4_PASS` | NOT RUN | Requires independently derived charged, CP-odd and CP-even mass matrices and physical-state map. |
| `PHASE_9_PASS` | NOT RUN | Requires exact charged-Higgs trilinear derivation. |
