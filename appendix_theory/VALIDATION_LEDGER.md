# Validation ledger — theory appendix

Issue: `#81`  
Ledger rule: `VERIFIED` requires source evidence plus the independent validation required by the mission. Source quotation alone is not sufficient when the claim is designated for derivation.

## C1 — Exact project-compatible scalar potential

Claim: The complete CP-conserving scalar potential used by the project has the Davidson–Haber/Branco operator normalization and the active 2HDMC implementation uses the same convention.  
Status: **PARTIALLY VERIFIED**  
Epistemic class: `SOURCE`, `TRANSLATED`  
Assumptions: CP-conserving real basis for the eventual project form.  
Primary source: DH05 `hbasis.tex`, Eq. `pot`; BFLRS11 `PhysRep_large.tex`, Eq. `2_VH1`.  
Source convention: `+m11^2`, `+m22^2`, `-[m12^2 Phi1†Phi2+h.c.]`; `lambda1,2,5` carry `1/2`.  
Independent derivation: not yet performed.  
Project translation: operator normalization selected; active 2HDMC implementation cross-check deliberately deferred until analytic derivation.  
Dimensional check: quadratic coefficients mass dimension 2; quartics dimension 0.  
Limiting checks: not yet performed.  
Implementation check: not yet performed.  
Numerical check: not yet performed.  
Notes: GHOO18 Eq. `Eq:pot` uses a different quadratic normalization and is not copied symbol-for-symbol.

## C2 — Definition and rotation of physical scalar states

Claim: The project physical states can be obtained from the generic doublets with a fixed `beta` rotation in charged/CP-odd sectors and an `alpha` rotation in the CP-even sector.  
Status: **PARTIALLY VERIFIED**  
Epistemic class: `SOURCE`  
Assumptions: neutral CP-conserving vacuum.  
Primary source: DH05 Eqs. `gpm`, `goldn`, `scalareigenstates`; BFLRS11 simple scalar-sector discussion; GHOO18 Eqs. `Eq:R-def`–`Rmatrix`.  
Source convention: DH05 CP-even signs are selected as the reference convention for later derivation.  
Independent derivation: not yet performed from the quadratic potential.  
Project translation: `h` and `phi` are not yet mapped to source `h/H/H_i`; mapping will be by couplings/mixing, not mass ordering alone.  
Dimensional check: fields dimension 1.  
Limiting checks: not yet performed.  
Implementation check: not yet performed.  
Numerical check: not yet performed.  
Notes: BFLRS11's displayed simple `h,H` fields are both global-sign reversals of DH05's displayed states.

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
Source convention: DH05 also writes the coefficient as `-M12^2`, with `Y3=-M12^2` at `chi=0`; BFLRS11 barred notation therefore has `bar m12^2=+bar lambda6 v^2/2`.  
Independent derivation: required directly from the rotated potential.  
Project translation: pending until Higgs-basis construction.  
Dimensional check: both sides mass dimension 2.  
Limiting checks: `Z6=0 -> Y3=0` is source-consistent but not yet used as proof of alignment.  
Implementation check: pending.  
Numerical check: pending.  
Notes: source agreement is not promoted to VERIFIED before re-derivation.

## C6 — Relation between exact alignment and Z6

Claim: Exact alignment is equivalent to `Z6=0` under the project assumptions.  
Status: **PARTIALLY VERIFIED**  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Assumptions: CP-conserving neutral sector; non-degenerate state identification must be handled carefully; decoupling limit wording must be separated from finite-mass exact alignment.  
Primary source: GHOO18 Appendix `app:al`; DH05 CP-even Higgs-basis mass matrix Eq. `massmhh`.  
Source convention: Higgs-basis `Z6/Lambda6`.  
Independent derivation: pending diagonalization.  
Project translation: pending.  
Dimensional check: `Z6` dimensionless.  
Limiting checks: pending.  
Implementation check: pending.  
Numerical check: pending.  
Notes: mission explicitly requires derivation rather than quotation.

## C7 — Exact phi H+ H- coupling in the Higgs basis

Claim: Exact Higgs-basis expression and exact-alignment reduction of the project `phi H+H-` trilinear.  
Status: **NOT VERIFIED**  
Epistemic class: `OPEN-QUESTION`  
Assumptions: CP conserving; state sign fixed.  
Primary source: GHOO18 cubic-potential section; GHOO18 approximate-alignment table for `q_i`.  
Source convention: `q_i` is explicitly a coefficient in `V`, not directly a Lagrangian coefficient or Feynman rule.  
Independent derivation: required by expanding the Higgs-basis potential.  
Project translation: pending.  
Dimensional check: trilinear coefficient mass dimension 1.  
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
Primary source: none used as a shortcut; exact source rotation exists in DH05.  
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
Dimensional check: expected mass dimension 1 for a trilinear coefficient.  
Limiting checks: pending.  
Implementation check: pending.  
Numerical check: pending.  
Notes: high-risk claim remains unverified.

## C11 — Re-expression using X=lambda6 tan(beta)

Claim: Determine whether the derived trilinear may be written approximately as `-v X cot(beta)`.  
Status: **NOT VERIFIED**  
Epistemic class: `PROJECT-DEFINITION`, `OPEN-QUESTION`  
Assumptions: only meaningful after C10; `X` itself is a project-defined derived coordinate.  
Primary source: none; `X` is not literature notation.  
Independent derivation: blocked on C10.  
Project translation: intentionally deferred.  
Dimensional check: `X` dimensionless.  
Limiting checks: pending.  
Implementation check: pending.  
Numerical check: pending.  
Notes: `X` is intentionally not used elsewhere in Phase 0.

## Phase gates

| Gate | Status | Evidence |
|---|---|---|
| `PHASE_0_PASS` | **PASS** | Source TeX inventory fixes basic field, potential, VEV, mixing, Higgs-basis and trilinear-object conventions sufficiently to start from the doublets. |
| `PHASE_2_PASS` | NOT RUN | Requires exact CP-conserving project potential plus 2HDMC convention cross-check. |
| `PHASE_4_PASS` | NOT RUN | Requires independently derived mass matrices and physical state map. |
| `PHASE_9_PASS` | NOT RUN | Requires exact charged-Higgs trilinear derivation. |
