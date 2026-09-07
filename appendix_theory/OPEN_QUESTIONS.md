# Open questions — theory appendix

Only unresolved items are recorded here. A question is removed only when its resolution is entered in `VALIDATION_LEDGER.md` with evidence.

Resolved entries:
- `Q0.1` active 2HDMC generic-potential convention: closed with C1.
- `Q0.5` project `h,phi` state-sign map: closed with C2 after the CP-even Hessian and mixing directions were derived. On the branch connected to `sin(beta-alpha)=1`, the frozen convention is `h=h_DH`, `phi=H_DH`, hence `h=rho_v`, `phi=-rho_perp` at exact alignment. The gauge-coupling interpretation remains a separate C4 question.
- Type-I project fermion modifier: closed with C3 in Phase 5. The project convention gives `kappa_f^phi=-cot(beta)` for all charged fermions; this result no longer depends on interpreting a GHOO table sign.

## Q0.2 — BFLRS11 displayed CP-even signs versus gauge-coupling prose

Status: `OPEN, NON-BLOCKING`  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Observation: the early displayed states satisfy `h_B=-h_DH`, `H_B=-H_DH`, while nearby phenomenology prose quotes the standard positive `sin(beta-alpha)` light-state gauge modifier.  
What Phase 4 added: the project now uses the DH05 field signs explicitly, so no BFLRS11 displayed sign is imported.  
Resolution path: derive gauge couplings from the kinetic terms in Phase 6 and compare only after applying the field-sign translation.  
Blocks: final source-to-source gauge-sign map, not the project derivation.

## Q0.3 — Meaning of the GHOO18 displayed Yukawa sign

Status: `OPEN, NON-BLOCKING`  
Epistemic class: `SOURCE`, `TRANSLATED`, `OPEN-QUESTION`  
Question: for a fully source-to-source sign map, how does each GHOO18 heavy-state label and displayed `bar f f H_j : ...` coefficient map to the project `phi` after all GHOO field signs are frozen?  
What Phase 5 resolved: C3 no longer depends on this question. Starting from the Type-I gauge-invariant Yukawa Lagrangian and the Phase-4 project scalar rotation gives `kappa_f^phi=-cot(beta)` independently. GHOO18 confirms the `Phi2` projection and `cot(beta)` dependence.  
Resolution path: if an exact GHOO coefficient-level comparison is needed in the final appendix, map its neutral rotation matrix `R` and any global heavy-state sign before comparing signed coefficients.  
Blocks: only the optional exact GHOO-to-project signed coefficient map; **does not block C3 or later project derivations**.

## Q0.4 — GHOO18 approximate-alignment prose typo

Status: `OPEN, NON-BLOCKING`  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Observation: Appendix `app:al` contains a prose statement with `|Z6| >> 1` in a context where the adjacent equations/footnote require small mixing controlled by `Z6`.  
Resolution path: ignore that sentence as evidence; derive the Higgs-basis mass matrix and alignment condition directly.  
Blocks: none if C6 is derived independently.

## Q0.6 — Residual Higgs-basis rephasing/sign

Status: `OPEN BY DESIGN`  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Question: which residual sign/rephasing of the second Higgs-basis doublet is frozen for the project?  
Phase-4 constraint: the generic-basis physical convention now fixes `phi=-rho_perp` at exact alignment.  
Resolution path: Phase 7 must define exactly `H1=Phi1 c_beta+Phi2 s_beta`, `H2=-Phi1 s_beta+Phi2 c_beta`; then the relation between `phi` and `sqrt(2) Re H2^0` is fixed rather than chosen later.  
Blocks: absolute signs in C5 and C7–C11.

## Q3.1 — Local quadratic minimum versus global electroweak vacuum

Status: `OPEN BY DESIGN, NON-BLOCKING FOR COUPLING DERIVATIONS`  
Epistemic class: `DERIVED`, `OPEN-QUESTION`  
What Phase 4 resolved: the physical scalar Hessians and their eigenvalues are now explicit. Positive `m_H+^2`, `m_A^2`, and CP-even eigenvalues are the corresponding local quadratic conditions around the neutral stationary point.  
What remains: those conditions do not prove that the point is the global electroweak minimum rather than a metastable/local neutral minimum or that no deeper charge/CP-breaking configuration exists.  
Resolution path: keep global-vacuum/stability checks separate from the appendix coupling derivation and tie them to the validated 2HDMC theory predicates when needed.  
Blocks: calling an arbitrary stationary point the global physical vacuum; does not block formal mass/coupling algebra.

## Q4.1 — Exactly degenerate CP-even eigenvalues

Status: `OPEN, NON-BLOCKING FOR CURRENT NON-DEGENERATE PROJECT BRANCH`  
Epistemic class: `DERIVED`, `OPEN-QUESTION`  
Question: how should the project state map be phrased at an exact CP-even degeneracy where `alpha` is not uniquely fixed by diagonalization?  
What is established: away from degeneracy, the DH05 `alpha` branch and the vacuum-aligned/orthogonal projections uniquely fix the state directions up to the already-frozen signs.  
Resolution path: when writing the final appendix, state the non-degeneracy condition on any equivalence between a unique mixing angle and state identity; do not infer a unique `alpha` at exact degeneracy.  
Blocks: only exact-degeneracy wording.

## Q4.2 — Early BFLRS11 restricted mass-term normalization

Status: `SOURCE-INTERNAL DISCREPANCY, NON-BLOCKING`  
Epistemic class: `SOURCE`, `DERIVED`  
Observation: the early pedagogical restricted subsection around source lines 317–330 displays charged/pseudoscalar “mass terms” with factors that do not match the Hessian of the potential printed immediately above.  
Independent evidence: the later general scalar-sector section defines mass matrices as second derivatives of `V` and agrees with our derivation; DH05, GHOO18 rotations, the CAS audit, and active 2HDMC are also consistent with the derived formulas.  
Resolution policy: do not use the early restricted formulas as normalization evidence; cite the later general section or DH05 for mass formulas.  
Blocks: nothing downstream.