# Open questions — theory appendix

Only unresolved items are recorded here. A question is removed only when its resolution is entered in `VALIDATION_LEDGER.md` with evidence.

Resolved entries:
- `Q0.1` active 2HDMC generic-potential convention: closed with C1.
- `Q0.5` project `h,phi` state-sign map: closed with C2 after the CP-even Hessian and mixing directions were derived. On the branch connected to `sin(beta-alpha)=1`, the frozen convention is `h=h_DH`, `phi=H_DH`, hence `h=rho_v`, `phi=-rho_perp` at exact alignment.
- `C3` project Type-I exact-alignment Yukawa sign: closed in Phase 5 from the explicit `Phi2` Yukawa Lagrangian. The project result is `kappa_f^phi=-cot(beta)` for `u,d,l`.

## Q0.2 — BFLRS11 displayed CP-even signs versus gauge-coupling prose

Status: `OPEN, NON-BLOCKING`  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Observation: the early displayed states satisfy `h_B=-h_DH`, `H_B=-H_DH`, while nearby phenomenology prose quotes the standard positive `sin(beta-alpha)` light-state gauge modifier.  
What Phases 4–5 added: the project now uses the DH05 field signs explicitly, so no BFLRS11 displayed sign is imported into either the state map or Yukawa derivation.  
Resolution path: derive gauge couplings from the kinetic terms in Phase 6 and compare only after applying the field-sign translation.  
Blocks: final source-to-source gauge-sign map, not the project derivation.

## Q0.3 — Semantic status of the GHOO18 displayed Yukawa “couplings”

Status: `OPEN, NON-BLOCKING AFTER PHASE 5`  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Question: in Appendix `Yuk_Type_I`, is each displayed `bar f f H_j : ...` entry a coefficient in `-L_Y`, in `L_Y`, or a vertex factor stripped of `i`?  
What Phase 5 resolved: this semantic ambiguity no longer affects C3 because the project Yukawa modifiers were derived independently from the explicit gauge-invariant Lagrangian and the frozen Phase-4 scalar signs. The matched GHOO18 AL field satisfies `H2_GHOO=-phi_project`, which explains the apparent sign difference.  
Resolution path if a literal GHOO18 vertex citation is needed later: trace its table notation back to the immediately preceding explicit `-L_Y` definition and state the extracted Feynman-rule convention separately.  
Blocks: nothing in the current project derivation.

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

## Q5.1 — Tree-level versus loop-renormalized Yukawa modifiers

Status: `OPEN BY SCOPE, NON-BLOCKING`  
Epistemic class: `DERIVED`, `OPEN-QUESTION`  
What Phase 5 establishes: the exact tree-level Type-I relations `kappa_f^h=c_alpha/s_beta` and `kappa_f^phi=s_alpha/s_beta`, hence `kappa_f^phi=-cot(beta)` in exact alignment.  
What is not claimed: a loop-renormalized effective Yukawa modifier or scheme-dependent higher-order correction.  
Resolution policy: remain at tree level unless a later phenomenology section explicitly requires higher-order Yukawa matching.  
Blocks: nothing in the current theory appendix chain.
