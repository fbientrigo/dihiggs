# Open questions — theory appendix

Only unresolved or source-internal issues are kept here. A project physics claim is removed as a blocker once it has been independently derived and promoted in `VALIDATION_LEDGER.md`.

Resolved project claims through Phase 10:
- active 2HDMC generic-potential convention: C1;
- project `h,phi` sign/state map: C2;
- Type-I exact-alignment `kappa_f^phi=-cot(beta)`: C3;
- exact-alignment gauge result and tree-level `AVV=0`: C4/C4b;
- residual project Higgs-basis sign `H2=-s_beta Phi1+c_beta Phi2`: frozen in Phase 7;
- Higgs-basis stationarity `Y3=-Z6v^2/2`: C5;
- exact alignment iff `Z6=0`: C6;
- exact generic-to-Higgs `Z7` map: C8;
- `lambda7=0`, large-`tan beta` at fixed `lambda6`: `Z7 -> -lambda6`: C9;
- exact `phi H+H-` potential/Lagrangian/Feynman-rule map: C7;
- active `gamma gamma` charged-scalar loop object/sign mapping: C12;
- active `Z gamma` charged-scalar object/sign mapping: C13 (external normalization comparison remains Q10.1).

## Q0.2 — BFLRS11 displayed CP-even signs versus nearby coupling prose

Status: `SOURCE-INTERNAL DISCREPANCY, NON-BLOCKING`  
Epistemic class: `SOURCE`, `DERIVED`  
Observation: an early simple scalar-state display uses global negatives of the DH05 `h,H` fields, while nearby phenomenology prose quotes the standard positive gauge modifiers.  
Project resolution: Phases 4–7 freeze field signs and derive gauge/Higgs-basis relations independently.  
Policy: retain only as a source-translation warning.

## Q0.3 — GHOO18 Yukawa table-sign semantics

Status: `SOURCE-TRANSLATION QUESTION, NON-BLOCKING`  
Epistemic class: `SOURCE`, `DERIVED`  
Question: whether each compact `bar f f H_j : ...` entry is a coefficient in `-L_Y`, `L_Y`, or a stripped vertex convention.  
Project resolution: Phase 5 derives the Type-I modifiers independently.  
Policy: use GHOO18's explicit Lagrangian, not compact table labels, for sign translation.

## Q0.4 — GHOO18 approximate-alignment prose typo

Status: `SOURCE-INTERNAL TYPO, NON-BLOCKING`  
Epistemic class: `SOURCE`, `DERIVED`  
Observation: a prose statement uses `|Z6| >> 1` in a context whose equations require small mixing.  
Project resolution: Phase 7 derives the CP-even Higgs-basis matrix and proves alignment iff `Z6=0`; the prose sentence is not evidence.

## Q3.1 — Local quadratic minimum versus global electroweak vacuum

Status: `OPEN BY DESIGN, NON-BLOCKING FOR FORMAL COUPLING DERIVATIONS`  
Epistemic class: `DERIVED`, `OPEN-QUESTION`  
Established: neutral stationarity and all physical scalar Hessians.  
Not established: that every stationary point used abstractly is the global electroweak minimum.  
Resolution path: keep global-vacuum/stability validation separate from formal coupling algebra and tie it to validated model-point predicates when needed.

## Q4.1 — Exactly degenerate CP-even eigenvalues

Status: `CLARIFIED, NON-BLOCKING`  
Epistemic class: `DERIVED`  
Phase-7 result: the Higgs-basis matrix has off-diagonal `Z6v^2`; exact degeneracy of a real symmetric 2x2 matrix requires the discriminant to vanish and therefore also `Z6=0`.  
What remains non-unique: the mixing-angle/eigenvector labels inside an exactly degenerate subspace.  
Policy: `alignment iff Z6=0` may be stated geometrically as the VEV direction being a mass eigenvector; avoid claiming a unique `alpha` at exact degeneracy.

## Q4.2 — Early BFLRS11 restricted mass-term normalization

Status: `SOURCE-INTERNAL DISCREPANCY, NON-BLOCKING`  
Epistemic class: `SOURCE`, `DERIVED`  
Observation: an early pedagogical restricted subsection has charged/pseudoscalar factors that do not match the Hessian of its printed potential.  
Independent evidence: the later general Hessian, DH05, symbolic audits and active 2HDMC agree with the project derivation.  
Policy: do not use the early restricted formulas as normalization evidence.

## Q6.1 — 2HDMC trilinear-local `Z6/Z7` sign layer

Status: `RESOLVED`  
Epistemic class: `DERIVED`, `IMPLEMENTATION-CHECKED`  
Resolution: Phase 7 proved `get_param_higgs` returns the project/DH `Z7` sign. Phase 9 then showed that `get_qki` uses second components `(-cba,+sba)` for `(h,H)` while `get_coupling_hhh` defines local `Z7=-l7`. Combining these gives exactly the project Feynman rules `-i v(Z3 sba+Z7 cba)` for `h` and `-i v(Z3 cba-Z7 sba)` for `phi=H`.  
Conclusion: the extra minus is a local implementation convention compensating the `qki` second-component sign; it is not a different Higgs-basis transformation.  
Blocks: none.

## Q9.1 — Historical project shorthand `g_phiH+H-=vZ7`

Status: `PARTLY RESOLVED; NUMERICAL APPROXIMATION REQUIRES REVALIDATION`  
Epistemic class: `PROJECT-DEFINITION`, `DERIVED`, `IMPLEMENTATION-CHECKED`, `OPEN-QUESTION`  
Observation: earlier project material calls `vZ7~-vXcot(beta)` the charged-Higgs trilinear `g`. Phase 9 established that `vZ7` is exactly the literal interaction-Lagrangian coefficient `C_L`, while `C_V=-vZ7`. Phase 10 established that the active photonic widths consume the Feynman rule derived from `C_V`, so existing 2HDMC widths are not invalidated by this naming.  
New correction: at fixed `X=lambda6 tan(beta)`, exact alignment and `lambda7=0`, `Z7=-(X+lambda1-lambda2)cot(beta)+O(cot^3 beta)`. Thus the additional approximation `vZ7~-vXcot(beta)` is not generic.  
Resolution path: for actual paper points compute `R_X=Z7/[-X cot(beta)]` and `Delta12=lambda1-lambda2`; only retain the historical approximation where its numerical accuracy is demonstrated.  
Blocks: semantic migration plus any interpretation that attributes the leading loop coefficient to `X` alone.

## Q10.1 — Active 2HDMC `Z gamma` normalization versus Anatomy II/HDECAY

Status: `OPEN IMPLEMENTATION-SOURCE NORMALIZATION QUESTION, NON-BLOCKING FOR OBJECT/SIGN MAPPING`  
Epistemic class: `IMPLEMENTATION-CHECKED`, `SOURCE`, `OPEN-QUESTION`  
Observation: active `DecayTable::hZga` uses the charged-scalar factor `(2*cW-1/cW)` and explicitly comments that this normalization gives the HDECAY result but is not consistent with Anatomy II Eqs. 2.23/2.33.  
Established: which trilinear object and sign enter the active project calculation is known exactly.  
Not established: an independent derivation of the normalization mismatch and which external convention should be regarded as canonical for a standalone appendix formula.  
Resolution path: reproduce the `Z H+H-` gauge vertex and the scalar loop diagram independently, then compare operator normalization with HDECAY and Djouadi before changing any code or published formula.  
Blocks: only an externally normalized standalone `Z gamma` formula; not existing 2HDMC-generated widths.

## Q6.2 — Tree-level gauge decoupling versus loop-induced photonic amplitudes

Status: `RESOLVED CONCEPTUALLY, NON-BLOCKING`  
Epistemic class: `DERIVED`  
Clarification: `kappa_V^phi=0` in exact alignment removes tree-level linear `phi WW/ZZ`; it does not force `Gamma(phi->gamma gamma)` or `Gamma(phi->Z gamma)` to vanish.
