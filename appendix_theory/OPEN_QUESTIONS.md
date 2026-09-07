# Open questions — theory appendix

Only unresolved or source-internal issues are kept here. A project physics claim is removed as a blocker once it has been independently derived and promoted in `VALIDATION_LEDGER.md`.

Resolved project claims through Phase 7:
- active 2HDMC generic-potential convention: C1;
- project `h,phi` sign/state map: C2;
- Type-I exact-alignment `kappa_f^phi=-cot(beta)`: C3;
- exact-alignment gauge result and tree-level `AVV=0`: C4/C4b;
- residual project Higgs-basis sign `H2=-s_beta Phi1+c_beta Phi2`: frozen in Phase 7;
- Higgs-basis stationarity `Y3=-Z6v^2/2`: C5;
- exact alignment iff `Z6=0`: C6;
- exact generic-to-Higgs `Z7` map: C8;
- `lambda7=0`, large-`tan beta`: `Z7 -> -lambda6`: C9.

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

Status: `OPEN, BLOCKING IMPLEMENTATION-LEVEL TRILINEAR SIGN VALIDATION ONLY`  
Epistemic class: `IMPLEMENTATION-CHECKED`, `OPEN-QUESTION`  
Phase-7 resolution of the basis part: `THDM::get_param_higgs` returns quantities named `Lambda6,Lambda7` that match the project/DH `Z6,Z7` formulas **including sign** under `H2=-s Phi1+c Phi2`.  
Remaining observation: `THDM::get_coupling_hhh` receives those values as `l6,l7` and then explicitly defines local `Z6=-l6`, `Z7=-l7`.  
Conclusion allowed now: this extra minus sign is **not** the generic-to-Higgs-basis transformation.  
Still unknown: the exact convention layer represented by those local signs in the routine cited to hep-ph/0602242.  
Resolution path:
1. derive `rho_v H+H-` and `rho_perp H+H-` coefficients directly from the verified Higgs-basis potential;
2. convert to `h,phi`, keeping `phi=-rho_perp` at alignment;
3. distinguish potential coefficient, `L_int` coefficient and Feynman rule;
4. only then map the 2HDMC routine term by term.  
Blocks: implementation-level validation of C7/C10/C11; does not block the analytic trilinear derivation.

## Q6.2 — Tree-level gauge decoupling versus loop-induced photonic amplitudes

Status: `RESOLVED CONCEPTUALLY, NON-BLOCKING`  
Epistemic class: `DERIVED`  
Clarification: `kappa_V^phi=0` in exact alignment removes tree-level linear `phi WW/ZZ`; it does not force `Gamma(phi->gamma gamma)` or `Gamma(phi->Z gamma)` to vanish.
