# Open questions — theory appendix

Only unresolved or source-internal issues are kept here. A project physics claim is removed as a blocker once it has been independently derived and promoted in `VALIDATION_LEDGER.md`.

Resolved project claims so far:
- active 2HDMC generic-potential convention: closed with C1;
- project `h,phi` sign/state map: closed with C2;
- Type-I exact-alignment `kappa_f^phi=-cot(beta)`: closed with C3;
- exact-alignment gauge result `kappa_V^h=1`, `kappa_V^phi=0`: closed with C4;
- tree-level `AVV=0`: closed with C4b.

## Q0.2 — BFLRS11 displayed CP-even signs versus nearby coupling prose

Status: `SOURCE-INTERNAL DISCREPANCY, NON-BLOCKING`  
Epistemic class: `SOURCE`, `DERIVED`  
Observation: an early simple scalar-state display uses global negatives of the DH05 `h,H` fields, while nearby phenomenology prose quotes the standard positive `sin(beta-alpha)` / `cos(beta-alpha)` gauge modifiers.  
What resolves the project side: Phase 4 froze DH field signs and Phase 6 derived gauge couplings directly from the kinetic terms. Project signs no longer depend on reconciling these two source snippets.  
Policy: retain the discrepancy as a source-translation note; never mix a displayed field sign from one convention with a coupling table from another without transforming all odd-field couplings coherently.  
Blocks: nothing in the project derivation.

## Q0.3 — GHOO18 Yukawa table-sign semantics

Status: `SOURCE-TRANSLATION QUESTION, NON-BLOCKING`  
Epistemic class: `SOURCE`, `DERIVED`  
Question: whether every displayed `bar f f H_j : ...` entry in the GHOO18 appendix is best read as a coefficient in `-L_Y`, `L_Y`, or a stripped vertex convention.  
What resolves the project side: Phase 5 derived the Type-I modifiers from the explicit project `-L_Y` and Phase-4 scalar signs, so C3 no longer depends on this table-label semantics.  
Policy: use GHOO18's explicit Lagrangian when a source-to-source sign translation is needed; do not use the compact table to choose a project sign.  
Blocks: nothing in C3.

## Q0.4 — GHOO18 approximate-alignment prose typo

Status: `OPEN, NON-BLOCKING`  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Observation: the approximate-alignment discussion contains a prose statement with `|Z6| >> 1` in a context where the adjacent equations/footnote require small mixing controlled by `Z6`.  
Policy: ignore that sentence as evidence. Phase 7 will derive the Higgs-basis mass matrix directly and decide the alignment criterion algebraically.  
Blocks: nothing if C6 is independently derived.

## Q0.6 — Residual Higgs-basis sign/rephasing

Status: `OPEN BY DESIGN, BLOCKING PHASE 7/9 SIGN CLAIMS`  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Question: how the residual real transformation `H2 -> -H2` is frozen in the project Higgs basis.  
Why it matters: `Y3`, `Z6`, `Z7`, the neutral field in `H2`, and every scalar interaction containing an odd number of `H2` fields flip sign together.  
Constraint inherited from Phase 4: at exact alignment the project state is `phi=-rho_perp`.  
Resolution path: Phase 7 must **define**

`H1=c_beta Phi1+s_beta Phi2`,  
`H2=-s_beta Phi1+c_beta Phi2`

before expanding the potential, and then retain that sign without later rephasing.  
Blocks: absolute signs in C5 and C7–C11.

## Q3.1 — Local quadratic minimum versus global electroweak vacuum

Status: `OPEN BY DESIGN, NON-BLOCKING FOR FORMAL COUPLING DERIVATIONS`  
Epistemic class: `DERIVED`, `OPEN-QUESTION`  
Established: Phase 3 gives a neutral stationary point; Phase 4 gives its physical scalar Hessians. Positive physical squared masses are the corresponding local quadratic conditions.  
Not established: that the point is the global electroweak minimum rather than a local/metastable neutral minimum, or that no deeper competing configuration exists.  
Resolution path: keep global-vacuum/stability validation conceptually separate from the analytic coupling appendix and tie it to validated theory predicates/model-point checks when required.  
Blocks: calling an arbitrary stationary point the global physical vacuum.

## Q4.1 — Exactly degenerate CP-even eigenvalues

Status: `OPEN, NON-BLOCKING FOR CURRENT NON-DEGENERATE BRANCH`  
Epistemic class: `DERIVED`, `OPEN-QUESTION`  
Question: how to phrase alignment/state identity at an exact CP-even degeneracy where diagonalization does not select a unique `alpha`.  
Established away from degeneracy: the DH alpha branch and vacuum-aligned/orthogonal projections fix the state directions up to the already-frozen field signs.  
Resolution path: state the non-degeneracy condition whenever a unique mixing angle is used; treat exact degeneracy separately in the final appendix wording.  
Blocks: only exact-degeneracy statements, including a naive unrestricted `alignment iff Z6=0` wording.

## Q4.2 — Early BFLRS11 restricted mass-term normalization

Status: `SOURCE-INTERNAL DISCREPANCY, NON-BLOCKING`  
Epistemic class: `SOURCE`, `DERIVED`  
Observation: an early pedagogical restricted subsection displays charged/pseudoscalar “mass terms” with factors that do not match the Hessian of the potential printed immediately above.  
Independent evidence: the later general scalar-sector Hessian agrees with our derivation; DH05, GHOO18 rotations, the symbolic audit and active 2HDMC are also consistent.  
Policy: do not use the early restricted formulas as normalization evidence.  
Blocks: nothing downstream.

## Q6.1 — 2HDMC Higgs-basis `Lambda6/Lambda7` versus trilinear-local `Z6/Z7` sign layer

Status: `OPEN, BLOCKING IMPLEMENTATION TRANSLATION FOR C5/C7–C11`  
Epistemic class: `IMPLEMENTATION-CHECKED`, `OPEN-QUESTION`  
Observation: active `THDM::get_param_higgs` computes/returns quantities named `Lambda6,Lambda7`. In `THDM::get_coupling_hhh`, the returned values are stored as local `l6,l7` and then the code explicitly defines `Z6=-l6`, `Z7=-l7` before constructing scalar trilinears.  
Why this matters: a naive statement that “2HDMC Lambda7 is the GHOO/DH Z7” can therefore carry the wrong sign in the scalar-coupling path.  
What is **not** being concluded yet: whether this is merely an internal convention inherited from the cited Feynman-rule convention, a Higgs-basis `H2` sign choice, or another naming translation.  
Resolution path:
1. derive the project Higgs-basis potential from the explicit `H1,H2` rotation;
2. derive project `Z6,Z7` algebraically;
3. evaluate 2HDMC `get_param_higgs` on the same convention;
4. inspect the convention used by `get_coupling_hhh` and map its local signs operator-by-operator;
5. only then use the implementation as evidence for trilinear signs.  
Blocks: implementation-level validation of C5 and C7–C11; does not block the analytic Phase-7 derivation itself.

## Q6.2 — Tree-level gauge decoupling versus loop-induced photonic amplitudes

Status: `RESOLVED CONCEPTUALLY, NON-BLOCKING`  
Epistemic class: `DERIVED`  
Clarification: `kappa_V^phi=0` in exact alignment means no tree-level linear `phi WW` or `phi ZZ` vertex. It does **not** imply `Gamma(phi->gamma gamma)=0` or `Gamma(phi->Z gamma)=0`; those amplitudes can receive fermion and charged-scalar loop contributions.  
Resolution path: preserve this distinction explicitly when Phase 10 constructs loop amplitudes.  
Blocks: nothing now.
