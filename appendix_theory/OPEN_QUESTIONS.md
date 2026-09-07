# Open questions — theory appendix

Only unresolved items are recorded here. A question is removed only when its resolution is entered in `VALIDATION_LEDGER.md` with evidence.

`Q0.1` (active 2HDMC generic-potential convention) was removed after C1 was promoted to `VERIFIED`: the analytic DH05/BFLRS11 potential was frozen first, then the active 2HDMC `set_param_gen` stationarity path was independently checked against the derived equation.

## Q0.2 — BFLRS11 displayed CP-even signs versus gauge-coupling prose

Status: `OPEN, NON-BLOCKING`  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Observation: the displayed simple-chapter fields satisfy `h_B=-h_DH`, `H_B=-H_DH`, while nearby prose quotes the standard positive `sin(beta-alpha)` modifier for the light state.  
Interpretation not yet promoted: this is consistent with scalar-field sign freedom only if every interaction sign is transformed coherently.  
Resolution path: derive gauge couplings from the kinetic terms in the selected DH05 field convention; do not import the BFLRS11 prose sign directly.  
Blocks: nothing downstream if the project uses one explicit convention consistently.

## Q0.3 — Meaning of the GHOO18 Yukawa “coupling” sign

Status: `OPEN`  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Question: In Appendix `Yuk_Type_I`, should each displayed `bar f f H_j : ...` entry be interpreted as a coefficient in `-L_Y`, in `L_Y`, or as a vertex factor stripped of `i`?  
Resolution path: do not use the sign as project evidence. Start from GHOO18's explicit `-L_Y=...` source equation and independently expand/rotate in Phase 5.  
Blocks: `C3 VERIFIED`.

## Q0.4 — GHOO18 source typo in approximate-alignment discussion

Status: `OPEN, NON-BLOCKING`  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Observation: the source text in Appendix `app:al` contains a sentence stating approximate alignment without decoupling with `|Z_6| >> 1`, while the immediately following footnote and surrounding equations require a small mixing controlled by `Z_6` relative to mass splittings.  
Resolution path: rely on the explicit mass matrix and independently derive the alignment condition; do not use this prose sentence as evidence.  
Blocks: must be avoided when validating C6.

## Q0.5 — Project `h,phi` to literature `h,H,H_i` map

Status: `OPEN BY DESIGN`  
Epistemic class: `PROJECT-DEFINITION`, `OPEN-QUESTION`  
Question: Which source mass eigenstate corresponds exactly to project `h` and project `phi` after the project mixing convention is fixed?  
Known project definitions: project `h` is the SM-like CP-even state near 125.13 GeV; project `phi` is the additional CP-even state.  
Resolution path: derive the CP-even matrix and gauge projection, then map by couplings/mixing. Do not map solely by numerical mass ordering.  
Blocks: C2–C4 and downstream trilinears.

## Q0.6 — Residual Higgs-basis rephasing/sign

Status: `OPEN BY DESIGN`  
Epistemic class: `SOURCE`, `OPEN-QUESTION`  
Question: Which residual sign/rephasing of the second Higgs-basis doublet will be frozen for the project (`H_2 -> -H_2` is still allowed in a real CP-conserving Higgs basis)?  
Why it matters: `Y_3`, `Z_6`, `Z_7`, the non-SM CP-even field, and odd-`H_2` trilinears change sign together under this convention change.  
Resolution path: Phase 7 must define the Higgs-basis field rotation explicitly as `H1=Phi1 c_beta+Phi2 s_beta`, `H2=-Phi1 s_beta+Phi2 c_beta`, then retain that sign throughout.  
Blocks: absolute signs in C5, C7–C11 if the field rotation is not frozen.

## Q3.1 — Stationary point versus physical vacuum

Status: `OPEN BY DESIGN`  
Epistemic class: `DERIVED`, `OPEN-QUESTION`  
Question: Do the stationary conditions derived in Phase 3 correspond to the desired local/global electroweak minimum for a given model point?  
What is established: the two neutral tadpoles vanish at the stationary point. This alone does not establish positive physical scalar masses or global-minimum status.  
Resolution path: derive the charged, CP-odd and CP-even Hessians in Phase 4, identify Goldstone zero modes, and keep boundedness/global-vacuum tests conceptually separate.  
Blocks: interpretation of a stationary point as a physical vacuum; does not invalidate the Phase-3 algebra.
