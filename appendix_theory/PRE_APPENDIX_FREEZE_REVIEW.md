# Pre-appendix freeze review — Issue #81

Status: **THEORY CORE READY FOR MANUSCRIPT IMPLEMENTATION WITH EXPLICIT CAVEATS**

This file is the final audit boundary before converting `THEORY_APPENDIX_MASTER.md` into polished paper text. It is intentionally stricter than the planned appendix: it records what may be stated as an exact result, what is conditional, and what must not be promoted.

## Executive decision

The first-principles chain is internally consistent through the physical couplings required by issue #81. The generic potential, vacuum conditions, mass/state conventions, Type-I Yukawas, gauge couplings, Higgs-basis map, charged-Higgs trilinear, diphoton charged-scalar loop mapping, and physical `h phi phi` trilinear have all been independently derived before implementation checks.

The main correction produced by the audit concerns the project coordinate

\[
X\equiv\lambda_6\tan\beta.
\]

`X` is useful for organizing scan families but **is not a sufficient one-dimensional description of the charged-Higgs trilinear or photonic transition**. In exact alignment with `lambda7=0`, define

\[
Q\equiv(m_\phi^2-M^2)\tan^2\beta.
\]

Then the exact physical-branch relation is

\[
\boxed{
Z_7=
\left(\frac{X}{2}-\frac{Q}{v^2}\right)\cot\beta
+\left(\frac{X}{2}+\frac{Q}{v^2}\right)\cot^3\beta
}.
\]

Thus the old one-variable statement `Z7 ~ -X cot(beta)` is not a generic fixed-`X` result. It may only be used as an empirical approximation for a restricted family after its numerical error is demonstrated.

## Freeze matrix

| Topic | Freeze status | Manuscript policy |
|---|---|---|
| CP-conserving generic 2HDM potential | **EXACT** | Include explicitly with signs/factors. |
| `v1,v2,v,tan(beta)` and `M^2=m12^2/(s_beta c_beta)` | **EXACT** | Define once; never conflate `M^2`, `m12^2`, `m22^2`. |
| Physical states `h,phi,A,H+/-` | **EXACT** | Keep the frozen DH/project field signs. |
| Exact alignment | **EXACT** | `h=rho_v`, `phi=-rho_perp`; `Z6=0`. |
| Type-I `kappa_f^phi=-cot(beta)` | **EXACT on stated branch** | Include. |
| Tree-level `kappa_V^phi=0` | **EXACT on stated branch** | Include; say tree-level. |
| `Y3=-Z6 v^2/2` | **EXACT in declared Higgs-basis convention** | Include. |
| Exact generic-to-Higgs `Z7` | **EXACT** | Use for physics claims. |
| `Z7 -> -lambda6` at large `tan(beta)` | **CONTROLLED ASYMPTOTIC** | Valid only when `lambda6` is held fixed. |
| `Z7 ~ -X cot(beta)` at fixed `X` | **NOT GENERIC** | Do not state as a model result. |
| Exact `(X,Q)` relation | **EXACT on exact-alignment, lambda7=0 physical branch** | May be used to explain why `X` is incomplete; `Q` need not become a headline variable. |
| `phi H+H-` | **EXACT** | Explicitly state whether coefficient is in `V`, `L_int`, or the Feynman rule. |
| `phi -> gamma gamma` charged-Higgs dependence | **EXACT object/sign mapping** | Include minimally using exact `Z7` or `C_V`. |
| `phi -> Z gamma` | **2HDMC object/sign mapping exact; external normalization open** | Prefer omission of detailed formula unless needed. |
| Physical `h phi phi` | **EXACT** | Include as physical 2HDM production coupling. |
| `g_hphiphi^phys=(mh^2+2mphi^2-2M^2)/v` | **EXACT in exact alignment** | Include with `L_int=-(1/2)g h phi^2`. |
| PI simplified `8mphi^2/v` | **PROJECT EFFECTIVE PRESCRIPTION** | Include only if manuscript discusses simplified production; visibly separate from physical coupling. |
| `EF_m2=8(mphi^2-m_2^2)/v` | **PROJECT ESTIMATE** | Exclude from theory appendix unless manuscript uses it; `m_2^2` formal meaning remains separate. |
| Existing 2HDMC photonic widths | **IMPLEMENTATION-CHECKED** | Not invalidated by old `g=vZ7` naming. |

## Physical `h phi phi` result required by issue #81

Direct Higgs-basis expansion gives

\[
V\supset\frac12v(Z_3+Z_4+Z_5)h\phi^2.
\]

With

\[
\mathcal L_{\rm int}\supset-\frac12g_{h\phi\phi}^{\rm phys}h\phi^2,
\]

one obtains

\[
\boxed{g_{h\phi\phi}^{\rm phys}=v(Z_3+Z_4+Z_5)}.
\]

Exact alignment further gives

\[
\boxed{
g_{h\phi\phi}^{\rm phys}
=\frac{m_h^2+2m_\phi^2-2M^2}{v}
}.
\]

A versioned valid 2HDMC benchmark (`H2scan_mH150_tb300000`) predicts `63.59142520025 GeV`; direct `THDM::get_coupling_hhh(1,2,2)` returns the Feynman rule `-i 63.59142520076 GeV`. The difference is below `1e-9 GeV`.

The simplified PI production prescription is a different object:

\[
\boxed{g_{h\phi\phi}^{\rm PI}=\frac{8m_\phi^2}{v}}.
\]

It must never be presented as the physical 2HDM trilinear or as a limit automatically implied by the physical formula.

## Numerical warning showing why `X` is not enough

For the same versioned valid benchmark,

\[
X=3\times10^{-5},\qquad
Q/v^2\simeq0.742.
\]

The exact two-coordinate expression gives

\[
Z_7\simeq-2.47\times10^{-6},
\]

whereas the historical fixed-`X` approximation gives

\[
-X\cot\beta=-10^{-10}.
\]

The ratio in magnitude is about `2.5e4`. This is a benchmark observation, not a universal factor; its purpose is to show that the missing parameter dependence can dominate in a valid project point.

## Remaining non-blockers

The following do not prevent writing the theory appendix if handled as stated:

- global-vacuum uniqueness is not derived in the appendix; theory-valid benchmark selection remains a separate model-point predicate;
- exact CP-even degeneracy makes the mixing-angle label non-unique, although the Higgs-basis alignment condition remains well defined;
- active 2HDMC `Z gamma` normalization has a documented HDECAY/Anatomy-II convention discrepancy, so no precision standalone `Z gamma` normalization should be asserted here;
- the exploratory `m_2^2` used in `EF_m2` is not identified with project `M^2` and should stay out unless needed by the manuscript.

## Reference completeness

The manuscript reference set should include:

- Davidson & Haber, hep-ph/0504050, for the convention-connected generic/Higgs-basis derivation;
- Branco et al., arXiv:1106.0034, for standard 2HDM/Type-I context;
- Grzadkowski, Haber, Ogreid & Osland, arXiv:1808.01472, for alignment-limit Higgs-basis/cubic conventions;
- Eriksson, Rathsman & Stål, arXiv:0902.0851 / CPC 181 (2010) 189–205, as the canonical 2HDMC calculator reference.

The active vendored 2HDMC source is implementation evidence, not a substitute for the published calculator citation.

## Recommended appendix structure

1. Define the generic potential and vacuum conventions.
2. Define physical states and exact alignment.
3. Give Type-I fermion and gauge modifiers.
4. Introduce the Higgs basis only to the depth required to define `Z6,Z7` and the relevant trilinears.
5. Give `phi H+H-` with one declared interaction convention.
6. Give the minimal charged-scalar `gamma gamma` amplitude dependence on the exact `Z7` and `m_H+`.
7. Give the physical `h phi phi` coupling and then, in a visibly separate paragraph, the effective MadGraph prescription if the production comparison is discussed.
8. Mention `X` only as a project-defined scan coordinate and explicitly state that it is not a sufficient one-dimensional model parameter.

## Final readiness decision

**GO for appendix implementation**, with two hard restrictions:

1. do not use `Z7 ~ -X cot(beta)` as a generic fixed-`X` model equation;
2. do not merge the physical `h phi phi` coupling with the PI simplified MadGraph prescription.

If these restrictions are respected, no unresolved sign or normalization issue blocks the main appendix. The only remaining open normalization issue (`Z gamma`) can be omitted or clearly quarantined without weakening the claims needed by issue #81.
