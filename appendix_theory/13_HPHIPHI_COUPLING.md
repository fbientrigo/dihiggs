# Phase 13 — Physical `h phi phi` coupling and separation from effective MadGraph prescriptions

## What is established

In exact alignment,

\[
h=\rho_v,\qquad \phi=-\rho_\perp,
\]

so the neutral CP-even Higgs-basis fields may be written

\[
H_1^0=\frac{v+h}{\sqrt2},\qquad
H_2^0=-\frac{\phi}{\sqrt2}
\]

when Goldstones and the CP-odd field are omitted for this cubic extraction.

The physical 2HDM potential contains

\[
\boxed{
V\supset\frac12v(Z_3+Z_4+Z_5)\,h\phi^2
}.
\]

Define the dimension-one physical trilinear by

\[
\boxed{
\mathcal L_{\rm int}\supset-\frac12g_{h\phi\phi}^{\rm phys}\,h\phi^2
}.
\]

Then

\[
\boxed{
g_{h\phi\phi}^{\rm phys}=v(Z_3+Z_4+Z_5)
}.
\]

Because the two `phi` fields are identical, the Feynman rule is

\[
\boxed{h\phi\phi:\quad -i g_{h\phi\phi}^{\rm phys}}.
\]

Using exact alignment and the previously derived mass relations,

\[
\boxed{
g_{h\phi\phi}^{\rm phys}
=\frac{m_h^2+2m_\phi^2-2M^2}{v}
}.
\]

This is the physical 2HDM coupling on the exact-alignment branch. It is not the simplified production prescription used in the stand-in MadGraph model.

## Derivation

Define

\[
A=H_1^\dagger H_1,\qquad
B=H_2^\dagger H_2,\qquad
C=H_1^\dagger H_2.
\]

For the neutral CP-even fields in exact alignment,

\[
A=\frac{(v+h)^2}{2},\qquad
B=\frac{\phi^2}{2},\qquad
C=-\frac{(v+h)\phi}{2}.
\]

Now inspect the Higgs-basis potential operator by operator.

From `Z3 AB`,

\[
Z_3AB\supset\frac12vZ_3h\phi^2.
\]

From `Z4 C C^dagger`, with `C` real in this restricted field slice,

\[
Z_4CC^\dagger\supset\frac12vZ_4h\phi^2.
\]

For real `Z5`,

\[
\frac12Z_5C^2+\mathrm{h.c.}=Z_5C^2
\supset\frac12vZ_5h\phi^2.
\]

`Z6` vanishes in exact alignment. `Z7` produces odd powers of `phi` such as `phi^3`, not `h phi^2`. Quadratic terms cannot generate this cubic monomial.

Therefore

\[
V\supset\frac12vZ_{345}h\phi^2,
\qquad
Z_{345}\equiv Z_3+Z_4+Z_5.
\]

The factor `1/2` is essential: differentiating twice with respect to the identical `phi` fields removes it, giving the vertex `-i v Z345`.

### Reduction to masses and `M^2`

Phase 7 gives, in exact alignment,

\[
Z_1v^2=m_h^2.
\]

From the generic-to-Higgs-basis map and the stationarity conditions one finds

\[
Y_2-
\left(M^2-\frac12Z_1v^2\right)
\propto Z_6.
\]

Hence exact alignment (`Z6=0`) gives

\[
\boxed{Y_2=M^2-\frac12m_h^2}.
\]

The orthogonal CP-even mass is

\[
m_\phi^2=Y_2+\frac12(Z_3+Z_4+Z_5)v^2.
\]

Solving for `Z345`,

\[
Z_3+Z_4+Z_5
=\frac{m_h^2+2m_\phi^2-2M^2}{v^2}.
\]

Thus

\[
\boxed{
g_{h\phi\phi}^{\rm phys}
=\frac{m_h^2+2m_\phi^2-2M^2}{v}}.
\]

This derivation uses only the already-frozen 2HDM potential, exact alignment, and the definition `M^2=m12^2/(s_beta c_beta)`.

## Source cross-check

GHOO18 explicitly states that its cubic quantities are coefficients of the potential. In exact alignment it gives, for the SM-like state `H1` and an orthogonal neutral state `Hj`,

\[
C_V(H_1H_jH_j)=\frac{q_1}{2}+\frac{M_j^2-M_{H^\pm}^2}{v},
\]

where `q1=C_V(H1 H+H-)`. Phase 9 gives `q1=vZ3`, while the Higgs-basis mass relations give

\[
M_j^2-M_{H^\pm}^2=\frac12(Z_4+Z_5)v^2
\]

for the CP-even orthogonal state. Therefore the source expression becomes exactly

\[
\frac12v(Z_3+Z_4+Z_5),
\]

matching the independent extraction above.

## Numerical 2HDMC spot-check

The versioned benchmark `benchmarks/H2scan_mH150_tb300000_production_coupling.json` records a direct call

`THDM::get_coupling_hhh(1,2,2,c)`

at

- `m_h=125.13 GeV`,
- `m_phi=150 GeV`,
- exact alignment,
- `M^2=2.24999999995003345e4 GeV^2`.

The analytic formula predicts

\[
g_{h\phi\phi}^{\rm phys}
=\frac{m_h^2+2m_\phi^2-2M^2}{v}
=63.5914252\ \mathrm{GeV}
\]

using the benchmark SM input for `v`.

2HDMC returns

\[
c=-i\,63.5914252007597\ \mathrm{GeV}.
\]

The agreement is at floating-point precision. This satisfies the issue-81 requirement for a numerical spot-check of the stated coupling mapping against a valid 2HDMC point.

## Separation from effective production prescriptions

The simplified MadGraph/UFO production model is a different theory object. The original PI prescription used in that model is

\[
\boxed{
g_{h\phi\phi}^{\rm PI}=\frac{8m_\phi^2}{v}}
\]

with the UFO encoding the corresponding vertex convention. This is a project-defined effective production prescription, not the physical 2HDM trilinear above.

A second exploratory prescription used in coupling-comparison work is

\[
\boxed{
g_{h\phi\phi}^{\rm EF\_m2}=\frac{8(m_\phi^2-m_2^2)}{v}},
\]

where the selected calculator pathway supplies `m_2^2` but the formula itself is project-defined. It must not be identified with `M^2` or with a 2HDMC definition unless that separate notation problem is resolved.

The appendix should only include an effective prescription if it is actually used in the manuscript. It must be visually separated from the physical 2HDM equation.

## Convention map

- `g_phys` above is defined by `L_int=-(1/2) g_phys h phi^2` and has mass dimension one.
- The physical Feynman rule is `-i g_phys`.
- The PI simplified prescription is an external model input for production and is not obtained by taking a limit of the physical equation.
- Equality between the two numerical couplings would require a particular choice of `M^2`; it is not a model identity.

## What remains uncertain

The physical exact-alignment coupling is closed. The simplified PI prescription is also well-defined as a project prescription. The formal interpretation of the separate `m_2^2` used by `EF_m2` remains outside the physical 2HDM derivation and should not be promoted in the appendix unless the manuscript actually needs it.

## Next smallest validation

When writing the appendix, present `g_phys` first as the physical model result. Put the simplified MadGraph prescription in a clearly labeled paragraph or table row titled “effective production prescription,” not next to the physical equation without qualification.
