# Phase 11 — Project coordinate `X=lambda6 tan(beta)`: exact fixed-`X` audit

## What is established

The project-defined coordinate

\[
\boxed{X\equiv\lambda_6\tan\beta}
\]

is useful for organizing scan families, but it is not a standard 2HDM invariant and it is not, by itself, a complete coordinate for the charged-Higgs trilinear or photonic loop coefficient.

In the exact-alignment, CP-conserving project branch with `lambda7=0`, define

\[
\boxed{t\equiv\tan\beta},\qquad
\boxed{Q\equiv(m_\phi^2-M^2)t^2}.
\]

Then the generic quartics reconstructed from the physical branch satisfy exactly

\[
\boxed{\lambda_1=\frac{m_h^2}{v^2}+\frac{Q}{v^2}-\frac32X},
\]

\[
\boxed{\lambda_2=\frac{m_h^2}{v^2}+\frac{Q}{v^2t^4}+\frac{X}{2t^4}}.
\]

Using the exact-alignment expression for `Z7`, one obtains the exact two-coordinate form

\[
\boxed{
Z_7=
\left(\frac{X}{2}-\frac{Q}{v^2}\right)\frac1t
+
\left(\frac{X}{2}+\frac{Q}{v^2}\right)\frac1{t^3}
}.
\]

Therefore, at large `tan(beta)` with `X,Q` finite,

\[
\boxed{
Z_7=
\left(\frac{X}{2}-\frac{Q}{v^2}\right)\cot\beta
+\mathcal O(\cot^3\beta)
}.
\]

This proves that `X` alone does not generically determine the leading charged-Higgs trilinear. A second independent combination, equivalently `Q` or `lambda1-lambda2`, contributes at the same order.

## Derivation

Phase 4/7 give the exact-alignment physical inversion of the generic potential. With

\[
M^2=\frac{m_{12}^2}{s_\beta c_\beta},\qquad
\lambda_6=\frac Xt,\qquad \lambda_7=0,
\]

the diagonal quartics become

\[
\lambda_1
=\frac{m_h^2}{v^2}
+\frac{(m_\phi^2-M^2)t^2}{v^2}
-\frac32X,
\]

\[
\lambda_2
=\frac{m_h^2}{v^2}
+\frac{m_\phi^2-M^2}{v^2t^2}
+\frac{X}{2t^4}.
\]

Defining `Q=(m_phi^2-M^2)t^2` gives the two boxed expressions above.

Phase 10 derived, in exact alignment and `lambda7=0`,

\[
Z_7=-\frac{t(\lambda_1-\lambda_2)+\lambda_6(t^2-1)}{1+t^2}.
\]

Substituting the physical-branch `lambda1-lambda2` and `lambda6=X/t`, and simplifying exactly, gives

\[
Z_7=
\frac{-2Qt^2+2Q+Xv^2t^2+Xv^2}{2t^3v^2},
\]

which is identical to

\[
Z_7=
\left(\frac{X}{2}-\frac{Q}{v^2}\right)t^{-1}
+
\left(\frac{X}{2}+\frac{Q}{v^2}\right)t^{-3}.
\]

No asymptotic expansion was used to obtain this formula.

## Convention map

`X` and `Q` are project-facing coordinates. Neither is a standard basis invariant of a general 2HDM. They are meaningful only after the Type-I/generic basis, exact-alignment branch, and the definition of `M^2` have been frozen.

The earlier fixed-`lambda6` asymptotic statement

\[
Z_7\to-\lambda_6
\]

remains correct when `lambda6` is held fixed. It cannot be converted into a fixed-`X` statement by replacing `lambda6` with `X/t` after taking the limit.

## Numerical check on a versioned valid 2HDMC point

The repository benchmark `H2scan_mH150_tb300000_production_coupling.json` records

- `m_phi=150 GeV`,
- `tan(beta)=300000`,
- `lambda6=1e-10`, hence `X=3e-5`,
- `M2=2.24999999995003345e4 GeV^2`,
- `construction_theory_status=VALIDATED`.

For this point,

\[
Q=(m_\phi^2-M^2)t^2\simeq4.497\times10^4\ \mathrm{GeV}^2,
\]

so

\[
\frac{Q}{v^2}\simeq0.742,
\]

whereas `X=3e-5`. The exact two-coordinate formula therefore gives

\[
Z_7\simeq-2.47\times10^{-6},
\]

while the historical approximation `-X/t` would give `-1e-10`. The discrepancy is about `2.5e4` in magnitude.

This is a numerical observation for this benchmark, not a universal ratio. Its role is to demonstrate that the extra coordinate is phenomenologically relevant in an actual validated project point.

## What was checked against the source

- The exact `Z7(lambda_i,beta)` relation was independently derived in Phase 7 and checked against DH05 and active 2HDMC.
- The exact-alignment `Z6=0` reduction was independently derived in Phase 7.
- The physical-input inversion used above matches active `THDM::set_param_phys`.
- The benchmark values and direct 2HDMC construction status come from the versioned repository artifact named above.
- Symbolic algebra independently reproduces the exact `X,Q` expression.

## What remains uncertain

Nothing about the algebraic conclusion that `X` alone is incomplete. What remains phenomenological is whether a particular restricted benchmark family happens to enforce a numerical relation between `Q` and `X`. That must be demonstrated from the actual points, not assumed.

## Next smallest validation

For any manuscript plot organized by `X`, persist `Q`, exact `Z7`, and the ratio between exact `Z7` and any displayed `X`-only approximation. Do not promote an `X`-only law unless its numerical error is explicitly shown over the plotted family.
