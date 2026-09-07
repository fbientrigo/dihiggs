# Phase 9 — Charged-Higgs trilinear extracted directly from the verified Higgs-basis potential

## What is established

The frozen Higgs basis is
\[
H_1=\begin{pmatrix}G^+\\(v+\rho_v+iG^0)/\sqrt2\end{pmatrix},\qquad
H_2=\begin{pmatrix}H^+\\(\rho_\perp+iA)/\sqrt2\end{pmatrix},
\]
with
\[
\rho_v=s_{\beta-\alpha}h+c_{\beta-\alpha}\phi,\qquad
\rho_\perp=c_{\beta-\alpha}h-s_{\beta-\alpha}\phi.
\]
At exact alignment, `sba=1,cba=0`, so `h=rho_v` and `phi=-rho_perp`.

The only Higgs-basis operators capable of producing a neutral CP-even field times `H+H-` at cubic order are
\[
Z_3(H_1^\dagger H_1)(H_2^\dagger H_2)
\]
and
\[
\{Z_7(H_2^\dagger H_2)H_1^\dagger H_2+\mathrm{h.c.}\}.
\]
Direct expansion gives
\[
\boxed{V\supset v\,(Z_3\rho_v+Z_7\rho_\perp)H^+H^-}.
\]
Therefore, before alignment,
\[
\boxed{C_V^{hH^+H^-}=v(Z_3s_{\beta-\alpha}+Z_7c_{\beta-\alpha})},
\]
\[
\boxed{C_V^{\phi H^+H^-}=v(Z_3c_{\beta-\alpha}-Z_7s_{\beta-\alpha})}.
\]
At exact alignment,
\[
\boxed{C_V^{hH^+H^-}=vZ_3},\qquad
\boxed{C_V^{\phi H^+H^-}=-vZ_7}.
\]

Here `C_V` means the coefficient of the displayed monomial in the scalar potential.

Because `L_int=-V_int`, the literal coefficient multiplying `phi H+H-` in the interaction Lagrangian is
\[
\boxed{C_{\mathcal L}^{\phi H^+H^-}=+vZ_7}
\]
in exact alignment.

If instead one defines a dimension-one coupling by
\[
\boxed{\mathcal L_{\rm int}\supset-g_{\phi H^+H^-}\,\phi H^+H^-},
\]
then
\[
\boxed{g_{\phi H^+H^-}=C_V^{\phi H^+H^-}=-vZ_7}.
\]

Since `phi`, `H+`, and `H-` are distinct external fields, there is no extra identical-particle factorial. The Feynman rule is
\[
\boxed{\phi H^+H^-:\quad -iC_V^{\phi H^+H^-}=+ivZ_7}
\]
in exact alignment.

## Derivation

Define
\[
A_1\equiv H_1^\dagger H_1,\qquad A_2\equiv H_2^\dagger H_2,\qquad C\equiv H_1^\dagger H_2.
\]
Keeping only pieces that can contain one neutral CP-even field and `H+H-`,
\[
A_1=\cdots+v\rho_v+\cdots,
\]
\[
A_2=H^-H^++\cdots,
\]
\[
C=\cdots+\frac v2(\rho_\perp+iA)+\cdots,
\qquad
C^\dagger=\cdots+\frac v2(\rho_\perp-iA)+\cdots.
\]
Thus
\[
Z_3A_1A_2\supset vZ_3\rho_vH^+H^-.
\]
For real `Z7`,
\[
Z_7A_2C+\mathrm{h.c.}=Z_7A_2(C+C^\dagger),
\]
and
\[
C+C^\dagger\supset v\rho_\perp.
\]
Hence
\[
Z_7A_2(C+C^\dagger)\supset vZ_7\rho_\perp H^+H^-.
\]
No other Higgs-basis term contributes at this field order: `Z2` has no VEV insertion; `Z4,Z5,Z6` require Goldstone/other charged fields to generate charged-scalar factors; quadratic terms do not generate a cubic scalar interaction.

Substitute the already-derived physical rotation:
\[
\rho_v=s h+c\phi,\qquad \rho_\perp=c h-s\phi,
\]
where `s=s_(beta-alpha)`, `c=c_(beta-alpha)`. Then
\[
V\supset v[(Z_3s+Z_7c)h+(Z_3c-Z_7s)\phi]H^+H^-.
\]
This proves the general result and its alignment limit without using a coupling table.

## Large-tan(beta) limit and lambda6

Phase 7 established, for `lambda7=0`,
\[
Z_7=-\lambda_6+\frac{\lambda_{345}-\lambda_1}{\tan\beta}+\mathcal O(\tan^{-2}\beta).
\]
Therefore the different trilinear objects behave as
\[
\boxed{C_V^{\phi H^+H^-}= -vZ_7
=+v\lambda_6-\frac{v(\lambda_{345}-\lambda_1)}{\tan\beta}+\cdots},
\]
\[
\boxed{C_{\mathcal L}^{\phi H^+H^-}=+vZ_7
=-v\lambda_6+\frac{v(\lambda_{345}-\lambda_1)}{\tan\beta}+\cdots},
\]
and, under `L_int=-g phi H+H-`,
\[
\boxed{g_{\phi H^+H^-}\simeq +v\lambda_6}.
\]
The Feynman rule tends to
\[
\boxed{-iv\lambda_6}.
\]

With the later project coordinate `X=lambda6 tan(beta)`, so `lambda6=X cot(beta)`,
\[
\boxed{C_V^{\phi H^+H^-}\simeq +vX\cot\beta},
\]
\[
\boxed{C_{\mathcal L}^{\phi H^+H^-}\simeq -vX\cot\beta},
\]
\[
\boxed{g_{\phi H^+H^-}\simeq +vX\cot\beta}\quad\text{if }\mathcal L_{\rm int}=-g\phi H^+H^-,
\]
and
\[
\boxed{\text{Feynman rule}\simeq-i vX\cot\beta}.
\]

## Convention map

GHOO18 explicitly states that its `q_i` is the coefficient of `H_i H+H-` in the potential. Its cubic-coupling appendix further states that potential coefficients become Feynman rules by multiplying by `-i` plus combinatorial factors for identical fields. For `H_i H+H-`, the three fields are distinct, so no extra factorial appears. This agrees with the object separation above.

Active 2HDMC implements
`c=-i v Re(q_{i1} Z3_local + q_{i2} Z7_local)` in `get_coupling_hhh`, after setting `Z7_local=-Lambda7_returned`. Phase 7 proved `Lambda7_returned=Z7_project`. Meanwhile `get_qki` uses second components `(-cba,+sba)` for `(h,H)`. Combining both sign layers gives
\[
-i v(Z_3s+Z_7c)
\]
for `h` and
\[
-i v(Z_3c-Z_7s)
\]
for `H=phi`, exactly the Feynman rules derived from the project potential. Thus the Phase-6 `Z7_local=-l7` mystery is resolved: it compensates the sign convention used in the second `qki` component and is not a different Higgs-basis `Z7`.

## Project-convention conflict discovered

Earlier project material used the shorthand
\[
g_{\phi H^+H^-}=vZ_7\simeq-vX\cot\beta.
\]
With the now-frozen physical state `phi=-rho_perp`, this equality is **not** the coefficient of `phi H+H-` in the potential and is **not** the `g` defined by `L_int=-g phi H+H-`.

It is, however, exactly the literal coefficient in `L_int`:
\[
C_{\mathcal L}^{\phi H^+H^-}=vZ_7.
\]
Therefore the old numerical expression can be retained only if its symbol is explicitly defined as the coefficient multiplying the monomial directly in `L_int`, not as the potential coefficient or as a `-g` convention. The master document must preserve this distinction.

## What was checked against the source

1. Direct symbolic differentiation of the full Higgs-basis potential gives
`d^3 V/(d rho_v dH+ dH-)=v Z3` and
`d^3 V/(d rho_perp dH+ dH-)=v Z7`.
2. Physical-state substitution gives `v(Z3 cba-Z7 sba)` for `phi`.
3. GHOO18 confirms that `q_i` denotes the coefficient in the potential and that the Feynman rule is obtained with `-i` (plus factorials only for identical particles).
4. Active 2HDMC reproduces exactly `-i` times the derived potential coefficient once its `qki` and local `Z7=-l7` conventions are combined.

## What remains uncertain

The analytic trilinear and the 2HDMC Feynman-rule translation are closed. What remains is a **project naming decision**: older plots/notes that call `vZ7` the coupling `g_phiH+H-` must be relabeled or explicitly defined as a Lagrangian monomial coefficient. No physics formula should be changed without tracking which convention those downstream loop-amplitude formulas assumed.

## Next smallest validation

Audit the charged-scalar contribution to `phi -> gamma gamma` and `phi -> Z gamma` starting from a declared interaction-Lagrangian convention. Determine whether the loop formula consumes `C_V`, `C_L`, `g` defined by `L=-g phiH+H-`, or the raw Feynman rule. Only then update any prior `g=vZ7` plotting convention.
