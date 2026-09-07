# Phase 10 — Audit of the charged-scalar loop in `phi -> gamma gamma` and `phi -> Z gamma`

This phase answers one narrow question left open by Phase 9:

> Which trilinear object is actually consumed by the photonic loop implementation, and what sign enters the interference?

The answer is different from an unqualified statement such as `g_phiH+H-=v Z7`. The active project scans call 2HDMC `DecayTable`, which consumes the **Feynman-rule object** returned by `THDM::get_coupling_hhh`, not a hand-inserted Lagrangian coefficient.

A second issue was discovered while auditing the large-`tan beta` rewrite: the approximation `Z7 ~ -lambda6` is a large-`tan beta` expansion at fixed `lambda6`; it cannot be converted without qualification into `Z7 ~ -X cot(beta)` when `X=lambda6 tan(beta)` is itself held fixed.

---

## What is established

Use the Phase-9 convention

\[
V\supset C_V^S\,S H^+H^-,
\qquad
\mathcal L_{\rm int}\supset-C_V^S\,S H^+H^-.
\]

Thus the Feynman rule is

\[
\boxed{\Gamma_{S H^+H^-}^{\rm FR}=-i C_V^S}.
\]

For the project state `phi`,

\[
C_V^\phi=v(Z_3c_{\beta-\alpha}-Z_7s_{\beta-\alpha}),
\]

and in exact alignment

\[
\boxed{C_V^\phi=-vZ_7}.
\]

The literal coefficient in `L_int` is therefore

\[
\boxed{C_{\mathcal L}^\phi=+vZ_7}.
\]

### Actual project code path

The di-Higgs evaluators/scans obtain photonic widths through

- `DecayTable::get_gamma_hgaga(2)` for `phi -> gamma gamma`;
- `DecayTable::get_gamma_hZga(2)` for `phi -> Z gamma`.

No separate project formula inserts `vZ7` in these calls. Consequently, existing widths produced by this path inherit the internal 2HDMC Feynman-rule convention and are not changed by how we subsequently name `C_V`, `C_L`, or `g` in prose.

---

# Part I — `phi -> gamma gamma`

## Derivation from the active implementation

In `DecayTable::hgaga`, 2HDMC obtains

```cpp
model.get_coupling_hhh(h,4,4,g_hhchc);
...
S_sum = S_sum
      + g_hww*F_1(tau_W)
      + g_hhchc/v*v2/(2.*pow(mHp,2))*F_0(tau_Hp);
```

Since `v2=v^2`, the charged-scalar contribution is

\[
S_{H^\pm}^{\gamma\gamma}
=
\Gamma_{S H^+H^-}^{\rm FR}
\frac{v}{2m_{H^\pm}^2}F_0.
\]

Using

\[
\Gamma_{S H^+H^-}^{\rm FR}=-iC_V^S,
\]

one obtains

\[
\boxed{
S_{H^\pm}^{\gamma\gamma}
=-i\left[
\frac{C_V^S v}{2m_{H^\pm}^2}A_0(\tau_{H^\pm})
\right]
}
\]

because the code states that its `F_0` has the same sign as the Djouadi spin-zero form factor.

The fermion vertex returned by 2HDMC is

\[
-i\frac{m_f}{v}\kappa_f^S,
\]

and its `F_sf` is one half of the conventional spin-1/2 form factor, compensated by an explicit factor 2 in the fermion sum. The W vertex is `+i 2m_W^2 kappa_V/v`; the code defines `F_1` with the opposite sign to Djouadi. Therefore every CP-even term shares one irrelevant global factor `-i` and the reduced amplitude is

\[
\boxed{
\widehat{\mathcal A}_{\gamma\gamma}^S
=
\sum_f N_cQ_f^2\kappa_f^S A_{1/2}(\tau_f)
+\kappa_V^S A_1(\tau_W)
+\frac{C_V^S v}{2m_{H^\pm}^2}A_0(\tau_{H^\pm})
}.
\]

The width is

\[
\boxed{
\Gamma(S\to\gamma\gamma)
=
\frac{\alpha^2m_S^3}{256\pi^3v^2}
\left|\widehat{\mathcal A}_{\gamma\gamma}^S\right|^2
}
\]

for a CP-even state in the present normalization.

This structure agrees with Djouadi, Anatomy II, Eq. (2.23): the charged-scalar term is proportional to the trilinear, `1/m_H+^2`, and the spin-zero form factor. Djouadi defines the scalar trilinear from derivatives of the potential and normalizes the displayed dimensionless `lambda_(H H+H-)` to a dimension-one Feynman-rule unit.

## Exact-alignment project result

Phase 5 and Phase 6 give

\[
\kappa_f^\phi=-\cot\beta,
\qquad
\kappa_V^\phi=0,
\]

while Phase 9 gives

\[
C_V^\phi=-vZ_7.
\]

Hence

\[
\boxed{
\widehat{\mathcal A}_{\gamma\gamma}^\phi
=
-\cot\beta\sum_fN_cQ_f^2A_{1/2}(\tau_f)
-
\frac{v^2Z_7}{2m_{H^\pm}^2}A_0(\tau_{H^\pm})
}.
\]

If one insists on using the historical object

\[
g_{\rm old}\equiv vZ_7=C_{\mathcal L}^\phi,
\]

then the same amplitude is

\[
\boxed{
\widehat{\mathcal A}_{\gamma\gamma}^\phi
=
-\cot\beta\sum_fN_cQ_f^2A_{1/2}
-
\frac{g_{\rm old}v}{2m_{H^\pm}^2}A_0
}.
\]

Therefore **the historical `vZ7` object enters the reduced gamma-gamma amplitude with a minus sign**.

---

# Part II — `phi -> Z gamma`

## What the active code actually computes

For `Z gamma`, the active 2HDMC code uses

```cpp
S_sum = S_sum
      - g_hww*FW(tau_W,lambda_W)
      - (2.*ctw-1./ctw)
        *g_hhchc/v*v2/(2.*pow(mHp,2))
        *FHp(tau_Hp,lambda_Hp);
```

with

\[
F_{H^\pm}=I_1(\tau,\lambda),
\]

and

\[
K_Z\equiv 2c_W-\frac1{c_W}
=\frac{2c_W^2-1}{c_W}.
\]

Substituting the Feynman rule `g_hhchc=-iC_V^S`, the charged-scalar code contribution is

\[
S_{H^\pm}^{Z\gamma}
=+i K_Z\frac{C_V^S v}{2m_{H^\pm}^2}I_1
=-i\left[
-K_Z\frac{C_V^S v}{2m_{H^\pm}^2}I_1
\right].
\]

Thus, after factoring the same global `-i` used for the fermion/W pieces, the active 2HDMC reduced charged-scalar term is

\[
\boxed{
\widehat{\mathcal A}_{H^\pm}^{Z\gamma,\,2HDMC}
=-K_Z\frac{C_V^S v}{2m_{H^\pm}^2}I_1(\tau_{H^\pm},\lambda_{H^\pm})
}.
\]

At exact alignment,

\[
C_V^\phi=-vZ_7,
\]

so

\[
\boxed{
\widehat{\mathcal A}_{H^\pm}^{Z\gamma,\,2HDMC}
=+K_Z\frac{v^2Z_7}{2m_{H^\pm}^2}I_1
}.
\]

In terms of the historical `g_old=vZ7=C_L`,

\[
\boxed{
\widehat{\mathcal A}_{H^\pm}^{Z\gamma,\,2HDMC}
=+K_Z\frac{g_{\rm old}v}{2m_{H^\pm}^2}I_1
}.
\]

Therefore the same historical `vZ7` object enters the active 2HDMC `Z gamma` reduced amplitude with a **plus** sign and the explicit `K_Z` factor.

## Source/implementation discrepancy that must remain explicit

Djouadi Anatomy II Eq. (2.33) contains the charged-Higgs `Z gamma` contribution and Eq. (2.34) defines the reduced `Z H+H-` coupling. However, the active 2HDMC source itself contains this warning:

```cpp
// Charged Higgs contribution above differs with a factor ctw compared to HDECAY.
// The normalisation below gives same result but is not consistent with formulas 2.23 and 2.33 in Anatomy II
```

Therefore two distinct claims must be separated:

1. **VERIFIED for the project implementation:** which trilinear object and sign the active `DecayTable` consumes.
2. **NOT independently closed here:** why the active `Z gamma` normalization chosen to match HDECAY differs from the normalization written in Anatomy II.

The project must not rewrite the active `Z gamma` formula to the Anatomy-II expression merely for aesthetic consistency. A separate HDECAY or diagrammatic cross-check is required before changing that implementation.

---

# Part III — Large-`tan beta` audit: fixed `lambda6` is not fixed `X`

Phase 7 derived, for `lambda7=0`,

\[
Z_7=-\frac{\lambda_6t^4+(\lambda_1-\lambda_{345})t^3-3\lambda_6t^2+(\lambda_{345}-\lambda_2)t}{(1+t^2)^2},
\qquad t\equiv\tan\beta.
\]

### Expansion at fixed `lambda6`

For fixed `lambda6`,

\[
\boxed{
Z_7=-\lambda_6+
\frac{\lambda_{345}-\lambda_1}{t}
+\frac{5\lambda_6}{t^2}
+\mathcal O(t^{-3})
}.
\]

This is the limit in which the statement `Z7 -> -lambda6` is valid.

### Expansion at fixed `X=lambda6 t`

If instead the project coordinate

\[
X\equiv\lambda_6t
\]

is held fixed, then `lambda6=X/t` must be substituted **before** taking the asymptotic limit. One obtains

\[
\boxed{
Z_7=
\frac{\lambda_{345}-\lambda_1-X}{t}
+\frac{5X+2\lambda_1+\lambda_2-3\lambda_{345}}{t^3}
+\mathcal O(t^{-5})
}.
\]

Thus the inference

\[
Z_7\simeq-\lambda_6
\quad\Rightarrow\quad
Z_7\simeq-X\cot\beta
\]

is **not generally valid at fixed `X`**.

### Exact-alignment simplification

For `lambda7=0`, imposing the exact Phase-7 condition `Z6=0` permits elimination of `lambda345`. The exact result simplifies to

\[
\boxed{
Z_7=-\frac{t(\lambda_1-\lambda_2)+\lambda_6(t^2-1)}{1+t^2}
}.
\]

At fixed `X=lambda6 t`,

\[
\boxed{
Z_7=-(\lambda_1-\lambda_2+X)\cot\beta
+\mathcal O(\cot^3\beta)
}.
\]

Therefore the old approximation

\[
Z_7\simeq-X\cot\beta
\]

requires an additional condition such as

\[
|\lambda_1-\lambda_2|\ll|X|
\]

(or a numerically equivalent cancellation). That condition has **not** been established generically.

This downgrades the previous unqualified `C11` statement. `X` still controls one contribution, but `X` alone does not generically determine the leading charged-Higgs trilinear at fixed `X`.

---

# Part IV — Correct leading loop amplitudes at fixed `X`

Assuming exact alignment, `lambda7=0`, finite `lambda1,lambda2`, and fixed `X`, define

\[
\Delta_X\equiv X+\lambda_1-\lambda_2.
\]

Then

\[
C_V^\phi=-vZ_7
= v\,\Delta_X\cot\beta+\mathcal O(\cot^3\beta).
\]

Hence the gamma-gamma reduced amplitude behaves as

\[
\boxed{
\widehat{\mathcal A}_{\gamma\gamma}^\phi
=\cot\beta\left[
-\sum_fN_cQ_f^2A_{1/2}
+\frac{v^2\Delta_X}{2m_{H^\pm}^2}A_0
\right]
+\mathcal O(\cot^3\beta)
}.
\]

For the active 2HDMC `Z gamma` convention,

\[
\boxed{
\widehat{\mathcal A}_{Z\gamma}^{\phi,\,2HDMC}
=\cot\beta\left[
-\mathcal F_f^{Z\gamma}
-K_Z\frac{v^2\Delta_X}{2m_{H^\pm}^2}I_1
\right]
+\mathcal O(\cot^3\beta)
}.
\]

Here `F_f^(Zgamma)` denotes the fermion form-factor sum before multiplying by `kappa_f=-cot(beta)`.

The important robust consequence is that both amplitudes can still factorize as `cot(beta)` at fixed `X` under the stated assumptions, so widths scale approximately as `cot^2(beta)`. What is **not** robust is the claim that the coefficient of that scaling depends only on `X`.

---

## What was checked against the source

1. Active project scans call `DecayTable::get_gamma_hgaga(2)` and `get_gamma_hZga(2)` rather than inserting a hand-defined `g_phiH+H-`.
2. `DecayTable::hgaga` consumes the output of `THDM::get_coupling_hhh` directly and multiplies it by `v/(2m_H+^2) F_0`.
3. `THDM::get_coupling_hhh` was already established in Phase 9 to return the Feynman-rule object `-i C_V`.
4. `DecayTable::F_0` is explicitly documented in code as having the same sign as Djouadi; `F_1` is explicitly documented as having the opposite sign, which produces a common overall phase in the reduced amplitude.
5. Djouadi Anatomy II Eq. (2.23) has the same charged-scalar `gamma gamma` mass suppression and spin-zero form-factor structure.
6. Djouadi Anatomy II Eqs. (2.33)–(2.35) give the charged-scalar `Z gamma` structure and `I_1`; active 2HDMC explicitly records a normalization difference relative to Anatomy II and HDECAY.
7. The fixed-`lambda6` and fixed-`X` asymptotic expansions were checked independently with symbolic algebra.

## What remains uncertain

- The **project/2HDMC object mapping and signs are closed** for both photonic channels.
- Existing 2HDMC-generated photonic widths are not invalidated by the historical `g=vZ7` naming ambiguity.
- The active 2HDMC `Z gamma` normalization versus Anatomy II/HDECAY remains a separate implementation-source audit item. Do not modify the code until it is independently reproduced.
- The approximation `g_old=vZ7 ~ -vX cot(beta)` is not generically established at fixed `X`; it must be checked point-by-point through `lambda1-lambda2` (in exact alignment) or through the exact `Z7` expression.

## Next smallest validation

For the actual benchmark/scanned points, compute and store

\[
R_X\equiv\frac{Z_7}{-X\cot\beta}
\]

and

\[
\Delta_{12}\equiv\lambda_1-\lambda_2.
\]

This immediately tests whether the historical `-X cot(beta)` approximation is numerically valid in the region actually used by the paper, rather than assuming it from the wrong asymptotic limit.
