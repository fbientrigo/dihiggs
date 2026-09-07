# 00 — Source index and Phase-0 evidence inventory

Issue: `fbientrigo/dihiggs#81`  
Mission phase: `PHASE 0 — SOURCE AND CONVENTION INVENTORY`  
Status: `PHASE_0_PASS` for proceeding to field definitions; no downstream coupling claim is verified by this gate.

## Evidence rule

The working priority is

`source TeX > rendered equation in PDF > prose in the same paper > review statement > implementation > general model knowledge`.

No project scan, MadGraph result, or 2HDMC output is used here to establish a theoretical formula.

## Source packs inspected

| Ref. | Primary source | Uploaded archive SHA256 | TeX/source file inspected | Main use in this appendix |
|---|---|---|---|---|
| DH05 | S. Davidson, H. E. Haber, *Basis-independent methods for the two-Higgs-doublet model*, `hep-ph/0504050v5` | `b6f133cc7875c71e3952386f56721f261aa51ccc9cfa1c8ee8aeeb17fbb8db01` | `hbasis.tex` | Generic potential, Higgs-basis rotation, Higgs-basis invariants `Y_i,Z_i`, CP-conserving CP-even convention |
| BFLRS11 | G. C. Branco et al., *Theory and phenomenology of two-Higgs-doublet models*, `arXiv:1106.0034v3` | `6e3693b8fc46a8f26dc546a4389ab3f4f5ead704fbcf0ee5e61bd8524ed19481` | `PhysRep_large.tex` | Standard generic-potential convention, Type-I convention, pedagogical field decomposition, Higgs-basis cross-check |
| GHOO18 | B. Grzadkowski, H. E. Haber, O. M. Ogreid, P. Osland, *Heavy Higgs boson decays in the alignment limit of the 2HDM*, `arXiv:1808.01472v4` | `03c4b3f6a32c30beb4f62645ea941c0735c38fd09b5d99f215e8ce229d3d874b` | `paper_heavyhiggs_jhep_revised3.tex` | Alignment-limit convention, Higgs-basis potential, `Y_3,Z_6,Z_7`, scalar-cubic convention, Type-I specialization |

## Exact source anchors

### DH05 — `hbasis.tex`

- `[SOURCE]` Two identical complex `SU(2)_L` scalar doublets with hypercharge `Y=1`: lines 457–460 and 2187–2191.
- `[SOURCE]` Generic scalar potential: Eq. `\label{pot}`, lines 2193–2207.
- `[SOURCE]` VEV convention with possible phase `xi`: Eq. `\label{potmin}`, lines 2225–2234.
- `[SOURCE]` `s_beta=v_2/v`, `c_beta=v_1/v`, `tan beta=v_2/v_1`: Eq. `\label{tanbdef}`, lines 2252–2259.
- `[SOURCE]` Goldstone directions: Eqs. `\label{gpm}` and `\label{goldn}`, lines 2270–2280.
- `[SOURCE]` Higgs-basis field rotation: Eq. `\label{higgsbasis}`, lines 2313–2316.
- `[SOURCE]` Higgs-basis potential with `M_ij^2, Lambda_i`: Eq. `\label{pothbasis}`, lines 2340–2354.
- `[SOURCE]` Exact generic-to-Higgs-basis quartic relations: Eqs. `\label{Lam1def}`–`\label{Lam7def}`, lines 2366–2410.
- `[SOURCE]` Higgs-basis stationarity in the `M_ij^2,Lambda_i` notation: Eq. `\label{hbasismincond}`, lines 2435–2444.
- `[SOURCE]` Basis-independent definitions of `Y_i,Z_i`: Eqs. `\label{yvv}`–`\label{zvv7}`, lines 2613–2642.
- `[SOURCE]` Their relation to Higgs-basis coefficients: Eq. `\label{hbasisinv}`, lines 2661–2666.
- `[SOURCE]` Invariant stationarity relation `Y_3=-Z_6 v^2/2`: text before Eq. `\label{potmininv2}`, lines 870–880.
- `[SOURCE]` CP-conserving CP-even rotation and Higgs-basis projection: Eqs. `\label{scalareigenstates}`, `\label{hbasis}`, `\label{Hbasis}`, lines 1546–1574.
- `[SOURCE]` Gauge-coupling angle dependence: Eq. `\label{littletable}`, lines 1607–1631.

### BFLRS11 — `PhysRep_large.tex`

- `[SOURCE]` Generic potential in “notation 1”: Eq. `\label{2_VH1}`, lines 6249–6280.
- `[SOURCE]` The authors explicitly state that notation 1 follows Davidson and Haber and warn of sign/factor/conjugation differences across literature: lines 6282–6294.
- `[SOURCE]` CP-conserving pedagogical field decomposition `Phi_a=(phi_a^+,(v_a+rho_a+i eta_a)/sqrt2)^T`: lines 298–306.
- `[SOURCE]` Charged, pseudoscalar and CP-even mass-sector setup in the restricted real `lambda_6=lambda_7=0` presentation: lines 307–379.
- `[SOURCE]` `tan beta=v_2/v_1` and Higgs-basis rotation: lines 371–391.
- `[SOURCE]` Type-I convention: all `u_R,d_R,e_R` couple to `Phi_2`: lines 530–545 and Table `\label{tab:3_models}`.
- `[SOURCE]` CP-even state convention in the phenomenology chapter: lines 566–577.
- `[SOURCE]` Yukawa interaction normalization and Type-I modifiers: Eq. `\label{Eq:Yukawa}` and Table `\label{tab:3_couplings}`, lines 713–765.
- `[SOURCE]` Higgs-basis rotation for a complex neutral vacuum: Eqs. `\label{2_eq:HBT}`–`\label{2_eq:higbas}`, lines 8838–8911.
- `[SOURCE]` Higgs-basis stationarity in barred `m^2,lambda` notation: Eq. `\label{2_stationarity_HB}`, lines 8965–8983.
- `[SOURCE]` Trilinear potential written explicitly as `V_3`: Eqs. `\label{nhydp}` and `\label{mbkpq}`, lines 12446–12485.

### GHOO18 — `paper_heavyhiggs_jhep_revised3.tex`

- `[SOURCE]` General 2HDM67 scalar potential: Eq. `\label{Eq:pot}`, section 2, lines 246–262.
- `[SOURCE]` Generic field/VEV parameterization with phases `xi_j`: Eq. `\label{vevs}`, lines 270–278.
- `[SOURCE]` Goldstone/charged rotation and neutral `R` matrix: section 2, lines 282–333.
- `[SOURCE]` Ordered neutral states `M_1<=M_2<=M_3` and `H_i=R_{ij} eta_j`: Eqs. `\label{Eq:R-def}`–`\label{Rmatrix}`.
- `[SOURCE]` Exact-alignment definition through gauge coupling `e_1=v`, with `alpha_1=beta`, `alpha_2=0`: Eqs. `\label{Eq:e_1}`, `\label{Eq:h1sm-limit}`.
- `[SOURCE]` Higgs-basis field rotation: Eq. in Appendix `\label{HiggsBasis}`, lines 1566–1575.
- `[SOURCE]` Higgs-basis potential with `+ [Y_3 H_1^\dagger H_2+H.c.]`: Eq. `\label{higgspot}`, lines 1580–1590.
- `[SOURCE]` Higgs-basis stationarity `Y_1=-Z_1v^2/2`, `Y_3=-Z_6v^2/2`: Eq. `\label{yz}`, lines 1593–1596.
- `[SOURCE]` Higgs-basis charged mass and neutral mass matrix: Eqs. `\label{mch}`, `\label{mtwo}`, lines 1610–1632.
- `[SOURCE]` Scalar trilinear convention: section `\label{sect:cubic_couplings}`, especially footnote at line 1399 — listed quantities are **coefficients of the potential**; convert to Feynman rules by multiplying by `-i` and the relevant identical-particle combinatorial factor.
- `[SOURCE]` `H_i H^+H^- : q_i` is therefore a potential coefficient: Eq. line 1417.
- `[SOURCE]` Type-I condition `eta_1^{u,d,l}=0` and `tan beta=v_2/v_1`: Eqs. `\label{rhou_I}`–`\label{rhol_I}`, lines 2092–2118.
- `[SOURCE]` In exact alignment with `alpha_3=0`, `H_2` is CP-even and the displayed neutral Type-I coefficient is proportional to `1/tan beta`: lines 2123–2136. The precise project sign is deliberately not promoted in Phase 0.

## Source-to-source compatibility established in Phase 0

1. `[SOURCE][TRANSLATED]` DH05 Eq. `pot` and BFLRS11 Eq. `2_VH1` are the same scalar-potential convention term by term: positive `m_11^2,m_22^2`, negative off-diagonal bilinear `-[m_12^2 Phi_1^dagger Phi_2+h.c.]`, factors `1/2` on `lambda_1,lambda_2,lambda_5`, and unit coefficients on `lambda_3,lambda_4,lambda_6,lambda_7` before adding the Hermitian conjugate.
2. `[SOURCE][TRANSLATED]` GHOO18 uses different symbols/normalization for the quadratic terms. Its Eq. `Eq:pot` has `-1/2` multiplying the entire quadratic bracket. Therefore its `m_ij^2` symbols must not be copied into the Davidson–Haber convention without an explicit translation.
3. `[SOURCE][TRANSLATED]` DH05 Higgs-basis potential uses `-[M_12^2 H_1^dagger H_2+h.c.]` while GHOO18 uses `+[Y_3 H_1^dagger H_2+h.c.]`. DH05 Eq. `hbasisinv` states `Y_3=-M_12^2 e^{-2i chi}`. Thus the differing stationarity signs in the two coefficient notations are a convention translation, not a physical contradiction.
4. `[SOURCE]` All three sources define `tan beta=v_2/v_1` in the relevant real/Type-I basis. DH05 also stresses that `tan beta` is basis-dependent in a fully general 2HDM and only becomes physically meaningful after extra structure such as a Type-I/II Yukawa convention fixes a preferred basis.

## Phase-0 gate

`PHASE_0_PASS = YES` for proceeding to Phase 1.

Reason: the field normalization, generic-potential convention, VEV convention, `beta` definition, principal CP-even rotation conventions, Higgs-basis field rotation, and the distinction between potential coefficients and Feynman rules are identifiable from source TeX.

This gate does **not** verify C1–C11. In particular, it does not yet prove the project sign of the Type-I non-SM Yukawa coupling, the alignment implication for `Z_6`, the `phi H^+H^-` trilinear, any large-`tan beta` limit, or any loop formula.
