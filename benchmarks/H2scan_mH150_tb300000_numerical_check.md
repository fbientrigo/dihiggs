# Independent numerical-stability check: `H2scan_mH150_tb300000`

Canonical construction: `THDM::set_param_phys_lam1` (lambda1-input inversion + round-trip check), as run by `dihiggs/app/Lambda1EvaluatorV2` (schema `dihiggs.lambda1.v2`).

Independent construction: `THDM::set_param_phys` called directly (never via `set_param_phys_lam1`, never performing its lambda1 round-trip) with an m12_2 value computed in `benchmarks/check_H2scan_mH150_tb300000.py` from the identical closed-form inversion formula 2HDMC's `set_param_phys_lam1` uses (2hdmc/src/THDM.cpp), fed to the tiny standalone C++ program `benchmarks/check_H2scan_mH150_tb300000.cpp`. No 2HDMC source was modified; the Yukawa convention (`set_yukawas_type(1)`) is unchanged.

## Step 1: reproduction

m12_2 recomputed faithfully in double precision (matching 2HDMC's own `beta_loc=atan(tan_beta); tb=tan(beta_loc)` round-trip, not a mathematically-simplified shortcut): `7.49999999996479594e-02` GeV^2. Constructing directly via `set_param_phys` with this value reproduces the canonical row (max relative difference over total_width_gev/ctau_mm/BRs: **0.000e+00**, i.e. exact to double-precision noise). This confirms the two code paths agree when fed the same m12_2 -- as expected, since `set_param_phys_lam1` itself is exactly this same formula followed by this same `set_param_phys` call.

An earlier version of this check used a mathematically-equivalent but numerically-different shortcut (`tb = tan_beta` directly instead of `tb = tan(atan(tan_beta))`), which differs from 2HDMC's actual computation by a relative ~2.9e-8 in `tb`. Because of the sensitivity below, that tiny difference alone shifted the reconstructed lambda1 from 1.0000009 to 1.956 -- a reminder that at this point, a mathematically-exact reformulation is not numerically equivalent to 2HDMC's own evaluation order, and only an exact reproduction of the formula is a meaningful baseline for the sensitivity probe below.

## Step 2: sensitivity probe (the actual stability check)

The round-trip formula is exactly linear in m12_2 at fixed masses/sba/tan_beta/lambda6/lambda7: d(lambda1_reconstructed)/d(m12_2) = tan_beta / (v^2 cos^2(beta)) = **4.453640e+11** GeV^-2 at this point -- changing m12_2 by ~2e-12 GeV^2 changes the reconstructed lambda1 by 1. A double has ~1.7e-17 absolute resolution near m12_2=0.075 GeV^2, so a floor of order 1e-7 to 1e-6 on `lambda1_abs_residual` is essentially unavoidable in double precision -- this is not a fixable bug in `set_param_phys_lam1`.

A small, bounded probe (perturbing m12_2 only, around the faithfully-reproduced value, for this exact point -- not a scan of nearby mass or tan_beta points) maps the theory-valid window:

| m12_2 offset (GeV^2) | theory_ok | lambda1_reconstructed | total_width_gev | ctau_mm | br_bb | br_gammagamma | br_Zgamma |
|---:|:-:|---:|---:|---:|---:|---:|---:|
| -1e-11 | 0 | 5.45364237810470964e+00 | 4.56084859220861037e-14 | 4.32654091471261903e+00 | 7.56793352287929655e-01 | 2.00296893055586169e-04 | 1.21832726046464183e-05 |
| -1e-12 | 1 | 1.44536831902384866e+00 | 4.56114928624164453e-14 | 4.32625568725017740e+00 | 7.56743460604678009e-01 | 2.58827206243751942e-04 | 1.95640773959976963e-05 |
| -2e-13 | 1 | 1.08908089204948388e+00 | 4.56117805399259659e-14 | 4.32622840117524365e+00 | 7.56738687753830908e-01 | 2.64392296963877709e-04 | 2.03043247334055615e-05 |
| 0.0 | 1 | 1.00000093418204150e+00 | 4.56118529862185007e-14 | 4.32622152973311191e+00 | 7.56737485808578692e-01 | 2.65792941770435500e-04 | 2.04915508605281493e-05 |
| 2e-13 | 1 | 9.10926377063833281e-01 | 4.56119256362645807e-14 | 4.32621463898712566e+00 | 7.56736280486718926e-01 | 2.67197199753718554e-04 | 2.06796245921116123e-05 |
| 1e-12 | 1 | 5.54638950089468286e-01 | 4.56122183089207189e-14 | 4.32618687965472848e+00 | 7.56731424854244050e-01 | 2.72851038474875257e-04 | 2.14404854052144614e-05 |
| 1e-11 | 0 | -3.45363510899139214e+00 | 4.56157403739331286e-14 | 4.32585284777623524e+00 | 7.56672996361025363e-01 | 3.40533440573350223e-04 | 3.09470435606505975e-05 |

Offsets beyond about +/-1e-11 GeV^2 fail `theory_ok` (positivity/unitarity/perturbativity) -- the theory-valid window in m12_2 at this exact (m_H2=150 GeV, tan_beta=3e5) point is only **a few times 1e-12 GeV^2 wide**, over which `lambda1_reconstructed` itself still ranges from about 0.55 to 1.45 (it is not pinned to 1.0 by theory-validity alone).

Within that theory-valid window, the relative spread of each required quantity is:

| Field | Relative spread across theory-valid window |
|---|---:|
| total_width_gev | 0.0016% |
| ctau_mm | 0.0016% |
| br_bb | 0.0016% |
| br_gammagamma | 5.4182% |
| br_Zgamma | 9.5911% |
| br_tautau | 0.0016% |
| br_gg | 0.0016% |

**total_width_gev, ctau_mm, br_bb, br_tautau and br_gg are stable to 0.0016% across the entire theory-valid window** (they are dominated by tree-level Yukawa-driven fermionic widths, which do not depend on lambda1/lambda3). `br_gammagamma` and `br_Zgamma` -- loop-induced, charged-Higgs-mediated, subdominant channels that are not part of the proposed DV+jets (H2->bb) recast -- vary by **5.4182%** and **9.5911%** respectively across the same window, because these loop amplitudes are directly sensitive to lambda3 (which shares lambda1's m12_2-linear dependence and is comparably ill-determined here).

## Operational classification: `NUMERICALLY_UNRESOLVED`

br_Zgamma varies by 9.5911% across the theory-valid m12_2 window (|offset| <= 1e-12 GeV^2), while total_width_gev/ctau_mm/br_bb vary by only 0.0016% over the same window.

This is an operational criterion for the first recast (usability of width/ctau/BR values for an initial detector-sensitivity study), not a claim about final publication precision. **Unstable quantity: `br_gammagamma` and `br_Zgamma` only. `total_width_gev`, `ctau_mm`, and `br_bb` -- the quantities the proposed H2->bb DV+jets recast actually consumes -- are numerically robust at this point**, independent of whether the model is built via `set_param_phys_lam1` or directly via `set_param_phys`.
