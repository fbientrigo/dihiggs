# Canonical evaluator contract v2

Status: authoritative for the maintained 2HDMC evaluation core.

## Scope and provenance

The repository commit recorded in each output row and orchestration manifest is
the authoritative code version for that result. This contract was frozen from
the operational v2 core beginning at commit
`9bc300299208dab1ce0b9f7e7510c3b92b8979f4`; later commits supersede that
baseline only when the output provenance says so.

The evaluators use the repository's patched 2HDMC checkout under `2hdmc/`.
The 2HDMC source provenance is the repository commit of that checkout (the
baseline is `05a6217a949ff947abec2c43cac04f6ba340be8b`); builds link
`2hdmc/lib/lib2HDMC.a`. Output metadata records the evaluator commit and dirty
state. Units are GeV for masses, GeV² for `M²`, `m12_sq`, and widths in GeV,
and millimetres for `ctau_mm`.

Acceptance is layered. Construction is the exact 2HDMC parameter-set return
value. Theory flags report positivity, unitarity, perturbativity, stability,
and the theory-only compatibility flag. `TRIPLE_OK`/`triple_ok_legacy` is
theory-only; `theory_ok_v1` currently equals the three theory predicates. The
M² producer leaves experimental fields unevaluated.

## Canonical lambda1 producer

| Field | Contract |
|---|---|
| Executable | `dihiggs/app/Lambda1EvaluatorV2` |
| Source | `dihiggs/src/Lambda1EvaluatorV2.cpp` |
| Schema | `dihiggs.lambda1.v2` |
| Coordinate | explicit `lambda1_target` |
| API | `THDM::set_param_phys_lam1` |

The input CSV header is exact and fixed:

```text
point_id,mh_gev,mH_gev,mA_gev,mHp_gev,sin_beta_minus_alpha,tan_beta,lambda1_target,lambda6_input,lambda7_input
```

Inputs are persisted as Float64 values with round-trip-safe formatting, while
the raw input lexeme is retained in the output. There is one output row per
attempted input row, including malformed, construction-failing,
theory-rejected, and accepted rows. Rows contain reconstructed quartics and
`m12_sq`, theory flags, partial widths, branching ratios, `width_ok`, and
`ctau_mm`.

## Canonical M² producer

| Field | Contract |
|---|---|
| Executable | `dihiggs/app/DihiggsPointV2Evaluator` |
| Source | `dihiggs/src/DihiggsPointV2Evaluator.cpp` |
| Schema | `dihiggs.point.v2` |
| Coordinates | rectangular `(mH, M²)` |
| API | `THDM::set_param_phys` |

The physical convention is

```text
M² = m12_sq / (sin(beta) * cos(beta))
m12_sq = M² * sin(beta) * cos(beta)
```

`lambda1` is reconstructed output, never a fixed input. The default
`mh = 125.20 GeV` is explicit in the CLI metadata and manifest provenance;
`mHp` and Yukawa type are explicit named inputs. Every attempted grid point
gets a row, including construction failures. Experimental fields remain
unevaluated.

### Canonical `h-H2-H2` observable

`DihiggsPointV2Evaluator` owns the production-coupling observable used by the
LLP downstream stack. For every point that completes 2HDMC construction and
numerical reconstruction, the evaluator calls

```text
THDM::get_coupling_hhh(1, 2, 2, c)
```

on the same `THDM` object used to compute the point's widths and branching
ratios, and exports

```text
g_hH2H2_GeV = abs(Im(c))
```

The frozen convention is

```text
2HDMC: c = -i*g
UFO:   GHphiphi = Im(c) = -g_hH2H2_GeV
```

There is no factorial or symmetry rescaling. The field is a derived observable,
not a scan coordinate, so adding it does not change `point_id`. It remains
`nan` for rows that fail before the coupling can be evaluated. Theory rejection
alone does not mask the observable: a successfully constructed numerical point
still carries its model coupling even if a later theory predicate fails.

### Explicit top-pair width (Gate A, high-mass factory)

Through commit `9f80196`, `width_bb_GeV` and `width_cc_GeV` were exported
explicitly but `H2 -> t tbar` was not; any top-pair contribution above
threshold was silently folded into `width_unaccounted_GeV`. This is corrected:
the evaluator now exports `width_tt_GeV` and `br_tt` from
`DecayTable::get_gamma_huu(2, 3, 3)` on the same `DecayTable` object as the
other partial widths, positioned immediately after `width_cc_GeV` /
`br_cc` in the CSV. This is a breaking column-order change to the
`dihiggs.point.v2` CSV (nine selected widths become ten); consumers must key
by header name, not position.

2HDMC's `get_gamma_huu` for `u1=u2=3` does not impose a hard step at
`M_H2 = 2*m_t^pole` (172.5 GeV pole mass, so 345 GeV); it uses a running-mass
threshold treatment that yields a small nonzero width below the naive
kinematic threshold and a rapid rise through it (empirically: exactly zero at
250 GeV, ~1.5e-7 GeV at 300 GeV, ~9e-5 GeV at 350 GeV, dominant by 400 GeV for
tan_beta=50 Type-I anchors). This is 2HDMC's own physical treatment, not an
artifact of this change; `width_tt_GeV` reports whatever value the shared
`DecayTable` object returns. See `docs/HIGH_MASS_H2_CONTRACT.md` for the full
high-mass audit.

The benchmark `H2scan_mH150_tb300000` freezes the cross-contract anchor

```text
mH = 150 GeV
g_hH2H2_GeV = 63.5914252007596588 GeV
ctau_mm = 4.32622152973311191
br_bb = 0.756737485808578692
```

and the point-v2 tests require the canonical producer to reproduce that anchor
from the same direct `set_param_phys` construction. Downstream repositories
must consume `g_hH2H2_GeV`; they must not rederive the coupling convention.

## Experimental M² tracker

`Phys_M2BandTracker` is a boundary-search helper, not a canonical rectangular
point producer. It is validated only inside bounded pilot domains. Its
intervals and points are non-canonical boundary artifacts and are not final
scientific boundary evidence. It makes no fixed-input-lambda1 claim.

## Legacy compatibility

`PhysScanWithFixings` is the legacy fixed-precision lambda1 path. Its
historical output loses lifetime and branching-ratio information and it is
retained for replay/compatibility only. `PhysScanWithFixings` is not permitted
for new LLP lifetime production. New work must use the v2 producers above.
