# First H2 LLP recast candidate — selection report

**Verdict: `NO_VALID_EXISTING_H2_BENCHMARK`.** No existing, post-Yukawa-correction
H2 point in this repository satisfies the LLP-benchmark eligibility gates. See
`benchmarks/FIRST_H2_RECAST_CANDIDATE.json` for the machine-readable record.

- Base commit evaluated: `27817ab156c23546117c93f1584dd4aa766f4850` (`main`)
- Yukawa-order fix commit (the validity boundary): `6bfad7662fd87750d838bf2fe0bd7ac00ee2326a`
  ("fix: initialize canonical evaluator Yukawas after construction", 2026-07-16 22:02:42 -0400)
- Branch: `analysis/select-first-h2-benchmark`, worked in the pre-existing worktree
  `/home/fabi/atlas_dihiggs/_worktrees/select-first-h2-benchmark`

## What was inspected

1. **`docs/pilots/dihiggs_point_v2_v1/`** — the `dihiggs.point.v2` engineering
   pilot: `L01_accepted_anchor.csv`, `L05_theory_rejected_anchor.csv`,
   `L06_llp_anchor.csv`, `ordering_boundary.csv`, `construction_failure.csv`.
2. **`docs/pilots/lambda1_v2_yukawa_fix_v1/`** — the same anchors re-evaluated
   through the independent `dihiggs.lambda1.v2` evaluator
   (`Lambda1EvaluatorV2`), used here only as a cross-check. Widths, BRs, and
   `ctau_mm` agree with (1) to the reported precision for `L01` and `L06`
   (verified below). The `ordering` row is **not** a like-for-like
   reproduction: `dihiggs.point.v2` takes `M2`/`m12_sq` inputs and produces
   two `ordering_boundary` rows, neither `theory_ok`; `dihiggs.lambda1.v2`
   takes `lambda1_target` directly and produces one `ordering` row that *is*
   `theory_ok=1`. This is a real difference in which physical point each
   evaluator's boundary test lands on, not a data-quality problem — but it
   means the two files must not be pooled when counting eligible rows (see
   filtering counts below, which use `dihiggs.point.v2` as the canonical
   basis).
3. **`docs/verification/dihiggs_point_v2_verification_v1.{md,json}`** and
   **`docs/verification/lambda1_v2_yukawa_fix_v1.{md,json}`** — the frozen
   verification reports, including an explicit before/after Yukawa-correction
   table for L01 and L06.
4. **The data lake** at `data_lake_dir` (`project_config.json`,
   `/mnt/c/Users/Asus/cern_db/dihiggs_lake`) — 119 `campaign=*` directories.
5. **`docs/audits/closure_2026-07/REPOSITORY_CLOSURE_INVENTORY_2026-07.{md,json}`**
   and **`/home/fabi/atlas_dihiggs/CLOSURE_HANDOFF_2026-07-18.md`** — confirmed
   this is precisely the deferred item tracked as GitHub issue
   [`fbientrigo/dihiggs#63`](https://github.com/fbientrigo/dihiggs/issues/63),
   "Search for a valid LLP benchmark point to replace L06," opened during the
   July closure and explicitly deferred because a genuine replacement
   "requires a parameter-space search... outside this closure's scope," with
   the directive that "no replacement point be invented during closure."
6. **`hep_cross/contracts/model_point_to_llp_recast_v1.yaml`** — the downstream
   contract an eventual benchmark must satisfy (`sigma_production_fb` among
   other required columns, not present in any pilot row).
7. **`docs/migration/data_contract_v0.1_draft.md`** — the `dihiggs.lambda1.v2`
   executable contract (input header, `Lambda1EvaluatorV2` semantics), used to
   specify the follow-up scan below.

Repository-wide `git log --all --grep=yukawa -i` and `git diff --stat
6bfad766..HEAD -- '*.csv' '*.json'` were used to confirm the Yukawa-fix commit
set and that no CSV/JSON data besides the five pilot anchors and their
verification reports was produced after it.

## Filtering counts

Counting basis: `docs/pilots/dihiggs_point_v2_v1/*.csv` (5 files, 6 rows total
— `ordering_boundary.csv` has 2 rows), the canonical `dihiggs.point.v2`
producer. The `dihiggs.lambda1.v2` cross-check is reported separately above
because its `ordering` row is a different physical point (see note above),
not because its result is discarded.

| Gate | Count |
|---|---:|
| Post-Yukawa-fix H2 candidate rows found (`dihiggs.point.v2`) | 6 (across 5 named anchors) |
| `construction_ok = 1` | 4 (`L01`, `L05`, `L06`, one of two `ordering_boundary` rows) |
| `theory_ok_v1 = 1` after construction | 2 (`L01`, `L06`) |
| Excluded: is `L06` (forbidden by task instructions / issue #63) | 1 |
| **Eligible after exclusions** | **1 (`L01`)** |
| Eligible **and** LLP-scale lifetime (`1 mm <= ctau_mm <= ~10 m`) | **0** |

`L05` fails at the unitarity gate (`theory_ok_v1 = 0`). Both `ordering_boundary`
rows fail theory validity or construction (`mh_gt_mH` construction rejection;
`check_positivity_false`) — neither reaches the eligible set regardless of the
lambda1_v2 cross-check discrepancy noted above. `construction_failure` fails
construction outright.

The entire `data_lake_dir` (119 `campaign=*` directories, most recent write
2026-05-20) predates the Yukawa fix by roughly two months and was excluded
wholesale — every width, branching ratio, and `ctau_mm` value it contains is a
pre-correction (invalid) quantity under this task's rules.

## ctau invariant check

`hep_cross/contracts/model_point_to_llp_recast_v1.yaml` declares
`ctau_mm == hbar_c / total_width_GeV` (`hbar_c = 1.973269804e-13 GeV·mm`) as a
required invariant. `benchmarks/verify_pilot_ctau_invariant.py` checks this
against every `width_ok=1` row in both pilot files:

```
$ python3 benchmarks/verify_pilot_ctau_invariant.py
docs/pilots/dihiggs_point_v2_v1/L01_accepted_anchor.csv: 1 width_ok=1 row(s) checked
docs/pilots/dihiggs_point_v2_v1/L05_theory_rejected_anchor.csv: 1 width_ok=1 row(s) checked
docs/pilots/dihiggs_point_v2_v1/L06_llp_anchor.csv: 1 width_ok=1 row(s) checked
docs/pilots/dihiggs_point_v2_v1/ordering_boundary.csv: 1 width_ok=1 row(s) checked
docs/pilots/dihiggs_point_v2_v1/construction_failure.csv: 0 width_ok=1 row(s) checked
docs/pilots/lambda1_v2_yukawa_fix_v1/lambda1_v2_pilot.csv: 4 width_ok=1 row(s) checked
PASS: ctau_mm == hbar_c_GeV_mm / total_width for all width_ok=1 pilot rows
```

All 8 checked rows satisfy the invariant to floating-point precision (rel_err
`0.0`), so "reproducible `ctau`" (one of the minimum validity requirements) is
confirmed for every row that reports a width, including `L01`.

## Why no point survives

`L01_accepted_anchor` is the only row that reaches the final eligibility gate
(construction- and theory-valid, not `L06`, not a boundary/failure case), so
it is the only one that needs a full criterion-by-criterion record:

| Minimum validity requirement | `L01` | Note |
|---|---|---|
| Generated after Yukawa-order correction | ✓ | `producer_commit = 6bfad766...`, listed in `docs/verification/dihiggs_point_v2_verification_v1.md`'s before/after table |
| Successful model construction | ✓ | `construction_ok = 1` |
| Theory-valid per canonical gates | ✓ | `theory_ok_v1 = 1` (positivity, unitarity, perturbativity, stability all `1`) |
| Unambiguous H2 identity and mass | ✓ | `mH_input_GeV = 130.0` (2HDMC scalar 2, per `docs/migration/data_contract_v0.1_draft.md`) |
| Finite positive total width | ✓ | `total_width_GeV = 5.441e-06` |
| Reproducible `ctau` | ✓ | verified above, rel_err `0.0` |
| BRs from the same corrected point | ✓ | same row, same `producer_commit` |
| Production params / source row recoverable | ✓ | full row in `docs/pilots/dihiggs_point_v2_v1/L01_accepted_anchor.csv` |
| Producer commit / provenance recorded | ✓ | `producer_commit`, `producer_dirty=no`, `evaluator_api` all present |
| **LLP-relevant lifetime** | **✗** | `ctau_mm = 3.63e-08` — decays effectively promptly, ~7 orders of magnitude below any ATLAS DV scale |
| **No known invalid/provisional status** | **✗** | the pilot set is scoped "engineering validation only; no campaign-level scientific conclusion" (`docs/verification/dihiggs_point_v2_verification_v1.md`); it was never intended as a physics-conclusion-bearing scan point |

`L01` independently fails two of the criteria — a short-lifetime physics
failure and a provisional-status/provenance-scope failure — either of which
alone would disqualify it.

- **`L06_llp_anchor`** (`m_H2 = 200 GeV`, Type I, alignment `sin(β-α)=1`,
  `tan β = 10⁴`) is theory-valid and has `br_bb = 0.675`, but is the point this
  task explicitly forbids from being reused, and on its own corrected numbers
  it no longer qualifies anyway: `ctau_mm` collapsed from `1.79 × 10⁴ mm`
  (pre-fix) to `3.40 × 10⁻³ mm` (post-fix) — three orders of magnitude below
  even the innermost ATLAS displaced-vertex scale.
- **`L01_accepted_anchor`** (`m_H2 = 130 GeV`, near-alignment `sin(β-α)=0.999`,
  `tan β = 50`) is theory-valid and construction-valid, with a reasonable
  `br_bb = 0.31`, but its couplings are essentially unsuppressed at this
  `tan β`, giving `ctau_mm = 3.63 × 10⁻⁸ mm` — a promptly-decaying particle,
  not an LLP candidate by any margin.
- **`L05_theory_rejected_anchor`** fails the unitarity gate outright
  (`theory_ok_v1 = 0`); no physics interpretation is available.
- **`ordering_boundary`** and **`construction_failure`** are deliberately
  constructed edge cases for exercising rejection code paths (degenerate
  `mH == mh`, negative reconstructed mass²), not physical scan points; the
  verification docs describe the whole pilot set as "engineering validation
  only; no campaign-level scientific conclusion."

No point in the corrected data simultaneously satisfies theory validity,
finite positive width, and an LLP-relevant lifetime. This matches the prior
finding recorded in issue #63 and in `paper/notas/CLAIM_LEDGER.md` (claim
`C010`, `SUPERSEDED`).

## What the recast agent needs next

A small, decision-oriented follow-up scan (not run by this task) is required
before a first H2 recast point exists. See the specification below.

## Follow-up scan specification

**Goal:** find the first `dihiggs.lambda1.v2`-valid H2 point with
LLP-relevant lifetime, reusing L06's suppression strategy (Type I, exact
alignment, near-zero `λ6`/`λ7`, large `tan β`) but pushed far enough to survive
the corrected Yukawa installation order.

- **Evaluator / contract:** `dihiggs/app/Lambda1EvaluatorV2`, input header
  exactly `point_id,mh_gev,mH_gev,mA_gev,mHp_gev,sin_beta_minus_alpha,tan_beta,lambda1_target,lambda6_input,lambda7_input`
  (`docs/migration/data_contract_v0.1_draft.md`).
- **Fixed:** `mh_gev = 125.13`, `sin_beta_minus_alpha = 1.0` (exact alignment),
  `lambda1_target = 1.0`, `lambda6_input = 1e-10`, `lambda7_input = 0.0`
  (Type I fermiophobic-heavy-Higgs configuration, as in `L06`).
- **Why `tan_beta` and how far to push it:** at exact alignment in Type I,
  `H2 -> WW/ZZ` vanishes and every open fermionic channel scales as
  `1/tan^2(beta)`, so `Gamma_total ~ 1/tan^2(beta)` in this configuration.
  Anchoring on `L06`'s measured `Gamma = 5.808e-11 GeV` at `tan_beta = 1e4`:
  `ctau_mm = 1 mm` extrapolates to `tan_beta ~ 1.72e5`, and `ctau_mm = 10 mm`
  to `tan_beta ~ 5.43e5` (`Gamma_target = hbar_c / ctau_target`,
  `tan_beta_new = tan_beta_L06 * sqrt(Gamma_L06 / Gamma_target)`). This is a
  falsifiable prediction the scan should confirm or refute, not an assumption
  to build the benchmark on.
- **Swept, bounded:**
  - `tan_beta` in `{1e4, 3e4, 1e5, 3e5, 1e6}` (5 values; brackets the
    `~1.7e5`-`5.4e5` window predicted above on both sides).
  - `mH_gev` (= `m_H2`) in `{150, 200, 250}` GeV (3 values; stays inside a
    mass range plausibly reachable by the validated DV+jets recast).
  - `mA_gev = mHp_gev = mH_gev + 300` GeV for each row (keeps the same
    moderate custodial splitting as `L06`, avoids introducing a second free
    dimension).
- **Maximum points:** 15 (5 × 3 grid), run as a single small CSV through the
  existing evaluator — no campaign harness, no adaptive search.
- **Required output columns:** the full `dihiggs.lambda1.v2` row (all
  `construction_ok`, `theory_ok`/`triple_ok`, `width_ok`, `total_width_gev`,
  `ctau_mm`, all nine `br_*` columns, `lambda1_abs_residual`,
  `lambda1_residual_warning`, and full input/reconstructed-parameter
  provenance) — i.e., run it unmodified through the existing evaluator, do not
  hand-pick a subset of fields.
- **Early-stop condition:** stop at the first row with `theory_ok = 1`,
  `width_ok = 1`, `1 mm <= ctau_mm <= 1e4 mm`, `br_bb > 0.05`, **and**
  `lambda1_residual_warning = 0`. The residual-warning requirement matters
  here: `L06` already carries `lambda1_residual_warning = 1`
  (`lambda1_abs_residual = 2.06e-09`) at `tan_beta = 1e4`, and pushing to
  `1e5`-`1e6` drives `cos(beta)` toward `~1e-5`-`1e-6`, which will plausibly
  worsen the `m12_sq`/`M2` round-trip further. A row that satisfies the
  lifetime/BR window but still carries the residual warning should be flagged
  for numerical-stability review, not silently accepted as the benchmark. If
  no row in this 15-point grid qualifies cleanly, report
  `NO_VALID_EXISTING_H2_BENCHMARK` again with the extended grid results and
  escalate to PI for a wider search (this task's mandate stops here).
- **Not included in this step:** `sigma_production_fb` (production cross
  section) — required by
  `hep_cross/contracts/model_point_to_llp_recast_v1.yaml` before the point can
  be packaged for `llp_recast`, but out of scope for benchmark *selection*
  (needs MadGraph, explicitly a non-goal of this task).

## Files added

- `benchmarks/FIRST_H2_RECAST_CANDIDATE.json`
- `benchmarks/FIRST_H2_RECAST_CANDIDATE.md` (this file)
- `benchmarks/verify_pilot_ctau_invariant.py` (focused schema/provenance
  consistency check, run above; no new dependencies)

No existing files were modified; no C++ or physics code was touched; no scan
was executed.

## Unresolved decisions for PI input

1. Whether to authorize running the 15-point follow-up scan specified above.
2. Whether the `tan_beta` upper bound (`1e6`) and mass grid (`150–250 GeV`)
   are the right first region to probe, or whether a different suppression
   mechanism (e.g. moving off exact alignment, or a different mass splitting)
   should be tried first.
3. Whether `1 mm <= ctau_mm <= 1e4 mm` is the right operational LLP-scale
   window for "population of the recast selections," or whether the
   `llp_recast`/DV+jets validation contract already fixes a narrower
   acceptance window that should be used instead.
