# First H2 LLP recast candidate — selection report

**The selected point is valid for the first `H2 -> bb` DV+jets recast because
its width, lifetime and `BRbb` are numerically stable. It remains a
provisional global benchmark because loop-induced channels and the
production normalization are unresolved.**

**Overall verdict: `PROVISIONAL_NUMERICAL_H2_BENCHMARK`.** No *existing*
post-Yukawa-correction H2 point in this repository satisfied the LLP-benchmark
eligibility gates (see "Part 1" below). A bounded 15-point follow-up scan
(Part 2) was then executed with the canonical corrected evaluator and found
six points with theory-valid, finite-width, LLP-scale lifetimes — but every
one of them carries `lambda1_residual_warning = 1`, so none qualifies as a
*clean* candidate. The best of the six, `H2scan_mH150_tb300000`, is recorded
as a provisional numerical candidate.

An independent numerical-stability check (see "Part 2 — Numerical-stability
check" below) has since been completed with a **channel-scoped** result, not
a single blanket verdict:

| Scope | Classification |
|---|---|
| Overall benchmark | `PROVISIONAL_NUMERICAL_H2_BENCHMARK` |
| `H2 -> bb` DV+jets recast inputs (`total_width_gev`, `ctau_mm`, `br_bb`) | **`VALID_FOR_FIRST_BB_RECAST`** |
| `H2 -> gamma gamma` | `NUMERICALLY_UNRESOLVED` |
| `H2 -> Z gamma` | `NUMERICALLY_UNRESOLVED` |

The `H2->bb` DV+jets recast consumes only `total_width_gev`, `ctau_mm`, and
`br_bb`, which are confirmed stable to ~0.0016% across the *entire*
theory-valid `m12_2` window at this point — well below the percent level. The
loop-induced `br_gammagamma`/`br_Zgamma` channels are not equally stable
(5.4%/9.6% variation across the same window) and remain
`NUMERICALLY_UNRESOLVED`; they are not used by the `bb` recast. The benchmark
stays **provisional at the global level** — not because the numerical warning
blocks the `bb` recast (it does not), but because `lambda1_residual_warning =
1`, the loop-induced BRs remain numerically sensitive, and the physical 2HDM
production cross section (`sigma_production_fb`) is still pending. This point
is **not** being relabeled as publication-ready or as a fully validated 2HDM
benchmark. See `benchmarks/FIRST_H2_RECAST_CANDIDATE.json` for the
machine-readable record (including the `channel_validity` fields) and
`benchmarks/first_h2_bounded_scan.csv` /
`first_h2_bounded_scan_manifest.json` for the full scan output.

- Base commit evaluated (Part 1, existing-data search): `27817ab156c23546117c93f1584dd4aa766f4850` (`main`)
- Scan producer commit (Part 2): `38c9cbd73d388e82ee4c4f3dcf8b25e4894ea151`, `evaluator_dirty=no`
- Yukawa-order fix commit (the validity boundary): `6bfad7662fd87750d838bf2fe0bd7ac00ee2326a`
  ("fix: initialize canonical evaluator Yukawas after construction", 2026-07-16 22:02:42 -0400)
- Branch: `analysis/select-first-h2-benchmark`, worked in the pre-existing worktree
  `/home/fabi/atlas_dihiggs/_worktrees/select-first-h2-benchmark`

# Part 1 — existing-data search (no valid candidate)

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
against every `width_ok=1` row in the pilot files and (after Part 2) the
bounded scan output:

```
$ python3 benchmarks/verify_pilot_ctau_invariant.py
docs/pilots/dihiggs_point_v2_v1/L01_accepted_anchor.csv: 1 width_ok=1 row(s) checked
docs/pilots/dihiggs_point_v2_v1/L05_theory_rejected_anchor.csv: 1 width_ok=1 row(s) checked
docs/pilots/dihiggs_point_v2_v1/L06_llp_anchor.csv: 1 width_ok=1 row(s) checked
docs/pilots/dihiggs_point_v2_v1/ordering_boundary.csv: 1 width_ok=1 row(s) checked
docs/pilots/dihiggs_point_v2_v1/construction_failure.csv: 0 width_ok=1 row(s) checked
docs/pilots/lambda1_v2_yukawa_fix_v1/lambda1_v2_pilot.csv: 4 width_ok=1 row(s) checked
benchmarks/first_h2_bounded_scan.csv: 15 width_ok=1 row(s) checked
PASS: ctau_mm == hbar_c_GeV_mm / total_width for all width_ok=1 pilot rows
```

All 23 checked rows satisfy the invariant to floating-point precision
(rel_err `0.0`), so "reproducible `ctau`" (one of the minimum validity
requirements) is confirmed for every row that reports a width, including `L01`
and every scan row.

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

# Part 2 — bounded follow-up scan (executed)

Part 1 found no valid existing candidate, so the smallest justified follow-up
scan was specified and then executed with the canonical corrected evaluator
(`dihiggs/app/Lambda1EvaluatorV2`, `dihiggs.lambda1.v2` schema), built locally
in this worktree from unmodified source
(`git diff --stat dihiggs/src/Lambda1EvaluatorV2.cpp dihiggs/src/ReplaySafeOutput.cpp`
is empty). No 2HDMC/HiggsTools physics code or Yukawa convention was changed.

## Scan specification (as run)

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
- **Runner:** `benchmarks/run_first_h2_bounded_scan.py`, resumable and
  interrupt-safe (rewrites the canonical CSV after every point, skips
  `point_id`s already present, refuses to mix rows from a different producer
  commit or dirty state).

## Scan results

All 15 grid points were evaluated (`stop_reason: exhausted_grid` — no *clean*
candidate was found, so the scan ran to completion rather than stopping
early). Every point is `construction_ok=1` and `theory_ok=1`
(positivity/unitarity/perturbativity all hold across this entire alignment-
limit, large-`tan_beta` region); every point has `width_ok=1` and
`br_bb > 0.5`. The only thing separating "provisional" from "rejected" here is
the `1 mm <= ctau_mm <= 1e4 mm` window and `lambda1_residual_warning`:

| point_id | m_H2 (GeV) | tan β | ctau_mm | br_bb | lambda1_abs_residual | lambda1_residual_warning | classification |
|---|---:|---:|---:|---:|---:|---:|---|
| H2scan_mH150_tb300000 | 150 | 3×10⁵ | 4.326 | 0.7567 | 9.342e-07 | 1 | **provisional** |
| H2scan_mH200_tb300000 | 200 | 3×10⁵ | 3.058 | 0.6755 | 9.342e-07 | 1 | provisional |
| H2scan_mH250_tb300000 | 250 | 3×10⁵ | 2.181 | 0.5776 | 9.342e-07 | 1 | provisional |
| H2scan_mH150_tb100000 | 150 | 1×10⁵ | 0.481 | 0.7567 | 2.700e-07 | 1 | rejected (ctau below 1 mm) |
| H2scan_mH200_tb100000 | 200 | 1×10⁵ | 0.340 | 0.6755 | 2.700e-07 | 1 | rejected (ctau below 1 mm) |
| H2scan_mH250_tb100000 | 250 | 1×10⁵ | 0.242 | 0.5776 | 2.700e-07 | 1 | rejected (ctau below 1 mm) |
| H2scan_mH150_tb1000000 | 150 | 1×10⁶ | 48.070 | 0.7567 | 8.750e-06 | 1 | provisional |
| H2scan_mH200_tb1000000 | 200 | 1×10⁶ | 33.976 | 0.6755 | 5.126e-05 | 1 | provisional |
| H2scan_mH250_tb1000000 | 250 | 1×10⁶ | 24.231 | 0.5776 | 5.126e-05 | 1 | provisional |
| H2scan_mH150_tb30000 | 150 | 3×10⁴ | 0.0433 | 0.7567 | 1.160e-08 | 1 | rejected (ctau below 1 mm) |
| H2scan_mH200_tb30000 | 200 | 3×10⁴ | 0.0306 | 0.6755 | 4.241e-08 | 1 | rejected (ctau below 1 mm) |
| H2scan_mH250_tb30000 | 250 | 3×10⁴ | 0.0218 | 0.5776 | 4.241e-08 | 1 | rejected (ctau below 1 mm) |
| H2scan_mH150_tb10000 | 150 | 1×10⁴ | 0.00481 | 0.7567 | 8.062e-09 | 1 | rejected (ctau below 1 mm) |
| H2scan_mH200_tb10000 | 200 | 1×10⁴ | 0.00340 | 0.6755 | 2.061e-09 | 1 | rejected (ctau below 1 mm) |
| H2scan_mH250_tb10000 | 250 | 1×10⁴ | 0.00242 | 0.5776 | 2.061e-09 | 1 | rejected (ctau below 1 mm) |

Full per-row provenance (all input parameters, reconstructed λ's, all nine
partial widths and BRs, producer commit/dirty state) is in
`benchmarks/first_h2_bounded_scan.csv`; run parameters and classification
counts (6 provisional, 0 clean, 9 rejected) are in
`benchmarks/first_h2_bounded_scan_manifest.json`.

**Lifetime scaling.** The predicted `1/tan^2(beta)` scaling from the scan
specification holds well: going from `tan_beta=1e4` to `3e4` (factor 3) moves
`ctau_mm` by very close to `3^2=9`x at fixed mass (e.g. `mH=200`:
`3.397e-3 -> 3.058e-2`, a 9.0x increase); `1e4 -> 1e5` (factor 10) gives
`3.397e-3 -> 0.340`, a 100x increase; `1e4 -> 3e5` gives `3.397e-3 -> 3.058`,
a 900x increase; `1e4 -> 1e6` gives `3.397e-3 -> 33.98`, a ~10,000x increase —
all within the predicted `tan_beta^2` law. The earlier estimate
(`tan_beta ~ 1.7e5` for `ctau=1mm`) is confirmed to fall between the scanned
`1e5` (still sub-mm, rejected) and `3e5` (already 2-4 mm) grid points.

**Numerical stability.** `lambda1_abs_residual` grows sharply with `tan_beta`:
`~2e-9` at `1e4`, `~1e-8` at `3e4`, `~3e-7` at `1e5`, `~9e-7` at `3e5`, and
`5e-5`-`9e-6` at `1e6`. All 15 points exceed the evaluator's built-in
`1e-9` warning threshold, so **every** point in the entire grid — not just the
six lifetime-window points — carries `lambda1_residual_warning=1`. This
confirms the concern flagged in the original scan specification: pushing
`tan_beta` this high measurably degrades the `set_param_phys_lam1` round-trip
(`m12_sq`/`lambda1` reconstruction), and no point in this grid is clean by the
numerical-stability criterion alone.

## Provisional candidate: `H2scan_mH150_tb300000`

Selected as the best of the six lifetime-window points per the task's ranking
order (LLP-plausible lifetime; DV+jets-compatible `H2->bb`; plausible mass;
complete provenance; stable numerics; minimal extra work): it is the first
point reached in scan priority order, has the highest `br_bb` (0.757) and the
smallest `lambda1_abs_residual` (9.34e-07) among the three `tan_beta=3e5`
points, and its `ctau_mm=4.33` sits comfortably inside the 1–1e4 mm window
rather than near either edge.

| Field | Value |
|---|---|
| `m_H2_GeV` (`mH_input_gev`) | `150.0` |
| `mA_input_gev` / `mHp_input_gev` | `450.0` / `450.0` |
| `sin_beta_minus_alpha_input` | `1.0` |
| `tan_beta_input` | `3.0e5` |
| `lambda1_target` | `1.0` |
| `lambda6_input` / `lambda7_input` | `1.0e-10` / `0.0` |
| `construction_ok` / `theory_ok` (`triple_ok`) | `1` / `1` |
| `yukawa_type_installed` | `1` (Type I) |
| `total_width_gev` | `4.56118529862185007e-14` |
| `ctau_mm` | `4.32622152973311191e+00` |
| `br_bb` | `7.56737485808578692e-01` |
| `br_gammagamma` | `2.65792941770435500e-04` |
| `br_Zgamma` | `2.04915508605281493e-05` |
| `br_tautau` / `br_cc` / `br_gg` | `0.0757` / `0.0339` / `0.1329` |
| `br_WW` / `br_ZZ` | `0.0` / `0.0` (vanish at exact alignment) |
| `lambda1_abs_residual` | `9.34182041500974947e-07` |
| `lambda1_residual_warning` | `1` |
| `sigma_production_fb` | `pending` — not computed; requires MadGraph/UFO (out of scope here) |
| `producer_commit` / `evaluator_dirty` | `38c9cbd73d388e82ee4c4f3dcf8b25e4894ea151` / `no` |
| `evaluator_schema` / `evaluator_api` | `dihiggs.lambda1.v2` / `THDM::set_param_phys_lam1+2HDMC::DecayTable` |
| Source row | `benchmarks/first_h2_bounded_scan.csv`, `point_id=H2scan_mH150_tb300000` |

**Why physically interesting:** it is a Type I, exact-alignment (`sin(β-α)=1`)
H2 with a detector-scale lifetime (`ctau_mm ≈ 4.3 mm`) and a dominant,
DV+jets-plausible `H2 -> bb` branching ratio (0.757), reached by the same
`cot^2(beta)` suppression mechanism already validated on `L06` but pushed far
enough in `tan_beta` to survive the post-Yukawa-fix corrected widths. Unlike
`L06`, its post-fix lifetime is genuinely in the LLP range rather than having
collapsed to micron scale.

**Why it is not cleanly accepted:** `lambda1_residual_warning=1` with
`lambda1_abs_residual=9.34e-07`, roughly 900x the evaluator's `1e-9` warning
threshold. This flags degraded numerical accuracy in the
`set_param_phys_lam1` round-trip (the reconstructed `m12_sq`/`lambda1` do not
close to the same precision as at moderate `tan_beta`), not a theory or
construction failure — `theory_ok=1` and `construction_ok=1` are unaffected.
Promoting a point with this residual to "clean" would risk overstating the
precision of its width/lifetime/BR values downstream.

## Numerical-stability check: channel-scoped classification

The recheck flagged above has been completed. Method: `THDM::set_param_phys`
was called directly (bypassing `set_param_phys_lam1` and its lambda1
round-trip entirely) with an `m12_2` recomputed in double precision from
2HDMC's own closed-form lambda1-inversion formula
(`2hdmc/src/THDM.cpp::set_param_phys_lam1`), via a small standalone C++
program, `benchmarks/check_H2scan_mH150_tb300000.cpp`, driven by
`benchmarks/check_H2scan_mH150_tb300000.py`. This reproduces the canonical
row bit-for-bit (0.0 relative difference) when fed the identical `m12_2` —
confirming the reimplementation is faithful — then probes `m12_2` at offsets
of order `1e-12 GeV^2` (the scale set by the analytic sensitivity
`d(lambda1)/d(m12_2) = tan_beta/(v^2 cos^2(beta)) = 4.4536e11 GeV^-2` at this
point) to map the theory-valid window and the stability of the physical
outputs within it. No 2HDMC source was modified, the Yukawa convention was
unchanged, and no other mass or `tan_beta` point was scanned.

Findings:

- The theory-valid (`positivity_ok`/`unitarity_ok`/`perturbativity_ok`)
  window in `m12_2` at this exact point is only a few times `1e-12 GeV^2`
  wide; `lambda1_reconstructed` itself ranges from about 0.55 to 1.45 across
  that window (it is not pinned to 1.0 by theory-validity alone) — this,
  amplified by the huge sensitivity coefficient above, is the root cause of
  `lambda1_residual_warning=1`, and is not a fixable bug in
  `set_param_phys_lam1`.
- Within that theory-valid window, **`total_width_gev`, `ctau_mm`, `br_bb`,
  `br_tautau`, and `br_gg` vary by only ~0.0016%** — they are dominated by
  tree-level Yukawa-driven fermionic widths, which do not depend on
  `lambda1`/`lambda3`.
- **`br_gammagamma` (5.42%) and `br_Zgamma` (9.59%) vary by more than the 5%
  `NUMERICALLY_UNRESOLVED` threshold** across the same window — these
  loop-induced, charged-Higgs-mediated channels are directly sensitive to
  `lambda3`, which shares `lambda1`'s `m12_2`-linear dependence.

**Unstable quantities: `br_gammagamma`, `br_Zgamma` only.** `total_width_gev`,
`ctau_mm`, and `br_bb` — the only quantities the proposed `H2->bb` DV+jets
recast consumes — are independently confirmed numerically robust at this
point, regardless of construction path. See
`benchmarks/H2scan_mH150_tb300000_numerical_check.md` for the full derivation
and probe table, and `benchmarks/FIRST_H2_RECAST_CANDIDATE.json`'s
`numerical_stability_check` object for the machine-readable record.

The channel-scoped result of this check:

| Scope | Classification | Basis |
|---|---|---|
| Overall benchmark | `PROVISIONAL_NUMERICAL_H2_BENCHMARK` | `lambda1_residual_warning=1`; loop-induced BRs numerically sensitive; `sigma_production_fb` pending |
| `H2 -> bb` DV+jets recast inputs | `VALID_FOR_FIRST_BB_RECAST` | `total_width_gev`/`ctau_mm`/`br_bb` stable to ~0.0016% across the entire theory-valid `m12_2` window |
| `H2 -> gamma gamma` | `NUMERICALLY_UNRESOLVED` | `br_gammagamma` varies 5.42% across the same window |
| `H2 -> Z gamma` | `NUMERICALLY_UNRESOLVED` | `br_Zgamma` varies 9.59% across the same window |

The candidate is **not** promoted to a fully validated, publication-ready
2HDM benchmark; the global status remains `PROVISIONAL_NUMERICAL_H2_BENCHMARK`.
That global caveat does **not** block the first `H2->bb` DV+jets recast: its
inputs (`total_width_gev`/`ctau_mm`/`br_bb`) are independently confirmed
numerically robust and are classified `VALID_FOR_FIRST_BB_RECAST` above.
Downstream work that depends on `br_gammagamma`/`br_Zgamma` should not
proceed on this point until those channels are resolved; downstream work
that depends only on `total_width_gev`/`ctau_mm`/`br_bb` (i.e. the `bb`
recast) is not blocked by that open item.

## Next required step

1. Proceed with the first `H2->bb` DV+jets recast using
   `total_width_gev`/`ctau_mm`/`br_bb` from this point — these inputs are
   `VALID_FOR_FIRST_BB_RECAST` (see channel-scoped table above).
2. Generate `sigma_production_fb` via MadGraph/UFO for this exact
   `(m_H2, tan_beta, sin(β-α), Type I)` point and package the result against
   `hep_cross/contracts/model_point_to_llp_recast_v1.yaml`. Not performed here
   (out of scope for this task).
3. Separately, PI decision on whether/how to pursue resolving the
   `br_gammagamma`/`br_Zgamma` instability before any downstream work that
   depends on those loop-induced channels, or a fully validated,
   publication-ready benchmark (see "Unresolved decisions" below).

## Files added

- `benchmarks/FIRST_H2_RECAST_CANDIDATE.json`
- `benchmarks/FIRST_H2_RECAST_CANDIDATE.md` (this file)
- `benchmarks/verify_pilot_ctau_invariant.py` (extended to also check the new
  scan CSV)
- `benchmarks/run_first_h2_bounded_scan.py` (deterministic, resumable 15-point
  scan runner)
- `benchmarks/first_h2_bounded_scan.csv` (full per-point evaluator output, 15
  rows)
- `benchmarks/first_h2_bounded_scan_manifest.json` (grid declaration,
  producer provenance, classification counts, stop reason)
- `benchmarks/check_H2scan_mH150_tb300000.cpp` (one-use standalone C++
  program: constructs the model via `THDM::set_param_phys` directly, given an
  externally supplied `m12_2`)
- `benchmarks/check_H2scan_mH150_tb300000.py` (one-use driver: recomputes
  `m12_2` from 2HDMC's closed-form formula, compiles/runs the `.cpp` checker,
  probes the theory-valid window, writes the two artifacts below)
- `benchmarks/H2scan_mH150_tb300000_numerical_check.csv` /
  `.md` (numerical-stability check results: reproduction check, sensitivity
  probe table, classification)

No existing C++ or physics code was touched; no Yukawa convention changed; no
MadGraph/Pythia generation was run; `dihiggs/app/Lambda1EvaluatorV2` and
`2hdmc/lib/lib2HDMC.a` were built locally from unmodified source (build
artifacts, gitignored, not committed).

## Unresolved decisions for PI input

1. ~~Whether the `br_gammagamma`/`br_Zgamma` instability blocks use of
   `total_width_gev`/`ctau_mm`/`br_bb` for the first recast.~~ **Resolved in
   this update:** it does not. The classification is now channel-scoped
   (see table above) — `bb_dvjets` is `VALID_FOR_FIRST_BB_RECAST` and is not
   gated by the `gammagamma`/`Zgamma` instability. What remains open for PI
   input is whether/how to pursue resolving `br_gammagamma`/`br_Zgamma`
   before any downstream work that specifically depends on those loop-induced
   channels, or before this point could be considered a fully validated,
   publication-ready 2HDM benchmark.
2. Whether to search other suppression mechanisms (e.g. moving slightly off
   exact alignment, or a different mass splitting) instead of pushing
   `tan_beta` further, given that **every** point in this 15-point grid
   already exceeds the evaluator's residual-warning threshold, and that the
   underlying cause (a ~1e-12 GeV^2-wide theory-valid `m12_2` window at large
   `tan_beta`) is intrinsic to the `lambda1`-input parameterization at this
   corner of parameter space, not a fixable implementation bug.
3. Whether `1 mm <= ctau_mm <= 1e4 mm` is the right operational LLP-scale
   window for "population of the recast selections," or whether the
   `llp_recast`/DV+jets validation contract already fixes a narrower
   acceptance window that should be used instead.
