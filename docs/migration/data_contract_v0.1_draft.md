# Canonical DiHiggs point-data contract — v0.1 DRAFT

Status: **DRAFT, not authoritative.** This is the candidate canonical data
contract for `dihiggs-core-v2`. It is derived from **executable evidence** —
two compiled golden oracles and the actual consumer source — not from prose.
Nothing here changes any producer's behavior.

Machine-readable companion: [`field_crosswalk_v0.1.json`](field_crosswalk_v0.1.json)
(34 field entries; some cover a column family).

## Why v0.1 and not v0.0

`v0.0` was reserved for the case where a material **local-vs-origin physics
divergence** in `main_dihiggs` remained unresolved. It did not:

- local HEAD `c090b83a` is a **strict ancestor** of `origin/main` `c021e0a1`
  (0 ahead, 25 behind, merge-base == local HEAD);
- every lambda1-critical file is **byte-identical** between them (blob-hash
  verified: `THDM.{cpp,h}`, `PhysLam1Scan.cpp`, `PhysScanWithFixings.cpp`,
  `M2PointEvaluator.cpp`, `Constraints.cpp`, `DecayTable.cpp`);
- the 25 remote-only commits touch this path only additively.

So no field is blocked by a local-vs-remote divergence. Fields **are** blocked
by genuine cross-repository and scientific divergences; those are marked
`UNRESOLVED` and listed below.

## Evidence base

| Producer | Oracle | Commit |
|---|---|---|
| `boundary/src/evaluate_point.cpp` | `tests/golden/evaluate_point_v1/expected.csv`, SHA-256 `513cc9de…` | boundary `main` `ae9bb132` |
| `main_dihiggs` lambda1 path | `tests/golden/lambda1_v1/expected.csv` | `origin/main` `c021e0a1` |

Consumers inspected as source, not documentation:
`autoresearch/harness/bounded_adaptive_search.py`,
`autoresearch/harness/dihiggs_validators.py`,
`mlpython/lake_pipeline/consolidate_lake.py`,
`dihiggs/app/orchestrator/engines/lambda1.py`,
`hep_cross/contracts/model_point_to_llp_recast_v1.yaml`, `llp_recast`.

`conventions/physics_conventions.yaml` (`physics_conventions_v1`, present on
`origin/main`) already fixes `hbar_c_gev_mm = 1.973269804e-13`,
`c_mm_per_ns = 299.792458`, and the `h1/h2/h3/hc` scalar naming with PDG ids and
mass columns. **This contract extends that file; it does not contradict it.**

## Canonical naming rules

1. lowercase `snake_case`.
2. Explicit units wherever ambiguity is plausible (`_GeV`, `_GeV2`, `_mm`, `_fb`).
3. **No silent metre/millimetre conversion.** `ctau_mm_H2` is canonical;
   `ctau_m` and `tau_ns` are explicit conversion outputs, never unitless aliases.
4. **Float64** is the canonical persisted scientific type; CSV uses round-trip
   safe 17-significant-digit formatting.
5. Booleans are `0/1` with exact documented semantics — never floats.
6. Stable `point_id` on every record; explicit `schema_version`.
7. **Input, requested target and reconstructed value are distinct fields**
   whenever they can differ (`lambda1_target` vs `lambda1_reconstructed`;
   `m12_sq_input` vs `m12_sq_derived`).
8. Names are not changed for aesthetics — only to remove ambiguity or reduce
   migration cost.

## Theory semantics

Defined, with the current implementation behavior recorded honestly:

| Predicate | Definition | Status |
|---|---|---|
| `construction_ok` | the 2HDMC parameter-setting call returned true | **UNRESOLVED in main** — main cannot express `0`: it *drops* the row |
| `positivity_ok` | `Constraints::check_positivity()` | **UNRESOLVED** — aliased to stability (below) |
| `unitarity_ok` | `Constraints::check_unitarity()` | defined; no golden isolates it |
| `perturbativity_ok` | `Constraints::check_perturbativity()` | defined; boundary G04 isolates it |
| `stability_ok` | vacuum stability | **UNRESOLVED** — aliased, and never evaluated at all in main |
| `triple_ok` | `positivity_ok && unitarity_ok && perturbativity_ok` | defined; **not** paper-valid |
| `theory_ok` | `construction_ok && triple_ok && stability_ok && width_ok` (boundary's expression) | **UNRESOLVED** — main's M2 evaluator means something else |
| `width_ok` | `total_width_H2_GeV` finite and `> 0` | **UNRESOLVED** — uncomputable from main data in the long-lived regime |

### The positivity/stability alias (implementation behavior, not truth)

In **both** repositories' 2HDMC — boundary's stock tree *and* main's fork —

```cpp
bool Constraints::check_positivity() { return model.check_stability(); }
bool Constraints::check_stability()  { return model.check_stability(); }
```

So `positivity_ok == stability_ok` identically on every constructed row, and the
`stability_ok` term in boundary's `theory_ok` is currently **redundant**. This
is recorded as **implementation behavior, not canonical scientific truth.** It
is *not* a boundary-only quirk — the audit noted it for stock 2HDMC; this work
confirms the fork carries it too.

This draft therefore **does not define independent stability semantics**: there
is no executable implementation to support them, and audit blocker 5 makes it a
scientific decision. No golden point with `triple_ok=1, stability_ok=0` can
exist under either dependency; that combination is exercised only by a clearly
labeled synthetic boolean-logic test in boundary.

**`triple_ok` is not a synonym for paper-valid.** It omits stability and width
validity. No downstream filter may treat it as acceptance merely because the
current dependency makes the extra term redundant.

## Summary of divergences

| # | Field | boundary | main lambda1 | Class | Human decision? |
|---|---|---|---|---|---|
| 1 | `mh_GeV` | 125.09, emitted | **125.0**, not emitted | UNRESOLVED | **Yes** |
| 2 | `mH_GeV` | `mH` | `m_phi` | RENAME | No |
| 3 | `m12_sq_GeV2` | `m12_sq_input`/`_derived` | **`m12`** (holds m12²) | SEMANTIC_DRIFT | No |
| 4 | `lambda1_reconstructed` | **`lambda1`** | `computed_lam1` | SEMANTIC_DRIFT | No |
| 5 | `lambda1_target` | — | `lam1` | RENAME | No |
| 6 | `construction_ok` | `set_param_phys_ok`, row kept | **row dropped** | SEMANTIC_DRIFT | No |
| 7 | `stability_ok` | emitted (aliased) | **never evaluated** | UNRESOLVED | **Yes** |
| 8 | `theory_ok` | stability+width aware | different predicate | UNRESOLVED | **Yes** |
| 9 | `total_width_H2_GeV` | `total_width_H2` | `total_width`, **fixed/15** | UNRESOLVED | **Yes** |
| 10 | `ctau_mm_H2` | `ctau_mm_H2` (mm) | **not emitted**; autoresearch derives **metres** | UNIT_CONVERSION | No |
| 11 | `br_gammagamma_H2` | `br_gammagamma_H2` | `br_gaga`, **zeroed by a guard** | UNRESOLVED | **Yes** |
| 12 | `point_id` | `point_id` | **none** | UNRESOLVED | No |
| 13 | `mHp_GeV` | emitted | **implicit** (`mA` passed twice) | SEMANTIC_DRIFT | No |

### The three that matter most

**A. `total_width_H2_GeV` — the long-lived regime is destroyed.** main serializes
with `std::fixed << std::setprecision(15)`, so any width below `5e-16` GeV
becomes exactly `0.000000000000000`. Golden case `L06` has a **true** width of
`1.214398311035274e-17` GeV — a **true ctau of 16.2 m**, an excellent LLP
candidate — and emits `0`. `autoresearch` computes
`ctau_m = HBARC_GEV_M / total_width` behind a `v > 0.0` filter, so the point is
**silently dropped from the ctau metric**. Consequences:

- maximum observable `ctau_m` ≈ **0.197 m** (one serialization quantum);
- any true `ctau_m` **> ≈0.395 m** is invisible;
- `ctau_m` is quantized to `ħc/(k·1e-15)`;
- `L07` (the `GEMINI.md` "current best" point) serializes as
  `0.000000000000008` — **one significant digit**.

The adaptive search's own objective is capped and coarsened by a serialization
artifact, and its **best candidates are the ones most likely to disappear.**
`mlpython/lake_pipeline` then casts `total_width` to **Float32** on top of this.

Fixing this **will change historical campaign results** — hence a recorded
decision, not a silent repair.

**B. `br_gammagamma_H2` — the displaced regime is zeroed.** main computes
`br_gaga = (w_tot > 1e-15) ? w_gaga/w_tot : 0.0`. A width of `1e-15` GeV is
`ctau ≈ 0.197 m`, so the guard forces `br_gaga` to **exactly 0** for any scalar
with `ctau` above roughly **20 cm** — precisely the displaced-vertex regime
`llp_recast` exists to study. `L06`'s true `br_gaga` is `0.685`; the emitted
value is `0`, indistinguishable from a genuine zero.

**C. `mh_GeV` — 125.0 vs 125.09 is not a free relabel.** The construction guard
is `m_h > m_H`, so raising `mh` *moves the boundary*
(`scripts/lambda1_mh_eps_sensitivity.py`):

| `m_phi` | mh=125.0 | mh=125.09 |
|---|---|---|
| `125.0` | constructed | **CONSTRUCTION FAIL** |
| `125.05` | constructed | **CONSTRUCTION FAIL** |
| `125.09` | constructed | constructed |

Inside `m_phi ∈ [125.0, 125.09)` the change is **semantic**; elsewhere it is
numerical (`m12²` moves ~`8e-7` relative). Because failures are *dropped*
(divergence 6), that deletion would leave **no trace** in the CSV.

## Fields blocked on human scientific input

1. **`mh_GeV`** — 125.0 vs 125.09, and what happens to `m_phi ∈ [125.0, 125.09)`.
2. **`stability_ok` / `positivity_ok`** — synonyms, or is a distinct check
   required? (audit blocker 5)
3. **`theory_ok`** — which predicate is canonical, and its version.
4. **`total_width_H2_GeV` / `br_gammagamma_H2`** — accepting the fix means
   historical campaign results change; the ctau ceiling and the `1e-15` guard
   must be explicitly retired or retained.
5. **Paper coordinate** — `M` (boundary) vs `lambda1_target` (main). Core can
   support both `construction_mode`s, but the frozen paper dataset must choose
   and document one sampling measure (audit blocker 1).
6. **lambda1 round-trip residual** — should ~`5.8e-07` (≈580× EPS) reject a
   point? Today it is accepted silently and the warning never reaches the CSV.

## Known non-blocking gaps

- No golden point in **either** repository isolates a unitarity-only failure.
- `dihiggs/app/orchestrator/engines/lambda1.py:expected_csv_columns` declares
  `mphi` and `ctau_m`; the binary emits `m_phi` and **no lifetime column**. The
  declaration is **never executed** against the real binary (audit 4.10,
  confirmed).
- `main` cannot be rebuilt from `main`: `dihiggs/include/*.hpp` is untracked, so
  `make -C dihiggs app/PhysLam1Scan` fails on a clean checkout (audit
  blocker 10, confirmed).
- `THDM::EPS` `1e-9` (fork) vs `1e-12` (stock) showed **no observable
  difference** over the sampled region — a bounded empirical result, not a proof
  of global equivalence. Audit blocker 2 is narrowed, not closed.

## Next step

Migration PR 3 (core skeleton/schema) should consume this crosswalk and the two
golden oracles directly, generating the C++ and Python record definitions from a
single source. The `UNRESOLVED` fields above must be resolved by their owners
**before** any frozen paper dataset is published — not during it.
