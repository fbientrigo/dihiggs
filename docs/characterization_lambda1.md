# lambda1 construction golden characterization

Status: **characterization only** — this suite **freezes** the current behavior
of the lambda1-target construction path; it does **not** certify that behavior
as physically correct, and it repairs nothing. It exists so that the planned
selective migration (`SELECTIVE_MIGRATION_AUDIT.md` at the workspace root,
sections 7–8; migration PR 2) starts from a reproducible executable oracle
instead of prose contracts.

This is the `main_dihiggs` counterpart to boundary's
`docs/characterization_evaluate_point.md`. Boundary owns the physical-`M`
evaluator; this repository owns the **only** implemented `lambda1`-target
constructor, `THDM::set_param_phys_lam1`.

## Scope

Frozen here:

| Artifact | Contents |
|---|---|
| `tests/golden/lambda1_v1/cases.json` | 9 cases: CLI argument vectors plus per-case provenance and the behavior each one characterizes |
| `tests/golden/lambda1_v1/expected.csv` | the emitted CSV rows (7 rows from 9 cases — see row cardinality below), prefixed with `case_id` |
| `tests/golden/lambda1_v1/expected_markers.json` | per-case stdout/stderr observations: attempts, emitted rows, `TRIPLE_OK_POINTS`, warning counts |
| `tests/golden/lambda1_v1/manifest.json` | commit/dirty state, toolchain, dependency blob hashes, compiled constants, checksums, row counts |
| `tests/golden/lambda1_v1/mh_eps_sensitivity.json` | machine-readable output of the mh/EPS study |
| `tests/test_golden_lambda1.py` | comparator + reviewed assertions |
| `scripts/run_lambda1_characterization.sh` | the entry point |

Frozen path: `dihiggs/app/PhysLam1Scan` → `THDM::set_param_phys_lam1` →
`set_param_phys` → `get_param_gen` → `Constraints` → `DecayTable` →
`ParamUtils::write_csv_row`.

### Test vocabulary

`tests/test_golden_lambda1.py` names each test by the kind of claim it makes,
so the suite never conflates "this is what the code does" with "this is right":

- `behavior_*` — freezes current behavior; no correctness claim.
- `defect_*` — freezes a **real, characterized, unrepaired defect**. These tests
  pass by asserting the *broken* behavior and must be inverted when it is fixed.
- `invariant_*` — should hold for any correct implementation.
- `compat_*` — a property downstream consumers currently depend on.
- `opendecision_*` — an unresolved scientific/contract decision.

## Base selection and provenance

| | |
|---|---|
| Implementation base | `origin/main` @ `c021e0a1` |
| Local worktree HEAD at audit | `c090b83a` (0 ahead, 25 behind; merge-base == local HEAD) |
| Verdict | **NOT BLOCKED** — behaviorally equivalent for this path |

The local checkout is a strict ancestor of `origin/main`, and every
lambda1-critical file is **byte-identical** between them (verified by blob
hash): `2hdmc/src/THDM.{cpp,h}`, `dihiggs/src/PhysLam1Scan.cpp`,
`PhysScanWithFixings.cpp`, `M2PointEvaluator.cpp`, `2hdmc/src/Constraints.cpp`,
`2hdmc/src/DecayTable.cpp`. The 25 remote-only commits touch this path only
additively (a new `GenScanWithFixings` target and `gen_fixings` orchestrator
engine, an `mlpython` reorganization, `chris/` Stage-1 tools, a new
`conventions/physics_conventions.yaml`, and pinning of two floating
`GIT_TAG master` dependencies). The `dihiggs/Makefile` diff adds a target and
does not change `CXXFLAGS` or how `PhysLam1Scan` is built.

Per the base-selection rule, `origin/main` is therefore the implementation base
while the local HEAD retains provenance for existing campaigns. Both revisions
produce the same oracle for this path.

## How to run and regenerate

```bash
bash scripts/run_lambda1_characterization.sh          # build + run the suite
python3 scripts/lambda1_mh_eps_sensitivity.py         # reproduce the mh/EPS study
```

Regeneration is a deliberate, reviewed act — never a way to make a red test
green:

```bash
bash scripts/build_lambda1_characterization.sh
python3 scripts/generate_golden_lambda1.py
```

The generator runs **every case twice** and refuses to freeze output that is not
byte-identical across the two runs, so nondeterminism cannot be baked into the
oracle. Review the diff of all four golden files and record why the behavior
legitimately changed.

## Compiled constants (part of the oracle)

| Constant | Value | Note |
|---|---|---|
| `M_H` | **125.0** | `constexpr` in `dihiggs/src/PhysLam1Scan.cpp`. **Not** the paper's 125.09 |
| `THDM::EPS` | **1e-9** | `constexpr` in `2hdmc/src/THDM.h`. Stock 2HDMC uses `1e-12` |
| `yukawas_type` | 1 | hard-coded |
| `mHp` | `= mA` | `mA_fixed` is passed for both arguments |

Neither constant is changed by this work.

## Numerical comparison policy

The production CSV is written with `std::fixed << std::setprecision(15)` — 15
digits **after the decimal point**. This is an *absolute*-precision format, so
the comparison floor is absolute (`1e-15`) and identical for every column: two
values that serialize to the same 15-decimal string are indistinguishable in
this format, and a *relative* tolerance is meaningless for a column whose
serialized value is `0.000000000000000`.

Compared **exactly** (any difference is real drift): the three flags, and the
pass-through inputs and compile-time constants copied unmodified into the row
(`m_phi`, `mA`, `lambda6`, `lambda7`, `sin_ba`, `tan_beta`, `lam1`).

Repeated runs of the same binary must be **byte-identical**; the tolerance is
not slack for nondeterminism. Non-finite values are classified before any
tolerance arithmetic, and row shape is guarded explicitly so truncated or
overlong rows fail loudly instead of being absorbed by `zip`.

## What is frozen: findings

### 1. Construction failures are dropped (not preserved)

`PhysLam1Scan.cpp` executes `if (!pset) continue;`, so an attempted point that
fails construction **disappears entirely**. `L08` freezes this as an executable
fact: **3 attempted points → 2 emitted rows**. Boundary's contract is the
opposite (one row per attempted point, with `rejection_stage`/`reason`).

Worse, only one of the three failure classes is even visible:

| Case | Guard | Diagnostic |
|---|---|---|
| `L02` | `m_h > m_H` | stderr warning |
| `L03` | `tan_beta <= 0` | **completely silent** |
| `L04` | `abs(sba) > 1` | **completely silent** |

A theory *rejection*, by contrast, **is** emitted (`L05`) — so "absent from the
CSV" and "rejected by physics" are different things that this schema cannot
distinguish.

### 2. `std::fixed<<setprecision(15)` destroys the long-lived regime

This is the most consequential finding, because it lands exactly on the regime
the project exists to study.

For `L06` (small `lambda6`, large `tan_beta`):

| Quantity | True value (17 digits) | Emitted in CSV |
|---|---|---|
| `total_width` | `1.214398311035274e-17` GeV | `0.000000000000000` |
| `br_gaga` | `0.68512795134260707` | `0.000000000000000` |
| implied `ctau` | **16.25 m** | undefined (division by zero) |

Two independent mechanisms cause this:

1. **Serialization.** Fixed/15 cannot represent any width below `5e-16`; it
   writes zero. The smallest nonzero representable width is `1e-15` GeV.
2. **The `br_gaga` guard.** `br_gaga = (w_tot > 1e-15) ? w_gaga/w_tot : 0.0`
   forces the branching ratio to exactly `0` — and `w_tot = 1e-15` GeV
   corresponds to `ctau = ħc/Γ ≈ 0.197 m`, i.e. the guard zeroes `br_gaga` for
   any scalar with `ctau` above roughly **20 cm**, precisely the
   displaced-vertex regime `llp_recast` exists to study.

The consumer closes the loop: `autoresearch/harness/bounded_adaptive_search.py`
computes `ctau_m = HBARC_GEV_M / total_width` behind a `v > 0.0` filter (lines
277 and 792). A width that serialized to zero is therefore **silently dropped
from the ctau metric**.

Consequences, stated as characterization rather than repair:

- Maximum observable `ctau_m` ≈ **0.197 m** (one serialization quantum).
- Any point with true `ctau_m` **> ≈0.395 m** is invisible to the metric.
- `ctau_m` is quantized: the width can only be an integer multiple of `1e-15`.
- `L07` (the `GEMINI.md` "current best" point) has `total_width ≈ 8.006e-15` —
  serialized as `0.000000000000008`, i.e. **one significant digit**, so any
  `ctau` derived from it carries ~1 digit of precision.

The adaptive search's own objective is therefore capped and coarsened by a
serialization artifact, and its best candidates are the ones most likely to be
discarded. **Not repaired here** (mission stop condition: characterize existing
defects, do not fix them during characterization).

### 3. The lambda1 round-trip warning never rejects, and never reaches the CSV

`set_param_phys_lam1` recomputes lambda1 after construction, stores the residual
in `lam1_validation_*`, and warns on stderr when it exceeds `THDM::EPS`. But
`PhysLam1Scan` never reads those fields, so:

| Case | `|lam1 − computed_lam1|` | vs EPS=1e-9 | Outcome |
|---|---|---|---|
| `L01` | < 1e-12 | ok | silent |
| `L06` | `2.06e-09` | ~2× | stderr warning; **still triple-OK** |
| `L07` | `5.82e-07` | **~580×** | stderr warning; **still triple-OK** |

The `GEMINI.md` "current best state" point (`L07`) carries a lambda1
reconstruction error ~580× the fork's own threshold, is emitted with all flags
set, and is counted in `TRIPLE_OK_POINTS`. The warning exists only in stderr;
no CSV column records it. The residual *is* recomputable post-hoc from
`lam1` − `computed_lam1`, which is how the golden test asserts it.

### 4. Schema semantics

- The column named **`m12` actually holds m12 *squared*** (`m12_2_g` from
  `get_param_gen`); G06's value ≈ 334 GeV², not a mass. Audit section 4.4.
- **No lifetime column** is emitted (`ctau`/`ctau_m`/`ctau_mm` all absent), and
  the column is `m_phi`, not `mphi`. This contradicts
  `dihiggs/app/orchestrator/engines/lambda1.py:expected_csv_columns`, which
  claims both `mphi` and `ctau_m` — audit section 4.10 confirmed. That
  declaration is never executed against the real binary.
- **No `triple_ok`, `theory_ok` or `stability_ok` column exists.** `triple_ok`
  survives only as the stdout marker `TRIPLE_OK_POINTS`.
- `lam2` and `computed_lam2` are **the same value** (the source writes `lam2_g`
  twice, "no separate l2 calculation"), so the `target`/`computed` column
  pairing is asymmetric: `lam1` is a real target, `lam2` is not.

### 5. Stability is never evaluated, and positivity is aliased to it

`PhysLam1Scan` computes only positivity/unitarity/perturbativity; it never calls
`check_stability`. Independently, the **vendored fork carries the same alias
boundary has**:

```cpp
bool Constraints::check_positivity() { return model.check_stability(); }
bool Constraints::check_stability()  { return model.check_stability(); }
```

So the reported `positivity_ok` is really a stability result, and the alias is
**not** a boundary-only quirk — it is common to both repositories. Recorded as
implementation behavior, **not** canonical scientific truth. Whether positivity
and stability are intended synonyms remains a scientific-contract decision
(audit section 9, blocker 5); this suite does not resolve it.

## mh and EPS sensitivity

Reproduce with `python3 scripts/lambda1_mh_eps_sensitivity.py`. The study builds
probe variants **outside** the repository, rewriting exactly one constant each,
and verifies that an unmodified copy reproduces the production binary
byte-for-byte before drawing any conclusion. It is a narrow deterministic
comparison, not a scan.

### mh = 125.0 (compiled default) vs 125.09 (paper)

**Semantic inside a narrow band; numerical everywhere else.**

The construction guard is `m_h > m_H`, so raising `mh` *moves the boundary* and
turns previously-constructed points into (dropped, warned) construction
failures:

| `m_phi` | mh=125.0 | mh=125.09 | |
|---|---|---|---|
| `124.99999999999999` | fail | fail | |
| `125.0` | constructed | **CONSTRUCTION FAIL** | **semantic** |
| `125.05` | constructed | **CONSTRUCTION FAIL** | **semantic** |
| `125.09` | constructed | constructed | (boundary is strict `>`) |
| `125.2` | constructed | constructed | |

Away from that band, no flag flips on any probed point; the difference is
numerical: `m12²` moves by ~`8.2e-07` relative (G06) and ~`3.7e-07`
(`campaign_best`), with correspondingly small width/BR shifts. At exact
alignment with large `tan_beta` (`small_l6_large_tb`, `sin_ba = 1.0`) the `m12²`
delta is exactly **0** — the `m_h` term drops out of the reconstruction.

Practical consequence for migration: switching to the paper's 125.09 is not a
free relabel. It silently deletes attempted points in `m_phi ∈ [125.0, 125.09)`,
and because failures are dropped rather than recorded (finding 1), that deletion
leaves no trace in the CSV.

### THDM::EPS = 1e-9 (fork) vs 1e-12 (stock)

**No observable effect — diagnostic-only in everything sampled.**

EPS is used in four places. Three are cosmetic (`print_param_*` Z2-violation
warnings) or diagnostic (the lambda1 round-trip warning). The one that could
matter is inside `THDM::check_stability()`:

```cpp
if ((abs(lambda[6])<EPS)&&(abs(lambda[7])<EPS)) { /* shortcut branch */ }
```

Since `check_positivity` aliases to `check_stability` (finding 5), a change here
could in principle flip the reported `positivity_ok` for
`lambda6 ∈ [1e-12, 1e-9)` — the range where the two EPS values disagree about
which branch to take.

Executable result: **no difference was observed.** At `lambda6 = 1e-10` (which
does flip the branch), a 5016-row comparison across four configurations
(`mA ∈ {300,500}` × `tan_beta ∈ {50,10000}`, `m_phi` 130–500, `lambda1` −2…6)
produced **identical flags and identical rows**, and all three golden points are
**byte-identical** between the two builds — including the round-trip warning,
which fires or stays silent identically.

So on current evidence the fork's EPS change is **not semantic**: the shortcut
and full stability branches agree wherever they were compared. This is a bounded
empirical result over the sampled region, **not** a proof of global equivalence,
and it does not license removing or changing EPS. Audit blocker 2 ("the minimal
patch and its scientific effect must be isolated") is **narrowed, not closed**.

## Known limitations

- **Nine cases cannot certify a parameter space.** They pin the discrete
  decision structure (construct / drop / reject), the schema, and representative
  numerics of each class.
- **The oracle is the current implementation's output**, not independently
  verified theory values.
- **No golden case isolates unitarity or perturbativity alone.** `L05` freezes
  its full flag vector as-is; no coverage of a single-constraint failure is
  claimed. No synthetic point is fabricated to fill the gap.
- **Stability is uncovered by construction** — this path never evaluates it.
- The EPS conclusion is bounded by the sampled region (above).
- The suite says nothing about HB/HS, campaigns, MadGraph, the recast, or ML.

## UNRESOLVED / blocked

1. **The production build is broken from a clean checkout.**
   `make -C dihiggs app/PhysLam1Scan` fails with
   `fatal error: M2PointEvaluator.hpp: No such file or directory`, because
   `libdihiggs.a` compiles `COMMON_SRCS` (including `M2PointEvaluator.cpp`)
   while `dihiggs/include/*.hpp` is **untracked** — six M2 headers exist only in
   one developer's dirty worktree. `LDFLAGS` also unconditionally links
   `-lHiggsTools`, which `PhysLam1Scan.cpp` never uses. This suite therefore
   links `PhysLam1Scan` standalone against 2HDMC only — the same pattern the
   repository already applies to `GenScanWithFixings`, whose Makefile rule
   carries an equivalent comment. Audit blocker 10 is **confirmed**, not fixed.
   Until the headers are committed, `main` cannot be rebuilt from `main`.
2. **`mh = 125.0` vs the paper's `125.09`** — quantified above; the default is
   deliberately unchanged. The `m_phi ∈ [125.0, 125.09)` band must be an
   explicit migration decision.
3. **Whether the EPS fork change is required at all** — no effect observed, but
   the comparison is bounded.
4. **Whether the lambda1 round-trip residual should reject a point**, and
   whether `lam1_validation_*` belongs in the schema. Today a ~5.8e-07 residual
   is accepted silently.
5. **Whether `positivity` and `stability` are intended synonyms** (blocker 5,
   shared with boundary).
6. **Whether `triple_ok` (no stability) or boundary's stability-aware
   `theory_ok` is canonical.** They are not the same predicate.
7. **The Float64 loss in the long-lived regime** (finding 2) invalidates the
   `ctau` metric above ~0.2 m. Any migration that keeps this serialization
   inherits the ceiling; any migration that fixes it will **change historical
   campaign results**, so the fix needs an explicit, recorded decision.
