# Inter-repo communication: assessment & highest-ROI changes

_Assessment of how the three dihiggs repos interoperate, the maintainability
risks in that communication, and the changes made to harden it._

## The three repos and how they talk

| Repo | Pipeline role | Produces | Consumes |
|------|---------------|----------|----------|
| `dihiggs` | parameter scan + data lake | scan CSV/parquet, lake | — |
| `dihiggs_boundary` | theory + experiment boundary (2HDMC + HiggsTools) | `evaluate_point.csv` → `hbhs_enriched.csv` → `boundary_atlas.csv` | scan points |
| `dihiggs_hep_cross` | LLP recast / paper cross-checks | recast yields, figures | model-point CSV |

They are **fully decoupled at the code level**: no git submodules, no
cross-repo imports, no shared package. They communicate through **CSV files**
carried between stages, described by **prose "contracts"** in
`docs/contracts/*.md` (hep_cross) and `docs/*_contract.md` (boundary).

That decoupling is good for isolation. The risk is in the *communication
channel itself*: the CSV column schemas and the physics constants that must
agree across repos were maintained by hand, in prose and in copy-pasted
literals, with almost nothing checking that the sides still agree.

## Findings, ranked by ROI

| # | Finding | Impact | Cost | ROI |
|---|---------|--------|------|-----|
| 1 | CSV contracts were prose-only (except `dhb.schema`); no validator ran at a stage boundary | High — a renamed/dropped column produces silent wrong science | Low | **Highest** |
| 2 | The physics constant `HBAR_C_GEV_MM = 1.973269804e-13` and the `ctau = ħc/Γ` formula were copy-pasted into **8+ files across 2 repos** | High — a divergent edit means inconsistent lifetimes | Low | **Highest** |
| 3 | The C++↔Python HBHS schema sync was guarded only by a static fixture | Med-High | Low | High |
| 4 | `dihiggs_boundary` and `dihiggs_hep_cross` had **no CI** — their good tests never ran | Med-High | Low | High |
| 5 | Column-name drift across the boundary: `ctau_mm` (hep_cross) vs `ctau_mm_H2` (boundary), same physics | Med | Low | High |
| 6 | Main `dihiggs` repo is not an importable package, which is *why* constants got copy-pasted | Med | Med | Med |
| 7 | `2HDMC 1.8` and `HiggsTools v1.2` vendored twice, no shared version pin | Med | Med-High | Med |

## What changed (Phases A–C)

### Phase A — machine-readable, validated CSV contracts (findings 1, 3, 5)

- **`dihiggs_boundary`**
  - `python/dhb/contracts.py` — single source of truth for the three
    cross-stage column sets, built from `dhb.schema`. Emits a committed YAML
    mirror under `contracts/*.yaml` (`python -m dhb.contracts --emit`).
  - `python/dhb/validate.py` — reusable validator (`validate_csv` + a CLI,
    `python -m dhb.validate --contract … --input …`) that checks required
    columns **and** the numeric invariant `ctau_mm_H2 == ħc / total_width_H2`.
  - `enrich.py` / `atlas.py` now take their required-column lists from
    `dhb.contracts` instead of inline copies (with an assertion tying the two
    together), so the pipeline code and the contract cannot drift.
  - `tests/test_schema.py` gained a **non-circular** C++/Python drift guard:
    it reads the name arrays (`kNeutralNames`, `kEffFermionNames`,
    `kHhPairNames`) straight out of `src/evaluate_point.cpp` and asserts they
    equal the `dhb.schema` constants — catching a rename on either side
    without compiling the evaluator.
- **`dihiggs_hep_cross`**
  - `src/llp_recast/data_contract.py` — enforceable version of
    `docs/contracts/model_point_to_llp_recast_contract.md`: required columns,
    the `ctau_mm == ħc / total_width_GeV` rule the prose contract demands, the
    `BR_* ≤ 1` sum rule, and the allowed `beta_gamma_source` enum. Emits a
    committed `contracts/model_point_to_llp_recast_v1.yaml`.
- **Column-name drift (finding 5)** is recorded as an explicit `aliases`
  block in both contracts (`ctau_mm ↔ ctau_mm_H2`,
  `total_width_GeV ↔ total_width_H2`, `m_scalar_GeV ↔ mS_GeV`), so a future
  consumer can map the two naming schemes deterministically.

### Phase B — one values-only source for physics constants (findings 2, 6)

- `conventions/physics_conventions.yaml` — **byte-identical in all three
  repos** (md5 pinned; see the file header) — holds `hbar_c_gev_mm`,
  `c_mm_per_ns`, and the neutral/charged scalar naming + PDG map.
- Each repo reads it through a thin loader with a pinned-literal fallback
  (`dihiggs/mlpython/lake_pipeline/physics_conventions.py`,
  `dihiggs_hep_cross/src/llp_recast/constants.py`,
  `dihiggs_boundary/python/dhb/contracts.HBAR_C_GEV_MM`), and a test in each
  repo asserts the loaded value equals the pinned literal.
- The `HBAR_C_GEV_MM` literal was removed from **six** `dihiggs`
  lake-pipeline scripts (`extract_ctau_population.py`,
  `paper_like_three_br_ctau_lambda1.py`, `ctau_population_hist.py`,
  `make_focused_grid_slices.py`, `ctau_landscape_scan.py`,
  `make_br_ctau_panel_v2.py`), which now import the shared constant. No
  forced cross-repo package dependency was introduced (lightweight coupling).

### Phase C — CI for the previously untested repos (finding 4)

- `dihiggs_boundary/.github/workflows/ci.yml` and
  `dihiggs_hep_cross/.github/workflows/ci.yml` run `pytest` on push/PR, plus a
  step that re-emits each contract and fails if the committed YAML drifted
  from the code, plus (boundary) a validator smoke against the fixture. These
  need no 2HDMC/HiggsTools build.

## Phase D — external-tool version pinning & provenance (finding 7)

The original idea — collapse the twice-vendored `2HDMC` / `HiggsTools` trees
into one shared submodule — turned out to be **unsafe on inspection**:

- **`dihiggs/2hdmc` is a patched fork, not stock v1.8.** It adds a whole
  `set_param_phys_lam1()` API (λ1-basis input, m12² reconstruction +
  validation) and custom `Calc*` executables that the `dihiggs_boundary` copy
  lacks; the dihiggs λ1-scan workflow depends on them. The two copies are
  legitimately different code and must not be merged.
- **`dihiggs/higgstools` core is identical to boundary's**, except its
  `external/CMakeLists.txt` floated `json` and `eigen3` on `GIT_TAG master` — a
  reproducibility bug — where boundary pins `v3.10.5` / `3.4.0`.

So Phase D shipped a **provenance + version-pinning + drift-detection** layer
instead of a physical merge:

- `external_tools.lock.yaml` in each repo — a single diffable manifest of every
  external tool (`2HDMC 1.8`, `HiggsTools v1.2`, HB `v1.7`, HS `v1.1`), its
  upstream URL, vendored path, and `local_patches`. The dihiggs manifest records
  the 2HDMC λ1 fork explicitly so it is never mistaken for stock v1.8.
- `tests/test_external_tools.py` (identical in all three repos) — the repos are
  separate and can't see each other in CI, so cross-repo agreement is enforced
  by **shared pinned constants**: the tool versions and the byte-identical
  `conventions/physics_conventions.yaml` md5 (`a2fea4c8…`). Drift in any repo
  fails that repo's CI (the boundary/hep_cross workflows from Phase C run it).
- Fixed the floating-pin bug: `dihiggs/higgstools/external/CMakeLists.txt` now
  pins `json v3.10.5` / `eigen3 3.4.0`, matching boundary.

Still out of scope: a single *physical* shared conventions file would need a new
shared submodule repo wired into all three — the drift-guard above achieves the
same safety without it.

## How to verify

```bash
# boundary: unit tests + contract checks (no HiggsTools)
cd dihiggs_boundary && PYTHONPATH=python pytest -q -m "not integration"
PYTHONPATH=python python -m dhb.validate --contract evaluate_point_v1 \
  --input tests/fixtures/evaluate_point_sample.csv

# hep_cross: full suite + contract checks
cd dihiggs_hep_cross && PYTHONPATH=src pytest -q

# dihiggs: the shared-constant loader + external-tool drift guard
cd dihiggs && pytest tests/test_physics_conventions.py tests/test_external_tools.py -q

# external-tool versions + shared conventions md5 are pinned across all repos
md5sum */conventions/physics_conventions.yaml   # one shared hash

# no stray HBAR_C literal definitions remain in the lake pipeline
grep -rn "= 1.973269804e-13" dihiggs/mlpython/lake_pipeline/*.py   # only the pinned fallback
```
