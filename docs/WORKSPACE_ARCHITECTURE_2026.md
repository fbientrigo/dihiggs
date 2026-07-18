# Workspace Architecture — 2026

Describes how the repositories under the `atlas_dihiggs` workspace interact. The workspace root
itself is not a git repository (no durable version there); this document is the canonical copy,
with a convenience copy kept at the workspace root for local navigation only.

Repository integration described here does **not** imply completion of the physics program — it
describes what is built, tested, and interoperable as of the closure SHA below, not what physics
results have been produced or validated at scale.

## Repository roles

```text
dihiggs (fbientrigo/dihiggs)
    canonical 2HDMC point evaluation: Lambda1EvaluatorV2 (dihiggs.lambda1.v2),
    DihiggsPointV2Evaluator (dihiggs.point.v2), experimental Phys_M2BandTracker,
    legacy PhysScanWithFixings. Closed at 5b6ab63a3917967c638502d152ca63e1b51e4f11.

dihiggs_boundary (fbientrigo/dihiggs_boundary)
    theory/experimental boundary mapping via LHS sampling + HiggsBounds/HiggsSignals
    enrichment, over canonical dihiggs outputs. Must not redefine M², m12_sq, lambda1
    input semantics, theory_ok_v1, ctau_mm units, or Yukawa installation order.

dihiggs_hep_cross (fbientrigo/dihiggs_hep_cross)
    HEPData cross-section/recast bridge; packages surviving 2HDM points against the
    formal schema in contracts/model_point_to_llp_recast_v1.yaml for llp_recast.

dihiggs_llp_recast (fbientrigo/dihiggs_llp_recast)
    reproduces a public ATLAS displaced-vertex+jets recast, then reinterprets it for
    long-lived 2HDM scalars. Governed by its own AGENTS.md upstream-validation gate.
    Active, separate workstream.

UFO / MadGraph workspace (ufos/)
    model export (UFO), cross-section campaigns, event generation. Active, separate
    workstream; not a git repository in this workspace.

paper/
    manuscript and maintained notes; consumes frozen result artifacts. Active, separate
    workstream; not a git repository in this workspace.

wiki/, main_dihiggs/ (stale checkout), lambda1-*/dihiggs-point-v2 worktrees,
autoresearch/mlpython history
    legacy or archived; not part of the maintained closure surface.
```

## Dependency graph

```text
dihiggs (canonical point CSV + manifests: ctau_mm, theory_ok_v1, point IDs, provenance)
    │
    ├──► dihiggs_boundary   (boundary detection / interval analysis / visualisation)
    │
    ├──► dihiggs_hep_cross  (packages points against model_point_to_llp_recast_v1.yaml)
    │        │  benchmark-selection manifest (m_scalar_GeV, total_width_GeV, ctau_mm,
    │        │  sigma_production_fb, BR columns, beta_gamma_source, recast_channel_hint)
    │        ▼
    │   dihiggs_llp_recast  (DV+jets recast; own upstream-validation gate)
    │
    ├──► UFO / MadGraph workspace  (model export, cross sections, event generation)
    │
    └──► paper/  (frozen result artifacts only; L06 references require review — see issues)
```

## Data contracts between stages

| Contract | Producer | Consumer(s) | Notes |
|---|---|---|---|
| 2HDMC canonical point CSV | `dihiggs` (`Lambda1EvaluatorV2`/`DihiggsPointV2Evaluator`) | `dihiggs_boundary`, `dihiggs_hep_cross` | schemas `dihiggs.lambda1.v2` / `dihiggs.point.v2`, frozen in `docs/contracts/canonical_evaluators_v2.md` |
| Benchmark-selection manifest | `dihiggs_hep_cross` | `dihiggs_llp_recast` | `contracts/model_point_to_llp_recast_v1.yaml`; invariant `ctau_mm = ħc/Γ` |
| UFO parameter/model definition | UFO workspace | MadGraph | not yet formally interfaced with canonical `dihiggs` point IDs — see interface matrix |
| MadGraph cross-section/event metadata | MadGraph workspace | `dihiggs_llp_recast` | active workstream, out of this closure's scope |
| Recast input manifest | `dihiggs_llp_recast` upstream gate | recast analysis | governed by `llp_recast/AGENTS.md`; `NOT_RUN`/`BLOCKED`/`FAILED`/`PROVISIONAL`/`VALIDATED` vocabulary |
| Paper-facing frozen result artifact | any upstream stage | `paper/` | must be explicitly frozen/versioned before citation; no maintained paper-result claim from `dihiggs` alone |

## Interface matrix

| Repository | Role | Upstream dependency | Schema consumed | Executable assumed | Units assumed | Current compatibility | Blocking defect | Action now | Deferred issue | Owner |
|---|---|---|---|---|---|---|---|---|---|---|
| `dihiggs_boundary` | boundary analysis | independent 2HDMC construction (`src/evaluate_point.cpp`), not yet consuming `dihiggs` CSV output | own C++ schema; conventions (`M² = m12_sq / (sin β cos β)`, lambda1 reconstructed, Yukawa installed post-construction) verified equivalent to canonical | `build/bin/evaluate_point` | GeV, GeV², mm | **ALIGNED** — Yukawa order already correct (construction precedes `set_yukawas_type`), M²/m12_sq relationship matches canonical exactly, unit pytest 84 passed/1 skipped | none found | none | downstream-interface issue: wire `evaluate_point` to consume `dihiggs.point.v2` directly instead of independent construction (non-blocking, forward-looking) | Fabian |
| `dihiggs_hep_cross` | HEPData/recast bridge | none yet — no code consuming `dihiggs.point.v2`/`dihiggs.lambda1.v2` found | `contracts/model_point_to_llp_recast_v1.yaml` (not yet mapped from canonical output; requires case-renaming + `BR_hadronic_proxy` derivation + external `sigma_production_fb` from MadGraph) | none (CSV-only) | GeV, GeV, mm, fb | not yet exercised against v2 canonical output; no L06 references found | none found by static check | none | downstream-interface issue: publish canonical-CSV → contract mapping + one example manifest | Fabian |
| `dihiggs_llp_recast` | DV+jets recast | `dihiggs_hep_cross` manifest (not yet wired) | benchmark manifest + validation-state vocabulary | none (own upstream gate) | ctau in mm assumed | governed separately; no dihiggs v2 consumption and no L06 references found | none found | none — hands off | none required yet | Fabian |
| UFO/MadGraph workspace | model export, cross sections | none formal yet | none formal | none | none formal | pre-formal-interface stage; no L06 references found | none found (no interface exists to break) | none | downstream-interface issue: define canonical-point-ID handoff when this workstream needs it | Fabian |
| `paper/` | manuscript | frozen artifacts only | none | none | none | `notas/CLAIM_LEDGER.md` and `notas/VALIDATION_2026-07-17.md` **already document L06 as SUPERSEDED / retired as an LLP benchmark**, independently of this closure | none — already current | none (out of scope to edit; already correct) | none needed | — |

Silent-corruption-class incompatibilities are fixed on the `dihiggs`/`dihiggs_boundary` side or
explicitly block that interface; non-blocking, unexercised gaps are filed as issues per the table
above.
