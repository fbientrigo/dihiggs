# Physics authority discrepancy and supersession ledger

Audit date: 2026-09-07  
Authority: [`conventions/physics_conventions.yaml`](../conventions/physics_conventions.yaml)  
Scope: the five active repositories at the audited `main` commits below

This ledger records disagreement; it does not rewrite frozen artifacts. A row
marked `superseded_for_new_production` remains valid for exact historical replay
when its original convention and provenance are retained.

## Status vocabulary

| Status | Meaning |
|---|---|
| `authoritative` | Governs new production after merge |
| `aligned_non_authoritative` | Semantically agrees but does not own the quantity |
| `superseded_for_new_production` | Retained for replay; cannot govern a new run |
| `migration_required` | Active code/document must be updated before producing under v3 |
| `incorrect_claim` | Statement is factually wrong and must not be propagated |

## `fbientrigo/dihiggs`

Audited `main`: `2264ffe1b8d8f59bec69341d6aa1a917396ccb2b`

| File | Observation | Status under v3 | Required migration |
|---|---|---|---|
| `conventions/physics_conventions.yaml` | v1 omits `m_h`, uses ambiguous `H` for PDG 35, and describes three byte-identical copies as co-authorities | `authoritative` after this change | Replace with v3; `dihiggs` becomes the only owner |
| `external_tools.lock.yaml` | Pins the v1 MD5 and byte-identical-copy model | `migration_required` | Pin canonical v3 with SHA-256; downstream provenance moves to sidecars |
| `tests/test_external_tools.py` | Enforces the obsolete MD5/co-authority policy | `migration_required` | Test v3 authority and SHA-256 instead |
| `tests/test_physics_conventions.py` | Parses constants but does not test schema or physics semantics | `migration_required` | Add v3 schema and semantic invariants |
| `docs/inter_repo_communication_assessment.md` | Records the former three-copy/MD5 design as current | `superseded_for_new_production` | Add an authority notice; retain the assessment as design history |
| `README.md` | New-looking smoke/orchestrator examples use `125.13` | `migration_required` | Issue #84 must make active examples consume `125.20`; historical examples require labels |
| `docs/contracts/canonical_evaluators_v2.md` | Calls `125.13` the default/current convention | `migration_required` | Issue #84 must distinguish v2 replay from new production |
| `docs/HIGH_MASS_H2_CONTRACT.md` | Freezes `125.13` for its campaign and incorrectly calls `lambda6`, `lambda7`, and `m12_sq` soft breaking | mixed: historical mass; `incorrect_claim` for Z2 wording | Preserve campaign numbers; issue #84 corrects active wording and labels the campaign historical |
| `dihiggs/app/orchestrator/cli.py`, `manifest.py`, `engines/m2_tracker.py` | Carry active-looking `125.13` defaults | `migration_required` | Issue #84 removes or updates silent new-production defaults without changing frozen artifacts |
| `docs/campaigns/high_mass_h2_physical_point_scan_v2_mh12520/` | A 125.20 recalculation exists in repository history but is not present at audited `main` | non-authoritative evidence | Do not infer authority from an unmerged/newer artifact; promote only through normal review |

## `fbientrigo/dihiggs_ufo`

Audited `main`: `7dd02c14e2314f9d172a0782d7860826a52600d6`

| File | Observation | Status under v3 | Required migration |
|---|---|---|---|
| `docs/PHYSICAL_POINT_UFO_HANDOFF.md` | Correctly assigns stage ownership, but says `125.20` is the PDG 2026 value and claims byte-identical v2 copies in three repos | `incorrect_claim` plus `migration_required` | Cite the project decision, consume v3 provenance, and remove co-authority wording |
| `README.md` | Describes the repository as owning MadGraph integration and reports `Recast: NOT_RUN`, which is only local Pack status and can be mistaken for project-wide state | `migration_required` | Narrow the repository role to UFO/model implementation; link external stage authorities rather than restating them |
| `CURRENT_STATE.md` | Local Pack/recast status can be read as ecosystem status | `superseded_for_new_production` outside local Pack scope | Retain as local history/status; add scope and authority links |
| `pack_b/operator/build_model_derived.py` | Frozen builder uses historical `125.13` and labels it historical | `superseded_for_new_production` | Preserve bytes/behavior for replay; do not turn it into the generic point interface |
| repository root | No canonical convention cache/lock pair exists | `migration_required` | Either consume point-level v3 provenance or add a byte-identical cache plus lock sidecar |

## `fbientrigo/dihiggs_hep_cross`

Audited `main`: `be331768ab644dc71fa4e8c46703a9fc3f567696`

| File | Observation | Status under v3 | Required migration |
|---|---|---|---|
| `conventions/physics_conventions.yaml` | v2 contains `125.20` but presents itself as a source of truth, retains `H/35`, and falsely attributes `125.20` to the PDG 2026 listing | `incorrect_claim` and `migration_required` | Replace with byte-identical v3 cache and provenance sidecar |
| `docs/MH_CONVENTION.md` | Repeats the PDG attribution and three-copy authority model | `incorrect_claim` | Cite v3 project authority and distinguish external measurement from production choice |
| `docs/CURRENT_DATA_AUTHORITY.md` | Correctly owns current-result promotion, but falsely says `dihiggs/main` already establishes `125.20` and blurs data authority with physics authority | `migration_required` | Keep result authority local; defer physics semantics to `dihiggs` v3 |
| `state/CURRENT_CAMPAIGN.yaml` | Serializes `m_h_GeV` as a YAML number, cites the wrong authority chain, and registers only two of the three comparison prescriptions | `migration_required` | Use decimal text plus v3 commit/checksum and all required prescription labels |
| `scripts/run_physical_point_madgraph.py` | Correctly refuses a missing point mass; current coupling handoff uses the magnitude compatibility field | `aligned_non_authoritative` with pending sign gap | Preserve fail-closed mass behavior; signed physical vertex waits for issue #85 |

## `fbientrigo/dihiggs_llp_recast`

Audited `main`: `72b8638671b788c35260491c06b3ab406775b5ae`

| File | Observation | Status under v3 | Required migration |
|---|---|---|---|
| `README.md` | Correctly owns cutflow/acceptance/efficiency and rejects production-normalization inference; its 150 GeV benchmark uses historical upstream quantities | `aligned_non_authoritative` | Label incoming physics convention IDs and retain the benchmark as historical response evidence |
| repository root | No canonical convention cache/lock pair and no v3 point-provenance requirement were found | `migration_required` | Validate convention schema/commit/checksum on new sample manifests |
| `results/` and benchmark documentation | Frozen response artifacts must not become model-point or cross-section authority | `superseded_for_new_production` when reused outside their frozen sample | Preserve hashes; compose only with explicitly compatible upstream provenance |

## `fbientrigo/dihiggs_boundary`

Audited `main`: `32dd9df533cd4b4d613daaf516b6b25254bdaf25`

| File | Observation | Status under v3 | Required migration |
|---|---|---|---|
| `conventions/physics_conventions.yaml` | v1 omits `m_h`, uses `H/35`, and presents a byte-identical copy as authority | `migration_required` | Replace with byte-identical v3 cache and provenance sidecar |
| `docs/model_contract.md` | Presents `125.09`, exact alignment, `lambda7=0`, and `mHp=mA` as the active model contract | `superseded_for_new_production` | Preserve as legacy contract or revise to require explicit campaign values |
| `docs/evaluate_point_contract.md` | Hard-codes `125.09` and exact alignment | `superseded_for_new_production` | Restrict to legacy replay; new composition consumes v3-bearing inputs |
| `docs/hbhs_contract.md`, `python/dhb/runner.py` | Use a `125.09` HiggsSignals reference | `migration_required` for new production; historical for exact regression | Separate the statistical reference mass from model-point `m_h` and record both explicitly |
| `docs/atlas_schema.md` | Freezes `125.09`, exact alignment, `mHp=mA`, and `lambda7=0` as Atlas-v0 values | `superseded_for_new_production` | Keep v0 replayable; v3-aware schemas must carry explicit campaign choices |
| `README.md` | Current composition role is substantially aligned and correctly requires direct per-point MadGraph results | `aligned_non_authoritative` | Add v3 provenance checks during the downstream migration |

## Supersession rules

1. This v3 contract supersedes semantic claims, not artifact bytes.
2. `125.13`, `125.09`, and `125.0` remain reproducible only under their named
   historical convention IDs.
3. A file timestamp, branch name, local copy, or later commit cannot supersede
   the authority automatically.
4. `dihiggs_hep_cross` may be authoritative for current MadGraph results while
   remaining non-authoritative for the model convention; authority is per
   stage, not per repository age.
5. An aligned downstream statement still requires canonical commit and
   checksum provenance before it can govern a new handoff.
6. Any conflict not listed here is rejected and added to this ledger before a
   migration proceeds.

## Review boundary

No appendix material or appendix branch was inspected or modified for this
ledger. No evaluator, scan, generated dataset, benchmark, MadGraph run, recast
run, or limit calculation was changed.
