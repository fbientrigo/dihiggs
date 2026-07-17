# Operating status v2

Authoritative commit: to be filled with the final operational-v3 commit.
Initial branch commit: `9bc300299208dab1ce0b9f7e7510c3b92b8979f4`.
This status is an engineering readiness record, not a physics-result claim.

## Canonical components

- `Lambda1EvaluatorV2` → `dihiggs.lambda1.v2`, explicit `lambda1_target`.
- `DihiggsPointV2Evaluator` → `dihiggs.point.v2`, rectangular `(mH, M²)`.
- `dihiggs/app/orchestrator` exposes `lambda1_v2` and `m2` as canonical paths.

Both producers preserve attempted rows, record repository provenance, and emit
theory-only acceptance fields. Lifetime output is `ctau_mm`.

## Experimental components

`Phys_M2BandTracker` is a bounded-pilot boundary-search helper. Its intervals
are non-canonical artifacts and no experimental gate is evaluated or promoted.

## Legacy components

`PhysScanWithFixings` and the `lambda1_legacy`/`lambda1` orchestration path are
replay and compatibility components. Autoresearch and historical campaigns are
frozen. They are not new LLP lifetime production paths.

## Passing checks

The final handoff records the exact commands and results for the clean 2HDMC
build, all three executable builds, focused evaluator/orchestrator tests, root
pytest, deterministic bounded micro-pilots, row-preservation cases,
schema/row-count checks, diff whitespace, workflow YAML parsing, and clean
checkout reproduction where available.

## Optional dependencies

The maintained v2 2HDMC producers require a C++17 compiler, Make, GSL, and the
patched repository 2HDMC checkout. HiggsTools, HiggsBounds, HiggsSignals,
MadGraph, recasting tools, ML packages, and plotting packages are optional and
are not required for the canonical v2 smoke tests.

## Deferred components

Full campaigns, HB/HS datasets, STU production integration, recast integration,
model training, large plotting campaigns, and paper-result synchronization are
deferred. The maintained core makes no full-campaign claim, no STU production
gate claim, no HB/HS production-gate claim, no recast-integration claim, and no
maintained paper-result claim.
