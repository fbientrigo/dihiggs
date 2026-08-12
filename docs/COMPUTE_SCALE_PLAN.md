# Compute-scale plan — high-mass H2 point factory

Status: preparation only. No production campaign was launched and the
machine was not saturated under this task. This document specifies the
resume-safe runner contract that `docs/contracts/high_mass_campaign_template.yaml`
already references; it does not implement the runner.

## 1. Task granularity

One task = one `(m_H2, Delta_heavy, tan_beta, M2, lambda6, lambda7,
sin_beta_minus_alpha, yukawa_type, model_variant)` tuple, evaluated by a
single `DihiggsPointV2Evaluator` invocation with `n-mH=1, n-M2=1` (a
single-point "grid"). The evaluator already supports batched grids
(`n-mH>1`/`n-M2>1`) internally, but the campaign runner treats each physical
point as its own atomic task regardless, so that per-task status, logs, and
resume semantics are meaningful at the finest granularity the contract
requires (`task_key` in `high_mass_campaign_template.yaml`).

## 2. Atomic status

Each task's status lives in its own `status.json`, written via the same
write-temp-then-rename pattern `DihiggsPointV2Evaluator` already uses for its
CSV output (`c.output + ".tmp"` then `std::rename`). No shared/central status
database is touched mid-task; a crash mid-write leaves the previous valid
`status.json` (or none) in place, never a half-written file that a resume
pass could misread as `COMPLETE`.

States: `PENDING -> RUNNING -> {COMPLETE, FAILED}`. `RUNNING` additionally
records a heartbeat timestamp and the owning worker's PID; a runner sweep
that finds a `RUNNING` task whose heartbeat is older than a configurable
timeout (default: 10x the median observed per-task runtime, since
`DihiggsPointV2Evaluator` itself runs in well under a second per point —
see pilot timing below) marks it `FAILED` (stale) and makes it eligible for
re-queue.

## 3. Per-task artifacts

```text
<campaign_id>/<task_key>/status.json
<campaign_id>/<task_key>/run.log
<campaign_id>/<task_key>/manifest.json
<campaign_id>/<task_key>/point.csv        # raw dihiggs.point.v2 row
```

`task_key` is the sha256 of the canonical-precision CSV encoding of the
point's immutable physical coordinates plus `model_variant`, per
`high_mass_campaign_template.yaml`'s `task_key` section — not a grid index,
so re-running the same campaign config against a different mass-list
ordering, or merging two campaign configs that happen to share a point,
still resolves to the same task directory.

## 4. Seed provenance

`DihiggsPointV2Evaluator` itself is fully deterministic (no RNG; verified in
this task's pilot — `reproducibility.same_point_id_on_repeat` and
`identical_total_width`/`identical_ctau_mm` in
`docs/pilots/high_mass_h2_v1/pilot_validation.json`). No `seed` field is
needed for canonical-producer tasks; the manifest schema still reserves a
`seed` field (`null` for these tasks, never omitted) because later stages of
the same campaign machinery (MadGraph/Pythia event generation, out of scope
here) are stochastic and must record their seed for replay. Reserving the
field now avoids a schema-breaking addition later.

## 5. Resume without silent stale reuse

A task already `COMPLETE` is skipped on restart only if its stored
manifest's `evaluator_cli_argv` and `producer_commit` exactly match what the
current campaign config would produce for that `task_key`. Any mismatch
(code changed, or — should it ever happen — a `task_key` collision with
different construction arguments, which the hash construction is designed to
make practically impossible) causes the task to be treated as `FAILED`
(stale) and re-queued, never silently reused. This mirrors the existing
`dihiggs.point.v2.verification.v1` pilot's own provenance discipline
(`producer_commit`/`producer_dirty` recorded per row today).

## 6. Bounded workers

`DihiggsPointV2Evaluator` is single-process and single-threaded per
invocation (no internal `-j`/thread flag). "Workers" at the campaign level
means concurrent OS processes, each running one task to completion. A
`max_workers` config value bounds concurrency; this plan does not set a
specific number, since no campaign has been run under it. For sizing
guidance: each pilot point above completed in well under 0.1s of wall time
(single evaluator invocation, no I/O beyond a small CSV), so the bottleneck
for a large `mH2<=800/mA=mHp<=2000` campaign is expected to be task-management
overhead and disk I/O for the many small per-task files, not CPU time per
point — this favors a moderate worker count with efficient task-queue
polling over a very high worker count.

## 7. What this plan explicitly does not cover

- MadGraph/Pythia execution: those stages belong to `dihiggs_ufo` /
  downstream repositories (see `UFO_GENERICIZATION_REQUIREMENTS.md`), which
  are not modified by this task and are not yet resume-safe or per-point
  generic themselves (Gaps 2-3 in that report).
- Any specific `max_workers` value, cluster/scheduler integration, or actual
  campaign launch.

## 8. Sensitivity-analysis readiness (data shape, not analysis)

The per-task `point.csv` rows, once assembled into a campaign-level table,
are tidy (one row per point, one column per field) and directly usable for:

- mass effects (`m_H2_GeV`, `Delta_heavy_GeV` as explicit columns, not
  implicit in a filename),
- `tan_beta`, `M2_GeV2`, `lambda6` effects (already explicit input columns),
- `g_hH2H2_GeV`, `ctau_physical_mm`, all ten BRs, and `theory_ok`/
  `rejection_stage` as explicit columns for rejection-stage analysis.

Per the mission's explicit instruction, no SHAP or causal-attribution
modeling is performed or implied by this data shape. For future work:

- **Exact parameter interventions / controlled scans** (varying one
  coordinate at fixed others, as the pilot's `MASS_CONTROL` class already
  does) are primary causal evidence.
- **Sobol/fANOVA** global sensitivity analysis is appropriate once a
  `PHYSICAL_POINT_SCAN` dataset with a defined sampling distribution exists.
- **SHAP** may only be used as surrogate-model attribution on top of a
  validated predictive model of some campaign output (e.g. predicting
  `theory_ok` or `ctau_physical_mm` from input coordinates); any such model
  must report its own predictive validation (e.g. held-out accuracy/R^2)
  before any SHAP attribution from it is interpreted, and SHAP output must
  never be presented as a causal claim about the underlying 2HDM physics.
