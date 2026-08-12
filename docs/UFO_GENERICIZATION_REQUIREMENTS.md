# UFO genericization requirements — `dihiggs_ufo` (read-only assessment)

Status: assessment only. `dihiggs_ufo` was inspected read-only for this
report; no files in that repository were modified.

## Summary

`dihiggs_ufo`'s current Pack A/Pack AA/Pack B chain is a single-benchmark
pipeline demonstration, not a generic high-mass 2HDM factory, and it is
explicitly self-documented as such (`PHYSICAL_POINT_UFO_HANDOFF.md`: *"must
not be treated as a generic generator for arbitrary 2HDM scan points"*;
`PACK_AA_RESEARCHER_RESULTS_BRIEF.md`: results labeled
`PIPELINE_DEMONSTRATION_ONLY / NOT_MODEL_DERIVED`). A future high-mass UFO
task starts from a materially different position than "extend an existing
generic factory" — the first gap below (no 2HDM UFO model) blocks everything
downstream of it.

## Gap 1 — no 2HDM UFO model exists in this repo (blocking)

The only UFO model present
(`releases/pack_a/frozen/pi_ufo_baseline_v1_frozen_hotfix1.zip` →
`model/LLscalar_v3_UFO_runtime/`) is a FeynRules "LLscalar" toy model: SM
particles plus one extra neutral real scalar `h2` (PDG `9000006`) and one
unused second scalar `MP`/`9000005`. It has:

- no pseudoscalar `A`,
- no charged Higgs `H+`/`H-`,
- no `tan_beta`, `sin(beta-alpha)`, or `lambda1..lambda7` parameters,
- a mass card (`point_param_card.dat`) with only `MH=125, MP=120, Mh2=200` —
  no `MA`, no `MHp`.

Its own `KNOWN_LIMITATIONS.md` states: *"Not a complete 2HDM UFO; no
lambda1...lambda7, tan(beta), m12^2, or full mixing structure"* and *"2HDMC
is not executed for new points."*

**Requirement:** a genuine 2HDM UFO model (independent `mH2`, `mA`, `mHp`,
`tan_beta`, `sin(beta-alpha)` or an equivalent physical-basis parameterization,
matching the same convention `DihiggsPointV2Evaluator` uses) must be produced
or imported (e.g. via FeynRules' 2HDM implementation) before `m_A`/`m_Hp` up
to 2000 GeV can be represented in MadGraph at all. This is a prerequisite for
every other item below, not a parallel work item.

## Gap 2 — Pack A/AA production is single-point and single-particle

`pack_aa/python/pack_aa.py` does read `mass_GeV`, `ctau_mm`, and a `decays`
list from a YAML config — decay/lifetime is parameterized within the
one-scalar model — but:

- the LLP PDG code is hardcoded to `9000006` in the `.cmnd` writer and in
  validation (`llp.pdg must be 9000006`), not read from config;
- `PACK_A_SHA256` is a module-level constant pinning one frozen zip;
- there is no `mA`/`mHp` concept anywhere, since the underlying model lacks
  them (Gap 1).

`design/pack_aa/PACK_AA_CONFIG_SCHEMA.yaml` already sketches a more generic
per-point config shape, but is explicitly marked *"design artifact only: no
Pack AA implementation reads it"* — i.e. the generic schema is aspirational,
not real, today.

**Requirement:** once Gap 1 is closed, Pack A/AA needs a real per-point
config path (point_id-keyed, not a module constant), and the LLP PDG-code
assumption needs to generalize past a single fixed particle.

## Gap 3 — Pack B is fully frozen to one benchmark, with no resume semantics

`pack_b_generate.py` extracts exactly `build/point_000001/decay.slha` and
`build/point_000001/point.json` from one pinned zip
(`config.json["pack_a_sha256"]`); there is no `point_id`/output-path
templating, no loop over points, and it refuses to write into a non-empty
output directory (i.e. no resume-safe per-point campaign semantics — compare
this task's `COMPUTE_SCALE_PLAN.md`, which requires atomic per-task status
and resume without silent stale reuse). `pack_b_preflight.py` cross-checks a
hardcoded `point["point_id"] == "A_PI_NATIVE_200"`.
`build_model_derived.py` hardcodes `POINT_ID`, `MASS_GEV`, `CTAU_MM`,
`GH_PHI_PHI`, `BR_BB` as module-level constants — running a different point
requires editing source, not passing arguments.
`docs/PHYSICAL_POINT_UFO_HANDOFF.md` self-declares this scope limit.

**Requirement:** Pack B needs to become a per-point, argument-driven
(point_id in, artifacts out under a point_id-keyed path) stage with the same
atomic-status/resume contract as the campaign runner described in
`COMPUTE_SCALE_PLAN.md`, once Gaps 1-2 are closed.

## Gap 4 — no `model_variant`, no `BR_tt`/`BR_hh`, decay is forced

No occurrence of `BR_tt`, `BR_hh`, `model_variant`, `FACTORIZED_G_ONLY`, or
`PHYSICAL_DECAYS_NO_HEAVY_CASCADES` exists anywhere in `dihiggs_ufo`.
`build_model_derived.py`/`RUNTIME_CONTRACT.md` hardcode
`pythia_decay_ownership: "forced H2 -> bb for sample efficiency"`; validation
`.cmnd` files hardcode BR splits (e.g. a 0.2/0.8 γγ/bb demo split) rather
than deriving them from a per-point BR table. This is consistent with
Variant A (`FACTORIZED_G_ONLY`) semantics as a *concept*, but the concept
itself is not represented anywhere in the code — there is no field or flag
that records which variant a given Pack B run corresponds to.

**Requirement:** once a real per-point BR table exists (Gap 1/2), Pack B
needs an explicit `model_variant` input that either (a) forces one Pythia
decay channel from the point's stored BRs (Variant A) or (b) hands Pythia the
full BR table for physical decays (Variant B), and every output artifact
must record which mode produced it.

## Gap 5 — lifetime is per-run-hardcoded, not looked up per point

`ctau_mm` appears as a fixed value in `config.json` (`100.0`) and as a module
constant in `build_model_derived.py` (`CTAU_MM = 4.326...`), not looked up
dynamically from a point's `ctau_physical_mm` /
`ctau_response_mm`. This is the same physical-vs-response ambiguity this
contract's schema (`docs/contracts/high_mass_point_schema.yaml`) makes
explicit; `dihiggs_ufo` today has no mechanism to distinguish the two even in
principle, since it only ever handles one frozen point.

**Requirement:** Pack B's per-point config (Gap 3) must carry both
`ctau_physical_mm` and `ctau_response_mm` as distinct, independently
sourced fields once it is genericized, matching this contract's schema.

## What is already in reasonable shape

- Pack AA's YAML-driven decay/lifetime configuration for its one scalar is a
  reasonable starting shape for a future per-point config, once the model
  supports more than one heavy scalar.
- The project's own design docs (`PACK_AA_CONFIG_SCHEMA.yaml`,
  `PHYSICAL_POINT_UFO_HANDOFF.md`) already anticipate several of the gaps
  above and explicitly recommend against building a generic framework
  "before it is needed" — consistent with this task's instruction not to
  modify `dihiggs_ufo` now. This report exists so that future work has an
  accurate, evidence-based starting checklist rather than rediscovering the
  same gaps.

## Explicit non-goal of this report

This report does not propose an implementation or timeline for closing these
gaps, and no code in `dihiggs_ufo` was changed to produce it.
