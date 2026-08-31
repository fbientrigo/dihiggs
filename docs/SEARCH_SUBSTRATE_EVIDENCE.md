# Deterministic search-substrate evidence report

Status: reconnaissance completed before implementation. This is a compact
decision record for the deterministic benchmark-search substrate; it is not a
general repository audit and it does not authorize a search campaign.

## Repository state

- Checkout: `/home/fabian/atlas_dihiggs/dihiggs`
- Git branch: `feat/h2-highmass-continuation-mc`
- Git HEAD at reconnaissance: `ea2681069d2edbe9b059e7b2e17f3169d5a4653d`
- Pre-existing worktree change: untracked
  `scripts/run_diphoton_postfix_validation_v1.py` (preserved)
- Built canonical executable: `dihiggs/app/DihiggsPointV2Evaluator`

## Authoritative evaluator and physics layers

The maintained production path is `DihiggsPointV2Evaluator`, implementing
`dihiggs.point.v2` through `THDM::set_param_phys`, generic-parameter
round-trip reconstruction, the canonical `h-H2-H2` coupling, and the shared
2HDMC `DecayTable`. `dihiggs/include/YukawaType.hpp` installs and verifies the
requested Yukawa type only after successful construction. The evaluator
records construction, numerical reconstruction, theory predicates,
experimental unevaluated state, widths, branching ratios, and lifetime as
separate fields. `theory_ok_v1` is the positivity/unitarity/perturbativity
theory-only conjunction, not detector acceptance.

The active project convention is `m_h=125.20 GeV` from
`conventions/physics_conventions.yaml` and `docs/HIGH_MASS_H2_CONTRACT.md`.
The canonical contract explicitly distinguishes `M2` from `m12_sq`, and the
canonical lifetime constant is `1.973269804e-13 GeV mm`.

The historical L06 artifact is explicitly not used as a golden physics
fixture. The corrected verification documents show why: the pre-Yukawa-fix
L06 lifetime and branching ratios are qualitatively wrong. New substrate
tests use synthetic evaluator rows or current v2 execution, never that old
decay artifact.

## Existing family evidence

1. `docs/pilots/h2_qs_continuation_v1/clean_replay_v1.md` and its JSON report
   support the deterministic continuation invariants
   `Q=(mH2^2-M2)tan_beta^2` and `S=mA^2-mH2^2`, with
   `mA=mHp`, Type-I, alignment, and fixed `tan_beta/lambda6/lambda7`.
   The active-125.20 replay is theory-valid for 150, 175, and 200 GeV.
2. `runs/boundary_benchmark_discovery_20260821_v1/` contains a post-fix
   canonical-evaluator campaign whose labels identify an approximately
   balanced family near `x=lambda6*tan_beta≈1500`, a transition near 3000,
   and a diphoton-dominated family near 10000. The campaign also contains
   150/200/250-range mass studies and fixed-x lifetime decompositions.
3. `runs/diphoton_effective_coordinate_v1/summary/scientific_handoff.md`
   supports `x=lambda6*tan_beta` as an approximate composition coordinate
   under its fixed prescription and recommends x=1500/3000/10000 families.

## Missing authoritative decisions

The artifacts contain descriptive labels and example rows, but no explicit,
machine-readable, collaborator-approved numerical predicates defining
`mixed` and `photonic` (nor trusted active-125.20 fixtures for such predicates).
Therefore the classifier in the new substrate is an interface with an
explicit `BLOCKED_SCIENTIFIC_DECISION` state. It does not silently turn the
descriptive x=1500 and x=10000 labels into thresholds.

The requested preferred anchor interval `[500,1000] mm` is frozen as a search
requirement. The existing family evidence reports nearby lifetimes mostly
below that interval and does not justify a broader diagnostic envelope, so no
diagnostic envelope is frozen yet.

## Consequence for implementation

The substrate can be made reproducible and resumable now: proposals are
restricted to the canonical physical inputs, family validation derives the
150/200/250 nodes by deterministic Q,S continuation, and all physics decisions
come from the v2 evaluator. Archive promotion remains fail-closed while the
scientific classifier decision is blocked.
