# High-mass H2 2HDM point factory — scientific contract

Status: Gate A and Gate B closed. Bounded pilots pass (Section 6). No
production campaign, no MadGraph run, and no 800/2000 GeV full scan has been
launched under this contract.

Canonical producer: `DihiggsPointV2Evaluator` (`dihiggs/src/DihiggsPointV2Evaluator.cpp`),
schema `dihiggs.point.v2`. See `docs/contracts/canonical_evaluators_v2.md` for
the evaluator-level contract this document extends.

## 1. Physics target

```text
m_H2 <= 800 GeV
m_A  <= 2000 GeV
m_Hp <= 2000 GeV
hierarchy: m_Hp = m_A > m_H2 > m_h
```

Initial signal process: `pp -> H2 H2`. Heavy A/Hp production and cascade
feed-down (`H2 -> A Z`, `H2 -> Hp W`, `H2 -> A A`, `H2 -> Hp Hm`) are outside
scope until a later contract revision explicitly enables them (see
`docs/contracts/cascade_contract.yaml`).

Channels `bb, cc, tt, tautau, WW, ZZ, gg, gammagamma, Zgamma, hh` are
represented whenever kinematically open, per the evaluator's `DecayTable`
output.

## 2. Gate A — top-threshold audit (CLOSED)

**Finding.** Through commit `9f80196` (branch point of this work),
`DihiggsPointV2Evaluator` exported `width_bb_GeV` and `width_cc_GeV`
explicitly but not `width_tt_GeV`. Any `H2 -> t tbar` contribution above
threshold was silently absorbed into `width_unaccounted_GeV`. This blocked
high-mass physical-decay promotion per the gate rule.

**Fix.** Added `width_tt_GeV` and `br_tt` to `DihiggsPointV2Evaluator`,
computed via `DecayTable::get_gamma_huu(2, 3, 3)` on the same `DecayTable`
object used for the other nine selected widths (`bb, cc, tautau, WW, ZZ,
gammagamma, Zgamma, gg, hh`). The field is positioned immediately after
`width_cc_GeV` / `br_cc` in the CSV. `width_unaccounted_GeV` is computed
generically as `total_width - sum(all ten selected widths)`, so `tt` is now
included in the closure invariant automatically; no separate accounting path
was introduced.

**Verified invariants** (see `tests/test_dihiggs_point_v2.py`):

- `width_tt_GeV == 0` strictly below 2HDMC's internal top-pair threshold
  treatment (confirmed exactly zero at `m_H2 = 150, 200, 250 GeV` for the
  frozen anchor coordinates).
- `width_tt_GeV > 0` and dominant above threshold (confirmed `br_tt = 0.9668`
  at `m_H2 = 400 GeV`, anchor `mA=mHp=800, tan_beta=50, yukawa_type=1`).
- `width_unaccounted_GeV` stays at the pre-existing closure scale
  (`O(1e-9)` relative) with `tt` included in the explicit sum; it is not used
  to hide the `tt` channel.
- `sum(explicit partial widths) + width_unaccounted == total_width` to
  `1e-14` relative, unchanged in form from the pre-Gate-A invariant.
- `ctau_mm = hbar_c / total_width` is unaffected in form; its value changes
  wherever `tt` was previously (incorrectly) absent from `total_width`'s
  implicit accounting — it was not, since `total_width` always came from
  `DecayTable::get_gammatot_h(2)` directly, independent of which channels the
  evaluator chose to export. **Only the explicit/unaccounted split changes,
  not `total_width` or `ctau_mm` themselves.**

**Threshold behavior note.** 2HDMC does not impose a hard step at
`2 * m_t^pole = 345 GeV`. Empirically (Type-I anchor, `tan_beta=50`,
`mA=mHp=800`): `width_tt_GeV` is exactly `0` at `250 GeV`, `~1.5e-7 GeV` at
`300 GeV`, `~9e-5 GeV` at `350 GeV`, and dominant by `400 GeV`. This is
2HDMC's own running-mass/threshold treatment inside `get_gamma_huu`, not an
artifact introduced by this change. Any downstream consumer that assumed a
hard `2*m_t` cutoff must not.

**Schema impact.** This is a breaking column-order change to the
`dihiggs.point.v2` CSV: nine selected widths become ten. Consumers must key
columns by header name (`csv.DictReader` or equivalent), never by position.
No repository script or test in this codebase was found to depend on
position; all use name-keyed access (verified by repo-wide grep).

## 3. Gate B — frozen model assumptions

| Coordinate | Frozen value / convention | Source |
|---|---|---|
| `model_family` | General 2HDM, CP-conserving, Z2-softly-broken (via `lambda6`, `lambda7`, `m12_sq`) | `THDM::set_param_phys` physical-basis parameterization |
| Yukawa type | Type I (`yukawa_type=1`) is the default control point; Types I-IV (`1..4`, hep-ph/0504050 convention) are all supported and must be recorded per point, never assumed | `dihiggs/include/YukawaType.hpp`, `2hdmc/src/THDM.cpp:757-773` |
| `sin(beta-alpha)` | Scan/control coordinate, not fixed; benchmark anchors use `sin_ba in {0.999, 1.0}`. Record per point; do not silently default. | CLI `--sin-ba`, required argument |
| `m_h` | **Frozen at 125.13 GeV** for this campaign (canonical `dihiggs.point.v2` default, PDG-sourced: `https://pdg.lbl.gov/encoder_listings/s126.pdf`, per `docs/verification/dihiggs_point_v2_verification_v1.json`). The `125.09 GeV` value found in `docs/characterization_lambda1.md`, `scripts/generate_golden_lambda1.py`, and `tests/test_golden_lambda1.py` belongs to the older `lambda1`-characterization line of work and is **out of scope / not to be mixed** into this campaign. | See "Conflict resolution" below |
| `lambda7` | Scan/control coordinate, not fixed; pilot anchors use `0.0` or `1e-10` (soft-scale numerical stabilization, not a physical zero). Record per point. | CLI `--lambda7` |
| `mHp = mA` | Enforced by construction: both variants require `m_Hp = m_A`. The evaluator does not itself enforce equality (it accepts independent `--mA`/`--mHp`); campaign tooling (Section 5, `high_mass_campaign_template.yaml`) must always pass equal values, and `cascade_contract.yaml`'s fail-closed check on `H2_to_HpHm_open` / `H2_to_HpW_open` provides a second, independent guard. | Contract requirement, not a code invariant |
| 2HDMC | Vendored source, not a submodule. Version `1.8.0` (2020-02-10) per `2hdmc/README:159`, with a local patch aliasing `check_stability` to `check_positivity` (`kStabilityAlias` in the evaluator; see `docs/contracts/canonical_evaluators_v2.md`). Provenance is this repository's own commit (evaluator output already records `producer_commit`/`producer_dirty`); there is no separate upstream 2HDMC commit hash to track. | `2hdmc/README`, `DihiggsPointV2Evaluator.cpp` |
| SM inputs | `GF=1.16637e-5`, `MZ=91.15349`, `MW=80.36951`, `alpha_s=0.119`, `alpha_em=1/127.934`, `m_t^pole=172.5`, `m_b^pole=4.75`, `m_c^pole=1.42`, `m_tau=1.77684`, full CKM as listed in `2hdmc/src/SM.h:18-89`. These are 2HDMC 1.8.0 defaults, unmodified by this contract. | `2hdmc/src/SM.h` |

### Conflict resolution: `m_h = 125.09` vs `125.13`

Both values exist in this repository's history for different purposes:

- `125.09 GeV` — used by the `lambda1` scan-characterization line
  (`docs/characterization_lambda1.md`, `scripts/generate_golden_lambda1.py`,
  `scripts/lambda1_mh_eps_sensitivity.py`, `tests/test_golden_lambda1.py`).
  This line of work is **not** part of the high-mass H2 point factory and is
  left untouched.
- `125.13 GeV` — the canonical `dihiggs.point.v2` default, already frozen and
  verified in `docs/verification/dihiggs_point_v2_verification_v1.json`
  (`"mh_convention_GeV": 125.13`), and used throughout every point-v2 pilot
  anchor and regression test in this repository.

**Decision for this contract: `m_h = 125.13 GeV` for every high-mass point,
both `MASS_CONTROL` and `PHYSICAL_POINT_SCAN` classes, both study variants.**
Do not mix `125.09` results into any high-mass campaign plot, table, or
theory-validity pilot.

## 4. Two study variants

Both variants share the same canonical producer output; they differ only in
what downstream MadGraph/Pythia does with the point, and in which fields a
handoff consumer is permitted to treat as physically meaningful. See
`docs/contracts/high_mass_point_schema.yaml` for the `model_variant` field.

### Variant A — `FACTORIZED_G_ONLY`

- Purpose: clean production/kinematics + detector-response benchmarking.
- `pp -> H2 H2` is generated directly from `g_hH2H2_GeV`; H2 is stable in the
  LHE; a specific H2 decay may be forced in Pythia for response studies.
- Carries **both** `ctau_physical_mm` (from `total_width_H2`, informational
  only in this variant) **and** `ctau_response_mm` (whatever lifetime the
  response study independently scans). These must never be identified with
  each other, and a schema validator must reject any point where downstream
  metadata conflates the two field names.
- Any quoted physical yield must explicitly apply
  `sigma(pp->H2H2) * BR_channel^2 * Aeff_channel`; the canonical producer
  itself never emits a MadGraph cross section (see Section 7).

### Variant B — `PHYSICAL_DECAYS_NO_HEAVY_CASCADES`

- Purpose: richer simplified-2HDM signal with physical H2 lifetime and
  principal decays.
- Requires `mA = mHp > mH2 > mh` (checked by `cascade_contract.yaml`).
- `ctau_physical_mm = hbar_c / total_width_H2` only; never independently
  tuned.
- Stores all ten principal BRs listed in Section 1.
- `pp -> H2 H2` only; no A/Hp production, no cascade feed-down. Fail closed
  (Section 5) if any forbidden heavy-state channel opens.

## 5. Cascade-state contract (summary)

Full machine-readable definition: `docs/contracts/cascade_contract.yaml`.

Forbidden group (must be `false` for a valid `PHYSICAL_DECAYS_NO_HEAVY_CASCADES`
point; computed and recorded for every point regardless of variant):

```text
H2_to_AZ_open, H2_to_HpW_open, H2_to_AA_open, H2_to_HpHm_open
```

Physical group (must **not** be treated as forbidden; expected to be `true`
above their respective thresholds):

```text
H2_to_hh_open, H2_to_tt_open
```

Each flag is a pure kinematic (mass-threshold) boolean, independent of
whether 2HDMC's `DecayTable` reports a nonzero width for that channel at the
evaluator's current feature set (e.g. `H2 -> A Z` / `H2 -> Hp W` are not
computed by `DihiggsPointV2Evaluator` today — the flag exists precisely so a
future evaluator extension, or a hand-computed cross-check, can fail closed
without relying on an absent width field silently reading as zero).

## 6. Scan coordinates and staged campaign

Coordinates: `(m_H2, Delta_heavy)` with `m_A = m_Hp = m_H2 + Delta_heavy`,
never a rectangular `(m_H2, m_A)` grid (which would admit invalid-hierarchy
points). Both the transformed `(m_H2, Delta_heavy)` and physical
`(m_H2, m_A, m_Hp)` coordinates are stored per point
(`docs/contracts/high_mass_point_schema.yaml`).

Staged plan (see `docs/contracts/high_mass_campaign_template.yaml` for the
machine-readable version):

```text
Stage 0: reproduce existing 150 GeV anchor (mH2=150, Delta_heavy=300)
Stage 1: 130-300 GeV controlled pilot, mH2 x Delta_heavy in {150,300,500} GeV
Stage 2: deliberately cross tt (~345 GeV) and hh (~2*mh~250 GeV) thresholds
Stage 3: 300-500 GeV
Stage 4: 500-800 GeV
Heavy-state scale: progressively extend Delta_heavy so mA=mHp -> 2000 GeV
```

No stage beyond the bounded pilot (Section 7) has been executed under this
contract. Exact node values for Stage 1+ are deferred until a small
theory-validity pilot (Section 7) establishes the parameterization is
numerically usable across the full mass range — which it now has, for the
seven required pilot points.

Two experiment classes must never be mixed in one plot without labeling:

- `MASS_CONTROL`: vary `(m_H2, Delta_heavy)` holding `tan_beta, M2` (or
  `lambda6`) fixed wherever theory-valid.
- `PHYSICAL_POINT_SCAN`: allow `tan_beta, M2, lambda6` to vary as needed to
  find theory-valid points at each `(m_H2, Delta_heavy)`.

Every attempted point, including theory-rejected ones, is retained with its
`rejection_stage` / `rejection_reason` (already a `dihiggs.point.v2`
invariant — see `test_ordering_boundary_emits_success_and_failure_with_nan_masks`
and `test_construction_failure_preserves_row_without_decay_evaluation`).

## 7. Pilot results

Seven bounded pilot points were run per the contract's minimum list. Full
data: `docs/pilots/high_mass_h2_v1/pilot_points.csv` and
`pilot_validation.json` in the same directory, with checksums for every
artifact this task produced in `docs/pilots/high_mass_h2_v1/checksums.sha256`.
Summary:

Anchors P1-P6b hold `sin(beta-alpha)=0.999, tan_beta=50, M2=16721.68154468371
GeV^2, lambda6=0.1, lambda7=0, yukawa_type=1` fixed (MASS_CONTROL class),
varying only `(m_H2, Delta_heavy)`; P0 reproduces the pre-existing frozen
150 GeV benchmark exactly (`sin_ba=1.0, tan_beta=300000`, its own historical
anchor coordinates).

| Pilot | m_H2 | Delta_heavy | Purpose | theory_ok | width_tt_GeV | br_hh |
|---|---:|---:|---|:-:|---:|---:|
| P0_anchor_150 | 150 | 300 | existing 150 GeV anchor reproduction | 1 | 0 | 0 |
| P1_below_hh | 200 | 300 | below `hh` threshold (`2*125.13=250.26`) | 0 | 0 | 0 |
| P2_above_hh | 300 | 500 | above `hh` threshold | 0 | 2.35e-07 | 0.286 |
| P3_below_tt | 330 | 470 | below `tt` threshold region | 0 | 3.22e-06 | 0.310 |
| P4_above_tt | 400 | 400 | above `tt` threshold | 0 | 2.56e-03 | 0.318 |
| P5_near_800 | 800 | 1200 | near-800 GeV high point, mA=mHp=2000 | 0 | 2.03e-02 | 0.331 |
| P6a/P6b_heavy_sep | 300 | {200,900} | multiple heavy-state separations | 0 | 2.35e-07 | 0.286 |

Full data: `docs/pilots/high_mass_h2_v1/pilot_points.csv` /
`pilot_validation.json`. `theory_ok=0` for P1-P6b reflects that these
MASS_CONTROL anchors were not tuned for theory validity (that is
`PHYSICAL_POINT_SCAN`'s job, out of scope for this bounded pilot); per
Section 6, rejected points are retained as data, not treated as failures of
the pilot itself. Notably, `br_hh` is essentially zero at exact alignment
(`sin(beta-alpha)=1.0`, as used by the historical P0 anchor) but a visible
~30% once `sin(beta-alpha)=0.999`, confirming the evaluator correctly
exposes the alignment-limit suppression of the `H2->hh` vertex rather than
hiding it — this was discovered while building the pilot and is recorded
here as a physics-validity note for anyone reusing exact-alignment anchors.

All eight pilot points pass the required checks: same input produces the
same `point_id` (hash of physical coordinates only, independent of
`campaign_id`/`run_id`), worker-count independence is structural (the
evaluator is single-process, single-threaded per invocation — there is no
`workers=N` mode to diverge from; reproducibility is instead verified via
byte-identical repeated CLI invocation, which is the applicable form of the
same check), width decomposition closes to `1e-9` relative or better,
`ctau_mm` is internally consistent with `total_width_GeV`, the mass
hierarchy holds by construction (`mH2 < mA = mHp`, `mh < mH2`), cascade
flags are computed and consistent with the forced hierarchy, theory status
is recorded (not assumed), and `g_hH2H2_GeV` is finite and non-negative for
every constructed point. `width_tt_GeV` is confirmed exactly zero below
threshold and monotonically rising through `P1 < P3 < P4 < P5`.

## 8. Downstream interfaces (not modified in this task)

- `dihiggs_hep_cross/contracts/model_point_to_llp_recast_v1.yaml`: gap report
  in `docs/DOWNSTREAM_INTERFACE_GAP_REPORT.md`.
- `dihiggs_ufo` Pack B builder: requirements list in
  `docs/UFO_GENERICIZATION_REQUIREMENTS.md`.

Neither downstream repository was modified by this task.

## 9. Compute-scale preparation

See `docs/COMPUTE_SCALE_PLAN.md` and
`docs/contracts/high_mass_campaign_template.yaml`. No production campaign was
launched; the machine was not saturated.

## 10. Explicitly out of scope for this task

- MadGraph production of any kind.
- The full `mH2 <= 800 GeV` / `mA=mHp <= 2000 GeV` scan.
- Modifications to `dihiggs_boundary`, `dihiggs_hep_cross`, `dihiggs_llp_recast`,
  or `dihiggs_ufo`.
- Any SHAP/causal-attribution modeling (deferred per contract; see Section 9
  of `COMPUTE_SCALE_PLAN.md`).
