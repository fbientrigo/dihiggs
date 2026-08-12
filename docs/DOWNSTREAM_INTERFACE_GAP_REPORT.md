# Downstream interface gap report — `model_point_to_llp_recast_v1`

Status: assessment only. `dihiggs_hep_cross` was inspected read-only for this
report (`dihiggs_hep_cross/contracts/model_point_to_llp_recast_v1.yaml`); no
files in that repository were modified.

## Current contract (v1), verbatim structure

```yaml
name: model_point_to_llp_recast_v1
produced_by: a BSM scalar-pair scan (e.g. dihiggs / dihiggs_boundary)
consumed_by: src/llp_recast, scripts/03_paper_s_benchmark_proxy.py
required_columns:
- model
- point_id
- m_scalar_GeV
- total_width_GeV
- ctau_mm
- sigma_production_fb
- BR_bb
- BR_WW
- BR_ZZ
- BR_gg
- BR_tautau
- BR_hadronic_proxy
- beta_gamma_source
- recast_channel_hint
aliases:
  ctau_mm: [ctau_mm_H2, total_width_H2 -> total_width_GeV alias]
invariants:
- ctau_mm == hbar_c / total_width_GeV   (rel_tol 1e-4)
- sum(BR_bb..BR_tautau) <= 1            (abs_tol 1e-6)
```

## Gap analysis against the high-mass point factory

| Missing / ambiguous in v1 | Why the high-mass factory needs it | Contract source in this task |
|---|---|---|
| `m_A_GeV`, `m_Hp_GeV` | v1 only has a single `m_scalar_GeV`; there is no way to express the heavy-state scale at all, let alone the `mA=mHp` constraint or `Delta_heavy`. | `high_mass_point_schema.yaml` fields `m_A_GeV`, `m_Hp_GeV`, `Delta_heavy_GeV` |
| `g_hH2H2_GeV` | Not present. Downstream currently has no way to receive the canonical production coupling; a v1 consumer that needs cross section context has no path except rederiving it, which the canonical-evaluator contract explicitly forbids. | `canonical_evaluators_v2.md` "Downstream repositories must consume `g_hH2H2_GeV`" |
| `model_variant` | Not present. v1's single `ctau_mm` cannot distinguish a physical lifetime from a forced-decay detector-response lifetime; a v1 consumer has no way to know which one it is looking at. | `high_mass_point_schema.yaml` field `model_variant` |
| Physical vs. response ctau | v1 has exactly one lifetime field (`ctau_mm`, aliased to `ctau_mm_H2`). This conflates the two concepts this contract requires kept separate. A v1-shaped record cannot represent a Variant A point without silently overloading the meaning of `ctau_mm`. | `high_mass_point_schema.yaml` fields `ctau_physical_mm` / `ctau_response_mm`; `HIGH_MASS_H2_CONTRACT.md` section 4 |
| Decay ownership | v1 has no field recording which stage computed a BR or forced a decay (canonical calculator vs. downstream Pythia forcing). Combined with the ctau ambiguity above, a v1 record cannot express "this BR_bb is the physical branching ratio" vs. "this sample was forced to decay 100% to bb for a response study." | `HIGH_MASS_H2_CONTRACT.md` section 4 (Variant A semantics) |
| `BR_tt` | Not present at all. This is the exact channel Gate A closed in the canonical evaluator (Section 2 of `HIGH_MASS_H2_CONTRACT.md`); v1 predates that fix and has no column for it. Any high-mass point with `m_H2` near or above 345 GeV has physically significant `BR_tt` that v1 cannot carry. | Gate A, `high_mass_point_schema.yaml` field `BR_tt` |
| `BR_hh` | Not present. `H2 -> hh` opens at `m_H2 > 250.26 GeV`, well inside the high-mass factory's target range, and is explicitly required as a "physical, not forbidden" channel by the cascade contract. | `cascade_contract.yaml` flag `H2_to_hh_open`; `high_mass_point_schema.yaml` field `BR_hh` |
| `BR_cc`, `BR_gammagamma`, `BR_Zgamma` | Not present in v1's required-column list at all (only `BR_bb, BR_WW, BR_ZZ, BR_gg, BR_tautau, BR_hadronic_proxy` are listed). The `sum(BR_bb..BR_tautau) <= 1` invariant sums `BR_bb, BR_WW, BR_ZZ, BR_gg, BR_tautau` (correctly excluding the derived `BR_hadronic_proxy` from the sum), but has no room for `BR_cc`, `BR_gammagamma`, or `BR_Zgamma` at all — three of the ten channels this contract's canonical evaluator exports are simply absent from v1. | `high_mass_point_schema.yaml` — full ten-channel BR list |
| Production process identity | `produced_by` is a free-text description ("a BSM scalar-pair scan"), not a structured field identifying the process (`pp -> H2 H2` vs. a future `pp -> H2 H2`-plus-cascade or `pp -> A/Hp` process). The contract explicitly requires the initial signal process to remain `pp -> H2 H2` only; v1 has no field that could assert or check this. | `HIGH_MASS_H2_CONTRACT.md` sections 1, 4 |
| Cascade-state flags | Not present. v1 has no way to record that a point's `H2_to_AZ_open` / `H2_to_HpW_open` / `H2_to_AA_open` / `H2_to_HpHm_open` flags were checked and are false, which is exactly the auditable guarantee a `PHYSICAL_DECAYS_NO_HEAVY_CASCADES` consumer would need before trusting a v1-shaped record's BR table as complete. | `cascade_contract.yaml` |
| `sigma_production_fb` semantics | v1 requires this as a plain required column with no provenance field distinguishing "actually produced by MadGraph for this exact point" from a place-holder/estimated value. The canonical calculator contract explicitly forbids the calculator itself from emitting a cross section; v1's schema does not record which stage produced this number or whether it is trustworthy. | `HIGH_MASS_H2_CONTRACT.md` section 8 |
| `beta_gamma_source` / `recast_channel_hint` | Present in v1, orthogonal to this contract's scope (ATLAS-recast-specific derived quantities); no gap identified here — flagged only as fields a v2 contract should carry forward unchanged unless `dihiggs_llp_recast` requests otherwise. | — |

## What a future `model_point_to_llp_recast_v2` should add, at minimum

Per the mission's own list (reproduced here for traceability, all confirmed
as genuine gaps above): `mA`, `mHp`, `g_hH2H2`, `model_variant`,
physical-versus-response `ctau`, decay ownership, `BR_tt`, `BR_hh`,
production process identity, cascade-state flags.

## Explicit non-goal of this report

This report does not modify `dihiggs_hep_cross/contracts/model_point_to_llp_recast_v1.yaml`
and does not draft a `v2` file in that repository. It exists so a future,
explicitly-scoped `v2` contract task in `dihiggs_hep_cross` has an accurate,
evidence-based starting gap list.
