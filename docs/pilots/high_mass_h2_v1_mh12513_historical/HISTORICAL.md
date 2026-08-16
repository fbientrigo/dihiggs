# HISTORICAL — high-mass H2 bounded pilot at m_h = 125.13 GeV

Frozen snapshot of `docs/pilots/high_mass_h2_v1/` as produced under the
**superseded** mass convention `m_h = 125.13 GeV`.

- Convention: `m_h = 125.13 GeV` (PDG 2026 listing), superseded by the
  canonical `m_h = 125.20 GeV` in `conventions/physics_conventions.yaml`
  (`sm_like_higgs.m_h_GeV`).
- Every `mh_input_GeV` in `pilot_points.csv` reads `1.25129999999999995e+02`.
- `H2_to_hh_open` in this snapshot uses the old threshold `2*m_h = 250.26 GeV`
  (the canonical threshold is now `250.40 GeV`).
- `point_id` values here are **not** comparable with post-migration ids:
  `m_h` is hashed into the canonical point identity, so recalculating at
  125.20 re-ids every point.

Retained for exact regression against the pre-migration state. Do **not**
reinterpret these rows as newly calculated `m_h = 125.20` points, and do not
mix them into the same plot or table as post-migration results.

See `docs/HIGH_MASS_H2_CONTRACT.md`, "Conflict resolution", Decision v1
(superseded) and Decision v2 (active).
