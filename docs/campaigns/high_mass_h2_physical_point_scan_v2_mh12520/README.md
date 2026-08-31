# high_mass_h2_physical_point_scan_v2_mh12520

Recalculation of the 208 accepted high-mass H2 points under the **canonical**
SM-like Higgs mass convention `m_h = 125.20 GeV`
(`conventions/physics_conventions.yaml`, `sm_like_higgs.m_h_GeV`).

The v1 campaign (`../high_mass_h2_physical_point_scan_v1/`, `m_h = 125.13 GeV`)
is **untouched and retained as historical**. Do not mix rows from the two
directories in a single plot or table.

## Files

| File | Contents |
|---|---|
| `high_mass_valid_points_mh12520.csv` | 208 recalculated points, `dihiggs.high_mass_point.v1` schema |
| `mh_convention_delta_report.csv` | per-point per-field old -> new delta (4230 rows) |
| `point_id_map_v1_to_v2.csv` | old -> new `point_id` (all 208 changed) |
| `manifest.json` | provenance, flip counts, threshold move |

## How it was produced

Each accepted point's construction coordinates were replayed through the
canonical `DihiggsPointV2Evaluator` binary at `--mh 125.20`, one degenerate
1x1 grid per point. No physics is reimplemented; 2HDMC does all of it.

**Validation**: replaying the same pipeline at `--mh 125.13` reproduces
`high_mass_valid_points.csv` **byte-identically** (0 delta rows, 0 point_id
changes). So every difference recorded here is attributable to the mass
convention and nothing else.

## Results

- **All 208 `point_id` values changed.** `m_h` is hashed into the canonical
  point identity, so a point recalculated at a new convention is a *different
  point*, not a relabelling. Use `point_id_map_v1_to_v2.csv` to cross-reference.
- **No theory-validity flips.** `positivity_ok`, `unitarity_ok`,
  `perturbativity_ok`, `stability_ok`, `theory_ok`: 0 changes. All 208 remain
  theory-valid.
- **No cascade-flag flips.** The `H2 -> hh` threshold moved `2*m_h`:
  250.26 -> 250.40 GeV, but no point crosses it (the highest region is 250 GeV,
  which stays below), so `H2_to_hh_open` is `False` throughout and `BR_hh = 0`.
- **Reconstructed quartics shift**, as expected, and only those that depend on
  `m_h`:

  | quantity | max abs relative shift |
  |---|---|
  | `lambda1` | 3.97e-03 |
  | `lambda2` | 1.82e-03 |
  | `lambda3` | 4.69e-05 |
  | `lambda4`, `lambda5` | **0** (exactly unchanged -- no `m_h` dependence) |
  | `g_hH2H2_GeV` | 3.92e-03 (median 2.21e-04) |

- **Widths, BRs and `ctau` are physically unchanged**: every one shifts by
  <= 1.4e-05 and typically ~1e-15, i.e. floating-point noise. The single
  larger figure (`width_unaccounted_GeV`, 1.4e-05) is a residual of a
  near-cancelling sum, not a physical change.

## Limitation

This re-evaluates the previously **accepted** points only. It cannot discover
points newly admitted at `m_h = 125.20`; that would require replaying the full
attempted grid (55,080 points), which was deliberately not run.

## HB/HS status

`classification` is carried over from v1 as `HBHS_BLOCKED` and is **stale**:
the boundary blocker was removed in the same migration (see
`dihiggs_boundary/src/evaluate_point.cpp`, canonical v2 input schema). A tiny
3-point enrichment was run as proof; a full re-classification of these 208
points is a separate, not-yet-run step.
