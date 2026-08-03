# First H2 LLP recast candidate — corrected review status

## Verdict

`H2scan_mH150_tb300000` remains the **leading provisional candidate** from the
reported bounded scan, but it is **not yet a frozen benchmark and is not
cleared for downstream generation or recast**.

Two P1 review findings changed the admissible claim:

1. The previous numerical check perturbed `m12_2` at approximately
   `1e-12 GeV^2`. Near `m12_2 ≈ 0.075 GeV^2`, that is about five orders of
   magnitude larger than one double-precision ULP and changes reconstructed
   `lambda1` by order one. Those probes compare different physical models,
   not alternative numerical representations of the fixed
   `lambda1_target = 1` benchmark.
2. The scan runner excluded the entire `benchmarks/` directory from its git
   dirtiness check. Consequently, the v1 manifest could report
   `producer_dirty = no` even when the runner, grid, and constants were not in
   the recorded commit. The current scan cannot be reproduced from its
   manifest alone.

The previous channel-scoped labels — `VALID_FOR_FIRST_BB_RECAST` and
`NUMERICALLY_UNRESOLVED` for the loop channels — are therefore withdrawn until
the corrected checks are rerun.

## Reported candidate

| Field | Reported value |
|---|---:|
| `point_id` | `H2scan_mH150_tb300000` |
| `m_H2` | `150 GeV` |
| `mA = mHp` | `450 GeV` |
| `sin(beta-alpha)` | `1.0` |
| `tan_beta` | `3e5` |
| `lambda1_target` | `1.0` |
| `lambda6`, `lambda7` | `1e-10`, `0.0` |
| `construction_ok`, `theory_ok`, `width_ok` | `1`, `1`, `1` |
| `total_width_GeV` | `4.56118529862185e-14` |
| `ctau_mm` | `4.326221529733112` |
| `br_bb` | `0.7567374858085787` |
| `lambda1_abs_residual` | `9.34182041500975e-7` |
| `lambda1_residual_warning` | `1` |
| `sigma_production_fb` | pending |

These values are retained as the reported output of the existing CSV. They are
not promoted to a reproducible campaign result until the clean rerun below.

## Corrected numerical check

`benchmarks/check_H2scan_mH150_tb300000.py` now:

- reproduces the exact double-precision `m12_2` center used by
  `set_param_phys_lam1`;
- evaluates only the immediately adjacent representable doubles with
  `math.nextafter`;
- reports the local ULP, the half-ULP rounding bound, and the propagated
  `lambda1` error bound;
- classifies output stability only across the previous/center/next-double
  bracket.

A wider `m12_2` interval may still be useful as a **physical sensitivity
study**, but it must not be labeled numerical uncertainty of this fixed point.
The stale numerical-check CSV and Markdown outputs were removed; the corrected
runner must regenerate them.

Current channel status:

| Scope | Status |
|---|---|
| `H2 -> bb` DV+jets inputs | `PENDING_CORRECTED_ULP_RERUN` |
| `H2 -> gamma gamma` | `PENDING_CORRECTED_ULP_RERUN` |
| `H2 -> Z gamma` | `PENDING_CORRECTED_ULP_RERUN` |

## Corrected scan provenance

`benchmarks/run_first_h2_bounded_scan.py` now excludes only files it generates:

- `benchmarks/first_h2_bounded_scan.csv`;
- `benchmarks/first_h2_bounded_scan_manifest.json`;
- the two scratch CSV files.

The runner refuses to execute if any experiment-defining source or
configuration is dirty. Manifest v2 also records SHA-256 hashes for:

- the runner itself;
- `dihiggs/src/Lambda1EvaluatorV2.cpp`;
- `dihiggs/src/ReplaySafeOutput.cpp`.

The existing v1 manifest is explicitly marked
`INVALIDATED_PROVENANCE_REQUIRES_RERUN`. The evaluator commit recorded in the
CSV is still useful for identifying the evaluator source, but it does not
recover the missing runner state that determined the grid and selection.

## Required rerun sequence

From a clean commit containing the corrected scripts:

```bash
mv benchmarks/first_h2_bounded_scan.csv \
   benchmarks/first_h2_bounded_scan.pre_provenance_fix.csv
mv benchmarks/first_h2_bounded_scan_manifest.json \
   benchmarks/first_h2_bounded_scan_manifest.pre_provenance_fix.json

python3 benchmarks/run_first_h2_bounded_scan.py
python3 benchmarks/verify_pilot_ctau_invariant.py
python3 benchmarks/check_H2scan_mH150_tb300000.py
python3 -m json.tool benchmarks/first_h2_bounded_scan_manifest.json >/dev/null
git diff --check
```

After the rerun, update `FIRST_H2_RECAST_CANDIDATE.json` and this report from
the regenerated artifacts. Only then may a channel receive a validity label or
the point proceed to MadGraph/UFO normalization and the recast chain.

## Allowed current claim

The bounded scan reported six theory-valid, finite-width points in the
operational lifetime window and selected `H2scan_mH150_tb300000` as the leading
candidate. Because the scan provenance and numerical representation test are
not yet rerun under the corrected contracts, the point remains provisional and
unfrozen.
