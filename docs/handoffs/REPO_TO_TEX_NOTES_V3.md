# Repository-to-TeX notes handoff v3

This is the authority input for the later maintained-notes cleanup. It does
not edit or authorize edits to `/home/fabi/atlas_dihiggs/paper`.

## Exact code state

- Branch: `agent/dihiggs-operational-v3-20260717`
- Initial commit: `9bc300299208dab1ce0b9f7e7510c3b92b8979f4`
- Final implementation commit: `5835006`; the close-out documentation commit
  is the final repository HEAD reported in the delivery record.
- PR context: draft PR from this branch targeting
  `fix/lifetime-integrity-core-v2`, descended from PR #58 with PR #59 merged.
- Dirty state: clean after the final commit; generated binaries and pilot
  outputs are not part of the commit.
- Schemas: `dihiggs.lambda1.v2` and `dihiggs.point.v2`.
- Executables: `Lambda1EvaluatorV2`, `DihiggsPointV2Evaluator`, and the
  experimental `Phys_M2BandTracker`.

## Claims active TeX notes must use

- Two canonical row-preserving producers exist.
- Lambda1 and M² are distinct coordinates.
- `Lambda1EvaluatorV2` uses `THDM::set_param_phys_lam1`.
- `DihiggsPointV2Evaluator` uses `THDM::set_param_phys`.
- `Phys_M2BandTracker` is experimental and bounded-pilot only.
- `PhysScanWithFixings` is legacy/compatibility.
- `TRIPLE_OK` or `triple_ok_legacy` remains theory-only.
- `theory_ok_v1` currently equals the three theory predicates.
- Experimental gates remain unevaluated.
- Rows are preserved, including malformed, construction-failing, and rejected
  rows.
- Lifetime output is `ctau_mm` in millimetres.
- M² and `m12_sq` are distinct, with
  `m12_sq = M² * sin(beta) * cos(beta)`.
- The M² producer uses explicit default `mh = 125.20 GeV` with provenance.
- Lambda1 v2 receives `mh` explicitly in its input CSV.

## Mandatory TeX non-claims

Active notes must not claim:

- full scans were executed;
- final M² boundaries were validated globally;
- HB/HS are production gates;
- STU is a production gate;
- experimental acceptance was evaluated;
- LHC recasting is integrated;
- LHS is the current production method;
- autoresearch is canonical;
- `PhysScanWithFixings` is the preferred LLP producer.

## Exact replacement guidance

Replace the following old wording wherever it occurs in active notes:

| Old wording | Required replacement |
|---|---|
| “current production path is PhysScanWithFixings” | “canonical LLP production uses `Lambda1EvaluatorV2`; `PhysScanWithFixings` is legacy compatibility” |
| “M² engine uses fixed lambda1” | “the canonical M² producer reconstructs lambda1 from rectangular `(mH, M²)` points” |
| “M² equals m12_sq” | “`M² = m12_sq/(sin(beta)cos(beta))`; they are distinct GeV² quantities” |
| “failed points disappear” | “canonical producers preserve one output row per attempted input/grid point” |
| “lifetime is ctau_m” | “lifetime is `ctau_mm`, in millimetres” |
| “experimental viability follows from TRIPLE_OK” | “`TRIPLE_OK`/`theory_ok_v1` is theory-only; experimental fields remain unevaluated” |

The notes cleanup must retain the distinction between bounded engineering
pilots and scientific boundary evidence. No full campaign, new physics result,
HB/HS or STU gate, or recast result is implied by this handoff.
