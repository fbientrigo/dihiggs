# Clean replay of the Q,S continuation pilot

The original pilot in [`continuation_report.json`](continuation_report.json)
was run from dirty working trees in both `dihiggs` and `dihiggs_boundary` (that
provenance is recorded honestly there and was not rewritten). After Task 1 and
Task 2 were frozen into commits, the identical anchor → 175 → 200 GeV
sequence was re-run from clean checkouts to check reproducibility:

- `dihiggs` at `ce722f70ee5965d2f4642a192d38240f953e36de` (clean, 0 dirty files)
- `dihiggs_boundary` at `f498e4df43c2afdf0953d1c693ff2b26f9da1cb4` (clean of
  source; only pre-existing, unrelated build-output directories under `lib/`
  remained untracked)

Compared field-for-field between the two runs: `point_id`, every `theory`
flag (`construction_ok`, `numerical_ok`, `theory_ok_v1`,
`positivity_reported_ok`, `unitarity_ok`, `perturbativity_ok`,
`stability_reported_ok`, `width_ok`, `rejection_stage`, `rejection_reason`),
`lambda1..lambda5`, all `diagnostics` (`g_hH2H2_GeV`, `total_width_GeV`,
`ctau_mm`, `br_bb`, `br_hh`, `br_tt`), the `qs_residuals` block, and the
evaluator executable's sha256.

**Result: every compared field was byte-identical between the dirty-tree
pilot and the clean replay**, for all three points (150, 175, 200 GeV).

This does not change the pilot's verdict
(`DIRECT_QS_CONTINUATION_VALID_TO_200_GEV`, unchanged) — it confirms the
result does not depend on the uncommitted state either run happened to start
from.
