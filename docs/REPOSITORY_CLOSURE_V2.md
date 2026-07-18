# Repository Closure v2

Status: closed. Final `main` SHA: `5b6ab63a3917967c638502d152ca63e1b51e4f11` (merge of PR #58, which
itself absorbed PR #60). Closure performed 2026-07-17/18. This is an engineering readiness record,
not a physics-result claim.

## Supported workflows

- Clean build: `make -C 2hdmc && make -C dihiggs clean && make -C dihiggs`.
- Canonical lambda1 evaluation: `dihiggs/app/Lambda1EvaluatorV2 <input.csv> <output.csv>`, schema
  `dihiggs.lambda1.v2`, or via `python -m dihiggs.app.orchestrator --engine lambda1_v2 ...`.
- Canonical M² evaluation: `dihiggs/app/DihiggsPointV2Evaluator ...`, schema `dihiggs.point.v2`, or
  via `python -m dihiggs.app.orchestrator --engine m2 ...`.
- Bounded experimental boundary-search pilots via `dihiggs/app/Phys_M2BandTracker`
  (`--engine m2_tracker`) — non-canonical, intervals are not final scientific evidence.
- Legacy replay/compatibility via `PhysScanWithFixings` (`--engine lambda1_legacy`) — not permitted
  for new LLP lifetime production; historical output loses lifetime/BR information.
- Root `python -m pytest` collects cleanly: 579 passed, 89 skipped, 0 failed (reproduced from a
  fresh clean clone at the closure SHA; see `docs/audits/closure_2026-07/postmerge_verification.md`).

## Unsupported / deferred workflows

- Full parameter-space scans or campaigns.
- HiggsBounds/HiggsSignals production integration (theory-only flags exist; no experimental gate is
  evaluated by the canonical core).
- STU production integration.
- MadGraph cross-section generation, event generation, or LLP recast execution (these are separate,
  active workstreams in `hep_cross/`, `llp_recast/`, and the UFO/MadGraph workspace — see
  `docs/WORKSPACE_ARCHITECTURE_2026.md`).
- Model training, SHAP analysis, or other `mlpython/`/`autoresearch/` functionality — frozen
  (`docs/autoresearch_frozen.md`).
- Paper-result synchronization.

## Known limitations

- **L06 is not a valid LLP benchmark.** The Yukawa installation-order fix (commit `6bfad766`)
  changes widths, branching ratios, and lifetimes computed under the prior (buggy) order. L06's
  original characterization as a long-lived benchmark point predates the fix and is invalidated.
  This closure does not propose a replacement point (see issue: physics-recompute / LLP benchmark
  search).
- All widths, branching ratios, and lifetimes produced before the Yukawa-order fix should be
  treated as unverified and are candidates for recomputation (see issue: physics-recompute).
- `Phys_M2BandTracker` is validated only inside bounded pilot domains; its intervals are not final
  boundary evidence.
- `theory_ok_v1` currently equals the three theory predicates (positivity, unitarity,
  perturbativity/stability); no experimental viability is implied by any theory-only flag.
- Several pre-existing test files carry honest skips for components removed or never committed
  prior to this closure (`tests/test_recompute_readiness.py`, `tests/test_run_quarantine.py`,
  legacy comparison-fixture tests) — see the PR #60 skip-classification table in
  `docs/audits/closure_2026-07/pr60_merge_gate_review.md`. None hide currently-maintained behavior.

## Build and test commands

```bash
make -C 2hdmc
make -C dihiggs clean
make -C dihiggs
python -m pytest -q
```

Prerequisites: Linux/WSL, C++17 compiler, GNU Make, GSL. HiggsTools/HiggsBounds/HiggsSignals,
MadGraph, recasting tools, ML packages, and plotting packages are optional and not required for the
canonical v2 evaluators or their tests.

## Current schemas

- `dihiggs.lambda1.v2` — explicit `lambda1_target` input, `THDM::set_param_phys_lam1`. See
  `docs/contracts/canonical_evaluators_v2.md`.
- `dihiggs.point.v2` — rectangular `(mH, M²)` input, `THDM::set_param_phys`, lambda1 reconstructed.
  Same contract document.

## Downstream interface

- `hep_cross/` packages canonical output against
  `contracts/model_point_to_llp_recast_v1.yaml` (invariant `ctau_mm = ħc/Γ`; explicitly forbids
  conflating the paper's `lambda_eff`/`sin_theta` with 2HDM `lambda6`/`lambda7`/`sin(β-α)`).
- `dihiggs_boundary` consumes canonical CSV output for boundary detection and visualisation; must
  not redefine `M²`, `m12_sq`, lambda1 input semantics, `theory_ok_v1`, `ctau_mm` units, or the
  Yukawa installation order. See `docs/WORKSPACE_ARCHITECTURE_2026.md` for the full dependency
  graph and the interface matrix for current alignment status.

## Issues (deferred work)

See the issue ledger opened in `fbientrigo/dihiggs` (labels: `physics-recompute`,
`experimental-gates`, `boundary-validation`, `downstream-interface`, `legacy-cleanup`,
`documentation`, `CI`) for the itemized backlog. None of the deferred items block the workflows
listed as supported above.

## Reopening criteria

Reopen repository-core work on `dihiggs` when: a defect is found in a canonical evaluator's
construction, theory-predicate, width/BR, or lifetime computation; a downstream consumer
(`dihiggs_boundary`, `hep_cross`, `llp_recast`) reports a schema or unit mismatch that could
silently corrupt results; or a new canonical producer/schema version is needed. Routine campaign
work, ML/autoresearch restoration, and paper-result production do not require reopening this core —
they are separate workstreams that consume the canonical output as-is.
