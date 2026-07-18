# PR #60 Merge-Gate Review

Reviewer: Sonnet 5. Independent gate check: advisor (acting as Opus per plan). Head: `daf1ab885bd3050ba5d051dd92cf6cbba463a4bc`.

## Checklist disposition

1. **lambda1_v2 orchestration correctness** — PASS. `stable_point_id` is a deterministic FNV-1a64 hash over `%.17e`-formatted values (platform-stable, no PRNG/wallclock dependence). `format_float64` uses `.17g` (IEEE754 float64 round-trip safe, confirmed by test asserting exact round-trip of `1/10`). Cartesian ordering is a fixed nested loop (tan_beta → mH → lambda1), deterministic and tested exactly in `test_exact_header_cartesian_order_and_stable_ids`.
2. **Overwrite/--force/resume, partial-output/failed-run manifests, row-count/schema validation, provenance/dirty-state** — PASS with one P1 gap: `run_lambda1_v2`'s `subprocess.run(..., timeout=timeout_s)` is not wrapped, so a `TimeoutExpired` propagates with no manifest and no stdout/stderr log written — contradicts the "failed-run manifests" contract but only triggers when `--timeout` is set and exceeded. Not gate-blocking; filed as issue.
3. **Malformed-input row preservation** — PASS. `result.input = partial_input(fields, line)` is set before the try/catch in `Lambda1EvaluatorV2.cpp`, guaranteeing row identity survives parse/evaluate failures; PR #60's diff to this file is a readability/comment change only (the guarantee predates PR #60, confirmed already present).
4. **CLI defaults/aliases, M² rejects fixed lambda1, canonical/experimental separation** — PASS. `cli.py`'s `69375e88` commit renames engine choices to `lambda1_v2`/`lambda1_legacy`/`m2`/`m2_tracker`, defaults to `lambda1_v2`, adds `--mH-min/--mH-max/--n-mH` aliases, and both `m2` and `m2_tracker` explicitly reject `--lambda1` with a clear stderr message. Directly asserted in `test_m2_rejects_fixed_lambda1` and `test_parser_defaults_to_canonical_lambda1_v2_and_exposes_tracker`.
5. **Skip audit** (all skips added or altered by PR #60, per diff):

| Location | Reason text | Classification |
|---|---|---|
| `autoresearch/tests/test_alerting.py`, `test_bounded_adaptive_search.py`, `test_dihiggs_ucb_selection.py` | "autoresearch is frozen outside the canonical evaluation core" | VALID_FROZEN_SKIP — matches `docs/autoresearch_frozen.md`, consistent labeling |
| `mlpython/lake_pipeline/test_nulls.py` | `pytest.importorskip("polars", ...)` — was a hard `import polars` before | VALID_OPTIONAL_SKIP — converts a hard ImportError into an honest skip; pre-existing hardcoded local path (`/home/fabi/cern_db/...`) is a separate, pre-existing WSL-local-path concern, not introduced by this PR (flagged for Workstream 8) |
| `scripts/test_m2_band_tracker.py` | skip if experimental tracker binary not built | VALID_OPTIONAL_SKIP |
| `tests/test_group_width_errors_by_config.py`, `test_show_width_rel_stats.py`, `test_validate_colab.py` | skip if legacy comparison fixture missing | VALID_LEGACY_SKIP |
| `tests/test_recompute_readiness.py` | `importorskip("scripts.recompute_readiness", "frozen/quarantine component absent from this checkout")` | **VALID_LEGACY_SKIP** (see investigation below) |
| `tests/test_run_quarantine.py` | `importorskip("scripts.run_quarantine", "legacy component absent from this checkout")` | VALID_LEGACY_SKIP — module removed in commit `34ba02c7`, predates this PR chain entirely |

**`test_recompute_readiness.py` investigation** (1118 lines, 29 tests, the largest single skip by test count): `scripts.recompute_readiness` was never committed under that exact path in git history. The related module `scripts/audit_recomputed_campaigns.py` (added alongside the test in `ebe3a415`) does NOT export the imported symbols (`APPROVED_SCHEMA_WHITELIST`, etc.) — and was itself deleted in `34ba02c7`, the same commit that deleted `scripts/run_quarantine.py`. Both deletions predate PR #58/#60 entirely. An untracked, never-committed copy exists on disk at `main_dihiggs/scripts/archive/recompute_readiness.py` (a separate, stale local checkout) — evidence the component was built and used once but never merged into any branch. **Conclusion: this is genuine pre-existing legacy debt, not something PR #60 introduced or hides.** PR #60 only converts a hard top-level `ImportError` into an honest, labeled skip — a strict improvement over the base branch. Recorded for the legacy-cleanup issue (§9 item 8) with a pointer to the recoverable untracked copy, since silently restoring 262+1118 lines of unreviewed legacy code is out of this mission's bounded-fix scope.

6. **Documentation truthfulness** — PASS. `README.md`, `docs/OPERATING_STATUS_V2.md`, and `docs/contracts/canonical_evaluators_v2.md` are consistent with each other and with the mission's scientific contracts: M²/m12_sq distinction stated explicitly, `theory_ok_v1` scoped as theory-only, legacy/experimental/canonical clearly separated, and explicit non-claims section ("makes no full-campaign claim... no maintained paper-result claim").
7. **Test-coverage quality** — PASS, not smoke-only. `tests/test_orchestrator/test_lambda1_v2.py` (92 lines) asserts exact cartesian order, exact stable-point-ID value, exact float64 round-trip string, schema/row-count validation errors, and a full orchestration run (via a fake shell-script evaluator) verifying manifest provenance and row counts end-to-end. `tests/test_orchestrator/test_cli_contract.py` (36 lines) asserts CLI defaults, the M²-rejects-lambda1 error path, and exact dry-run input-CSV contents.

## Pytest tally caveat

Ran `python -m pytest -q -rs` in a **minimal venv** (pandas/numpy/pytest only, no C++ binaries built, no HiggsTools, no PyYAML): **551 passed, 117 skipped**. This is NOT the PR's claimed 579/89 — that comparison belongs to Phase 3/4's clean-checkout build with all binaries compiled. Most of the extra skips here (`DihiggsPointV2Evaluator`/`Lambda1EvaluatorV2`/tracker "not built", `yaml` import errors) are artifacts of this minimal environment and would pass/run in a full CI build. The categorical skip audit (above) is the Phase-1-relevant finding and holds regardless of environment.

## Verdict: **MERGE**

No defect found blocks merging PR #60. All P0-candidate items (recompute_readiness skip) resolved to VALID_LEGACY_SKIP with evidence, not MISSING_MAINTAINED_COMPONENT. One P1 (timeout/manifest gap) routed to the issue ledger, not gate-blocking.
