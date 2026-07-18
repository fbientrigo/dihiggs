# PR #58 Finalization Review (Phase 3)

Head reviewed: `bd4f6d495ee689bbf60fab2b004841eb19b04c0d` (23 commits ahead of `origin/main` @ `4d19bb97`).

## Full-diff checklist

- **lambda1 v2 / point v2 / M² v2**: reviewed in Phase 1 and confirmed working end-to-end in a fresh clean-checkout build (below).
- **Yukawa installation order (`6bfad766`)**: confirmed. Prior code called `model.set_yukawas_type(...)` *before* `model.set_param_phys(...)`; the fix installs the Yukawa type only after construction succeeds, and records `yukawa_type_installed` (read back from the model) as an explicit output field so consumers can verify it rather than trust the input echo. Applied consistently across `DihiggsPointV2Evaluator.cpp`, `Lambda1EvaluatorV2.cpp`, and the M² path. Adds `width_unaccounted_GeV` as a width-closure diagnostic (total width minus sum of tracked partial widths). Adds explicit `supported_yukawa_type()` validation replacing a bare int parse.
- **Lifetime `ctau_mm`**: unchanged formula (`kHbarCGeVmm / total_width`), correctly gated behind `width_ok`.
- **Row preservation / schema contracts / orchestrator**: confirmed in Phase 1.
- **Legacy compatibility (`PhysScanWithFixings`)**: untouched by this diff; clearly labeled non-canonical in docs.
- **CI/build logic**: `ci_pipeline.yml` (`push: main` only, Docker-based C++ unit/oracle suite) vs `python_guards.yml` (`push+PR`, lightweight guards) is an **intentional, in-file-documented split**, not an oversight — no fix needed.
- **Deletions** (`05a6217a` "remove generated artifacts"): reviewed name-by-name. All deletions are Sphinx `docs/_build/` output, stale committed 2HDMC/`chris/` compiled binaries (`2hdmc/working/ParamScan*`, `chris/diag_roundtrip`, etc.), and two dead archive notebooks/scripts. No maintained source lost — this is a hygiene fix (removing previously-committed binaries).
- **No committed binaries**: `git ls-files` under `dihiggs/app/`, `2hdmc/working/`, `chris/` at HEAD returns zero ELF files. The diff only *removes* previously-tracked binaries, adds none.

## Clean-checkout reproduction

Fresh `git clone` into `_verify/dihiggs-pr58-verify/` (new directory, no reuse of dirty worktrees), checked out to `bd4f6d49`.

```
make -C 2hdmc            → clean build, all legacy CLIs + lib2HDMC.a
make -C dihiggs clean
make -C dihiggs           → clean build, all canonical + legacy + experimental binaries linked
python -m pytest -q -rs   → 579 passed, 89 skipped, 0 failed
```

**579/89 matches the PR's original claim exactly** (my Phase-1 tallies of 551/117 and 562/106 were minimal-venv artifacts — no yaml/full toolchain — not a discrepancy; this fresh-clone run is authoritative).

**Deterministic micro-pilots** (bounded, 5–7 fixed points each, `OMP_NUM_THREADS=1` where applicable):
- `run_lambda1_v2_yukawa_pilot.py`: run-vs-run byte-identical; physics values (all columns except `evaluator_commit`/`evaluator_dirty`) bit-exact match against the committed golden pilot (0 mismatches / 5 rows × ~70 columns).
- `run_dihiggs_point_v2_pilot.py`: run-vs-run byte-identical across all 5 anchor cases; physics values (excluding `producer_commit`/`evaluator_dirty`) bit-exact match against committed golden output (0 mismatches).
- `tests/test_m2_tracker_intervals.py` + `scripts/test_m2_band_tracker.py` (experimental `Phys_M2BandTracker`): 4/4 passed.
- README command blocks (Lambda1 v2 smoke test, M² v2 smoke test, orchestrator `lambda1_v2` example) executed verbatim: all exit 0, correct row counts.

No full scans were run.

## Verdict: READY_TO_MERGE

All Phase 3 acceptance criteria met. CI reconfirmed green on the current head `bd4f6d49` (`guards`, `clean-lambda1-build`, GitGuardian — all SUCCESS, run completed 2026-07-17T23:33Z, post-dating the PR #60 merge). No unresolved review threads (draft PR, no reviewers assigned). Branch is 0 commits behind `main` (`main` unchanged throughout this session at `4d19bb97`). PR body already rewritten to describe the combined work (Phase 2). Non-claims preserved. Zero changes under paper/UFO/recast paths (confirmed via Phase 0 inventory — this PR only touches the `dihiggs` repo).

**Remaining gate: HD-1 — Fabian's explicit authorization is required before merging into `main`.** This is a human-only decision per the closure plan (`main` has no branch protection, so this merge is the sole safeguard before it's live). Per plan discipline, execution stops here pending that go-ahead; Phase 4 (merge + post-merge verification) is queued and ready to run immediately once authorized.
