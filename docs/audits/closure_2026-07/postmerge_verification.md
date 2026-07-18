# Post-Merge Verification — PR #58 → main

**Merged**: 2026-07-18T01:26:19Z. Merge commit / final `main` SHA: `5b6ab63a3917967c638502d152ca63e1b51e4f11`.
Authorized by Fabian (HD-1) after the Phase 3 READY_TO_MERGE report.

Verification performed from a **fresh clean clone** at `/home/fabi/atlas_dihiggs/_verify/dihiggs-main-5b6ab63/` (new directory, never reused from a dirty worktree).

| Step | Result |
|---|---|
| `make -C 2hdmc` | clean, all legacy CLIs + `lib2HDMC.a` built |
| `make -C dihiggs clean && make -C dihiggs` | clean, all canonical/legacy/experimental binaries built |
| `python -m pytest -q -rs` | **579 passed, 89 skipped, 0 failed** |
| `run_lambda1_v2_yukawa_pilot.py` | bit-exact match vs. committed golden pilot (0 mismatches, provenance fields excluded) |
| `run_dihiggs_point_v2_pilot.py` (5 anchor cases) | bit-exact match vs. committed golden pilots (0 mismatches, provenance fields excluded) |
| README Lambda1 v2 smoke test | exit 0, 1 row written |
| README M² v2 smoke test | exit 0, 1 row written |
| `git status` after all steps | clean except expected untracked build binaries (`dihiggs/app/{Lambda1EvaluatorV2,DihiggsPointV2Evaluator,Phys_M2BandTracker}` — matches PR body's stated "compiled v2 binary remains untracked") and pilot-regen provenance-field diffs (expected, embeds current commit) |

## Findings for follow-up (non-blocking)

- `docs/OPERATING_STATUS_V2.md` line 3 cites `5835006` as the "authoritative implementation commit" — stale as of the actual merge SHA `5b6ab63`. Queued for the Phase 5 closure-docs PR.

## Verdict: main is operationally verified at `5b6ab63a3917967c638502d152ca63e1b51e4f11`.

No incident. No hotfix needed.
