# Phase 6 — `dihiggs_boundary` Alignment Report

Checked against `origin/main` (`b3a70be`) of `fbientrigo/dihiggs_boundary`, no open PRs.

## Yukawa installation order (the named silent-corruption risk)

**Not incompatible.** `src/evaluate_point.cpp`'s main evaluation path already installs the Yukawa
type *after* physical-model construction: `model.set_param_phys(...)` (line 494) is checked for
success and returns early on failure (lines 506-509) before `model.set_yukawa_type(1)` is ever
called (line 512). This is the same ordering the canonical `dihiggs` repo adopted in its
`6bfad766` fix. A second, unrelated `set_yukawas_type(1)` call at line 160 operates on a disposable
`sm_like` reference object used only for Standard-Model coupling-ratio comparison inside
`compute_hbhs_block`, not on the physical model under evaluation — not the same bug class.

**Verdict: boundary does not have the pre-fix bug. No patch required.**

## M² / m12_sq convention

`r.m12_sq_input = r.M2 * denom` where `denom` is `sin(beta)*cos(beta)`, matching the canonical
`m12_sq = M² * sin(beta) * cos(beta)` relationship exactly. `lambda1` is populated only after
construction (`r.lambda1 = l1`, reconstructed output), never taken as a fixed input. No semantic
redefinition of the canonical conventions found.

## Schema / test health

- `python -m pytest` (unit only, `FakeHP` path, no HiggsTools): **84 passed, 1 skipped** — clean.
- Local checkout (`test/freeze-evaluate-point-v1`, `b4bf2ac`) is fully contained in `origin/main`
  (`b3a70be`, 1 ahead) — no divergent local work to reconcile.
- 4 untracked files present (audit notes, a plotting script) — left untouched, not part of this
  review's scope.

## Decision

No blocking incompatibility found. **No compatibility PR needed against `dihiggs_boundary` for this
closure.** `dihiggs_boundary` is classified **ALIGNED** with the canonical `dihiggs` v2 contracts as
of this check. Non-blocking, forward-looking gap: `evaluate_point.cpp` still constructs 2HDMC
points independently rather than consuming `dihiggs.point.v2` CSV output directly — this is an
architecturally valid choice (boundary's own C++ evaluator predates and is independent of the v2
schema) but means future canonical-core changes are not automatically inherited. Filed as a
non-blocking `downstream-interface`/`boundary-validation` issue, not a merge blocker.
