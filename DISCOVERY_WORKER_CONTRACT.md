# Discovery worker contract (READY_FOR_DISCOVERY_SEARCH)

You are a search worker. You propose *regions and parameters*; the frozen
deterministic evaluator decides all physics. Nothing you believe about a point
is evidence. Only `evaluation.json` / CLI output is evidence.

## Absolute prohibitions

- Do NOT edit anything under `search_substrate/` (frozen evaluator, gates, ctau,
  normalization, provenance, ledger). Do NOT edit `search_discovery/`.
- Do NOT edit the contract, the archive, or `ledger.jsonl` by hand.
- Do NOT run MadGraph, cross-section, detector-acceptance or event-yield work.
- Do NOT call a point "valid", "mixed" or "photonic" on your own judgement.
- Do NOT report a number you did not read out of tool output.

## Environment

```
cd /home/fabian/atlas_dihiggs/dihiggs
PY=/home/fabian/atlas_dihiggs/.venv/bin/python
RD=runs/discovery_v1
```

## Coordinates

The search space is `(mH_GeV, mA_GeV, tan_beta, X, Q)`. You never pass `M2` or
`lambda6`; they are derived:

```
lambda6 = X / tan_beta
M2_GeV2 = mH_GeV^2 - Q / tan_beta^2
mHp     = mA
```

Fixed and not yours to touch: `m_h = 125.20`, Yukawa type 1,
`sin(beta-alpha) = 1.0`, `lambda7 = 0`.

## Master calibration already established (deterministic, reproducible)

- `lambda1 = m_h^2/v^2 - 1.5 * X` exactly. Positivity needs `X < ~0.172`.
- Theory gates respond to the invariant `Q = (mH^2 - M2) * tan_beta^2`, not to
  `M2`. Unitarity/perturbativity degrade above `|Q| ~ 3e5`.
- `ctau ∝ tan_beta^2` at fixed `Q`. A theory-valid **mixed** anchor exists at
  `mH=200, mA=500, tan_beta=4.54e6, X=0.10, Q=0` giving `ctau = 700.0 mm`.
- The **photonic** pattern (`br_gammagamma + br_Zgamma ~ 1.0`) has so far only
  been observed where unitarity/perturbativity FAIL. Whether any theory-valid
  photonic point exists is an OPEN SEARCH QUESTION, not a settled negative.

Treat all of the above as telemetry to search against, never as a boundary you
must respect. The frozen gates decide.

## Commands (each one deduplicates, evaluates, ledgers and archives atomically)

```bash
$PY -m search_discovery envelope --run-dir $RD          # active tier-2 envelope
$PY -m search_discovery explore  --run-dir $RD --worker-id A1 --seed 11 --count 64 --mode sobol
$PY -m search_discovery refine   --run-dir $RD --worker-id C1 --parent <candidate_id> \
                                 --seed 21 --count 16 --radius 0.05 --target mixed
$PY -m search_discovery optimize --run-dir $RD --worker-id B1 --seed 31 --target photonic \
                                 --generations 12 --sigma0 0.25 [--start-candidate <id>]
$PY -m search_discovery validate-family --run-dir $RD --anchor-candidate <candidate_id>
$PY -m search_discovery summary  --run-dir $RD
```

`explore` modes: `sobol` (scrambled Sobol' coverage) and `random` (independent
log-uniform stream). Two workers in the same strategy MUST use different seeds
and should use different modes, so their trajectories stay independent.

`validate-family` re-derives the 150 / 200 / 250 GeV members by the frozen
`Q, S` continuation, holding `X` and `tan_beta` fixed, evaluates all three, and
offers the result to the Top-5 archive. It only accepts a `mH = 200` anchor
(the command rewrites `mH` to 200 and keeps your other coordinates).

## Output discipline

Report compactly. Per batch: attempts, valid fraction, provisional families
seen, best `ctau_mm`, best objective value, candidate IDs worth refining, and
whether you are still improving. Do NOT paste raw logs or full JSON dumps.
Quote candidate IDs so other stages can seed from them.

Your final message must be a short structured summary, not a narrative.

## Lineage selection

```bash
$PY -m search_discovery candidates --run-dir $RD --valid-only --target mixed --limit 20
$PY -m search_discovery candidates --run-dir $RD --family photonic --limit 20
```
Returns evaluated candidates best-first with their `candidate_id`, coordinates,
`ctau_mm`, provisional family and objective value. Use it to pick parents for
`refine` and start points for `optimize --start-candidate`.

## Anchor re-anchoring (IMPORTANT)

The scientific target is `ctau` **at the mH = 200 GeV anchor**, not at the mass a
candidate happened to be discovered at. `validate-family` re-derives the family
from a mH=200 anchor by the frozen `Q, S` continuation, so a candidate with
ctau = 968 mm at mH = 239 can have ctau_200 = 1263 mm and miss the window.

Two commands handle this:

```bash
# evaluate one explicit coordinate (global physics bounds still apply)
$PY -m search_discovery evaluate-point --run-dir $RD --mH 200 --mA 500 \
      --tan-beta 4.54e6 --X 0.10 --Q 0 --target mixed

# re-anchor a candidate to mH=200 and solve tan_beta for ctau_200 in [500,1000]
$PY -m search_discovery reanchor --run-dir $RD --candidate <candidate_id>
```

`reanchor` preserves `S = mA^2 - mH^2`, `Q` and `X`, sets `mH = 200`, and uses the
measured `ctau ~ tan_beta^2` scaling as its update step — but every iterate is
confirmed by the frozen evaluator. It reports the full trace. It typically
converges in 2 iterations. Use it BEFORE `validate-family` to produce finalists
that actually sit in the anchor window.

## Campaign checkpoint

```bash
$PY -m search_discovery checkpoint --run-dir $RD
```
Compact state: attempts and valid counts by strategy, valid fraction, duplicate
rate, gate-failure counts, best point per family, distinct basin count, archive
counts.

## Constrained continuation (feasible-manifold ascent)

For the "how far can a known-good mixed family be pushed toward higher photon
fraction without ever leaving validity" question, use `continue-lineage`. It
implements the FULL trust-region ascent algorithm as deterministic code — one
call performs many generations of: full 150/200/250 family evaluation, hard
feasibility check (all three masses theory-valid, same X, S > 501 GeV^2,
ctau_200 in [500,1000] mm via a re-solved tan_beta), step growth on a streak,
step shrink on failure, and branching when two different coordinates both
improve. **Do not manually perturb coordinates yourself for this task — call
the tool, it already does the algorithm.**

```bash
$PY -m search_discovery continue-lineage --run-dir $RD --lineage-id <ID> \
    --worker-id <ID> --seed-family <family_id> \
    --max-generations 30 --step0 0.10 --patience 6
# or continue from any already-evaluated candidate:
    --seed-candidate <candidate_id>
```

Returns JSON with `seed_pf200`, `best_pf200`, `enrichment_factor`, and one or
more `branches`, each with a `trace` of per-generation rows containing
`pf150/pf200/pf250`, `min_pf`, `mean_pf`, `ctau150/200/250`,
`theory150/200/250`, `X`, `tan_beta`, `lambda6`, `mA_GeV`, `Q_requested`,
`Q_effective_200`, `candidate_id`, `coordinates`, `direction`.

If the best branch is still improving when `max_generations` is reached (not
stopped by patience), you may extend it once:
`continue-lineage --seed-candidate <best final candidate_id> --lineage-id <ID>_ext ...`

Free coordinates are `mA_GeV`, `X`, `Q` only. `mH` is pinned to 200. `tan_beta`
is re-solved automatically, never a direction you choose.
