# SEARCH_INTERFACE.md

Stable handoff for the later Orca workers. This interface is proposal-only;
the workers cannot redefine physics, classification, scoring, or archive state.

The immutable contract is
`docs/contracts/llp_benchmark_search_v1.json` (`dihiggs.llp_benchmark_search.v1`).
The append-only run artifact is `<run-dir>/ledger.jsonl`; the derived archive is
`<run-dir>/archive.json`.

## Allowed operations

```bash
python -m search_substrate contract
python -m search_substrate request-candidates --count 8 --seed 123 --worker-id orca-1
python -m search_substrate perturb --parent parent.json --count 8 --seed 456 --radius 0.05
python -m search_substrate evaluate --candidate candidate.json --run-dir runs/llp_search_v1 --evaluator dihiggs/app/DihiggsPointV2Evaluator
python -m search_substrate validate-family --family anchor_candidate.json --run-dir runs/llp_search_v1 --evaluator dihiggs/app/DihiggsPointV2Evaluator
python -m search_substrate summary --run-dir runs/llp_search_v1
```

All commands use structured JSON on stdout. A proposal contains only metadata
plus `parameters`: `mH_GeV`, `mA_GeV`, `M2_GeV2`, `tan_beta`, and `lambda6`.
Aliases are canonicalized before identity hashing. `X` is always derived as
`tan_beta * lambda6`; workers may not submit or redefine it. `mHp`, `m_h`,
`lambda7`, alignment, Yukawa type, family labels, scores, and archive fields
are forbidden proposal inputs.

`validate-family` accepts only an anchor at `mH=200 GeV`. It derives the 150,
200, and 250 GeV members by the frozen Q,S continuation and checks the same X.
Every member is a separate evaluator attempt and separate ledger record.

## Read-only information available to workers

- `summary` archive entries and aggregate ledger counts;
- individual `evaluation.json` records under `<run-dir>/attempts/`;
- failure stage/reason statistics derived from the ledger.

## Forbidden operations

Workers must not edit the contract, evaluator, classifier, score, or archive;
write directly to `archive.json`; delete or rewrite `ledger.jsonl`; pass
experimental acceptance/recast quantities; use pre-Yukawa-fix/L06 decay data;
or launch MadGraph, production cross-section, detector, or multi-agent search
workflows through this substrate.

Archive promotion is code-owned and fail-closed. In v1 it is blocked because
the scientific owner has not supplied authoritative numerical mixed/photonic
predicates and active-125.20 fixtures.

