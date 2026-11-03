# Legacy and experimental components

The repository contains historical and experimental code because exact replay
and provenance remain important. Presence in the tree does **not** make a path
canonical for new production.

## Canonical

- `Lambda1EvaluatorV2`;
- `DihiggsPointV2Evaluator`;
- `python -m dihiggs.app.orchestrator --engine lambda1_v2`;
- `python -m dihiggs.app.orchestrator --engine m2`.

## Experimental

- `Phys_M2BandTracker` / `m2_tracker`: bounded-pilot boundary helper.
  Its intervals are not canonical point-production evidence.

## Replay or compatibility only

- `PhysLam1Scan`: frozen oracle behind `tests/test_golden_lambda1.py`;
- ML workflows that predate the v2 closure;
- frozen `autoresearch/` components;
- historical campaigns and quarantine scripts.

## Removed

`PhysScanWithFixings`, `PhysParamScan`, `GenScanWithFixings`, the
`lambda1_legacy`/`gen_fixings` orchestrator engines, `orchestrate_scans.py`,
the adaptive explorers and the christopher campaign scripts were removed.
Reproduce old campaigns with `git checkout legacy-lambda1-final`.

Do not use legacy artifacts to make new production claims unless a current
contract explicitly promotes that path.
