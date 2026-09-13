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

- `PhysScanWithFixings`;
- `lambda1_legacy` and compatibility alias `lambda1`;
- adaptive/ML workflows that predate the v2 closure;
- frozen `autoresearch/` components;
- historical campaigns and quarantine scripts.

Do not use legacy artifacts to make new production claims unless a current
contract explicitly promotes that path.
