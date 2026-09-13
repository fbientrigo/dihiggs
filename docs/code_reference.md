# Code reference

This page maps the maintained execution surface to source paths.

## Canonical C++ evaluators

- `dihiggs/src/Lambda1EvaluatorV2.cpp` → `dihiggs/app/Lambda1EvaluatorV2`
- `dihiggs/src/DihiggsPointV2Evaluator.cpp` → `dihiggs/app/DihiggsPointV2Evaluator`
- `dihiggs/src/Phys_M2BandTracker.cpp` → experimental bounded-pilot helper

## Python orchestration

Primary package:

- `dihiggs/app/orchestrator/__main__.py`
- `dihiggs/app/orchestrator/cli.py`
- `dihiggs/app/orchestrator/lambda1_v2.py`
- `dihiggs/app/orchestrator/engines/m2.py`
- `dihiggs/app/orchestrator/physics_authority.py`
- `dihiggs/app/orchestrator/runner.py`
- `dihiggs/app/orchestrator/manifest.py`

The package exports `ScanRunner`, `ScanGrid`, `FixedParams`,
`TaskSpec`, `TaskResult`, `Lambda1Engine`, `M2Engine`, and
`run_lambda1_v2`.

## Physics contract

- `conventions/physics_conventions.yaml`: normative machine contract;
- `docs/PROJECT_PHYSICS_AUTHORITY.md`: human-readable interpretation;
- `docs/contracts/canonical_evaluators_v2.md`: evaluator schema contract.

Browse the current source tree on
[GitHub](https://github.com/fbientrigo/dihiggs/tree/main).
