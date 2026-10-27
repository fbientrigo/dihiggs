# Quickstart

The current maintained core has two canonical production paths:

- `Lambda1EvaluatorV2`: explicit `(mH, lambda1_target)` inputs;
- `DihiggsPointV2Evaluator`: rectangular `(mH, M2)` inputs with reconstructed
  `lambda1`.

The active SM-like Higgs-mass convention for new production is **125.20 GeV**.
Do not silently reuse historical 125.13, 125.09, or 125.0 conventions.

## Build

From the repository root, build only the library and two maintained evaluators
first. This avoids the separate HiggsTools path required by legacy/HBHS
executables and does not remove or clean existing build outputs.

```bash
make -C 2hdmc lib -j2
make -C dihiggs app/Lambda1EvaluatorV2 app/DihiggsPointV2Evaluator -j2
```

The canonical binaries are written under `dihiggs/app/`.

The full `make -C dihiggs -j2` target builds all configured executables and
may require the separately installed HiggsTools library. Use it only when
those legacy/HBHS programs are needed; no clean target is required for the
canonical evaluator workflow.

## Lambda1 v2 through the orchestrator

```bash
python -m dihiggs.app.orchestrator \
  --engine lambda1_v2 \
  --campaign lambda1_docs_smoke \
  --outdir /tmp/dihiggs_output \
  --mH-min 130 --mH-max 130 --n-mH 1 \
  --axis-min 1.0 --axis-max 1.0 --n-axis 1 \
  --mA 300 --mHp 300 \
  --mh 125.20 \
  --sin-ba 0.995 \
  --lambda6 0.1 --lambda7 0.0 \
  --tanbeta 50
```

## M2 v2 through the orchestrator

```bash
python -m dihiggs.app.orchestrator \
  --engine m2 \
  --exec ./dihiggs/app/DihiggsPointV2Evaluator \
  --campaign m2_docs_smoke \
  --outdir /tmp/dihiggs_output \
  --mH-min 130 --mH-max 130 --n-mH 1 \
  --axis-min 15000 --axis-max 15000 --n-axis 1 \
  --mA 300 --mHp 300 \
  --mh 125.20 \
  --yukawa-type 1 \
  --sin-ba 0.995 \
  --lambda6 0.1 --lambda7 0.0 \
  --tanbeta 50
```

## Direct M2 evaluator smoke test

```bash
dihiggs/app/DihiggsPointV2Evaluator \
  --campaign-id smoke \
  --run-id m2 \
  --mh 125.20 \
  --mH-min 130 --mH-max 130 --n-mH 1 \
  --mA 300 --mHp 300 \
  --yukawa-type 1 \
  --sin-ba 0.995 \
  --tan-beta 50 \
  --M2-min 15000 --M2-max 15000 --n-M2 1 \
  --lambda6 0.1 --lambda7 0 \
  --output /tmp/dihiggs-point-v2.csv
```

## Validate the repository

```bash
python -m pytest -q
```

For exact schema semantics and current operational limits, read
[Canonical evaluators v2](contracts/canonical_evaluators_v2.md) and
[Operating status v2](OPERATING_STATUS_V2.md).
