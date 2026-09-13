# Quickstart

The current maintained core has two canonical production paths:

- `Lambda1EvaluatorV2`: explicit `(mH, lambda1_target)` inputs;
- `DihiggsPointV2Evaluator`: rectangular `(mH, M2)` inputs with reconstructed
  `lambda1`.

The active SM-like Higgs-mass convention for new production is **125.13 GeV**.
Do not silently reuse older 125.20, 125.09, or 125.0 conventions.

## Build

From the repository root:

```bash
make -C 2hdmc -j2
make -C dihiggs clean
make -C dihiggs -j2
```

The canonical binaries are written under `dihiggs/app/`.

## Lambda1 v2 through the orchestrator

```bash
python -m dihiggs.app.orchestrator \
  --engine lambda1_v2 \
  --campaign lambda1_docs_smoke \
  --outdir /tmp/dihiggs_output \
  --mH-min 130 --mH-max 130 --n-mH 1 \
  --axis-min 1.0 --axis-max 1.0 --n-axis 1 \
  --mA 300 --mHp 300 \
  --mh 125.13 \
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
  --mh 125.13 \
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
  --mh 125.13 \
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
