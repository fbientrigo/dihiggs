# 2HDMC / DiHiggs evaluation core

This repository is an executable 2HDMC evaluation and scan core for a generic
2HDM. It produces row-preserving theory evaluations, reconstructed quartics,
decay widths, branching ratios, and ctau_mm for bounded engineering pilots.
It is not a maintained ML or campaign-results repository.

## Canonical evaluators

- dihiggs/app/Lambda1EvaluatorV2 from dihiggs/src/Lambda1EvaluatorV2.cpp:
  explicit (mH, lambda1_target) inputs, schema dihiggs.lambda1.v2, and
  THDM::set_param_phys_lam1.
- dihiggs/app/DihiggsPointV2Evaluator from
  dihiggs/src/DihiggsPointV2Evaluator.cpp: rectangular (mH, M2) inputs,
  schema dihiggs.point.v2, and THDM::set_param_phys.
- dihiggs/app/Phys_M2BandTracker: experimental bounded-pilot boundary helper;
  its intervals are not canonical point-production evidence.

The frozen contract is docs/contracts/canonical_evaluators_v2.md.

## Canonical schemas

dihiggs.lambda1.v2 preserves one row for every attempted input, including
malformed, construction-failing, theory-rejected, and accepted rows. It stores
raw input lexemes, exact Float64 values, reconstructed quartics and m12_sq,
theory flags, selected widths, branching ratios, width_ok, and ctau_mm.

dihiggs.point.v2 preserves every rectangular grid point. It records
M2 = m12_sq / (sin(beta) * cos(beta)) and
m12_sq = M2 * sin(beta) * cos(beta). lambda1 is reconstructed output.
Experimental fields are deliberately unevaluated.

## Theory acceptance semantics

construction_ok is the exact 2HDMC parameter-construction result.
positivity, unitarity, and perturbativity are theory predicates;
TRIPLE_OK/triple_ok_legacy is theory-only. theory_ok_v1 currently equals
those three predicates. No experimental acceptance is inferred from these
flags.

## Build prerequisites

Linux or WSL, a C++17 compiler, GNU Make, GSL, and the repository's patched
2HDMC checkout under 2hdmc/ are required. HiggsTools is not needed for the
three v2 2HDMC-only executables.

## Build commands

    make -C 2hdmc -j2
    make -C dihiggs clean
    make -C dihiggs -j2

The v2 binaries are written under dihiggs/app/.

## Lambda1 v2 smoke test

    tmp=$(mktemp -d)
    printf '%s\n' 'point_id,mh_gev,mH_gev,mA_gev,mHp_gev,sin_beta_minus_alpha,tan_beta,lambda1_target,lambda6_input,lambda7_input' 'p0,125.13,130,300,300,0.995,50,1.0,0.1,0.0' > "$tmp/input.csv"
    dihiggs/app/Lambda1EvaluatorV2 "$tmp/input.csv" "$tmp/output.csv"

## M2 v2 smoke test

    dihiggs/app/DihiggsPointV2Evaluator --campaign-id smoke --run-id m2 --mh 125.13 --mH-min 130 --mH-max 130 --n-mH 1 --mA 300 --mHp 300 --yukawa-type 1 --sin-ba 0.995 --tan-beta 50 --M2-min 15000 --M2-max 15000 --n-M2 1 --lambda6 0.1 --lambda7 0 --output /tmp/dihiggs-point-v2.csv

## Orchestrator examples

Canonical lambda1 v2 orchestration generates one exact input CSV and validates
the output schema and row count:

    python -m dihiggs.app.orchestrator --engine lambda1_v2 --campaign lambda1_pilot --mH-min 130 --mH-max 290 --n-mH 3 --axis-min 0 --axis-max 12 --n-axis 4 --mA 300 --mHp 300 --mh 125.13 --sin-ba 0.995 --lambda6 0.1 --lambda7 0 --tanbeta 50

Canonical M2 orchestration uses named flags and reconstructs lambda1:

    python -m dihiggs.app.orchestrator --engine m2 --exec ./dihiggs/app/DihiggsPointV2Evaluator --mH-min 130 --mH-max 290 --n-mH 3 --axis-min 15000 --axis-max 18000 --n-axis 4 --mA 300 --mHp 300 --mh 125.13 --yukawa-type 1 --sin-ba 0.995 --lambda6 0.1 --lambda7 0 --tanbeta 50

m2_tracker is explicitly experimental and bounded-pilot only.
lambda1_legacy (or the compatibility alias lambda1) invokes
PhysScanWithFixings; it is replay-only and not an LLP production path.

## Output interpretation and lifetime units

Masses and widths use GeV; M2 and m12_sq use GeV2. ctau_mm is the proper
decay length in millimetres, derived from hbar*c = 1.973269804e-13 GeV mm
when width_ok is true. M2 must not be read as m12_sq.

## Tracker experimental status

Phys_M2BandTracker searches for intervals in bounded pilot domains. Its
interval output is a boundary-search artifact, not final global boundary
evidence, and it does not provide a fixed-input-lambda1 evaluation.

## Legacy and frozen components

PhysScanWithFixings, the autoresearch/ tree, ML notebooks, historical
campaigns, and quarantine/replay scripts remain available only where their
historical or compatibility role is explicit. Autoresearch is frozen and no
canonical code imports it.

## Tests

    python -m pytest -q

Focused C++ and orchestrator tests are documented in
docs/OPERATING_STATUS_V2.md. Optional, frozen, legacy, and quarantine tests
skip explicitly when their components are absent.

## Documentation map

- docs/contracts/canonical_evaluators_v2.md
- docs/OPERATING_STATUS_V2.md
- docs/autoresearch_frozen.md
- docs/handoffs/REPO_TO_TEX_NOTES_V3.md
- docs/verification/

## Current non-claims

The maintained core currently claims no LHS production, SHAP, Bayesian
optimization, active learning, surrogate training, complete seven-dimensional
coverage, experimental acceptance, full campaign result, HB/HS gate, STU gate,
recast integration, or maintained paper result.
