# T2 → T4 campaign status

Campaign directory: `dihiggs/campaigns/t2_t4_20261008T043858Z/`.
The T1 reproduction uses `mh=125.13 GeV` and `v=246 GeV`; the active
`mh=125.20 GeV` convention has a separate 90-point table. No benchmark is
classified as physically validated. The 300 T1-convention and 90 active-
convention points are `CANDIDATE_VALID` under the recorded tree-level
certificates; this label does not mean benchmark certification.

The work is partial. The pilot-wide actual-calculator comparison completed all 390 points through legacy, explicit, and analytic routes at 50/80/120-digit input precision. Analytic constructor agreement passed all 1170 route/precision records; legacy and explicit routes failed their recorded per-observable family-agreement checks, with exact reasons and residuals preserved. The independent native diphoton component reconstruction passed; the leading-order appendix versus running-mass/QCD width comparison remains `INCONCLUSIVE` as a scheme diagnostic. Full stationary-point enumeration and certified boundary refinement remain incomplete. The additive `--engine t4` CLI supports legacy, explicit, and analytic modes while preserving existing engine defaults and schemas. Nondegenerate reversed CP-even states are mapped to calculator indices with raw/mapped phase provenance; exactly degenerate states remain `INCONCLUSIVE` with widths `NOT_TESTED`. Displaced-photon selection and experimental exclusions remain `NOT_TESTED`.

The precision audit includes stable forward construction and arbitrary-
precision inversion at 50/80/120 digits, naive binary64 inversion, decimal
serialization loss, a lambda7 perturbation diagnostic, and three actual
2HDMC routes at one high-tanβ anchor. Its measured numerical envelopes are
anchor-specific; they do not cover route agreement across the pilot. General
BFB uses exact rational witnesses and a sufficient certificate for the
generated aligned family, but does not certify every arbitrary quartic
potential. Competing-extrema enumeration, exact treatment of all degenerate
cases, and physical-boundary refinement remain uncompleted, so candidate
labels are limited to the checks recorded for each sampled point.

The seven deliverables are:

1. `reproducibility_manifest.json` — source/artifact checksums, repository
   revisions and dirty-diff fingerprints, runtime dependencies and build
   executable hashes.
2. `numerical_precision_audit.json` — decimal constructor precision,
   complete scalar-scattering spectra at 50/80/120 digits, convergence and
   backend route comparisons; records the separate active-mass table.
3. `physical_constraints_report.json` and
   `physical_constraints_report.jsonl` — status summary and complete
   per-point constraints, values, and failure reasons.
4. `pilot_scan_tables.csv` — T1-convention pilot; the separate active table is
   `pilot_scan_tables_mh125_20.csv` with detailed rows in
   `active_mh125_20_constraints.jsonl`.
5. `physical_boundary_analysis.json` — sampled status transitions and
   limitations on refinement/certification.
6. `validated_benchmarks.json` — intentionally empty because no point meets
   the complete certification and phenomenology requirements.
7. `T4_integration_report.md` — constructor route behavior, backend
   comparisons, decay diagnostics and remaining limitations. The machine
   readable measurements are in `T4_integration_report.json`.

Intermediate resumable/checkpoint files include the point batches under
`t2_lightweight/`, `physical_constraints_checkpoint.jsonl`, and the
`*_pre_*` trajectory snapshots. `baseline_checksums.json` records the
checksummed T1 and source baseline preserved before these edits.

## Reproduction commands

From `dihiggs/`, using the isolated campaign Python environment:

```sh
/tmp/dihiggs-t2-t4-venv/bin/python campaigns/t2_t4_20261008T043858Z/scripts/run_t2_lightweight.py
/tmp/dihiggs-t2-t4-venv/bin/python campaigns/t2_t4_20261008T043858Z/scripts/run_t2_constraints.py
/tmp/dihiggs-t2-t4-venv/bin/python -m pytest -q tests/test_aligned_lambda7.py tests/test_aligned_constraints.py tests/test_analytic_bridge.py tests/test_t2_lightweight_checkpoint.py tests/test_t2_active_table.py
make -C campaigns/t2_t4_20261008T043858Z/build_analytic/2hdmc analytic-bridge
/tmp/dihiggs-t2-t4-venv/bin/python campaigns/t2_t4_20261008T043858Z/scripts/run_t4_bridge.py
```

Latest verification: the focused Python suite reports 26 passed; the generic
cache regression passed; and `CalcRoundTrip` passed all four legacy cases.
The T4 runner invoked all nine route/precision combinations: analytic passed
all three, while explicit and legacy failed for the exact reasons retained in
`T4_integration_report.json`.

The isolated C++ regressions were compiled against the campaign library:

```sh
g++ -std=c++11 -Wall -O3 -fopenmp -I campaigns/t2_t4_20261008T043858Z/build_analytic/2hdmc/src campaigns/t2_t4_20261008T043858Z/scripts/test_gen_cache.cpp -L campaigns/t2_t4_20261008T043858Z/build_analytic/2hdmc/lib -l2HDMC -lgsl -lgslcblas -lm -o /tmp/t2_t4_test_gen_cache
/tmp/t2_t4_test_gen_cache
g++ -std=c++11 -Wall -O3 -fopenmp -I campaigns/t2_t4_20261008T043858Z/build_analytic/2hdmc/src campaigns/t2_t4_20261008T043858Z/build_analytic/2hdmc/src/CalcRoundTrip.cpp -L campaigns/t2_t4_20261008T043858Z/build_analytic/2hdmc/lib -l2HDMC -lgsl -lgslcblas -lm -o /tmp/t2_t4_CalcRoundTrip
/tmp/t2_t4_CalcRoundTrip
```

The canonical `2hdmc/lib/lib2HDMC.a` and additive `2hdmc/analytic-bridge` were rebuilt with all five library translation units for ABI consistency. The existing `DihiggsPointV2Evaluator` and `Lambda1EvaluatorV2` executables retained their baseline SHA256 values; existing legacy calculator executables were not relinked. The campaign-local relinked evaluator was compared at historical 450/500 GeV and an active-mass anchor, with exact CSV equality.


The completed actual-calculator comparison is recorded in `calculator_scan_agreement.jsonl` (3510 records: 390 points × 3 precisions × 3 routes), with replay identity/checksum in `calculator_scan_checkpoint.json` and aggregation in `calculator_scan_summary.json`. The T1 and active-mass slices remain separate. The scan's constructor-agreement `FAIL` statuses report residuals to the stable-family target; they do not mean the point is physically invalid. Analytic states use the requested CP-even eigenstate and phase map; explicit and legacy failures retain their measured soft-term/quartic/scale residuals. Earlier passes and interrupted pre-audit attempts are preserved under `snapshots/` or labeled `.pre_`/`.aborted_`; only the current identity-matching output is summarized. Certified extrema enumeration and boundary refinement remain incomplete.

The locked v4 scan identity is `cfe361dffbdd5a8eac7cabc8755a00944edb9ae2264be7fd494225587cc3ca71`; the previous pre-lock v3 dataset remains separately checksummed under `snapshots/pre_scan_lock_v3_20261008T/`. The current runner acquires the lock before loading or mutating scan files and rejects a concurrent writer with the owner PID and scan identity. A fresh 390-point v4 pass and a strict no-new-calculator-call replay both completed successfully.

The reproducible v4 scan and strict replay commands (run separately; both acquire the same nonblocking lock) are:

```sh
prlimit --as=9663676416 --cpu=900 -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    /tmp/dihiggs-t2-t4-venv/bin/python campaigns/t2_t4_20261008T043858Z/scripts/run_t4_scan.py --batch-size 10
```

The first invocation produced the current v4 records; the second validated all record identities, inputs, schemas, and payload checksums and emitted the summary without additional calculator calls. The lock helper's concurrent-writer behavior is covered by `test_scan_lock_blocks_a_concurrent_writer_and_records_owner`.
