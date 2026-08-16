# Scan Harness Architecture — v2 (Modular Orchestrator)

> **Location**: `dihiggs/app/orchestrator/`
> **Supersedes**: `dihiggs/app/orchestrate_scans.py` (legacy monolith, kept for reference)
> **Standard library only** — no pandas, Ray, Optuna, MLflow, rich.

---

## 1. Purpose

The modular orchestrator (`dihiggs/app/orchestrator/`) manages 2HDM parameter
scans against compiled C++ binaries.  It supports three physics engines with
explicitly distinct axis semantics, preventing the historical confusion between
`lambda1`, `M2`, `m12`, and `m12_sq`.

---

## 2. File Map

```
dihiggs/app/orchestrator/
├── __init__.py        Public API re-exports
├── __main__.py        python -m dihiggs.app.orchestrator
├── cli.py             Argument parser + main() entrypoint
├── runner.py          ScanRunner — main execution loop
├── engines/
│   ├── __init__.py
│   ├── base.py        EngineAdapter protocol + ScanAxis enum
│   ├── lambda1.py     PhysScanWithFixings adapter
│   ├── m2.py          Phys_M2BoundaryScan adapter
│   └── gen_fixings.py GenScanWithFixings adapter (Stage 2 calibration)
├── grid.py            ScanGrid dataclass + grid_signature()
├── models.py          FixedParams, TaskSpec, TaskResult dataclasses
├── layout.py          Data-lake path construction
├── manifest.py        run_manifest.json + scan_meta.json writers
├── resume.py          Resume / skip logic + cross-run CSV reuse
├── parsing.py         parse_cpp_stats() — Total Attempts, TRIPLE_OK_POINTS
└── io_utils.py        Atomic JSON write, JSONL append, git metadata
```

Tests live in `tests/test_orchestrator/`.

---

## 3. Engine Adapters

Each engine implements the `EngineAdapter` protocol:

| Method | Purpose |
|--------|---------|
| `engine_name` | Stable identifier included in grid signature |
| `scan_axis` | `ScanAxis.LAMBDA1`, `ScanAxis.M2`, or `ScanAxis.GEN_FIXINGS` |
| `executable_basename` | Default binary name |
| `build_command(exe, grid, fixed, csv)` | Build the CLI command (13-token positional for lambda1/m2; named flags for gen_fixings) |
| `axis_metadata()` | Authoritative dict: units, labels, m12/m12_sq note |
| `expected_csv_columns()` | Partial column list for validation |

### Lambda1Engine (PhysScanWithFixings)

```
./PhysScanWithFixings  mphi_min mphi_max N_mphi  lam1_min lam1_max N_lam1  mA sin_ba tan_beta lambda6 lambda7  output.csv
```

- Second triplet (positions 4–6): **lambda_1** (dimensionless)
- `fixed.lambda1` must be **None** (lambda1 is the scan axis)

### M2Engine (DihiggsPointV2Evaluator)

```
./DihiggsPointV2Evaluator --campaign-id ID --run-id ID --mh 125.20 \
  --mH-min MIN --mH-max MAX --n-mH N --mA MASS --mHp MASS --yukawa-type 1 \
  --sin-ba SBA --tan-beta TB --M2-min MIN --M2-max MAX --n-M2 N \
  --lambda6 L6 --lambda7 L7 --output output.csv
```

- M² is defined as `m12_sq/(sin(beta)*cos(beta))` in **GeV²**.
- lambda1 is reconstructed output only.

### GenFixingsEngine (GenScanWithFixings, "Stage 2" calibration)

```
./GenScanWithFixings --bronze-csv=<shard.csv> --output-csv=<validated.csv> \
    [--calibration-n=50] [--calibration-frac=0.10] [--rng-seed=0]
```

`Phys_M2BandTracker` remains an experimental, non-canonical boundary helper.
Its branch selection and refinement are validated only within the bounded pilot
domains recorded in `docs/verification/dihiggs_point_v2_verification_v1.*`.

- `scan_axis` = `ScanAxis.GEN_FIXINGS`: there is **no generated (m_phi, axis)
  grid**. The `ScanGrid` passed to `build_command` is a required-but-unused
  placeholder (its `mphi_*`/`axis_*` fields are ignored).
- The scan domain is a pre-generated **bronze shard CSV**
  (`chris/CalcLambda1ScanFixings` output, `fixed.bronze_shard_csv`), one
  shard = one `(tan_beta, m_A, lambda6, lambda7, sin_ba)` fixing tuple.
- `cli.py` derives `mA`/`sin_ba`/`lambda6`/`lambda7`/`tan_beta` from the
  bronze shard's first data row (not from `--mA`/`--sin-ba`/etc.), so
  run-dir naming and grid signatures can't silently diverge from the data
  being processed. `--mA`/`--sin-ba`/`--lambda6`/`--lambda7`/`--tanbeta` are
  ignored for this engine.
- Each row is calibrated against 2HDMC via a `±calibration-frac`
  `N=calibration-n` random search (`GenScanPointEvaluator`), producing a
  "validated/silver" CSV (legacy 29-column schema + calibration/stability/
  chris cross-check diagnostics).

---

## 4. Axis Semantics Contract

> **Do not conflate lambda1, M2, m12, and m12_sq.**

| Quantity | Symbol | Units | Engine |
|----------|--------|-------|--------|
| Quartic coupling | lambda_1 | dimensionless | Lambda1Engine (scan axis) |
| Soft-breaking parameter | M² = m12_sq/(sin(beta)cos(beta)) | GeV² | M2Engine (scan axis) |
| Diagnostic coupling | lambda_1 | dimensionless | M2Engine (reconstructed output) |

**Historical ambiguity**: Some CSV files produced by `Phys_M2BoundaryScan` use
the column name `m12` while actually storing **GeV^2** values (`m12_sq`).  The
`axis_metadata()` returned by each engine documents this unambiguously and is
written into every `scan_meta.json` and `run_manifest.json`.

### M2 vs m12_sq contract

* M2 is NOT identical to m12_sq.
* The exact relation is:
  `M2 = m12_sq / (sin_beta * cos_beta)`
  `m12_sq = M2 * sin_beta * cos_beta`
* Ensure that you interpret historical CSV columns correctly!

**Rule**: Never infer axis units from CSV column names alone.  Always read
`axis_metadata` from the scan metadata.

### gen_fixings: no generated axis

`ScanAxis.GEN_FIXINGS` is **not a numeric scan axis** — `axis_metadata()`
reports `axis_units: "n/a"`. The "axis" is which bronze shard CSV
(`fixed.bronze_shard_csv`) is being calibrated; the per-row `(tan_beta, m_A,
lambda6, lambda7, sin_ba)` values come from the shard itself, not from a
generated grid.

---

## 5. Grid Signature

```python
from dihiggs.app.orchestrator.grid import grid_signature
sig = grid_signature(engine_name, grid, fixed)   # 16 hex chars
```

The engine name is **included in the hash** to prevent resume collisions
between lambda1 and M2 scans that share the same numeric ranges.

```
SHA-256({engine, mphi_min, mphi_max, n_mphi,
         axis_min, axis_max, n_axis,
         mA, sin_ba, tan_beta, lambda6, lambda7,
         [lambda1_fixed],
         [bronze_shard_csv, calibration_n, calibration_frac]})[:16]
```

For `gen_fixings`, the bronze shard path and calibration config are folded
into the hash so distinct shards/configs don't collide on the dedup key.

---

## 6. Data-Lake Folder Layout

```
<outdir>/<lake_name>/
    campaign=<campaign>/
        fixed_sinba=<...>_l6=<...>_l7=<...>_mA=<...>/   ← fixed_dir
            <run_name>/                                   ← run_dir
                run_manifest.json
                orchestrator.log
                task_summary.jsonl
                tb_<NNNNN>/                              ← per-tanbeta
                    scan_tb_<tag>.csv
                    scan_meta.json
                    stdout.log   (on failure)
                    stderr.log   (on failure)
```

Preserved 1:1 from the legacy orchestrator.

---

## 7. Preserved Legacy Features

| Feature | Module |
|---------|--------|
| campaign/run/fixed-param hierarchy | `layout.py` |
| `run_manifest.json` (initial + summary) | `manifest.py` |
| `task_summary.jsonl` (one record per task) | `runner.py` |
| `scan_meta.json` per task | `manifest.py` |
| Dry-run (manifests written, no subprocess) | `runner.py` |
| Force overwrite | `runner.py` |
| Resume by exact grid signature | `resume.py` |
| Cross-run CSV reuse (`--resume-scope=fixed`) | `resume.py` |
| Materialize previous matching CSV | `resume.py` |
| OMP_NUM_THREADS injection | `runner.py` |
| stdout/stderr capture on failure | `runner.py` |
| `Total Attempts` + `TRIPLE_OK_POINTS` parsing | `parsing.py` |
| Atomic JSON writes (tmp+rename) | `io_utils.py` |
| Git metadata in manifests | `io_utils.py` |
| `tb_NNNNN` folder naming | `layout.py` |

---

## 8. New Features

| Feature | Where |
|---------|-------|
| Engine adapter protocol (`EngineAdapter`) | `engines/base.py` |
| `ScanAxis` enum (LAMBDA1 / M2 / GEN_FIXINGS) | `engines/base.py` |
| Engine-specific grid signatures | `grid.py` |
| Engine-specific command builders | `engines/lambda1.py`, `engines/m2.py`, `engines/gen_fixings.py` |
| Engine-specific `axis_metadata()` in every JSON | `manifest.py` |
| M2 engine support | `engines/m2.py` |
| gen_fixings ("Stage 2" calibration) engine support | `engines/gen_fixings.py` |
| `--engine lambda1|m2|m2_tracker|gen_fixings` CLI flag | `cli.py` |
| Generic `--axis-min/max/n-axis` flags | `cli.py` |
| Legacy `--lam1-*` flags kept as aliases | `cli.py` |
| `--bronze-csv/--calibration-n/--calibration-frac/--rng-seed` flags | `cli.py` |
| Fake-executable tests (no C++ compile needed) | `tests/test_orchestrator/` |

---

## 9. CLI Quick Reference

```bash
# Lambda1 scan (backward compatible)
python -m dihiggs.app.orchestrator \
    --exec ./PhysScanWithFixings \
    --campaign scan_999 \
    --sin-ba 0.995 --lambda6 0.10 --lambda7 0.00 --mA 300 \
    --tanbeta 10000,15000,20000 \
    --lam1-min 0 --lam1-max 12 --n-lam1 666

# M2 boundary scan (new)
python -m dihiggs.app.orchestrator \
    --engine m2 \
    --exec ./Phys_M2BoundaryScan \
    --campaign m2_scan_001 \
    --sin-ba 1.0 --lambda6 0.001 --lambda7 0.0 --mA 500 \
    --lambda1 1.0 \
    --tanbeta 50.0 \
    --axis-min 0 --axis-max 500000 --n-axis 50

# gen_fixings calibration (Stage 2: validate a bronze shard against 2HDMC)
python -m dihiggs.app.orchestrator \
    --engine gen_fixings \
    --exec ./GenScanWithFixings \
    --campaign gen_fixings_001 \
    --bronze-csv chris/bronze/tb_10000_mA_300/shard.csv \
    --calibration-n 50 --calibration-frac 0.10 --rng-seed 0
# Note: --mA/--sin-ba/--lambda6/--lambda7/--tanbeta are ignored here -- they
# are derived from the bronze shard's first row instead.

# Dry-run
python -m dihiggs.app.orchestrator --dry-run --engine lambda1 ...

# Force overwrite
python -m dihiggs.app.orchestrator --force ...
```

---

## 10. Backward Compatibility Note

The legacy `dihiggs/app/orchestrate_scans.py` is **not deleted**.  It remains
as-is for the autoresearch harness (`orchestrator_adapter.py`) which invokes it
via CLI subprocess.  Migration of autoresearch to the new orchestrator is a
separate task.

---

## 11. What Is Intentionally Not Here

- Adaptive search logic (handled by `autoresearch/harness/`)
- `autoresearch` modifications
- Heavy dependencies (pandas, Ray, Optuna, MLflow, rich)
- Per-tanbeta `--n-lam1-map` override (deferred to v1.1)

## 12. Known Gaps

- `gen_fixings` has no `tests/test_orchestrator/test_command_gen_fixings.py`
  yet, unlike `lambda1`/`m2` (see `test_command_lambda1.py`,
  `test_command_m2.py`). Follow-up work.
