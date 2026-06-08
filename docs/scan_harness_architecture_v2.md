# Scan Harness Architecture — v2 (Modular Orchestrator)

> **Location**: `dihiggs/app/orchestrator/`
> **Supersedes**: `dihiggs/app/orchestrate_scans.py` (legacy monolith, kept for reference)
> **Standard library only** — no pandas, Ray, Optuna, MLflow, rich.

---

## 1. Purpose

The modular orchestrator (`dihiggs/app/orchestrator/`) manages 2HDM parameter
scans against compiled C++ binaries.  It supports two physics engines with
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
│   └── m2.py          Phys_M2BoundaryScan adapter
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
| `scan_axis` | `ScanAxis.LAMBDA1` or `ScanAxis.M2` |
| `executable_basename` | Default binary name |
| `build_command(exe, grid, fixed, csv)` | Build the 13-token CLI command |
| `axis_metadata()` | Authoritative dict: units, labels, m12/m12_sq note |
| `expected_csv_columns()` | Partial column list for validation |

### Lambda1Engine (PhysScanWithFixings)

```
./PhysScanWithFixings  mphi_min mphi_max N_mphi  lam1_min lam1_max N_lam1  mA sin_ba tan_beta lambda6 lambda7  output.csv
```

- Second triplet (positions 4–6): **lambda_1** (dimensionless)
- `fixed.lambda1` must be **None** (lambda1 is the scan axis)

### M2Engine (Phys_M2BoundaryScan)

```
./Phys_M2BoundaryScan  mphi_min mphi_max N_mphi  M2_min M2_max N_M2  mA sin_ba tan_beta lambda6 lambda7  output.csv
```

- Second triplet (positions 4–6): **M^2 = m12_sq** (units: **GeV^2**)
- `fixed.lambda1` must be **set** (lambda1 is a fixed constant here)

---

## 4. Axis Semantics Contract

> **Do not conflate lambda1, M2, m12, and m12_sq.**

| Quantity | Symbol | Units | Engine |
|----------|--------|-------|--------|
| Quartic coupling | lambda_1 | dimensionless | Lambda1Engine (scan axis) |
| Soft-breaking parameter | M^2 = m12_sq | GeV^2 | M2Engine (scan axis) |
| Fixed coupling | lambda_1 | dimensionless | M2Engine (fixed constant) |

**Historical ambiguity**: Some CSV files produced by `Phys_M2BoundaryScan` use
the column name `m12` while actually storing **GeV^2** values (`m12_sq`).  The
`axis_metadata()` returned by each engine documents this unambiguously and is
written into every `scan_meta.json` and `run_manifest.json`.

**Rule**: Never infer axis units from CSV column names alone.  Always read
`axis_metadata` from the scan metadata.

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
         [lambda1_fixed]})[:16]
```

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
| `ScanAxis` enum (LAMBDA1 / M2) | `engines/base.py` |
| Engine-specific grid signatures | `grid.py` |
| Engine-specific command builders | `engines/lambda1.py`, `engines/m2.py` |
| Engine-specific `axis_metadata()` in every JSON | `manifest.py` |
| M2 engine support | `engines/m2.py` |
| `--engine lambda1|m2` CLI flag | `cli.py` |
| Generic `--axis-min/max/n-axis` flags | `cli.py` |
| Legacy `--lam1-*` flags kept as aliases | `cli.py` |
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
