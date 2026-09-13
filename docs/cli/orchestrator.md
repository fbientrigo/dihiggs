# Orchestrator CLI

The supported entry point is:

```bash
python -m dihiggs.app.orchestrator
```

The old `dihiggs/app/orchestrate_scans.py` wrapper remains only for compatibility.

## Engine status

| Engine | Status | Meaning |
|---|---|---|
| `lambda1_v2` | canonical | `Lambda1EvaluatorV2`, explicit `lambda1_target` |
| `m2` | canonical | `DihiggsPointV2Evaluator`, explicit `M2`, reconstructed `lambda1` |
| `m2_tracker` | experimental | bounded-pilot `Phys_M2BandTracker` boundary helper |
| `lambda1_legacy` / `lambda1` | replay only | compatibility adapter for `PhysScanWithFixings` |
| `gen_fixings` | calibration / compatibility | fixed-input calibration path |

New LLP lifetime production must use a canonical v2 producer.

## Common runtime flags

- `--engine`: select the execution path.
- `--exec`: override the selected C++ executable.
- `--threads N`: set an explicit OpenMP thread cap.
- `--all-cores`: disable the automatic CPU-headroom cap.
- `--dry-run`: create layout/manifests without executing C++.
- `--force`: overwrite existing output where supported.
- `--timeout`: per-task subprocess timeout.
- `--campaign`, `--outdir`: run identity and output root.

By default the orchestrator leaves CPU headroom instead of pinning every logical
core. This is intentional: long scans must not starve the host or SSH service.

## Physics coordinates

The canonical paths distinguish the two scan axes explicitly:

- `lambda1_v2`: the axis is `lambda1_target`;
- `m2`: the axis is `M2 = m12_sq / (sin(beta) cos(beta))`.

`M2`, `m12_sq`, and `m22_sq` are different quantities and must not be
interchanged by name similarity.

For current examples, use the repository [README](https://github.com/fbientrigo/dihiggs/blob/main/README.md)
and [Quickstart](../quickstart.md).
