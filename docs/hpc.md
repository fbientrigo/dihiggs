# Long-run and HPC operation

The orchestrator is designed to leave CPU headroom by default. On a host that
also provides SSH or interactive services, this is safer than consuming every
logical core.

## CPU control

Use:

```bash
python -m dihiggs.app.orchestrator --threads N ...
```

Only use `--all-cores` on a dedicated batch node where starving interactive
services is acceptable.

## Dry run first

Before a long campaign, render the planned layout and commands without executing
the C++ evaluator:

```bash
python -m dihiggs.app.orchestrator --dry-run ...
```

Inspect the produced manifest and command lines before the full run.

## Reproducibility

For every production campaign preserve:

- repository commit SHA;
- engine and schema version;
- physics convention ID;
- exact decimal Higgs mass text;
- fixed parameters and scan axes;
- run manifest;
- output CSVs and validation records.

Do not mix frozen historical artifacts with newly produced rows under a new
physics convention.
