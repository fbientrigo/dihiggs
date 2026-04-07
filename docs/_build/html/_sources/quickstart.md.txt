# Quickstart

This guide shows the fastest way to run the scan orchestrator.

## Minimal Setup

1. Create and activate a Python environment.
2. Ensure the C++ executable exists (default: `./PhysScanWithFixings`).
3. Choose an output root outside the repo (CERNBox-synced path recommended).

Example:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r docs/requirements.txt
```

## Dry Run (No C++ execution)

This creates folders and `run_manifest.json`, prints planned commands, and does not run the binary.

```bash
python dihiggs/app/orchestrate_scans.py \
  --exec ./PhysScanWithFixings \
  --campaign scan_docs \
  --tanbeta 10000 \
  --n-lam1 10 \
  --dry-run
```

## Tiny Real Run

This executes one small point set for quick validation.

```bash
python dihiggs/app/orchestrate_scans.py \
  --exec ./PhysScanWithFixings \
  --campaign scan_docs \
  --tanbeta 10000 \
  --n-lam1 10 \
  --mphi-min 130 --mphi-max 132 --n-mphi 2 \
  --lam1-min 0 --lam1-max 1
```

Outputs are written under:

```text
<outdir>/<lake-name>/campaign=<campaign>/fixed_.../run_.../
```
