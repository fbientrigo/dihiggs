# dihiggs

## PhysLam1Scan Executable Usage

`PhysLam1Scan` is a parameter scanning tool for Two-Higgs-Doublet Model (THDM) physics calculations. It performs a two-dimensional scan over m_φ and λ₁ parameters, evaluating physics constraints (positivity, unitarity, perturbativity) and computing decay widths and branching ratios.

### Prerequisites

Before building, ensure you have:
- **C++ compiler**: g++ with C++17 support
- **Build tool**: GNU make
- **2HDMC library**: Installed and accessible (typically in system library paths)
- **HiggsTools**: Required physics calculation library
- **Standard libraries**: Boost (optional, for advanced features)

### Build Command

Build the entire project, including `PhysLam1Scan`:

```bash
make -C dihiggs
```

This produces the executable at `dihiggs/app/PhysLam1Scan`.

### CLI Arguments

`PhysLam1Scan` accepts 12 positional command-line arguments:

| Argument # | Name | Type | Description |
|------------|------|------|-------------|
| 1 | `mphi_min` | double | Minimum m_φ value (GeV) |
| 2 | `mphi_max` | double | Maximum m_φ value (GeV) |
| 3 | `N_mphi` | int | Number of m_φ scan points |
| 4 | `lam1_min` | double | Minimum λ₁ value |
| 5 | `lam1_max` | double | Maximum λ₁ value |
| 6 | `N_lam1` | int | Number of λ₁ scan points |
| 7 | `mA` | double | Pseudoscalar mass (GeV); note: mₐ = m_H⁺ |
| 8 | `sin_ba` | double | sin(β-α) value, must be in [-1, 1] |
| 9 | `tan_beta` | double | tan(β) value, must be > 0 |
| 10 | `lambda6` | double | λ₆ coupling parameter |
| 11 | `lambda7` | double | λ₇ coupling parameter |
| 12 | `output.csv` | string | Path to output CSV file |

### Stdout Markers

`PhysLam1Scan` outputs two key markers to stdout:

- **`Total Attempts: <int>`** — Total number of parameter points evaluated during the scan
- **`TRIPLE_OK_POINTS <int>`** — Number of points passing all three physics constraints (positivity, unitarity, and perturbativity)

Example output:
```
Total Attempts: 15
TRIPLE_OK_POINTS 15
```

### Output CSV Schema

The output CSV contains 29 columns:

| Column # | Name | Type | Description |
|----------|------|------|-------------|
| 1 | `m_phi` | float | Scanned m_φ value (GeV) |
| 2 | `mA` | float | Pseudoscalar mass (GeV) |
| 3 | `alpha` | float | Mixing angle α (radians) |
| 4 | `beta` | float | Mixing angle β (radians) |
| 5 | `lambda6` | float | λ₆ coupling |
| 6 | `lambda7` | float | λ₇ coupling |
| 7 | `m12` | float | Mass parameter m₁₂² (GeV²) |
| 8 | `sin_ba` | float | sin(β-α) value |
| 9 | `tan_beta` | float | tan(β) value |
| 10 | `positivity_ok` | int | 1 if positivity constraint satisfied, 0 otherwise |
| 11 | `unitarity_ok` | int | 1 if unitarity constraint satisfied, 0 otherwise |
| 12 | `perturbativity_ok` | int | 1 if perturbativity constraint satisfied, 0 otherwise |
| 13 | `width_bb` | float | Decay width to b-quarks (GeV) |
| 14 | `width_tautau` | float | Decay width to τ leptons (GeV) |
| 15 | `width_WW` | float | Decay width to W bosons (GeV) |
| 16 | `width_ZZ` | float | Decay width to Z bosons (GeV) |
| 17 | `width_gaga` | float | Decay width to photons (GeV) |
| 18 | `width_Zga` | float | Decay width to Z and photon (GeV) |
| 19 | `width_gg` | float | Decay width to gluons (GeV) |
| 20 | `width_hh` | float | Decay width to Higgs pairs (GeV) |
| 21 | `total_width` | float | Total decay width (GeV) |
| 22 | `br_gaga` | float | Branching ratio to photons |
| 23 | `lam1` | float | Input λ₁ value |
| 24 | `computed_lam1` | float | Computed λ₁ from model (verification) |
| 25 | `lam2` | float | λ₂ coupling |
| 26 | `computed_lam2` | float | Computed λ₂ from model (verification) |
| 27 | `lam3` | float | λ₃ coupling |
| 28 | `lam4` | float | λ₄ coupling |
| 29 | `lam5` | float | λ₅ coupling |

### Command Examples

#### Example 1: Single-Point Verification

Scan a single m_φ and λ₁ point to verify correct execution:

```bash
OMP_NUM_THREADS=1 dihiggs/app/PhysLam1Scan 130 130 1 0.1 0.1 1 300 0.999 50 0.1 0.0 /tmp/physlam1scan_verify.csv
```

**Parameters:**
- m_φ: 130 GeV (fixed)
- λ₁: 0.1 (fixed)
- m_A: 300 GeV
- sin(β-α): 0.999
- tan(β): 50
- λ₆, λ₇: 0.1, 0.0

**Expected output:**
```
Total Attempts: 1
TRIPLE_OK_POINTS 1
```

**Output file:** 1 CSV row containing physics constraints and decay widths for the single point.

#### Example 2: Multi-Point Parameter Scan

Perform a 3×5 parameter scan (15 total points):

```bash
OMP_NUM_THREADS=1 dihiggs/app/PhysLam1Scan 130 135 3 0.1 0.5 5 300 0.999 50 0.1 0.0 /tmp/physlam1scan_scan.csv
```

**Parameters:**
- m_φ: scan from 130 to 135 GeV in 3 points (130, 132.5, 135 GeV)
- λ₁: scan from 0.1 to 0.5 in 5 points (0.1, 0.2, 0.3, 0.4, 0.5)
- m_A: 300 GeV
- sin(β-α): 0.999
- tan(β): 50
- λ₆, λ₇: 0.1, 0.0

**Expected output:**
```
Total Attempts: 15
TRIPLE_OK_POINTS 15
```

**Output file:** 15 CSV rows, one for each parameter combination, all passing physics constraints.

## Lam1 Adaptive Backend Usage

The lam1 adaptive backend is configured in:

- `autoresearch/configs/dihiggs_explorers_lam1.json`
- `dihiggs/app/adaptive_explorer_lam1.py`

This backend delegates each proposal to `dihiggs/app/PhysLam1Scan` and writes:

- `events.jsonl` in the selected output directory
- `checkpoints/<arm-id>/adaptive_state.json`
- `checkpoints/<arm-id>/run_*/results.csv`

### Backend preflight check

```bash
python3 -c "import json,sys; sys.path.insert(0,'autoresearch'); from harness.dihiggs_preflight import run_all_preflight_checks; cfg=json.load(open('autoresearch/configs/dihiggs_explorers_lam1.json')); print(run_all_preflight_checks(cfg))"
```

### One-round adaptive smoke run

```bash
python3 -c "import json,sys; sys.path.insert(0,'autoresearch'); from harness.dihiggs_runner import DiHiggsRunner; cfg=json.load(open('autoresearch/configs/dihiggs_explorers_lam1.json')); out='/tmp/physlam1_backend_readme_smoke'; r=DiHiggsRunner(config=cfg,outdir=out); print(r.run_single_round('adaptive-smoke'))"
```

Use `OMP_NUM_THREADS=1` for deterministic and thread-safe execution of this backend.
