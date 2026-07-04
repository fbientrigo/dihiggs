# 🚀 Quick Run Guide: Parquet Analysis Pipeline
**Last Updated:** April 2026  
**Target Parquet:** `./temp_subspace_l6_tb_high.parquet`  
**Working Directory:** `<repo>/mlpython/lake_pipeline/`

---

## 📋 Table of Contents
1. [Dataset Overview](#dataset-overview)
2. [Prerequisites](#prerequisites)
3. [Quick Start (Sequential Workflow)](#quick-start-sequential-workflow)
4. [Individual Script Reference](#individual-script-reference)
5. [Parallel Coordinates (New)](#parallel-coordinates-new)
6. [Output Directories](#output-directories)
7. [Troubleshooting](#troubleshooting)

---

## 📊 Dataset Overview

The parquet file `temp_subspace_l6_tb_high.parquet` contains:
- **~1M+ parameter points** from campaigns: `scan_l6_tb_high` and `scan_l6_tb_xhigh`
- **Filtered by:** `triple_ok` (perturbativity × unitarity × positivity)
- **Key Columns:**
  - Masses: `mA`, `m_phi`
  - Couplings: `lam1`, `lambda6`, `tan_beta`
  - Widths: `total_width`, decay channels
  - Branching Ratios: `br_gaga` (photon pair), others
  - Flags: `positivity_ok`, `unitarity_ok`, `perturbativity_ok`

---

## ⚙️ Prerequisites

### 1. Environment Setup
```bash
# Activate your Python environment
source "$HIGGS_ENV_ACTIVATE"  # optional: your Python 3.12 venv

# Verify required packages
pip list | grep -E "polars|pandas|matplotlib|numpy|pyarrow"
```

**Required Packages:**
- `polars` ≥ 0.19
- `pandas`
- `matplotlib`
- `numpy`
- `pyarrow` (for parquet I/O)

### 2. File Paths
All commands assume you're in: `<repo>/mlpython/lake_pipeline/`

```bash
cd <repo>/mlpython/lake_pipeline
```

---

## 🔄 Quick Start (Sequential Workflow)

Run this sequence to generate all plots and analyses **in recommended order**:

### Phase 1: Exploratory Data Analysis (Notebook)
```bash
# Option A: Use Jupyter in terminal
jupyter notebook EDA_subSpace.ipynb

# Option B: View with VSCode (interactive execution)
code EDA_subSpace.ipynb
```
**Produces:** Statistical summaries, density plots, phase-space visualizations  
**Time:** ~1-2 minutes  
**Output:** In-notebook plots + `temp_subspace_l6_tb_high.parquet` (if regenerating)

---

### Phase 2: Photon Pair Decay Analysis (ctau vs BR)
```bash
python ctau_vs_br_multislice_patched.py \
    --input ./temp_subspace_l6_tb_high.parquet \
    --color-by both \
    --br-col br_gaga
```

**Parameters:**
- `--input`: Path to parquet (default: `temp_subspace.parquet`)
- `--color-by`: `lambda6` | `tan_beta` | `both` (generates separate variants)
- `--br-col`: Branching ratio column (default: `br_gaga`)
- `--apply-phys-filter`: Skip this flag (your parquet is already `triple_ok` filtered)
- `--max-combinations`: Optional safety limit (e.g., `10` for testing)

**Produces:**
- PNG/PDF plots: `ctau_plots/ctau_vs_*.{png,pdf}`
- Summary JSON: `ctau_plots/ctau_summary.json`

**Time:** 5-15 minutes  
**Memory:** ~1-2 GB (streaming Polars engine)

---

### Phase 3: Mass vs Branching Ratio (m_φ vs BR)
```bash
python paper_like_mphi_vs_br_patched.py \
    --input ./temp_subspace_l6_tb_high.parquet \
    --mphi-col m_phi \
    --br-col br_gaga
```

**Parameters:**
- `--input`: Path to parquet
- `--mphi-col`: Mass column name (default: `m_phi`)
- `--br-col`: Branching ratio column (default: `br_gaga`)
- `--apply-phys-filter`: Skip (already filtered)
- `--max-slices`: Optional limit on (mA, lambda6) combinations

**Produces:**
- PNG/PDF plots: `paper_plots_mphi_br/mphi_vs_*.{png,pdf}`
- Summary JSON: `paper_plots_mphi_br/summary.json`

**Time:** 5-15 minutes  
**Memory:** ~1-2 GB

---

### Phase 4: Parallel Coordinates (Raw + Log, color by ctau)
```bash
python parallel_coordinates_roi.py \
        --input ./temp_subspace_l6_tb_high.parquet \
        --plot-modes both \
        --axes tan_beta,lambda6,m_phi,lam1,total_width,br_gaga \
        --color-by ctau
```

**Qué hace ahora:**
- El plot global se colorea por `ctau` (la información de `ctau` va en color).
- Por defecto `ctau` no se grafica como eje en coordenadas paralelas.
- Genera dos versiones sin escala global normalizada:
    - `raw`: valores numéricos directos
    - `log`: valores en log10 para columnas definidas en `--log-cols`
- Dibuja límites numéricos por eje (mínimo y máximo visibles) para identificar regiones.

**Salida:** `parallel_plots/` + `parallel_summary.json`

---

## 📚 Individual Script Reference

### 1. `ctau_vs_br_multislice_patched.py`
**Purpose:** Analyze lifetime vs decay width for different parameter slices

| Aspect | Details |
|--------|---------|
| **Input** | Parquet with `total_width`, branching ratios, `mA`, `lambda6`, `tan_beta` |
| **Computation** | $c\tau = \frac{\hbar c}{\Gamma}$ (lifetime in mm) |
| **Variants** | Generated for each `color_mode`: color curves by λ₆ or tan β |
| **Output Dir** | `ctau_plots/` |
| **Key Outputs** | PNG/PDF plots per (mA, λ₆) slice with tan β family curves |

**Example with Testing:**
```bash
python ctau_vs_br_multislice_patched.py \
    --input ./temp_subspace_l6_tb_high.parquet \
    --color-by both \
    --max-combinations 5  # Test with only 5 (mA, λ₆) combinations first
```

---

### 2. `paper_like_mphi_vs_br_patched.py`
**Purpose:** Generate publication-quality plots of scalar mass vs branching ratio

| Aspect | Details |
|--------|---------|
| **Input** | Parquet with `m_phi`, branching ratios, `mA`, `lambda6`, `tan_beta` |
| **Style** | Paper-friendly serif font, log scales, error bands |
| **Variants** | One plot per (mA, λ₆) pair, with tan β family |
| **Output Dir** | `paper_plots_mphi_br/` |
| **Key Outputs** | High-res PNG/PDF suitable for publications |

An example of execution for the paper like
```bash
# Pipeline-compatible family mode (recomendado)
python paper_like_mphi_vs_br_patched.py --input ./temp_subspace_l6_tb_high.parquet --mphi-col m_phi --br-col br_gaga

# Fixed-cut legacy mode (solo si quieres un corte puntual)
python paper_like_mphi_vs_br_patched.py --input temp_subspace.parquet --sin-ba 1 --tan-beta FIXED_VALUE --mA FIXED_VALUE --lambda6 FIXED_VALUE --lambda7 FIXED_VALUE
```
---

## Parallel Coordinates (New)

Script: `parallel_coordinates_roi.py`

### Objetivos
- Ver flujo multidimensional en el espacio paramétrico.
- Resaltar ROI: `ctau < threshold` y `br_gaga > threshold`.
- Colorear el conjunto completo por `ctau` para evitar meter `ctau` como eje.
- Generar versiones `raw` y `log` sin escala global [0,1] visible.

### Selección modular de columnas (incluir / quitar)

Definir columnas explícitas (orden importa):
```bash
--axes tan_beta,lambda6,m_phi,lam1,total_width,br_gaga
```

Excluir columnas después de `--axes`:
```bash
--exclude-axes lam1,total_width
```

Ejemplo completo (quitar `lam1`):
```bash
python parallel_coordinates_roi.py \
    --input ./temp_subspace_l6_tb_high.parquet \
    --plot-modes both \
    --axes tan_beta,lambda6,m_phi,lam1,total_width,br_gaga \
    --exclude-axes lam1 \
    --color-by ctau
```

### Modo numérico sin normalización global

- `--plot-modes raw`: versión con valores directos.
- `--plot-modes log`: versión log10 para columnas en `--log-cols`.
- `--plot-modes both`: genera ambas.

Ejemplos:
```bash
# Solo versión raw
python parallel_coordinates_roi.py --input ./temp_subspace_l6_tb_high.parquet --plot-modes raw

# Solo versión log
python parallel_coordinates_roi.py --input ./temp_subspace_l6_tb_high.parquet --plot-modes log
```

### ROI (bajo ctau, alto BR)

Umbrales manuales:
```bash
python parallel_coordinates_roi.py \
    --input ./temp_subspace_l6_tb_high.parquet \
    --ctau-threshold 1e-3 \
    --br-threshold 1e-2
```

Umbrales dinámicos (por cuantiles, default):
- `ctau < q20`
- `br_gaga > q80`

---

## 📁 Output Directories

After running the full pipeline, your directory structure will look like:

```
lake_pipeline/
├── temp_subspace_l6_tb_high.parquet          # Your input data
├── ctau_plots/
│   ├── ctau_vs_*.png                         # Lifetime vs BR plots
│   ├── ctau_vs_*.pdf
│   └── ctau_summary.json                     # Metadata + plot inventory
├── paper_plots_mphi_br/
│   ├── mphi_vs_*.png                         # Publication plots
│   ├── mphi_vs_*.pdf
│   └── summary.json
├── subspace_comparisons_logs/
│   ├── compare_*.png                         # Comparison histograms
│   └── execution_log.txt
├── parallel_plots/
│   ├── parallel_full_colorby_ctau_raw.png
│   ├── parallel_full_colorby_ctau_log.png
│   ├── parallel_roi_raw.png
│   ├── parallel_roi_log.png
│   └── parallel_summary.json
└── plots/                                     # (Pre-existing directory)
```

---

## 🔧 Troubleshooting

### Issue: `ModuleNotFoundError: No module named 'polars'`
**Solution:**
```bash
pip install polars pyarrow --upgrade
```

### Issue: Memory Error ("MemoryError" or "out of memory")
**Solution:** Use `--max-combinations` or `--max-slices` to process subsets:
```bash
python ctau_vs_br_multislice_patched.py \
    --input ./temp_subspace_l6_tb_high.parquet \
    --color-by lambda6 \
    --max-combinations 20
```

### Issue: "Could not resolve column 'br_gaga'"
**Solution:** Check available columns:
```bash
python -c "
import polars as pl
lf = pl.scan_parquet('./temp_subspace_l6_tb_high.parquet')
print(lf.collect_schema().names())
"
```

### Issue: Plots not generating (empty output directories)
**Solution:** 
1. Verify the parquet isn't empty:
   ```bash
   python -c "
   import polars as pl
   lf = pl.scan_parquet('./temp_subspace_l6_tb_high.parquet')
   print('Rows:', lf.select(pl.len()).collect().item())
   print('Schema:', lf.collect_schema().names()[:5], '...')
   "
   ```

2. Check script logs (they print to stdout with `[INFO]`, `[WARNING]` prefixes)

3. Try with reduced data:
   ```bash
   python ctau_vs_br_multislice_patched.py \
       --input ./temp_subspace_l6_tb_high.parquet \
       --color-by lambda6 \
       --max-combinations 3
   ```

---

## 📖 Recommended Workflow Summary

```
┌─────────────────────────────────────────────────────────────┐
│ 1. Prepare Environment                  (~2 min)           │
│    source "$HIGGS_ENV_ACTIVATE"  # optional: your Python 3.12 venv           │
│    cd <repo>/mlpython/lake_pipeline       │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 2. Run EDA Notebook                     (~2 min)            │
│    jupyter notebook EDA_subSpace.ipynb                       │
│    → Review data statistics & phase-space plots             │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 3. Generate ctau vs BR Plots            (~10 min)           │
│    python ctau_vs_br_multislice_patched.py --color-by both  │
│    → Outputs: ctau_plots/                                   │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│ 4. Generate m_φ vs BR Plots             (~10 min)           │
│    python paper_like_mphi_vs_br_patched.py                  │
│    → Outputs: paper_plots_mphi_br/                          │
└─────────────────────────────────────────────────────────────┘
                            ↓
               ↓
┌─────────────────────────────────────────────────────────────┐
│ 6. Parallel Coordinates (New)         (~5-12 min)          │
│    python parallel_coordinates_roi.py --plot-modes both     │
│    → Outputs: parallel_plots/ (raw + log + ROI)             │
└─────────────────────────────────────────────────────────────┘
                            ↓
                    ✅ COMPLETE!
     Total time: ~35-50 minutes
         All results in timestamped directories
```

---

## 💡 Quick Tips

- **Test first:** Use `--max-combinations 3` to verify setup before full runs
- **Monitor progress:** All scripts print `[INFO]` / `[WARNING]` messages to stdout
- **Resume workflows:** If a script fails, fix the issue and re-run (outputs are idempotent)
- **Batch processing:** Run Phase 2 & 3 simultaneously in separate terminals for 50% speedup
- **Check outputs:** Look for `.json` summary files in each output directory for statistics

---

## 📞 Support

For detailed parameter documentation, see individual script headers:
```bash
head -30 ctau_vs_br_multislice_patched.py
head -30 paper_like_mphi_vs_br_patched.py
```

