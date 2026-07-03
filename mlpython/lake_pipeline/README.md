# 📚 Parquet Analysis Pipeline - Documentation Index

**Location:** `/home/fabi/wt_dihiggs_exploratory/mlpython/2603/`  
**Data:** `temp_subspace_l6_tb_high.parquet` (502,250 rows, 29 columns)

---

## 🎯 Start Here

### For The Impatient 
```bash
cd /home/fabi/wt_dihiggs_exploratory/mlpython/2603
source /home/fabi/higgs_env_py312/bin/activate
./run_pipeline.sh check        # Verify setup
./run_pipeline.sh all          # Run everything
```

### For The Detail-Oriented
1. Read **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** (5 min) - All commands at a glance
2. Read **[QUICK_RUN_GUIDE.md](QUICK_RUN_GUIDE.md)** (10 min) - Full workflow explanation  
3. Read **[PARQUET_SCHEMA.md](PARQUET_SCHEMA.md)** (5 min) - Understand your data

---

## 📖 Documentation Files

### 1. **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** ⭐ START HERE
**Purpose:** Fast lookup of all commands  
**Length:** ~5 minutes  
**Contains:**
- One-liners for every task
- Common issues & fixes
- Output interpretation
- Data inspection snippets

**Best for:** Running specific analyses, troubleshooting

---

### 2. **[QUICK_RUN_GUIDE.md](QUICK_RUN_GUIDE.md)** 📋 MAIN GUIDE
**Purpose:** Comprehensive workflow explanation  
**Length:** ~10-15 minutes  
**Contains:**
- Dataset overview
- Prerequisites & setup
- Sequential workflow (4 phases)
- Individual script documentation
- Troubleshooting guide
- Recommended workflow summary

**Best for:** Understanding the full pipeline, first-time setup

---

### 3. **[PARQUET_SCHEMA.md](PARQUET_SCHEMA.md)** 🔍 DATA REFERENCE
**Purpose:** Data schema and column documentation  
**Length:** ~5 minutes  
**Contains:**
- Complete column listing (all 29 columns)
- Column purposes and types
- Grouping by physical meaning
- Query examples
- Data characteristics
- Notes for data users

**Best for:** Understanding data structure, writing custom queries

---

## 🚀 Quick Start Workflows

### Scenario 1: "Just run everything!"
```bash
./run_pipeline.sh all
```
**Time:** ~30-40 minutes  
**Output:** All 3 plot types (ctau, m_phi, comparison)

---

### Scenario 2: "I want just ctau plots"
```bash
./run_pipeline.sh ctau
```
**Time:** ~10 minutes  
**Output:** `ctau_plots/` directory with PNG/PDF files

---

### Scenario 3: "Let me test first"
```bash
python ctau_vs_br_multislice_patched.py \
    --input ./temp_subspace_l6_tb_high.parquet \
    --color-by lambda6 \
    --max-combinations 3
```
**Time:** ~1-2 minutes  
**Use:** Verify setup before full run

---

### Scenario 4: "I want publication-quality plots"
```bash
./run_pipeline.sh mphi
```
**Time:** ~10 minutes  
**Output:** `paper_plots_mphi_br/` directory with high-res PNG/PDF

---

### Scenario 5: "Check my data filtering"
```bash
./run_pipeline.sh compare
```
**Time:** ~5 minutes  
**Output:** `subspace_comparisons_logs/` with comparison histograms

---

## 🛠️ Available Tools

### The Runner Script
**File:** `run_pipeline.sh` ✓ Executable  
**Usage:** `./run_pipeline.sh [STAGE]`

| Command | Effect |
|---------|--------|
| `./run_pipeline.sh help` | Show all options |
| `./run_pipeline.sh check` | Verify environment (fast) |
| `./run_pipeline.sh ctau` | Generate ctau plots |
| `./run_pipeline.sh mphi` | Generate m_phi plots |
| `./run_pipeline.sh compare` | Generate comparison plots |
| `./run_pipeline.sh all` | Run all stages |

---

### Direct Python Scripts

| Script | Purpose | Time |
|--------|---------|------|
| `ctau_vs_br_multislice_patched.py` | Lifetime vs decay plots | 5-15 min |
| `paper_like_mphi_vs_br_patched.py` | Mass vs BR plots (paper-quality) | 5-15 min |
| `subspace_comparator.py` | Compare data distributions | 5-10 min |
| `EDA_subSpace.ipynb` | Exploratory analysis (Jupyter) | 1-2 min |

---

### Python Utilities

| File | Purpose |
|------|---------|
| `test_nulls.py` | Check for null values in flags |
| `lake_explorer.py` | Explore data lake structure |
| `subset_unique_values_to_json.py` | Export unique parameter values |

---

## 📊 Output Directories

After running the pipeline:

```
2603/
├── ctau_plots/                    # Lifetime vs branching ratio
│   ├── ctau_vs_br_gaga_*.png
│   ├── ctau_vs_br_gaga_*.pdf
│   └── ctau_summary.json
│
├── paper_plots_mphi_br/           # Mass vs branching ratio (pub-quality)
│   ├── mphi_vs_*.png
│   ├── mphi_vs_*.pdf
│   └── summary.json
│
└── subspace_comparisons_logs/     # Dataset comparisons
    ├── compare_*.png
    └── execution_log.txt
```

---

## ✅ Checklist: Getting Started

- [ ] Read **QUICK_REFERENCE.md** (5 min)
- [ ] Run `./run_pipeline.sh check` (1 min)  
- [ ] Run `python ctau_vs_br_multislice_patched.py ... --max-combinations 2` (2 min)
- [ ] Review outputs in `ctau_plots/`
- [ ] Run `./run_pipeline.sh all` (30 min)
- [ ] Review all outputs

---

## 🔧 Environment

**Python:** 3.12  
**Environment:** `/home/fabi/higgs_env_py312/`

**Key Dependencies:**
- polars ≥ 0.19
- pandas
- matplotlib
- numpy
- pyarrow

**Verify:**
```bash
./run_pipeline.sh check
```

---

## 📞 Quick Help

### "Where do I start?"
→ Read **QUICK_REFERENCE.md**, then run `./run_pipeline.sh check`

### "How do I run just ctau plots?"
→ Read **QUICK_REFERENCE.md** section "Individual Phases", then:
```bash
./run_pipeline.sh ctau
```

### "What columns are in my data?"
→ Read **PARQUET_SCHEMA.md**, or run:
```bash
python -c "import polars as pl; lf = pl.scan_parquet('./temp_subspace_l6_tb_high.parquet'); print('\n'.join(lf.collect_schema().names()))"
```

### "I got an error about missing columns"
→ See **PARQUET_SCHEMA.md** "Complete Column Listing" or run:
```bash
./run_pipeline.sh check
```

### "How long will this take?"
→ See **QUICK_RUN_GUIDE.md** "Quick Start (Sequential Workflow)" or **QUICK_REFERENCE.md** "Time Estimate"

### "Can I run things in parallel?"
→ See **QUICK_REFERENCE.md** "Batch Processing (Background)"

### "What if a script fails?"
→ See **QUICK_RUN_GUIDE.md** "Troubleshooting" section

---

## 🎯 Common Workflows

### Complete Analysis (All Plots)
```bash
./run_pipeline.sh all
# Total time: 30-40 minutes
# Generates: 3 directories with PNG/PDF plots + metadata
```

### Fast Testing
```bash
python ctau_vs_br_multislice_patched.py \
    --input ./temp_subspace_l6_tb_high.parquet \
    --color-by lambda6 \
    --max-combinations 3
# Total time: 2-3 minutes
# Tests: Basic setup and data flow
```

### Publication Plotting
```bash
./run_pipeline.sh mphi
# Total time: 10-15 minutes
# Generates: `paper_plots_mphi_br/` with high-res PNG/PDF
```

---

## 📈 Next Steps After Analysis

1. **Review Plots:** Open PNG files in your image viewer
2. **Check Metadata:** Read `.json` summary files for statistics
3. **Customize:** Edit scripts' parameters for different analyses
4. **Export Results:** Copy plots to your paper/presentation

Example metadata file:
```bash
cat ctau_plots/ctau_summary.json | python -m json.tool
```

---

## 📝 File Organization

```
2603/
├── Documentation/                 # ← YOU ARE HERE
│   ├── README.md                 # This file
│   ├── QUICK_REFERENCE.md        # Fast lookup
│   ├── QUICK_RUN_GUIDE.md        # Full guide
│   └── PARQUET_SCHEMA.md         # Data reference
│
├── Scripts/
│   ├── run_pipeline.sh           # Main runner ✓
│   ├── ctau_vs_br_multislice_patched.py
│   ├── paper_like_mphi_vs_br_patched.py
│   ├── subspace_comparator.py
│   └── [utilities]
│
├── Data/
│   └── temp_subspace_l6_tb_high.parquet  # Your data
│
├── Notebooks/
│   ├── EDA_subSpace.ipynb
│   └── quick_view.ipynb
│
└── Outputs/
    ├── ctau_plots/              # Generated  
    ├── paper_plots_mphi_br/     # Generated
    └── subspace_comparisons_logs/  # Generated
```

---

## ⚡ Pro Tips

1. **Run in background:** `nohup ./run_pipeline.sh all > pipeline.log 2>&1 &`
2. **Monitor progress:** `tail -f pipeline.log`
3. **Quick column check:** `python -c "..."`  (See PARQUET_SCHEMA.md)
4. **Test before full run:** Always use `--max-combinations` first
5. **Save disk space:** Delete old plot directories before re-running

---

## 🆘 Troubleshooting Quick Links

| Problem | Solution |
|---------|----------|
| Import errors | See QUICK_RUN_GUIDE.md → Troubleshooting |
| Memory error | Use `--max-combinations 10` |
| No plots generated | Check logs (printed to stdout), see QUICK_RUN_GUIDE.md |
| "Column not found" | Check PARQUET_SCHEMA.md or run `./run_pipeline.sh check` |
| Notebook not opening | Install Jupyter: `pip install jupyter` |

---

## 📞 Document Index TL;DR

| Need | Read | Time |
|------|------|------|
| Quick commands | QUICK_REFERENCE.md | 5 min |
| Full walkthrough | QUICK_RUN_GUIDE.md | 15 min |
| Understand data | PARQUET_SCHEMA.md | 5 min |
| Run everything | `./run_pipeline.sh all` | 40 min |

---

**Last Updated:** April 7, 2026  
**Data Version:** temp_subspace_l6_tb_high.parquet (502,250 rows)  
**Environment:** Python 3.12 + higgs_env_py312

