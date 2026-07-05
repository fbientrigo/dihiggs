# 📋 Quick Reference Card
**Parquet Analysis Pipeline** | `<repo>/mlpython/lake_pipeline/`

---

## 🚀 Quick Commands

### Full Pipeline (All Phases)
```bash
cd <repo>/mlpython/lake_pipeline
source "$HIGGS_ENV_ACTIVATE"  # optional: your Python 3.12 venv
./run_pipeline.sh all
```

### Individual Phases
```bash
./run_pipeline.sh check        # Verify setup
./run_pipeline.sh eda          # Open notebook (manual)
./run_pipeline.sh ctau         # Generate ctau plots (~10 min)
./run_pipeline.sh mphi         # Generate m_phi plots (~10 min)
./run_pipeline.sh compare      # Generate comparison plots (~5 min)
```

---

## 🎯 Direct Python Commands

### ctau vs BR Plots (Full)
```bash
python ctau_vs_br_multislice_patched.py \
    --input ./temp_subspace_l6_tb_high.parquet \
    --color-by both
```

### ctau vs BR Plots (Testing with 5 combinations)
```bash
python ctau_vs_br_multislice_patched.py \
    --input ./temp_subspace_l6_tb_high.parquet \
    --color-by lambda6 \
    --max-combinations 5
```

### m_phi vs BR Plots (Full)
```bash
python paper_like_mphi_vs_br_patched.py \
    --input ./temp_subspace_l6_tb_high.parquet
```

### m_phi vs BR Plots (Testing with 3 slices)
```bash
python paper_like_mphi_vs_br_patched.py \
    --input ./temp_subspace_l6_tb_high.parquet \
    --max-slices 3
```

## 📁 Output Structure

After running pipeline:
```
lake_pipeline/
├── ctau_plots/
│   ├── ctau_vs_br_gaga_*.png      # Main outputs
│   ├── ctau_vs_br_gaga_*.pdf
│   └── ctau_summary.json           # Metadata
├── paper_plots_mphi_br/
│   ├── mphi_vs_*.png
│   ├── mphi_vs_*.pdf
│   └── summary.json
└── subspace_comparisons_logs/
    ├── compare_*.png
    └── execution_log.txt
```

---

## 🔍 Data Inspection

### List Available Columns
```bash
python -c "
import polars as pl
lf = pl.scan_parquet('./temp_subspace_l6_tb_high.parquet')
print('Columns:', lf.collect_schema().names())
"
```

### Check Row Count
```bash
python -c "
import polars as pl
lf = pl.scan_parquet('./temp_subspace_l6_tb_high.parquet')
print('Rows:', lf.select(pl.len()).collect().item())
"
```

### Quick Statistics
```bash
python -c "
import polars as pl
lf = pl.scan_parquet('./temp_subspace_l6_tb_high.parquet')
stats = lf.select(['mA', 'lambda6', 'tan_beta']).describe().collect()
print(stats)
"
```

---

## ⚙️ Configuration

### Key Parameter Explained

| Parameter | Script | Effect |
|-----------|--------|--------|
| `--input` | All | Path to parquet file |
| `--color-by` | ctau | Generate plots with `lambda6` or `tan_beta` families |
| `--br-col` | mphi, ctau | Which branching ratio column (`br_gaga`, etc.) |
| `--apply-phys-filter` | All | Skip (your data is already `triple_ok` filtered) |
| `--max-combinations` | ctau | Limit tested (mA, λ₆) pairs for testing |
| `--max-slices` | mphi | Limit tested (mA, λ₆) pairs for testing |
| `--logy` | compare | Log scale for y-axis |
| `--logx-width` | compare | Log scale for x-axis on width columns |

---

## 🐛 Common Issues & Fixes

| Issue | Fix |
|-------|-----|
| `ModuleNotFoundError: polars` | `pip install polars --upgrade` |
| Out of memory | Use `--max-combinations 10` to process smaller subset |
| "Column not found" | Verify column name exists: `python -c "import polars as pl; lf = pl.scan_parquet(...); print(lf.collect_schema().names())"` |
| No plots generated | Check logs (printed to stdout with `[INFO]`, `[WARNING]`) and verify data isn't empty |
| Notebook won't open | Install Jupyter: `pip install jupyter` |

---

## 📊 Output Interpretation

### ctau_summary.json
```json
{
  "input_file": "temp_subspace_l6_tb_high.parquet",
  "valid_ctau_rows": 1234567,
  "variants": [
    {
      "mode": "lambda6",
      "plots_generated": 42,
      "discarded_combinations": 5
    }
  ]
}
```

### paper_plots_mphi_br/summary.json
```json
{
  "total_combinations_investigated": 50,
  "plots_generated": [
    {
      "mA": 200,
      "lambda6": 0.5,
      "tan_beta_curves": [1.0, 5.0, 10.0],
      "points_in_slice": 5000
    }
  ]
}
```

---

## 🎓 Workflow Recommendations

### For First Time
```bash
./run_pipeline.sh check           # Verify setup
./run_pipeline.sh ctau --max-combinations 2  # Quick test
```

### For Production Run
```bash
./run_pipeline.sh all             # Full pipeline
```

### For Debugging
```bash
# Terminal 1: Run ctau
python ctau_vs_br_multislice_patched.py --max-combinations 3

# Terminal 2: Monitor output
watch -n 1 "ls -ltrh ctau_plots/ | tail -5"
```

### For Batch Processing (Background)
```bash
# Run in background
nohup ./run_pipeline.sh all > pipeline.log 2>&1 &

# Monitor progress
tail -f pipeline.log
```

---

## 📞 Quick Help

```bash
# Show script help
./run_pipeline.sh help

# Check environment status
./run_pipeline.sh check

# View script details
head -50 ctau_vs_br_multislice_patched.py
head -50 paper_like_mphi_vs_br_patched.py
```

---

## 📝 Notes

- **Default Parquet:** `./temp_subspace_l6_tb_high.parquet` (must be in working directory)
- **Data Already Filtered:** Your parquet is pre-filtered to `triple_ok`, omit `--apply-phys-filter`
- **Memory Usage:** ~1-2GB for full dataset (Polars streaming engine optimizes this)
- **Time Estimate:** Full pipeline ~30-40 minutes
- **Output Format:** PNG (300 dpi) + PDF for publications
