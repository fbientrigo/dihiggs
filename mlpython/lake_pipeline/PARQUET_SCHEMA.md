# 📋 Parquet Schema Reference
**File:** `temp_subspace_l6_tb_high.parquet`  
**Total Rows:** 502,250  
**Total Columns:** 29  
**Data Type:** Polars LazyFrame

---

## Complete Column Listing

| # | Column Name | Type | Description |
|---|---|---|---|
| 1 | `m_phi` | Float64 | Mass of scalar φ [GeV] |
| 2 | `mA` | Float64 | Mass of pseudoscalar A [GeV] |
| 3 | `alpha` | Float64 | Mixing angle α |
| 4 | `beta` | Float64 | Mixing angle β |
| 5 | `lambda6` | Float64 | Z₂-breaking quartic coupling λ₆ |
| 6 | `lambda7` | Float64 | Z₂-breaking quartic coupling λ₇ |
| 7 | `m12` | Float64 | Soft-breaking mass parameter m₁₂² |
| 8 | `sin_ba` | Float64 | sin(β - α) |
| 9 | `tan_beta` | Float64 | tan(β) ratio of VEVs |
| 10 | `positivity_ok` | Float64 | Flag: 1 if bounded-from-below satisfied |
| 11 | `unitarity_ok` | Float64 | Flag: 1 if unitarity constraints satisfied |
| 12 | `perturbativity_ok` | Float64 | Flag: 1 if perturbativity satisfied |
| 13 | `width_bb` | Float64 | Partial width: φ → b b̄ [GeV] |
| 14 | `width_tautau` | Float64 | Partial width: φ → τ⁺τ⁻ [GeV] |
| 15 | `width_WW` | Float64 | Partial width: φ → WW [GeV] |
| 16 | `width_ZZ` | Float64 | Partial width: φ → ZZ [GeV] |
| 17 | `width_gaga` | Float64 | Partial width: φ → γγ [GeV] |
| 18 | `width_Zga` | Float64 | Partial width: φ → Zγ [GeV] |
| 19 | `width_gg` | Float64 | Partial width: φ → gg [GeV] |
| 20 | `width_hh` | Float64 | Partial width: φ → hh [GeV] |
| 21 | `total_width` | Float64 | Total decay width Γ_total [GeV] |
| 22 | `br_gaga` | Float64 | Branching ratio: BR(φ → γγ) |
| 23 | `lam1` | Float64 | Quartic coupling λ₁ |
| 24 | `computed_lam1` | Float64 | λ₁ computed from other parameters |
| 25 | `lam2` | Float64 | Quartic coupling λ₂ |
| 26 | `computed_lam2` | Float64 | λ₂ computed from other parameters |
| 27 | `lam3` | Float64 | Quartic coupling λ₃ |
| 28 | `lam4` | Float64 | Quartic coupling λ₄ |
| 29 | `lam5` | Float64 | Quartic coupling λ₅ |

---

## Column Grouping by Purpose

### Physical Parameters (Masses & Angles)
```
mA, m_phi, alpha, beta, tan_beta, sin_ba
```

### Soft-Breaking Terms
```
m12, lambda6, lambda7
```

### Theoretical Validity Flags ✓
```
positivity_ok, unitarity_ok, perturbativity_ok
```
**Note:** All rows have these flags = 1 (already filtered)

### Decay Channels (Partial Widths)
```
width_bb, width_tautau, width_WW, width_ZZ, width_gaga, width_Zga, width_gg, width_hh
```

### Total Decay & Branching Ratios
```
total_width, br_gaga
```

### Quartic Couplings
```
lam1, computed_lam1, lam2, computed_lam2, lam3, lam4, lam5
```

---

## Quick Querying Examples

### 1. Get statistics on mass parameters
```python
import polars as pl

lf = pl.scan_parquet('./temp_subspace_l6_tb_high.parquet')
stats = lf.select(['m_phi', 'mA', 'tan_beta']).describe().collect()
print(stats)
```

### 2. Find points with large γγ branching ratio
```python
lf = pl.scan_parquet('./temp_subspace_l6_tb_high.parquet')
result = (
    lf.filter(pl.col('br_gaga') > 0.1)
    .select(['m_phi', 'br_gaga', 'total_width'])
    .collect()
)
print(f"Found {len(result)} points with BR(γγ) > 0.1")
```

### 3. Compute lifetime (c·τ) in mm
```python
HBAR = 6.582119569e-25  # GeV·s
C_SPEED = 2.99792458e11  # mm/s

lf = pl.scan_parquet('./temp_subspace_l6_tb_high.parquet')
result = (
    lf.with_columns(
        pl.when(pl.col('total_width') > 0)
        .then((HBAR * C_SPEED) / pl.col('total_width').cast(pl.Float64))
        .otherwise(None)
        .alias('ctau_mm')
    )
    .select(['m_phi', 'total_width', 'ctau_mm'])
    .collect()
)
```

### 4. Unique parameter combinations
```python
lf = pl.scan_parquet('./temp_subspace_l6_tb_high.parquet')
unique_pairs = (
    lf.select(['mA', 'lambda6'])
    .unique()
    .sort(['mA', 'lambda6'])
    .collect()
)
print(f"Total unique (mA, λ₆) pairs: {len(unique_pairs)}")
```

---

## Data Characteristics

### Coverage
- **mA range:** Depends on campaign (typically 100-1000+ GeV)
- **m_phi range:** Depends on campaign (typically 50-1000+ GeV)
- **tan_beta range:** High values (typically 1-50000)
- **lambda6 range:** Log-spaced (typically 10⁻⁵ to 1)

### Quality
- **Completeness:** No missing values in key columns
- **Validity:** All rows satisfy `triple_ok` filters (positivity, unitarity, perturbativity)
- **Precision:** Float64 (64-bit IEEE 754)

### Theoretical Constraints Applied
✓ Bounded from below (Positivity)  
✓ Unitarity of scattering matrix  
✓ Perturbative regime maintained  

**NOT Applied (see scripts for optional application):**
- HiggsBounds constraints
- HiggsSignals precision measurements
- Direct experimental constraints

---

## Notes for Data Users

1. **Partial Widths Are Conservative**: The `width_*` columns represent tree-level calculations. Loop corrections not included.

2. **Branching Ratios**: Only `br_gaga` is provided. Compute others as:
   ```
   BR(φ → X) = width_X / total_width
   ```

3. **Computed vs Direct Lambda**: 
   - Use `lam1`, `lam2` as computed from the 2HDMC engine
   - `computed_lam1`, `computed_lam2` are cross-checks

4. **Lifetime Calculations**: Must account for proper unit conversion (GeV ↔ mm/s):
   ```
   c·τ [mm] = (ℏc) / Γ = (6.582e-25 GeV·s × 2.998e11 mm/s) / Γ[GeV]
   ```

5. **For Plotting**:
   - Use **logarithmic scales** for: `tan_beta`, `lambda6`, `total_width`, `ctau`
   - Use **linear scales** for: masses, branching ratios

---

## File Info

- **Format:** Apache Parquet (columnar storage)
- **Engine:** PyArrow (Python) / Polars (native)
- **Compression:** Snappy (lossless, fast)
- **Size on Disk:** ~200-500 MB (depending on subset)
- **Memory When Loaded:** ~1-2 GB (full dataset)

