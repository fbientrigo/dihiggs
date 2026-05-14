# scripts/: active Colab-vs-2HDMC analyzer

Scope
- `scripts/` is the active Colab-vs-2HDMC analyzer directory.
- This workflow is not lambda drift centered.
- `warning_flag` is an audit signal, not a physics veto.

Environment and output policy
- Activate: `source ~/higgs_env_py312/bin/activate`
- Run from repo root: `/home/fabi/wt_dihiggs_exploratory`
- Output root: `scripts/out/`
- Preserved and protected:
  - `scripts/out/current_baseline/` (unchanged)
  - `scripts/out/refactor_colab_compare/` (unchanged)

Active files in scripts/
- `validate_colab.py`
- `show_width_rel_stats.py`
- `group_width_errors_by_config.py`
- `build_clean_analysis_notebook.py`
- `analyze_colab_vs_2hdmc_clean.ipynb`
- `README.md`
- `__init__.py`
- `out/`

Canonical commands
```bash
source ~/higgs_env_py312/bin/activate
cd /home/fabi/wt_dihiggs_exploratory

python3 scripts/validate_colab.py \
  --out-csv scripts/out/refactor_colab_compare/real_run/colab_vs_2hdmc_merged_comparison.csv \
  --top-n 20

python3 scripts/show_width_rel_stats.py \
  scripts/out/refactor_colab_compare/real_run/colab_vs_2hdmc_merged_comparison.csv \
  --sort-by mean_rel \
  --max-point \
  --max-point-cols point_id,mH,lambda1,lambda6,tan_beta,rel_err_width_gg,rel_err_width_gaga,rel_err_width_total \
  --out-csv scripts/out/refactor_colab_compare/real_run/width_rel_stats.csv

python3 scripts/group_width_errors_by_config.py \
  scripts/out/refactor_colab_compare/real_run/colab_vs_2hdmc_merged_comparison.csv \
  --group-cols lambda1,lambda6,tan_beta \
  --out-csv scripts/out/refactor_colab_compare/real_run/grouped_error_by_l1_l6_tb.csv
```

Where files moved
- Maintenance script moved to:
  - `tools/recompute_audit/audit_recomputed_campaigns.py`
- Benchmark script moved to:
  - `tools/lake_benchmark/benchmark_lam1_from_official_lake.py`
- Legacy drift/lambda1 scripts and notebook moved to:
  - `legacy/colab_compare/stat_drift_helpers.py`
  - `legacy/colab_compare/export_lambda1_compact_report.py`
  - `legacy/colab_compare/analyze_colab_vs_2hdmc.ipynb`
- Legacy drift/lambda1 artifacts moved to:
  - `legacy/colab_compare/artifacts/`

Archive policy
- Loose artifacts and accidental nested output were archived under:
  - `scripts/out/archive/loose_artifacts_before_cleanup/`
  - `scripts/out/archive/accidental_nested_scripts_out/`

Notes
- No formulas or numerical definitions were changed in cleanup.
- `validate_colab.py` numerical definitions were not modified.

Fixed-lambda1 Christopher-style control campaign
- Purpose: reproduce a controlled Christopher-style setup with fixed lambda1 and fixed lambda6 per curve, sweeping only m_phi for fixed tan_beta, mA(=mHp in scanner), lambda7, and sin_ba.
- This is a control campaign, not an adaptive search.
- Uses existing `dihiggs/app/orchestrate_scans.py` + `PhysScanWithFixings`; no new scan framework.

Commands
```bash
source ~/higgs_env_py312/bin/activate
cd /home/fabi/wt_dihiggs_exploratory

# Dry-run (prints six orchestrator commands)
python3 scripts/run_christopher_fixed_lam1_campaign.py --dry-run

# Smoke run (tiny)
python3 scripts/run_christopher_fixed_lam1_campaign.py \
  --campaign christopher_fixed_lam1_2026apr_smoke \
  --n-mphi 5 \
  --lambda1-values 0.4 \
  --lambda6-values 0.01 \
  --force

# Full run
python3 scripts/run_christopher_fixed_lam1_campaign.py \
  --campaign christopher_fixed_lam1_2026apr \
  --threads 8 \
  --force

# Summary
python3 scripts/summarize_christopher_fixed_lam1_campaign.py \
  --campaign christopher_fixed_lam1_2026apr
```

Outputs
- Wrapper outputs:
  - `scripts/out/<campaign>/wrapper_manifest.json`
  - `scripts/out/<campaign>/wrapper.log`
- Lake campaign outputs:
  - `/mnt/c/Users/Asus/cern_db/dihiggs_lake/campaign=<campaign>/...`
  - includes `run_manifest.json`, `orchestrator.log`, `task_summary.jsonl`, `tb_*/scan_tb_*.csv`, `tb_*/scan_meta.json`
- Summary outputs:
  - `scripts/out/<campaign>/summary.csv`
  - `scripts/out/<campaign>/summary.md`

Comparing Christopher Colab, recomputed Colab, and fixed-lambda1 campaign
- Purpose:
  - Compare physical masks and channel curves between:
    1) original Christopher/Colab analytical points (as represented in the merged artifact),
    2) recomputed/validated Colab-vs-2HDMC points,
    3) fixed-lambda1 control campaign curves.
- Inputs:
  - merged artifact:
    - `scripts/out/refactor_colab_compare/real_run/colab_vs_2hdmc_merged_comparison.csv`
    - fallback: `scripts/out/current_baseline/colab_vs_2hdmc_merged_comparison.csv`
  - fixed campaign root:
    - `/mnt/c/Users/Asus/cern_db/dihiggs_lake/campaign=christopher_fixed_lam1_2026apr`
  - fixed campaign summary:
    - `scripts/out/christopher_fixed_lam1_2026apr/summary.csv`
- Command:
```bash
source ~/higgs_env_py312/bin/activate
cd /home/fabi/wt_dihiggs_exploratory
python3 scripts/compare_christopher_fixed_lam1.py --make-plots
```
- Outputs:
  - `scripts/out/christopher_fixed_lam1_2026apr/comparison/comparison_summary_by_group.csv`
  - `scripts/out/christopher_fixed_lam1_2026apr/comparison/comparison_summary_by_group.md`
  - `scripts/out/christopher_fixed_lam1_2026apr/comparison/matched_points.csv`
  - `scripts/out/christopher_fixed_lam1_2026apr/comparison/mask_confusion_by_group.csv`
  - `scripts/out/christopher_fixed_lam1_2026apr/comparison/channel_diff_by_group.csv`
  - `scripts/out/christopher_fixed_lam1_2026apr/comparison/comparison_report.md`
  - `scripts/out/christopher_fixed_lam1_2026apr/comparison/figures/*.png`
- Interpretation caveat:
  - This compares fixed-lambda1 control curves against the Colab-vs-2HDMC merged artifact. It does not prove global lambda1 irrelevance.
