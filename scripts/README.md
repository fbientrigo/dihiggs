# scripts/: study and campaign scripts

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
