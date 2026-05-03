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
