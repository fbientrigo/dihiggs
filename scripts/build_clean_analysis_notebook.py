#!/usr/bin/env python3
from pathlib import Path
import nbformat as nbf


def build_notebook(path: Path):
    nb = nbf.v4.new_notebook()
    cells = []

    cells.append(nbf.v4.new_markdown_cell("""# Colab vs 2HDMC clean discrepancy analysis

Purpose: compare Colab analytical points vs 2HDMC validation outputs using stable artifacts.

This notebook is **not** a lambda1 drift analysis.
Main focus: channel-level and physical-configuration-level discrepancy patterns."""))

    cells.append(nbf.v4.new_code_cell("""from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style='whitegrid')
BASE = Path('scripts/out/refactor_colab_compare/real_run')
FIG_DIR = Path('scripts/out/refactor_colab_compare/figures')
FIG_DIR.mkdir(parents=True, exist_ok=True)

merged = pd.read_csv(BASE/'colab_vs_2hdmc_merged_comparison.csv')
width_stats = pd.read_csv(BASE/'width_rel_stats.csv')
grouped = pd.read_csv(BASE/'grouped_error_by_l1_l6_tb.csv')

print('merged rows:', len(merged))
print('set_param_ok:', int(merged['set_param_ok'].sum()) if 'set_param_ok' in merged.columns else 'n/a')
print('triple_ok:', int(merged['triple_ok'].sum()) if 'triple_ok' in merged.columns else 'n/a')
print('warning_flag:', int(merged['warning_flag'].sum()) if 'warning_flag' in merged.columns else 'n/a')

required = ['point_id','lambda1','lambda6','tan_beta','rel_err_width_gg','rel_err_width_gaga','rel_err_width_total']
missing = [c for c in required if c not in merged.columns]
print('missing required columns:', missing)
print('warning_flag caveat: likely strict lambda1 round-trip threshold; not automatic physics veto.')

df = merged[merged['triple_ok'] == True].copy() if 'triple_ok' in merged.columns else merged.copy()"""))

    cells.append(nbf.v4.new_code_cell("""# Global channel discrepancy overview
ordered = width_stats.sort_values('median_rel', ascending=False)

plt.figure(figsize=(8,4))
sns.barplot(data=ordered, x='channel', y='median_rel', color='#4C78A8')
plt.title('Median relative width error by channel (triple_ok)')
plt.tight_layout()
plt.savefig(FIG_DIR/'channel_median_rel_errors.png', dpi=160)

plt.figure(figsize=(8,4))
sns.barplot(data=ordered, x='channel', y='max_rel', color='#F58518')
plt.title('Max relative width error by channel (triple_ok)')
plt.tight_layout()
plt.savefig(FIG_DIR/'channel_max_rel_errors.png', dpi=160)

ordered[['channel','mean_rel','median_rel','max_rel']]"""))

    cells.append(nbf.v4.new_code_cell("""# Physical configuration summary
show = grouped.sort_values('median_rel_total', ascending=False).copy()
show[['lambda1','lambda6','tan_beta','n','median_rel_gg','median_rel_gaga','median_rel_total']]"""))

    cells.append(nbf.v4.new_code_cell("""plt.figure(figsize=(7,4))
sns.barplot(data=show, x='lambda6', y='median_rel_total', hue='lambda1')
plt.title('Grouped median rel error: total width')
plt.tight_layout()
plt.savefig(FIG_DIR/'grouped_config_median_rel_total.png', dpi=160)

fig, axes = plt.subplots(1,2, figsize=(12,4), sharex=True)
sns.barplot(data=show, x='lambda6', y='median_rel_gg', hue='lambda1', ax=axes[0])
axes[0].set_title('median_rel_gg')
sns.barplot(data=show, x='lambda6', y='median_rel_gaga', hue='lambda1', ax=axes[1])
axes[1].set_title('median_rel_gaga')
for ax in axes: ax.legend(title='lambda1')
plt.tight_layout()
plt.savefig(FIG_DIR/'grouped_config_median_rel_gg_gaga.png', dpi=160)"""))

    cells.append(nbf.v4.new_code_cell("""# Channel contribution and control channels
long_cols = ['rel_err_width_bb','rel_err_width_cc','rel_err_width_tautau','rel_err_width_gg','rel_err_width_gaga','rel_err_width_Zga','rel_err_width_total']
avail = [c for c in long_cols if c in df.columns]
long = df[avail].melt(var_name='channel', value_name='rel_err')
long['channel'] = long['channel'].str.replace('rel_err_width_','', regex=False)

plt.figure(figsize=(10,4))
sns.boxplot(data=long, x='channel', y='rel_err')
plt.ylim(0, np.nanpercentile(long['rel_err'], 99.5))
plt.title('Relative error distributions by channel (triple_ok)')
plt.tight_layout()

print('tautau median:', float(df['rel_err_width_tautau'].median()) if 'rel_err_width_tautau' in df.columns else 'n/a')
print('gg median:', float(df['rel_err_width_gg'].median()) if 'rel_err_width_gg' in df.columns else 'n/a')"""))

    cells.append(nbf.v4.new_code_cell("""# mH dependence
plt.figure(figsize=(8,4))
sns.scatterplot(data=df, x='mH', y='rel_err_width_gg', hue='lambda6', style='lambda1', alpha=0.5, s=20)
plt.title('mH vs rel_err_width_gg')
plt.tight_layout()
plt.savefig(FIG_DIR/'mH_vs_rel_err_gg.png', dpi=160)

plt.figure(figsize=(8,4))
sns.scatterplot(data=df, x='mH', y='rel_err_width_gaga', hue='lambda6', style='lambda1', alpha=0.5, s=20)
plt.title('mH vs rel_err_width_gaga')
plt.tight_layout()
plt.savefig(FIG_DIR/'mH_vs_rel_err_gaga.png', dpi=160)"""))

    cells.append(nbf.v4.new_code_cell("""# Parallel coordinates / hyperparameter population view
from pandas.plotting import parallel_coordinates
cols = ['mH','lambda1','lambda6','tan_beta','rel_err_width_gg','rel_err_width_gaga','rel_err_width_total']
sub = df[cols].copy().sample(min(200, len(df)), random_state=42)
sub['color_bin'] = pd.qcut(sub['rel_err_width_total'], q=4, labels=['q1','q2','q3','q4'])

scaled = sub.copy()
for c in cols:
    v = scaled[c].astype(float)
    denom = (v.max() - v.min())
    scaled[c] = 0.0 if denom == 0 else (v - v.min()) / denom

plt.figure(figsize=(12,5))
parallel_coordinates(scaled[['color_bin']+cols], 'color_bin', colormap='viridis', alpha=0.2)
plt.title('Parallel coordinates (scaled)')
plt.tight_layout()
plt.savefig(FIG_DIR/'parallel_hyperparams_errors.png', dpi=160)"""))

    cells.append(nbf.v4.new_code_cell("""# Top problematic points (compact)
top_cols = ['point_id','mH','lambda1','lambda6','tan_beta','rel_err_width_gg','rel_err_width_gaga','rel_err_width_Zga','rel_err_width_total','width_gg','width_gg_colab','total_width','total_width_colab']
top_cols = [c for c in top_cols if c in df.columns]
rank = df.sort_values('rel_err_width_total', ascending=False)[top_cols].head(20)
rank"""))

    cells.append(nbf.v4.new_markdown_cell("""## Final interpretation
- gg discrepancy is systematic and dominant.
- tautau behaves as a stable control channel.
- gaga and total-width discrepancies are configuration/lambda6-sensitive.
- No lambda1-drift claim is made here.
- warning_flag should be audited separately before any veto use."""))

    nb.cells = cells
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(nbf.writes(nb), encoding='utf-8')


if __name__ == '__main__':
    build_notebook(Path('scripts/analyze_colab_vs_2hdmc_clean.ipynb'))
    print('Wrote scripts/analyze_colab_vs_2hdmc_clean.ipynb')
