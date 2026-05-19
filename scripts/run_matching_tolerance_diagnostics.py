from pathlib import Path
import pandas as pd
import numpy as np
import importlib.util

root = Path('/home/fabi/wt_dihiggs_exploratory')
script_path = root / 'scripts/compare_christopher_fixed_lam1.py'
spec = importlib.util.spec_from_file_location('cmp_mod', script_path)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

merged_path = mod._pick_merged_path(mod.DEFAULT_MERGED_CSV)
merged_df = mod._normalize_merged(pd.read_csv(merged_path))
fixed_df = mod._load_fixed_campaign(mod.DEFAULT_CAMPAIGN_ROOT)

tols = [1e-6, 1e-5, 1e-4, 1e-3]
param_tol = 1e-8

c = merged_df.copy().reset_index(drop=True)
f = fixed_df.copy().reset_index(drop=True)
c['_key'] = mod._group_key(c, param_tol)
f['_key'] = mod._group_key(f, param_tol)

rows = []

def mask_counts(m):
    if m.empty:
        return (0, 0, 0, 0)
    return (int(m['both_ok'].sum()), int(m['only_colab_ok'].sum()), int(m['only_fixed_ok'].sum()), int(m['both_fail'].sum()))

for tol in tols:
    m = mod.match_points(merged_df, fixed_df, mphi_tol=tol, param_tol=param_tol)
    n_colab = len(merged_df)
    n_fixed = len(fixed_df)
    n_matched = len(m)
    unmatched_colab = n_colab - n_matched
    unmatched_fixed = n_fixed - n_matched

    amb_colab = 0
    amb_fixed = 0
    inter_keys = sorted(set(c['_key']).intersection(set(f['_key'])))
    for k in inter_keys:
        cg = c[c['_key'] == k].sort_values('m_phi')
        fg = f[f['_key'] == k].sort_values('m_phi')
        cvals = cg['m_phi'].to_numpy(dtype=float)
        fvals = fg['m_phi'].to_numpy(dtype=float)
        if len(cvals) and len(fvals):
            dmat_cf = np.abs(cvals[:, None] - fvals[None, :])
            amb_colab += int(np.sum((dmat_cf <= tol).sum(axis=1) > 1))
            dmat_fc = np.abs(fvals[:, None] - cvals[None, :])
            amb_fixed += int(np.sum((dmat_fc <= tol).sum(axis=1) > 1))

    dup_colab = int(m.duplicated(subset=['lambda1', 'lambda6', 'tan_beta', 'm_phi_colab']).sum()) if not m.empty else 0
    dup_fixed = int(m.duplicated(subset=['lambda1', 'lambda6', 'tan_beta', 'm_phi_fixed']).sum()) if not m.empty else 0

    both_ok, only_colab_ok, only_fixed_ok, both_fail = mask_counts(m)

    rows.append({
        'scope': 'global', 'mphi_tol': tol,
        'lambda1': np.nan, 'lambda6': np.nan, 'tan_beta': np.nan,
        'n_colab': n_colab, 'n_fixed': n_fixed, 'n_matched': n_matched,
        'unmatched_colab': unmatched_colab, 'unmatched_fixed': unmatched_fixed,
        'both_ok': both_ok, 'only_colab_ok': only_colab_ok, 'only_fixed_ok': only_fixed_ok, 'both_fail': both_fail,
        'ambiguous_colab_points': amb_colab, 'ambiguous_fixed_points': amb_fixed,
        'duplicate_matched_colab_rows': dup_colab, 'duplicate_matched_fixed_rows': dup_fixed,
    })

    if not m.empty:
        for (l1, l6, tb), g in m.groupby(['lambda1', 'lambda6', 'tan_beta'], dropna=False):
            mask_c = (merged_df['lambda1'] == l1) & (merged_df['lambda6'] == l6) & (merged_df['tan_beta'] == tb)
            mask_f = (fixed_df['lambda1'] == l1) & (fixed_df['lambda6'] == l6) & (fixed_df['tan_beta'] == tb)
            rows.append({
                'scope': 'group', 'mphi_tol': tol,
                'lambda1': l1, 'lambda6': l6, 'tan_beta': tb,
                'n_colab': int(mask_c.sum()),
                'n_fixed': int(mask_f.sum()),
                'n_matched': len(g),
                'unmatched_colab': int(mask_c.sum()) - len(g),
                'unmatched_fixed': int(mask_f.sum()) - len(g),
                'both_ok': int(g['both_ok'].sum()),
                'only_colab_ok': int(g['only_colab_ok'].sum()),
                'only_fixed_ok': int(g['only_fixed_ok'].sum()),
                'both_fail': int(g['both_fail'].sum()),
                'ambiguous_colab_points': np.nan,
                'ambiguous_fixed_points': np.nan,
                'duplicate_matched_colab_rows': int(g.duplicated(subset=['m_phi_colab']).sum()),
                'duplicate_matched_fixed_rows': int(g.duplicated(subset=['m_phi_fixed']).sum()),
            })

diag = pd.DataFrame(rows)
out_dir = root / 'scripts/out/christopher_fixed_lam1_2026apr/comparison'
out_dir.mkdir(parents=True, exist_ok=True)
csv_path = out_dir / 'matching_tolerance_diagnostics.csv'
diag.to_csv(csv_path, index=False)

lines = []
lines.append('# Matching tolerance diagnostics')
lines.append('')
lines.append(f'- merged_csv: {merged_path}')
lines.append(f'- campaign_root: {mod.DEFAULT_CAMPAIGN_ROOT}')
lines.append('- tolerances tested: 1e-6, 1e-5, 1e-4, 1e-3')
lines.append('')
for tol in tols:
    g = diag[(diag.scope == 'global') & (diag.mphi_tol == tol)].iloc[0]
    lines.append(f'## mphi_tol={tol:g}')
    lines.append(f"- n_colab={int(g.n_colab)}")
    lines.append(f"- n_fixed={int(g.n_fixed)}")
    lines.append(f"- n_matched={int(g.n_matched)}")
    lines.append(f"- unmatched_colab={int(g.unmatched_colab)}")
    lines.append(f"- unmatched_fixed={int(g.unmatched_fixed)}")
    lines.append(f"- mask_confusion: both_ok={int(g.both_ok)}, only_colab_ok={int(g.only_colab_ok)}, only_fixed_ok={int(g.only_fixed_ok)}, both_fail={int(g.both_fail)}")
    lines.append(f"- ambiguity: ambiguous_colab_points={int(g.ambiguous_colab_points)}, ambiguous_fixed_points={int(g.ambiguous_fixed_points)}")
    lines.append(f"- duplicate matches in matched output: duplicate_colab_rows={int(g.duplicate_matched_colab_rows)}, duplicate_fixed_rows={int(g.duplicate_matched_fixed_rows)}")
    lines.append('')
    sub = diag[(diag.scope == 'group') & (diag.mphi_tol == tol)].copy()
    if sub.empty:
        lines.append('(no matched groups)')
        lines.append('')
    else:
        sub = sub[['lambda1', 'lambda6', 'tan_beta', 'n_colab', 'n_fixed', 'n_matched', 'unmatched_colab', 'unmatched_fixed', 'both_ok', 'only_colab_ok', 'only_fixed_ok', 'both_fail']].sort_values(['lambda1', 'lambda6', 'tan_beta'])
        try:
            lines.append(sub.to_markdown(index=False))
        except Exception:
            cols = list(sub.columns)
            header = '| ' + ' | '.join(cols) + ' |'
            sep = '| ' + ' | '.join(['---'] * len(cols)) + ' |'
            body = [header, sep]
            for _, rr in sub.iterrows():
                body.append('| ' + ' | '.join(str(rr[c]) for c in cols) + ' |')
            lines.extend(body)
        lines.append('')

lines.append('## Recommendation')
safe = []
for tol in tols:
    g = diag[(diag.scope == 'global') & (diag.mphi_tol == tol)].iloc[0]
    if int(g.ambiguous_colab_points) == 0 and int(g.ambiguous_fixed_points) == 0 and int(g.duplicate_matched_colab_rows) == 0 and int(g.duplicate_matched_fixed_rows) == 0:
        safe.append((tol, int(g.n_matched)))
if safe:
    best = sorted(safe, key=lambda x: x[0])[0]
    lines.append(f'- Smallest safe tolerance with no ambiguity/duplicates: mphi_tol={best[0]:g} (n_matched={best[1]}).')
else:
    lines.append('- No tested tolerance met the no-ambiguity criterion.')

md_path = out_dir / 'matching_tolerance_diagnostics.md'
md_path.write_text('\n'.join(lines) + '\n', encoding='utf-8')

print(csv_path)
print(md_path)
print(diag[diag.scope == 'global'][['mphi_tol', 'n_colab', 'n_fixed', 'n_matched', 'unmatched_colab', 'unmatched_fixed', 'both_ok', 'only_colab_ok', 'only_fixed_ok', 'both_fail', 'ambiguous_colab_points', 'ambiguous_fixed_points', 'duplicate_matched_colab_rows', 'duplicate_matched_fixed_rows']].to_string(index=False))