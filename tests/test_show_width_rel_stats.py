import subprocess
import sys
from pathlib import Path
import pandas as pd

FIX = Path('tests/fixtures')


def _make_merged(tmp_path):
    from scripts.validate_colab import build_merged_comparison
    df = build_merged_comparison(FIX/'colab_points_minimal.csv', FIX/'calc_points_minimal.csv')
    out = tmp_path/'merged.csv'
    df.to_csv(out, index=False)
    return out


def test_filters_triple_ok_when_present(tmp_path):
    merged = _make_merged(tmp_path)
    cmd = [sys.executable, 'scripts/show_width_rel_stats.py', str(merged), '--out-csv', str(tmp_path/'summary.csv')]
    p = subprocess.run(cmd, capture_output=True, text=True)
    assert p.returncode == 0
    s = pd.read_csv(tmp_path/'summary.csv')
    assert (s['n_finite'] == 1).all()


def test_max_point_cols_limits_output(tmp_path):
    merged = _make_merged(tmp_path)
    cmd = [sys.executable, 'scripts/show_width_rel_stats.py', str(merged), '--max-point', '--max-point-cols', 'point_id,rel_err_width_gg']
    p = subprocess.run(cmd, capture_output=True, text=True)
    assert p.returncode == 0
    assert 'point_id' in p.stdout and 'rel_err_width_gg' in p.stdout
    assert 'total_width_colab' not in p.stdout


def test_out_csv_writes_compact_summary(tmp_path):
    merged = _make_merged(tmp_path)
    out_csv = tmp_path/'summary.csv'
    p = subprocess.run([sys.executable, 'scripts/show_width_rel_stats.py', str(merged), '--out-csv', str(out_csv)], capture_output=True, text=True)
    assert p.returncode == 0
    s = pd.read_csv(out_csv)
    assert {'channel','mean_rel','median_rel','max_rel'}.issubset(set(s.columns))
