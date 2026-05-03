import pandas as pd
from pathlib import Path

from scripts.validate_colab import build_merged_comparison
from scripts.group_width_errors_by_config import build_grouped_report

FIX = Path('tests/fixtures')


def test_groups_and_median_and_n():
    df = build_merged_comparison(FIX/'colab_points_minimal.csv', FIX/'calc_points_minimal.csv')
    out = build_grouped_report(df, ['lambda1','lambda6','tan_beta'])
    assert {'lambda1','lambda6','tan_beta','n','median_rel_gg'}.issubset(set(out.columns))
    assert out['n'].sum() == 1  # triple_ok filter retains only p1


def test_skips_missing_optional_channels(tmp_path):
    df = build_merged_comparison(FIX/'colab_points_minimal.csv', FIX/'calc_points_minimal.csv')
    df = df.drop(columns=['rel_err_width_cc'])
    out = build_grouped_report(df, ['lambda1','lambda6','tan_beta'])
    assert 'median_rel_cc' not in out.columns
