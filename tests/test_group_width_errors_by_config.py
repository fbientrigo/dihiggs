from pathlib import Path
import pytest

pd = pytest.importorskip(
    "pandas",
    reason="historical comparison test requires pandas; not part of canonical v2 core",
)

from scripts.validate_colab import build_merged_comparison
from scripts.group_width_errors_by_config import build_grouped_report

FIX = Path('tests/fixtures')
pytestmark = pytest.mark.skipif(
    not (FIX / 'colab_points_minimal.csv').is_file(),
    reason='legacy comparison fixtures are not part of the maintained checkout',
)


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
