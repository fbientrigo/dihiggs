from pathlib import Path
import pytest

np = pytest.importorskip(
    "numpy",
    reason="historical comparison test requires numpy; not part of canonical v2 core",
)
pd = pytest.importorskip(
    "pandas",
    reason="historical comparison test requires pandas; not part of canonical v2 core",
)

from scripts.validate_colab import build_merged_comparison, safe_rel_err

FIX = Path('tests/fixtures')
pytestmark = pytest.mark.skipif(
    not (FIX / 'colab_points_minimal.csv').is_file(),
    reason='legacy comparison fixtures are not part of the maintained checkout',
)


def test_successful_merge_one_to_one():
    df = build_merged_comparison(FIX/'colab_points_minimal.csv', FIX/'calc_points_minimal.csv')
    assert len(df) == 2
    assert 'rel_err_width_gg' in df.columns


def test_duplicate_point_id_fails_clearly(tmp_path):
    colab = pd.read_csv(FIX/'colab_points_minimal.csv')
    colab = pd.concat([colab, colab.iloc[[0]]], ignore_index=True)
    cpath = tmp_path/'colab_dup.csv'
    colab.to_csv(cpath, index=False)
    with pytest.raises(ValueError, match='duplicate point_id'):
        build_merged_comparison(cpath, FIX/'calc_points_minimal.csv')


def test_missing_required_column_fails_clearly(tmp_path):
    colab = pd.read_csv(FIX/'colab_points_minimal.csv').drop(columns=['width_gg_colab'])
    cpath = tmp_path/'colab_missing.csv'
    colab.to_csv(cpath, index=False)
    with pytest.raises(ValueError, match='missing required columns'):
        build_merged_comparison(cpath, FIX/'calc_points_minimal.csv')


def test_rel_error_nan_for_zero_reference_not_inf():
    out = safe_rel_err(pd.Series([1.0, 0.0]), pd.Series([0.0, 0.0]))
    assert np.isnan(out.iloc[0]) and np.isnan(out.iloc[1])


def test_optional_width_cc_handling_works(tmp_path):
    calc = pd.read_csv(FIX/'calc_points_minimal.csv').drop(columns=['width_cc'])
    p = tmp_path/'calc_no_cc.csv'
    calc.to_csv(p, index=False)
    df = build_merged_comparison(FIX/'colab_points_minimal.csv', p)
    assert 'rel_err_width_cc' not in df.columns
