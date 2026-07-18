import sys
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.validate_colab import safe_rel_err, safe_fraction, summarize_comparison


def test_safe_rel_err_nan_when_reference_zero():
    a = pd.Series([1.0, 2.0, 0.0])
    b = pd.Series([1.0, 0.0, 0.0])
    out = safe_rel_err(a, b)
    assert np.isclose(out.iloc[0], 0.0)
    assert np.isnan(out.iloc[1])
    assert np.isnan(out.iloc[2])


def test_safe_fraction_nan_when_denominator_zero():
    num = pd.Series([1.0, 2.0, 3.0])
    den = pd.Series([2.0, 0.0, -0.0])
    out = safe_fraction(num, den)
    assert np.isclose(out.iloc[0], 0.5)
    assert np.isnan(out.iloc[1])
    assert np.isnan(out.iloc[2])


def test_summarize_comparison_has_expected_counts():
    df = pd.DataFrame(
        {
            "set_param_ok": [1, 1, 0],
            "triple_ok": [1, 0, 0],
            "warning_flag": [1, 0, 1],
            "abs_err_br_bb": [0.1, 0.2, 0.3],
            "abs_err_br_gaga": [0.01, 0.02, 0.03],
            "abs_err_br_Zga": [0.03, 0.04, 0.05],
            "rel_err_br_bb": [0.1, 0.2, 0.3],
            "rel_err_br_gaga": [0.1, 0.2, 0.3],
            "rel_err_br_Zga": [0.1, 0.2, 0.3],
            "abs_err_width_bb": [1e-3, 2e-3, 3e-3],
            "abs_err_width_tautau": [1e-4, 2e-4, 3e-4],
            "abs_err_width_gg": [1e-5, 2e-5, 3e-5],
            "abs_err_width_gaga": [1e-6, 2e-6, 3e-6],
            "abs_err_width_Zga": [1e-6, 2e-6, 3e-6],
            "abs_err_width_total": [1e-2, 2e-2, 3e-2],
            "rel_err_width_total": [0.01, 0.02, 0.03],
        }
    )
    s = summarize_comparison(df)
    assert int(s["n_points"]) == 3
    assert int(s["n_set_param_ok"]) == 2
    assert int(s["n_triple_ok"]) == 1
    assert int(s["n_warning_flag"]) == 2
