import pytest

pd = pytest.importorskip(
    "pandas",
    reason="historical comparison test requires pandas; not part of canonical v2 core",
)

from scripts.compare_christopher_fixed_lam1 import (
    _normalize_merged,
    build_group_summary,
    build_mask_confusion_by_group,
    match_points,
)


def _mk_colab_df():
    return pd.DataFrame(
        [
            {
                "mH": 125.0,
                "lambda1": 0.4,
                "lambda6": 0.01,
                "tan_beta": 10000.0,
                "triple_ok": 1,
                "width_bb": 10.0,
                "width_cc": 1.0,
                "width_tautau": 2.0,
                "width_gg": 3.0,
                "width_gaga": 4.0,
                "width_Zga": 5.0,
                "total_width": 30.0,
                "br_gaga_2hdmc": 0.2,
            },
            {
                "mH": 126.0,
                "lambda1": 0.4,
                "lambda6": 0.01,
                "tan_beta": 10000.0,
                "triple_ok": 0,
                "width_bb": 11.0,
                "width_cc": 1.1,
                "width_tautau": 2.1,
                "width_gg": 3.1,
                "width_gaga": 4.1,
                "width_Zga": 5.1,
                "total_width": 31.0,
                "br_gaga_2hdmc": 0.21,
            },
        ]
    )


def _mk_fixed_df(include_cc: bool = False):
    rows = [
        {
            "m_phi": 125.0,
            "lam1": 0.4,
            "lambda6": 0.01,
            "tan_beta": 10000.0,
            "positivity_ok": 1,
            "unitarity_ok": 1,
            "perturbativity_ok": 1,
            "width_bb": 10.0,
            "width_tautau": 2.0,
            "width_gg": 3.0,
            "width_gaga": 4.0,
            "width_Zga": 5.0,
            "total_width": 30.0,
            "br_gaga": 0.2,
        },
        {
            "m_phi": 126.0,
            "lam1": 0.4,
            "lambda6": 0.01,
            "tan_beta": 10000.0,
            "positivity_ok": 1,
            "unitarity_ok": 1,
            "perturbativity_ok": 1,
            "width_bb": 11.0,
            "width_tautau": 2.1,
            "width_gg": 3.1,
            "width_gaga": 4.1,
            "width_Zga": 5.1,
            "total_width": 31.0,
            "br_gaga": 0.21,
        },
    ]
    if include_cc:
        for r in rows:
            r["width_cc"] = 1.0
    return pd.DataFrame(rows)


def _normalize_fixed_for_match(df: pd.DataFrame) -> pd.DataFrame:
    out = df.rename(columns={"lam1": "lambda1"}).copy()
    out["triple_ok_fixed_campaign"] = out["positivity_ok"].astype(bool) & out["unitarity_ok"].astype(bool) & out["perturbativity_ok"].astype(bool)
    return out


def test_exact_key_matching():
    colab = _normalize_merged(_mk_colab_df())
    fixed = _normalize_fixed_for_match(_mk_fixed_df())
    matched = match_points(colab, fixed, mphi_tol=1e-6, param_tol=1e-8)
    assert len(matched) == 2
    assert matched["m_phi_abs_diff"].max() == 0.0


def test_nearest_mphi_matching_within_tolerance():
    colab = _normalize_merged(_mk_colab_df())
    fixed = _normalize_fixed_for_match(_mk_fixed_df())
    fixed.loc[fixed["m_phi"] == 126.0, "m_phi"] = 126.0 + 5e-7
    matched = match_points(colab, fixed, mphi_tol=1e-6, param_tol=1e-8)
    assert len(matched) == 2
    assert (matched["m_phi_abs_diff"] <= 1e-6).all()


def test_physical_mask_confusion_counts():
    colab = _normalize_merged(_mk_colab_df())
    fixed = _normalize_fixed_for_match(_mk_fixed_df())
    matched = match_points(colab, fixed, mphi_tol=1e-6, param_tol=1e-8)
    # row1 both ok, row2 only fixed ok
    assert int(matched["both_ok"].sum()) == 1
    assert int(matched["only_fixed_ok"].sum()) == 1
    assert int(matched["only_colab_ok"].sum()) == 0
    assert int(matched["both_fail"].sum()) == 0

    conf = build_mask_confusion_by_group(matched)
    assert int(conf["both_ok"].sum()) == 1
    assert int(conf["only_fixed_ok"].sum()) == 1


def test_missing_optional_width_cc_does_not_crash():
    colab = _normalize_merged(_mk_colab_df())
    fixed = _normalize_fixed_for_match(_mk_fixed_df(include_cc=False))
    matched = match_points(colab, fixed, mphi_tol=1e-6, param_tol=1e-8)
    assert len(matched) == 2
    assert "width_cc_colab" in matched.columns
    assert "width_cc_fixed" in matched.columns


def test_group_summary_metrics_compute_correctly():
    colab = _normalize_merged(_mk_colab_df())
    fixed = _normalize_fixed_for_match(_mk_fixed_df())
    matched = match_points(colab, fixed, mphi_tol=1e-6, param_tol=1e-8)
    summary = build_group_summary(colab, fixed, matched)
    assert len(summary) == 1
    r = summary.iloc[0]
    assert int(r["n_colab"]) == 2
    assert int(r["n_fixed"]) == 2
    assert int(r["n_matched"]) == 2
    assert int(r["n_triple_ok_colab"]) == 1
    assert int(r["n_triple_ok_fixed"]) == 2
    assert int(r["n_both_ok"]) == 1
    assert int(r["n_only_fixed_ok"]) == 1
    assert abs(float(r["fraction_mask_agreement"]) - 0.5) < 1e-12
