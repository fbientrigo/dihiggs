from __future__ import annotations

import math

import pytest

from autoresearch.harness.dihiggs_physics_score import (
    finite_float,
    quantile,
    read_physics_csv,
    score_csv,
    score_rows,
    score_run_dir,
    score_run_dir_checked,
    score_task_summary,
)


def test_finite_float_accepts_int_float_and_scientific_string() -> None:
    assert finite_float(3) == 3.0
    assert finite_float(2.5) == 2.5
    assert finite_float("1e-3") == 1e-3


def test_finite_float_rejects_bool_none_empty_nan_inf_and_invalid() -> None:
    assert finite_float(True) is None
    assert finite_float(None) is None
    assert finite_float("") is None
    assert finite_float("nan") is None
    assert finite_float(float("nan")) is None
    assert finite_float("inf") is None
    assert finite_float(float("inf")) is None
    assert finite_float("abc") is None


def test_quantile_empty_returns_none() -> None:
    assert quantile([], 0.5) is None


def test_quantile_q0_is_min_q1_is_max() -> None:
    vals = [3.0, 1.0, 2.0]
    assert quantile(vals, 0.0) == 1.0
    assert quantile(vals, 1.0) == 3.0


def test_quantile_median_and_q95_are_deterministic() -> None:
    vals = [1.0, 2.0, 3.0, 4.0]
    assert quantile(vals, 0.5) == 2.5
    assert quantile(vals, 0.95) == pytest.approx(3.85)


def test_quantile_rejects_q_out_of_range() -> None:
    with pytest.raises(ValueError):
        quantile([1.0], -0.01)
    with pytest.raises(ValueError):
        quantile([1.0], 1.01)


def test_read_physics_csv_missing_raises(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        read_physics_csv(tmp_path / "missing.csv")


def test_read_physics_csv_header_only_returns_empty(tmp_path) -> None:
    path = tmp_path / "h.csv"
    path.write_text("a,b\n", encoding="utf-8")
    assert read_physics_csv(path) == []


def test_read_physics_csv_simple_rows_load(tmp_path) -> None:
    path = tmp_path / "x.csv"
    path.write_text("m_phi,br_gaga\n125,1e-4\n126,2e-4\n", encoding="utf-8")
    rows = read_physics_csv(path)
    assert len(rows) == 2
    assert rows[0]["m_phi"] == "125"


def test_score_rows_empty_not_rankable_with_warnings() -> None:
    result = score_rows([])
    assert result["rankable"] is False
    assert result["objective_score"] is None
    assert result["warnings"]


def test_score_rows_counts_triple_ok_and_derives_br_bb_and_total_width_rules() -> None:
    rows = [
        {
            "positivity_ok": 1,
            "unitarity_ok": 1,
            "perturbativity_ok": 1,
            "br_gaga": "0.01",
            "width_bb": "2.0",
            "total_width": "10.0",
        },
        {
            "positivity_ok": 1,
            "unitarity_ok": 0,
            "perturbativity_ok": 1,
            "br_gaga": "bad",
            "width_bb": "3.0",
            "total_width": "0.0",
            "br_bb": "0.4",
        },
        {
            "positivity_ok": 1,
            "unitarity_ok": 1,
            "perturbativity_ok": 1,
            "br_gaga": "0.03",
            "width_bb": "3.0",
            "total_width": "12.0",
        },
    ]

    result = score_rows(rows)
    assert result["csv_rows"] == 3
    assert result["triple_ok_count"] == 2
    assert result["triple_ok_rate_over_csv_rows"] == pytest.approx(2 / 3)

    assert result["br_gaga_count"] == 2
    assert result["br_gaga_max"] == 0.03

    assert result["br_bb_count"] == 3
    assert result["br_bb_min"] == 0.2
    assert result["total_width_count"] == 2
    assert result["total_width_median"] == 11.0

    assert result["rankable"] is True
    assert result["objective_score"] is not None


def test_score_rows_objective_none_when_not_rankable() -> None:
    rows = [{"br_gaga": "bad", "width_bb": "bad", "total_width": "bad"}]
    result = score_rows(rows)
    assert result["rankable"] is False
    assert result["objective_score"] is None


def test_score_task_summary_missing_returns_zeroes_and_warning(tmp_path) -> None:
    result = score_task_summary(tmp_path / "task_summary.jsonl")
    assert result["tasks_total"] == 0
    assert result["attempts_total"] == 0
    assert result["warnings"]


def test_score_task_summary_aggregates_done_skip_fail_and_timeout(tmp_path) -> None:
    path = tmp_path / "task_summary.jsonl"
    path.write_text(
        "\n".join(
            [
                '{"status": "done", "attempts": 2, "triple_ok_points": 1}',
                '{"status": "skipped"}',
                '{"status": "failed", "attempts": "3", "triple_ok_points": "2"}',
                '{"status": "timeout", "metrics": {"attempts": "1", "triple_ok_points": "1"}}',
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    result = score_task_summary(path)
    assert result["tasks_total"] == 4
    assert result["tasks_done"] == 1
    assert result["tasks_failed"] == 1
    assert result["tasks_skipped"] == 1
    assert result["tasks_timeout"] == 1
    assert result["attempts_total"] == 6
    assert result["triple_ok_total"] == 4
    assert result["triple_ok_rate_over_attempts"] == pytest.approx(4 / 6)


def test_score_task_summary_real_orchestrator_records(tmp_path) -> None:
    # Records as actually emitted by dihiggs/app/orchestrator/runner.py:
    # status is one of done|fail|crash|skip and attempts use total_attempts.
    path = tmp_path / "task_summary.jsonl"
    path.write_text(
        "\n".join(
            [
                '{"task_id": "tb_10000", "status": "done", "total_attempts": 100, "triple_ok_points": 42, "returncode": 0}',
                '{"task_id": "tb_20000", "status": "skip", "total_attempts": 0, "triple_ok_points": 0}',
                '{"task_id": "tb_30000", "status": "fail", "total_attempts": 7, "triple_ok_points": 0, "returncode": 1}',
                '{"task_id": "tb_40000", "status": "crash", "total_attempts": 3, "triple_ok_points": 0, "returncode": 134}',
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    result = score_task_summary(path)
    assert result["tasks_total"] == 4
    assert result["tasks_done"] == 1
    assert result["tasks_skipped"] == 1
    # fail + crash are both failures.
    assert result["tasks_failed"] == 2
    assert result["attempts_total"] == 110
    assert result["triple_ok_total"] == 42


def test_score_task_summary_metrics_total_attempts_fallback(tmp_path) -> None:
    path = tmp_path / "task_summary.jsonl"
    path.write_text(
        '{"status": "done", "metrics": {"total_attempts": 5, "triple_ok_points": 2}}\n',
        encoding="utf-8",
    )
    result = score_task_summary(path)
    assert result["attempts_total"] == 5
    assert result["triple_ok_total"] == 2


def test_score_task_summary_nonzero_returncode_counts_as_failure(tmp_path) -> None:
    # No recognized status string, but a nonzero exit code must still count.
    path = tmp_path / "task_summary.jsonl"
    path.write_text('{"task_id": "x", "returncode": 2, "total_attempts": 4}\n', encoding="utf-8")
    result = score_task_summary(path)
    assert result["tasks_failed"] == 1
    assert result["attempts_total"] == 4


def test_score_task_summary_total_attempts_takes_precedence(tmp_path) -> None:
    path = tmp_path / "task_summary.jsonl"
    path.write_text(
        '{"status": "done", "total_attempts": 9, "attempts": 1, "triple_ok_points": 3}\n',
        encoding="utf-8",
    )
    result = score_task_summary(path)
    assert result["attempts_total"] == 9


def test_score_task_summary_invalid_json_raises_line_number(tmp_path) -> None:
    path = tmp_path / "task_summary.jsonl"
    path.write_text('{"status": "done"}\n{bad json}\n', encoding="utf-8")
    with pytest.raises(ValueError, match="line 2"):
        score_task_summary(path)


def test_score_task_summary_string_numeric_triple_ok_points_parsed(tmp_path) -> None:
    path = tmp_path / "task_summary.jsonl"
    path.write_text('{"status": "done", "attempts": "2", "triple_ok_points": "1.0"}\n', encoding="utf-8")
    result = score_task_summary(path)
    assert result["attempts_total"] == 2
    assert result["triple_ok_total"] == 1


def test_score_csv_includes_path_and_matches_score_rows(tmp_path) -> None:
    path = tmp_path / "scan_tb_1.csv"
    path.write_text(
        "positivity_ok,unitarity_ok,perturbativity_ok,br_gaga,width_bb,total_width\n"
        "1,1,1,0.01,1.0,10.0\n",
        encoding="utf-8",
    )
    csv_result = score_csv(path)
    rows_result = score_rows(read_physics_csv(path))

    assert csv_result["csv_path"] == str(path)
    for key, value in rows_result.items():
        assert csv_result[key] == value


def test_score_run_dir_synthetic_with_two_csvs_and_summary(tmp_path) -> None:
    run_dir = tmp_path / "run"
    tb1 = run_dir / "tb_1000"
    tb2 = run_dir / "tb_2000"
    tb1.mkdir(parents=True)
    tb2.mkdir(parents=True)

    (run_dir / "task_summary.jsonl").write_text(
        '{"status": "done", "attempts": 2, "triple_ok_points": 1}\n',
        encoding="utf-8",
    )

    (tb1 / "scan_tb_1000.csv").write_text(
        "positivity_ok,unitarity_ok,perturbativity_ok,br_gaga,width_bb,total_width\n"
        "1,1,1,0.02,2.0,10.0\n",
        encoding="utf-8",
    )
    (tb2 / "scan_tb_2000.csv").write_text(
        "positivity_ok,unitarity_ok,perturbativity_ok,br_gaga,width_bb,total_width\n"
        "1,0,1,0.01,1.0,5.0\n",
        encoding="utf-8",
    )

    result = score_run_dir(run_dir)
    assert result["run_dir"] == str(run_dir)
    assert result["rankable"] is True
    assert len(result["csv_metrics"]) == 2
    agg = result["aggregate_metrics"]
    assert agg["csv_rows"] == 2
    assert agg["triple_ok_count"] == 1


def test_score_run_dir_no_csvs_is_not_rankable(tmp_path) -> None:
    run_dir = tmp_path / "run_empty"
    run_dir.mkdir()
    (run_dir / "task_summary.jsonl").write_text("", encoding="utf-8")

    result = score_run_dir(run_dir)
    assert result["csv_metrics"] == []
    assert result["rankable"] is False
    assert result["aggregate_metrics"]["csv_rows"] == 0


def _minimal_valid_scan_csv(path) -> None:
    path.write_text(
        "m_phi,lambda6,lambda7,tan_beta,lam1,positivity_ok,unitarity_ok,perturbativity_ok,width_bb,width_gaga,width_Zga,width_gg,width_WW,width_ZZ,width_hh,total_width,br_gaga\n"
        "250,0.0012,0.0,30000,6.0,1,1,1,1.0,1e-5,1e-5,0.1,0.1,0.1,0.1,10.0,0.02\n",
        encoding="utf-8",
    )


def test_score_run_dir_checked_rejects_dry_run_default(tmp_path) -> None:
    run_dir = tmp_path / "dry"
    run_dir.mkdir()
    (run_dir / "run_manifest.json").write_text('{"runtime":{"dry_run":true}}', encoding="utf-8")
    out = score_run_dir_checked(run_dir)
    assert out["rankable"] is False
    assert out["health"]["dry_run"] is True
    assert "run_dir_not_scoreable" in out["aggregate_metrics"]["warnings"]


def test_score_run_dir_checked_allows_warning_when_enabled(tmp_path) -> None:
    run_dir = tmp_path / "warn_ok"
    tb = run_dir / "tb_1000"
    tb.mkdir(parents=True)
    _minimal_valid_scan_csv(tb / "scan_tb_1000.csv")
    out = score_run_dir_checked(run_dir, allow_warnings=True)
    assert out["health"]["health_status"] == "warning"
    assert out["rankable"] is True


def test_score_run_dir_checked_blocks_warning_when_disabled(tmp_path) -> None:
    run_dir = tmp_path / "warn_block"
    tb = run_dir / "tb_1000"
    tb.mkdir(parents=True)
    _minimal_valid_scan_csv(tb / "scan_tb_1000.csv")
    out = score_run_dir_checked(run_dir, allow_warnings=False)
    assert out["rankable"] is False
    assert "run_dir_health_warning_blocked" in out["aggregate_metrics"]["warnings"]


def test_score_run_dir_checked_attaches_health_on_success(tmp_path) -> None:
    run_dir = tmp_path / "ok"
    tb = run_dir / "tb_1000"
    tb.mkdir(parents=True)
    (run_dir / "task_summary.jsonl").write_text('{"status":"done"}\n', encoding="utf-8")
    _minimal_valid_scan_csv(tb / "scan_tb_1000.csv")
    out = score_run_dir_checked(run_dir)
    assert isinstance(out.get("health"), dict)
    assert out["health"]["is_scoreable"] is True
    assert out["rankable"] is True


def test_score_run_dir_checked_allow_failed_true_proceeds(tmp_path) -> None:
    run_dir = tmp_path / "failed_but_forced"
    tb = run_dir / "tb_1000"
    tb.mkdir(parents=True)
    (tb / "scan_tb_1000.csv").write_text(
        "m_phi,lambda6,lambda7,tan_beta,lam1,positivity_ok,unitarity_ok,perturbativity_ok,width_bb,width_gaga,width_Zga,width_gg,width_WW,width_ZZ,width_hh,total_width,br_gaga\n",
        encoding="utf-8",
    )
    out = score_run_dir_checked(run_dir, allow_failed=True)
    assert isinstance(out.get("health"), dict)
    assert out["health"]["is_scoreable"] is False
    assert out["rankable"] is False
