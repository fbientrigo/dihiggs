from __future__ import annotations

import copy
import csv
import json
import subprocess
import sys

import pytest

from autoresearch.harness.dihiggs_decision_report import (
    build_decision_report,
    build_decision_row,
    classify_decision,
    compute_derived_decision_metrics,
    load_registry_state,
    map_scores_by_run_dir,
    write_decision_report_json,
    write_decision_report_markdown,
    write_decision_rows_csv,
)


def _proposal(pid: str, run_dir: str | None, generation: int = 0) -> dict[str, object]:
    return {
        "schema_version": "1.0",
        "proposal_id": pid,
        "parent_id": None,
        "generation": generation,
        "operator": "seed",
        "bounds": {"lambda6": [1e-6, 0.5]},
        "fixed": {"mh": 125.0},
        "resolution": {"lambda6": 8},
        "budget": {"requested_points": 100, "max_walltime_sec": 60.0},
        "status": "planned",
        "created_utc": "2026-05-03T22:00:00Z",
        "updated_utc": "2026-05-03T22:00:00Z",
        "run_dir": run_dir,
        "metrics": {},
        "error": None,
    }


def _score(run_dir: str, **aggregate: object) -> dict[str, object]:
    return {
        "run_dir": run_dir,
        "rankable": bool(aggregate.get("rankable", True)),
        "aggregate_metrics": dict(aggregate),
        "task_summary_metrics": {
            "attempts_total": aggregate.get("attempts_total"),
            "triple_ok_rate_over_attempts": aggregate.get("triple_ok_rate_over_attempts"),
        },
    }


def test_classify_decision_not_rankable() -> None:
    out = classify_decision({"rankable": False})
    assert out["decision_label"] == "not_rankable"
    assert out["triggered_rules"]


def test_classify_decision_few_rows_promising_needs_more_points() -> None:
    out = classify_decision(
        {
            "rankable": True,
            "csv_rows": 5,
            "br_gaga_q95": 0.02,
            "br_bb_median": 0.05,
            "triple_ok_rate": 0.02,
        }
    )
    assert out["decision_label"] == "needs_more_points"
    assert out["triggered_rules"]


def test_classify_decision_high_br_gaga_suspicious_width_sanity_check() -> None:
    out = classify_decision(
        {
            "rankable": True,
            "csv_rows": 100,
            "br_gaga_q95": 0.02,
            "br_bb_median": 0.05,
            "triple_ok_rate": 0.02,
            "total_width_median": 1e-20,
        }
    )
    assert out["decision_label"] == "sanity_check"
    assert out["triggered_rules"]


def test_classify_decision_high_br_gaga_low_br_bb_acceptable_triple_ok_exploit() -> None:
    out = classify_decision(
        {
            "rankable": True,
            "csv_rows": 100,
            "br_gaga_q95": 0.02,
            "br_bb_median": 0.05,
            "triple_ok_rate": 0.02,
            "total_width_median": 1.0,
        }
    )
    assert out["decision_label"] == "exploit"
    assert out["triggered_rules"]


def test_classify_decision_weak_metrics_discard() -> None:
    out = classify_decision(
        {
            "rankable": True,
            "csv_rows": 100,
            "br_gaga_q95": 0.001,
            "br_bb_median": 0.2,
            "triple_ok_rate": 0.0,
        }
    )
    assert out["decision_label"] == "discard"
    assert out["triggered_rules"]


def test_compute_derived_metrics_prefers_attempts_rate() -> None:
    out = compute_derived_decision_metrics(
        {
            "csv_rows": 10,
            "attempts_total": 20,
            "triple_ok_rate_over_attempts": 0.2,
            "triple_ok_rate_over_csv_rows": 0.8,
            "br_gaga_q95": 0.05,
            "br_gaga_median": 0.02,
            "br_bb_median": 0.1,
        }
    )
    assert out["feasibility_score"] == 0.2
    assert out["physics_score"] == 0.05 - 0.1 * 0.1
    assert out["robustness_gap_br_gaga"] == pytest.approx(0.03)
    assert out["efficiency_score"] == 0.5
    assert out["novelty_score"] is None
    assert out["coverage_warning"] == "not_computed_yet"


def test_compute_derived_metrics_falls_back_to_csv_rows_rate() -> None:
    out = compute_derived_decision_metrics(
        {
            "csv_rows": 10,
            "triple_ok_rate_over_csv_rows": 0.3,
        }
    )
    assert out["feasibility_score"] == 0.3


def test_build_decision_row_required_fields_and_no_mutation() -> None:
    proposal = _proposal("p1", "/tmp/run1")
    before = copy.deepcopy(proposal)
    row = build_decision_row(
        proposal,
        _score(
            "/tmp/run1",
            rankable=True,
            csv_rows=100,
            br_gaga_q95=0.02,
            br_gaga_median=0.01,
            br_bb_median=0.05,
            total_width_median=1.0,
            objective_score=0.015,
            triple_ok_rate_over_csv_rows=0.02,
            attempts_total=200,
            triple_ok_rate_over_attempts=0.015,
        ),
    )
    required = {
        "proposal_id",
        "parent_id",
        "generation",
        "operator",
        "status",
        "run_dir",
        "bounds",
        "requested_points",
        "max_walltime_sec",
        "rankable",
        "csv_rows",
        "triple_ok_rate",
        "br_gaga_q95",
        "br_gaga_median",
        "br_bb_median",
        "total_width_median",
        "objective_score",
        "feasibility_score",
        "physics_score",
        "robustness_gap_br_gaga",
        "efficiency_score",
        "novelty_score",
        "coverage_warning",
        "health_status",
        "health_scoreable",
        "health_recommended_action",
        "dry_run",
        "decision_label",
        "decision_reason",
        "triggered_rules",
        "warnings",
        "recommended_next_action",
    }
    assert required.issubset(set(row.keys()))
    assert proposal == before


def test_build_decision_row_missing_score_is_not_rankable() -> None:
    row = build_decision_row(_proposal("p2", None), score=None)
    assert row["decision_label"] == "not_rankable"
    assert row["recommended_next_action"] == "obtain valid scored outputs first"
    assert row["health_status"] is None
    assert row["health_scoreable"] is None
    assert row["health_recommended_action"] is None
    assert row["dry_run"] is None


def test_build_decision_row_includes_health_fields() -> None:
    proposal = _proposal("p3", "/tmp/run3")
    score = _score(
        "/tmp/run3",
        rankable=True,
        csv_rows=50,
        br_gaga_q95=0.01,
        br_gaga_median=0.005,
        br_bb_median=0.4,
        total_width_median=1.0,
    )
    score["health"] = {
        "health_status": "warning",
        "is_scoreable": True,
        "recommended_action": "inspect_task_summary",
        "dry_run": False,
    }
    row = build_decision_row(proposal, score=score)
    assert row["health_status"] == "warning"
    assert row["health_scoreable"] is True
    assert row["health_recommended_action"] == "inspect_task_summary"
    assert row["dry_run"] is False


def test_build_decision_report_counts_sorts_and_best_ids() -> None:
    proposals = {
        "a": _proposal("a", "/runs/a"),
        "b": _proposal("b", "/runs/b"),
        "c": _proposal("c", "/runs/c"),
    }
    scores = map_scores_by_run_dir(
        [
            _score(
                "/runs/a",
                rankable=True,
                csv_rows=100,
                br_gaga_q95=0.02,
                br_gaga_median=0.01,
                br_bb_median=0.05,
                total_width_median=1.0,
                objective_score=0.2,
                triple_ok_rate_over_csv_rows=0.02,
            ),
            _score(
                "/runs/b",
                rankable=True,
                csv_rows=100,
                br_gaga_q95=0.001,
                br_gaga_median=0.001,
                br_bb_median=0.3,
                total_width_median=1.0,
                objective_score=-0.1,
                triple_ok_rate_over_csv_rows=0.0,
            ),
            _score("/runs/c", rankable=False, csv_rows=0),
        ]
    )
    report = build_decision_report(proposals, scores)
    summary = report["summary"]
    assert summary["n_proposals"] == 3
    assert summary["n_rankable"] == 2
    assert summary["n_score_rankable"] == 2
    assert summary["n_decision_rankable"] == 2
    assert summary["n_exploit"] == 1
    assert summary["n_discard"] == 1
    assert summary["n_not_rankable"] == 1
    assert summary["best_proposal_id_by_physics_score"] == "a"
    assert summary["best_proposal_id_by_objective_score"] == "a"
    labels = [row["decision_label"] for row in report["rows"]]
    assert labels == ["exploit", "discard", "not_rankable"]


def test_build_decision_report_empty_mapping() -> None:
    report = build_decision_report({}, {})
    assert report["summary"]["n_proposals"] == 0
    assert report["rows"] == []


def test_writers_json_markdown_csv(tmp_path) -> None:
    report = {
        "schema_version": "x",
        "summary": {"n_proposals": 1},
        "rows": [
            {
                "proposal_id": "p1",
                "decision_label": "explore",
                "triggered_rules": ["rule:mixed"],
                "bounds": {"x": [0, 1]},
            }
        ],
    }
    json_path = tmp_path / "report.json"
    md_path = tmp_path / "report.md"
    csv_path = tmp_path / "rows.csv"

    write_decision_report_json(report, json_path)
    loaded = json.loads(json_path.read_text(encoding="utf-8"))
    assert loaded["schema_version"] == "x"

    write_decision_report_markdown(report, md_path)
    text = md_path.read_text(encoding="utf-8")
    assert "# DiHiggs Autoresearch Decision Report" in text
    assert "## Summary" in text
    assert "## Top Candidates" in text
    assert "## Sanity Checks" in text

    write_decision_rows_csv(report, csv_path)
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    assert "triggered_rules" in rows[0]
    assert rows[0]["triggered_rules"] == '["rule:mixed"]'
    assert rows[0]["bounds"] == '{"x":[0,1]}'


def test_build_decision_report_health_summary_counts() -> None:
    proposals = {
        "ok": _proposal("ok", "/runs/ok"),
        "warn": _proposal("warn", "/runs/warn"),
        "fail": _proposal("fail", "/runs/fail"),
    }
    s_ok = _score("/runs/ok", rankable=True, csv_rows=10, br_gaga_q95=0.01, br_bb_median=0.4)
    s_ok["health"] = {"health_status": "ok", "is_scoreable": True, "recommended_action": "score", "dry_run": False}
    s_warn = _score("/runs/warn", rankable=True, csv_rows=10, br_gaga_q95=0.01, br_bb_median=0.4)
    s_warn["health"] = {"health_status": "warning", "is_scoreable": True, "recommended_action": "inspect_task_summary", "dry_run": False}
    s_fail = _score("/runs/fail", rankable=False, csv_rows=0)
    s_fail["health"] = {"health_status": "failed", "is_scoreable": False, "recommended_action": "reject_dry_run", "dry_run": True}

    report = build_decision_report(proposals, map_scores_by_run_dir([s_ok, s_warn, s_fail]))
    summary = report["summary"]
    assert summary["n_health_ok"] == 1
    assert summary["n_health_warning"] == 1
    assert summary["n_health_failed"] == 1
    assert summary["n_health_scoreable"] == 2
    assert summary["n_health_dry_run"] == 1


def test_load_registry_state_and_cli(tmp_path) -> None:
    registry = tmp_path / "registry.jsonl"
    run_dir = tmp_path / "run1"
    tb = run_dir / "tb_1000"
    tb.mkdir(parents=True)

    proposal = _proposal("pcli", str(run_dir))
    registry.write_text(json.dumps({"event_type": "create_proposal", "proposal": proposal}) + "\n", encoding="utf-8")

    (run_dir / "task_summary.jsonl").write_text(
        '{"status":"done","attempts":10,"triple_ok_points":1}\n',
        encoding="utf-8",
    )
    (tb / "scan_tb_1000.csv").write_text(
        "m_phi,lambda6,lambda7,tan_beta,lam1,positivity_ok,unitarity_ok,perturbativity_ok,width_bb,width_gaga,width_Zga,width_gg,width_WW,width_ZZ,width_hh,total_width,br_gaga\n"
        "250,0.0006,0.0,30000,6.0,1,1,1,1.0,1e-5,1e-5,0.1,0.1,0.1,0.1,10.0,0.02\n",
        encoding="utf-8",
    )

    state = load_registry_state(registry)
    assert "pcli" in state

    out_json = tmp_path / "out.json"
    out_md = tmp_path / "out.md"
    out_csv = tmp_path / "out.csv"

    cmd = [
        sys.executable,
        "-m",
        "autoresearch.harness.dihiggs_decision_report",
        "--registry",
        str(registry),
        "--run-dir",
        str(run_dir),
        "--out-json",
        str(out_json),
        "--out-md",
        str(out_md),
        "--out-csv",
        str(out_csv),
    ]
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    assert result.returncode == 0
    assert out_json.exists()
    assert out_md.exists()
    assert out_csv.exists()


def test_cli_reports_dry_run_as_non_rankable_with_health_counts(tmp_path) -> None:
    registry = tmp_path / "registry.jsonl"
    run_dir = tmp_path / "run_dry"
    run_dir.mkdir(parents=True)
    (run_dir / "run_manifest.json").write_text('{"runtime":{"dry_run":true}}', encoding="utf-8")

    proposal = _proposal("pdry", str(run_dir))
    registry.write_text(json.dumps({"event_type": "create_proposal", "proposal": proposal}) + "\n", encoding="utf-8")

    out_json = tmp_path / "dry.json"
    cmd = [
        sys.executable,
        "-m",
        "autoresearch.harness.dihiggs_decision_report",
        "--registry",
        str(registry),
        "--run-dir",
        str(run_dir),
        "--out-json",
        str(out_json),
    ]
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    assert result.returncode == 0
    payload = json.loads(out_json.read_text(encoding="utf-8"))
    row = payload["rows"][0]
    assert row["rankable"] is False
    assert row["health_status"] == "failed"
    assert row["health_scoreable"] is False
    assert row["dry_run"] is True
    assert payload["summary"]["n_health_failed"] == 1
    assert payload["summary"]["n_health_dry_run"] == 1
