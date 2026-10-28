from __future__ import annotations

import json
from pathlib import Path

import pytest

from dihiggs.app.orchestrator import analytic_bridge
from dihiggs.app.orchestrator.aligned_lambda7 import build_aligned_lambda7_family
from dihiggs.app.orchestrator.cli import build_parser, main


def _point_file(tmp_path: Path) -> Path:
    point = build_aligned_lambda7_family(
        tan_beta="100", lambda1_ref="1.2", mh_GeV="125.13",
        mphi_GeV="150", mA_GeV="200", mHp_GeV="200", v_GeV="246",
        X_lambda6_tanbeta="1",
    )
    path = tmp_path / "point.json"
    path.write_text(json.dumps(point.as_dict()))
    return path


def test_t4_cli_defaults_preserve_old_engine_and_legacy_mode():
    args = build_parser().parse_args(["--engine", "t4", "--point-json", "p", "--output", "o"])
    assert args.engine == "t4"
    assert args.t4_mode == "legacy"
    assert build_parser().parse_args([]).engine == "lambda1_v2"


def test_t4_cli_dry_run_validates_and_never_calls_backend(tmp_path, monkeypatch, capsys):
    point = _point_file(tmp_path)
    output = tmp_path / "out/result.json"

    def backend_must_not_run(*args, **kwargs):
        raise AssertionError("dry-run executed the calculator")

    monkeypatch.setattr(analytic_bridge, "evaluate_t4", backend_must_not_run)
    assert main(["--engine", "t4", "--t4-mode", "analytic", "--point-json", str(point),
                 "--output", str(output), "--dry-run"]) == 0
    result = json.loads(output.read_text())
    assert result["constructor_status"] == "NOT_TESTED"
    assert result["dry_run"] is True
    assert result["requested_mode"] == "analytic"
    assert "not executed" in result["dry_run_reason"]
    assert "NOT_TESTED" in capsys.readouterr().out


def test_t4_cli_dry_run_still_rejects_tampered_point(tmp_path, capsys):
    point = _point_file(tmp_path)
    data = json.loads(point.read_text())
    data["values"]["Z7"] = "0"
    point.write_text(json.dumps(data))
    assert main(["--engine", "t4", "--point-json", str(point),
                 "--output", str(tmp_path / "out.json"), "--dry-run"]) == 2
    assert "derived field Z7" in capsys.readouterr().err


def test_t4_cli_dispatches_validated_analytic_point_to_backend(tmp_path):
    point = _point_file(tmp_path)
    output = tmp_path / "out/result.json"
    bridge = Path(__file__).resolve().parents[2] / "campaigns/t2_t4_20261008T043858Z/scripts/analytic_bridge"
    assert main(["--engine", "t4", "--t4-mode", "analytic", "--point-json", str(point),
                 "--output", str(output), "--exec", str(bridge)]) == 0
    result = json.loads(output.read_text())
    assert result["requested_mode"] == "analytic"
    assert result["constructor_status"] == "INCONCLUSIVE"
    assert result["backend_constructor_status"] == "PASS"
    assert result["family_constructor_agreement_status"] == "NOT_TESTED"
    assert result["mode"] == result["requested_mode"]
