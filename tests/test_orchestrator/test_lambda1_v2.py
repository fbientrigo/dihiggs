from __future__ import annotations

import json
import stat
import subprocess
from pathlib import Path

import pytest

from dihiggs.app.orchestrator.lambda1_v2 import (
    INPUT_HEADER,
    Lambda1Fixed,
    SCHEMA_VERSION,
    build_command,
    cartesian_rows,
    format_float64,
    grid_values,
    run_lambda1_v2,
    stable_point_id,
    validate_output,
    write_input_csv,
)


def rows() -> list[dict[str, str]]:
    return cartesian_rows(
        fixed=Lambda1Fixed(125.13, 300.0, 310.0, 0.995, 0.1, 0.0),
        mH_values=[130.0, 290.0], lambda1_values=[0.0, 12.0],
        tan_beta_values=[50.0],
    )


def test_exact_header_cartesian_order_and_stable_ids(tmp_path: Path) -> None:
    generated = rows()
    assert list(generated[0]) == list(INPUT_HEADER)
    assert len(generated) == 4
    assert [row["mH_gev"] for row in generated] == ["130", "130", "290", "290"]
    assert [row["lambda1_target"] for row in generated] == ["0", "12", "0", "12"]
    assert generated == cartesian_rows(
        fixed=Lambda1Fixed(125.13, 300.0, 310.0, 0.995, 0.1, 0.0),
        mH_values=[130.0, 290.0], lambda1_values=[0.0, 12.0], tan_beta_values=[50.0],
    )
    assert generated[0]["point_id"] == stable_point_id(
        (125.13, 130.0, 300.0, 310.0, 0.995, 50.0, 0.0, 0.1, 0.0)
    )
    input_csv = tmp_path / "input.csv"
    write_input_csv(input_csv, generated)
    assert input_csv.read_text().splitlines()[0] == ",".join(INPUT_HEADER)


def test_float64_format_is_round_trip_safe() -> None:
    value = 1.0 / 10.0
    assert float(format_float64(value)) == value
    assert format_float64(value) == "0.10000000000000001"
    assert grid_values(0.0, 1.0, 3) == [0.0, 0.5, 1.0]


def test_command_uses_only_input_and_output_paths() -> None:
    assert build_command(Path("/bin/Lambda1EvaluatorV2"), Path("in.csv"), Path("out.csv")) == [
        "/bin/Lambda1EvaluatorV2", "in.csv", "out.csv"
    ]


def test_schema_and_row_count_validation(tmp_path: Path) -> None:
    output = tmp_path / "output.csv"
    output.write_text("schema_version,point_id\ndihiggs.lambda1.v2,p0\n", encoding="utf-8")
    assert validate_output(output, 1)["schema"] == SCHEMA_VERSION
    with pytest.raises(ValueError, match="row-count"):
        validate_output(output, 2)


def test_run_preserves_manifest_provenance_and_row_count(tmp_path: Path) -> None:
    executable = tmp_path / "fake-evaluator"
    executable.write_text(
        "#!/bin/sh\n"
        "printf 'schema_version,point_id\\n' > \"$2\"\n"
        "tail -n +2 \"$1\" | cut -d, -f1 | while read -r id; do "
        "printf 'dihiggs.lambda1.v2,%s\\n' \"$id\" >> \"$2\"; done\n",
        encoding="utf-8",
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    manifest = run_lambda1_v2(
        executable=executable, rows=rows(), outdir=tmp_path / "out",
        campaign="pilot", run_name="run-1", repo_root=tmp_path,
    )
    assert manifest["status"] == "complete"
    assert manifest["schema"] == SCHEMA_VERSION
    assert manifest["row_count"] == 4
    assert manifest["requested_row_count"] == 4
    assert len(manifest["input_sha256"]) == 64
    assert manifest["mass_convention"]["convention_id"] == "mh_125p13_pdg_2026"
    assert manifest["mass_convention"]["m_h_GeV_text"] == "125.13"
    assert manifest["mass_convention"]["schema_version"] == "physics_conventions_v3"
    assert len(manifest["mass_convention"]["source_sha256"]) == 64
    saved = json.loads((tmp_path / "out/campaign=pilot/run-1/run_manifest.json").read_text())
    assert saved["command"][0] == str(executable)
    assert saved["twohdmc_provenance"] == "unknown"


def test_timeout_preserves_logs_and_writes_timed_out_manifest(tmp_path: Path) -> None:
    executable = tmp_path / "slow-evaluator"
    executable.write_text(
        "#!/bin/sh\n"
        "printf 'partial stdout\\n'\n"
        "printf 'partial stderr\\n' >&2\n"
        "sleep 2\n",
        encoding="utf-8",
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)

    with pytest.raises(subprocess.TimeoutExpired):
        run_lambda1_v2(
            executable=executable, rows=rows(), outdir=tmp_path / "out",
            campaign="pilot", run_name="timeout", repo_root=tmp_path, timeout_s=0.1,
        )

    run_dir = tmp_path / "out/campaign=pilot/timeout"
    manifest = json.loads((run_dir / "run_manifest.json").read_text())
    assert manifest["status"] == "timed_out"
    assert manifest["returncode"] is None
    assert manifest["timeout_s"] == 0.1
    assert "partial stdout" in (run_dir / "stdout.log").read_text()
    assert "partial stderr" in (run_dir / "stderr.log").read_text()
    assert "output_sha256" not in manifest


def test_new_lambda1_production_rejects_an_unlabelled_historical_mass(tmp_path: Path) -> None:
    historical_rows = cartesian_rows(
        fixed=Lambda1Fixed(125.09, 300.0, 310.0, 0.995, 0.1, 0.0),
        mH_values=[130.0], lambda1_values=[1.0], tan_beta_values=[50.0],
    )
    with pytest.raises(RuntimeError, match="historical replay requires"):
        run_lambda1_v2(
            executable=tmp_path / "unused",
            rows=historical_rows,
            outdir=tmp_path / "out",
            campaign="new-production",
            run_name="historical-override",
            repo_root=tmp_path,
        )
    assert not (tmp_path / "out/campaign=new-production/historical-override/input_v2.csv").exists()
