import json

from dihiggs.app.orchestrator.resume import find_previous_csv, should_skip


def _write_meta(tb_dir, *, event, grid_signature):
    (tb_dir / "scan_meta.json").write_text(
        json.dumps({"event": event, "grid_signature": grid_signature}),
        encoding="utf-8",
    )


def test_done_signature_without_output_does_not_skip(tmp_path):
    output_csv = tmp_path / "tb_100" / "scan_tb_100.csv"

    skip, reason = should_skip(
        output_csv=output_csv,
        grid_sig="sig-current",
        force=False,
        done_signatures={"sig-current"},
    )

    assert skip is False
    assert reason == ""


def test_existing_csv_requires_done_task_metadata(tmp_path):
    tb_dir = tmp_path / "tb_100"
    tb_dir.mkdir()
    output_csv = tb_dir / "scan_tb_100.csv"
    output_csv.write_text("header\nrow\n", encoding="utf-8")
    _write_meta(tb_dir, event="fail", grid_signature="sig-current")

    skip, reason = should_skip(
        output_csv=output_csv,
        grid_sig="sig-current",
        force=False,
        done_signatures={"sig-current"},
    )

    assert skip is False
    assert reason == "stale_csv_meta_missing_or_not_done"


def test_cross_run_reuse_does_not_fallback_from_failed_task_to_run_manifest(tmp_path):
    fixed_dir = tmp_path / "fixed"
    previous_run = fixed_dir / "run_old"
    previous_tb = previous_run / "tb_100"
    previous_tb.mkdir(parents=True)
    previous_csv = previous_tb / "scan_tb_100.csv"
    previous_csv.write_text("header\npartial-row\n", encoding="utf-8")
    _write_meta(previous_tb, event="fail", grid_signature="sig-current")

    # A run-level signature must not override explicit failed per-task metadata.
    (previous_run / "run_manifest.json").write_text(
        json.dumps({"scan_grid": {"grid_signature": "sig-current"}}),
        encoding="utf-8",
    )

    current_run = fixed_dir / "run_current"
    current_run.mkdir(parents=True)

    found = find_previous_csv(
        fixed_dir=fixed_dir,
        current_run_dir=current_run,
        tb_tag="100",
        desired_grid_sig="sig-current",
    )

    assert found is None


def test_cross_run_reuse_accepts_explicit_done_task(tmp_path):
    fixed_dir = tmp_path / "fixed"
    previous_run = fixed_dir / "run_old"
    previous_tb = previous_run / "tb_100"
    previous_tb.mkdir(parents=True)
    previous_csv = previous_tb / "scan_tb_100.csv"
    previous_csv.write_text("header\nrow\n", encoding="utf-8")
    _write_meta(previous_tb, event="done", grid_signature="sig-current")

    current_run = fixed_dir / "run_current"
    current_run.mkdir(parents=True)

    found = find_previous_csv(
        fixed_dir=fixed_dir,
        current_run_dir=current_run,
        tb_tag="100",
        desired_grid_sig="sig-current",
    )

    assert found == previous_csv
