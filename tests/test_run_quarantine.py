from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
import sys
import threading
import time
from typing import cast

import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from scripts.run_quarantine import (
    OUTPUT_FIELDNAMES,
    build_command,
    build_manifest,
    compute_manifest_sha256,
    list_ready_campaigns,
    main,
    normalize_campaign_name,
    read_csv_rows,
    run_manifest_row,
)


def make_input_row(
    *,
    readiness_status: str = "READY",
    source_campaign: str = "campaign=target",
    source_run: str = "run-a",
    filename: str = "a.csv",
    source_row_index: str = "0",
    source_row_sha256: str = "sha-a0",
    run_id: str = "run-id-a",
    tb_folder: str = "tb-000",
    m_phi: str = "125.0",
    mA: str = "300.0",
    alpha: str = "0.2",
    beta: str = "0.6",
    sin_ba: str = "0.999",
    tan_beta: str = "50.0",
    lambda6: str = "0.1",
    lambda7: str = "0.0",
    m12: str = "1000.0",
    explorer_guess: str = "lam1_explorer",
    lam1: str = "0.25",
) -> dict[str, str]:
    return {
        "readiness_status": readiness_status,
        "source_campaign": source_campaign,
        "source_run": source_run,
        "filename": filename,
        "source_row_index": source_row_index,
        "source_row_sha256": source_row_sha256,
        "run_id": run_id,
        "tb_folder": tb_folder,
        "m_phi": m_phi,
        "mA": mA,
        "alpha": alpha,
        "beta": beta,
        "sin_ba": sin_ba,
        "tan_beta": tan_beta,
        "lambda6": lambda6,
        "lambda7": lambda7,
        "m12": m12,
        "explorer_guess": explorer_guess,
        "lam1": lam1,
    }


def write_input_csv(path: Path, rows: list[dict[str, str]]) -> None:
    base_fieldnames = [
        "readiness_status",
        "source_campaign",
        "source_run",
        "filename",
        "source_row_index",
        "source_row_sha256",
        "run_id",
        "tb_folder",
        "m_phi",
        "mA",
        "alpha",
        "beta",
        "sin_ba",
        "tan_beta",
        "lambda6",
        "lambda7",
        "m12",
        "explorer_guess",
        "lam1",
    ]
    extras = [key for row in rows for key in row if key not in base_fieldnames]
    fieldnames = [*base_fieldnames, *list(dict.fromkeys(extras))]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_output_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), [dict(row) for row in reader]


def logical_output_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    ignored_fields = {
        "worker_subset",
        "started_at_epoch",
        "finished_at_epoch",
        "elapsed_sec",
    }
    return [
        {field: value for field, value in row.items() if field not in ignored_fields}
        for row in rows
    ]


def test_build_manifest_filters_ready_campaign_sorts_limits_files_subsets_and_keeps_physics() -> None:
    manifest = build_manifest(
        rows=[
            make_input_row(source_campaign="campaign=other", source_run="run-b", filename="b.csv", source_row_index="2", source_row_sha256="zzz"),
            make_input_row(source_run="run-b", filename="b.csv", source_row_index="5", source_row_sha256="bbb", m_phi="126.0", alpha="0.3"),
            make_input_row(readiness_status="BLOCKED", source_run="run-a", filename="a.csv", source_row_index="1", source_row_sha256="aaa"),
            make_input_row(source_run="run-a", filename="a.csv", source_row_index="3", source_row_sha256="ccc", m_phi="127.0", alpha="0.4"),
            make_input_row(source_run="run-a", filename="a.csv", source_row_index="1", source_row_sha256="aaa"),
            make_input_row(source_run="run-c", filename="c.csv", source_row_index="9", source_row_sha256="ddd"),
        ],
        campaign="campaign=target",
        max_files=2,
    )

    subset_manifest = [row for index, row in enumerate(manifest) if index % 2 == 1]

    assert [
        (
            row["source_campaign"],
            row["source_run"],
            row["filename"],
            row["source_row_index"],
            row["source_row_sha256"],
        )
        for row in subset_manifest
    ] == [
        ("campaign=target", "run-a", "a.csv", "3", "ccc"),
    ]
    assert subset_manifest[0]["m_phi"] == "127.0"
    assert subset_manifest[0]["alpha"] == "0.4"
    assert subset_manifest[0]["beta"] == "0.6"
    assert subset_manifest[0]["m12"] == "1000.0"
    assert subset_manifest[0]["explorer_guess"] == "lam1_explorer"


def test_build_manifest_deduplicates_and_hash_is_stable_across_input_order() -> None:
    logical_rows = [
        make_input_row(source_run="run-b", filename="b.csv", source_row_index="2", source_row_sha256="sha-b2"),
        make_input_row(source_run="run-a", filename="a.csv", source_row_index="1", source_row_sha256="sha-a1"),
        make_input_row(source_run="run-c", filename="c.csv", source_row_index="9", source_row_sha256="sha-c9"),
        make_input_row(source_run="run-a", filename="a.csv", source_row_index="3", source_row_sha256="sha-a3"),
    ]
    rows_a = [
        {**logical_rows[0], "source_quarantine_path": "/lake-a/quarantine/run-b/b.csv"},
        {**logical_rows[1], "source_quarantine_path": "/lake-a/quarantine/run-a/a.csv"},
        {**logical_rows[2], "source_quarantine_path": "/lake-a/quarantine/run-c/c.csv"},
        {**logical_rows[1], "source_quarantine_path": "/lake-a/duplicate-prefix/run-a/a.csv"},
        {**logical_rows[3], "source_quarantine_path": "/lake-a/quarantine/run-a/a.csv"},
    ]
    rows_b = [
        {**logical_rows[3], "source_quarantine_path": "/other-root/run-a/a.csv"},
        {**logical_rows[1], "source_quarantine_path": "/other-root/run-a/a.csv"},
        {**logical_rows[2], "source_quarantine_path": "/other-root/run-c/c.csv"},
        {**logical_rows[0], "source_quarantine_path": "/other-root/run-b/b.csv"},
        {**logical_rows[1], "source_quarantine_path": "/other-root/dupe/run-a/a.csv"},
    ]

    manifest_a = build_manifest(rows=rows_a, campaign="campaign=target", max_files=None)
    manifest_b = build_manifest(rows=rows_b, campaign="campaign=target", max_files=None)

    assert manifest_a == manifest_b
    assert [
        (
            row["source_campaign"],
            row["source_run"],
            row["filename"],
            row["source_row_index"],
            row["source_row_sha256"],
        )
        for row in manifest_a
    ] == [
        ("campaign=target", "run-a", "a.csv", "1", "sha-a1"),
        ("campaign=target", "run-a", "a.csv", "3", "sha-a3"),
        ("campaign=target", "run-b", "b.csv", "2", "sha-b2"),
        ("campaign=target", "run-c", "c.csv", "9", "sha-c9"),
    ]
    expected_hash = hashlib.sha256(
        json.dumps(manifest_a, ensure_ascii=False, separators=(",", ":"), sort_keys=True).encode("utf-8")
    ).hexdigest()
    assert compute_manifest_sha256(manifest_a) == expected_hash
    assert compute_manifest_sha256(manifest_b) == expected_hash


def test_build_manifest_keeps_rows_distinct_by_run_id_and_tb_folder_but_limits_by_filename() -> None:
    manifest = build_manifest(
        rows=[
            make_input_row(run_id="run-id-a", tb_folder="tb-001", source_row_index="1", source_row_sha256="sha-a1"),
            make_input_row(run_id="run-id-a", tb_folder="tb-002", source_row_index="1", source_row_sha256="sha-a1"),
            make_input_row(run_id="run-id-b", tb_folder="tb-001", source_row_index="1", source_row_sha256="sha-a1"),
            make_input_row(source_run="run-b", filename="b.csv", source_row_index="2", source_row_sha256="sha-b2"),
        ],
        campaign="campaign=target",
        max_files=1,
    )

    assert [
        (
            row["source_campaign"],
            row["source_run"],
            row["filename"],
            row["source_row_index"],
            row["source_row_sha256"],
            row["run_id"],
            row["tb_folder"],
        )
        for row in manifest
    ] == [
        ("campaign=target", "run-a", "a.csv", "1", "sha-a1", "run-id-a", "tb-001"),
        ("campaign=target", "run-a", "a.csv", "1", "sha-a1", "run-id-a", "tb-002"),
        ("campaign=target", "run-a", "a.csv", "1", "sha-a1", "run-id-b", "tb-001"),
    ]


def test_compute_manifest_sha256_ignores_dictionary_insertion_order() -> None:
    ordered_row = build_manifest(
        rows=[make_input_row(source_row_index="1", source_row_sha256="sha-a1")],
        campaign="campaign=target",
        max_files=None,
    )[0]
    reversed_row = {key: ordered_row[key] for key in reversed(list(ordered_row.keys()))}

    assert compute_manifest_sha256([ordered_row]) == compute_manifest_sha256([reversed_row])


def test_build_manifest_allows_identical_duplicates_but_rejects_conflicts() -> None:
    identical_rows = [
        make_input_row(source_row_index="1", source_row_sha256="sha-a1", m_phi="125.0", lam1="0.25"),
        make_input_row(source_row_index="1", source_row_sha256="sha-a1", m_phi="125.0", lam1="0.25"),
    ]

    manifest = build_manifest(rows=identical_rows, campaign="campaign=target", max_files=None)

    assert len(manifest) == 1
    assert manifest[0]["source_row_sha256"] == "sha-a1"
    assert manifest[0]["m_phi"] == "125.0"
    assert manifest[0]["lam1"] == "0.25"

    conflicting_rows = [
        make_input_row(source_row_index="1", source_row_sha256="sha-a1", m_phi="125.0", lam1="0.25"),
        make_input_row(source_row_index="1", source_row_sha256="sha-a1", m_phi="126.0", lam1="0.25"),
    ]

    with pytest.raises(ValueError, match=r"campaign=target.*run-a.*a\.csv.*1.*sha-a1"):
        _ = build_manifest(rows=conflicting_rows, campaign="campaign=target", max_files=None)


def test_build_manifest_skips_malformed_unrelated_rows() -> None:
    manifest = build_manifest(
        rows=[
            make_input_row(source_campaign="broken-campaign-path", source_run="run-b", filename="b.csv", source_row_index="2"),
            make_input_row(source_campaign="campaign=target", source_run="run-a", filename="a.csv", source_row_index="1", source_row_sha256="sha-a1"),
        ],
        campaign="campaign=target",
        max_files=None,
    )

    assert [
        (row["source_campaign"], row["source_run"], row["filename"], row["source_row_index"])
        for row in manifest
    ] == [
        ("campaign=target", "run-a", "a.csv", "1"),
    ]


def test_list_ready_campaigns_skips_malformed_rows() -> None:
    campaigns = list_ready_campaigns(
        [
            make_input_row(source_campaign="broken-campaign-path", source_row_index="0"),
            make_input_row(source_campaign="/lake-a/scan/campaign=beta/run-1", source_row_index="1", source_row_sha256="sha-a1"),
            make_input_row(source_campaign="/lake-b/scan/campaign=alpha/run-2", source_row_index="2", source_row_sha256="sha-a2"),
        ]
    )

    assert campaigns == ["campaign=alpha", "campaign=beta"]


def test_normalize_campaign_name_accepts_windows_style_paths() -> None:
    assert normalize_campaign_name(r"C:\lake\campaign=alpha\run-1") == "campaign=alpha"


def test_normalize_campaign_name_rejects_paths_missing_campaign_token() -> None:
    with pytest.raises(ValueError, match="campaign="):
        _ = normalize_campaign_name("/lake/no-campaign-token/run-1")


def test_compute_manifest_sha256_changes_when_physics_payload_changes() -> None:
    manifest_a = build_manifest(
        rows=[make_input_row(source_row_index="1", source_row_sha256="sha-a1", m_phi="125.0", lam1="0.25")],
        campaign="campaign=target",
        max_files=None,
    )
    manifest_b = build_manifest(
        rows=[make_input_row(source_row_index="1", source_row_sha256="sha-a1", m_phi="126.0", lam1="0.25")],
        campaign="campaign=target",
        max_files=None,
    )

    assert [manifest_a[0][field] for field in ("source_campaign", "source_run", "filename", "source_row_index", "source_row_sha256")] == [
        manifest_b[0][field] for field in ("source_campaign", "source_run", "filename", "source_row_index", "source_row_sha256")
    ]
    assert compute_manifest_sha256(manifest_a) != compute_manifest_sha256(manifest_b)


def test_build_command_uses_physical_calculator_arguments() -> None:
    row = make_input_row(
        m_phi="130.0",
        mA="310.0",
        sin_ba="0.998",
        tan_beta="45.0",
        lambda6="0.2",
        lambda7="-0.1",
        lam1="1.5",
    )

    command = build_command(Path("PhysLam1Scan"), row, Path("/tmp/point.csv"))

    assert command == [
        "PhysLam1Scan",
        "130.0",
        "130.0",
        "1",
        "1.5",
        "1.5",
        "1",
        "310.0",
        "0.998",
        "45.0",
        "0.2",
        "-0.1",
        "/tmp/point.csv",
    ]


def test_build_command_rejects_unsupported_explorer() -> None:
    row = make_input_row(explorer_guess="tb_scan_explorer")

    with pytest.raises(ValueError, match="tb_scan_explorer"):
        _ = build_command(Path("PhysLam1Scan"), row, Path("/tmp/point.csv"))


@pytest.mark.parametrize(
    ("fieldname", "overrides"),
    [
        ("m_phi", {"m_phi": ""}),
        ("mA", {"mA": ""}),
        ("sin_ba", {"sin_ba": ""}),
        ("tan_beta", {"tan_beta": ""}),
        ("lambda6", {"lambda6": ""}),
        ("lambda7", {"lambda7": ""}),
        ("lam1", {"lam1": "", "legacy_lam1": "", "computed_lam1": "", "legacy_computed_lam1": ""}),
    ],
)
def test_build_manifest_rejects_selected_rows_missing_required_fields(
    fieldname: str,
    overrides: dict[str, str],
) -> None:
    row = {**make_input_row(source_row_index="1", source_row_sha256="sha-a1"), **overrides}

    with pytest.raises(ValueError, match=fieldname):
        _ = build_manifest(rows=[row], campaign="campaign=target", max_files=None)


@pytest.mark.parametrize(
    ("fallback_field", "fallback_value"),
    [
        ("legacy_lam1", "1.25"),
        ("computed_lam1", "1.5"),
        ("legacy_computed_lam1", "1.75"),
    ],
)
def test_build_command_uses_lam1_fallback_fields(fallback_field: str, fallback_value: str) -> None:
    row = make_input_row(lam1="")
    row[fallback_field] = fallback_value

    command = build_command(Path("PhysLam1Scan"), row, Path("/tmp/point.csv"))

    assert command[4] == fallback_value
    assert command[5] == fallback_value


@pytest.mark.parametrize("fallback_field", ["legacy_lam1", "computed_lam1", "legacy_computed_lam1"])
def test_build_manifest_accepts_selected_rows_with_lam1_fallback_fields(fallback_field: str) -> None:
    row = make_input_row(source_row_index="1", source_row_sha256="sha-a1", lam1="")
    row[fallback_field] = "0.75"

    manifest = build_manifest(rows=[row], campaign="campaign=target", max_files=None)

    assert len(manifest) == 1
    assert manifest[0][fallback_field] == "0.75"


def test_run_manifest_row_streams_stdout_and_parses_real_markers(monkeypatch: pytest.MonkeyPatch) -> None:
    progress_events: list[tuple[float, float]] = []
    sleep_calls = 0

    class FakeStream:
        def __init__(self, lines: list[str]) -> None:
            self._lines: list[str] = list(lines)
            self.drained: bool = False

        def readline(self) -> str:
            if self._lines:
                return self._lines.pop(0)
            self.drained = True
            return ""

        def close(self) -> None:
            return None

    class FakeProcess:
        def __init__(self) -> None:
            self.stdout: FakeStream = FakeStream(["warming up\n", "Total Attempts: 9\n", "TRIPLE_OK_POINTS 4\n"])
            self.stderr: FakeStream = FakeStream(["note: still healthy\n"])
            self.returncode: int | None = None

        def poll(self) -> int | None:
            if self.stdout.drained and self.stderr.drained:
                self.returncode = 0
            return self.returncode

    process = FakeProcess()
    time_values = iter([100.0, 100.2, 100.4, 100.6, 100.8, 101.0, 101.2, 101.4])
    perf_values = iter([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])

    def fake_sleep(_: float) -> None:
        nonlocal sleep_calls
        sleep_calls += 1
        if sleep_calls > 5:
            raise AssertionError("poll loop did not finish after pipes drained")

    def fake_popen(*args: object, **kwargs: object) -> FakeProcess:
        _ = args
        _ = kwargs
        return process

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen)
    monkeypatch.setattr("scripts.run_quarantine.time.time", lambda: next(time_values))
    monkeypatch.setattr("scripts.run_quarantine.time.perf_counter", lambda: next(perf_values))
    monkeypatch.setattr("scripts.run_quarantine.time.sleep", fake_sleep)

    output_row = run_manifest_row(
        Path("PhysLam1Scan"),
        make_input_row(),
        "0/1",
        progress_callback=lambda now_epoch, perf_elapsed_sec: progress_events.append((now_epoch, perf_elapsed_sec)),
        progress_every_sec=0,
    )

    assert output_row["run_status"] == "OK"
    assert output_row["calculator_returncode"] == "0"
    assert output_row["total_attempts"] == "9"
    assert output_row["triple_ok_points"] == "4"
    assert progress_events


def test_run_manifest_row_uses_real_subprocess_stub(tmp_path: Path) -> None:
    stub = tmp_path / "phys_stub.py"
    stub.write_text(
        "\n".join(
            [
                f"#!{sys.executable}",
                "import pathlib",
                "import sys",
                "import time",
                "pathlib.Path(sys.argv[-1]).write_text('result\\n', encoding='utf-8')",
                "print('boot', flush=True)",
                "time.sleep(0.05)",
                "print('Total Attempts: 4', flush=True)",
                "time.sleep(0.05)",
                "print('TRIPLE_OK_POINTS 2', flush=True)",
                "print('stub stderr note', file=sys.stderr, flush=True)",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    stub.chmod(0o755)

    progress_events: list[tuple[float, float]] = []
    output_row = run_manifest_row(
        stub,
        make_input_row(),
        "0/1",
        progress_callback=lambda now_epoch, perf_elapsed_sec: progress_events.append((now_epoch, perf_elapsed_sec)),
        progress_every_sec=0.01,
    )

    assert output_row["run_status"] == "OK"
    assert output_row["calculator_returncode"] == "0"
    assert output_row["total_attempts"] == "4"
    assert output_row["triple_ok_points"] == "2"
    assert progress_events


def test_run_manifest_row_pins_omp_threads_in_subprocess_env(monkeypatch: pytest.MonkeyPatch) -> None:
    popen_env: dict[str, str] = {}

    class FakeProcess:
        def __init__(self) -> None:
            self.stdout: io.StringIO = io.StringIO("Total Attempts: 1\nTRIPLE_OK_POINTS 2\n")
            self.stderr: io.StringIO = io.StringIO("")
            self.returncode: int | None = 0

        def poll(self) -> int:
            return 0

    def fake_popen(*args: object, **kwargs: object) -> FakeProcess:
        nonlocal popen_env
        _ = args
        popen_env = cast(dict[str, str], kwargs.get("env"))
        return FakeProcess()

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen)

    output_row = run_manifest_row(Path("PhysLam1Scan"), make_input_row(), "0/1")

    assert output_row["run_status"] == "OK"
    assert popen_env["OMP_NUM_THREADS"] == "1"


def test_run_manifest_row_returns_error_for_unsupported_explorer() -> None:
    output_row = run_manifest_row(
        Path("PhysLam1Scan"),
        make_input_row(explorer_guess="tb_scan_explorer"),
        "0/1",
    )

    assert output_row["run_status"] == "ERROR"
    assert output_row["error_code"] == "ValueError"
    assert "tb_scan_explorer" in cast(str, output_row["error_detail"])


def test_run_manifest_row_times_out_hung_child_process(monkeypatch: pytest.MonkeyPatch) -> None:
    class FakeStream:
        def readline(self) -> str:
            return ""

        def close(self) -> None:
            return None

    class FakeProcess:
        def __init__(self) -> None:
            self.stdout: FakeStream = FakeStream()
            self.stderr: FakeStream = FakeStream()
            self.returncode: int | None = None
            self.terminated = False

        def poll(self) -> int | None:
            return self.returncode

        def terminate(self) -> None:
            self.terminated = True
            self.returncode = -15

        def wait(self, timeout: float | None = None) -> int:
            _ = timeout
            return cast(int, self.returncode)

    process = FakeProcess()
    perf_values = iter([0.0, 0.1, 0.35, 0.45])
    time_values = iter([100.0, 100.1, 100.35, 100.45])

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", lambda *args, **kwargs: process)
    monkeypatch.setattr("scripts.run_quarantine.time.perf_counter", lambda: next(perf_values))
    monkeypatch.setattr("scripts.run_quarantine.time.time", lambda: next(time_values))
    monkeypatch.setattr("scripts.run_quarantine.time.sleep", lambda _: None)

    output_row = run_manifest_row(
        Path("PhysLam1Scan"),
        make_input_row(),
        "0/1",
        process_timeout_sec=0.3,
    )

    assert process.terminated is True
    assert output_row["run_status"] == "ERROR"
    assert output_row["error_code"] == "TimeoutError"
    assert "0.3s" in cast(str, output_row["error_detail"])


def test_run_manifest_row_does_not_deadlock_when_child_exits_but_pipe_stays_open(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    entered_blocking_read = threading.Event()
    release_pipe = threading.Event()

    class BlockingStream:
        def __init__(self, initial_line: str) -> None:
            self._initial_line = initial_line
            self._served_initial_line = False

        def readline(self) -> str:
            if not self._served_initial_line:
                self._served_initial_line = True
                return self._initial_line
            entered_blocking_read.set()
            release_pipe.wait(timeout=5.0)
            return ""

        def close(self) -> None:
            return None

    class FakeProcess:
        def __init__(self) -> None:
            self.stdout: BlockingStream = BlockingStream("Total Attempts: 1\nTRIPLE_OK_POINTS 2\n")
            self.stderr: BlockingStream = BlockingStream("")
            self.returncode: int | None = 0

        def poll(self) -> int:
            return 0

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", lambda *args, **kwargs: FakeProcess())
    monkeypatch.setattr("scripts.run_quarantine.STREAM_JOIN_TIMEOUT_SEC", 0.01)

    started = time.perf_counter()
    output_row = run_manifest_row(Path("PhysLam1Scan"), make_input_row(), "0/1")
    elapsed_sec = time.perf_counter() - started
    release_pipe.set()

    assert entered_blocking_read.is_set() is True
    assert elapsed_sec < 0.5
    assert output_row["run_status"] == "OK"
    assert output_row["calculator_returncode"] == "0"


def test_run_manifest_row_keeps_late_draining_output_after_child_exit(monkeypatch: pytest.MonkeyPatch) -> None:
    class LateStream:
        def __init__(self, chunks: list[tuple[float, str]]) -> None:
            self._chunks = list(chunks)

        def readline(self) -> str:
            if not self._chunks:
                return ""
            delay_sec, chunk = self._chunks.pop(0)
            if delay_sec > 0:
                time.sleep(delay_sec)
            return chunk

        def close(self) -> None:
            return None

    class FakeProcess:
        def __init__(self) -> None:
            self.stdout: LateStream = LateStream([(0.0, "Total Attempts: 11\n"), (0.02, "")])
            self.stderr: LateStream = LateStream([(0.02, "late stderr marker\n"), (0.0, "")])
            self.returncode: int | None = 7

        def poll(self) -> int:
            return 7

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", lambda *args, **kwargs: FakeProcess())
    monkeypatch.setattr("scripts.run_quarantine.STREAM_JOIN_TIMEOUT_SEC", 0.05)
    monkeypatch.setattr("scripts.run_quarantine.STREAM_JOIN_SLICE_SEC", 0.005)
    monkeypatch.setattr("scripts.run_quarantine.STREAM_DRAIN_QUIET_SEC", 0.005)

    output_row = run_manifest_row(Path("PhysLam1Scan"), make_input_row(), "0/1")

    assert output_row["run_status"] == "ERROR"
    assert output_row["calculator_returncode"] == "7"
    assert output_row["total_attempts"] == "11"
    assert output_row["error_detail"] == "late stderr marker"


def test_main_runs_jobs_appends_results_and_writes_status_json(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_row_index="0", source_row_sha256="sha-a0", m_phi="130.0", lam1="0.1"),
            make_input_row(source_row_index="1", source_row_sha256="sha-a1", m_phi="131.0", lam1="0.2"),
            make_input_row(source_campaign="campaign=other", source_run="run-z", filename="z.csv", source_row_index="9", source_row_sha256="sha-z9"),
        ],
    )
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDNAMES)
        writer.writeheader()
        writer.writerow(
            {
                "source_campaign": "campaign=existing",
                "source_run": "run-old",
                "filename": "old.csv",
                "source_row_index": "99",
                "source_row_sha256": "old-sha",
                "run_status": "OK",
                "error_code": "",
                "error_detail": "",
                "worker_subset": "0/1",
                "started_at_epoch": "1.0",
                "finished_at_epoch": "2.0",
                "elapsed_sec": "1.0",
                "calculator_returncode": "0",
                "total_attempts": "1",
                "triple_ok_points": "5",
            }
        )

    call_log: list[list[str]] = []

    class FakeProcess:
        def __init__(self, stdout_text: str) -> None:
            self.stdout: io.StringIO = io.StringIO(stdout_text)
            self.stderr: io.StringIO = io.StringIO("")
            self.returncode: int | None = 0

        def poll(self) -> int:
            return 0

    def fake_popen(cmd: list[str], **_: object) -> FakeProcess:
        call_log.append(cmd)
        return FakeProcess("Total Attempts: 1\nTRIPLE_OK_POINTS 7\n")

    time_values = iter([100.0, 101.0, 110.0, 111.5, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0])
    perf_values = iter([0.0, 1.0, 1.0, 2.5, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0])

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen)
    monkeypatch.setattr("scripts.run_quarantine.time.time", lambda: next(time_values))
    monkeypatch.setattr("scripts.run_quarantine.time.perf_counter", lambda: next(perf_values))

    expected_manifest_hash = compute_manifest_sha256(
        build_manifest(rows=read_csv_rows(input_csv), campaign="campaign=target", max_files=None)
    )

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
            "--progress-every-sec",
            "0.001",
        ]
    )

    assert exit_code == 0
    assert len(call_log) == 2
    assert call_log[0][:12] == [
        "PhysLam1Scan",
        "130.0",
        "130.0",
        "1",
        "0.1",
        "0.1",
        "1",
        "300.0",
        "0.999",
        "50.0",
        "0.1",
        "0.0",
    ]
    assert call_log[1][:12] == [
        "PhysLam1Scan",
        "131.0",
        "131.0",
        "1",
        "0.2",
        "0.2",
        "1",
        "300.0",
        "0.999",
        "50.0",
        "0.1",
        "0.0",
    ]

    fieldnames, rows = read_output_csv(output_csv)
    assert fieldnames == OUTPUT_FIELDNAMES
    assert len(rows) == 3
    assert rows[1]["run_status"] == "OK"
    assert rows[1]["run_id"] == "run-id-a"
    assert rows[1]["tb_folder"] == "tb-000"
    assert rows[1]["calculator_returncode"] == "0"
    assert rows[1]["total_attempts"] == "1"
    assert rows[1]["triple_ok_points"] == "7"
    assert rows[2]["source_row_index"] == "1"

    status = cast(dict[str, object], json.loads(status_json.read_text(encoding="utf-8")))
    assert status["worker_subset"] == "0/1"
    assert status["processed_rows"] == 2
    assert status["ok_rows"] == 2
    assert status["error_rows"] == 0
    assert status["done"] is True
    assert status["alive"] is False
    assert status["total_rows"] == 2
    assert status["subset_total_rows"] == 2
    assert status["campaign_total_rows"] == 2
    assert status["manifest_sha256"] == expected_manifest_hash
    assert status["rows_per_sec"] == pytest.approx(0.2)

    stdout = capsys.readouterr().out
    assert "alive=True" in stdout
    assert "alive=False" in stdout
    assert "processed=2/2" in stdout
    assert "campaign_total=2" in stdout
    assert "rows_per_sec=0.2" in stdout
    assert f"manifest_sha256={expected_manifest_hash}" in stdout


def test_main_reports_live_rows_per_sec_for_inflight_row(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    write_input_csv(
        input_csv,
        [make_input_row(source_row_index="0", source_row_sha256="sha-a0")],
    )

    class FakeProcess:
        def __init__(self) -> None:
            self.stdout: io.StringIO = io.StringIO("Total Attempts: 1\nTRIPLE_OK_POINTS 1\n")
            self.stderr: io.StringIO = io.StringIO("")
            self.returncode: int | None = None
            self.poll_calls = 0

        def poll(self) -> int | None:
            self.poll_calls += 1
            if self.poll_calls >= 4:
                self.returncode = 0
            return self.returncode

    time_values = iter([100.0, 101.0, 102.0, 103.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0])
    perf_values = iter([0.0, 0.4, 1.2, 2.2, 2.4, 3.0, 3.4, 3.8, 4.2, 4.6])

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", lambda *args, **kwargs: FakeProcess())
    monkeypatch.setattr("scripts.run_quarantine.time.time", lambda: next(time_values))
    monkeypatch.setattr("scripts.run_quarantine.time.perf_counter", lambda: next(perf_values))
    monkeypatch.setattr("scripts.run_quarantine.time.sleep", lambda _: None)

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--binary",
            "PhysLam1Scan",
            "--progress-every-sec",
            "1.5",
        ]
    )

    assert exit_code == 0
    stdout = capsys.readouterr().out
    assert "rows_per_sec=0.2" in stdout
    assert "campaign_total=1" in stdout


def test_main_continues_after_nonzero_exit_and_parses_attempt_count(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_row_index="0", source_row_sha256="sha-a0", lam1="0.3"),
            make_input_row(source_row_index="1", source_row_sha256="sha-a1", lam1="0.4"),
        ],
    )

    responses = iter([
        (0, "Total Attempts: 3\nTRIPLE_OK_POINTS 3\n", ""),
        (9, "Total Attempts: 9\n", "calculator exploded\n"),
    ])

    class FakeProcess:
        def __init__(self, returncode: int, stdout_text: str, stderr_text: str) -> None:
            self.stdout: io.StringIO = io.StringIO(stdout_text)
            self.stderr: io.StringIO = io.StringIO(stderr_text)
            self.returncode: int | None = returncode

        def poll(self) -> int:
            return cast(int, self.returncode)

    def fake_popen(*args: object, **kwargs: object) -> FakeProcess:
        _ = args
        _ = kwargs
        returncode, stdout_text, stderr_text = next(responses)
        return FakeProcess(returncode, stdout_text, stderr_text)

    time_values = iter([200.0, 201.0, 202.0, 203.0, 204.0, 205.0, 206.0, 207.0, 208.0, 209.0])
    perf_values = iter([0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0])

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen)
    monkeypatch.setattr("scripts.run_quarantine.time.time", lambda: next(time_values))
    monkeypatch.setattr("scripts.run_quarantine.time.perf_counter", lambda: next(perf_values))

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--binary",
            "PhysLam1Scan",
            "--progress-every-sec",
            "0.001",
        ]
    )

    assert exit_code == 0

    _, rows = read_output_csv(output_csv)
    assert len(rows) == 2
    assert rows[0]["run_status"] == "OK"
    assert rows[0]["triple_ok_points"] == "3"
    assert rows[0]["total_attempts"] == "3"
    assert rows[1]["run_status"] == "ERROR"
    assert rows[1]["error_code"] == "NONZERO_EXIT"
    assert rows[1]["error_detail"] == "calculator exploded"
    assert rows[1]["calculator_returncode"] == "9"
    assert rows[1]["total_attempts"] == "9"
    assert rows[1]["triple_ok_points"] == ""

    stdout = capsys.readouterr().out
    assert "errors=1" in stdout


def test_main_returns_nonzero_when_all_selected_rows_fail_to_execute(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_row_index="0", source_row_sha256="sha-a0"),
            make_input_row(source_row_index="1", source_row_sha256="sha-a1"),
        ],
    )

    def fake_popen(*args: object, **kwargs: object) -> None:
        _ = args
        _ = kwargs
        raise FileNotFoundError("missing PhysLam1Scan")

    time_values = iter([300.0, 301.0, 302.0, 303.0, 304.0, 305.0, 306.0, 307.0, 308.0, 309.0])
    perf_values = iter([0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen)
    monkeypatch.setattr("scripts.run_quarantine.time.time", lambda: next(time_values))
    monkeypatch.setattr("scripts.run_quarantine.time.perf_counter", lambda: next(perf_values))

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
        ]
    )

    assert exit_code == 1
    _, rows = read_output_csv(output_csv)
    assert len(rows) == 2
    assert [row["run_status"] for row in rows] == ["ERROR", "ERROR"]
    assert [row["error_code"] for row in rows] == ["FileNotFoundError", "FileNotFoundError"]
    status = cast(dict[str, object], json.loads(status_json.read_text(encoding="utf-8")))
    assert status["ok_rows"] == 0
    assert status["error_rows"] == 2
    assert status["done"] is True
    assert status["subset_total_rows"] == 2
    assert status["campaign_total_rows"] == 2
    assert status["rows_per_sec"] == pytest.approx(0.2)
    assert "errors=2" in capsys.readouterr().out


def test_main_fails_fast_when_campaign_has_no_ready_rows(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(
        input_csv,
        [make_input_row(source_campaign="campaign=other")],
    )

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=missing",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
        ]
    )

    assert exit_code != 0
    assert not output_csv.exists()
    assert not status_json.exists()
    assert "0 rows" in capsys.readouterr().err


def test_main_returns_success_for_empty_subset_of_existing_campaign(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(
        input_csv,
        [make_input_row(source_row_index="0")],
    )

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
            "--subset",
            "1/2",
        ]
    )

    assert exit_code == 0
    assert not output_csv.exists()
    status = cast(dict[str, object], json.loads(status_json.read_text(encoding="utf-8")))
    assert status["total_rows"] == 0
    assert status["subset_total_rows"] == 0
    assert status["campaign_total_rows"] == 1
    assert status["processed_rows"] == 0
    assert status["done"] is True
    assert status["alive"] is False
    stdout = capsys.readouterr().out
    assert "received 0 rows" in stdout


def test_main_treats_max_files_zero_as_empty_success(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_row_index="0", source_row_sha256="sha-a0"),
            make_input_row(source_row_index="1", source_row_sha256="sha-a1"),
        ],
    )

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
            "--max-files",
            "0",
        ]
    )

    assert exit_code == 0
    assert not output_csv.exists()
    status = cast(dict[str, object], json.loads(status_json.read_text(encoding="utf-8")))
    assert status["total_rows"] == 0
    assert status["subset_total_rows"] == 0
    assert status["campaign_total_rows"] == 2
    assert status["processed_rows"] == 0
    assert status["done"] is True
    assert status["alive"] is False
    captured = capsys.readouterr()
    assert "--max-files 0 selected no files" in captured.out
    assert captured.err == ""


def test_main_rejects_missing_campaign_when_max_files_zero(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(
        input_csv,
        [make_input_row(source_campaign="campaign=other")],
    )

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=missing",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
            "--max-files",
            "0",
        ]
    )

    assert exit_code == 1
    assert not output_csv.exists()
    assert not status_json.exists()
    assert "0 rows" in capsys.readouterr().err


def test_main_subset_isolated_from_unselected_unsupported_explorer_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_run="run-a", filename="a.csv", source_row_index="0", source_row_sha256="sha-a0", lam1="0.1"),
            make_input_row(
                source_run="run-b",
                filename="b.csv",
                source_row_index="1",
                source_row_sha256="sha-b1",
                run_id="run-id-b",
                tb_folder="tb-001",
                explorer_guess="tb_scan_explorer",
                lam1="0.2",
            ),
        ],
    )

    class FakeProcess:
        def __init__(self) -> None:
            self.stdout: io.StringIO = io.StringIO("Total Attempts: 1\nTRIPLE_OK_POINTS 1\n")
            self.stderr: io.StringIO = io.StringIO("")
            self.returncode: int | None = 0

        def poll(self) -> int:
            return 0

    popen_calls = 0

    def fake_popen(*args: object, **kwargs: object) -> FakeProcess:
        nonlocal popen_calls
        popen_calls += 1
        _ = args
        _ = kwargs
        return FakeProcess()

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen)

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
            "--subset",
            "0/2",
        ]
    )

    assert exit_code == 0
    assert popen_calls == 1
    _, rows = read_output_csv(output_csv)
    assert len(rows) == 1
    assert rows[0]["source_row_sha256"] == "sha-a0"
    assert cast(dict[str, object], json.loads(status_json.read_text(encoding="utf-8")))["manifest_sha256"]
    captured = capsys.readouterr()
    assert captured.err == ""


def test_main_fails_when_selected_subset_contains_unsupported_explorer(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_row_index="0", source_row_sha256="sha-a0", explorer_guess="lam1_explorer"),
            make_input_row(
                source_run="run-b",
                filename="b.csv",
                source_row_index="1",
                source_row_sha256="sha-a1",
                run_id="run-id-b",
                tb_folder="tb-001",
                explorer_guess="tb_scan_explorer",
            ),
        ],
    )
    popen_calls = 0

    def fake_popen(*args: object, **kwargs: object) -> None:
        nonlocal popen_calls
        popen_calls += 1
        _ = args
        _ = kwargs
        raise AssertionError("Popen must not run when selected subset validation fails")

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen)

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
            "--subset",
            "1/2",
        ]
    )

    assert exit_code == 1
    assert popen_calls == 0
    assert not output_csv.exists()
    assert not status_json.exists()
    stderr = capsys.readouterr().err
    assert "Unsupported explorer manifest" in stderr
    assert "tb_scan_explorer" in stderr


def test_main_skips_unselected_unsupported_explorer_rows_when_max_files_limits_subset(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_row_index="0", source_row_sha256="sha-a0", explorer_guess="lam1_explorer"),
            make_input_row(
                source_run="run-b",
                filename="b.csv",
                source_row_index="1",
                source_row_sha256="sha-a1",
                run_id="run-id-b",
                tb_folder="tb-001",
                explorer_guess="tb_scan_explorer",
            ),
        ],
    )
    popen_calls = 0

    class FakeProcess:
        def __init__(self) -> None:
            self.stdout: io.StringIO = io.StringIO("Total Attempts: 1\nTRIPLE_OK_POINTS 1\n")
            self.stderr: io.StringIO = io.StringIO("")
            self.returncode: int | None = 0

        def poll(self) -> int:
            return 0

    def fake_popen(*args: object, **kwargs: object) -> FakeProcess:
        nonlocal popen_calls
        popen_calls += 1
        _ = args
        _ = kwargs
        return FakeProcess()

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen)

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
            "--max-files",
            "1",
        ]
    )

    assert exit_code == 0
    assert popen_calls == 1
    assert output_csv.exists()
    assert status_json.exists()
    captured = capsys.readouterr()
    assert captured.err == ""


def test_main_skips_unselected_conflicting_duplicate_rows_when_max_files_limits_subset(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_run="run-a", filename="a.csv", source_row_index="0", source_row_sha256="sha-a0", m_phi="125.0"),
            make_input_row(source_run="run-z", filename="z.csv", source_row_index="9", source_row_sha256="sha-z9", m_phi="126.0"),
            make_input_row(source_run="run-z", filename="z.csv", source_row_index="9", source_row_sha256="sha-z9", m_phi="127.0"),
        ],
    )

    popen_calls = 0

    class FakeProcess:
        def __init__(self) -> None:
            self.stdout: io.StringIO = io.StringIO("Total Attempts: 1\nTRIPLE_OK_POINTS 1\n")
            self.stderr: io.StringIO = io.StringIO("")
            self.returncode: int | None = 0

        def poll(self) -> int:
            return 0

    def fake_popen(*args: object, **kwargs: object) -> FakeProcess:
        nonlocal popen_calls
        popen_calls += 1
        _ = args
        _ = kwargs
        return FakeProcess()

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen)

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
            "--max-files",
            "1",
        ]
    )

    assert exit_code == 0
    assert popen_calls == 1
    assert output_csv.exists()
    assert status_json.exists()
    captured = capsys.readouterr()
    assert captured.err == ""



def test_main_fails_cleanly_for_conflicting_duplicate_manifest_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_row_index="1", source_row_sha256="sha-a1", m_phi="125.0", lam1="0.25"),
            make_input_row(source_row_index="1", source_row_sha256="sha-a1", m_phi="126.0", lam1="0.25"),
        ],
    )
    popen_calls = 0
    append_calls = 0
    status_calls = 0

    def fake_popen(*args: object, **kwargs: object) -> None:
        nonlocal popen_calls
        popen_calls += 1
        _ = args
        _ = kwargs
        raise AssertionError("Popen must not run for conflicting duplicates")

    def fake_append_output_row(*args: object, **kwargs: object) -> None:
        nonlocal append_calls
        append_calls += 1
        _ = args
        _ = kwargs
        raise AssertionError("output must not be written for conflicting duplicates")

    def fake_write_status_json(*args: object, **kwargs: object) -> None:
        nonlocal status_calls
        status_calls += 1
        _ = args
        _ = kwargs
        raise AssertionError("status must not be written for conflicting duplicates")

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen)
    monkeypatch.setattr("scripts.run_quarantine.append_output_row", fake_append_output_row)
    monkeypatch.setattr("scripts.run_quarantine.write_status_json", fake_write_status_json)

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
        ]
    )

    assert exit_code != 0
    assert popen_calls == 0
    assert append_calls == 0
    assert status_calls == 0
    assert not output_csv.exists()
    assert not status_json.exists()
    stderr = capsys.readouterr().err
    assert "conflicting duplicate manifest rows" in stderr
    assert "Traceback" not in stderr


def test_main_rejects_negative_max_files_without_traceback(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(input_csv, [make_input_row()])

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
            "--max-files",
            "-1",
        ]
    )

    assert exit_code != 0
    assert not output_csv.exists()
    assert not status_json.exists()
    stderr = capsys.readouterr().err
    assert "max-files" in stderr
    assert "non-negative" in stderr
    assert "Traceback" not in stderr


def test_main_rejects_nonpositive_progress_every_sec_without_traceback(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(input_csv, [make_input_row()])

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
            "--progress-every-sec",
            "0",
        ]
    )

    assert exit_code != 0
    assert not output_csv.exists()
    assert not status_json.exists()
    stderr = capsys.readouterr().err
    assert "progress-every-sec" in stderr
    assert "positive" in stderr
    assert "Traceback" not in stderr


def test_main_rejects_nonpositive_process_timeout_sec_without_traceback(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(input_csv, [make_input_row()])

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
            "--process-timeout-sec",
            "0",
        ]
    )

    assert exit_code != 0
    assert not output_csv.exists()
    assert not status_json.exists()
    stderr = capsys.readouterr().err
    assert "process-timeout-sec" in stderr
    assert "positive" in stderr
    assert "Traceback" not in stderr


def test_main_rejects_invalid_subset_without_traceback(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    status_json = tmp_path / "status.json"
    write_input_csv(input_csv, [make_input_row()])

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--status-json",
            str(status_json),
            "--binary",
            "PhysLam1Scan",
            "--subset",
            "2/2",
        ]
    )

    assert exit_code != 0
    assert not output_csv.exists()
    assert not status_json.exists()
    stderr = capsys.readouterr().err
    assert "subset" in stderr
    assert "0 <= K < N" in stderr
    assert "Traceback" not in stderr


def test_main_lists_unique_ready_campaigns_normalized_from_absolute_paths(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_campaign="/lake-a/scan/campaign=beta/run-1"),
            make_input_row(source_campaign="/lake-b/other/campaign=alpha/run-2", source_row_index="1", source_row_sha256="sha-a1"),
            make_input_row(source_campaign="/lake-c/other/campaign=beta/run-3", source_row_index="2", source_row_sha256="sha-a2"),
            make_input_row(
                source_campaign="/lake-d/other/campaign=gamma/run-4",
                source_row_index="3",
                source_row_sha256="sha-a3",
                readiness_status="BLOCKED",
            ),
        ],
    )

    exit_code = main([
        "--input-csv",
        str(input_csv),
        "--list-campaigns",
    ])

    assert exit_code == 0
    captured = capsys.readouterr()
    assert captured.out.splitlines() == ["1 campaign=beta", "2 campaign=alpha"]
    assert captured.err == ""


def test_main_list_campaigns_skips_malformed_rows(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_campaign="broken-campaign-path", source_row_index="0"),
            make_input_row(source_campaign="/lake-a/scan/campaign=beta/run-1", source_row_index="1", source_row_sha256="sha-a1"),
        ],
    )

    exit_code = main([
        "--input-csv",
        str(input_csv),
        "--list-campaigns",
    ])

    assert exit_code == 0
    captured = capsys.readouterr()
    assert captured.out.splitlines() == ["1 campaign=beta"]
    assert "Skipping malformed row:" in captured.err


def test_main_campaign_can_be_selected_by_numeric_index(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    write_input_csv(
        input_csv,
        [
            make_input_row(
                source_campaign="/lake-a/scan/campaign=older/run-1",
                source_row_index="0",
                source_row_sha256="sha-old",
            ),
            make_input_row(
                source_campaign="/lake-b/scan/campaign=newer/run-1",
                source_row_index="1",
                source_row_sha256="sha-new",
            ),
        ],
    )

    class FakeProcess:
        def __init__(self) -> None:
            self.stdout: io.StringIO = io.StringIO("Total Attempts: 1\nTRIPLE_OK_POINTS 1\n")
            self.stderr: io.StringIO = io.StringIO("")
            self.returncode: int | None = 0

        def poll(self) -> int:
            return 0

    def fake_popen(*args: object, **kwargs: object) -> FakeProcess:
        _ = args
        _ = kwargs
        return FakeProcess()

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen)

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "1",
            "--output-csv",
            str(output_csv),
            "--binary",
            "PhysLam1Scan",
            "--max-files",
            "1",
        ]
    )

    assert exit_code == 0
    assert output_csv.exists()
    header, rows = read_output_csv(output_csv)
    assert header == OUTPUT_FIELDNAMES
    assert {row["source_campaign"] for row in rows} == {"campaign=newer"}


def test_main_campaign_with_max_files_skips_unrelated_malformed_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_campaign="broken-campaign-path", source_run="run-z", filename="z.csv", source_row_index="9", source_row_sha256="sha-z9"),
            make_input_row(source_campaign="campaign=target", source_run="run-a", filename="a.csv", source_row_index="0", source_row_sha256="sha-a0"),
        ],
    )

    popen_calls = 0

    class FakeProcess:
        def __init__(self) -> None:
            self.stdout: io.StringIO = io.StringIO("Total Attempts: 1\nTRIPLE_OK_POINTS 1\n")
            self.stderr: io.StringIO = io.StringIO("")
            self.returncode: int | None = 0

        def poll(self) -> int:
            return 0

    def fake_popen(*args: object, **kwargs: object) -> FakeProcess:
        nonlocal popen_calls
        popen_calls += 1
        _ = args
        _ = kwargs
        return FakeProcess()

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen)

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--binary",
            "PhysLam1Scan",
            "--max-files",
            "1",
        ]
    )

    assert exit_code == 0
    assert popen_calls == 1
    assert output_csv.exists()
    captured = capsys.readouterr()
    assert "Skipping malformed row:" in captured.err


def test_main_uses_full_manifest_hash_for_all_subsets(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_a = tmp_path / "out-a.csv"
    output_b = tmp_path / "out-b.csv"
    status_a = tmp_path / "status-a.json"
    status_b = tmp_path / "status-b.json"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_row_index="0", source_row_sha256="sha-a0", lam1="0.1"),
            make_input_row(source_row_index="1", source_row_sha256="sha-a1", lam1="0.2"),
        ],
    )

    class FakeProcess:
        def __init__(self) -> None:
            self.stdout: io.StringIO = io.StringIO("Total Attempts: 1\nTRIPLE_OK_POINTS 1\n")
            self.stderr: io.StringIO = io.StringIO("")
            self.returncode: int | None = 0

        def poll(self) -> int:
            return 0

    def fake_popen_subset(*args: object, **kwargs: object) -> FakeProcess:
        _ = args
        _ = kwargs
        return FakeProcess()

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen_subset)

    full_manifest = build_manifest(rows=read_csv_rows(input_csv), campaign="campaign=target", max_files=None)
    expected_hash = compute_manifest_sha256(full_manifest)

    assert main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_a),
            "--status-json",
            str(status_a),
            "--binary",
            "PhysLam1Scan",
            "--subset",
            "0/2",
        ]
    ) == 0
    assert main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_b),
            "--status-json",
            str(status_b),
            "--binary",
            "PhysLam1Scan",
            "--subset",
            "1/2",
        ]
    ) == 0

    payload_a = cast(dict[str, object], json.loads(status_a.read_text(encoding="utf-8")))
    payload_b = cast(dict[str, object], json.loads(status_b.read_text(encoding="utf-8")))
    assert payload_a["manifest_sha256"] == expected_hash
    assert payload_b["manifest_sha256"] == expected_hash
    assert payload_a["campaign_total_rows"] == 2
    assert payload_b["campaign_total_rows"] == 2
    assert payload_a["subset_total_rows"] == 1
    assert payload_b["subset_total_rows"] == 1


def test_main_rejects_existing_output_with_mismatched_header(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    output_csv = tmp_path / "out.csv"
    write_input_csv(input_csv, [make_input_row(source_row_index="0", source_row_sha256="sha-a0")])
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["bad", "header"])
        writer.writeheader()
        writer.writerow({"bad": "1", "header": "2"})

    popen_calls = 0

    def fake_popen(*args: object, **kwargs: object) -> None:
        nonlocal popen_calls
        popen_calls += 1
        _ = args
        _ = kwargs
        raise AssertionError("Popen must not run when output schema is invalid")

    monkeypatch.setattr("scripts.run_quarantine.subprocess.Popen", fake_popen)

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--binary",
            "PhysLam1Scan",
        ]
    )

    assert exit_code == 1
    assert popen_calls == 0
    stderr = capsys.readouterr().err
    assert "existing output CSV header mismatch" in stderr
    _, rows = read_output_csv(output_csv)
    assert rows == [{"bad": "1", "header": "2"}]


def test_main_reports_missing_input_csv_cleanly(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "priority" / "missing.csv"
    output_csv = tmp_path / "out.csv"

    exit_code = main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_csv),
            "--binary",
            "PhysLam1Scan",
        ]
    )

    assert exit_code == 1
    stderr = capsys.readouterr().err
    assert f"Failed to read input CSV '{input_csv}'" in stderr
    assert "Traceback" not in stderr


def test_main_subset_outputs_append_to_full_run_exactly_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "priority" / "recompute_inputs_v1.csv"
    full_output = tmp_path / "full.csv"
    output_a = tmp_path / "out-a.csv"
    output_b = tmp_path / "out-b.csv"
    write_input_csv(
        input_csv,
        [
            make_input_row(source_run="run-a", filename="a.csv", source_row_index="0", source_row_sha256="sha-a0", lam1="0.1"),
            make_input_row(source_run="run-a", filename="a.csv", source_row_index="1", source_row_sha256="sha-a1", lam1="0.2"),
            make_input_row(source_run="run-b", filename="b.csv", source_row_index="2", source_row_sha256="sha-b2", lam1="0.3"),
            make_input_row(source_run="run-b", filename="b.csv", source_row_index="2", source_row_sha256="sha-b2", lam1="0.3"),
        ],
    )

    def fake_run_manifest_row(
        binary: Path,
        row: dict[str, str],
        worker_subset: str,
        progress_callback: object = None,
        progress_every_sec: float = 5.0,
        process_timeout_sec: float = 600.0,
    ) -> dict[str, str]:
        _ = binary
        _ = progress_callback
        _ = progress_every_sec
        _ = process_timeout_sec
        return {
            **row,
            "run_status": "OK",
            "error_code": "",
            "error_detail": "",
            "worker_subset": worker_subset,
            "started_at_epoch": "1.0",
            "finished_at_epoch": "2.0",
            "elapsed_sec": "1.0",
            "calculator_returncode": "0",
            "total_attempts": row["source_row_index"],
            "triple_ok_points": row["source_row_index"],
        }

    monkeypatch.setattr("scripts.run_quarantine.run_manifest_row", fake_run_manifest_row)

    assert main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(full_output),
            "--binary",
            "PhysLam1Scan",
        ]
    ) == 0
    assert main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_a),
            "--binary",
            "PhysLam1Scan",
            "--subset",
            "0/2",
        ]
    ) == 0
    assert main(
        [
            "--input-csv",
            str(input_csv),
            "--campaign",
            "campaign=target",
            "--output-csv",
            str(output_b),
            "--binary",
            "PhysLam1Scan",
            "--subset",
            "1/2",
        ]
    ) == 0

    _, full_rows = read_output_csv(full_output)
    _, rows_a = read_output_csv(output_a)
    _, rows_b = read_output_csv(output_b)

    combined_rows = logical_output_rows(rows_a) + logical_output_rows(rows_b)
    assert len(combined_rows) == len(logical_output_rows(full_rows))
    assert sorted(combined_rows, key=lambda row: (row["source_run"], row["filename"], row["source_row_index"])) == sorted(
        logical_output_rows(full_rows),
        key=lambda row: (row["source_run"], row["filename"], row["source_row_index"]),
    )
