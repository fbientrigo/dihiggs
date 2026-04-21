#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import threading
import time
from collections import defaultdict, deque
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, cast

DEFAULT_SOURCE_ROOT = Path("/mnt/c/Users/Asus/cern_db/dihiggs_lake")
DEFAULT_BINARY = Path("2hdmc/CalcLam1Scan")
DEFAULT_PROCESS_TIMEOUT_SEC = 3600.0
POLL_INTERVAL_SEC = 0.1
STREAM_JOIN_TIMEOUT_SEC = 1.0
STREAM_JOIN_SLICE_SEC = 0.05
STREAM_DRAIN_QUIET_SEC = 0.05
RUNS_PER_SEC_WINDOW = 10
TOTAL_ATTEMPTS_PATTERN = re.compile(r"Total Attempts:\s*([0-9]+)")
TRIPLE_OK_POINTS_PATTERN = re.compile(r"TRIPLE_OK_POINTS\s+([0-9]+)")
RECOMPUTED_DIRNAME = "recomputed"
SCAN_GLOB = "tb_*/scan_tb_*.csv"
REQUIRED_NUMERIC_FIELDS = (
    "m_phi",
    "mA",
    "sin_ba",
    "tan_beta",
    "lambda6",
    "lambda7",
)
LAM1_FALLBACK_FIELDS = (
    "lam1",
    "legacy_lam1",
    "computed_lam1",
    "legacy_computed_lam1",
)
RECOMPUTE_INPUT_FIELDNAMES = [
    "point_id",
    "m_phi",
    "mA",
    "sin_ba",
    "tan_beta",
    "lambda6",
    "lambda7",
    "lam1",
]
STATUS_FIELDNAMES = [
    "source_campaign",
    "source_run_relative",
    "source_csv_relative",
    "destination_csv",
    "worker_subset",
    "source_row_count",
    "output_row_count",
    "run_status",
    "error_code",
    "error_detail",
    "started_at_epoch",
    "finished_at_epoch",
    "elapsed_sec",
    "calculator_returncode",
    "total_attempts",
    "triple_ok_points",
    "source_sha256",
    "destination_sha256",
]


@dataclass(frozen=True)
class ScanFileWorkItem:
    source_csv: Path
    source_root: Path

    @property
    def relative_csv(self) -> Path:
        return self.source_csv.relative_to(self.source_root)

    @property
    def source_campaign(self) -> str:
        return self.relative_csv.parts[0]

    @property
    def tb_dir(self) -> Path:
        return self.source_csv.parent

    @property
    def run_dir(self) -> Path:
        return self.tb_dir.parent

    @property
    def fixed_dir(self) -> Path:
        return self.run_dir.parent

    @property
    def relative_run_dir(self) -> Path:
        return self.run_dir.relative_to(self.source_root)


@dataclass(frozen=True)
class RunGroup:
    run_dir: Path
    source_root: Path
    scan_files: tuple[ScanFileWorkItem, ...]

    @property
    def relative_run_dir(self) -> Path:
        return self.run_dir.relative_to(self.source_root)

    @property
    def source_campaign(self) -> str:
        return self.relative_run_dir.parts[0]


def normalize_value(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def parse_subset(value: str) -> tuple[int, int]:
    match = re.fullmatch(r"(\d+)/(\d+)", value.strip())
    if match is None:
        raise argparse.ArgumentTypeError("subset must use K/N format")
    subset_index = int(match.group(1))
    subset_count = int(match.group(2))
    if subset_count <= 0:
        raise argparse.ArgumentTypeError("subset denominator must be positive")
    if subset_index < 0 or subset_index >= subset_count:
        raise argparse.ArgumentTypeError("subset numerator must satisfy 0 <= K < N")
    return (subset_index, subset_count)


def parse_non_negative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return parsed


def parse_positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def parse_positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def is_under(path: Path, parent: Path) -> bool:
    try:
        _ = path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def count_csv_data_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        _ = next(reader, None)
        return sum(1 for _ in reader)


def append_jsonl(path: Path, record: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(dict(record), ensure_ascii=False) + "\n")


def append_status_row(path: Path, row: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    need_header = not path.exists() or path.stat().st_size == 0
    if not need_header:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader, [])
        if header != STATUS_FIELDNAMES:
            raise ValueError(
                "existing status CSV header mismatch: "
                + f"expected {STATUS_FIELDNAMES!r}, got {header!r}"
            )
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDNAMES)
        if need_header:
            writer.writeheader()
        writer.writerow({field: normalize_value(row.get(field)) for field in STATUS_FIELDNAMES})


def write_json_atomic(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    tmp_path.write_text(json.dumps(dict(payload), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    tmp_path.replace(path)


def list_campaign_dirs(source_root: Path) -> list[Path]:
    campaign_dirs: list[Path] = []
    for entry in source_root.iterdir():
        if not entry.is_dir():
            continue
        if entry.name == RECOMPUTED_DIRNAME:
            continue
        if entry.name.startswith("campaign="):
            campaign_dirs.append(entry)
    campaign_dirs.sort(key=lambda path: (path.stat().st_mtime, path.name), reverse=True)
    return campaign_dirs


def resolve_campaign_argument(campaign_arg: str, source_root: Path) -> str:
    candidate = campaign_arg.strip()
    campaigns = list_campaign_dirs(source_root)
    if re.fullmatch(r"\d+", candidate):
        index = int(candidate)
        if index < 1 or index > len(campaigns):
            raise ValueError(f"Campaign index {index} out of range (1..{len(campaigns)}).")
        return campaigns[index - 1].name
    if not candidate.startswith("campaign="):
        candidate = f"campaign={candidate}"
    return candidate


def discover_run_groups(source_root: Path, campaign_name: str) -> list[RunGroup]:
    campaign_dir = source_root / campaign_name
    if not campaign_dir.exists() or not campaign_dir.is_dir():
        raise ValueError(f"Campaign directory not found: {campaign_dir}")

    grouped: dict[Path, list[ScanFileWorkItem]] = defaultdict(list)
    for source_csv in sorted(campaign_dir.glob(f"**/{SCAN_GLOB}")):
        source_csv = source_csv.resolve()
        if not source_csv.is_file():
            continue
        if is_under(source_csv, source_root / RECOMPUTED_DIRNAME):
            continue
        work_item = ScanFileWorkItem(source_csv=source_csv, source_root=source_root)
        grouped[work_item.run_dir.resolve()].append(work_item)

    run_groups: list[RunGroup] = []
    for run_dir, work_items in grouped.items():
        ordered = tuple(sorted(work_items, key=lambda item: str(item.relative_csv)))
        run_groups.append(RunGroup(run_dir=run_dir, source_root=source_root, scan_files=ordered))

    run_groups.sort(key=lambda group: str(group.relative_run_dir))
    return run_groups


def split_run_groups_by_subset(groups: Sequence[RunGroup], subset: tuple[int, int]) -> list[RunGroup]:
    subset_index, subset_count = subset
    return [group for index, group in enumerate(groups) if index % subset_count == subset_index]


def stop_process(process: object) -> None:
    terminate = getattr(process, "terminate", None)
    kill = getattr(process, "kill", None)
    wait = getattr(process, "wait", None)
    if callable(terminate):
        try:
            _ = terminate()
        except OSError:
            pass
    if callable(wait):
        try:
            _ = wait(timeout=1.0)
            return
        except (OSError, subprocess.TimeoutExpired):
            pass
    if callable(kill):
        try:
            _ = kill()
        except OSError:
            pass
    if callable(wait):
        try:
            _ = wait(timeout=1.0)
        except (OSError, subprocess.TimeoutExpired):
            pass


def join_stream_threads(
    threads: Sequence[threading.Thread],
    sinks: Sequence[list[str]],
    *,
    total_timeout_sec: float,
    join_slice_sec: float = STREAM_JOIN_SLICE_SEC,
    quiet_sec: float = STREAM_DRAIN_QUIET_SEC,
) -> None:
    deadline = time.perf_counter() + max(0.0, total_timeout_sec)
    last_change_perf = time.perf_counter()
    previous_lengths = [len(sink) for sink in sinks]
    while True:
        alive_threads = [thread for thread in threads if thread.is_alive()]
        current_lengths = [len(sink) for sink in sinks]
        if current_lengths != previous_lengths:
            previous_lengths = current_lengths
            last_change_perf = time.perf_counter()
        now_perf = time.perf_counter()
        if not alive_threads:
            return
        if now_perf >= deadline and now_perf - last_change_perf >= quiet_sec:
            return
        join_timeout = min(join_slice_sec, max(0.0, deadline - now_perf))
        for thread in alive_threads:
            thread.join(timeout=join_timeout)


def extract_total_attempts(stdout_text: str) -> str:
    match = TOTAL_ATTEMPTS_PATTERN.search(stdout_text)
    if match is None:
        return ""
    return match.group(1)


def extract_triple_ok_points(stdout_text: str) -> str:
    marker_match = TRIPLE_OK_POINTS_PATTERN.search(stdout_text)
    if marker_match is not None:
        return marker_match.group(1)
    stripped = stdout_text.strip()
    if not stripped:
        return ""
    try:
        payload = cast(object, json.loads(stripped))
    except json.JSONDecodeError:
        payload = None
    if isinstance(payload, dict) and "triple_ok_points" in payload:
        payload_dict = cast(Mapping[str, object], payload)
        return normalize_value(payload_dict.get("triple_ok_points"))
    match = re.search(r"triple_ok_points\s*[:=]\s*([0-9]+)", stdout_text, flags=re.IGNORECASE)
    if match is None:
        return ""
    return match.group(1)


def log_message(message: str, log_path: Path) -> None:
    timestamp = time.strftime("[%Y-%m-%d %H:%M:%S]")
    full_line = f"{timestamp} {message}"
    print(full_line, flush=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(full_line + "\n")


def pick_first_present(row: Mapping[str, str], fieldnames: Sequence[str]) -> str:
    for fieldname in fieldnames:
        value = normalize_value(row.get(fieldname))
        if value:
            return value
    raise KeyError(fieldnames[0])


def parse_finite_float(raw_value: str, fieldname: str, source_csv: Path, row_index: int) -> float:
    try:
        parsed = float(raw_value)
    except ValueError as exc:
        raise ValueError(
            f"row {row_index} in {source_csv} has invalid float for {fieldname!r}: {raw_value!r}"
        ) from exc
    if not math.isfinite(parsed):
        raise ValueError(
            f"row {row_index} in {source_csv} has non-finite value for {fieldname!r}: {raw_value!r}"
        )
    return parsed



# Filtering helpers
def is_triple_ok_row(row: Mapping[str, str]) -> bool:
    def is_one(val: object) -> bool:
        if val is None:
            return False
        try:
            return float(str(val).strip()) == 1.0
        except Exception:
            return False

    return (
        is_one(row.get("positivity_ok"))
        and is_one(row.get("unitarity_ok"))
        and is_one(row.get("perturbativity_ok"))
    )

def build_recompute_input_from_source_filtered(
    source_csv: Path, input_csv_path: Path, source_filter: str
) -> tuple[int, int]:
    """
    Returns: (total_source_rows, selected_rows)
    """
    with source_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames or []
        missing_required = [field for field in REQUIRED_NUMERIC_FIELDS if field not in header]
        if missing_required:
            raise ValueError(
                f"source CSV {source_csv} is missing required columns: {missing_required!r}"
            )
        if not any(field in header for field in LAM1_FALLBACK_FIELDS):
            raise ValueError(
                f"source CSV {source_csv} is missing any usable lam1 column from {LAM1_FALLBACK_FIELDS!r}"
            )
        if source_filter == "triple_ok":
            for col in ("positivity_ok", "unitarity_ok", "perturbativity_ok"):
                if col not in header:
                    raise ValueError(f"source CSV {source_csv} missing required column '{col}' for triple_ok filtering")
        input_csv_path.parent.mkdir(parents=True, exist_ok=True)
        total_rows = 0
        selected_rows = 0
        with input_csv_path.open("w", encoding="utf-8", newline="") as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=RECOMPUTE_INPUT_FIELDNAMES)
            writer.writeheader()
            for row_index, row in enumerate(reader):
                total_rows += 1
                if source_filter == "triple_ok" and not is_triple_ok_row(row):
                    continue
                prepared: dict[str, str] = {"point_id": str(row_index)}
                for fieldname in REQUIRED_NUMERIC_FIELDS:
                    raw_value = normalize_value(row.get(fieldname))
                    if not raw_value:
                        raise ValueError(
                            f"row {row_index} in {source_csv} is missing required field {fieldname!r}"
                        )
                    parsed = parse_finite_float(raw_value, fieldname, source_csv, row_index)
                    prepared[fieldname] = repr(parsed)
                lam1_value = pick_first_present(row, LAM1_FALLBACK_FIELDS)
                parsed_lam1 = parse_finite_float(lam1_value, "lam1", source_csv, row_index)
                prepared["lam1"] = repr(parsed_lam1)
                writer.writerow(prepared)
                selected_rows += 1
    return total_rows, selected_rows


def build_command(binary: Path, input_csv_path: Path, output_csv_path: Path, yukawas_type: int) -> list[str]:
    return [
        str(binary),
        "quarantine",
        str(input_csv_path),
        str(output_csv_path),
        str(yukawas_type),
    ]


def read_scan_meta_event(meta_path: Path) -> str | None:
    try:
        if not meta_path.exists():
            return None
        payload = json.loads(meta_path.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            event = payload.get("event")
            if isinstance(event, str) and event:
                return event
        return None
    except Exception:
        return None


def make_run_manifest_payload(
    *,
    source_root: Path,
    recomputed_root: Path,
    source_run_dir: Path,
    mirror_run_dir: Path,
    binary: Path,
    yukawas_type: int,
    omp_num_threads: int | None,
    worker_subset: str,
    dry_run: bool,
    scan_files: Sequence[ScanFileWorkItem],
    summary: Mapping[str, object],
) -> dict[str, object]:
    return {
        "recompute_mode": True,
        "created_at_epoch": time.time(),
        "source_root": str(source_root),
        "recomputed_root": str(recomputed_root),
        "source_run_dir": str(source_run_dir),
        "mirror_run_dir": str(mirror_run_dir),
        "source_campaign": scan_files[0].source_campaign if scan_files else "",
        "source_run_relative": str(source_run_dir.relative_to(source_root)),
        "fixed_dir_relative": str(source_run_dir.parent.relative_to(source_root)),
        "binary": str(binary),
        "yukawas_type": yukawas_type,
        "omp_num_threads": omp_num_threads,
        "worker_subset": worker_subset,
        "dry_run": dry_run,
        "source_scan_files": [str(item.relative_csv) for item in scan_files],
        "summary": dict(summary),
    }


def invoke_batch_recompute(
    *,
    binary: Path,
    source_csv: Path,
    destination_csv: Path,
    destination_tb_dir: Path,
    yukawas_type: int,
    omp_num_threads: int | None,
    process_timeout_sec: float,
    progress_callback: Callable[[float, float], None] | None,
    progress_every_sec: float,
) -> dict[str, object]:
    # started_at_epoch = time.time()
    # started_perf = time.perf_counter()
    # source_sha256 = sha256_file(source_csv)

    # with tempfile.TemporaryDirectory(prefix="run_quarantine_batch_") as temp_dir_str:
    #     temp_dir = Path(temp_dir_str)
    #     input_csv_path = temp_dir / "input_points.csv"
    #     source_row_count = build_recompute_input_from_source(source_csv, input_csv_path)
    #     output_csv_path = temp_dir / destination_csv.name
    #     command = build_command(binary=binary, input_csv_path=input_csv_path, output_csv_path=output_csv_path, yukawas_type=yukawas_type)
    started_at_epoch = time.time()
    started_perf = time.perf_counter()
    source_sha256 = sha256_file(source_csv)

    # IMPORTANTE:
    # el temp debe vivir en el mismo filesystem que destination_csv para que
    # os.replace(...) no falle con "Invalid cross-device link".
    destination_tb_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix=".run_quarantine_batch_",
        dir=str(destination_tb_dir),
    ) as temp_dir_str:
        temp_dir = Path(temp_dir_str)
        input_csv_path = temp_dir / "input_points.csv"
        source_filter = globals().get("CURRENT_SOURCE_FILTER", "triple_ok")
        total_source_rows, selected_rows = build_recompute_input_from_source_filtered(
            source_csv, input_csv_path, source_filter
        )
        output_csv_path = temp_dir / destination_csv.name
        command = build_command(
            binary=binary,
            input_csv_path=input_csv_path,
            output_csv_path=output_csv_path,
            yukawas_type=yukawas_type,
        )
#---------
        process_env = dict(os.environ)
        if omp_num_threads is not None:
            process_env["OMP_NUM_THREADS"] = str(omp_num_threads)

        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            env=process_env,
        )
        stdout_lines: list[str] = []
        stderr_lines: list[str] = []

        def drain_stream(stream: object, sink: list[str]) -> None:
            readline = getattr(stream, "readline", None)
            close = getattr(stream, "close", None)
            try:
                if not callable(readline):
                    return
                while True:
                    chunk = readline()
                    if not isinstance(chunk, str) or chunk == "":
                        break
                    sink.append(chunk)
            finally:
                if callable(close):
                    _ = close()

        threads = [
            threading.Thread(target=drain_stream, args=(process.stdout, stdout_lines), daemon=True),
            threading.Thread(target=drain_stream, args=(process.stderr, stderr_lines), daemon=True),
        ]
        for thread in threads:
            thread.start()

        last_progress_perf = started_perf
        while True:
            returncode = process.poll()
            if returncode is not None:
                break
            now_epoch = time.time()
            now_perf = time.perf_counter()
            if now_perf - started_perf >= process_timeout_sec:
                stop_process(process)
                raise TimeoutError(
                    f"child process exceeded timeout after {process_timeout_sec:g}s for {source_csv}"
                )
            if progress_callback is not None and (
                progress_every_sec <= 0 or now_perf - last_progress_perf >= progress_every_sec
            ):
                progress_callback(now_epoch, now_perf - started_perf)
                last_progress_perf = now_perf
            time.sleep(POLL_INTERVAL_SEC)

        join_stream_threads(threads, (stdout_lines, stderr_lines), total_timeout_sec=STREAM_JOIN_TIMEOUT_SEC)
        stdout_text = "".join(stdout_lines)
        stderr_text = "".join(stderr_lines)

        destination_tb_dir.mkdir(parents=True, exist_ok=True)
        (destination_tb_dir / "stdout.log").write_text(stdout_text, encoding="utf-8")
        (destination_tb_dir / "stderr.log").write_text(stderr_text, encoding="utf-8")

        finished_at_epoch = time.time()
        elapsed_sec = time.perf_counter() - started_perf
        total_attempts = extract_total_attempts(stdout_text)
        triple_ok_points = extract_triple_ok_points(stdout_text)
        process_returncode = getattr(process, "returncode", None)
        returncode = returncode if process_returncode is None else process_returncode

        output_row_count = 0
        destination_sha = ""
        error_code = ""
        error_detail = ""
        run_status = "OK"
        if returncode == 0:
            if not output_csv_path.exists():
                run_status = "ERROR"
                error_code = "MISSING_OUTPUT_CSV"
                error_detail = f"calculator returned 0 but did not create output file for {source_csv}"
            else:
                output_row_count = count_csv_data_rows(output_csv_path)
                if output_row_count != selected_rows:
                    run_status = "ERROR"
                    error_code = "OUTPUT_ROW_COUNT_MISMATCH"
                    error_detail = (
                        f"selected row count {selected_rows} != output row count {output_row_count} for {source_csv}"
                    )
                else:
                    destination_csv.parent.mkdir(parents=True, exist_ok=True)
                    output_csv_path.replace(destination_csv)
                    destination_sha = sha256_file(destination_csv)
        else:
            run_status = "ERROR"
            error_code = "NONZERO_EXIT"
            error_detail = stderr_text.strip() or stdout_text.strip()

        return {
            "source_row_count_total": total_source_rows,
            "source_row_count_selected": selected_rows,
            "output_row_count": output_row_count,
            "started_at_epoch": f"{started_at_epoch:.6f}",
            "finished_at_epoch": f"{finished_at_epoch:.6f}",
            "elapsed_sec": f"{max(0.0, elapsed_sec):.6f}",
            "calculator_returncode": str(returncode),
            "total_attempts": total_attempts,
            "triple_ok_points": triple_ok_points,
            "run_status": run_status,
            "error_code": error_code,
            "error_detail": error_detail,
            "source_sha256": source_sha256,
            "destination_sha256": destination_sha,
        }


def make_status_payload(
    *,
    worker_subset: str,
    subset_total_runs: int,
    campaign_total_runs: int,
    processed_runs: int,
    processed_files: int,
    ok_files: int,
    skipped_files: int,
    error_files: int,
    manifest_sha256: str,
    started_at_epoch: float,
    now_epoch: float,
    runs_per_sec: float,
    alive: bool,
    done: bool,
    current_run_relative: str,
) -> dict[str, object]:
    return {
        "alive": alive,
        "done": done,
        "worker_subset": worker_subset,
        "subset_total_runs": subset_total_runs,
        "campaign_total_runs": campaign_total_runs,
        "processed_runs": processed_runs,
        "processed_files": processed_files,
        "ok_files": ok_files,
        "skipped_files": skipped_files,
        "error_files": error_files,
        "manifest_sha256": manifest_sha256,
        "started_at_epoch": started_at_epoch,
        "last_update_epoch": now_epoch,
        "elapsed_sec": round(max(0.0, now_epoch - started_at_epoch), 6),
        "runs_per_sec": round(max(0.0, runs_per_sec), 6),
        "current_run_relative": current_run_relative,
    }


def emit_progress(status: Mapping[str, object]) -> None:
    print(
        " ".join(
            [
                f"alive={normalize_value(status.get('alive'))}",
                f"worker_subset={normalize_value(status.get('worker_subset'))}",
                f"runs={normalize_value(status.get('processed_runs'))}/{normalize_value(status.get('subset_total_runs'))}",
                f"campaign_runs={normalize_value(status.get('campaign_total_runs'))}",
                f"files={normalize_value(status.get('processed_files'))}",
                f"ok={normalize_value(status.get('ok_files'))}",
                f"skipped={normalize_value(status.get('skipped_files'))}",
                f"errors={normalize_value(status.get('error_files'))}",
                f"runs_per_sec={normalize_value(status.get('runs_per_sec'))}",
                f"current_run={normalize_value(status.get('current_run_relative'))}",
                f"manifest_sha256={normalize_value(status.get('manifest_sha256'))}",
            ]
        ),
        flush=True,
    )


def write_status_json(path: Path | None, payload: Mapping[str, object]) -> None:
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    tmp_path.write_text(json.dumps(dict(payload), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    tmp_path.replace(path)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Recompute live dihiggs_lake campaign CSVs into a mirrored recomputed/ tree, "
            "preserving the original campaign/fixed/run/tb/scan layout."
        )
    )
    source_group = parser.add_mutually_exclusive_group(required=True)
    _ = source_group.add_argument(
        "--campaign",
        help=(
            "Campaign name (campaign=...) or numeric index from --list-campaigns output. "
            "Index 1 is newest; largest index is oldest."
        ),
    )
    _ = source_group.add_argument("--list-campaigns", action="store_true")
    _ = parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    _ = parser.add_argument(
        "--recomputed-root",
        type=Path,
        default=None,
        help="Destination mirror root. Default: <source-root>/recomputed",
    )
    _ = parser.add_argument("--subset", type=parse_subset, default=(0, 1))
    _ = parser.add_argument(
        "--max-runs",
        "--max-files",
        dest="max_runs",
        type=parse_non_negative_int,
        default=None,
        help="Limit the number of source run directories selected before subset splitting.",
    )
    _ = parser.add_argument("--status-csv", type=Path, default=None)
    _ = parser.add_argument("--status-json", type=Path, default=None)
    _ = parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    _ = parser.add_argument("--yukawas-type", type=int, default=1)
    _ = parser.add_argument("--omp-num-threads", type=parse_positive_int, default=None)
    _ = parser.add_argument("--process-timeout-sec", type=parse_positive_float, default=DEFAULT_PROCESS_TIMEOUT_SEC)
    _ = parser.add_argument("--progress-every-sec", type=parse_positive_float, default=5.0)
    _ = parser.add_argument("--force", action="store_true")
    _ = parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:

    args = parse_args(argv)
    global CURRENT_SOURCE_FILTER
    CURRENT_SOURCE_FILTER = getattr(args, "source_filter", "triple_ok")

    source_root = cast(Path, args.source_root).expanduser().resolve()
    campaign_arg = cast(str | None, args.campaign)
    list_campaigns = cast(bool, args.list_campaigns)
    recomputed_root_arg = cast(Path | None, args.recomputed_root)
    recomputed_root = (
        recomputed_root_arg.expanduser().resolve()
        if recomputed_root_arg is not None
        else (source_root / RECOMPUTED_DIRNAME).resolve()
    )

    if recomputed_root == source_root:
        print("Recomputed root must not be the same as source root.", file=sys.stderr)
        return 1

    if list_campaigns:
        campaigns = list_campaign_dirs(source_root)
        for index, campaign_dir in enumerate(campaigns, start=1):
            print(f"{index} {campaign_dir.name}")
        return 0

    if campaign_arg is None:
        print("--campaign required unless --list-campaigns is used.", file=sys.stderr)
        return 2

    try:
        campaign_name = resolve_campaign_argument(campaign_arg, source_root)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    try:
        campaign_groups = discover_run_groups(source_root, campaign_name)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    if not campaign_groups:
        print(f"Campaign '{campaign_name}' contains no {SCAN_GLOB} files.", file=sys.stderr)
        return 1

    max_runs = cast(int | None, args.max_runs)
    if max_runs is not None:
        campaign_groups = campaign_groups[:max_runs]

    subset = cast(tuple[int, int], args.subset)
    selected_groups = split_run_groups_by_subset(campaign_groups, subset)
    manifest_sha256 = hashlib.sha256(
        json.dumps([str(group.relative_run_dir) for group in selected_groups], ensure_ascii=False).encode("utf-8")
    ).hexdigest()
    worker_subset = f"{subset[0]}/{subset[1]}"

    status_json = cast(Path | None, args.status_json)
    status_csv = cast(Path | None, args.status_csv)
    binary = cast(Path, args.binary).expanduser().resolve()
    if not binary.exists() and not cast(bool, args.dry_run):
        print(f"Binary not found: {binary}", file=sys.stderr)
        return 1

    progress_every_sec = cast(float, args.progress_every_sec)
    process_timeout_sec = cast(float, args.process_timeout_sec)
    yukawas_type = cast(int, args.yukawas_type)
    omp_num_threads = cast(int | None, args.omp_num_threads)
    dry_run = cast(bool, args.dry_run)
    force = cast(bool, args.force)

    started_at_epoch = time.time()
    recent_completion_epochs: deque[float] = deque(maxlen=RUNS_PER_SEC_WINDOW)
    current_run_started_at_epoch: float | None = None
    current_run_relative = ""
    processed_runs = 0
    processed_files = 0
    ok_files = 0
    skipped_files = 0
    error_files = 0
    subset_total_runs = len(selected_groups)
    campaign_total_runs = len(campaign_groups)

    def current_runs_per_sec(now_epoch: float) -> float:
        if current_run_started_at_epoch is not None:
            current_run_elapsed_sec = now_epoch - current_run_started_at_epoch
            if current_run_elapsed_sec > 0:
                return 1.0 / current_run_elapsed_sec
        while recent_completion_epochs and now_epoch < recent_completion_epochs[0]:
            _ = recent_completion_epochs.popleft()
        if len(recent_completion_epochs) < 2:
            return 0.0
        elapsed_sec = now_epoch - recent_completion_epochs[0]
        if elapsed_sec <= 0:
            return 0.0
        return (len(recent_completion_epochs) - 1) / elapsed_sec

    initial_status = make_status_payload(
        worker_subset=worker_subset,
        subset_total_runs=subset_total_runs,
        campaign_total_runs=campaign_total_runs,
        processed_runs=processed_runs,
        processed_files=processed_files,
        ok_files=ok_files,
        skipped_files=skipped_files,
        error_files=error_files,
        manifest_sha256=manifest_sha256,
        started_at_epoch=started_at_epoch,
        now_epoch=started_at_epoch,
        runs_per_sec=0.0,
        alive=True,
        done=False,
        current_run_relative=current_run_relative,
    )
    write_status_json(status_json, initial_status)
    emit_progress(initial_status)

    if not selected_groups:
        finished_status = make_status_payload(
            worker_subset=worker_subset,
            subset_total_runs=0,
            campaign_total_runs=campaign_total_runs,
            processed_runs=0,
            processed_files=0,
            ok_files=0,
            skipped_files=0,
            error_files=0,
            manifest_sha256=manifest_sha256,
            started_at_epoch=started_at_epoch,
            now_epoch=time.time(),
            runs_per_sec=0.0,
            alive=False,
            done=True,
            current_run_relative="",
        )
        write_status_json(status_json, finished_status)
        emit_progress(finished_status)
        return 0

    def emit_live_status(now_epoch: float, _perf_elapsed_sec: float) -> None:
        status = make_status_payload(
            worker_subset=worker_subset,
            subset_total_runs=subset_total_runs,
            campaign_total_runs=campaign_total_runs,
            processed_runs=processed_runs,
            processed_files=processed_files,
            ok_files=ok_files,
            skipped_files=skipped_files,
            error_files=error_files,
            manifest_sha256=manifest_sha256,
            started_at_epoch=started_at_epoch,
            now_epoch=now_epoch,
            runs_per_sec=current_runs_per_sec(now_epoch),
            alive=True,
            done=False,
            current_run_relative=current_run_relative,
        )
        write_status_json(status_json, status)
        emit_progress(status)

    for group in selected_groups:
        current_run_started_at_epoch = time.time()
        current_run_relative = str(group.relative_run_dir)
        mirror_run_dir = recomputed_root / group.relative_run_dir
        mirror_run_dir.mkdir(parents=True, exist_ok=True)
        run_log_path = mirror_run_dir / "orchestrator.log"
        run_summary_path = mirror_run_dir / "task_summary.jsonl"
        run_manifest_path = mirror_run_dir / "run_manifest.json"

        run_summary = {
            "tasks_total": len(group.scan_files),
            "tasks_ok": 0,
            "tasks_skipped": 0,
            "tasks_failed": 0,
            "started_at_epoch": current_run_started_at_epoch,
        }
        write_json_atomic(
            run_manifest_path,
            make_run_manifest_payload(
                source_root=source_root,
                recomputed_root=recomputed_root,
                source_run_dir=group.run_dir,
                mirror_run_dir=mirror_run_dir,
                binary=binary,
                yukawas_type=yukawas_type,
                omp_num_threads=omp_num_threads,
                worker_subset=worker_subset,
                dry_run=dry_run,
                scan_files=group.scan_files,
                summary=run_summary,
            ),
        )
        log_message(f"[INIT] source_run={group.relative_run_dir}", run_log_path)
        log_message(f"[CONF] dry_run={dry_run} force={force} omp_num_threads={omp_num_threads}", run_log_path)

        for scan_file in group.scan_files:
            processed_files += 1
            destination_csv = recomputed_root / scan_file.relative_csv
            destination_tb_dir = destination_csv.parent
            scan_meta_path = destination_tb_dir / "scan_meta.json"
            source_csv_relative = str(scan_file.relative_csv)
            source_run_relative = str(group.relative_run_dir)
            file_started_at_epoch = time.time()

            source_row_count = ""
            output_row_count = ""
            run_status = ""
            error_code = ""
            error_detail = ""
            finished_at_epoch = ""
            elapsed_sec = ""
            calculator_returncode = ""
            total_attempts = ""
            triple_ok_points = ""
            source_sha256 = ""
            destination_sha256 = ""

            try:
                if not scan_file.source_csv.exists():
                    raise FileNotFoundError(f"source CSV not found: {scan_file.source_csv}")

                source_sha256 = sha256_file(scan_file.source_csv)

                # Inicializar SIEMPRE antes de cualquier posible excepción
                source_row_count_total = 0
                source_row_count_selected = 0
                source_row_count = "0"

                # Obtener conteos con un temp Path válido
                # with tempfile.TemporaryDirectory(prefix=".rowcount_probe_", dir=str(destination_tb_dir)) as probe_dir_str:
                #     probe_dir = Path(probe_dir_str)
                #     probe_input_csv = probe_dir / "probe_input.csv"
                #     total_source_rows, selected_rows = build_recompute_input_from_source_filtered(
                #         scan_file.source_csv,
                #         probe_input_csv,
                #         CURRENT_SOURCE_FILTER,
                #     )
                with tempfile.TemporaryDirectory(prefix=".rowcount_probe_") as probe_dir_str:
                    probe_dir = Path(probe_dir_str)
                    probe_input_csv = probe_dir / "probe_input.csv"
                    total_source_rows, selected_rows = build_recompute_input_from_source_filtered(
                        scan_file.source_csv,
                        probe_input_csv,
                        CURRENT_SOURCE_FILTER,
                    )

                source_row_count_total = int(total_source_rows)
                source_row_count_selected = int(selected_rows)
                source_row_count = str(source_row_count_total)

                existing_event = read_scan_meta_event(scan_meta_path)
                meta_payload = None
                if scan_meta_path.exists():
                    try:
                        meta_payload = json.loads(scan_meta_path.read_text(encoding="utf-8"))
                    except Exception:
                        meta_payload = None
                meta_source_filter = meta_payload.get("source_filter") if meta_payload else None
                meta_output_row_count = meta_payload.get("output_row_count") if meta_payload else None
                if (
                    destination_csv.exists()
                    and existing_event == "done"
                    and not force
                    and meta_source_filter == CURRENT_SOURCE_FILTER
                    and str(meta_output_row_count) == str(source_row_count_selected)
                ):
                    destination_row_count = count_csv_data_rows(destination_csv)
                    if destination_row_count == source_row_count_selected:
                        run_status = "SKIPPED"
                        output_row_count = str(destination_row_count)
                        destination_sha256 = sha256_file(destination_csv)
                        error_detail = "destination exists, scan_meta.json indicates done, and selected row counts match"
                        run_summary["tasks_skipped"] = int(run_summary["tasks_skipped"]) + 1
                        skipped_files += 1
                        append_jsonl(
                            run_summary_path,
                            {
                                "event": "skip",
                                "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                "source_run_relative": source_run_relative,
                                "source_csv_relative": source_csv_relative,
                                "destination_csv": str(destination_csv),
                                "worker_subset": worker_subset,
                                "source_row_count_total": source_row_count_total,
                                "source_row_count_selected": source_row_count_selected,
                                "output_row_count": destination_row_count,
                                "source_filter": CURRENT_SOURCE_FILTER,
                                "reason": "existing_done_output",
                            },
                        )
                        log_message(f"[SKIP] {source_csv_relative}", run_log_path)
                        continue
                    else:
                        log_message(
                            f"[WARN] stale recomputed file row-count mismatch for {source_csv_relative}: "
                            f"selected={source_row_count_selected} dest={destination_row_count}; recomputing",
                            run_log_path,
                        )
                        existing_event = None

                if not run_status:
                    if source_row_count_selected == 0:
                        # No triple_ok points: skip, do not call binary, write meta and summary
                        run_status = "SKIPPED"
                        error_detail = "No triple_ok points selected; skipping recomputation."
                        run_summary["tasks_skipped"] = int(run_summary["tasks_skipped"]) + 1
                        skipped_files += 1
                        destination_tb_dir.mkdir(parents=True, exist_ok=True)
                        write_json_atomic(
                            scan_meta_path,
                            {
                                "event": "skip_empty_selection",
                                "source_csv": str(scan_file.source_csv),
                                "destination_csv": str(destination_csv),
                                "source_csv_relative": source_csv_relative,
                                "source_run_relative": source_run_relative,
                                "source_row_count_total": source_row_count_total,
                                "source_row_count_selected": source_row_count_selected,
                                "output_row_count": 0,
                                "source_filter": CURRENT_SOURCE_FILTER,
                                "binary": str(binary),
                                "yukawas_type": yukawas_type,
                                "omp_num_threads": omp_num_threads,
                            },
                        )
                        append_jsonl(
                            run_summary_path,
                            {
                                "event": "skip_empty_selection",
                                "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                "source_run_relative": source_run_relative,
                                "source_csv_relative": source_csv_relative,
                                "destination_csv": str(destination_csv),
                                "worker_subset": worker_subset,
                                "source_row_count_total": source_row_count_total,
                                "source_row_count_selected": source_row_count_selected,
                                "output_row_count": 0,
                                "source_filter": CURRENT_SOURCE_FILTER,
                                "reason": "no triple_ok points",
                            },
                        )
                        log_message(f"[SKIP] {source_csv_relative} (no triple_ok points)", run_log_path)
                        continue
                    if dry_run:
                        run_status = "SKIPPED"
                        error_detail = "dry-run: command was not executed"
                        run_summary["tasks_skipped"] = int(run_summary["tasks_skipped"]) + 1
                        skipped_files += 1
                        destination_tb_dir.mkdir(parents=True, exist_ok=True)
                        write_json_atomic(
                            scan_meta_path,
                            {
                                "event": "dry_run",
                                "source_csv": str(scan_file.source_csv),
                                "destination_csv": str(destination_csv),
                                "source_csv_relative": source_csv_relative,
                                "source_run_relative": source_run_relative,
                                "source_row_count_total": source_row_count_total,
                                "source_row_count_selected": source_row_count_selected,
                                "output_row_count": 0,
                                "source_filter": CURRENT_SOURCE_FILTER,
                                "binary": str(binary),
                                "yukawas_type": yukawas_type,
                                "omp_num_threads": omp_num_threads,
                            },
                        )
                        append_jsonl(
                            run_summary_path,
                            {
                                "event": "dry_run",
                                "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                "source_run_relative": source_run_relative,
                                "source_csv_relative": source_csv_relative,
                                "destination_csv": str(destination_csv),
                                "worker_subset": worker_subset,
                                "source_row_count_total": source_row_count_total,
                                "source_row_count_selected": source_row_count_selected,
                                "output_row_count": 0,
                                "source_filter": CURRENT_SOURCE_FILTER,
                            },
                        )
                        log_message(f"[DRY ] {source_csv_relative}", run_log_path)
                        continue
                    # Normal recompute
                    result = invoke_batch_recompute(
                        binary=binary,
                        source_csv=scan_file.source_csv,
                        destination_csv=destination_csv,
                        destination_tb_dir=destination_tb_dir,
                        yukawas_type=yukawas_type,
                        omp_num_threads=omp_num_threads,
                        process_timeout_sec=process_timeout_sec,
                        progress_callback=emit_live_status,
                        progress_every_sec=progress_every_sec,
                    )
                    run_status = normalize_value(result.get("run_status"))
                    error_code = normalize_value(result.get("error_code"))
                    error_detail = normalize_value(result.get("error_detail"))
                    finished_at_epoch = normalize_value(result.get("finished_at_epoch"))
                    elapsed_sec = normalize_value(result.get("elapsed_sec"))
                    calculator_returncode = normalize_value(result.get("calculator_returncode"))
                    total_attempts = normalize_value(result.get("total_attempts"))
                    triple_ok_points = normalize_value(result.get("triple_ok_points"))
                    source_sha256 = normalize_value(result.get("source_sha256")) or source_sha256
                    destination_sha256 = normalize_value(result.get("destination_sha256"))
                    output_row_count = normalize_value(result.get("output_row_count"))
                    source_row_count_total = normalize_value(result.get("source_row_count_total"))
                    source_row_count_selected = normalize_value(result.get("source_row_count_selected"))

                    event_name = "done" if run_status == "OK" else "fail"
                    event_payload = {
                        "event": event_name,
                        "source_csv": str(scan_file.source_csv),
                        "destination_csv": str(destination_csv),
                        "source_csv_relative": source_csv_relative,
                        "source_run_relative": source_run_relative,
                        "source_row_count_total": int(source_row_count_total or "0"),
                        "source_row_count_selected": int(source_row_count_selected or "0"),
                        "output_row_count": int(output_row_count or "0"),
                        "source_filter": CURRENT_SOURCE_FILTER,
                        "calculator_returncode": calculator_returncode,
                        "total_attempts": total_attempts,
                        "triple_ok_points": triple_ok_points,
                        "source_sha256": source_sha256,
                        "destination_sha256": destination_sha256,
                        "error_code": error_code,
                        "error_detail": error_detail,
                        "binary": str(binary),
                        "yukawas_type": yukawas_type,
                        "omp_num_threads": omp_num_threads,
                    }
                    write_json_atomic(scan_meta_path, event_payload)
                    append_jsonl(
                        run_summary_path,
                        {
                            "event": event_name,
                            "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                            "source_run_relative": source_run_relative,
                            "source_csv_relative": source_csv_relative,
                            "destination_csv": str(destination_csv),
                            "worker_subset": worker_subset,
                            "source_row_count_total": int(source_row_count_total or "0"),
                            "source_row_count_selected": int(source_row_count_selected or "0"),
                            "output_row_count": int(output_row_count or "0"),
                            "source_filter": CURRENT_SOURCE_FILTER,
                            "calculator_returncode": calculator_returncode,
                            "total_attempts": total_attempts,
                            "triple_ok_points": triple_ok_points,
                            "error_code": error_code,
                            "error_detail": error_detail,
                        },
                    )
                    if run_status == "OK":
                        ok_files += 1
                        run_summary["tasks_ok"] = int(run_summary["tasks_ok"]) + 1
                        log_message(
                            f"[DONE] {source_csv_relative} total={source_row_count_total} selected={source_row_count_selected} triple_ok={triple_ok_points}",
                            run_log_path,
                        )
                    else:
                        error_files += 1
                        run_summary["tasks_failed"] = int(run_summary["tasks_failed"]) + 1
                        log_message(
                            f"[FAIL] {source_csv_relative} error={error_code or 'UNKNOWN'} detail={error_detail}",
                            run_log_path,
                        )
            except Exception as exc:
                run_status = "ERROR"
                error_code = exc.__class__.__name__
                error_detail = str(exc)
                error_files += 1
                run_summary["tasks_failed"] = int(run_summary["tasks_failed"]) + 1
                destination_tb_dir.mkdir(parents=True, exist_ok=True)
                write_json_atomic(
                    scan_meta_path,
                    {
                        "event": "crash",
                        "source_csv": str(scan_file.source_csv),
                        "destination_csv": str(destination_csv),
                        "source_csv_relative": source_csv_relative,
                        "source_run_relative": source_run_relative,
                        "source_row_count": int(source_row_count or "0"),
                        "error_code": error_code,
                        "error_detail": error_detail,
                        "binary": str(binary),
                        "yukawas_type": yukawas_type,
                        "omp_num_threads": omp_num_threads,
                    },
                )
                append_jsonl(
                    run_summary_path,
                    {
                        "event": "crash",
                        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                        "source_run_relative": source_run_relative,
                        "source_csv_relative": source_csv_relative,
                        "destination_csv": str(destination_csv),
                        "worker_subset": worker_subset,
                        "source_row_count": int(source_row_count or "0"),
                        "error_code": error_code,
                        "error_detail": error_detail,
                    },
                )
                log_message(
                    f"[CRSH] {source_csv_relative} error={error_code} detail={error_detail}",
                    run_log_path,
                )
            finally:
                file_finished_at_epoch = time.time()
                if not finished_at_epoch:
                    finished_at_epoch = f"{file_finished_at_epoch:.6f}"
                if not elapsed_sec:
                    elapsed_sec = f"{max(0.0, file_finished_at_epoch - file_started_at_epoch):.6f}"
                status_row = {
                    "source_campaign": group.source_campaign,
                    "source_run_relative": source_run_relative,
                    "source_csv_relative": source_csv_relative,
                    "destination_csv": str(destination_csv),
                    "worker_subset": worker_subset,
                    "source_row_count": source_row_count_total,
                    "output_row_count": output_row_count,
                    "run_status": run_status,
                    "error_code": error_code,
                    "error_detail": error_detail,
                    "started_at_epoch": f"{file_started_at_epoch:.6f}",
                    "finished_at_epoch": finished_at_epoch,
                    "elapsed_sec": elapsed_sec,
                    "calculator_returncode": calculator_returncode,
                    "total_attempts": total_attempts,
                    "triple_ok_points": triple_ok_points,
                    "source_sha256": source_sha256,
                    "destination_sha256": destination_sha256,
                    "source_row_count_selected": source_row_count_selected,
                    "source_filter": CURRENT_SOURCE_FILTER,
                }
                if status_csv is not None:
                    append_status_row(status_csv, status_row)

        processed_runs += 1
        recent_completion_epochs.append(time.time())
        run_summary["finished_at_epoch"] = time.time()
        write_json_atomic(
            run_manifest_path,
            make_run_manifest_payload(
                source_root=source_root,
                recomputed_root=recomputed_root,
                source_run_dir=group.run_dir,
                mirror_run_dir=mirror_run_dir,
                binary=binary,
                yukawas_type=yukawas_type,
                omp_num_threads=omp_num_threads,
                worker_subset=worker_subset,
                dry_run=dry_run,
                scan_files=group.scan_files,
                summary=run_summary,
            ),
        )
        current_run_started_at_epoch = None
        emit_live_status(time.time(), 0.0)

    final_status = make_status_payload(
        worker_subset=worker_subset,
        subset_total_runs=subset_total_runs,
        campaign_total_runs=campaign_total_runs,
        processed_runs=processed_runs,
        processed_files=processed_files,
        ok_files=ok_files,
        skipped_files=skipped_files,
        error_files=error_files,
        manifest_sha256=manifest_sha256,
        started_at_epoch=started_at_epoch,
        now_epoch=time.time(),
        runs_per_sec=current_runs_per_sec(time.time()),
        alive=False,
        done=True,
        current_run_relative=current_run_relative,
    )
    write_status_json(status_json, final_status)
    emit_progress(final_status)

    if error_files > 0:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
