from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from collections.abc import Iterable
from typing import cast


JsonDict = dict[str, object]
PathLike = str | Path


_RUN_DIR_RE = re.compile(r"\[PATH\]\s+run_dir\s*=\s*(?P<path>.+?)\s*$")


def parse_task_summary_jsonl(path: PathLike) -> list[JsonDict]:
    p = Path(path)
    try:
        if not p.exists():
            return []
        records: list[JsonDict] = []
        with p.open("r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                try:
                    obj = cast(object, json.loads(s))
                except Exception:
                    continue
                if not isinstance(obj, dict):
                    continue

                obj_dict = cast(dict[object, object], obj)
                rec: JsonDict = {}
                for k, v in obj_dict.items():
                    if isinstance(k, str):
                        rec[k] = v
                records.append(rec)
        return records
    except Exception:
        return []


def parse_triple_ok_points(value: object) -> int:
    if value is None:
        return 0

    if isinstance(value, bool):
        return 0

    if isinstance(value, int):
        return max(0, value)

    if isinstance(value, float):
        try:
            return max(0, int(value))
        except Exception:
            return 0

    if isinstance(value, str):
        s = value.strip()
        if not s:
            return 0

        m = re.search(r"-?\d+", s)
        if m:
            try:
                return max(0, int(m.group(0)))
            except Exception:
                return 0
        try:
            return max(0, int(float(s)))
        except Exception:
            return 0

    return 0


def _parse_attempts(value: object) -> int:
    if value is None:
        return 0
    if isinstance(value, bool):
        return 0
    if isinstance(value, int):
        return max(0, value)
    if isinstance(value, float):
        try:
            return max(0, int(value))
        except Exception:
            return 0
    if isinstance(value, str):
        s = value.strip()
        if not s:
            return 0
        m = re.search(r"-?\d+", s)
        if not m:
            return 0
        try:
            return max(0, int(m.group(0)))
        except Exception:
            return 0
    return 0


def _is_one_flag(value: object) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return False
    if isinstance(value, int):
        return value == 1
    if isinstance(value, float):
        return value == 1.0
    if isinstance(value, str):
        s = value.strip()
        if not s:
            return False
        try:
            return float(s) == 1.0
        except Exception:
            return False
    return False


def strict_valid_count_from_csv(path: str | Path) -> int:
    try:
        p = Path(path)
    except Exception:
        return 0

    try:
        if not p.exists() or not p.is_file():
            return 0
        count = 0
        with p.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                pos = row.get("positivity_ok")
                uni = row.get("unitarity_ok")
                per = row.get("perturbativity_ok")
                if _is_one_flag(pos) and _is_one_flag(uni) and _is_one_flag(per):
                    count += 1
        return count
    except Exception:
        return 0


def _parse_axis_count_from_grid_signature(grid_signature: object, axis: str) -> int:
    if not isinstance(grid_signature, str):
        return 0
    axis_esc = re.escape(axis)
    m = re.search(rf"(?:^|\|){axis_esc}=[^|]*-N(?P<n>\d+)(?:\||$)", grid_signature)
    if not m:
        return 0
    try:
        return max(0, int(m.group("n")))
    except Exception:
        return 0


def _parse_record_path(record: JsonDict, primary_key: str, fallback_key: str) -> Path | None:
    raw_primary = record.get(primary_key)
    raw_fallback = record.get(fallback_key)
    raw = raw_primary if raw_primary else raw_fallback
    if raw is None:
        return None
    if isinstance(raw, Path):
        return raw
    if isinstance(raw, str):
        s = raw.strip()
        if not s:
            return None
        return Path(s)
    return None


def _triple_ok_points_from_scan_meta(csv_path: Path) -> int | None:
    try:
        meta_path = Path(csv_path).parent / "scan_meta.json"
        if not meta_path.exists():
            return None
        with meta_path.open("r", encoding="utf-8") as f:
            obj = cast(object, json.load(f))
        if not isinstance(obj, dict):
            return None
        d = cast(dict[object, object], obj)
        if "triple_ok_points" not in d:
            return None
        return parse_triple_ok_points(d.get("triple_ok_points"))
    except Exception:
        return None


def successes_trials_from_task_record(record: dict[str, object]) -> tuple[int, int]:
    rec = record
    event = rec.get("event")

    n_lam1_effective = _parse_attempts(rec.get("n_lam1_effective"))
    n_mphi = _parse_axis_count_from_grid_signature(rec.get("grid_signature"), "mphi")
    inferred_trials = n_mphi * n_lam1_effective

    if event == "skip":
        csv_path = _parse_record_path(rec, "previous_csv", "output_csv")
        if csv_path is None:
            return 0, inferred_trials
        meta_successes = _triple_ok_points_from_scan_meta(csv_path)
        if meta_successes is not None:
            return int(meta_successes), inferred_trials
        successes = strict_valid_count_from_csv(csv_path)
        return successes, inferred_trials

    if event == "done":
        successes = parse_triple_ok_points(rec.get("triple_ok_points"))
        trials = _parse_attempts(rec.get("attempts"))
        if trials <= 0:
            trials = inferred_trials
        return successes, trials

    return 0, 0


def summarize_task_summary(records: Iterable[JsonDict]) -> tuple[int, int]:
    attempts_total = 0
    triple_ok_total = 0

    for rec in records:
        try:
            if rec.get("event") != "done":
                continue
            attempts_total += _parse_attempts(rec.get("attempts"))
            triple_ok_total += parse_triple_ok_points(rec.get("triple_ok_points"))
        except Exception:
            continue

    return attempts_total, triple_ok_total


def parse_run_dir_from_orchestrator_output(text: str) -> Path | None:
    if not text:
        return None
    for line in text.splitlines():
        m = _RUN_DIR_RE.search(line)
        if not m:
            continue
        raw = (m.group("path") or "").strip()
        if not raw:
            return None
        try:
            return Path(raw).expanduser()
        except Exception:
            return None
    return None


def read_omp_num_threads_from_manifest(run_dir: Path) -> int | None:
    try:
        mani_path = Path(run_dir) / "run_manifest.json"
        if not mani_path.exists():
            return None
        with mani_path.open("r", encoding="utf-8") as f:
            mani_obj = cast(object, json.load(f))
        if not isinstance(mani_obj, dict):
            return None

        mani_dict = cast(dict[object, object], mani_obj)
        mani: JsonDict = {}
        for k, v in mani_dict.items():
            if isinstance(k, str):
                mani[k] = v

        runtime_any = mani.get("runtime")
        if not isinstance(runtime_any, dict):
            return None

        runtime_dict = cast(dict[object, object], runtime_any)
        runtime: JsonDict = {}
        for k, v in runtime_dict.items():
            if isinstance(k, str):
                runtime[k] = v

        v = runtime.get("omp_num_threads")
        if v is None:
            return None
        if isinstance(v, bool):
            return None
        if isinstance(v, int):
            return v
        if isinstance(v, float):
            try:
                return int(v)
            except Exception:
                return None
        if isinstance(v, str):
            s = v.strip()
            if not s:
                return None
            m = re.search(r"-?\d+", s)
            if not m:
                return None
            try:
                return int(m.group(0))
            except Exception:
                return None
        return None
    except Exception:
        return None
