from __future__ import annotations

import json
import re
from pathlib import Path
from collections.abc import Iterable
from typing import cast


JsonDict = dict[str, object]
PathLike = str | Path


_RUN_DIR_RE = re.compile(r"^\s*\[PATH\]\s+run_dir\s*=\s*(?P<path>.+?)\s*$")


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
        m = _RUN_DIR_RE.match(line)
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
