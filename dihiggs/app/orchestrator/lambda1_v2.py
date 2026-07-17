"""Minimal orchestration for the canonical ``dihiggs.lambda1.v2`` producer."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from dihiggs.app.orchestrator.io_utils import detect_git_info

INPUT_HEADER = (
    "point_id", "mh_gev", "mH_gev", "mA_gev", "mHp_gev",
    "sin_beta_minus_alpha", "tan_beta", "lambda1_target",
    "lambda6_input", "lambda7_input",
)
SCHEMA_VERSION = "dihiggs.lambda1.v2"


def format_float64(value: float) -> str:
    """Return a decimal lexeme that round-trips through IEEE Float64."""
    return format(value, ".17g")


def grid_values(low: float, high: float, count: int) -> list[float]:
    if count < 1 or (count == 1 and low != high):
        raise ValueError("grid requires count >= 1; one point requires equal bounds")
    if count == 1:
        return [low]
    values = [low + (high - low) * i / (count - 1) for i in range(count)]
    values[-1] = high
    return values


def _fnv1a64(text: str) -> int:
    value = 14695981039346656037
    for byte in text.encode():
        value = (value ^ byte) * 1099511628211 & 0xFFFFFFFFFFFFFFFF
    return value


def stable_point_id(values: Sequence[float]) -> str:
    canonical = "".join(f"{value:.17e}," for value in values)
    return f"lambda1_{_fnv1a64(canonical):016x}"


@dataclass(frozen=True)
class Lambda1Fixed:
    mh: float
    mA: float
    mHp: float
    sin_ba: float
    lambda6: float
    lambda7: float


def cartesian_rows(
    *,
    fixed: Lambda1Fixed,
    mH_values: Iterable[float],
    lambda1_values: Iterable[float],
    tan_beta_values: Iterable[float],
) -> list[dict[str, str]]:
    mH_values = list(mH_values)
    lambda1_values = list(lambda1_values)
    tan_beta_values = list(tan_beta_values)
    rows: list[dict[str, str]] = []
    for tan_beta in tan_beta_values:
        for mH in mH_values:
            for lambda1 in lambda1_values:
                values = (
                    fixed.mh, mH, fixed.mA, fixed.mHp, fixed.sin_ba,
                    tan_beta, lambda1, fixed.lambda6, fixed.lambda7,
                )
                rows.append({
                    "point_id": stable_point_id(values),
                    **dict(zip(INPUT_HEADER[1:], map(format_float64, values))),
                })
    return rows


def write_input_csv(path: Path, rows: Sequence[dict[str, str]]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=INPUT_HEADER, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return sha256(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_command(executable: Path, input_csv: Path, output_csv: Path) -> list[str]:
    return [str(executable), str(input_csv), str(output_csv)]


def validate_output(path: Path, expected_rows: int) -> dict[str, object]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames or []
        rows = list(reader)
    if not columns or columns[0] != "schema_version":
        raise ValueError("lambda1 v2 output is missing schema_version")
    if len(rows) != expected_rows:
        raise ValueError(f"lambda1 v2 row-count mismatch: expected {expected_rows}, got {len(rows)}")
    bad_schema = next((row.get("schema_version") for row in rows if row.get("schema_version") != SCHEMA_VERSION), None)
    if bad_schema is not None:
        raise ValueError(f"lambda1 v2 schema mismatch: expected {SCHEMA_VERSION}, got {bad_schema}")
    return {"schema": SCHEMA_VERSION, "row_count": len(rows), "columns": columns}


def _twohdmc_provenance(repo_root: Path) -> str:
    try:
        result = subprocess.run(
            ["git", "log", "-1", "--format=%H", "--", "2hdmc"],
            cwd=repo_root, capture_output=True, text=True, check=True,
        )
        return result.stdout.strip() or "unknown"
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def run_lambda1_v2(
    *,
    executable: Path,
    rows: Sequence[dict[str, str]],
    outdir: Path,
    campaign: str,
    run_name: str,
    repo_root: Path,
    dry_run: bool = False,
    force: bool = False,
    timeout_s: float | None = None,
) -> dict[str, object]:
    run_dir = outdir / f"campaign={campaign}" / run_name
    input_csv, output_csv, manifest_path = (
        run_dir / "input_v2.csv", run_dir / "output_v2.csv", run_dir / "run_manifest.json"
    )
    if run_dir.exists() and not force and manifest_path.exists():
        raise FileExistsError(f"run already exists: {run_dir}; use --force to overwrite")
    input_sha = write_input_csv(input_csv, rows)
    git = detect_git_info(repo_root)
    command = build_command(executable, input_csv, output_csv)
    manifest: dict[str, object] = {
        "schema_version": "orchestrator.lambda1_v2",
        "engine": "lambda1_v2",
        "producer": "Lambda1EvaluatorV2",
        "producer_schema": SCHEMA_VERSION,
        "run_name": run_name,
        "campaign": campaign,
        "command": command,
        "input_csv": str(input_csv),
        "output_csv": str(output_csv),
        "input_sha256": input_sha,
        "requested_row_count": len(rows),
        "git": git,
        "repository_commit": git.get("commit"),
        "dirty_state": git.get("is_dirty"),
        "twohdmc_provenance": _twohdmc_provenance(repo_root),
        "status": "dry_run" if dry_run else "planned",
    }
    run_dir.mkdir(parents=True, exist_ok=True)
    if dry_run:
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        return manifest

    environment = {
        **dict(os.environ),
        "DIHIGGS_GIT_COMMIT": str(git.get("commit") or "unknown"),
        "DIHIGGS_GIT_DIRTY": str(git.get("is_dirty") or "unknown"),
    }
    completed = subprocess.run(
        command, cwd=repo_root, env=environment, capture_output=True, text=True,
        timeout=timeout_s, check=False,
    )
    (run_dir / "stdout.log").write_text(completed.stdout, encoding="utf-8")
    (run_dir / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    manifest["returncode"] = completed.returncode
    if completed.returncode != 0:
        manifest["status"] = "failed"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        raise RuntimeError(f"Lambda1EvaluatorV2 failed with return code {completed.returncode}")

    try:
        metadata = validate_output(output_csv, len(rows))
    except (OSError, ValueError) as error:
        manifest["status"] = "failed"
        manifest["validation_error"] = str(error)
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        raise
    manifest.update(metadata)
    manifest["output_sha256"] = sha256(output_csv)
    manifest["status"] = "complete"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest
