from __future__ import annotations

import argparse
import hashlib
import json
import os
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

DEFAULT_DATALAKE_ROOT = Path("/mnt/c/Users/Asus/cern_db/dihiggs_lake")
REQUIRED_ARTIFACTS: tuple[str, ...] = ("events.jsonl", "task_summary.jsonl")
OPTIONAL_ARTIFACT_NAMES: tuple[str, ...] = ("orchestrator.log",)
STATUS_SUFFIXES: tuple[str, ...] = ("status.json", "state.json")
MAX_JSONL_INSPECTION_LINES = 2000
MAX_LOG_INSPECTION_LINES = 2000


@dataclass(frozen=True)
class ManifestRow:
    campaign_dir: str
    artifact_type: str
    artifact_name: str
    source_path: str
    source_mtime: str | None
    parser_status: str
    records_read: int
    records_parsed: int
    records_quarantined: int
    detail: str
    run_dir_fingerprint: str


def _mtime_to_iso(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat().replace("+00:00", "Z")


def _fingerprint_path(path: Path) -> str:
    return hashlib.sha256(str(path).encode("utf-8")).hexdigest()[:16]


def _row(
    *,
    campaign_dir: Path,
    artifact_type: str,
    artifact_name: str,
    source_path: Path,
    source_mtime: str | None,
    parser_status: str,
    records_read: int,
    records_parsed: int,
    records_quarantined: int,
    detail: str,
) -> ManifestRow:
    return ManifestRow(
        campaign_dir=str(campaign_dir),
        artifact_type=artifact_type,
        artifact_name=artifact_name,
        source_path=str(source_path),
        source_mtime=source_mtime,
        parser_status=parser_status,
        records_read=records_read,
        records_parsed=records_parsed,
        records_quarantined=records_quarantined,
        detail=detail,
        run_dir_fingerprint=_fingerprint_path(campaign_dir),
    )


def _scan_jsonl_artifact(campaign_dir: Path, artifact_path: Path, artifact_type: str) -> ManifestRow:
    if not artifact_path.exists():
        return _row(
            campaign_dir=campaign_dir,
            artifact_type=artifact_type,
            artifact_name=artifact_path.name,
            source_path=artifact_path,
            source_mtime=None,
            parser_status="missing",
            records_read=0,
            records_parsed=0,
            records_quarantined=0,
            detail=f"required artifact missing: {artifact_path.name}",
        )

    records_read = 0
    records_parsed = 0
    records_quarantined = 0
    detail = "parsed successfully"
    parser_status = "success"

    try:
        with artifact_path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                if records_read >= MAX_JSONL_INSPECTION_LINES:
                    if parser_status == "success":
                        detail = f"inspection truncated at {MAX_JSONL_INSPECTION_LINES} lines"
                    else:
                        detail = f"{detail}; inspection truncated at {MAX_JSONL_INSPECTION_LINES} lines"
                    break
                records_read += 1
                line = raw_line.strip()
                if not line:
                    continue
                try:
                    record = json.loads(line)
                except json.JSONDecodeError as exc:
                    records_quarantined += 1
                    parser_status = "quarantined_malformed"
                    detail = f"malformed JSONL record: {exc}"
                    continue

                if isinstance(record, dict):
                    records_parsed += 1
                else:
                    records_quarantined += 1
                    parser_status = "quarantined_schema"
                    detail = "non-mapping JSONL record encountered"
    except OSError as exc:
        return _row(
            campaign_dir=campaign_dir,
            artifact_type=artifact_type,
            artifact_name=artifact_path.name,
            source_path=artifact_path,
            source_mtime=None,
            parser_status="quarantined_unreadable",
            records_read=records_read,
            records_parsed=records_parsed,
            records_quarantined=records_quarantined,
            detail=str(exc),
        )

    return _row(
        campaign_dir=campaign_dir,
        artifact_type=artifact_type,
        artifact_name=artifact_path.name,
        source_path=artifact_path,
        source_mtime=_mtime_to_iso(artifact_path),
        parser_status=parser_status,
        records_read=records_read,
        records_parsed=records_parsed,
        records_quarantined=records_quarantined,
        detail=detail,
    )


def _scan_log_artifact(campaign_dir: Path, artifact_path: Path) -> ManifestRow:
    records_read = 0
    try:
        with artifact_path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                if records_read >= MAX_LOG_INSPECTION_LINES:
                    break
                records_read += 1
                raw_line.strip()
    except OSError as exc:
        return _row(
            campaign_dir=campaign_dir,
            artifact_type="log",
            artifact_name=artifact_path.name,
            source_path=artifact_path,
            source_mtime=None,
            parser_status="quarantined_unreadable",
            records_read=records_read,
            records_parsed=0,
            records_quarantined=0,
            detail=str(exc),
        )

    return _row(
        campaign_dir=campaign_dir,
        artifact_type="log",
        artifact_name=artifact_path.name,
        source_path=artifact_path,
        source_mtime=_mtime_to_iso(artifact_path),
        parser_status="success",
        records_read=records_read,
        records_parsed=records_read,
        records_quarantined=0,
        detail=(
            f"inspection truncated at {MAX_LOG_INSPECTION_LINES} lines"
            if records_read >= MAX_LOG_INSPECTION_LINES
            else "parsed successfully"
        ),
    )


def _scan_status_artifact(campaign_dir: Path, artifact_path: Path) -> ManifestRow:
    try:
        text = artifact_path.read_text(encoding="utf-8")
    except OSError as exc:
        return _row(
            campaign_dir=campaign_dir,
            artifact_type="status",
            artifact_name=artifact_path.name,
            source_path=artifact_path,
            source_mtime=None,
            parser_status="quarantined_unreadable",
            records_read=0,
            records_parsed=0,
            records_quarantined=0,
            detail=str(exc),
        )

    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        return _row(
            campaign_dir=campaign_dir,
            artifact_type="status",
            artifact_name=artifact_path.name,
            source_path=artifact_path,
            source_mtime=_mtime_to_iso(artifact_path),
            parser_status="quarantined_malformed",
            records_read=1,
            records_parsed=0,
            records_quarantined=1,
            detail=f"malformed status JSON: {exc}",
        )

    detail = "parsed successfully" if isinstance(payload, (dict, list)) else "scalar status payload"
    return _row(
        campaign_dir=campaign_dir,
        artifact_type="status",
        artifact_name=artifact_path.name,
        source_path=artifact_path,
        source_mtime=_mtime_to_iso(artifact_path),
        parser_status="success",
        records_read=1,
        records_parsed=1,
        records_quarantined=0,
        detail=detail,
    )


def discover_campaign_dirs(root: str | Path) -> list[Path]:
    root_path = Path(root)
    if not root_path.exists() or not root_path.is_dir():
        return []

    return sorted(
        [entry for entry in root_path.iterdir() if entry.is_dir()],
        key=lambda item: str(item),
    )


def scan_campaign_artifacts(campaign_dir: str | Path) -> list[dict[str, Any]]:
    campaign_path = Path(campaign_dir)
    rows: list[ManifestRow] = []
    discovered_required: dict[str, list[Path]] = {name: [] for name in REQUIRED_ARTIFACTS}
    optional_paths: list[Path] = []
    status_paths: list[Path] = []

    for dirpath, _, filenames in os.walk(campaign_path, topdown=True):
        current_dir = Path(dirpath)
        for filename in sorted(filenames):
            file_path = current_dir / filename
            if filename in discovered_required:
                discovered_required[filename].append(file_path)
            elif filename in OPTIONAL_ARTIFACT_NAMES:
                optional_paths.append(file_path)
            elif any(filename.endswith(suffix) for suffix in STATUS_SUFFIXES):
                status_paths.append(file_path)

    for artifact_name, artifact_type in (("events.jsonl", "events"), ("task_summary.jsonl", "task_summary")):
        paths = sorted(discovered_required[artifact_name], key=lambda item: str(item))
        if not paths:
            rows.append(_scan_jsonl_artifact(campaign_path, campaign_path / artifact_name, artifact_type))
            continue
        for artifact_path in paths:
            rows.append(_scan_jsonl_artifact(campaign_path, artifact_path, artifact_type))

    for artifact_path in sorted(optional_paths, key=lambda item: str(item)):
        rows.append(_scan_log_artifact(campaign_path, artifact_path))

    for artifact_path in sorted(status_paths, key=lambda item: str(item)):
        rows.append(_scan_status_artifact(campaign_path, artifact_path))

    return [asdict(row) for row in rows]


def scan_datalake(root: str | Path = DEFAULT_DATALAKE_ROOT) -> dict[str, Any]:
    root_path = Path(root)
    if not root_path.exists() or not root_path.is_dir():
        return {
            "root_path": str(root_path),
            "status": "root_not_found",
            "campaign_count": 0,
            "manifest_row_count": 0,
            "parser_status_counts": {"root_not_found": 1},
            "manifest": [],
        }

    campaigns = discover_campaign_dirs(root_path)
    manifest: list[dict[str, Any]] = []
    for campaign_dir in campaigns:
        manifest.extend(scan_campaign_artifacts(campaign_dir))

    parser_status_counts: dict[str, int] = {}
    for row in manifest:
        status = str(row.get("parser_status", "unknown"))
        parser_status_counts[status] = parser_status_counts.get(status, 0) + 1

    status = "success"
    if any(key.startswith("quarantined") for key in parser_status_counts) or parser_status_counts.get("missing"):
        status = "success_with_quarantine"

    return {
        "root_path": str(root_path),
        "status": status,
        "campaign_count": len(campaigns),
        "manifest_row_count": len(manifest),
        "parser_status_counts": dict(sorted(parser_status_counts.items())),
        "manifest": manifest,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Scan datalake artifacts into a normalized manifest.")
    parser.add_argument("root", nargs="?", default=str(DEFAULT_DATALAKE_ROOT))
    parser.add_argument("--pretty", action="store_true")
    args = parser.parse_args(argv)

    report = scan_datalake(args.root)
    if args.pretty:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(json.dumps(report, sort_keys=True))
    return 0 if report["status"] != "root_not_found" else 1


if __name__ == "__main__":
    raise SystemExit(main())
