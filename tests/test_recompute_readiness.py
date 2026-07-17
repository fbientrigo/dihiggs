from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys
from typing import cast

import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

pytest.importorskip(
    "scripts.recompute_readiness",
    reason="recompute-readiness is a frozen/quarantine component absent from this checkout",
)

from scripts.recompute_readiness import (
    APPROVED_SCHEMA_WHITELIST,
    BLOCKED_OUTPUT_COLUMNS,
    BLOCKED_ROWS_OUTPUT_RELATIVE_PATH,
    BLOCKER_REASON_CODES,
    CONTRACT_OUTPUT_RELATIVE_PATH,
    DEFAULT_QUARANTINE_ROOT,
    EXPLORER_CANONICAL_MAPPINGS,
    EXPLORER_MISMATCH,
    FILE_STATUS_BLOCKED_ALL,
    FILE_STATUS_EMPTY_DATA,
    FILE_STATUS_READY_ALL,
    FILE_STATUS_READY_PARTIAL,
    FILE_EMPTY,
    FILE_HEADER_DUPLICATE,
    FILE_HEADER_INVALID,
    FILE_NOT_FOUND,
    FILE_SCHEMA_MISMATCH,
    FILE_UNREADABLE,
    INVENTORY_LOOKUP_MISSING,
    OPTIONAL_AUDIT_COLUMNS,
    READY_LEGACY_AUDIT_COLUMNS,
    READY_OUTPUT_FIELDNAMES,
    READY_ROWS_OUTPUT_RELATIVE_PATH,
    REQUIRED_PHYSICS_INPUTS,
    REQUIRED_PROVENANCE_INPUTS,
    ROW_INCONSISTENT_ANGLE_REPRESENTATION,
    ROW_INVALID_FLOAT,
    ROW_MISSING_REQUIRED_FIELD,
    ROW_NONFINITE_REQUIRED_VALUE,
    SUPPORTED_EXPLORERS,
    UNSUPPORTED_EXPLORER,
    UNSUPPORTED_SCHEMA,
    build_file_gate_records,
    build_row_readiness_records,
    build_recompute_contract_v1,
    blocker_reason_for_schema,
    compute_schema_signature,
    compute_canonical_input_sha256,
    compute_source_row_sha256,
    format_float_canonical,
    normalize_value,

    contract_output_path,
    parse_args,
    write_recompute_artifacts,
    write_recompute_contract_v1,
)

APPROVED_HEADER_COLUMNS = [
    "m_phi",
    "mA",
    "alpha",
    "beta",
    "lambda6",
    "lambda7",
    "m12",
    "sin_ba",
    "tan_beta",
    "positivity_ok",
    "unitarity_ok",
    "perturbativity_ok",
    "width_bb",
    "width_tautau",
    "width_WW",
    "width_ZZ",
    "width_gaga",
    "width_Zga",
    "width_gg",
    "width_hh",
    "total_width",
    "br_gaga",
    "lam1",
    "computed_lam1",
    "lam2",
    "computed_lam2",
    "lam3",
    "lam4",
    "lam5",
]
APPROVED_SCHEMA_SIGNATURE = compute_schema_signature(APPROVED_HEADER_COLUMNS)
ALT_SCHEMA_SIGNATURE = compute_schema_signature([*APPROVED_HEADER_COLUMNS, "extra_col"])


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_header_only_csv(path: Path, header_columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    _ = path.write_text(",".join(header_columns) + "\n", encoding="utf-8")


def write_source_csv(path: Path, rows: list[dict[str, str]]) -> None:
    write_csv(path, APPROVED_HEADER_COLUMNS, rows)


def build_candidate_row(path: Path, *, schema_signature: str, explorer_guess: str) -> dict[str, str]:
    path_str = str(path)
    return {
        "path": path_str,
        "campaign": "campaign=test",
        "run_id": "run-test",
        "schema_signature": schema_signature,
        "explorer_guess": explorer_guess,
        "recompute_priority": "HIGH",
        "source_quarantine_path": path_str,
        "source_campaign": "campaign=test",
        "source_run": "run_dir=test",
        "recomputed_at": "",
        "calculator_git_sha": "",
    }


def build_inventory_row(path: Path, *, schema_signature: str, explorer_guess: str) -> dict[str, str]:
    path_str = str(path)
    return {
        "full_path": path_str,
        "relative_path": path.name,
        "campaign": "campaign=test",
        "run_dir": "run_dir=test",
        "run_id": "run-test",
        "tb_folder": "tb_10000",
        "filename": path.name,
        "size_bytes": "123",
        "mtime_epoch": "1710000000.0",
        "mtime_iso": "2026-04-17T00:00:00+00:00",
        "header_columns": json.dumps(APPROVED_HEADER_COLUMNS),
        "num_columns": str(len(APPROVED_HEADER_COLUMNS)),
        "schema_signature": schema_signature,
        "explorer_guess": explorer_guess,
        "copy_status": "unchanged",
        "mapping_error": "",
        "original_source_path": f"/lake/{path.name}",
    }


def write_candidate_and_inventory_inputs(
    quarantine_root: Path,
    *,
    candidate_rows: list[dict[str, str]],
    inventory_rows: list[dict[str, str]],
) -> None:
    write_csv(
        quarantine_root / "priority" / "recompute_candidates.csv",
        [
            "path",
            "campaign",
            "run_id",
            "schema_signature",
            "explorer_guess",
            "recompute_priority",
            "source_quarantine_path",
            "source_campaign",
            "source_run",
            "recomputed_at",
            "calculator_git_sha",
        ],
        candidate_rows,
    )
    write_csv(
        quarantine_root / "inventories" / "csv_inventory.csv",
        [
            "full_path",
            "relative_path",
            "campaign",
            "run_dir",
            "run_id",
            "tb_folder",
            "filename",
            "size_bytes",
            "mtime_epoch",
            "mtime_iso",
            "header_columns",
            "num_columns",
            "schema_signature",
            "explorer_guess",
            "copy_status",
            "mapping_error",
            "original_source_path",
        ],
        inventory_rows,
    )


def build_source_data_row(
    *,
    m_phi: str = "125.0",
    mA: str = "300.0",
    alpha: str = "0.2",
    beta: str = "0.6",
    lambda6: str = "0.01",
    lambda7: str = "0.02",
    m12: str = "1000.0",
    sin_ba: str | None = None,
    tan_beta: str | None = None,
    positivity_ok: str = "1",
    unitarity_ok: str = "1",
    perturbativity_ok: str = "1",
    width_bb: str = "0.1",
    width_tautau: str = "0.2",
    width_WW: str = "0.3",
    width_ZZ: str = "0.4",
    width_gaga: str = "0.5",
    width_Zga: str = "0.6",
    width_gg: str = "0.7",
    width_hh: str = "0.8",
    total_width: str = "1.9",
    br_gaga: str = "0.01",
    lam1: str = "1.1",
    computed_lam1: str = "1.1",
    lam2: str = "1.2",
    computed_lam2: str = "1.2",
    lam3: str = "1.3",
    lam4: str = "1.4",
    lam5: str = "1.5",
) -> dict[str, str]:
    alpha_value = float(alpha)
    beta_value = float(beta)
    return {
        "m_phi": m_phi,
        "mA": mA,
        "alpha": alpha,
        "beta": beta,
        "lambda6": lambda6,
        "lambda7": lambda7,
        "m12": m12,
        "sin_ba": sin_ba or format_float_canonical(math.sin(beta_value - alpha_value)),
        "tan_beta": tan_beta or format_float_canonical(math.tan(beta_value)),
        "positivity_ok": positivity_ok,
        "unitarity_ok": unitarity_ok,
        "perturbativity_ok": perturbativity_ok,
        "width_bb": width_bb,
        "width_tautau": width_tautau,
        "width_WW": width_WW,
        "width_ZZ": width_ZZ,
        "width_gaga": width_gaga,
        "width_Zga": width_Zga,
        "width_gg": width_gg,
        "width_hh": width_hh,
        "total_width": total_width,
        "br_gaga": br_gaga,
        "lam1": lam1,
        "computed_lam1": computed_lam1,
        "lam2": lam2,
        "computed_lam2": computed_lam2,
        "lam3": lam3,
        "lam4": lam4,
        "lam5": lam5,
    }


def row_records_for(file_record: dict[str, object]) -> list[dict[str, object]]:
    return cast(list[dict[str, object]], file_record["rows"])


def read_csv_output(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return cast(list[str], reader.fieldnames), [dict(row) for row in reader]


def test_contract_scope_is_locked() -> None:
    contract = build_recompute_contract_v1()

    assert SUPPORTED_EXPLORERS == ("tb_scan_explorer", "lam1_explorer")
    assert contract["supported_explorers"] == ["tb_scan_explorer", "lam1_explorer"]
    assert APPROVED_SCHEMA_WHITELIST == (
        "770adc466bf92cc9a65f39b84dd8bdd34be917553a11fff896aa2638bd6326bf",
    )
    assert contract["approved_schema_whitelist"] == [
        "770adc466bf92cc9a65f39b84dd8bdd34be917553a11fff896aa2638bd6326bf"
    ]


def test_contract_cli_accepts_quarantine_root_argument(tmp_path: Path) -> None:
    args = parse_args(["--quarantine-root", str(tmp_path)])
    assert cast(Path, args.quarantine_root) == tmp_path

    default_args = parse_args([])
    assert cast(Path, default_args.quarantine_root) == DEFAULT_QUARANTINE_ROOT


def test_contract_categories_keep_optional_audit_columns_separate() -> None:
    contract = build_recompute_contract_v1()
    categories = cast(dict[str, list[str]], contract["canonical_categories"])

    assert categories["required_physics_inputs"] == [
        "m_phi",
        "mA",
        "alpha",
        "beta",
        "sin_ba",
        "tan_beta",
        "lambda6",
        "lambda7",
        "m12",
    ]
    assert categories["required_provenance_inputs"] == [
        "source_quarantine_path",
        "original_source_path",
        "source_campaign",
        "source_run",
        "run_id",
        "tb_folder",
        "filename",
        "schema_signature",
        "explorer_guess",
        "recompute_priority",
        "source_row_index",
        "source_row_sha256",
        "canonical_input_sha256",
    ]
    assert categories["optional_audit_columns"] == list(OPTIONAL_AUDIT_COLUMNS)

    required_union = set(REQUIRED_PHYSICS_INPUTS) | set(REQUIRED_PROVENANCE_INPUTS)
    assert required_union.isdisjoint(OPTIONAL_AUDIT_COLUMNS)
    assert "positivity_ok" in OPTIONAL_AUDIT_COLUMNS
    assert "br_gaga" in OPTIONAL_AUDIT_COLUMNS
    assert "total_width" in OPTIONAL_AUDIT_COLUMNS


def test_contract_exposes_fixed_artifact_paths_and_columns() -> None:
    contract = build_recompute_contract_v1()

    assert contract["artifact_paths"] == {
        "contract": CONTRACT_OUTPUT_RELATIVE_PATH.as_posix(),
        "ready_rows": READY_ROWS_OUTPUT_RELATIVE_PATH.as_posix(),
        "blocked_rows": BLOCKED_ROWS_OUTPUT_RELATIVE_PATH.as_posix(),
    }
    assert contract["artifact_schemas"] == {
        "ready_rows": list(READY_OUTPUT_FIELDNAMES),
        "blocked_rows": list(BLOCKED_OUTPUT_COLUMNS),
    }


def test_reason_codes_match_contract_taxonomy() -> None:
    assert BLOCKER_REASON_CODES == (
        "INVENTORY_LOOKUP_MISSING",
        "FILE_NOT_FOUND",
        "FILE_UNREADABLE",
        "FILE_EMPTY",
        "FILE_HEADER_INVALID",
        "FILE_HEADER_DUPLICATE",
        "FILE_SCHEMA_MISMATCH",
        "UNSUPPORTED_SCHEMA",
        "UNSUPPORTED_EXPLORER",
        "EXPLORER_MISMATCH",
        "ROW_MISSING_REQUIRED_FIELD",
        "ROW_INVALID_FLOAT",
        "ROW_NONFINITE_REQUIRED_VALUE",
        "ROW_INCONSISTENT_ANGLE_REPRESENTATION",
        "ROW_DUPLICATE_SOURCE_ID",
    )


def test_unsupported_schema_maps_to_reason_code() -> None:
    assert (
        blocker_reason_for_schema(
            "e21e4f3bcc27d7889acd93c98890f0bdcfde7cf010753ee0ff04e9cb683292f4"
        )
        == UNSUPPORTED_SCHEMA
    )
    assert (
        blocker_reason_for_schema(
            "770adc466bf92cc9a65f39b84dd8bdd34be917553a11fff896aa2638bd6326bf"
        )
        is None
    )


def test_contract_writer_is_deterministic_and_has_no_timestamps(tmp_path: Path) -> None:
    expected_path = tmp_path / CONTRACT_OUTPUT_RELATIVE_PATH

    assert contract_output_path(tmp_path) == expected_path

    written_path = write_recompute_contract_v1(tmp_path)
    assert written_path == expected_path

    payload = cast(dict[str, object], json.loads(written_path.read_text(encoding="utf-8")))
    assert "created_at" not in payload
    assert "timestamp" not in payload
    assert payload == build_recompute_contract_v1()

def test_stable_hashes_and_normalization() -> None:
    # Test float formatting
    assert format_float_canonical(1.23456789012345678) == "1.2345678901234567"
    assert format_float_canonical(float("nan")) == "nan"
    assert format_float_canonical(float("inf")) == "inf"
    assert format_float_canonical(float("-inf")) == "-inf"

    # Test normalization
    assert normalize_value("  trimmed  ") == "trimmed"
    assert normalize_value(1.23456789012345678) == "1.2345678901234567"
    assert normalize_value(None) == ""

    # Test source row hashing (deterministic with order)
    row = {"a": 1.0, "b": " val "}
    h1 = compute_source_row_sha256(row, ["a", "b"])
    h2 = compute_source_row_sha256(row, ["b", "a"])
    assert h1 != h2
    assert h1 == compute_source_row_sha256({"a": 1.0, "b": "val"}, ["a", "b"])

    # Test canonical hashing (uses REQUIRED_PHYSICS_INPUTS order)
    physics_row = {k: 1.0 for k in REQUIRED_PHYSICS_INPUTS}
    ch1 = compute_canonical_input_sha256(physics_row)
    # Different order in dict should NOT change hash
    shuffled_row = {k: 1.0 for k in reversed(list(REQUIRED_PHYSICS_INPUTS))}
    ch2 = compute_canonical_input_sha256(shuffled_row)
    assert ch1 == ch2


def test_duplicate_physics_preserved_by_identity() -> None:
    # Same physics, different provenance
    physics = {k: 1.0 for k in REQUIRED_PHYSICS_INPUTS}
    row1 = {**physics, "original_source_path": "file1.csv", "source_row_index": 0}
    row2 = {**physics, "original_source_path": "file2.csv", "source_row_index": 0}

    # Canonical hashes must be identical
    assert compute_canonical_input_sha256(row1) == compute_canonical_input_sha256(row2)

    # Source row hashes must be different (provenance differs)
    cols = ["original_source_path", "source_row_index"] + list(REQUIRED_PHYSICS_INPUTS)
    h1 = compute_source_row_sha256(row1, cols)
    h2 = compute_source_row_sha256(row2, cols)
    assert h1 != h2


def test_inventory_join_carries_metadata_and_sorts_by_path(tmp_path: Path) -> None:
    a_path = tmp_path / "raw" / "a.csv"
    b_path = tmp_path / "raw" / "b.csv"
    write_header_only_csv(a_path, APPROVED_HEADER_COLUMNS)
    write_header_only_csv(b_path, APPROVED_HEADER_COLUMNS)

    candidate_rows = [
        build_candidate_row(b_path, schema_signature=APPROVED_SCHEMA_SIGNATURE, explorer_guess="tb_scan_explorer"),
        build_candidate_row(a_path, schema_signature=APPROVED_SCHEMA_SIGNATURE, explorer_guess="tb_scan_explorer"),
    ]
    inventory_rows = [
        build_inventory_row(b_path, schema_signature=APPROVED_SCHEMA_SIGNATURE, explorer_guess="tb_scan_explorer"),
        build_inventory_row(a_path, schema_signature=APPROVED_SCHEMA_SIGNATURE, explorer_guess="tb_scan_explorer"),
    ]
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=candidate_rows,
        inventory_rows=inventory_rows,
    )

    records = build_file_gate_records(tmp_path)

    assert [record["path"] for record in records] == [str(a_path), str(b_path)]
    assert records[0]["original_source_path"] == "/lake/a.csv"
    assert records[0]["copy_status"] == "unchanged"
    assert records[0]["mtime_epoch"] == "1710000000.0"
    assert records[0]["size_bytes"] == "123"
    assert records[0]["campaign"] == "campaign=test"
    assert records[0]["run_dir"] == "run_dir=test"
    assert records[0]["run_id"] == "run-test"
    assert records[0]["tb_folder"] == "tb_10000"
    assert records[0]["filename"] == "a.csv"
    assert records[0]["blocker_reason"] is None


def test_inventory_join_requires_exact_path_equality(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "candidate.csv"
    write_header_only_csv(csv_path, APPROVED_HEADER_COLUMNS)
    candidate_row = build_candidate_row(
        csv_path,
        schema_signature=APPROVED_SCHEMA_SIGNATURE,
        explorer_guess="tb_scan_explorer",
    )
    candidate_row["source_quarantine_path"] = str(tmp_path / "raw" / "different.csv")
    inventory_row = build_inventory_row(
        csv_path,
        schema_signature=APPROVED_SCHEMA_SIGNATURE,
        explorer_guess="tb_scan_explorer",
    )
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[candidate_row],
        inventory_rows=[inventory_row],
    )

    records = build_file_gate_records(tmp_path)

    assert records[0]["blocker_reason"] == INVENTORY_LOOKUP_MISSING


def test_file_gate_missing_inventory_maps_reason(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "candidate.csv"
    write_header_only_csv(csv_path, APPROVED_HEADER_COLUMNS)
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
        inventory_rows=[],
    )

    records = build_file_gate_records(tmp_path)
    assert records[0]["blocker_reason"] == INVENTORY_LOOKUP_MISSING


def test_file_gate_missing_file_maps_reason(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "missing.csv"
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
        inventory_rows=[
            build_inventory_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
    )

    records = build_file_gate_records(tmp_path)
    assert records[0]["blocker_reason"] == FILE_NOT_FOUND


def test_file_gate_unreadable_file_maps_reason(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "unreadable.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    _ = csv_path.write_bytes(b"\xff\xfe\xfd")
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
        inventory_rows=[
            build_inventory_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
    )

    records = build_file_gate_records(tmp_path)
    assert records[0]["blocker_reason"] == FILE_UNREADABLE


def test_file_gate_empty_file_maps_reason(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "empty.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    _ = csv_path.write_text("", encoding="utf-8")
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
        inventory_rows=[
            build_inventory_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
    )

    records = build_file_gate_records(tmp_path)
    assert records[0]["blocker_reason"] == FILE_EMPTY


def test_file_gate_blank_header_column_maps_reason(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "invalid_header.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    _ = csv_path.write_text("m_phi,,alpha\n", encoding="utf-8")
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
        inventory_rows=[
            build_inventory_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
    )

    records = build_file_gate_records(tmp_path)
    assert records[0]["blocker_reason"] == FILE_HEADER_INVALID


def test_file_gate_duplicate_header_maps_reason(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "duplicate_header.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    _ = csv_path.write_text("m_phi,m_phi,alpha\n", encoding="utf-8")
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
        inventory_rows=[
            build_inventory_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
    )

    records = build_file_gate_records(tmp_path)
    assert records[0]["blocker_reason"] == FILE_HEADER_DUPLICATE


def test_schema_mismatch_detects_candidate_inventory_observed_disagreement(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "schema_mismatch.csv"
    write_header_only_csv(csv_path, APPROVED_HEADER_COLUMNS)
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                csv_path,
                schema_signature=ALT_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
        inventory_rows=[
            build_inventory_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
    )

    records = build_file_gate_records(tmp_path)
    assert records[0]["observed_schema_signature"] == APPROVED_SCHEMA_SIGNATURE
    assert records[0]["blocker_reason"] == FILE_SCHEMA_MISMATCH


def test_file_gate_unsupported_schema_maps_reason(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "unsupported_schema.csv"
    write_header_only_csv(csv_path, APPROVED_HEADER_COLUMNS + ["extra_col"])
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                csv_path,
                schema_signature=ALT_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
        inventory_rows=[
            build_inventory_row(
                csv_path,
                schema_signature=ALT_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
    )

    records = build_file_gate_records(tmp_path)
    assert records[0]["blocker_reason"] == UNSUPPORTED_SCHEMA


def test_unsupported_explorer_maps_reason(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "unsupported_explorer.csv"
    write_header_only_csv(csv_path, APPROVED_HEADER_COLUMNS)
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="mystery_explorer",
            )
        ],
        inventory_rows=[
            build_inventory_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="mystery_explorer",
            )
        ],
    )

    records = build_file_gate_records(tmp_path)
    assert records[0]["blocker_reason"] == UNSUPPORTED_EXPLORER


def test_file_gate_explorer_mismatch_maps_reason(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "explorer_mismatch.csv"
    write_header_only_csv(csv_path, APPROVED_HEADER_COLUMNS)
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
        inventory_rows=[
            build_inventory_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="lam1_explorer",
            )
        ],
    )

    records = build_file_gate_records(tmp_path)
    assert records[0]["blocker_reason"] == EXPLORER_MISMATCH


def test_row_ready_identity_mapping_and_legacy_fields(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "ready.csv"
    source_row = build_source_data_row(width_bb="0.123")
    write_source_csv(csv_path, [source_row])
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
        inventory_rows=[
            build_inventory_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
    )

    records = build_row_readiness_records(tmp_path)

    assert EXPLORER_CANONICAL_MAPPINGS["tb_scan_explorer"] == {
        field: field for field in REQUIRED_PHYSICS_INPUTS
    }
    assert EXPLORER_CANONICAL_MAPPINGS["lam1_explorer"] == {
        field: field for field in REQUIRED_PHYSICS_INPUTS
    }
    assert records[0]["file_status"] == FILE_STATUS_READY_ALL
    row = row_records_for(records[0])[0]
    assert row["source_row_index"] == 1
    assert row["blocker_reason"] is None
    assert row["m_phi"] == 125.0
    assert row["legacy_width_bb"] == "0.123"
    assert row["legacy_positivity_ok"] == "1"
    assert row["canonical_input_sha256"] == compute_canonical_input_sha256(
        {field: row[field] for field in REQUIRED_PHYSICS_INPUTS}
    )


def test_row_blocked_missing_required_field(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "missing_required.csv"
    write_source_csv(csv_path, [build_source_data_row(m12="")])
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
        inventory_rows=[
            build_inventory_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
    )

    records = build_row_readiness_records(tmp_path)

    assert records[0]["file_status"] == FILE_STATUS_BLOCKED_ALL
    row = row_records_for(records[0])[0]
    assert row["source_row_index"] == 1
    assert row["blocker_reason"] == ROW_MISSING_REQUIRED_FIELD
    assert row["canonical_input_sha256"] == ""


def test_row_blocked_invalid_float_and_nonfinite(tmp_path: Path) -> None:
    invalid_path = tmp_path / "raw" / "invalid_float.csv"
    nonfinite_path = tmp_path / "raw" / "nonfinite.csv"
    write_source_csv(invalid_path, [build_source_data_row(lambda6="not-a-float")])
    write_source_csv(nonfinite_path, [build_source_data_row(lambda7="nan")])
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                invalid_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            ),
            build_candidate_row(
                nonfinite_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="lam1_explorer",
            ),
        ],
        inventory_rows=[
            build_inventory_row(
                invalid_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            ),
            build_inventory_row(
                nonfinite_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="lam1_explorer",
            ),
        ],
    )

    records = build_row_readiness_records(tmp_path)

    assert [record["path"] for record in records] == [str(invalid_path), str(nonfinite_path)]
    assert row_records_for(records[0])[0]["blocker_reason"] == ROW_INVALID_FLOAT
    assert row_records_for(records[1])[0]["blocker_reason"] == ROW_NONFINITE_REQUIRED_VALUE


def test_mixed_quality_rows_drive_ready_partial_and_row_order(tmp_path: Path) -> None:
    csv_path = tmp_path / "raw" / "mixed.csv"
    write_source_csv(
        csv_path,
        [
            build_source_data_row(m_phi="126.0"),
            build_source_data_row(mA="oops"),
            build_source_data_row(m_phi="127.0"),
        ],
    )
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
        inventory_rows=[
            build_inventory_row(
                csv_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
    )

    records = build_row_readiness_records(tmp_path)

    assert records[0]["file_status"] == FILE_STATUS_READY_PARTIAL
    rows = row_records_for(records[0])
    assert [row["source_row_index"] for row in rows] == [1, 2, 3]
    assert [row["blocker_reason"] for row in rows] == [
        None,
        ROW_INVALID_FLOAT,
        None,
    ]


def test_angle_consistency_blocker_and_empty_data_status(tmp_path: Path) -> None:
    inconsistent_path = tmp_path / "raw" / "angle_inconsistent.csv"
    empty_path = tmp_path / "raw" / "header_only.csv"
    write_source_csv(
        inconsistent_path,
        [build_source_data_row(sin_ba="0.999999", tan_beta="0.888888")],
    )
    write_header_only_csv(empty_path, APPROVED_HEADER_COLUMNS)
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                inconsistent_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            ),
            build_candidate_row(
                empty_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            ),
        ],
        inventory_rows=[
            build_inventory_row(
                inconsistent_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            ),
            build_inventory_row(
                empty_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            ),
        ],
    )

    records = build_row_readiness_records(tmp_path)

    assert records[0]["file_status"] == FILE_STATUS_BLOCKED_ALL
    inconsistent_row = row_records_for(records[0])[0]
    assert inconsistent_row["blocker_reason"] == ROW_INCONSISTENT_ANGLE_REPRESENTATION
    assert inconsistent_row["canonical_input_sha256"]
    assert records[1]["file_status"] == FILE_STATUS_EMPTY_DATA
    assert row_records_for(records[1]) == []


def test_artifact_columns_and_row_reconciliation(tmp_path: Path) -> None:
    ready_path = tmp_path / "raw" / "ready.csv"
    mixed_path = tmp_path / "raw" / "mixed.csv"
    write_source_csv(ready_path, [build_source_data_row(width_bb="0.123")])
    write_source_csv(
        mixed_path,
        [
            build_source_data_row(m_phi="126.0"),
            build_source_data_row(m12=""),
        ],
    )
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                mixed_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            ),
            build_candidate_row(
                ready_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            ),
        ],
        inventory_rows=[
            build_inventory_row(
                mixed_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            ),
            build_inventory_row(
                ready_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            ),
        ],
    )

    outputs = write_recompute_artifacts(tmp_path)
    assert outputs["contract"] == tmp_path / CONTRACT_OUTPUT_RELATIVE_PATH
    assert outputs["ready_rows"] == tmp_path / READY_ROWS_OUTPUT_RELATIVE_PATH
    assert outputs["blocked_rows"] == tmp_path / BLOCKED_ROWS_OUTPUT_RELATIVE_PATH

    ready_columns, ready_rows = read_csv_output(outputs["ready_rows"])
    blocked_columns, blocked_rows = read_csv_output(outputs["blocked_rows"])

    assert ready_columns == list(READY_OUTPUT_FIELDNAMES)
    assert blocked_columns == list(BLOCKED_OUTPUT_COLUMNS)

    assert [(row["source_quarantine_path"], row["source_row_index"]) for row in ready_rows] == [
        (str(mixed_path), "1"),
        (str(ready_path), "1"),
    ]
    assert [row["readiness_status"] for row in ready_rows] == ["READY", "READY"]
    assert [row["file_readiness_status"] for row in ready_rows] == [FILE_STATUS_READY_PARTIAL, FILE_STATUS_READY_ALL]
    assert READY_LEGACY_AUDIT_COLUMNS == tuple(sorted(READY_LEGACY_AUDIT_COLUMNS))
    assert ready_rows[0]["source_row_sha256"]
    assert ready_rows[0]["canonical_input_sha256"]
    assert ready_rows[1]["source_row_sha256"]
    assert ready_rows[1]["canonical_input_sha256"]
    assert ready_rows[0]["legacy_width_bb"] == "0.1"
    assert ready_rows[1]["legacy_width_bb"] == "0.123"

    assert len(blocked_rows) == 1
    blocked_row = blocked_rows[0]
    assert blocked_row["readiness_scope"] == "ROW"
    assert blocked_row["file_readiness_status"] == FILE_STATUS_READY_PARTIAL
    assert blocked_row["row_readiness_status"] == "BLOCKED"
    assert blocked_row["blocked_reason_code"] == ROW_MISSING_REQUIRED_FIELD
    assert blocked_row["blocked_reason_detail"] == "m12"
    assert blocked_row["missing_required_fields"] == "m12"
    assert blocked_row["candidate_schema_signature"] == APPROVED_SCHEMA_SIGNATURE
    assert blocked_row["observed_schema_signature"] == APPROVED_SCHEMA_SIGNATURE
    assert blocked_row["candidate_explorer_guess"] == "tb_scan_explorer"
    assert blocked_row["inventory_explorer_guess"] == "tb_scan_explorer"


def test_blocked_rows_keep_file_scope_blockers_with_blank_row_identity(tmp_path: Path) -> None:
    missing_path = tmp_path / "raw" / "missing.csv"
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                missing_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
        inventory_rows=[
            build_inventory_row(
                missing_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            )
        ],
    )

    outputs = write_recompute_artifacts(tmp_path)
    ready_columns, ready_rows = read_csv_output(outputs["ready_rows"])
    blocked_columns, blocked_rows = read_csv_output(outputs["blocked_rows"])

    assert ready_columns == list(READY_OUTPUT_FIELDNAMES)
    assert ready_rows == []
    assert blocked_columns == list(BLOCKED_OUTPUT_COLUMNS)
    assert len(blocked_rows) == 1

    blocked_row = blocked_rows[0]
    assert blocked_row["readiness_scope"] == "FILE"
    assert blocked_row["file_readiness_status"] == FILE_STATUS_BLOCKED_ALL
    assert blocked_row["row_readiness_status"] == ""
    assert blocked_row["blocked_reason_code"] == FILE_NOT_FOUND
    assert blocked_row["blocked_reason_detail"] == f"missing file: {missing_path}"
    assert blocked_row["source_row_index"] == ""
    assert blocked_row["source_row_sha256"] == ""
    assert blocked_row["canonical_input_sha256"] == ""
    assert blocked_row["source_quarantine_path"] == str(missing_path)
    assert blocked_row["original_source_path"] == "/lake/missing.csv"


def test_deterministic_outputs_across_repeated_runs(tmp_path: Path) -> None:
    ready_path = tmp_path / "raw" / "repeatable.csv"
    blocked_path = tmp_path / "raw" / "repeatable_blocked.csv"
    write_source_csv(ready_path, [build_source_data_row(), build_source_data_row(m_phi="126.0")])
    write_source_csv(blocked_path, [build_source_data_row(lambda6="oops")])
    write_candidate_and_inventory_inputs(
        tmp_path,
        candidate_rows=[
            build_candidate_row(
                blocked_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="lam1_explorer",
            ),
            build_candidate_row(
                ready_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            ),
        ],
        inventory_rows=[
            build_inventory_row(
                blocked_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="lam1_explorer",
            ),
            build_inventory_row(
                ready_path,
                schema_signature=APPROVED_SCHEMA_SIGNATURE,
                explorer_guess="tb_scan_explorer",
            ),
        ],
    )

    first = write_recompute_artifacts(tmp_path)
    first_contract = first["contract"].read_text(encoding="utf-8")
    first_ready = first["ready_rows"].read_text(encoding="utf-8")
    first_blocked = first["blocked_rows"].read_text(encoding="utf-8")

    second = write_recompute_artifacts(tmp_path)
    assert second["contract"].read_text(encoding="utf-8") == first_contract
    assert second["ready_rows"].read_text(encoding="utf-8") == first_ready
    assert second["blocked_rows"].read_text(encoding="utf-8") == first_blocked
