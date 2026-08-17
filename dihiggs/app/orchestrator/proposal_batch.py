"""proposal_batch.py — thin `proposals.csv -> DihiggsPointV2Evaluator` handoff.

Implements the socket described in GitHub issue #72: an arbitrary batch of
proposed 2HDM points is submitted, one at a time, to the canonical
``DihiggsPointV2Evaluator`` and every attempt (malformed, evaluator-error, or
evaluated) is preserved as one row of ``attempts.csv``.

This module owns none of the physics. It never reconstructs the 2HDM, never
recomputes a width/BR/coupling, and never decides which proposals are
"interesting" — see docs/contracts/adaptive_proposal_batch_v1.yaml for the
input contract this module validates against, and
docs/contracts/canonical_evaluators_v2.md for the evaluator it wraps.

Input contract
--------------
``proposals.csv`` header: the columns in ``REQUIRED_COLUMNS`` (any order),
plus the optional ``mh_GeV`` column. No other columns are permitted. See the
YAML contract for full field semantics.

Output: ``dihiggs.adaptive_attempt.v1``
----------------------------------------
``attempts.csv`` has one row per input proposal row, in the same order.
Columns are ``ATTEMPT_ENVELOPE_COLUMNS`` followed by the canonical
``dihiggs.point.v2`` columns, copied verbatim (as strings, unparsed) from the
evaluator's own output row. Rows that never reached the evaluator (MALFORMED)
or whose evaluator invocation could not be trusted (EVALUATOR_ERROR) carry
empty strings in every canonical column — the canonical schema itself is
never mutated or reformatted.

Validation split
----------------
This module checks *formatting* validity of a proposal (does it parse, is
yukawa_type in {1,2,3,4}, does mh_GeV match the pinned convention, is
proposal_id unique). It deliberately does NOT check *physical* validity
(tan_beta > 0, |sin_ba| <= 1, mh < mH, ...) — those are 2HDM theory
predicates owned exclusively by ``THDM::set_param_phys`` inside
``DihiggsPointV2Evaluator``. A row that is well-formed here but physically
invalid still reaches the evaluator subprocess and is recorded however the
evaluator reports it (construction_ok=0 as a normal EVALUATED row, or a
nonzero evaluator exit as EVALUATOR_ERROR for CLI-level guards such as
``invalid_physics_input``).
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from dihiggs.app.orchestrator.io_utils import detect_git_info, safe_write_json, utc_now_iso
from dihiggs.app.orchestrator.lambda1_v2 import _twohdmc_provenance, format_float64, sha256

PROPOSAL_SCHEMA_VERSION = "dihiggs.adaptive_proposal_batch.v1"
ATTEMPT_SCHEMA_VERSION = "dihiggs.adaptive_attempt.v1"

REQUIRED_COLUMNS = (
    "proposal_id", "mH_GeV", "mA_GeV", "mHp_GeV", "M2_GeV2", "tan_beta",
    "lambda6", "lambda7", "sin_beta_minus_alpha", "yukawa_type",
)
OPTIONAL_COLUMNS = ("mh_GeV",)
ALL_INPUT_COLUMNS = frozenset(REQUIRED_COLUMNS) | frozenset(OPTIONAL_COLUMNS)

# Active project convention (docs/HIGH_MASS_H2_CONTRACT.md,
# "Active project/high-mass convention: m_h = 125.20 GeV"). Never overridable from this wrapper's CLI.
MH_CONVENTION_GEV = 125.20
MH_CONVENTION_SOURCE = "docs/HIGH_MASS_H2_CONTRACT.md: Active project convention m_h = 125.20 GeV"
SUPPORTED_YUKAWA_TYPES = (1, 2, 3, 4)

ATTEMPT_ENVELOPE_COLUMNS = (
    "row_index", "proposal_id", "attempt_status", "attempt_stage",
    "attempt_reason", "mh_GeV_used",
)

# The exact dihiggs.point.v2 header (tests/test_dihiggs_point_v2.py::HEADER).
# Used verbatim as the attempts.csv canonical-column block ONLY when no
# proposal in the batch ever reached a successful evaluator invocation (so
# there is no observed header to copy); otherwise the real header wins.
CANONICAL_HEADER_FALLBACK = (
    "schema_version,producer,producer_commit,producer_dirty,evaluator_api,campaign_id,run_id,point_id,"
    "yukawa_type,yukawa_type_installed,mh_input_GeV,mH_input_GeV,mA_input_GeV,mHp_input_GeV,sin_beta_minus_alpha_input,"
    "tan_beta_input,beta_input_rad,M2_input_GeV2,m12_sq_input_GeV2,lambda6_input,lambda7_input,"
    "lambda1_reconstructed,lambda2_reconstructed,lambda3_reconstructed,lambda4_reconstructed,"
    "lambda5_reconstructed,lambda6_reconstructed,lambda7_reconstructed,tan_beta_reconstructed,"
    "m12_sq_reconstructed_GeV2,M2_reconstructed_GeV2,construction_ok,numerical_ok,rejection_stage,rejection_reason,"
    "positivity_reported_ok,unitarity_ok,perturbativity_ok,stability_reported_ok,stability_dependency_alias,"
    "triple_ok_legacy,theory_ok_v1,experimental_evaluated,experimental_ok,g_hH2H2_GeV,width_bb_GeV,width_cc_GeV,"
    "width_tt_GeV,width_tautau_GeV,width_WW_GeV,width_ZZ_GeV,width_gammagamma_GeV,width_Zgamma_GeV,width_gg_GeV,"
    "width_hh_GeV,total_width_GeV,width_unaccounted_GeV,br_bb,br_cc,br_tt,br_tautau,br_WW,br_ZZ,br_gammagamma,br_Zgamma,br_gg,"
    "br_hh,width_ok,ctau_mm"
).split(",")

# Fail-closed presence check (Task 3): a header missing any of these can
# never be recorded as a successful batch, even if the subprocess exit code
# was 0. g_hH2H2_GeV is explicitly required per the task instructions.
REQUIRED_CANONICAL_COLUMNS = (
    "schema_version", "producer", "producer_commit", "producer_dirty",
    "point_id", "construction_ok", "rejection_stage", "rejection_reason",
    "theory_ok_v1", "g_hH2H2_GeV", "total_width_GeV", "ctau_mm",
)

MALFORMED = "MALFORMED"
EVALUATOR_ERROR = "EVALUATOR_ERROR"
EVALUATED = "EVALUATED"


@dataclass
class Attempt:
    row_index: int
    proposal_id: str
    attempt_status: str
    attempt_stage: str
    attempt_reason: str
    mh_GeV_used: str
    canonical: Dict[str, str] = field(default_factory=dict)


def _parse_float(raw: Optional[str]) -> float:
    if raw is None:
        raise ValueError("missing")
    value = float(raw.strip())
    if value != value or value in (float("inf"), float("-inf")):  # NaN/inf
        raise ValueError("non_finite")
    return value


def _row_violations(row_index: int, raw: Dict[str, str], dupe_ids: frozenset) -> List[str]:
    """Formatting-only contract checks. Never a physics predicate."""
    violations: List[str] = []

    proposal_id = (raw.get("proposal_id") or "").strip()
    if not proposal_id:
        violations.append("missing_proposal_id")
    elif proposal_id in dupe_ids:
        violations.append("duplicate_proposal_id")

    for column in ("mH_GeV", "mA_GeV", "mHp_GeV", "M2_GeV2", "tan_beta", "lambda6", "lambda7", "sin_beta_minus_alpha"):
        try:
            _parse_float(raw.get(column))
        except ValueError:
            violations.append(f"invalid_float:{column}")

    yukawa_raw = raw.get("yukawa_type")
    try:
        yukawa_value = int(str(yukawa_raw).strip())
    except (TypeError, ValueError):
        violations.append("invalid_int:yukawa_type")
    else:
        if yukawa_value not in SUPPORTED_YUKAWA_TYPES:
            violations.append(f"unsupported_yukawa_type:{yukawa_value}")

    mh_raw = (raw.get("mh_GeV") or "").strip() if "mh_GeV" in raw else ""
    if mh_raw:
        try:
            mh_value = _parse_float(mh_raw)
        except ValueError:
            violations.append("invalid_float:mh_GeV")
        else:
            if mh_value != MH_CONVENTION_GEV:
                violations.append(f"mh_convention_mismatch:{mh_value!r}")

    return violations


def read_proposals(path: Path) -> tuple[List[Dict[str, str]], List[str]]:
    """Read proposals.csv. Returns (rows, header_violations).

    header_violations is non-empty when the header itself breaks the
    contract (missing a required column, or an unrecognized extra column);
    in that case every row is unusable regardless of its own content.
    """
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames or []
        rows = list(reader)

    header_set = set(header)
    missing = sorted(set(REQUIRED_COLUMNS) - header_set)
    unknown = sorted(header_set - ALL_INPUT_COLUMNS)
    header_violations: List[str] = []
    if missing:
        header_violations.append(f"missing_required_column:{','.join(missing)}")
    if unknown:
        header_violations.append(f"unknown_column:{','.join(unknown)}")
    return rows, header_violations


def build_command(executable: Path, output: Path, *, campaign_id: str, proposal_id: str,
                   mh: float, mH: float, mA: float, mHp: float, sba: float, tan_beta: float,
                   M2: float, lambda6: float, lambda7: float, yukawa_type: int) -> List[str]:
    return [
        str(executable),
        "--campaign-id", campaign_id, "--run-id", proposal_id,
        "--mh", format_float64(mh),
        "--mH-min", format_float64(mH), "--mH-max", format_float64(mH), "--n-mH", "1",
        "--mA", format_float64(mA), "--mHp", format_float64(mHp),
        "--yukawa-type", str(yukawa_type),
        "--sin-ba", format_float64(sba), "--tan-beta", format_float64(tan_beta),
        "--M2-min", format_float64(M2), "--M2-max", format_float64(M2), "--n-M2", "1",
        "--lambda6", format_float64(lambda6), "--lambda7", format_float64(lambda7),
        "--output", str(output),
    ]


def _evaluate_one(
    *, row_index: int, raw: Dict[str, str], executable: Path, run_dir: Path,
    campaign_id: str, repo_root: Path, git: Dict[str, Optional[str]],
    timeout_s: Optional[float], expected_canonical_columns: Optional[List[str]],
    keep_raw: bool = False,
) -> tuple[Attempt, Optional[List[str]]]:
    """Evaluate one well-formed proposal row.

    Returns (attempt, observed_canonical_header). observed_canonical_header
    is non-None only on a trustworthy EVALUATED outcome, so the caller can
    freeze the batch's canonical column list from the first such row.
    """
    proposal_id = raw["proposal_id"].strip()
    mH = _parse_float(raw["mH_GeV"])
    mA = _parse_float(raw["mA_GeV"])
    mHp = _parse_float(raw["mHp_GeV"])
    M2 = _parse_float(raw["M2_GeV2"])
    tan_beta = _parse_float(raw["tan_beta"])
    lambda6 = _parse_float(raw["lambda6"])
    lambda7 = _parse_float(raw["lambda7"])
    sba = _parse_float(raw["sin_beta_minus_alpha"])
    yukawa_type = int(raw["yukawa_type"].strip())
    mh_raw = (raw.get("mh_GeV") or "").strip() if "mh_GeV" in raw else ""
    mh = _parse_float(mh_raw) if mh_raw else MH_CONVENTION_GEV
    mh_used = format_float64(mh)

    output = run_dir / f"_raw_{row_index:04d}.csv"
    command = build_command(
        executable, output, campaign_id=campaign_id, proposal_id=proposal_id,
        mh=mh, mH=mH, mA=mA, mHp=mHp, sba=sba, tan_beta=tan_beta, M2=M2,
        lambda6=lambda6, lambda7=lambda7, yukawa_type=yukawa_type,
    )
    environment = {
        **dict(os.environ),
        "DIHIGGS_GIT_COMMIT": str(git.get("commit") or "unknown"),
        "DIHIGGS_GIT_DIRTY": str(git.get("is_dirty") or "unknown"),
    }

    def _error(stage: str, reason: str) -> tuple[Attempt, None]:
        return Attempt(row_index, proposal_id, EVALUATOR_ERROR, stage, reason, mh_used), None

    try:
        completed = subprocess.run(
            command, cwd=repo_root, env=environment, capture_output=True, text=True,
            timeout=timeout_s, check=False,
        )
    except subprocess.TimeoutExpired:
        return _error("evaluator_invocation", f"evaluator_timeout:{timeout_s}s")

    if completed.returncode != 0:
        stderr = completed.stderr.strip()
        last_line = stderr.splitlines()[-1] if stderr else f"nonzero_exit_code:{completed.returncode}"
        return _error("evaluator_invocation", f"evaluator_nonzero_exit:{last_line}")

    try:
        if not output.is_file():
            return _error("evaluator_output", "missing_output_file")
        if output.stat().st_size == 0:
            return _error("evaluator_output", "empty_output_file")
        with output.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            observed_header = reader.fieldnames or []
            data_rows = list(reader)
    except (OSError, csv.Error) as error:
        return _error("evaluator_output", f"output_unreadable:{error}")
    finally:
        if not keep_raw:
            try:
                output.unlink(missing_ok=True)
            except OSError:
                pass

    if not observed_header:
        return _error("evaluator_output", "output_header_missing")
    missing_required = [c for c in REQUIRED_CANONICAL_COLUMNS if c not in observed_header]
    if missing_required:
        return _error("evaluator_output", f"output_missing_canonical_columns:{','.join(missing_required)}")
    if expected_canonical_columns is not None and observed_header != expected_canonical_columns:
        return _error("evaluator_output", "canonical_header_mismatch")
    if len(data_rows) != 1:
        return _error("evaluator_output", f"unexpected_row_count:{len(data_rows)}")

    canonical = data_rows[0]
    if not canonical.get("point_id"):
        return _error("evaluator_output", "evaluated_row_missing_point_id")
    if canonical.get("schema_version") != "dihiggs.point.v2":
        return _error("evaluator_output", f"evaluated_row_schema_mismatch:{canonical.get('schema_version')}")

    attempt = Attempt(row_index, proposal_id, EVALUATED, "", "", mh_used, canonical)
    return attempt, observed_header


def _tally(attempts: Sequence[Attempt]) -> Dict[str, int]:
    counts = {"attempted": len(attempts), "malformed": 0, "evaluator_error": 0,
              "evaluated": 0, "successful_construction": 0, "theory_valid": 0, "rejected": 0}
    for attempt in attempts:
        if attempt.attempt_status == MALFORMED:
            counts["malformed"] += 1
        elif attempt.attempt_status == EVALUATOR_ERROR:
            counts["evaluator_error"] += 1
        elif attempt.attempt_status == EVALUATED:
            counts["evaluated"] += 1
            construction_ok = attempt.canonical.get("construction_ok") == "1"
            if construction_ok:
                counts["successful_construction"] += 1
            try:
                theory_ok = float(attempt.canonical.get("theory_ok_v1", "nan")) == 1.0
            except ValueError:
                theory_ok = False
            if theory_ok:
                counts["theory_valid"] += 1
            else:
                counts["rejected"] += 1
    return counts


def write_attempts_csv(path: Path, attempts: Sequence[Attempt], canonical_columns: List[str]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(ATTEMPT_ENVELOPE_COLUMNS) + canonical_columns
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for attempt in attempts:
            row = {
                "row_index": attempt.row_index, "proposal_id": attempt.proposal_id,
                "attempt_status": attempt.attempt_status, "attempt_stage": attempt.attempt_stage,
                "attempt_reason": attempt.attempt_reason, "mh_GeV_used": attempt.mh_GeV_used,
            }
            for column in canonical_columns:
                row[column] = attempt.canonical.get(column, "")
            writer.writerow(row)
    tmp.replace(path)
    return sha256(path)


def validate_attempts_file(path: Path, expected_count: int, canonical_columns: List[str]) -> None:
    """Fail-closed structural check. Raises ValueError on any integrity problem.

    Never trusts subprocess return code == 0 by itself: re-reads and
    re-checks the file this module itself just wrote.
    """
    if not path.is_file():
        raise ValueError("attempts_output_missing")
    if expected_count > 0 and path.stat().st_size == 0:
        raise ValueError("attempts_output_empty")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames or []
        rows = list(reader)
    if not header:
        raise ValueError("attempts_header_missing")
    missing = [c for c in REQUIRED_CANONICAL_COLUMNS if c not in header]
    if missing:
        raise ValueError(f"attempts_missing_canonical_columns:{','.join(missing)}")
    if len(rows) != expected_count:
        raise ValueError(f"attempts_row_count_mismatch:expected={expected_count},got={len(rows)}")
    for row in rows:
        if row.get("attempt_status") == EVALUATED:
            if not row.get("point_id"):
                raise ValueError(f"evaluated_row_missing_point_id:proposal_id={row.get('proposal_id')}")
            if row.get("schema_version") != "dihiggs.point.v2":
                raise ValueError(f"evaluated_row_schema_mismatch:proposal_id={row.get('proposal_id')}")


def _executable_hash(executable: Path) -> Optional[str]:
    return sha256(executable) if executable.is_file() else None


def run_proposal_batch(
    *, executable: Path, proposals_csv: Path, outdir: Path, campaign_id: str, run_name: str,
    repo_root: Path, timeout_s: Optional[float] = None, keep_raw: bool = False,
) -> Dict[str, Any]:
    """Run one proposal batch end to end. Writes attempts.csv and
    batch_manifest.json under outdir/campaign=<campaign_id>/<run_name>/.

    Raises RuntimeError (after writing a status="failed" manifest) if the
    batch cannot be trusted: any structural integrity problem, or any
    MALFORMED / EVALUATOR_ERROR attempt anywhere in the batch. A batch
    containing only well-formed, cleanly-evaluated attempts succeeds even
    when individual points are theory-rejected — that is a normal physics
    outcome, not a failure of this interface.
    """
    run_dir = outdir / f"campaign={campaign_id}" / run_name
    input_csv = run_dir / "proposals.csv"
    attempts_csv = run_dir / "attempts.csv"
    manifest_path = run_dir / "batch_manifest.json"

    run_dir.mkdir(parents=True, exist_ok=True)
    input_bytes = proposals_csv.read_bytes()
    input_csv.write_bytes(input_bytes)
    input_sha = sha256(input_csv)

    git = detect_git_info(repo_root)
    manifest: Dict[str, Any] = {
        "schema_version": "orchestrator.proposal_batch.v1",
        "proposal_contract_version": PROPOSAL_SCHEMA_VERSION,
        "attempt_schema_version": ATTEMPT_SCHEMA_VERSION,
        "created_utc": utc_now_iso(),
        "campaign_id": campaign_id,
        "run_name": run_name,
        "producer": "DihiggsPointV2Evaluator",
        "wrapper": "dihiggs.app.orchestrator.proposal_batch",
        "input_csv": str(input_csv),
        "attempts_csv": str(attempts_csv),
        "input_sha256": input_sha,
        "repository_commit": git.get("commit"),
        "repository_dirty": git.get("is_dirty"),
        "evaluator_executable": str(executable),
        "evaluator_executable_sha256": _executable_hash(executable),
        "evaluator_source_sha256": sha256(repo_root / "dihiggs/src/DihiggsPointV2Evaluator.cpp")
        if (repo_root / "dihiggs/src/DihiggsPointV2Evaluator.cpp").is_file() else None,
        "twohdmc_provenance": _twohdmc_provenance(repo_root),
        "mh_GeV_pin": MH_CONVENTION_GEV,
        "mh_GeV_pin_source": MH_CONVENTION_SOURCE,
        "status": "running",
    }

    rows, header_violations = read_proposals(proposals_csv)
    manifest["requested_row_count"] = len(rows)

    # Duplicate proposal_id detection spans the whole batch, independent of
    # per-row parse success.
    ids = [(r.get("proposal_id") or "").strip() for r in rows]
    dupe_ids = frozenset(pid for pid in set(ids) if pid and ids.count(pid) > 1)

    attempts: List[Attempt] = []
    canonical_columns: Optional[List[str]] = None

    for index, raw in enumerate(rows, start=1):
        proposal_id = (raw.get("proposal_id") or "").strip()
        mh_raw = (raw.get("mh_GeV") or "").strip() if "mh_GeV" in raw else ""
        mh_used_display = mh_raw if mh_raw else format_float64(MH_CONVENTION_GEV)

        if header_violations:
            attempts.append(Attempt(
                index, proposal_id, MALFORMED, "proposal_contract",
                ";".join(header_violations), mh_used_display,
            ))
            continue

        violations = _row_violations(index, raw, dupe_ids)
        if violations:
            attempts.append(Attempt(
                index, proposal_id, MALFORMED, "proposal_contract",
                ";".join(violations), mh_used_display,
            ))
            continue

        attempt, observed_header = _evaluate_one(
            row_index=index, raw=raw, executable=executable, run_dir=run_dir,
            campaign_id=campaign_id, repo_root=repo_root, git=git,
            timeout_s=timeout_s, expected_canonical_columns=canonical_columns,
            keep_raw=keep_raw,
        )
        attempts.append(attempt)
        if observed_header is not None and canonical_columns is None:
            canonical_columns = observed_header

    if canonical_columns is None:
        canonical_columns = list(CANONICAL_HEADER_FALLBACK)

    output_sha = write_attempts_csv(attempts_csv, attempts, canonical_columns)
    manifest["canonical_columns"] = canonical_columns
    manifest["output_sha256"] = output_sha

    counts = _tally(attempts)
    manifest["counts"] = counts

    try:
        validate_attempts_file(attempts_csv, len(rows), canonical_columns)
        if counts["malformed"] or counts["evaluator_error"]:
            raise ValueError(
                f"attempt_contract_violations_present:malformed={counts['malformed']},"
                f"evaluator_error={counts['evaluator_error']}"
            )
    except ValueError as error:
        manifest["status"] = "failed"
        manifest["validation_error"] = str(error)
        safe_write_json(manifest_path, manifest)
        raise RuntimeError(f"proposal batch failed validation: {error}") from error

    manifest["status"] = "complete"
    safe_write_json(manifest_path, manifest)
    return manifest
