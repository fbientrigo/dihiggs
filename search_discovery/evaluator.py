"""Discovery-layer evaluation.

The frozen substrate is NOT modified and NOT reimplemented. Everything that
carries physics is imported directly from search_substrate so there is exactly
one definition of it in the repository:

    HBAR_C_GEV_MM      frozen ctau constant
    _validity          frozen construction/numerical/theory/width/lifetime gates
    repository_identity, evaluator_identity   frozen provenance stamps
    FIXED              frozen model settings (m_h, yukawa type, alignment, l7)

The only thing this module does differently from search_substrate.evaluator is
the *bounds* applied during normalization: the frozen contract's v1 search
envelope is replaced by the three-tier discovery bounds model. Identity hashing,
argv construction, Yukawa initialization, gate evaluation, ctau computation and
provenance are unchanged. tests/test_discovery_layer.py pins that equivalence.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import subprocess
import time
from pathlib import Path
from typing import Any

from search_substrate.contract import FIXED
from search_substrate.evaluator import HBAR_C_GEV_MM, _float_or_none, _validity
from search_substrate.provenance import evaluator_identity, repository_identity

from .bounds import Envelope, BoundsError, check_global, derive

SCHEMA_VERSION = "dihiggs.llp_discovery.v1"
REQUIRED_COORDS = ("mH_GeV", "mA_GeV", "tan_beta", "X", "Q")


class DiscoveryError(ValueError):
    pass


def _write_json(path: Path, document: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def normalize(proposal: dict[str, Any], envelope: Envelope | None = None) -> dict[str, Any]:
    """Validate a discovery proposal and derive frozen-evaluator parameters.

    candidate_id is computed from the *physical* parameters with exactly the
    construction used by search_substrate.contract.normalize_proposal, so a point
    reached by the discovery layer and the same point reached by the frozen
    substrate carry the same identity and deduplicate against one another.
    """
    if not isinstance(proposal, dict):
        raise DiscoveryError("proposal_must_be_object")
    allowed = {"proposal_id", "strategy_id", "worker_id", "generation", "parent_ids",
               "random_seed", "rationale", "coordinates"}
    unknown = sorted(set(proposal) - allowed)
    if unknown:
        raise DiscoveryError("unknown_proposal_fields:" + ",".join(unknown))
    proposal_id = proposal.get("proposal_id")
    if not isinstance(proposal_id, str) or not proposal_id or "," in proposal_id:
        raise DiscoveryError("invalid_proposal_id")
    coords = proposal.get("coordinates")
    if not isinstance(coords, dict):
        raise DiscoveryError("coordinates_must_be_object")
    extra = sorted(set(coords) - set(REQUIRED_COORDS))
    if extra:
        raise DiscoveryError("forbidden_coordinate:" + ",".join(extra))
    missing = [name for name in REQUIRED_COORDS if name not in coords]
    if missing:
        raise DiscoveryError("missing_coordinates:" + ",".join(missing))
    values = {}
    for name in REQUIRED_COORDS:
        raw = coords[name]
        if isinstance(raw, bool) or not isinstance(raw, (int, float)) or not math.isfinite(float(raw)):
            raise DiscoveryError(f"nonfinite_or_wrong_type:{name}")
        values[name] = float(raw)
    try:
        check_global(**values)
    except BoundsError as error:
        raise DiscoveryError(f"global_bounds:{error}") from error
    if envelope is not None:
        for name in REQUIRED_COORDS:
            low, high = getattr(envelope, name)
            if not low <= values[name] <= high:
                raise DiscoveryError(f"outside_active_envelope:{name}:{values[name]}")
    parameters = derive(**values)
    physical = {**parameters, **FIXED, "mHp_GeV": parameters["mA_GeV"]}
    canonical = json.dumps(physical, sort_keys=True, separators=(",", ":"), allow_nan=False)
    digest = hashlib.sha256(canonical.encode()).hexdigest()
    return {
        "schema_version": SCHEMA_VERSION,
        "proposal_id": proposal_id,
        "strategy_id": str(proposal.get("strategy_id", "unspecified")),
        "worker_id": str(proposal.get("worker_id", "unspecified")),
        "generation": int(proposal.get("generation", 0)),
        "parent_ids": list(proposal.get("parent_ids", []) or []),
        "random_seed": proposal.get("random_seed"),
        "rationale": str(proposal.get("rationale", "")),
        "coordinates": values,
        "parameters": parameters,
        "fixed_model_settings": dict(FIXED),
        "derived": {"X": values["X"], "Q": values["Q"],
                    "S": values["mA_GeV"] ** 2 - values["mH_GeV"] ** 2,
                    "mHp_GeV": parameters["mA_GeV"]},
        "candidate_id": "candidate_" + digest[:16],
        "physical_identity": digest,
        "canonical_physical_json": canonical,
        "envelope_id": envelope.envelope_id if envelope is not None else None,
    }


# Raw observables persisted for every terminated evaluation, per campaign brief.
RAW_OBSERVABLE_FIELDS = (
    "br_bb", "br_cc", "br_tt", "br_tautau", "br_gg", "br_gammagamma", "br_Zgamma",
    "br_WW", "br_ZZ", "br_hh", "total_width_GeV", "width_unaccounted_GeV", "ctau_mm",
    "width_bb_GeV", "width_gg_GeV", "width_gammagamma_GeV", "width_Zgamma_GeV",
    "width_tautau_GeV", "width_hh_GeV", "g_hH2H2_GeV",
    "lambda1_reconstructed", "lambda2_reconstructed", "lambda3_reconstructed",
    "lambda4_reconstructed", "lambda5_reconstructed",
    "m12_sq_input_GeV2", "M2_reconstructed_GeV2", "tan_beta_reconstructed",
)
THEORY_FLAG_FIELDS = (
    "construction_ok", "numerical_ok", "positivity_reported_ok", "unitarity_ok",
    "perturbativity_ok", "stability_reported_ok", "theory_ok_v1",
    "rejection_stage", "rejection_reason",
)


class DiscoveryEvaluator:
    """Runs the frozen DihiggsPointV2Evaluator binary. Physics is untouched."""

    def __init__(self, root: Path, executable: Path):
        self.root = Path(root).resolve()
        self.executable = Path(executable).resolve()

    def build_command(self, normalized: dict[str, Any], output_csv: Path) -> list[str]:
        p = normalized["parameters"]
        f = normalized["fixed_model_settings"]
        return [
            str(self.executable), "--campaign-id", "llp_benchmark_search_v1",
            "--run-id", normalized["proposal_id"], "--mh", str(f["m_h_GeV"]),
            "--mH-min", str(p["mH_GeV"]), "--mH-max", str(p["mH_GeV"]), "--n-mH", "1",
            "--mA", str(p["mA_GeV"]), "--mHp", str(normalized["derived"]["mHp_GeV"]),
            "--yukawa-type", str(f["yukawa_type"]),
            "--sin-ba", str(f["sin_beta_minus_alpha"]),
            "--tan-beta", str(p["tan_beta"]), "--M2-min", str(p["M2_GeV2"]),
            "--M2-max", str(p["M2_GeV2"]), "--n-M2", "1", "--lambda6", str(p["lambda6"]),
            "--lambda7", str(f["lambda7"]), "--output", str(output_csv),
        ]

    def evaluate(self, proposal: dict[str, Any], attempt_dir: Path,
                 envelope: Envelope | None = None) -> dict[str, Any]:
        started = time.time()
        attempt_dir.mkdir(parents=True, exist_ok=True)
        try:
            normalized = normalize(proposal, envelope)
        except DiscoveryError as error:
            result = {
                "schema_version": SCHEMA_VERSION, "status": "FAILED",
                "attempt_id": attempt_dir.name, "failure_stage": "normalization",
                "failure_reason": str(error), "source_candidate": proposal,
                "provenance": {**repository_identity(self.root),
                               **evaluator_identity(self.root, self.executable)},
                "runtime": {"started_unix": started, "finished_unix": time.time()},
            }
            _write_json(attempt_dir / "evaluation.json", result)
            return result
        output_csv = attempt_dir / "canonical_point.csv"
        command = self.build_command(normalized, output_csv)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        env["LC_ALL"] = "C"
        try:
            completed = subprocess.run(command, cwd=self.root, env=env,
                                       capture_output=True, text=True, check=False)
            (attempt_dir / "evaluator.stdout").write_text(completed.stdout, encoding="utf-8")
            (attempt_dir / "evaluator.stderr").write_text(completed.stderr, encoding="utf-8")
        except OSError as error:
            return self._failed(normalized, "evaluator_execution", repr(error), started, command, attempt_dir)
        if completed.returncode != 0:
            return self._failed(normalized, "evaluator_execution",
                                completed.stderr.strip() or f"returncode={completed.returncode}",
                                started, command, attempt_dir)
        try:
            with output_csv.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            if len(rows) != 1:
                raise ValueError(f"expected_one_canonical_row:got_{len(rows)}")
            row = rows[0]
            if row.get("schema_version") != "dihiggs.point.v2" or row.get("producer") != "DihiggsPointV2Evaluator":
                raise ValueError("unexpected_canonical_evaluator_schema")
        except (OSError, ValueError, csv.Error) as error:
            return self._failed(normalized, "canonical_output", str(error), started, command, attempt_dir)

        total_width = _float_or_none(row.get("total_width_GeV"))
        ctau = _float_or_none(row.get("ctau_mm"))
        expected_ctau = HBAR_C_GEV_MM / total_width if total_width and total_width > 0 else None
        ctau_ok = (ctau is not None and expected_ctau is not None
                   and abs(ctau - expected_ctau) <= 1e-10 * max(1.0, abs(ctau)))
        valid, gates = _validity(row, ctau_ok)

        observables = {name: _float_or_none(row.get(name)) for name in RAW_OBSERVABLE_FIELDS}
        theory_flags = {name: row.get(name) for name in THEORY_FLAG_FIELDS}
        from .families import provisional_family  # local import avoids cycle
        provisional = provisional_family(observables, valid)

        result = {
            "schema_version": SCHEMA_VERSION, "status": "TERMINATED",
            "attempt_id": attempt_dir.name,
            "candidate_id": normalized["candidate_id"],
            "proposal_id": normalized["proposal_id"],
            "strategy_id": normalized["strategy_id"], "worker_id": normalized["worker_id"],
            "generation": normalized["generation"], "parent_ids": normalized["parent_ids"],
            "envelope_id": normalized["envelope_id"],
            "coordinates": normalized["coordinates"],
            "parameters": normalized["parameters"],
            "derived": normalized["derived"],
            "fixed_model_settings": normalized["fixed_model_settings"],
            "canonical_point_id": row.get("point_id"),
            "gates": gates, "validity_gate": valid,
            "theory_flags": theory_flags,
            "observables": observables,
            "ctau_mm": ctau, "ctau_recomputed_mm": expected_ctau,
            "ctau_consistency_ok": ctau_ok,
            "provisional_family": provisional,
            "command": command,
            "provenance": {**repository_identity(self.root),
                           **evaluator_identity(self.root, self.executable),
                           "source_lineage": {"proposal_id": normalized["proposal_id"],
                                              "parent_ids": normalized["parent_ids"]}},
            "runtime": {"started_unix": started, "finished_unix": time.time(),
                        "elapsed_s": time.time() - started},
            "artifacts": {"canonical_csv": str(output_csv)},
        }
        _write_json(attempt_dir / "evaluation.json", result)
        return result

    def _failed(self, normalized: dict[str, Any], stage: str, reason: str, started: float,
                command: list[str], attempt_dir: Path) -> dict[str, Any]:
        result = {
            "schema_version": SCHEMA_VERSION, "status": "FAILED", "attempt_id": attempt_dir.name,
            "candidate_id": normalized.get("candidate_id"), "proposal_id": normalized.get("proposal_id"),
            "source_candidate": normalized, "failure_stage": stage, "failure_reason": reason,
            "command": command,
            "provenance": {**repository_identity(self.root),
                           **evaluator_identity(self.root, self.executable)},
            "runtime": {"started_unix": started, "finished_unix": time.time()},
        }
        _write_json(attempt_dir / "evaluation.json", result)
        return result
