from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import time
from pathlib import Path
from typing import Any

from .classification import classify_observables
from .contract import ContractError, normalize_proposal
from .provenance import evaluator_identity, repository_identity
from .score import score_components


HBAR_C_GEV_MM = 1.973269804e-13


def _float_or_none(value: str | None) -> float | None:
    if value is None or value.strip().lower() in {"nan", "inf", "+inf", "-inf"}:
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _write_json(path: Path, document: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _validity(row: dict[str, str], ctau_ok: bool) -> tuple[bool, dict[str, bool]]:
    fields = {
        "construction": row.get("construction_ok") == "1",
        "numerical": row.get("numerical_ok") == "1.00000000000000000e+00",
        "positivity": row.get("positivity_reported_ok") == "1.00000000000000000e+00",
        "unitarity": row.get("unitarity_ok") == "1.00000000000000000e+00",
        "perturbativity": row.get("perturbativity_ok") == "1.00000000000000000e+00",
        "theory": row.get("theory_ok_v1") == "1.00000000000000000e+00",
        "width": row.get("width_ok") == "1.00000000000000000e+00" and _float_or_none(row.get("total_width_GeV")) is not None,
        "lifetime": ctau_ok,
    }
    return all(fields.values()), fields


class CanonicalEvaluator:
    def __init__(self, root: Path, executable: Path):
        self.root = root.resolve()
        self.executable = executable.resolve()

    def evaluate(self, proposal: dict[str, Any], attempt_dir: Path) -> dict[str, Any]:
        started = time.time()
        attempt_dir.mkdir(parents=True, exist_ok=True)
        try:
            normalized = normalize_proposal(proposal)
        except ContractError as error:
            result = {
                "schema_version": "dihiggs.llp.evaluation.v1", "status": "FAILED",
                "attempt_id": attempt_dir.name,
                "failure_stage": "normalization", "failure_reason": str(error),
                "source_candidate": proposal,
                "provenance": {**repository_identity(self.root), **evaluator_identity(self.root, self.executable),
                                "source_lineage": {"proposal_id": proposal.get("proposal_id"), "parent_ids": proposal.get("parent_ids", [])}},
                "runtime": {"started_unix": started, "finished_unix": time.time()},
            }
            _write_json(attempt_dir / "evaluation.json", result)
            return result
        output_csv = attempt_dir / "canonical_point.csv"
        stdout_path, stderr_path = attempt_dir / "evaluator.stdout", attempt_dir / "evaluator.stderr"
        p = normalized["parameters"]
        command = [
            str(self.executable), "--campaign-id", "llp_benchmark_search_v1",
            "--run-id", normalized["proposal_id"], "--mh", str(normalized["fixed_model_settings"]["m_h_GeV"]),
            "--mH-min", str(p["mH_GeV"]), "--mH-max", str(p["mH_GeV"]), "--n-mH", "1",
            "--mA", str(p["mA_GeV"]), "--mHp", str(normalized["derived"]["mHp_GeV"]),
            "--yukawa-type", str(normalized["fixed_model_settings"]["yukawa_type"]),
            "--sin-ba", str(normalized["fixed_model_settings"]["sin_beta_minus_alpha"]),
            "--tan-beta", str(p["tan_beta"]), "--M2-min", str(p["M2_GeV2"]),
            "--M2-max", str(p["M2_GeV2"]), "--n-M2", "1", "--lambda6", str(p["lambda6"]),
            "--lambda7", str(normalized["fixed_model_settings"]["lambda7"]), "--output", str(output_csv),
        ]
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        env["LC_ALL"] = "C"
        try:
            completed = subprocess.run(command, cwd=self.root, env=env, capture_output=True, text=True, check=False)
            stdout_path.write_text(completed.stdout, encoding="utf-8")
            stderr_path.write_text(completed.stderr, encoding="utf-8")
        except OSError as error:
            return self._failed(normalized, "evaluator_execution", repr(error), started, command, attempt_dir)
        if completed.returncode != 0:
            return self._failed(normalized, "evaluator_execution", completed.stderr.strip() or f"returncode={completed.returncode}", started, command, attempt_dir)
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
        ctau_ok = ctau is not None and expected_ctau is not None and abs(ctau - expected_ctau) <= 1e-10 * max(1.0, abs(ctau))
        valid, gates = _validity(row, ctau_ok)
        observables = {name: _float_or_none(row.get(name)) for name in (
            "br_bb", "br_cc", "br_tautau", "br_gammagamma", "br_Zgamma", "br_gg", "br_hh", "ctau_mm"
        )}
        classification = classify_observables(observables)
        result = {
            "schema_version": "dihiggs.llp.evaluation.v1", "status": "TERMINATED", "attempt_id": attempt_dir.name,
            "candidate_id": normalized["candidate_id"], "proposal_id": normalized["proposal_id"],
            "source_candidate": normalized, "canonical_normalized_inputs": normalized["parameters"],
            "fixed_model_settings": normalized["fixed_model_settings"], "derived_inputs": normalized["derived"],
            "X": normalized["derived"]["X"], "canonical_point_id": row.get("point_id"),
            "gates": gates, "validity_gate": valid, "ctau_mm": ctau,
            "ctau_recomputed_mm": expected_ctau, "ctau_consistency_ok": ctau_ok,
            "classification": classification,
            "score_components": score_components(valid=valid, classification=classification, ctau_mm=ctau),
            "canonical_evaluator_row": row,
            "command": command,
            "provenance": {**repository_identity(self.root), **evaluator_identity(self.root, self.executable),
                            "source_lineage": {"proposal_id": normalized["proposal_id"], "parent_ids": normalized["parent_ids"]}},
            "runtime": {"started_unix": started, "finished_unix": time.time(), "elapsed_s": time.time() - started},
            "runtime_and_timestamps": {"started_unix": started, "finished_unix": time.time()},
            "artifacts": {"canonical_csv": str(output_csv), "stdout": str(stdout_path), "stderr": str(stderr_path)},
        }
        _write_json(attempt_dir / "evaluation.json", result)
        return result

    def _failed(self, normalized: dict[str, Any], stage: str, reason: str, started: float,
                command: list[str], attempt_dir: Path) -> dict[str, Any]:
        result = {
            "schema_version": "dihiggs.llp.evaluation.v1", "status": "FAILED", "attempt_id": attempt_dir.name,
            "candidate_id": normalized.get("candidate_id"), "proposal_id": normalized.get("proposal_id"),
            "source_candidate": normalized, "failure_stage": stage, "failure_reason": reason,
            "command": command, "provenance": {**repository_identity(self.root), **evaluator_identity(self.root, self.executable),
                                                  "source_lineage": {"proposal_id": normalized.get("proposal_id"), "parent_ids": normalized.get("parent_ids", [])}},
            "runtime": {"started_unix": started, "finished_unix": time.time(), "elapsed_s": time.time() - started},
            "runtime_and_timestamps": {"started_unix": started, "finished_unix": time.time()},
        }
        _write_json(attempt_dir / "evaluation.json", result)
        return result


def family_proposals(anchor: dict[str, Any]) -> list[dict[str, Any]]:
    normalized = normalize_proposal(anchor)
    p = normalized["parameters"]
    if p["mH_GeV"] != 200.0:
        raise ContractError("family_anchor_mass_must_be_200_GeV")
    q = (p["mH_GeV"] ** 2 - p["M2_GeV2"]) * p["tan_beta"] ** 2
    s = p["mA_GeV"] ** 2 - p["mH_GeV"] ** 2
    x = normalized["derived"]["X"]
    output = []
    for mass in (150.0, 200.0, 250.0):
        m_a = math.sqrt(mass ** 2 + s)
        m2 = mass ** 2 - q / p["tan_beta"] ** 2
        output.append({
            "proposal_id": f"{normalized['proposal_id']}_family_mH{int(mass)}",
            "strategy_id": "deterministic_QS_family_validation", "worker_id": normalized["worker_id"],
            "generation": normalized["generation"], "parent_ids": [normalized["proposal_id"]],
            "random_seed": normalized["random_seed"], "rationale": "derived by frozen Q,S continuation",
            "parameters": {"mH_GeV": mass, "mA_GeV": m_a, "M2_GeV2": m2, "tan_beta": p["tan_beta"], "lambda6": x / p["tan_beta"]},
        })
    return output
