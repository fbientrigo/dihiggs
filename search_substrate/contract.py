from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "docs/contracts/llp_benchmark_search_v1.json"
SCHEMA_VERSION = "dihiggs.llp_benchmark_search.v1"
FIXED = {
    "m_h_GeV": 125.20,
    "yukawa_type": 1,
    "sin_beta_minus_alpha": 1.0,
    "lambda7": 0.0,
}
REQUIRED = ("mH_GeV", "mA_GeV", "M2_GeV2", "tan_beta", "lambda6")
ALIASES = {
    "mH": "mH_GeV", "mH_gev": "mH_GeV", "mA": "mA_GeV", "mA_gev": "mA_GeV",
    "M2": "M2_GeV2", "M2_gev2": "M2_GeV2", "tb": "tan_beta", "tanbeta": "tan_beta",
    "l6": "lambda6",
}
BOUNDS = {
    "mH_GeV": (150.0, 250.0), "mA_GeV": (450.0, 550.0),
    "M2_GeV2": (20000.0, 63000.0), "tan_beta": (30000.0, 3000000.0),
    "lambda6": (0.0001, 0.05),
}


class ContractError(ValueError):
    pass


def load_contract() -> dict[str, Any]:
    document = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    if document.get("schema_version") != SCHEMA_VERSION:
        raise ContractError("contract_schema_version_mismatch")
    return document


def _number(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ContractError(f"unit_or_type_error:{name}")
    value = float(value)
    if not math.isfinite(value):
        raise ContractError(f"nonfinite:{name}")
    return value


def _canonical_parameters(raw: dict[str, Any]) -> dict[str, float]:
    if not isinstance(raw, dict):
        raise ContractError("parameters_must_be_object")
    out: dict[str, float] = {}
    for key, value in raw.items():
        canonical = ALIASES.get(key, key)
        if canonical not in BOUNDS:
            raise ContractError(f"unknown_or_forbidden_parameter:{key}")
        parsed = _number(value, key)
        if canonical in out and out[canonical] != parsed:
            raise ContractError(f"ambiguous_parameter_alias:{canonical}")
        out[canonical] = parsed
    missing = [name for name in REQUIRED if name not in out]
    if missing:
        raise ContractError("missing_parameters:" + ",".join(missing))
    for name, (low, high) in BOUNDS.items():
        if not low <= out[name] <= high:
            raise ContractError(f"out_of_range:{name}:{out[name]}")
    x = out["tan_beta"] * out["lambda6"]
    if not 200.0 <= x <= 100000.0:
        raise ContractError(f"out_of_range:X:{x}")
    if out["mA_GeV"] <= out["mH_GeV"]:
        raise ContractError("invalid_hierarchy:mA_not_above_mH")
    return out


def normalize_proposal(proposal: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(proposal, dict):
        raise ContractError("proposal_must_be_object")
    allowed = {"proposal_id", "strategy_id", "worker_id", "generation", "parent_ids", "random_seed", "rationale", "parameters"}
    unknown = sorted(set(proposal) - allowed)
    if unknown:
        raise ContractError("unknown_proposal_fields:" + ",".join(unknown))
    proposal_id = proposal.get("proposal_id")
    if not isinstance(proposal_id, str) or not proposal_id or "," in proposal_id:
        raise ContractError("invalid_proposal_id")
    parameters = _canonical_parameters(proposal.get("parameters", {}))
    parent_ids = proposal.get("parent_ids", [])
    if not isinstance(parent_ids, list) or not all(isinstance(item, str) and item for item in parent_ids):
        raise ContractError("invalid_parent_ids")
    generation = proposal.get("generation", 0)
    if isinstance(generation, bool) or not isinstance(generation, int) or generation < 0:
        raise ContractError("invalid_generation")
    seed = proposal.get("random_seed")
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
        raise ContractError("invalid_random_seed")
    normalized = {
        "proposal_id": proposal_id,
        "strategy_id": str(proposal.get("strategy_id", "unspecified")),
        "worker_id": str(proposal.get("worker_id", "unspecified")),
        "generation": generation,
        "parent_ids": list(parent_ids),
        "random_seed": seed,
        "rationale": str(proposal.get("rationale", "")),
        "parameters": parameters,
        "fixed_model_settings": dict(FIXED),
        "derived": {"X": parameters["tan_beta"] * parameters["lambda6"], "mHp_GeV": parameters["mA_GeV"]},
        "schema_version": SCHEMA_VERSION,
    }
    physical = {
        **parameters,
        **FIXED,
        "mHp_GeV": parameters["mA_GeV"],
    }
    canonical = json.dumps(physical, sort_keys=True, separators=(",", ":"), allow_nan=False)
    normalized["candidate_id"] = "candidate_" + hashlib.sha256(canonical.encode()).hexdigest()[:16]
    normalized["physical_identity"] = hashlib.sha256(canonical.encode()).hexdigest()
    normalized["canonical_physical_json"] = canonical
    return normalized


def canonical_proposal_json(normalized: dict[str, Any]) -> str:
    return json.dumps(normalized, sort_keys=True, separators=(",", ":"), allow_nan=False)

