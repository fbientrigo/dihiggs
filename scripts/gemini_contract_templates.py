from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from typing import Any


@dataclass
class TemplateMaterialization:
    contract: dict[str, Any]
    real_points: int
    template_name: str


TEMPLATES: dict[str, dict[str, Any]] = {
    "ctau_mphi_mA_refine": {
        "mA_values": [450.0, 475.0, 500.0],
        "lambda6_local_multipliers": [0.97, 1.0, 1.03],
        "tan_beta_local_multipliers": [0.98, 1.0, 1.02],
        "mphi_min": 180.0,
        "mphi_max": 220.0,
        "n_mphi": 40,
        "n_lam1": 1,
    },
    "mphi_lower_edge_probe": {
        "mA_values": [500.0],
        "lambda6_local_multipliers": [0.97, 1.0, 1.03],
        "tan_beta_local_multipliers": [0.96, 0.98, 1.0, 1.02, 1.04],
        "mphi_min": 160.0,
        "mphi_max": 220.0,
        "n_mphi": 80,
        "n_lam1": 1,
    },
    "mA_upper_edge_probe": {
        "mA_values": [450.0, 475.0, 500.0],
        "lambda6_local_multipliers": [1.0],
        "tan_beta_local_multipliers": [0.98, 1.0, 1.02],
        "mphi_min": 180.0,
        "mphi_max": 220.0,
        "n_mphi": 120,
        "n_lam1": 1,
    },
    "control_box": {
        "mA_values": [500.0],
        "lambda6_local_multipliers": [1.0],
        "tan_beta_local_multipliers": [1.0],
        "mphi_min": 195.0,
        "mphi_max": 205.0,
        "n_mphi": 1000,
        "n_lam1": 1,
    },
}

ALLOWED_OVERRIDES = {
    "mphi_min",
    "mphi_max",
    "n_mphi",
    "mA_values",
    "tan_beta_center",
    "lambda6_center",
    "tan_beta_local_multipliers",
    "lambda6_local_multipliers",
}

FORBIDDEN_OVERRIDES = {
    "lambda1",
    "n_lam1",
    "scoring",
    "triple_ok",
    "envelope",
    "command",
    "shell",
    "bash",
    "exec_path",
    "force",
}


def compute_template_points(spec: dict[str, Any]) -> int:
    return int(len(spec["mA_values"]) * len(spec["lambda6_local_multipliers"]) * len(spec["tan_beta_local_multipliers"]) * int(spec["n_mphi"]) * int(spec.get("n_lam1", 1)))


def _inside_range(vals: list[float], lo: float, hi: float) -> bool:
    return all(lo <= float(v) <= hi for v in vals)


def materialize_template(
    *,
    template_name: str,
    overrides: dict[str, Any] | None,
    envelope: dict[str, Any],
    min_points: int,
    max_points: int,
    iteration: int,
    max_total_points: int,
    tan_beta_center_default: float,
    lambda6_center_default: float,
    cycle_id: str | None = None,
    iteration_id: str | None = None,
) -> TemplateMaterialization:
    if template_name not in TEMPLATES:
        raise ValueError(f"unknown template: {template_name}")
    spec = dict(TEMPLATES[template_name])
    ov = dict(overrides or {})

    for k in ov:
        if k in FORBIDDEN_OVERRIDES:
            raise ValueError(f"forbidden override: {k}")
        if k not in ALLOWED_OVERRIDES:
            raise ValueError(f"override not allowed: {k}")

    spec.update({k: ov[k] for k in ov if k in ALLOWED_OVERRIDES})

    spec["mA_values"] = [float(v) for v in spec["mA_values"]]
    spec["lambda6_local_multipliers"] = [float(v) for v in spec["lambda6_local_multipliers"]]
    spec["tan_beta_local_multipliers"] = [float(v) for v in spec["tan_beta_local_multipliers"]]
    spec["mphi_min"] = float(spec["mphi_min"])
    spec["mphi_max"] = float(spec["mphi_max"])
    spec["n_mphi"] = int(spec["n_mphi"])
    spec["n_lam1"] = 1

    if spec["mphi_min"] > spec["mphi_max"]:
        raise ValueError("mphi_min > mphi_max")
    if not _inside_range(spec["mA_values"], float(envelope["mA"][0]), float(envelope["mA"][1])):
        raise ValueError("mA_values outside envelope")
    if not (float(envelope["mphi"][0]) <= spec["mphi_min"] <= float(envelope["mphi"][1])):
        raise ValueError("mphi_min outside envelope")
    if not (float(envelope["mphi"][0]) <= spec["mphi_max"] <= float(envelope["mphi"][1])):
        raise ValueError("mphi_max outside envelope")

    tan_beta_center = float(ov.get("tan_beta_center", tan_beta_center_default))
    lambda6_center = float(ov.get("lambda6_center", lambda6_center_default))
    if not (float(envelope["tan_beta"][0]) <= tan_beta_center <= float(envelope["tan_beta"][1])):
        raise ValueError("tan_beta_center outside envelope")
    if not (float(envelope["lambda6"][0]) <= lambda6_center <= float(envelope["lambda6"][1])):
        raise ValueError("lambda6_center outside envelope")

    points_per_mphi = int(len(spec["mA_values"]) * len(spec["lambda6_local_multipliers"]) * len(spec["tan_beta_local_multipliers"]))
    real_points = int(points_per_mphi * int(spec["n_mphi"]) * 1)
    if real_points < int(min_points):
        needed_n_mphi = max(1, (int(min_points) + points_per_mphi - 1) // points_per_mphi)
        spec["n_mphi"] = int(needed_n_mphi)
        real_points = int(points_per_mphi * int(spec["n_mphi"]) * 1)
    if real_points < int(min_points):
        raise ValueError("template real_points below min_points_per_subrun")
    if real_points > int(max_points):
        raise ValueError("template real_points above max_points_per_subrun")

    identity_payload = {
        "template_name": template_name,
        "mA_values": spec["mA_values"],
        "lambda6_center": lambda6_center,
        "tan_beta_center": tan_beta_center,
        "lambda6_local_multipliers": spec["lambda6_local_multipliers"],
        "tan_beta_local_multipliers": spec["tan_beta_local_multipliers"],
        "mphi_min": spec["mphi_min"],
        "mphi_max": spec["mphi_max"],
        "n_mphi": int(spec["n_mphi"]),
    }
    contract_hash_short = hashlib.sha1(json.dumps(identity_payload, sort_keys=True).encode("utf-8")).hexdigest()[:8]
    run_identity_suffix = f"h{contract_hash_short}"
    if cycle_id:
        run_identity_suffix = f"{run_identity_suffix}_{cycle_id}"
    if iteration_id:
        run_identity_suffix = f"{run_identity_suffix}_{iteration_id}"

    contract = {
        "search_envelope": {
            "sin_ba": [1.0, 1.0],
            "mA": [float(min(spec["mA_values"])), float(max(spec["mA_values"]))],
            "lambda7": [0.0, 0.0],
            "lambda6": [float(envelope["lambda6"][0]), float(envelope["lambda6"][1])],
            "tan_beta": [float(envelope["tan_beta"][0]), float(envelope["tan_beta"][1])],
            "mphi": [spec["mphi_min"], spec["mphi_max"]],
            "lambda1": [1.0, 1.0],
        },
        "budget": {
            "max_iterations": 1,
            "min_raw_points_per_subrun": int(min_points),
            "max_raw_points_per_subrun": int(max_points),
            "max_total_points": int(max_total_points),
            "max_per_subrun_minutes": 45,
            "min_remaining_minutes_to_start": 1,
            "max_failed_subcampaigns": 1,
            "preferred_min_n_lam1": 1,
            "production_target_n_lam1": 1,
            "min_n_mphi_when_varied": int(spec["n_mphi"]),
            "operational_smoke": False,
        },
        "runtime": {
            "exec_path": "/home/fabi/dihiggs/dihiggs/app/PhysScanWithFixings",
            "outdir": "/home/fabi/dihiggs/scripts/out",
            "lake_name": "dihiggs_lake",
            "campaign": f"gemini_worker_iter{iteration:02d}",
        },
        "runtime_thread_policy": {
            "initial_threads": 4,
            "min_threads": 1,
            "backoff_factor": 2,
            "max_retries_per_subcampaign": 2,
            "retry_on_gsl_abort": True,
        },
        "strategy": {
            "name": "template_grid_probe",
            "template_name": template_name,
            "lambda1_fixed": 1.0,
            "n_lam1": 1,
            "sin_ba_fixed": 1.0,
            "lambda7_fixed": 0.0,
            "mA_values": spec["mA_values"],
            "lambda6_center": lambda6_center,
            "lambda6_local_multiplier_values": spec["lambda6_local_multipliers"],
            "tan_beta_center": tan_beta_center,
            "tan_beta_local_multiplier_values": spec["tan_beta_local_multipliers"],
            "mphi_min": spec["mphi_min"],
            "mphi_max": spec["mphi_max"],
            "n_mphi": int(spec["n_mphi"]),
            "expected_raw_points": int(real_points),
            "run_identity_suffix": run_identity_suffix,
        },
        "template_metadata": {
            "template_name": template_name,
            "requested_mA_values": spec["mA_values"],
            "effective_n_mphi": int(spec["n_mphi"]),
            "template_real_points": int(real_points),
            "contract_hash_short": contract_hash_short,
            "cycle_id": cycle_id,
            "iteration_id": iteration_id,
            "run_identity_suffix": run_identity_suffix,
        },
    }

    return TemplateMaterialization(contract=contract, real_points=real_points, template_name=template_name)
