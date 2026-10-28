from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


# Vercel documents this TypeSafe-compatible base URL for existing System One clients.
DEFAULT_BASE_URL = "https://ai-gateway.vercel.sh/typesafe"
DEFAULT_PROXY_URL = "http://127.0.0.1:8765"

TEMPLATE_DESCRIPTIONS: dict[str, str] = {
    "ctau_mphi_mA_refine": (
        "Exploit locally around the current tan(beta)/lambda6 anchor while scanning "
        "m_phi and several m_A values. Prefer when recent runs show useful valid points "
        "and the current best region still appears improvable."
    ),
    "mphi_lower_edge_probe": (
        "Probe the lower m_phi edge with a wider local tan(beta) bracket. Prefer when "
        "good points accumulate near the low-m_phi boundary or coverage there is weak."
    ),
    "mA_upper_edge_probe": (
        "Probe the upper m_A edge while keeping lambda6 narrow. Prefer when evidence "
        "suggests improvement toward high m_A or that edge remains underexplored."
    ),
    "sensitivity_probe": (
        "Run a controlled local sensitivity study around the anchor. Prefer when the "
        "next useful step is to learn which local parameter direction matters rather "
        "than immediately maximize the objective."
    ),
    "control_box": (
        "Run a dense, narrow control region. Prefer when a stable reference/control "
        "measurement is more valuable than further broad exploration."
    ),
    "stop": (
        "Stop this autonomous cycle. Prefer when the point budget is exhausted, the "
        "target has been reached, repeated runs show saturation, or evidence is too "
        "weak to justify another bounded experiment."
    ),
}


class JevError(RuntimeError):
    pass


@dataclass(frozen=True)
class JevDecision:
    worker_action: dict[str, Any]
    request_payload: dict[str, Any]
    raw_response: dict[str, Any]

    @property
    def selected(self) -> str:
        return str(self.worker_action.get("template_name") or self.worker_action.get("action"))


def build_scientific_state(config: Mapping[str, Any], worker_state: Mapping[str, Any]) -> dict[str, Any]:
    """Build the compact, auditable state for one bounded DiHiggs policy decision."""
    executed = int(worker_state.get("executed_points_total", 0) or 0)
    max_total = int(config["max_points_total"])

    # Explicit unit resolution: canonical policy unit is mm. When *_mm is present it is
    # always authoritative; *_m is always derived from it (mm/1000), never read as stored,
    # so a stale/mismatched *_m value can never silently override the mm value (factor-1000 guard).
    if "best_ctau_mm" in worker_state:
        best_ctau_mm = float(worker_state["best_ctau_mm"])
        best_ctau_m = best_ctau_mm / 1000.0
    elif "best_ctau_m" in worker_state:
        best_ctau_m = float(worker_state["best_ctau_m"])
        best_ctau_mm = best_ctau_m * 1000.0
    elif "baseline_best_ctau_mm" in config:
        best_ctau_mm = float(config["baseline_best_ctau_mm"])
        best_ctau_m = best_ctau_mm / 1000.0
    else:
        best_ctau_m = float(config["baseline_best_ctau_m"])
        best_ctau_mm = best_ctau_m * 1000.0

    if "target_ctau_mm" in config:
        target_ctau_mm = float(config["target_ctau_mm"])
        target_ctau_m = target_ctau_mm / 1000.0
    elif "target_ctau_m" in config:
        target_ctau_m = float(config["target_ctau_m"])
        target_ctau_mm = target_ctau_m * 1000.0
    else:
        target_ctau_mm = 0.0
        target_ctau_m = 0.0

    history = worker_state.get("history", [])
    if not isinstance(history, list):
        history = []

    last_summary = worker_state.get("last_iteration_summary_for_prompt")
    if not isinstance(last_summary, dict):
        last_summary = None

    return {
        "domain": "bounded DiHiggs / 2HDM parameter-space search",
        "objective": {
            "primary": "increase best valid ctau_mm within the configured bounded search",
            "unit": "mm",
            "best_ctau_mm": best_ctau_mm,
            "target_ctau_mm": target_ctau_mm,
            "best_ctau_m": best_ctau_m,
            "target_ctau_m": target_ctau_m,
            "target_fraction_reached": (best_ctau_mm / target_ctau_mm) if target_ctau_mm > 0 else None,
        },
        "budget": {
            "executed_points": executed,
            "max_points_total": max_total,
            "remaining_points": max(0, max_total - executed),
            "iterations_completed": int(worker_state.get("iterations_completed", 0) or 0),
            "max_iterations": int(config["max_iterations"]),
        },
        "hard_constraints": {
            "search_envelope": config["envelope"],
            "min_points_per_subrun": int(config["min_points_per_subrun"]),
            "max_points_per_subrun": int(config["max_points_per_subrun"]),
            "lambda1_policy": "fixed at 1.0 for this worker unless a separately audited mode enables variation",
            "contract_policy": "only allowlisted templates; materialization and plancheck must pass before execution",
            "forbidden_actions": list(config.get("forbidden_actions", [])),
        },
        "search_context": {
            "last_uncontaminated_iteration": last_summary,
            "recent_history": history[-12:],
            "no_improvement_count": int(worker_state.get("no_improvement_count", 0) or 0),
            "rejected_points_total": int(worker_state.get("rejected_points_total", 0) or 0),
            "contaminated_summary_count": int(worker_state.get("contaminated_summary_count", 0) or 0),
        },
        "available_operators": TEMPLATE_DESCRIPTIONS,
        "decision_rule": (
            "Choose the next bounded experiment that is most informative/useful per point, balancing "
            "objective improvement, valid-point yield, information gain, and avoiding repeated unproductive probes."
        ),
    }


def build_jev_request(
    config: Mapping[str, Any],
    worker_state: Mapping[str, Any],
    *,
    model: str | None = None,
) -> dict[str, Any]:
    """Build a TypeSafe System One request.

    Vercel's TypeSafe-compatible endpoint preserves the System One request shape and
    does not require the native AI Gateway model slug in the payload. ``model`` is
    accepted only for compatibility/testing and omitted when unset.
    """
    state = build_scientific_state(config, worker_state)
    payload: dict[str, Any] = {
        "state": state,
        "questions": {
            "next_operator": {
                "type": "choice",
                "instructions": (
                    "Choose exactly one next operator for the bounded scientific search. Use only the evidence "
                    "in state. Prefer expected scientific value per evaluated point, not novelty for its own sake."
                ),
                "criteria": dict(TEMPLATE_DESCRIPTIONS),
            },
            "local_saturation": {
                "type": "score",
                "instructions": "Rate evidence that the currently exploited local region is saturated.",
                "criteria": [
                    "0: no evidence of saturation; local exploitation is still strongly justified",
                    "1: weak evidence of saturation",
                    "2: mixed/ambiguous evidence",
                    "3: substantial evidence of saturation",
                    "4: strong evidence; another similar local probe is unlikely to help",
                ],
            },
            "exploitation_value": {
                "type": "score",
                "instructions": "Rate the expected scientific value of another local exploit/refine run.",
                "criteria": [
                    "0: very low value", "1: low value", "2: uncertain/moderate value",
                    "3: high value", "4: very high value",
                ],
            },
            "diversification_value": {
                "type": "score",
                "instructions": "Rate the expected scientific value of probing a different edge/direction.",
                "criteria": [
                    "0: very low value", "1: low value", "2: uncertain/moderate value",
                    "3: high value", "4: very high value",
                ],
            },
            "frontier_escalation": {
                "type": "noul",
                "instructions": (
                    "Should this state be escalated to a more capable generative/reasoning model because the "
                    "allowlisted operators are insufficient or the evidence is unusually ambiguous?"
                ),
                "criteria": {
                    "true": "A new strategy or deeper scientific reasoning is likely needed.",
                    "false": "One of the existing bounded operators is sufficient for the next step.",
                },
            },
        },
    }
    if model:
        payload["model"] = model
    return payload


def _request_json(
    url: str,
    payload: Mapping[str, Any],
    *,
    api_key: str | None,
    timeout_s: float,
) -> dict[str, Any]:
    data = json.dumps(payload, separators=(",", ":")).encode("utf-8")
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json",
        "User-Agent": "dihiggs-jev-policy/0.2",
    }
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"
    req = urllib.request.Request(url, data=data, method="POST", headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            body = resp.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise JevError(f"Jev HTTP {exc.code}: {body[:1000]}") from exc
    except urllib.error.URLError as exc:
        raise JevError(f"Jev connection failed: {exc.reason}") from exc

    try:
        obj = json.loads(body)
    except json.JSONDecodeError as exc:
        raise JevError(f"Jev returned non-JSON: {body[:1000]}") from exc
    if not isinstance(obj, dict):
        raise JevError("Jev response must be a JSON object")
    return obj


def call_jev(
    payload: Mapping[str, Any],
    *,
    api_key: str | None = None,
    base_url: str | None = None,
    proxy_url: str | None = None,
    timeout_s: float = 30.0,
    retries: int = 1,
) -> dict[str, Any]:
    """Call Jev through either the local credential-isolating proxy or Vercel directly.

    Secure mode: set JEV_PROXY_URL=http://127.0.0.1:8765. The client then has no
    credential at all. Direct mode uses Vercel's standard AI_GATEWAY_API_KEY.
    """
    selected_proxy = proxy_url if proxy_url is not None else os.environ.get("JEV_PROXY_URL")
    if selected_proxy:
        url = f"{selected_proxy.rstrip('/')}/v1/systemone"
        key = None
    else:
        key = api_key or os.environ.get("AI_GATEWAY_API_KEY")
        if not key:
            raise JevError(
                "No credential path configured. Prefer JEV_PROXY_URL for secure mode; "
                "otherwise set Vercel's standard AI_GATEWAY_API_KEY."
            )
        selected_base = base_url or os.environ.get("JEV_BASE_URL", DEFAULT_BASE_URL)
        url = f"{selected_base.rstrip('/')}/v1/systemone"

    last: Exception | None = None
    for attempt in range(retries + 1):
        try:
            return _request_json(url, payload, api_key=key, timeout_s=timeout_s)
        except JevError as exc:
            last = exc
            if attempt >= retries:
                break
            time.sleep(0.5 * (2**attempt))
    assert last is not None
    raise last


def response_to_worker_action(response: Mapping[str, Any]) -> dict[str, Any]:
    answers = response.get("answers")
    if not isinstance(answers, dict):
        raise JevError("Jev response missing answers")

    next_answer = answers.get("next_operator")
    if not isinstance(next_answer, dict) or next_answer.get("type") != "choice":
        raise JevError("Jev response missing choice answer 'next_operator'")

    selected = str(next_answer.get("choice", ""))
    if selected not in TEMPLATE_DESCRIPTIONS:
        raise JevError(f"Jev selected non-allowlisted operator: {selected!r}")

    probabilities = next_answer.get("probabilities", {})
    confidence = next_answer.get("confidence")
    rationale = (
        f"Jev structured choice={selected}; confidence={confidence}; "
        f"probabilities={json.dumps(probabilities, sort_keys=True)}"
    )

    if selected == "stop":
        return {
            "action": "stop",
            "rationale": rationale,
            "safety_checks": ["jev_choice_from_allowlist"],
            "stop_conditions": ["jev_selected_stop"],
        }

    return {
        "action": "run_template",
        "template_name": selected,
        "template_overrides": {},
        "rationale": rationale,
        "safety_checks": [
            "jev_choice_from_allowlist",
            "template_materialization_required",
            "plancheck_required_before_execution",
        ],
        "stop_conditions": ["budget_exhausted", "target_reached", "future_policy_stop"],
    }


def _extract_json_object(text: str) -> dict[str, Any]:
    """Parse a JSON object out of CLI output, tolerating markdown fences/prose."""
    stripped = text.strip()
    try:
        obj = json.loads(stripped)
        if isinstance(obj, dict):
            return obj
    except json.JSONDecodeError:
        pass
    match = re.search(r"\{.*\}", stripped, re.DOTALL)
    if match:
        obj = json.loads(match.group(0))
        if isinstance(obj, dict):
            return obj
    raise JevError(f"gemini CLI returned non-JSON-object output: {stripped[:1000]}")


GEMINI_CLI_BIN = "agy"  # Antigravity CLI; the old `gemini` CLI's free tier was retired.


def choose_with_gemini(
    config: Mapping[str, Any],
    worker_state: Mapping[str, Any],
    *,
    timeout_s: float = 180.0,
) -> JevDecision:
    """Answer the same System One payload Jev would get, via the local Antigravity (`agy`) CLI.

    Uses response_to_worker_action for parsing, so the allowlist enforcement is
    identical to the Jev path (a non-allowlisted choice raises JevError). This makes
    no Jev/Vercel request; it is a free `agy -p` call.
    """
    payload = build_jev_request(config, worker_state)
    prompt = (
        "Answer these questions about the given scientific state. Output ONLY a JSON "
        "object of this exact shape, no markdown, no prose:\n"
        '{"answers": {"next_operator": {"type": "choice", "choice": "<one allowed key>", '
        '"confidence": <0-1 float>, "probabilities": {"<key>": <0-1 float>, ...}}, '
        '"local_saturation": {"type": "score", "value": <int>}, '
        '"exploitation_value": {"type": "score", "value": <int>}, '
        '"diversification_value": {"type": "score", "value": <int>}, '
        '"frontier_escalation": {"type": "noul", "value": <true|false>}}}\n\n'
        f"Request:\n{json.dumps(payload, sort_keys=True)}\n"
    )
    proc = subprocess.run(
        [GEMINI_CLI_BIN, "-p", prompt], text=True, capture_output=True, timeout=timeout_s
    )
    if proc.returncode != 0:
        raise JevError(f"gemini CLI exited {proc.returncode}: {proc.stderr[-1000:]}")
    raw_response = _extract_json_object(proc.stdout)
    raw_response.setdefault("model", "gemini-cli")
    action = response_to_worker_action(raw_response)
    return JevDecision(worker_action=action, request_payload=payload, raw_response=raw_response)


def choose_with_jev(
    config: Mapping[str, Any],
    worker_state: Mapping[str, Any],
    *,
    model: str | None = None,
    api_key: str | None = None,
    base_url: str | None = None,
    proxy_url: str | None = None,
) -> JevDecision:
    payload = build_jev_request(config, worker_state, model=model)
    response = call_jev(payload, api_key=api_key, base_url=base_url, proxy_url=proxy_url)
    action = response_to_worker_action(response)
    return JevDecision(worker_action=action, request_payload=payload, raw_response=response)


def build_decision_event(
    decision: JevDecision,
    *,
    timestamp: str,
    iteration: int,
    iteration_id: str,
    cycle_id: str,
) -> dict[str, Any]:
    """Training-oriented event: state -> policy probabilities -> chosen action."""
    return {
        "schema_version": "dihiggs_policy_decision_v1",
        "timestamp": timestamp,
        "iteration": iteration,
        "iteration_id": iteration_id,
        "cycle_id": cycle_id,
        "backend": "jev_via_vercel_ai_gateway",
        "state": decision.request_payload["state"],
        "questions": decision.request_payload["questions"],
        "selected_action": decision.worker_action,
        "model": decision.raw_response.get("model", "typesafe-ai/jev"),
        "answers": decision.raw_response.get("answers"),
        "usage": decision.raw_response.get("usage"),
    }
