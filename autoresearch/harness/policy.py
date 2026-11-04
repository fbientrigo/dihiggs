"""Which arm next: selector (UCB1) -> abstain -> agy delegate -> rules. Smallest version of RPE ADR-018.

UCB1 is the default policy. It abstains only with `abstain:no_support`: every runnable action arm has run and its latest
attempt brought zero fresh rows, so its rewards say nothing. Only then is the local `agy` CLI asked, and it may only
answer one offered menu key. Invalid output, a timeout, a failed call or a spent call budget fall through to rules
(the least-pulled arm), which always answer. Every hop is recorded with its reason. The law is EXACT for UCB1 and rules
(deterministic: the chosen arm has probability 1) and UNKNOWN for agy. Arms whose strategy cannot read --explored are
masked (unavailable, not failed, never retried). An arm whose last 2 closed attempts had 0 fresh rows is masked as
saturated; when every runnable arm is, the decision is `stop:saturated` (no agy call, nothing launched) and the campaign
ends SATURATED.
"""
from __future__ import annotations

import json
import subprocess
from collections.abc import Mapping

from autoresearch.harness import boundary_actions
from autoresearch.harness.dihiggs_adaptation import select_arm_ucb1

PROMPT = """You choose the next exploration action for a physics parameter scan. Every action below has already run and
its latest run found no new points. Answer with exactly one action key from this list and nothing else.
{actions}"""


class Saturated(Exception):
    """Every runnable arm's last 2 closed attempts had 0 fresh rows: a clean terminal stop."""


def _runnable(arm: Mapping) -> bool:
    return "action" not in arm or boundary_actions.STRATEGIES[str(arm["action"]).split("|")[0]]["table"] is not None


def choose(config: Mapping, bandit_state, closed: list[dict], agy_calls: int) -> dict:
    """The decision payload: arm_id, action, law and the hops that led to it. `closed` = journal closed payloads."""
    arms = [a for a in config.get("arms", []) if isinstance(a, Mapping) and a.get("id")]
    hops = [{"hop": "masked", "arm_id": a["id"], "reason": "unavailable:cannot_read_explored"} for a in arms if not _runnable(a)]
    offered = [a for a in arms if _runnable(a)]
    if not offered:
        raise ValueError("no runnable arm: every arm's strategy cannot read --explored")

    fresh = {}
    for c in closed:
        if c.get("status") != "not_executed":
            fresh.setdefault(c["arm_id"], []).append(int(c.get("fresh_rows", 0)))
    saturated = [a for a in offered if len(fresh.get(a["id"], [])) >= 2 and not any(fresh[a["id"]][-2:])]
    hops += [{"hop": "masked", "arm_id": a["id"], "reason": "saturated:last_2_attempts_0_fresh_rows"} for a in saturated]
    offered = [a for a in offered if a not in saturated]
    if not offered:
        hops.append({"hop": "stop:saturated", "reason": "every runnable arm's last 2 closed attempts had 0 fresh rows"})
        return {"arm_id": None, "action": None, "law": "EXACT", "hops": hops, "stop": "saturated"}

    def done(arm, law, reason):
        hops.append({"hop": reason, "arm_id": arm["id"]})
        return {"arm_id": arm["id"], "action": arm.get("action"), "law": law, "hops": hops}

    latest = {c["arm_id"]: c for c in closed if c.get("status") != "not_executed"}
    actions = [a for a in offered if "action" in a]
    if not actions or any(a["id"] not in latest or latest[a["id"]].get("fresh_rows", 0) > 0 for a in actions):
        ucb = select_arm_ucb1(bandit_state, {**config, "arms": offered})
        return done(next(a for a in offered if a["id"] == ucb), "EXACT", "ucb1:selected")
    hops.append({"hop": "abstain:no_support", "reason": "every runnable action's latest attempt had 0 fresh rows"})

    policy = config.get("policy") or {}
    if agy_calls >= int(policy.get("agy_max_calls", 10)):
        hops.append({"hop": "agy:budget_exhausted", "calls": agy_calls})
    else:
        keys = {a["action"]: a for a in actions}
        stats = [{"key": k, **{f: latest[a["id"]].get(f) for f in ("fresh_rows", "seeded_rows", "valid_rows")}} for k, a in keys.items()]
        cmd = [*policy.get("agy_cmd", ["agy", "-p"]), PROMPT.format(actions=json.dumps(stats, indent=1))]
        timeout = float(policy.get("agy_timeout_sec", 120))
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=False)
            answer = out.stdout.strip().strip("`").strip()
            if out.returncode:
                hops.append({"hop": "agy:error", "returncode": out.returncode, "stderr": out.stderr[-500:]})
            elif answer in keys:
                return done(keys[answer], "UNKNOWN", "agy:selected")
            else:
                hops.append({"hop": "agy:invalid_output", "raw": out.stdout[:2000], "stderr": out.stderr[-500:]})
        except subprocess.TimeoutExpired:
            hops.append({"hop": "agy:timeout", "timeout_sec": timeout})
        except OSError as exc:
            hops.append({"hop": "agy:unavailable", "error": str(exc)})
    pulls = {arm_id: int(s.get("pulls", 0)) for arm_id, s in bandit_state.arm_stats.items()}
    return done(min(offered, key=lambda a: pulls.get(a["id"], 0)), "EXACT", "rules:selected")


def agy_calls(decisions: list[dict]) -> int:
    """Calls already spent (any agy hop except budget_exhausted), counted from the journal so a restart keeps them."""
    return sum(h["hop"].startswith("agy:") and h["hop"] != "agy:budget_exhausted" for d in decisions for h in d["hops"])
