"""The finite, registered action menu over dihiggs_boundary strategies.

An action is (strategy, tan_beta slice, regime). The lake is the memory: before a launch `lake explored --export`
for the action's exact regime feeds boundary's --explored, and afterwards the output is `lake commit`ted. autoresearch only picks a key from this menu (RPE ADR-020): the
regime (lambda7, mh, mass convention) is the only physics the action sets; every other number comes from the
boundary config, and a runtime override may never touch `physics`.
"""
from __future__ import annotations

import copy
import hashlib
import json
import subprocess
from pathlib import Path

import yaml

# strategy -> boundary module, its config (relative to the boundary checkout), whether it takes --tan-beta, and the
# lake table its points go to (None: it cannot read --explored yet, so launching it would repeat lake points)
STRATEGIES = {
    "adaptive_m2_mesh": {"module": "dhb.adaptive_m2_mesh", "config": "configs/adaptive_m2_boundary_v1.yaml", "slices": False, "table": None},
    "profiled_m2": {"module": "dhb.profiled_m2_scan", "config": "configs/profiled_m2_boundary_v1.yaml", "slices": True, "table": None},
    "survivor_v2": {"module": "dhb.survivor_adaptive_campaign", "config": "campaigns/survivor_adaptive_12520_tanbeta_v2/campaign.yaml",
                    "slices": True, "table": "survivor_points=raw_fresh_12520_tanbeta_v2_points.csv@compact"},
    # MH-MCMC / PSO: fresh continuous points every attempt (seeded per attempt), so a repeat never returns 0 rows
    **{sampler: {"module": "dhb.stochastic_search", "config": "campaigns/stochastic_12520_v1/campaign.yaml", "slices": True,
                 "table": "stochastic_points=raw_stochastic_points.csv@compact", "sampler": sampler} for sampler in ("mh", "pso")},
}
# lambda7 = 0 is superseded: the default menu offers lambda7 != 0 only
REGIMES = {
    "l7_0p01_mh125p20": {"lambda7": 0.01, "mh_GeV": 125.2, "mass_convention_id": "mh_125p20_active_2026"},
}


def _config(boundary_root: Path, strategy: str) -> tuple[Path, dict]:
    path = Path(boundary_root) / STRATEGIES[strategy]["config"]
    return path, yaml.safe_load(path.read_text(encoding="utf-8"))


def build_menu(boundary_root: str | Path) -> list[dict]:
    """Every (strategy, slice, regime); slices are read from boundary's own configs, so the digest covers them."""
    menu = []
    for strategy, spec in sorted(STRATEGIES.items()):
        slices = _config(Path(boundary_root), strategy)[1]["tan_beta_slices"] if spec["slices"] else [None]
        for tb in slices:
            for regime_id, regime in sorted(REGIMES.items()):
                key = f"{strategy}|tb={'all' if tb is None else repr(float(tb))}|{regime_id}"
                menu.append({"key": key, "strategy": strategy, "tan_beta": tb, "regime": regime})
    return menu


def digest(menu: list[dict]) -> str:
    return hashlib.sha256(json.dumps(menu, sort_keys=True).encode()).hexdigest()


def register(path: str | Path, menu: list[dict]) -> str:
    """Write the menu and its digest once; a campaign resumed with another menu is refused."""
    path, sha = Path(path), digest(menu)
    if path.exists():
        old = json.loads(path.read_text(encoding="utf-8"))["digest"]
        if old != sha:
            raise ValueError(f"action menu changed under this campaign: registered {old[:12]}, now {sha[:12]}")
    else:
        path.write_text(json.dumps({"digest": sha, "actions": menu}, indent=1), encoding="utf-8")
    return sha


def lookup(menu: list[dict], key: str) -> dict:
    for action in menu:
        if action["key"] == key:
            return action
    raise ValueError(f"action {key!r} is not in the registered menu")


def boundary_config(action: dict, boundary_root: str | Path, overrides: dict | None = None) -> dict:
    """The boundary config for this action: its regime on top of the strategy's config, plus non-physics overrides."""
    path, cfg = _config(Path(boundary_root), action["strategy"])
    cfg["physics"].update(action["regime"])
    # the copy lives outside boundary: its config-relative paths become absolute
    cfg["dihiggs_root"] = str((path.parent / cfg["dihiggs_root"]).resolve())
    for section, keys in (("seed", ("raw_csv", "manifest")), ("historical_hint", ("manifest", "state", "archive_manifest"))):
        for key in keys:
            if key in cfg.get(section, {}):
                cfg[section][key] = str((path.parent / cfg[section][key]).resolve())
    for section, values in (overrides or {}).items():
        if section == "physics":
            raise ValueError("runtime overrides may not set physics; the action's regime owns it")
        if isinstance(values, dict):
            cfg.setdefault(section, {}).update(copy.deepcopy(values))
        else:
            cfg[section] = copy.deepcopy(values)
    return cfg


def attempt_seed(tag: str) -> int:
    """The sampler seed of an attempt (`campaign/attempt`): a retry of the same attempt repeats it, the next one differs."""
    return int(hashlib.sha256(tag.encode()).hexdigest()[:8], 16)


def command(action: dict, config_path: str | Path, output_dir: str | Path, python: str, seed: int = 0) -> list[str]:
    spec = STRATEGIES[action["strategy"]]
    cmd = [python, "-m", spec["module"], "--config", str(config_path), "--output-dir", str(output_dir)]
    cmd += ["--tan-beta", repr(float(action["tan_beta"]))] if action["tan_beta"] is not None else []
    return cmd + (["--sampler", spec["sampler"], "--seed", str(seed)] if "sampler" in spec else [])


def explored_export(lake: list[str], regime: dict, path: str | Path) -> bool:
    """`lake explored --export` for this regime (boundary re-checks it exactly); False only for a lake with no points yet."""
    result = subprocess.run([*lake, "explored", "--lambda7", repr(regime["lambda7"]), "--mh", repr(regime["mh_GeV"]), "--export", str(path)],
                            capture_output=True, text=True, check=False)
    if result.returncode == 0:
        return True
    if "no point tables" in result.stderr or '"commits/*.json"' in result.stderr:  # no table / no commit at all
        return False
    raise RuntimeError(f"lake explored failed, not launching: {result.stderr.strip()}")


def lake_commit(lake: list[str], run_dir: str | Path, action: dict, tag: str, on_start=lambda pid: None) -> str:
    """`tag` (campaign/attempt) goes in the message so a restart can find a commit whose journal row was lost;
    `on_start(pid)` lets the caller journal the running commit, which outlives a killed daemon."""
    proc = subprocess.Popen([*lake, "commit", str(run_dir), "-m", f"autoresearch {tag} {action['key']}", "--authority", "candidate",
                             "--table", STRATEGIES[action["strategy"]]["table"], "--producer-manifest", str(Path(run_dir) / "manifest.json")],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    on_start(proc.pid)
    out, err = proc.communicate()
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.returncode, proc.args, out, err)
    return out.split()[1].rstrip(":")  # "committed <cid>: ..."


def cleanup(lake: list[str], cid: str | None, work: str | Path) -> int:
    """The lake as pipeline: once `lake verify` confirms commit `cid` remotely, delete the action's outputs under work/out
    except its manifest.json files (logs, boundary config and explored.csv are outside out/ and stay). Returns bytes freed;
    0 and nothing deleted without a verified commit."""
    if not cid or subprocess.run([*lake, "verify", cid], capture_output=True, text=True, check=False).stdout.strip() != "ok":
        return 0
    freed, out = 0, Path(work) / "out"
    for path in sorted(out.rglob("*"), reverse=True):  # children before their directory
        if path.is_file() and path.name != "manifest.json":
            freed += path.stat().st_size
            path.unlink()
        elif path.is_dir() and not any(path.iterdir()):
            path.rmdir()
    return freed


def find_commit(lake: list[str], tag: str) -> str | None:
    """The cid of a commit made for this attempt (`lake ls` message field), or None."""
    out = subprocess.run([*lake, "ls"], capture_output=True, text=True, check=True).stdout
    return next((line.split("\t")[0] for line in out.splitlines() if f"autoresearch {tag} " in line), None)
