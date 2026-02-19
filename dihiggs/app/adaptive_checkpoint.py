from __future__ import annotations

import csv
import json
import shlex
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TypeAlias, cast

JSONPrimitive: TypeAlias = None | bool | int | float | str
JSONValue: TypeAlias = JSONPrimitive | list["JSONValue"] | dict[str, "JSONValue"]


def _checkpoint_dir(path: str | Path) -> Path:
    checkpoint_dir = Path(path)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    return checkpoint_dir


def write_adaptive_state(checkpoint_dir: str | Path, state: Mapping[str, JSONValue]) -> Path:
    out_path = _checkpoint_dir(checkpoint_dir) / "adaptive_state.json"
    content = json.dumps(dict(state), sort_keys=True) + "\n"
    _ = out_path.write_text(content, encoding="utf-8")
    return out_path


def read_adaptive_state(checkpoint_dir: str | Path) -> dict[str, JSONValue]:
    in_path = Path(checkpoint_dir) / "adaptive_state.json"
    parsed: object = json.loads(in_path.read_text(encoding="utf-8"))  # pyright: ignore[reportAny]
    if not isinstance(parsed, dict):
        raise ValueError("adaptive_state.json must contain a JSON object")
    return cast(dict[str, JSONValue], parsed)


def write_proposals_csv(
    checkpoint_dir: str | Path,
    proposals: Sequence[Mapping[str, object]],
) -> Path:
    out_path = _checkpoint_dir(checkpoint_dir) / "proposals.csv"
    fieldnames = ["proposal_id", "status", "command"]

    with out_path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for proposal in proposals:
            status = str(proposal.get("status", ""))
            if status not in {"PENDING", "DONE"}:
                raise ValueError("proposal status must be PENDING or DONE")
            row = {
                "proposal_id": str(proposal.get("proposal_id", "")),
                "status": status,
                "command": str(proposal.get("command", "")),
            }
            writer.writerow(row)
    return out_path


def _command_to_line(command: str | Sequence[str]) -> str:
    if isinstance(command, str):
        return command
    return shlex.join([str(part) for part in command])


def write_commands_sh(
    checkpoint_dir: str | Path,
    commands: Sequence[str | Sequence[str]],
) -> Path:
    out_path = _checkpoint_dir(checkpoint_dir) / "commands.sh"
    lines = ["#!/usr/bin/env bash", "set -euo pipefail"]
    lines.extend(_command_to_line(command) for command in commands)
    _ = out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path


def load_replay_commands(checkpoint_dir: str | Path) -> list[str]:
    commands_path = Path(checkpoint_dir) / "commands.sh"
    commands: list[str] = []

    for raw_line in commands_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("#!"):
            continue
        if line.startswith("#"):
            continue
        if line.startswith("set "):
            continue
        commands.append(line)

    return commands
