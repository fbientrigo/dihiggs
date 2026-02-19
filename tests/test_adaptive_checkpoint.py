# pyright: reportMissingImports=false, reportUnknownVariableType=false
from __future__ import annotations

import json
from pathlib import Path


def test_adaptive_state_roundtrip_and_deterministic_json(tmp_path: Path):
    from dihiggs.app.adaptive_checkpoint import read_adaptive_state, write_adaptive_state

    checkpoint_dir = tmp_path / "ckpt"
    state = {
        "step": 3,
        "seed": 42,
        "pending": ["lam1_bin_000", "lam1_bin_001"],
    }

    _ = write_adaptive_state(checkpoint_dir, state)
    first_bytes = (checkpoint_dir / "adaptive_state.json").read_bytes()
    _ = write_adaptive_state(checkpoint_dir, state)
    second_bytes = (checkpoint_dir / "adaptive_state.json").read_bytes()

    assert first_bytes == second_bytes
    assert first_bytes.endswith(b"\n")
    assert first_bytes.decode("utf-8") == json.dumps(state, sort_keys=True) + "\n"
    assert read_adaptive_state(checkpoint_dir) == state


def test_write_proposals_csv_deterministic_header_and_row_order(tmp_path: Path):
    from dihiggs.app.adaptive_checkpoint import write_proposals_csv

    checkpoint_dir = tmp_path / "ckpt"
    proposals = [
        {"proposal_id": "p-001", "status": "PENDING", "command": "python run.py --id p-001"},
        {"proposal_id": "p-002", "status": "DONE", "command": "python run.py --id p-002"},
    ]

    _ = write_proposals_csv(checkpoint_dir, proposals)
    first = (checkpoint_dir / "proposals.csv").read_text(encoding="utf-8")
    _ = write_proposals_csv(checkpoint_dir, proposals)
    second = (checkpoint_dir / "proposals.csv").read_text(encoding="utf-8")

    assert first == second
    assert first == (
        "proposal_id,status,command\n"
        "p-001,PENDING,python run.py --id p-001\n"
        "p-002,DONE,python run.py --id p-002\n"
    )


def test_write_commands_sh_with_shebang_and_deterministic_content(tmp_path: Path):
    from dihiggs.app.adaptive_checkpoint import write_commands_sh

    checkpoint_dir = tmp_path / "ckpt"
    commands = [
        ["python", "scan.py", "--label", "A B"],
        "python scan.py --label C",
    ]

    _ = write_commands_sh(checkpoint_dir, commands)
    content = (checkpoint_dir / "commands.sh").read_text(encoding="utf-8")

    assert content == (
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "python scan.py --label 'A B'\n"
        "python scan.py --label C\n"
    )


def test_checkpoint_replay_commands(tmp_path: Path):
    from dihiggs.app.adaptive_checkpoint import load_replay_commands

    checkpoint_dir = tmp_path / "ckpt"
    checkpoint_dir.mkdir(parents=True)
    script = "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            "",
            "# comment",
            "python scan.py --id 1",
            "python scan.py --id 2 --flag",
            "",
        ]
    )
    _ = (checkpoint_dir / "commands.sh").write_text(script, encoding="utf-8")

    assert load_replay_commands(checkpoint_dir) == [
        "python scan.py --id 1",
        "python scan.py --id 2 --flag",
    ]
