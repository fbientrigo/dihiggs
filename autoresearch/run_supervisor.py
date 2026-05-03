#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .harness.campaign_supervisor import CampaignSupervisor
from .harness.dihiggs_preflight import run_all_preflight_checks
from .harness.status_dashboard import print_cli_status


def _load_config(config_path: Path) -> dict[str, object]:
    config_raw = json.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(config_raw, dict):
        raise ValueError(f"Config root must be a JSON object: {config_path}")
    config: dict[str, object] = dict(config_raw)
    runtime = config.get("runtime")
    if isinstance(runtime, dict):
        runtime["python_exe"] = sys.executable
    return config


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the autoresearch campaign supervisor")
    parser.add_argument("config", help="Path to the JSON config file")
    parser.add_argument("--preflight-only", action="store_true", help="Run only preflight checks and exit")
    parser.add_argument("--status", action="store_true", help="Print compact status after the run")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config = _load_config(config_path)
    paths = config.get("paths")
    outdir_value = paths.get("outdir") if isinstance(paths, dict) else None
    outdir = Path(str(outdir_value)).resolve() if outdir_value else config_path.parent.resolve()

    if args.preflight_only:
        print(json.dumps(run_all_preflight_checks(config), indent=2))
        return 0

    supervisor = CampaignSupervisor(config, outdir=outdir)
    terminal_state = supervisor.run()
    print(f"terminal_state={terminal_state}")
    print(f"outdir={outdir}")
    print(f"status_json={outdir / 'campaign_status.json'}")
    print(f"status_html={outdir / 'campaign_status.html'}")
    print(f"status_html={outdir} / 'campaign_status.html'")

    if args.status:
        print()
        print_cli_status(outdir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

import argparse
import json
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from autoresearch.harness.campaign_supervisor import CampaignSupervisor
from autoresearch.harness.dihiggs_preflight import run_all_preflight_checks
from autoresearch.harness.status_dashboard import print_cli_status


def _load_config(config_path: Path) -> dict[str, object]:
    config_raw = json.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(config_raw, dict):
        raise ValueError(f"Config root must be a JSON object: {config_path}")
    config: dict[str, object] = dict(config_raw)
    runtime = config.get("runtime")
    if isinstance(runtime, dict):
        runtime["python_exe"] = sys.executable
    return config


def _as_mapping(value: object) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the autoresearch campaign supervisor")
    parser.add_argument("config", help="Path to the JSON config file")
    parser.add_argument("--preflight-only", action="store_true", help="Run only preflight checks and exit")
    parser.add_argument("--status", action="store_true", help="Print compact status after the run")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config = _load_config(config_path)
    paths = _as_mapping(config.get("paths"))
    outdir_value = paths.get("outdir")
    outdir = Path(str(outdir_value)).resolve() if outdir_value else config_path.parent.resolve()

    if args.preflight_only:
        print(json.dumps(run_all_preflight_checks(config), indent=2))
        return 0

    supervisor = CampaignSupervisor(config, outdir=outdir)
    terminal_state = supervisor.run()
    print(f"terminal_state={terminal_state}")
    print(f"outdir={outdir}")
    print(f"status_json={outdir / 'campaign_status.json'}")
    print(f"status_html={outdir / 'campaign_status.html'}")

    if args.status:
        print()
        print_cli_status(outdir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
