#!/usr/bin/env python3
"""Bounded, resumable background harvesting daemon over the frozen LLP
discovery search substrate (search_substrate/ + search_discovery/).

    resume state -> propose -> dedup -> evaluate -> append evidence
    -> update archive -> checkpoint -> repeat

This is a thin CLI surface; all loop logic lives in search_discovery.daemon.
See README_DAEMON.md for operational instructions and IMPLEMENTATION_DECISION.md
for the architecture rationale.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from search_discovery import daemon as daemon_mod  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="benchmark_search_daemon",
        description="Resumable background harvesting daemon over the frozen LLP discovery substrate",
    )
    parser.add_argument("--config", required=True, help="path to a campaign YAML config")
    parser.add_argument("--resume", action="store_true", help="resume from the run-dir's daemon checkpoint")
    parser.add_argument("--dry-run", action="store_true",
                        help="validate config/paths/schemas and print the execution plan; invokes zero physics evaluations")
    parser.add_argument("--workers", type=int, default=None, help="override runtime.workers from the config")
    parser.add_argument("--run-dir", default=None, help="override run_dir from the config")
    parser.add_argument("--max-cycles", type=int, default=None,
                        help="stop after this many cycles regardless of budget (mainly for tests/smoke runs)")
    args = parser.parse_args(argv)

    try:
        config = daemon_mod.load_config(Path(args.config))
    except daemon_mod.DaemonConfigError as error:
        print(f"CONFIG_ERROR: {error}", file=sys.stderr)
        return 2

    run_dir = Path(args.run_dir) if args.run_dir else Path(config["run_dir"])
    if not run_dir.is_absolute():
        run_dir = daemon_mod.ROOT / run_dir

    try:
        daemon_mod.guard_run_dir(run_dir)
    except daemon_mod.DaemonConfigError as error:
        print(f"REFUSED: {error}", file=sys.stderr)
        return 2

    return daemon_mod.run(config, run_dir, resume=args.resume, workers_override=args.workers,
                          dry_run=args.dry_run, stop_after_cycles=args.max_cycles)


if __name__ == "__main__":
    raise SystemExit(main())
