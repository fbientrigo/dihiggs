#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass
from collections.abc import Iterable
from pathlib import Path
from types import ModuleType
from typing import Protocol, TypeAlias, cast

JSONPrimitive: TypeAlias = None | bool | int | float | str
JSONValue: TypeAlias = JSONPrimitive | list["JSONValue"] | dict[str, "JSONValue"]


def _load_module(name: str, path: Path) -> ModuleType:
    import importlib.util

    spec = importlib.util.spec_from_file_location(name, str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load module: {name} from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_APP_DIR = Path(__file__).resolve().parent
_ARTIFACTS = _load_module("adaptive_artifacts", _APP_DIR / "adaptive_artifacts.py")
_CHECKPOINT = _load_module("adaptive_checkpoint", _APP_DIR / "adaptive_checkpoint.py")
_POLICY = _load_module("adaptive_policy", _APP_DIR / "adaptive_policy.py")


class _Artifacts(Protocol):
    def parse_task_summary_jsonl(self, path: str | Path) -> list[dict[str, object]]: ...

    def summarize_task_summary(self, records: object) -> tuple[int, int]: ...

    def successes_trials_from_task_record(self, record: dict[str, object]) -> tuple[int, int]: ...

    def parse_run_dir_from_orchestrator_output(self, text: str) -> Path | None: ...

    def read_omp_num_threads_from_manifest(self, run_dir: Path) -> int | None: ...


class _Checkpoint(Protocol):
    def write_adaptive_state(self, checkpoint_dir: str | Path, state: object) -> Path: ...

    def write_proposals_csv(self, checkpoint_dir: str | Path, proposals: object) -> Path: ...

    def write_commands_sh(self, checkpoint_dir: str | Path, commands: object) -> Path: ...

    def load_replay_commands(self, checkpoint_dir: str | Path) -> list[str]: ...


class _Policy(Protocol):
    def make_lam1_bins(self, lam1_min: float, lam1_max: float, n_bins: int) -> list[dict[str, object]]: ...

    def allocate_budget(
        self,
        bins: object,
        successes_by_bin: dict[str, int],
        trials_by_bin: dict[str, int],
        total_budget: int,
        floor_points: int,
        alpha: float = 1.0,
        beta: float = 1.0,
    ) -> dict[str, int]: ...

    def allocate_budget_per_tanbeta(
        self,
        bins: object,
        successes_by_tb_by_bin: dict[str, dict[str, int]],
        trials_by_tb_by_bin: dict[str, dict[str, int]],
        total_budget_per_tb: int,
        floor_points: int,
        alpha: float = 1.0,
        beta: float = 1.0,
    ) -> dict[str, dict[str, int]]: ...


adaptive_artifacts = cast(_Artifacts, cast(object, _ARTIFACTS))
adaptive_checkpoint = cast(_Checkpoint, cast(object, _CHECKPOINT))
adaptive_policy = cast(_Policy, cast(object, _POLICY))


def _as_int(x: object) -> int:
    return int(cast(int | float | str, x))


def _as_float(x: object) -> float:
    return float(cast(int | float | str, x))


@dataclass(frozen=True)
class OrchestratorConfig:
    exec_path: str
    outdir: str
    lake_name: str
    campaign: str
    threads: int | None
    run_name_prefix: str
    orchestrator_dry_run: bool


def _iter_dir(checkpoint_root: Path, iter_index_1based: int) -> Path:
    return checkpoint_root / f"iter_{iter_index_1based:04d}"


def _orchestrator_path() -> Path:
    return Path(__file__).with_name("orchestrate_scans.py")


def _proposal_id(iter_index_1based: int, bin_index: int) -> str:
    return f"iter{iter_index_1based:04d}_bin{bin_index:03d}"


def _run_name(prefix: str, iter_index_1based: int, bin_index: int) -> str:
    return f"{prefix}_iter{iter_index_1based:04d}_bin{bin_index:03d}"


def _build_orchestrator_command(
    *,
    ocfg: OrchestratorConfig,
    lam1_min: float,
    lam1_max: float,
    n_lam1: int,
    run_name: str,
    n_lam1_map: str | None = None,
) -> list[str]:
    cmd: list[str] = [
        sys.executable,
        str(_orchestrator_path()),
        "--exec",
        ocfg.exec_path,
        "--outdir",
        ocfg.outdir,
        "--lake-name",
        ocfg.lake_name,
        "--campaign",
        ocfg.campaign,
        "--run-name",
        run_name,
        "--lam1-min",
        str(lam1_min),
        "--lam1-max",
        str(lam1_max),
        "--n-lam1",
        str(int(n_lam1)),
    ]

    if n_lam1_map is not None:
        cmd.extend(["--n-lam1-map", str(n_lam1_map)])

    if ocfg.threads is not None:
        cmd.extend(["--threads", str(int(ocfg.threads))])
    if ocfg.orchestrator_dry_run:
        cmd.append("--dry-run")
    return cmd


def _run_orchestrator(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


def _read_bin_metrics_from_run_dir(run_dir: Path) -> tuple[int, int]:
    records = adaptive_artifacts.parse_task_summary_jsonl(run_dir / "task_summary.jsonl")
    attempts_total, triple_ok_total = adaptive_artifacts.summarize_task_summary(records)
    return int(attempts_total), int(triple_ok_total)


def _tb_sort_key(tb_tag: str) -> tuple[int, float | str]:
    try:
        return (0, float(tb_tag))
    except ValueError:
        return (1, tb_tag)


def _sorted_tb_tags(tb_tags: Iterable[str]) -> list[str]:
    tags = [str(x) for x in tb_tags]
    tags.sort(key=_tb_sort_key)
    return tags


def _read_per_tb_bin_metrics_from_run_dir(run_dir: Path) -> tuple[dict[str, int], dict[str, int]]:
    records = adaptive_artifacts.parse_task_summary_jsonl(run_dir / "task_summary.jsonl")
    successes_by_tb: dict[str, int] = {}
    trials_by_tb: dict[str, int] = {}
    for rec in records:
        tb_tag = str(rec.get("tb_tag", ""))
        succ, trials = adaptive_artifacts.successes_trials_from_task_record(rec)
        successes_by_tb[tb_tag] = int(successes_by_tb.get(tb_tag, 0) + int(succ))
        trials_by_tb[tb_tag] = int(trials_by_tb.get(tb_tag, 0) + int(trials))
    return successes_by_tb, trials_by_tb


def _write_iteration_checkpoint(
    *,
    checkpoint_root: Path,
    iter_index_1based: int,
    state: dict[str, JSONValue],
    proposals: list[dict[str, object]],
    commands: list[list[str]],
) -> Path:
    out_dir = _iter_dir(checkpoint_root, iter_index_1based)
    _ = adaptive_checkpoint.write_adaptive_state(out_dir, state)
    _ = adaptive_checkpoint.write_proposals_csv(out_dir, proposals)
    _ = adaptive_checkpoint.write_commands_sh(out_dir, commands)
    return out_dir


def _do_replay(checkpoint_dir: Path, *, list_commands: bool) -> int:
    commands_path = checkpoint_dir / "commands.sh"
    if list_commands:
        _ = sys.stdout.buffer.write(commands_path.read_bytes())
        return 0

    commands = adaptive_checkpoint.load_replay_commands(checkpoint_dir)
    worst_rc = 0
    for line in commands:
        proc = subprocess.run(["bash", "-lc", line], check=False)
        worst_rc = max(worst_rc, int(proc.returncode))
    return int(worst_rc)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Adaptive lam1 explorer (CLI wrapper over orchestrate_scans.py)")

    _ = p.add_argument("--lam1-min", type=float, required=False, default=0.0)
    _ = p.add_argument("--lam1-max", type=float, required=False, default=12.0)
    _ = p.add_argument("--n-bins", type=int, required=False, default=12)
    _ = p.add_argument("--seed", type=int, required=False, default=0)
    _ = p.add_argument("--n-coarse", type=int, required=False, default=30)
    _ = p.add_argument("--total-budget", type=int, required=False, default=200)
    _ = p.add_argument("--floor-points", type=int, required=False, default=5)
    _ = p.add_argument("--n-iters", type=int, required=False, default=2)
    _ = p.add_argument("--checkpoint-root", type=str, required=False, default=None)

    _ = p.add_argument("--exec", dest="exec_path", type=str, required=False, default="./PhysScanWithFixings")
    _ = p.add_argument("--outdir", type=str, required=False, default="/mnt/c/Users/Asus/cern_db")
    _ = p.add_argument("--lake-name", type=str, required=False, default="dihiggs_lake")
    _ = p.add_argument("--campaign", type=str, required=False, default="scan")
    _ = p.add_argument("--threads", type=int, required=False, default=None)
    _ = p.add_argument("--run-name-prefix", type=str, required=False, default="adaptive")
    _ = p.add_argument("--orchestrator-dry-run", action="store_true")

    _ = p.add_argument("--replay", type=str, required=False, default=None, help="Replay commands from a checkpoint dir")
    _ = p.add_argument(
        "--list-commands",
        action="store_true",
        help="With --replay: print commands.sh byte-identical to stdout",
    )
    return p


class Args(Protocol):
    lam1_min: float
    lam1_max: float
    n_bins: int
    seed: int
    n_coarse: int
    total_budget: int
    floor_points: int
    n_iters: int
    checkpoint_root: str | None
    exec_path: str
    outdir: str
    lake_name: str
    campaign: str
    threads: int | None
    run_name_prefix: str
    orchestrator_dry_run: bool
    replay: str | None
    list_commands: bool


def main() -> int:
    args = cast(Args, cast(object, build_parser().parse_args()))

    if args.replay:
        return _do_replay(Path(args.replay), list_commands=bool(args.list_commands))

    if args.list_commands:
        raise SystemExit("--list-commands requires --replay")

    if args.checkpoint_root is None:
        raise SystemExit("--checkpoint-root is required (for writing iter_XXXX checkpoints)")

    if args.n_iters < 1:
        raise SystemExit("--n-iters must be >= 1")
    if args.n_bins < 1:
        raise SystemExit("--n-bins must be >= 1")
    if args.n_coarse < 1:
        raise SystemExit("--n-coarse must be >= 1")

    checkpoint_root = Path(args.checkpoint_root)
    checkpoint_root.mkdir(parents=True, exist_ok=True)

    ocfg = OrchestratorConfig(
        exec_path=str(args.exec_path),
        outdir=str(args.outdir),
        lake_name=str(args.lake_name),
        campaign=str(args.campaign),
        threads=(int(args.threads) if args.threads is not None else None),
        run_name_prefix=str(args.run_name_prefix),
        orchestrator_dry_run=bool(args.orchestrator_dry_run),
    )

    bins = adaptive_policy.make_lam1_bins(args.lam1_min, args.lam1_max, int(args.n_bins))
    rng_seed = int(args.seed)
    omp_env_threads = ocfg.threads if ocfg.threads is not None else os.environ.get("OMP_NUM_THREADS")

    successes_by_bin: dict[str, int] = {str(b["id"]): 0 for b in bins}
    trials_by_bin: dict[str, int] = {str(b["id"]): 0 for b in bins}

    successes_by_tb_by_bin: dict[str, dict[str, int]] = {}
    trials_by_tb_by_bin: dict[str, dict[str, int]] = {}
    bin_run_dirs: dict[str, str] = {}
    bin_run_names: dict[str, str] = {}

    bins_json: list[JSONValue] = []
    for b in bins:
        bins_json.append(
            {
                "id": str(b["id"]),
                "index": _as_int(b["index"]),
                "lam1_lo": _as_float(b["lam1_lo"]),
                "lam1_hi": _as_float(b["lam1_hi"]),
            }
        )

    for iter_index_1based in range(1, int(args.n_iters) + 1):
        if iter_index_1based == 1:
            budgets = {str(b["id"]): int(args.n_coarse) for b in bins}
            budgets_by_tb_by_bin: dict[str, dict[str, int]] | None = None
        else:
            if successes_by_tb_by_bin:
                tb_tags = _sorted_tb_tags(successes_by_tb_by_bin.keys())
                budgets_by_tb_by_bin = adaptive_policy.allocate_budget_per_tanbeta(
                    bins=bins,
                    successes_by_tb_by_bin=successes_by_tb_by_bin,
                    trials_by_tb_by_bin=trials_by_tb_by_bin,
                    total_budget_per_tb=int(args.total_budget),
                    floor_points=int(args.floor_points),
                )
                budgets = {}
                for b in bins:
                    bid = str(b["id"])
                    budgets[bid] = int(sum(int(budgets_by_tb_by_bin.get(tb, {}).get(bid, 0)) for tb in tb_tags))
            else:
                budgets_by_tb_by_bin = {}
                budgets = adaptive_policy.allocate_budget(
                    bins=bins,
                    successes_by_bin=successes_by_bin,
                    trials_by_bin=trials_by_bin,
                    total_budget=int(args.total_budget),
                    floor_points=int(args.floor_points),
                )

        proposals: list[dict[str, object]] = []
        commands: list[list[str]] = []
        per_iter_state_proposals: list[JSONValue] = []

        for b in bins:
            bid = str(b["id"])
            bidx = _as_int(b["index"])
            lam1_lo = _as_float(b["lam1_lo"])
            lam1_hi = _as_float(b["lam1_hi"])
            n_lam1 = int(budgets.get(bid, 0))
            if n_lam1 <= 0:
                continue

            proposal_id = _proposal_id(iter_index_1based, bidx)
            run_name = _run_name(ocfg.run_name_prefix, iter_index_1based, bidx)

            n_lam1_map: str | None = None
            if iter_index_1based >= 2 and budgets_by_tb_by_bin is not None:
                tb_tags = _sorted_tb_tags(budgets_by_tb_by_bin.keys())
                if tb_tags:
                    parts: list[str] = []
                    for tb_tag in tb_tags:
                        n_for_tb = int(budgets_by_tb_by_bin.get(tb_tag, {}).get(bid, 0))
                        parts.append(f"{tb_tag}:{n_for_tb}")
                    n_lam1_map = ",".join(parts)
                    n_lam1 = max(1, int(args.floor_points))

            cmd = _build_orchestrator_command(
                ocfg=ocfg,
                lam1_min=lam1_lo,
                lam1_max=lam1_hi,
                n_lam1=n_lam1,
                run_name=run_name,
                n_lam1_map=n_lam1_map,
            )
            commands.append(cmd)

            proposals.append(
                {
                    "proposal_id": proposal_id,
                    "status": "PENDING",
                    "command": shlex.join(cmd),
                }
            )

            result = _run_orchestrator(cmd)
            run_dir = adaptive_artifacts.parse_run_dir_from_orchestrator_output(
                (result.stdout or "") + "\n" + (result.stderr or "")
            )

            proposals[-1]["status"] = "DONE"

            per_prop_state: dict[str, JSONValue] = {
                "proposal_id": proposal_id,
                "status": "DONE",
                "command": shlex.join(cmd),
                "returncode": int(result.returncode),
                "run_name": run_name,
                "bin_id": bid,
                "bin_index": bidx,
                "lam1_min": lam1_lo,
                "lam1_max": lam1_hi,
                "n_lam1": n_lam1,
            }

            if run_dir is not None:
                run_dir_s = str(run_dir)
                per_prop_state["run_dir"] = run_dir_s
                bin_run_dirs[bid] = run_dir_s
                bin_run_names[bid] = run_name

                attempts_total, triple_ok_total = _read_bin_metrics_from_run_dir(run_dir)
                trials_by_bin[bid] = int(attempts_total)
                successes_by_bin[bid] = int(triple_ok_total)
                per_prop_state["attempts_total"] = int(attempts_total)
                per_prop_state["triple_ok_total"] = int(triple_ok_total)

                succ_by_tb, trials_by_tb = _read_per_tb_bin_metrics_from_run_dir(run_dir)
                per_prop_state["successes_by_tb"] = {str(k): int(v) for k, v in succ_by_tb.items()}
                per_prop_state["trials_by_tb"] = {str(k): int(v) for k, v in trials_by_tb.items()}

                for tb_tag, succ in succ_by_tb.items():
                    successes_by_tb_by_bin.setdefault(str(tb_tag), {})[bid] = int(succ)
                for tb_tag, tr in trials_by_tb.items():
                    trials_by_tb_by_bin.setdefault(str(tb_tag), {})[bid] = int(tr)

                omp_from_manifest = adaptive_artifacts.read_omp_num_threads_from_manifest(run_dir)
                if omp_from_manifest is not None:
                    per_prop_state["omp_num_threads_manifest"] = int(omp_from_manifest)
            else:
                per_prop_state["run_dir"] = None

            per_iter_state_proposals.append(per_prop_state)

        state: dict[str, JSONValue] = {
            "seed": rng_seed,
            "n_iters": int(args.n_iters),
            "iter_index": int(iter_index_1based),
            "lam1_min": float(args.lam1_min),
            "lam1_max": float(args.lam1_max),
            "n_bins": int(args.n_bins),
            "n_coarse": int(args.n_coarse),
            "total_budget": int(args.total_budget),
            "floor_points": int(args.floor_points),
            "omp_num_threads": (str(omp_env_threads) if omp_env_threads is not None else None),
            "orchestrator": {
                "exec": ocfg.exec_path,
                "outdir": ocfg.outdir,
                "lake_name": ocfg.lake_name,
                "campaign": ocfg.campaign,
                "threads": ocfg.threads,
                "run_name_prefix": ocfg.run_name_prefix,
                "dry_run": ocfg.orchestrator_dry_run,
            },
            "bins": bins_json,
            "budgets": {str(k): int(v) for k, v in budgets.items()},
            "budgets_by_tb_by_bin": (
                None
                if budgets_by_tb_by_bin is None
                else {
                    str(tb): {str(bid): int(n) for bid, n in per_bin.items()}
                    for tb, per_bin in budgets_by_tb_by_bin.items()
                }
            ),
            "bin_run_dirs": dict(bin_run_dirs),
            "bin_run_names": dict(bin_run_names),
            "trials_by_bin": dict(trials_by_bin),
            "successes_by_bin": dict(successes_by_bin),
            "trials_by_tb_by_bin": {
                str(tb): {str(bid): int(n) for bid, n in per_bin.items()}
                for tb, per_bin in trials_by_tb_by_bin.items()
            },
            "successes_by_tb_by_bin": {
                str(tb): {str(bid): int(n) for bid, n in per_bin.items()}
                for tb, per_bin in successes_by_tb_by_bin.items()
            },
            "proposals": per_iter_state_proposals,
        }

        _ = _write_iteration_checkpoint(
            checkpoint_root=checkpoint_root,
            iter_index_1based=iter_index_1based,
            state=state,
            proposals=proposals,
            commands=commands,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
