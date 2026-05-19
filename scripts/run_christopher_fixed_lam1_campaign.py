#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


DEFAULT_CAMPAIGN = "christopher_fixed_lam1_2026apr"
DEFAULT_LAMBDA1_VALUES = [0.4, 1.0, 12.566370614359172]
DEFAULT_LAMBDA6_VALUES = [0.01, 0.1]


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def format_tag(x: float, ndp: int = 8) -> str:
    return f"{x:.{ndp}f}".replace("-", "m").replace(".", "p")


def parse_float_csv(s: str) -> List[float]:
    out: List[float] = []
    for p in s.split(","):
        p = p.strip()
        if p:
            out.append(float(p))
    return out


def normalize_campaign_name(name: str) -> str:
    name = name.strip()
    if name.startswith("campaign="):
        return name[len("campaign="):]
    return name


def sanitize_for_path(s: str) -> str:
    return s.replace(".", "p").replace("-", "m")


def expected_run_dir(outdir: str, lake_name: str, campaign: str, sin_ba: float, lam6: float, lambda7: float, mA: float, run_name: str) -> Path:
    fixed = (
        f"fixed_sinba={sanitize_for_path(f'{sin_ba:.4f}')}"
        f"_l6={sanitize_for_path(f'{lam6:.4f}')}"
        f"_l7={sanitize_for_path(f'{lambda7:.4f}')}"
        f"_mA={sanitize_for_path(f'{mA:.1f}')}"
    )
    return Path(outdir) / lake_name / f"campaign={campaign}" / fixed / run_name


def read_orchestrator_status(run_dir: Path) -> Dict[str, int]:
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.exists():
        return {"tasks_failed": 1}
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        summary = payload.get("summary") or {}
        return {
            "tasks_failed": int(summary.get("tasks_failed", 0)),
            "tasks_ok": int(summary.get("tasks_ok", 0)),
            "tasks_total": int(summary.get("tasks_total", 0)),
        }
    except Exception:
        return {"tasks_failed": 1}


@dataclass(frozen=True)
class WrapperConfig:
    campaign: str
    outdir: str
    lake_name: str
    exec_path: str
    threads: int | None
    dry_run: bool
    force: bool
    mphi_min: float
    mphi_max: float
    n_mphi: int
    tanbeta: float
    mA: float
    sin_ba: float
    lambda7: float
    lambda1_values: List[float]
    lambda6_values: List[float]


def build_orchestrator_command(cfg: WrapperConfig, lam1: float, lam6: float, run_name: str, orchestrator_path: Path) -> List[str]:
    cmd: List[str] = [
        sys.executable,
        str(orchestrator_path),
        "--exec",
        cfg.exec_path,
        "--campaign",
        cfg.campaign,
        "--run-name",
        run_name,
        "--outdir",
        cfg.outdir,
        "--lake-name",
        cfg.lake_name,
        "--mphi-min",
        str(cfg.mphi_min),
        "--mphi-max",
        str(cfg.mphi_max),
        "--n-mphi",
        str(cfg.n_mphi),
        "--lam1-min",
        str(lam1),
        "--lam1-max",
        str(lam1),
        "--n-lam1",
        "1",
        "--mA",
        str(cfg.mA),
        "--sin-ba",
        str(cfg.sin_ba),
        "--lambda6",
        str(lam6),
        "--lambda7",
        str(cfg.lambda7),
        "--tanbeta",
        str(cfg.tanbeta),
    ]
    if cfg.threads is not None:
        cmd.extend(["--threads", str(cfg.threads)])
    if cfg.force:
        cmd.append("--force")
    if cfg.dry_run:
        cmd.append("--dry-run")
    return cmd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run fixed-lambda1 Christopher-style campaign through existing orchestrate_scans.py")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--force", action="store_true")
    p.add_argument("--threads", type=int, default=None)
    p.add_argument("--outdir", type=str, default="/mnt/c/Users/Asus/cern_db")
    p.add_argument("--lake-name", type=str, default="dihiggs_lake")

    repo_root = Path(__file__).resolve().parents[1]
    p.add_argument("--exec", dest="exec_path", type=str, default=str((repo_root / "dihiggs/app/PhysScanWithFixings").resolve()))
    p.add_argument("--campaign", type=str, default=DEFAULT_CAMPAIGN)

    p.add_argument("--lambda1-values", type=str, default=",".join(str(x) for x in DEFAULT_LAMBDA1_VALUES))
    p.add_argument("--lambda6-values", type=str, default=",".join(str(x) for x in DEFAULT_LAMBDA6_VALUES))

    p.add_argument("--mphi-min", type=float, default=125.0)
    p.add_argument("--mphi-max", type=float, default=300.0)
    p.add_argument("--n-mphi", type=int, default=400)

    p.add_argument("--tanbeta", type=float, default=10000.0)
    p.add_argument("--mA", type=float, default=300.0)
    p.add_argument("--sin-ba", type=float, default=1.0)
    p.add_argument("--lambda7", type=float, default=0.0)
    return p


def main() -> int:
    args = build_parser().parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    orchestrator_path = (repo_root / "dihiggs/app/orchestrate_scans.py").resolve()
    if not orchestrator_path.exists():
        print(f"[ERROR] Orchestrator not found: {orchestrator_path}", file=sys.stderr)
        return 2

    campaign = normalize_campaign_name(args.campaign)
    lambda1_values = parse_float_csv(args.lambda1_values)
    lambda6_values = parse_float_csv(args.lambda6_values)

    cfg = WrapperConfig(
        campaign=campaign,
        outdir=args.outdir,
        lake_name=args.lake_name,
        exec_path=args.exec_path,
        threads=args.threads,
        dry_run=args.dry_run,
        force=args.force,
        mphi_min=args.mphi_min,
        mphi_max=args.mphi_max,
        n_mphi=args.n_mphi,
        tanbeta=args.tanbeta,
        mA=args.mA,
        sin_ba=args.sin_ba,
        lambda7=args.lambda7,
        lambda1_values=lambda1_values,
        lambda6_values=lambda6_values,
    )

    out_base = repo_root / "scripts" / "out" / campaign
    out_base.mkdir(parents=True, exist_ok=True)
    wrapper_log = out_base / "wrapper.log"
    wrapper_manifest = out_base / "wrapper_manifest.json"

    jobs = []
    for lam1 in cfg.lambda1_values:
        for lam6 in cfg.lambda6_values:
            run_name = (
                f"branchcurve_l6_{format_tag(lam6)}"
                f"_lam1_{format_tag(lam1)}"
                f"_tb_{format_tag(cfg.tanbeta)}"
                "_fixed_lambda1_christopher"
            )
            cmd = build_orchestrator_command(cfg, lam1=lam1, lam6=lam6, run_name=run_name, orchestrator_path=orchestrator_path)
            jobs.append({"lam1": lam1, "lam6": lam6, "run_name": run_name, "cmd": cmd})

    with wrapper_log.open("w", encoding="utf-8") as logf:
        logf.write(f"[{utc_now_iso()}] wrapper start\n")
        for i, job in enumerate(jobs, start=1):
            line = f"[{i}/{len(jobs)}] {' '.join(job['cmd'])}"
            print(line)
            logf.write(line + "\n")

    manifest = {
        "created_utc": utc_now_iso(),
        "campaign": cfg.campaign,
        "wrapper_log": str(wrapper_log),
        "defaults": {
            "mphi_min": cfg.mphi_min,
            "mphi_max": cfg.mphi_max,
            "n_mphi": cfg.n_mphi,
            "tanbeta": cfg.tanbeta,
            "mA": cfg.mA,
            "sin_ba": cfg.sin_ba,
            "lambda7": cfg.lambda7,
            "n_lam1": 1,
        },
        "jobs": [],
    }

    failures = 0
    for job in jobs:
        rc = 0
        orchestrator_summary: Dict[str, int] = {}
        if not cfg.dry_run:
            proc = subprocess.run(job["cmd"], text=True)
            rc = proc.returncode
            run_dir = expected_run_dir(
                outdir=cfg.outdir,
                lake_name=cfg.lake_name,
                campaign=cfg.campaign,
                sin_ba=cfg.sin_ba,
                lam6=job["lam6"],
                lambda7=cfg.lambda7,
                mA=cfg.mA,
                run_name=job["run_name"],
            )
            orchestrator_summary = read_orchestrator_status(run_dir)
            if orchestrator_summary.get("tasks_failed", 0) > 0:
                rc = 1
        status = "dry_run" if cfg.dry_run else ("ok" if rc == 0 else "fail")
        if rc != 0:
            failures += 1
        manifest["jobs"].append(
            {
                "lam1": job["lam1"],
                "lam6": job["lam6"],
                "run_name": job["run_name"],
                "cmd": job["cmd"],
                "returncode": rc,
                "status": status,
                "orchestrator_summary": orchestrator_summary,
            }
        )

    manifest["summary"] = {
        "total_jobs": len(jobs),
        "failed_jobs": failures,
        "status": "ok" if failures == 0 else "failed",
    }
    wrapper_manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    if failures:
        print(f"[ERROR] {failures} orchestrator command(s) failed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
