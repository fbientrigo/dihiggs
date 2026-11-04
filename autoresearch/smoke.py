"""#19 end-to-end smoke (`make autoresearch-smoke`): REAL boundary + a throwaway local-dir lake, 32 points, ~1 min.

1. build the evaluator  2. campaign runs boundary and commits to the lake  3. a second campaign repeats nothing
4. kill -9 mid-action, restart reconciles: one commit, nothing lost  5. the TUI shows budget left, actions, coverage.
6. an mh action brings fresh rows every attempt, so repeating it never ends the campaign SATURATED.
Passing = TESTED (not VALIDATED). Nothing touches the real lake or any physics number.
"""
from __future__ import annotations

import asyncio
import json
import os
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from autoresearch.harness.journal import Journal
from autoresearch.tui.app import AutoresearchApp
from autoresearch.tui.services.monitor import journal_text, journal_view

ROOT = Path(__file__).resolve().parents[1]  # dihiggs checkout
BOUNDARY = ROOT.parent / "dihiggs_boundary"
CONFIG = json.loads((ROOT / "autoresearch/configs/boundary_lambda7_smoke.json").read_text())
BUDGET = ["--cpu-hours", "0.6", "--disk-gb", "1", "--wall-hours", "0.5"]


def check(ok: bool, msg: str) -> None:
    print(("ok   " if ok else "FAIL ") + msg, flush=True)
    if not ok:
        raise SystemExit(1)


def supervisor(work: Path, name: str, env: dict, config: dict = CONFIG) -> subprocess.Popen:
    cfg = {**config, "paths": {**config["paths"], "outdir": str(work / name)}}
    (work / f"{name}.json").write_text(json.dumps(cfg))
    return subprocess.Popen([sys.executable, "-m", "autoresearch.run_supervisor", str(work / f"{name}.json"), *BUDGET],
                            cwd=ROOT, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def closed(work: Path, name: str) -> list[dict]:
    return [p for _, _, p in Journal(work / name / "journal.db").events("closed")]


def commits(env: dict) -> list[str]:
    out = subprocess.run([sys.executable, str(ROOT.parent / "lake/lake.py"), "ls"], env=env, capture_output=True, text=True, check=True).stdout
    return [line.split("\t")[0] for line in out.splitlines() if "autoresearch" in line]


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="autoresearch-smoke-") as tmp:
        work = Path(tmp)
        (work / "remote").mkdir()
        env = {**os.environ, "DIHIGGS_ROOT": str(ROOT), "LAKE_REMOTE": str(work / "remote"), "LAKE_CACHE": str(work / "cache")}

        build = subprocess.run(["bash", str(BOUNDARY / "scripts/build_evaluate_point.sh")], capture_output=True, text=True)
        check(build.returncode == 0 and (BOUNDARY / "build/bin/evaluate_point").exists(), "1 evaluator compiles (boundary scripts/build_evaluate_point.sh)")

        t0 = time.time()
        check(supervisor(work, "run1", env).wait(timeout=900) == 0, f"2 campaign ran boundary ({time.time() - t0:.0f}s)")
        c1 = closed(work, "run1")
        check(len(c1) == 1 and c1[0]["fresh_rows"] == 32 and c1[0]["status"] == "success", f"2 {c1[0]['fresh_rows']} fresh points computed")
        check(len(commits(env)) == 1 and c1[0]["lake_commit"] in commits(env), "2 result uploaded to the lake (1 commit)")
        kept = sorted(str(f.relative_to(work / "run1/actions/a0000")) for f in (work / "run1/actions/a0000/out").rglob("*") if f.is_file())
        check(kept == ["out/survivor_adaptive_12520_tanbeta_v2/manifest.json"], f"2 verified commit -> local outputs deleted, kept {kept}")

        check(supervisor(work, "run2", env).wait(timeout=900) == 0, "3 second campaign (same lake, empty journal) ran")
        c2 = closed(work, "run2")
        check(c2[0]["fresh_rows"] == 0 and c2[0]["seeded_rows"] > 0 and len(commits(env)) == 1,
              f"3 repeated nothing: fresh_rows=0, seeded_rows={c2[0]['seeded_rows']}, still 1 commit")

        (work / "remote2").mkdir()
        env2 = {**env, "LAKE_REMOTE": str(work / "remote2"), "LAKE_CACHE": str(work / "cache2")}
        proc = supervisor(work, "run3", env2)
        for _ in range(1200):
            db = work / "run3" / "journal.db"
            if db.exists() and Journal(db).get("started:a0000"):
                break
            time.sleep(0.05)
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait()
        check(Journal(work / "run3" / "journal.db").open_attempts() == ["a0000"], "4 daemon kill -9'd mid-action: attempt a0000 open in the journal")
        check(supervisor(work, "run3", env2).wait(timeout=900) == 0, "4 restart ran")
        j = Journal(work / "run3" / "journal.db")
        c3 = closed(work, "run3")
        check(not j.open_attempts() and c3[0]["fresh_rows"] == 32 and len(commits(env2)) == 1, "4 reconciled: attempt closed, 32 rows, exactly 1 commit, no relaunch")

        async def snapshot() -> str:
            app = AutoresearchApp()
            async with app.run_test(size=(140, 50)) as pilot:
                app.load_campaign(work / "run3")
                await pilot.pause()
                return str(app.query_one("#journal-info").render())

        text = asyncio.run(snapshot())
        check(text == journal_text(journal_view(work / "run3")) and "Budget remaining" in text and "Coverage: 1/" in text, "5 TUI shows budget, actions, coverage:")
        print("\n".join("     " + line for line in text.splitlines()))

        runtime = {**CONFIG["runtime"], "max_rounds": 2, "boundary_overrides": {**CONFIG["runtime"]["boundary_overrides"],
                                                                               "stochastic": {"evaluation_cap": 32, "epsilon": 0.05}}}
        mh = {**CONFIG, "runtime": runtime, "arms": [{"id": "mh-tb1e7-l7", "action": "mh|tb=10000000.0|l7_0p01_mh125p20", "disk_mb_per_action": 32}]}
        check(supervisor(work, "run4", env, mh).wait(timeout=900) == 0, "6 mh campaign ran")
        c4 = closed(work, "run4")
        state = json.loads((work / "run4" / "campaign_state.json").read_text())["campaign_state"]
        check([c["fresh_rows"] for c in c4] == [32, 32] and state != "SATURATED", f"6 mh repeated: fresh_rows {[c['fresh_rows'] for c in c4]}, state {state}")
    print("autoresearch-smoke TESTED")


if __name__ == "__main__":
    main()
