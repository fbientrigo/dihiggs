import json

import pytest
import yaml

from autoresearch.harness import boundary_actions as ba
from autoresearch.harness.dihiggs_runner import DiHiggsRunner


def _boundary(tmp_path):
    root = tmp_path / "dihiggs_boundary"
    for spec in ba.STRATEGIES.values():
        path = root / spec["config"]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump({"campaign_id": path.stem, "dihiggs_root": "../../dihiggs", "tan_beta_slices": [100.0, 1e7],
                                        "physics": {"lambda7": 0.0, "mh_GeV": 125.13}}))
    return root


def test_menu_registered_with_digest_and_out_of_menu_rejected(tmp_path):
    root = _boundary(tmp_path)
    menu = ba.build_menu(root)
    assert {a["strategy"] for a in menu} == set(ba.STRATEGIES) and all(a["regime"]["lambda7"] != 0 for a in menu)
    sha = ba.register(tmp_path / "menu.json", menu)
    assert json.loads((tmp_path / "menu.json").read_text())["digest"] == sha == ba.register(tmp_path / "menu.json", menu)
    with pytest.raises(ValueError, match="menu changed"):
        ba.register(tmp_path / "menu.json", menu[1:])
    with pytest.raises(ValueError, match="not in the registered menu"):
        ba.lookup(menu, "survivor_v2|tb=100.0|l7_0_mh125p20")

    cfg = {"campaign_id": "t", "paths": {"repo_root": str(tmp_path / "dihiggs")}, "arms": [{"id": "a", "action": "free_cmd"}]}
    with pytest.raises(ValueError, match="not in the registered menu"):
        DiHiggsRunner(cfg, str(tmp_path / "out"))


def test_lambda7_0p01_reaches_the_boundary_config(tmp_path):
    root = _boundary(tmp_path)
    action = ba.lookup(ba.build_menu(root), "survivor_v2|tb=10000000.0|l7_0p01_mh125p20")
    cfg = ba.boundary_config(action, root, {"mesh": {"root_cells": 1}})
    assert cfg["physics"]["lambda7"] == 0.01 and cfg["physics"]["mh_GeV"] == 125.2 and cfg["mesh"]["root_cells"] == 1
    assert ba.command(action, "c.yaml", "out", "py")[-2:] == ["--tan-beta", "10000000.0"]
    with pytest.raises(ValueError, match="may not set physics"):
        ba.boundary_config(action, root, {"physics": {"lambda7": 0.0}})


def test_action_reads_the_lake_first_and_commits_only_fresh_points(tmp_path, monkeypatch):
    import sys
    from autoresearch.tests.test_journal_budget import KEY, commits, launches, world

    cfg = world(tmp_path)
    runner = DiHiggsRunner(cfg, str(tmp_path / "out"))
    first = runner._run_action("a", cfg["arms"][0], False)
    argv = json.loads((tmp_path / "launches.jsonl").read_text().splitlines()[0])
    sent = yaml.safe_load(open(argv[argv.index("--config") + 1]))
    assert sent["physics"]["lambda7"] == 0.01 and "--explored" in argv and first["lake_commit"] == "cid0"
    assert first["seeded_rows"] == 5 and first["valid_rows"] == 1
    monkeypatch.setenv("FAKE_FRESH", "0")
    second = runner._run_action("a", cfg["arms"][0], False)
    assert second["fresh_rows"] == 0 and second["lake_commit"] is None and len(commits(tmp_path)) == 1 and launches(tmp_path) == 2
    with pytest.raises(ValueError, match="would repeat lake points"):
        runner._run_action("a", {"action": "adaptive_m2_mesh|tb=all|l7_0p01_mh125p20"}, False)

    regime = ba.REGIMES["l7_0p01_mh125p20"]
    fail = lambda msg: [sys.executable, "-c", f"import sys; sys.stderr.write({msg!r}); sys.exit(1)"]
    assert ba.explored_export(fail("no point tables in the lake"), regime, tmp_path / "e.csv") is False
    with pytest.raises(RuntimeError, match="not launching"):
        ba.explored_export(fail("rclone: network down"), regime, tmp_path / "e.csv")
