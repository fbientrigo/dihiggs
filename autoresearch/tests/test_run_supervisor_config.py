from __future__ import annotations

import json
from pathlib import Path

from .. import run_supervisor

CONFIGS = Path(run_supervisor.__file__).parent / "configs"


def test_repo_root_derived_from_config_location_and_env(tmp_path, monkeypatch) -> None:
    cfg = tmp_path / "dihiggs" / "autoresearch" / "configs" / "c.json"
    cfg.parent.mkdir(parents=True)
    cfg.write_text(json.dumps({"paths": {"outdir": "{repo_root}/runs/x"}}))

    monkeypatch.delenv("DIHIGGS_ROOT", raising=False)
    paths = run_supervisor._load_config(cfg)["paths"]
    assert paths["repo_root"] == str(tmp_path / "dihiggs")
    assert paths["outdir"] == str(tmp_path / "dihiggs" / "runs" / "x")

    monkeypatch.setenv("DIHIGGS_ROOT", "/elsewhere")
    assert run_supervisor._load_config(cfg)["paths"]["repo_root"] == "/elsewhere"


def test_shipped_configs_have_no_machine_paths() -> None:
    for cfg in CONFIGS.glob("*.json"):
        assert "/home/" not in cfg.read_text(), cfg.name
