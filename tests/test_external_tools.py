"""The authoritative conventions file must match the sha256 pinned in external_tools.lock.yaml.

Cross-repo tool-version agreement is checked by the root scripts/check_workspace.py.
"""

import hashlib
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
PINNED_CONVENTIONS_SHA256 = "4f3a938ed54dbc779d25f691ce1dd06fe17af69751105dc3b69ccb6b7d92fff6"


def test_conventions_sha256_matches_pin():
    conv = yaml.safe_load((ROOT / "external_tools.lock.yaml").read_text())["conventions"]
    actual = hashlib.sha256((ROOT / "conventions/physics_conventions.yaml").read_bytes()).hexdigest()
    assert (conv["sha256"], conv["schema_version"]) == (PINNED_CONVENTIONS_SHA256, "physics_conventions_v3")
    assert actual == PINNED_CONVENTIONS_SHA256, "conventions/physics_conventions.yaml drifted from its pin"
