"""Drift guard for external-tool versions and the physics authority file.

The repositories are separate and cannot see each other in CI. `dihiggs` owns
the canonical contract; downstream repositories verify an exact cached copy
against a source commit and SHA-256 sidecar.

  * EXPECTED_VERSIONS  -- every repo's external_tools.lock.yaml must declare
    these same tool versions.
  * PINNED_CONVENTIONS_SHA256 -- the authoritative file in this repository must
    match its manifest. Downstream migrations pin this value independently.
"""

import hashlib
import os

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST = os.path.join(REPO_ROOT, "external_tools.lock.yaml")
CONVENTIONS = os.path.join(REPO_ROOT, "conventions", "physics_conventions.yaml")

# --- shared cross-repo pins (keep identical in all three repos) ------------
EXPECTED_VERSIONS = {
    "2HDMC": "1.8",
    "HiggsTools": "v1.2",
    "HiggsBounds_dataset": "v1.7",
    "HiggsSignals_dataset": "v1.1",
}
PINNED_CONVENTIONS_SHA256 = (
    "f7cea46d28f328f9b66a5de3b9d76fbbb65948c384f813828fb1462dfe1ee147"
)
CONVENTIONS_SCHEMA_VERSION = "physics_conventions_v3"


def _load_manifest():
    yaml = pytest.importorskip("yaml")
    with open(MANIFEST) as fh:
        return yaml.safe_load(fh)


def test_manifest_present_and_parses():
    assert os.path.exists(MANIFEST), "missing external_tools.lock.yaml"
    manifest = _load_manifest()
    assert manifest["schema_version"] == "external_tools_v1"
    assert "tools" in manifest


def test_declared_versions_match_shared_pins():
    """Every tool version must equal the cross-repo pinned value; a bump in one
    repo without the others fails here."""
    tools = _load_manifest()["tools"]
    for name, expected in EXPECTED_VERSIONS.items():
        assert name in tools, "manifest missing tool %s" % name
        assert tools[name]["version"] == expected, (
            "%s version %r != pinned %r" % (name, tools[name]["version"], expected)
        )


def test_conventions_sha256_matches_pin():
    """The authoritative conventions file must match its SHA-256 pin."""
    assert os.path.exists(CONVENTIONS)
    with open(CONVENTIONS, "rb") as fh:
        actual = hashlib.sha256(fh.read()).hexdigest()
    assert actual == PINNED_CONVENTIONS_SHA256, (
        "conventions/physics_conventions.yaml sha256 %s != pinned %s "
        "(did the authority file drift?)" % (actual, PINNED_CONVENTIONS_SHA256)
    )


def test_manifest_conventions_block_matches_file():
    """The manifest's SHA-256/schema must match the actual authority file, so
    the manifest cannot silently lag the conventions file."""
    conv = _load_manifest()["conventions"]
    assert conv["authority_repository"] == "fbientrigo/dihiggs"
    assert conv["authority"] is True
    assert conv["sha256"] == PINNED_CONVENTIONS_SHA256
    assert conv["schema_version"] == CONVENTIONS_SCHEMA_VERSION
    with open(CONVENTIONS, "rb") as fh:
        actual = hashlib.sha256(fh.read()).hexdigest()
    assert conv["sha256"] == actual


def test_vendored_paths_exist_when_declared():
    """If a tool declares a vendored_path, it must exist in this repo (catches a
    renamed/moved vendored tree that would break the build)."""
    tools = _load_manifest()["tools"]
    for name, meta in tools.items():
        path = meta.get("vendored_path")
        if path:
            assert os.path.exists(os.path.join(REPO_ROOT, path)), (
                "%s vendored_path %r does not exist" % (name, path)
            )
