"""Drift guards for the canonical SM-like Higgs mass convention.

`mh = 125.13` used to be a bare literal in five orchestrator files, and they
disagreed: `engines/m2.py` hard-coded it into `axis_metadata()` while
`manifest.py` honoured an explicit `--mh`, so a run launched with `--mh X`
wrote a manifest claiming the canonical value. These tests pin the value to its
single owner and lock the disagreement shut.
"""

import os

import pytest

from dihiggs.app.orchestrator import conventions
from dihiggs.app.orchestrator.cli import build_parser
from dihiggs.app.orchestrator.engines.m2 import M2Engine


def _yaml_conventions():
    yaml = pytest.importorskip("yaml")
    assert os.path.exists(conventions.CONVENTIONS_PATH)
    with open(conventions.CONVENTIONS_PATH) as fh:
        return yaml.safe_load(fh)


def test_canonical_mh_matches_conventions_file():
    """The module must serve the YAML's value, not its pinned fallback."""
    conv = _yaml_conventions()
    assert conv["sm_like_higgs"]["m_h_GeV"] == conventions.M_H_GEV_TEXT
    assert float(conv["sm_like_higgs"]["m_h_GeV"]) == conventions.M_H_GEV


def test_pinned_fallback_agrees_with_loaded_value():
    """Fallback and file can never silently disagree."""
    assert conventions.M_H_GEV_TEXT == conventions._M_H_GEV_TEXT_PINNED


def test_mh_is_a_decimal_string_not_a_float():
    """float64 repr of 125.20 is 1.25200000000000003e+02, which is not the
    byte form UFO/SLHA cards need. The convention must carry the string."""
    conv = _yaml_conventions()
    assert isinstance(conv["sm_like_higgs"]["m_h_GeV"], str)
    assert "%.17e" % conventions.M_H_GEV != "1.25200000000000000e+02"


def test_convention_records_provenance():
    conv = _yaml_conventions()["sm_like_higgs"]
    assert conv["source"]
    assert conv["source_url"].startswith("https://")
    superseded = {entry["value"] for entry in conv["supersedes"]}
    assert {"125.13", "125.09"} <= superseded


def test_orchestrator_mh_default_is_the_canonical_value():
    """The CLI default must come from the conventions file, not a literal."""
    parser = build_parser()
    default = parser.get_default("mh")
    assert default == conventions.M_H_GEV


def test_axis_metadata_reports_the_canonical_convention():
    meta = M2Engine().axis_metadata()["mass_convention"]
    assert meta["mh_GeV"] == conventions.M_H_GEV
    assert meta["source"] == conventions.M_H_SOURCE


def test_explicit_mh_override_is_labelled_not_silently_canonical():
    """Regression: an explicit --mh must never be reported as the canonical
    convention. Previously engines/m2.py claimed 125.13 regardless."""
    block = conventions.mass_convention_block(125.13)
    assert block["mh_GeV"] == 125.13
    assert block["source"] == "explicit --mh override"
    assert block["canonical_mh_GeV"] == conventions.M_H_GEV
