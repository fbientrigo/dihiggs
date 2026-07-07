import importlib
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAKE_PIPELINE = os.path.join(REPO_ROOT, "mlpython", "lake_pipeline")
CONVENTIONS = os.path.join(REPO_ROOT, "conventions", "physics_conventions.yaml")

sys.path.insert(0, LAKE_PIPELINE)
physics_conventions = importlib.import_module("physics_conventions")


def test_conventions_file_present():
    assert os.path.exists(CONVENTIONS)


def test_pinned_values():
    assert physics_conventions.HBAR_C_GEV_MM == 1.973269804e-13
    assert physics_conventions.C_MM_PER_NS == 299.792458


def test_loaded_matches_pinned_fallback():
    assert physics_conventions.HBAR_C_GEV_MM == physics_conventions._HBAR_C_GEV_MM_PINNED
    assert physics_conventions.C_MM_PER_NS == physics_conventions._C_MM_PER_NS_PINNED


def test_loaded_matches_conventions_yaml():
    yaml = pytest.importorskip("yaml")
    with open(CONVENTIONS) as fh:
        conv = yaml.safe_load(fh)
    assert physics_conventions.HBAR_C_GEV_MM == float(conv["hbar_c_gev_mm"])
    assert physics_conventions.C_MM_PER_NS == float(conv["c_mm_per_ns"])


def test_ctau_helper():
    w = 2e-13
    assert physics_conventions.ctau_mm_from_width_gev(w) == (
        physics_conventions.HBAR_C_GEV_MM / w
    )
