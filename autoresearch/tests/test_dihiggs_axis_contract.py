from __future__ import annotations

import copy
import sys
from typing import cast

from autoresearch.harness.dihiggs_axis_contract import (
    bin_lam1,
    bin_mphi,
    decode_cell_id,
    encode_cell_id,
    substitute_placeholders,
)


def _test_config() -> dict[str, object]:
    return {
        "search": {
            "lam1": {"min": -20.0, "max": 20.0, "n_bins": 40},
            "mphi": {"min": 0.0, "max": 1000.0, "n_bins": 50},
        }
    }


def test_canonical_encoding_backward_compatible_bin_alias() -> None:
    assert encode_cell_id({"tb": 10000, "lam1_bin": 5}) == "tb=10000|bin=5"


def test_extended_encoding_with_additional_axes() -> None:
    assert (
        encode_cell_id({"tb": 10000, "lam1_bin": 5, "mphi_bin": 3})
        == "tb=10000|bin=5|mphi_bin=3"
    )


def test_roundtrip_encode_decode_encode_is_stable() -> None:
    original = "tb=10000|bin=5|explorer=adaptive|mphi_bin=3|strategy=adaptive-v1"
    assert encode_cell_id(decode_cell_id(original)) == original


def test_backward_compat_decode_old_cell_id() -> None:
    parsed = decode_cell_id("tb=5000|bin=7")
    assert parsed["tb"] == 5000
    assert parsed["lam1_bin"] == 7


def test_lam1_binning_boundaries_and_clamping() -> None:
    cfg = _test_config()
    assert bin_lam1(-20.0, cfg) == 0
    assert bin_lam1(20.0, cfg) == 39
    assert bin_lam1(-100.0, cfg) == 0
    assert bin_lam1(100.0, cfg) == 39
    assert bin_lam1(-19.0, cfg) == 1


def test_mphi_binning_boundaries_and_clamping() -> None:
    cfg = _test_config()
    assert bin_mphi(0.0, cfg) == 0
    assert bin_mphi(1000.0, cfg) == 49
    assert bin_mphi(-1.0, cfg) == 0
    assert bin_mphi(1001.0, cfg) == 49
    assert bin_mphi(20.0, cfg) == 1


def test_placeholder_substitution_supports_required_keys_and_recursion() -> None:
    cfg = {
        "runtime": {"python_exe": "{sys_executable}"},
        "paths": {
            "checkpoint_root": "{outdir}/checkpoints/{campaign}",
            "leaf": "{checkpoint_root}/{track_id}",
        },
        "cmd": [
            "{python}",
            "--checkpoint-root",
            "{checkpoint_root}",
            "--tb",
            "{tb}",
            "--iter",
            "{iter}",
            "--phys",
            "{phys_exec}",
            "--hb-env",
            "{hb_dataset_env}",
            "--hs-env",
            "{hs_dataset_env}",
            "--lake",
            "{lake_name}",
        ],
    }
    runtime_ctx = {
        "python": "{sys_executable}",
        "sys_executable": sys.executable,
        "campaign": "camp-1",
        "outdir": "/tmp/out",
        "checkpoint_root": "{outdir}/ckpt/{campaign}",
        "lake_name": "events.jsonl",
        "threads": 8,
        "tb": 10000,
        "iter": 12,
        "track_id": "track-abc",
        "phys_exec": "/opt/PhysScanWithFixings",
        "hb_dataset_env": "HB_DATASET",
        "hs_dataset_env": "HS_DATASET",
    }
    runtime_ctx_copy = copy.deepcopy(runtime_ctx)
    cfg_copy = copy.deepcopy(cfg)

    out = substitute_placeholders(cfg, runtime_ctx)
    runtime_out = cast(dict[str, object], out["runtime"])
    paths_out = cast(dict[str, object], out["paths"])
    cmd_out = cast(list[object], out["cmd"])

    assert runtime_out["python_exe"] == sys.executable
    assert paths_out["checkpoint_root"] == "/tmp/out/checkpoints/camp-1"
    assert paths_out["leaf"] == "/tmp/out/ckpt/camp-1/track-abc"
    assert cmd_out[-1] == "events.jsonl"

    assert cfg == cfg_copy
    assert runtime_ctx == runtime_ctx_copy


def test_placeholder_substitution_missing_key_raises_clear_error() -> None:
    cfg = {"cmd": ["{python}", "{missing_one}"]}
    runtime_ctx = {"python": sys.executable}

    try:
        _ = substitute_placeholders(cfg, runtime_ctx)
        raise AssertionError("Expected ValueError")
    except ValueError as exc:
        assert "missing_one" in str(exc)
