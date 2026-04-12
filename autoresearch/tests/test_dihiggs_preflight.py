from __future__ import annotations
# pyright: reportMissingImports=false, reportUnknownVariableType=false, reportUnknownMemberType=false, reportUnknownParameterType=false, reportMissingParameterType=false

from _pytest.monkeypatch import MonkeyPatch
from unittest.mock import patch

from ..harness.dihiggs_preflight import (
    check_2hdmc_lib,
    check_datasets,
    check_phys_exec,
    run_all_preflight_checks,
)


def _base_config() -> dict[str, object]:
    return {
        "paths": {"repo_root": "/repo"},
        "dihiggs": {
            "phys_exec": "/repo/dihiggs/app/PhysScanWithFixings",
            "hb_dataset_env": "HB_DATASET",
            "hs_dataset_env": "HS_DATASET",
        },
    }


def test_check_phys_exec_pass_when_binary_exists_and_executable() -> None:
    cfg = _base_config()
    with patch("os.path.exists", return_value=True), patch("os.access", return_value=True):
        result = check_phys_exec(cfg)

    assert result["check"] == "phys_exec"
    assert result["status"] == "pass"
    assert "executable" in result["message"]


def test_check_phys_exec_fail_when_binary_missing() -> None:
    cfg = _base_config()
    with patch("os.path.exists", return_value=False):
        result = check_phys_exec(cfg)

    assert result == {
        "check": "phys_exec",
        "status": "fail",
        "message": "PhysScanWithFixings not found at /repo/dihiggs/app/PhysScanWithFixings",
    }


def test_check_phys_exec_fail_when_not_executable() -> None:
    cfg = _base_config()
    with patch("os.path.exists", return_value=True), patch("os.access", return_value=False):
        result = check_phys_exec(cfg)

    assert result["check"] == "phys_exec"
    assert result["status"] == "fail"
    assert "not executable" in result["message"]


def test_check_phys_exec_exception_safe_returns_fail() -> None:
    cfg = _base_config()
    with patch("os.path.exists", side_effect=PermissionError("denied")):
        result = check_phys_exec(cfg)

    assert result["check"] == "phys_exec"
    assert result["status"] == "fail"
    assert "Unexpected error" in result["message"]


def test_check_datasets_pass_when_both_env_set_and_dirs_exist(monkeypatch: MonkeyPatch) -> None:
    cfg = _base_config()
    monkeypatch.setenv("HB_DATASET", "/datasets/HBDataset")
    monkeypatch.setenv("HS_DATASET", "/datasets/HSDataset")

    def exists_side_effect(path: object) -> bool:
        return isinstance(path, str) and path in {"/datasets/HBDataset", "/datasets/HSDataset"}

    with patch("os.path.exists", side_effect=exists_side_effect), patch("os.path.isdir", return_value=True):
        result = check_datasets(cfg)

    assert result["check"] == "datasets"
    assert result["status"] == "pass"


def test_check_datasets_warn_when_env_missing(monkeypatch: MonkeyPatch) -> None:
    cfg = _base_config()
    monkeypatch.delenv("HB_DATASET", raising=False)
    monkeypatch.delenv("HS_DATASET", raising=False)

    result = check_datasets(cfg)
    assert result["check"] == "datasets"
    assert result["status"] == "warn"
    assert "HB_DATASET" in result["message"]
    assert "HS_DATASET" in result["message"]


def test_check_datasets_fail_when_env_set_but_path_missing(monkeypatch: MonkeyPatch) -> None:
    cfg = _base_config()
    monkeypatch.setenv("HB_DATASET", "/datasets/HBDataset")
    monkeypatch.setenv("HS_DATASET", "/datasets/HSDataset")

    with patch("os.path.exists", return_value=False):
        result = check_datasets(cfg)

    assert result["check"] == "datasets"
    assert result["status"] == "fail"
    assert "path does not exist" in result["message"]


def test_check_datasets_fail_when_env_path_not_directory(monkeypatch: MonkeyPatch) -> None:
    cfg = _base_config()
    monkeypatch.setenv("HB_DATASET", "/datasets/HBDataset")
    monkeypatch.setenv("HS_DATASET", "/datasets/HSDataset")

    with patch("os.path.exists", return_value=True), patch("os.path.isdir", return_value=False):
        result = check_datasets(cfg)

    assert result["check"] == "datasets"
    assert result["status"] == "fail"
    assert "not a directory" in result["message"]


def test_check_datasets_missing_config_keys_returns_fail() -> None:
    result = check_datasets({"paths": {"repo_root": "/repo"}})
    assert result["check"] == "datasets"
    assert result["status"] == "fail"
    assert "Missing config key" in result["message"]


def test_check_2hdmc_lib_pass_when_any_library_exists() -> None:
    cfg = _base_config()

    def exists_side_effect(path: object) -> bool:
        return path == "/repo/2hdmc/lib/lib2HDMC.so"

    with patch("os.path.exists", side_effect=exists_side_effect):
        result = check_2hdmc_lib(cfg)

    assert result["check"] == "2hdmc_lib"
    assert result["status"] == "pass"
    assert "lib2HDMC.so" in result["message"]


def test_check_2hdmc_lib_warn_when_missing() -> None:
    cfg = _base_config()
    with patch("os.path.exists", return_value=False):
        result = check_2hdmc_lib(cfg)

    assert result["check"] == "2hdmc_lib"
    assert result["status"] == "warn"
    assert "not found" in result["message"]


def test_check_2hdmc_lib_missing_repo_root_returns_fail() -> None:
    result = check_2hdmc_lib({"dihiggs": {}})
    assert result["check"] == "2hdmc_lib"
    assert result["status"] == "fail"
    assert "Missing config key" in result["message"]


def test_run_all_preflight_checks_aggregates_with_fail_priority(monkeypatch: MonkeyPatch) -> None:
    cfg = _base_config()
    monkeypatch.delenv("HB_DATASET", raising=False)
    monkeypatch.delenv("HS_DATASET", raising=False)

    with patch("os.path.exists", return_value=False):
        result = run_all_preflight_checks(cfg)

    assert result["overall"] == "fail"
    checks = result["checks"]
    assert isinstance(checks, list)
    assert [check["check"] for check in checks] == ["phys_exec", "datasets", "2hdmc_lib"]


def test_run_all_preflight_checks_aggregates_warn_when_no_fail(monkeypatch: MonkeyPatch) -> None:
    cfg = _base_config()
    monkeypatch.delenv("HB_DATASET", raising=False)
    monkeypatch.delenv("HS_DATASET", raising=False)

    with patch("os.path.exists", return_value=True), patch("os.access", return_value=True):
        result = run_all_preflight_checks(cfg)

    assert result["overall"] == "warn"


def test_run_all_preflight_checks_pass_when_all_checks_pass(monkeypatch: MonkeyPatch) -> None:
    cfg = _base_config()
    monkeypatch.setenv("HB_DATASET", "/datasets/HBDataset")
    monkeypatch.setenv("HS_DATASET", "/datasets/HSDataset")

    def exists_side_effect(path: object) -> bool:
        if not isinstance(path, str):
            return False
        return path in {
            "/repo/dihiggs/app/PhysScanWithFixings",
            "/datasets/HBDataset",
            "/datasets/HSDataset",
            "/repo/2hdmc/lib/lib2HDMC.a",
        }

    def isdir_side_effect(path: object) -> bool:
        return isinstance(path, str) and path in {"/datasets/HBDataset", "/datasets/HSDataset"}

    with (
        patch("os.path.exists", side_effect=exists_side_effect),
        patch("os.path.isdir", side_effect=isdir_side_effect),
        patch("os.access", return_value=True),
    ):
        result = run_all_preflight_checks(cfg)

    assert result["overall"] == "pass"
    statuses = [entry["status"] for entry in result["checks"]]
    assert statuses == ["pass", "pass", "pass"]
