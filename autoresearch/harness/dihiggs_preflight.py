from __future__ import annotations

import os
from collections.abc import Mapping
from typing import cast


def _result(check: str, status: str, message: str) -> dict[str, str]:
    return {"check": check, "status": status, "message": message}


def _as_mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
    return cast(Mapping[str, object], value)


def check_phys_exec(config: dict[str, object]) -> dict[str, str]:
    try:
        dihiggs_cfg = _as_mapping(config["dihiggs"], "dihiggs")
        phys_exec_raw = dihiggs_cfg["phys_exec"]
        if not isinstance(phys_exec_raw, str) or not phys_exec_raw:
            return _result("phys_exec", "fail", "dihiggs.phys_exec must be a non-empty string")

        # Interpolate {repo_root} placeholder if present
        phys_exec = phys_exec_raw
        if "{repo_root}" in phys_exec:
            paths_cfg = _as_mapping(config.get("paths", {}), "paths")
            repo_root = paths_cfg.get("repo_root", "")
            if not isinstance(repo_root, str):
                return _result("phys_exec", "fail", "paths.repo_root must be a string for interpolation")
            phys_exec = phys_exec.format(repo_root=repo_root)
        
        if not os.path.exists(phys_exec):
            return _result("phys_exec", "fail", f"phys_exec not found at {phys_exec}")
        if not os.access(phys_exec, os.X_OK):
            return _result("phys_exec", "fail", f"phys_exec is not executable: {phys_exec}")
        return _result("phys_exec", "pass", f"phys_exec found and executable at {phys_exec}")
    except KeyError as exc:
        return _result("phys_exec", "fail", f"Missing config key: {exc.args[0]}")
    except Exception as exc:
        return _result("phys_exec", "fail", f"Unexpected error while checking phys_exec: {exc}")


def check_datasets(config: dict[str, object]) -> dict[str, str]:
    try:
        dihiggs_cfg = _as_mapping(config["dihiggs"], "dihiggs")
        hb_env_raw = dihiggs_cfg["hb_dataset_env"]
        hs_env_raw = dihiggs_cfg["hs_dataset_env"]

        if not isinstance(hb_env_raw, str) or not hb_env_raw:
            return _result("datasets", "fail", "dihiggs.hb_dataset_env must be a non-empty string")
        if not isinstance(hs_env_raw, str) or not hs_env_raw:
            return _result("datasets", "fail", "dihiggs.hs_dataset_env must be a non-empty string")

        hb_env = hb_env_raw
        hs_env = hs_env_raw
        hb_path = os.environ.get(hb_env)
        hs_path = os.environ.get(hs_env)

        missing_envs = [name for name, value in ((hb_env, hb_path), (hs_env, hs_path)) if not value]
        if missing_envs:
            missing_joined = ", ".join(missing_envs)
            return _result(
                "datasets",
                "warn",
                f"Dataset env var(s) not set: {missing_joined}; continuing preflight",
            )

        invalid_paths: list[str] = []
        assert hb_path is not None
        assert hs_path is not None
        for name, path in ((hb_env, hb_path), (hs_env, hs_path)):
            if not os.path.exists(path):
                invalid_paths.append(f"{name}={path} (path does not exist)")
            elif not os.path.isdir(path):
                invalid_paths.append(f"{name}={path} (not a directory)")

        if invalid_paths:
            return _result("datasets", "fail", "Invalid dataset path(s): " + "; ".join(invalid_paths))

        return _result(
            "datasets",
            "pass",
            f"Dataset env vars configured and directories exist: {hb_env}={hb_path}, {hs_env}={hs_path}",
        )
    except KeyError as exc:
        return _result("datasets", "fail", f"Missing config key: {exc.args[0]}")
    except Exception as exc:
        return _result("datasets", "fail", f"Unexpected error while checking datasets: {exc}")


def check_2hdmc_lib(config: dict[str, object]) -> dict[str, str]:
    try:
        paths_cfg = _as_mapping(config["paths"], "paths")
        repo_root_raw = paths_cfg["repo_root"]
        if not isinstance(repo_root_raw, str) or not repo_root_raw:
            return _result("2hdmc_lib", "fail", "paths.repo_root must be a non-empty string")

        lib_dir = os.path.join(repo_root_raw, "2hdmc", "lib")
        candidates = (
            os.path.join(lib_dir, "lib2HDMC.so"),
            os.path.join(lib_dir, "lib2HDMC.a"),
            os.path.join(lib_dir, "lib2HDMC.dylib"),
        )
        for candidate in candidates:
            if os.path.exists(candidate):
                return _result("2hdmc_lib", "pass", f"2HDMC library found at {candidate}")

        return _result(
            "2hdmc_lib",
            "warn",
            f"2HDMC library not found in {lib_dir} (checked lib2HDMC.so/.a/.dylib)",
        )
    except KeyError as exc:
        return _result("2hdmc_lib", "fail", f"Missing config key: {exc.args[0]}")
    except Exception as exc:
        return _result("2hdmc_lib", "fail", f"Unexpected error while checking 2hdmc lib: {exc}")


def run_all_preflight_checks(config: dict[str, object]) -> dict[str, object]:
    checks = [check_phys_exec(config), check_datasets(config), check_2hdmc_lib(config)]
    statuses = {check["status"] for check in checks}

    if "fail" in statuses:
        overall = "fail"
    elif "warn" in statuses:
        overall = "warn"
    else:
        overall = "pass"

    return {"checks": checks, "overall": overall}
