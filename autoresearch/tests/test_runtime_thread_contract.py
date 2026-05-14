from __future__ import annotations

from autoresearch.harness.campaign_supervisor import CampaignSupervisor


def test_autoscaling_does_not_mutate_cpp_omp_threads_without_explicit_code_path() -> None:
    supervisor = CampaignSupervisor.__new__(CampaignSupervisor)
    supervisor.config = {
        "runtime": {
            "threads": 4,
            "timeout_sec": 60,
            "cpp_omp_threads": 1,
        },
        "limits": {
            "max_new_run_dirs_per_round": 1,
            "max_new_run_dirs_per_arm_call": 1,
        },
    }

    supervisor._apply_autoscaling(
        {
            "action": "scale_up",
            "next_threads": 8,
            "next_timeout_sec": 120,
            "next_max_new_run_dirs_per_round": 2,
            "next_max_new_run_dirs_per_arm_call": 2,
        }
    )

    assert supervisor.config["runtime"]["threads"] == 8
    assert supervisor.config["runtime"]["cpp_omp_threads"] == 1
