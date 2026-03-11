import csv
from pathlib import Path
from typing import Callable, cast


def test_per_tanbeta_budget_allocation_with_doe_floor():
    from dihiggs.app import adaptive_policy

    alloc_any = getattr(adaptive_policy, "allocate_budget_per_tanbeta", None)
    assert alloc_any is not None, "missing adaptive_policy.allocate_budget_per_tanbeta (v2 API)"
    alloc = cast(
        Callable[..., dict[str, dict[str, int]]],
        alloc_any,
    )

    bins = adaptive_policy.make_lam1_bins(0.0, 1.0, 2)
    bin_ids = [b["id"] for b in bins]
    b0, b1 = bin_ids

    successes_by_tb_by_bin = {
        "10": {b0: 8, b1: 0},
        "2000": {b0: 0, b1: 8},
    }
    trials_by_tb_by_bin = {
        "10": {b0: 10, b1: 10},
        "2000": {b0: 10, b1: 10},
    }

    total_budget_per_tb = 5
    floor_points = 2

    budgets1 = alloc(
        bins=bins,
        successes_by_tb_by_bin=successes_by_tb_by_bin,
        trials_by_tb_by_bin=trials_by_tb_by_bin,
        total_budget_per_tb=total_budget_per_tb,
        floor_points=floor_points,
        alpha=1.0,
        beta=1.0,
    )
    budgets2 = alloc(
        bins=bins,
        successes_by_tb_by_bin=successes_by_tb_by_bin,
        trials_by_tb_by_bin=trials_by_tb_by_bin,
        total_budget_per_tb=total_budget_per_tb,
        floor_points=floor_points,
        alpha=1.0,
        beta=1.0,
    )

    assert budgets1 == budgets2
    assert set(budgets1.keys()) == {"10", "2000"}

    for tb_tag in ["10", "2000"]:
        tb_budgets = budgets1[tb_tag]
        assert set(tb_budgets.keys()) == set(bin_ids)
        assert sum(tb_budgets.values()) == total_budget_per_tb
        assert all(tb_budgets[bid] >= floor_points for bid in bin_ids)

    assert budgets1["10"][b0] >= budgets1["10"][b1]
    assert budgets1["2000"][b1] >= budgets1["2000"][b0]


def test_skip_event_recovers_valids_from_csv_for_budgeting(tmp_path: Path):
    from dihiggs.app import adaptive_artifacts

    st_any = getattr(adaptive_artifacts, "successes_trials_from_task_record", None)
    assert (
        st_any is not None
    ), "missing adaptive_artifacts.successes_trials_from_task_record (v2 API)"
    successes_trials_from_task_record = cast(
        Callable[[dict[str, object]], tuple[int, int]],
        st_any,
    )

    out_csv = tmp_path / "task_output.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["positivity_ok", "unitarity_ok", "perturbativity_ok", "dummy"],
        )
        w.writeheader()
        w.writerow(
            {
                "positivity_ok": 1,
                "unitarity_ok": 1,
                "perturbativity_ok": 1,
                "dummy": "a",
            }
        )
        w.writerow(
            {
                "positivity_ok": 1,
                "unitarity_ok": 1,
                "perturbativity_ok": 1,
                "dummy": "b",
            }
        )
        w.writerow(
            {
                "positivity_ok": 1,
                "unitarity_ok": 1,
                "perturbativity_ok": 1,
                "dummy": "c",
            }
        )
        w.writerow(
            {
                "positivity_ok": 1,
                "unitarity_ok": 1,
                "perturbativity_ok": 0,
                "dummy": "bad1",
            }
        )
        w.writerow(
            {
                "positivity_ok": 0,
                "unitarity_ok": 1,
                "perturbativity_ok": 1,
                "dummy": "bad2",
            }
        )

    record: dict[str, object] = {
        "event": "skip",
        "output_csv": str(out_csv),
        "previous_csv": None,
        "n_lam1_effective": 2,
        "grid_signature": "mphi=130.0000-290.0000-N15|lam1=0.0-0.5-N2",
        "tanbeta": 10.0,
        "tb_tag": "10",
    }

    successes, trials = successes_trials_from_task_record(record)
    assert successes == 3
    assert trials == 15 * 2


def test_skip_event_prefers_scan_meta_over_csv(tmp_path: Path):
    from dihiggs.app import adaptive_artifacts

    out_csv = tmp_path / "task_output.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["positivity_ok", "unitarity_ok", "perturbativity_ok", "dummy"],
        )
        w.writeheader()
        for i in range(3):
            w.writerow(
                {
                    "positivity_ok": 1,
                    "unitarity_ok": 1,
                    "perturbativity_ok": 1,
                    "dummy": f"ok{i}",
                }
            )

    meta_path = tmp_path / "scan_meta.json"
    meta_path.write_text(
        '{"event":"done","triple_ok_points":"7"}\n',
        encoding="utf-8",
    )

    record: dict[str, object] = {
        "event": "skip",
        "output_csv": str(out_csv),
        "previous_csv": None,
        "n_lam1_effective": 2,
        "grid_signature": "mphi=130.0000-290.0000-N15|lam1=0.0-0.5-N2",
        "tanbeta": 10.0,
        "tb_tag": "10",
    }

    successes, trials = adaptive_artifacts.successes_trials_from_task_record(record)
    assert successes == 7
    assert trials == 15 * 2


def test_parse_run_dir_from_timestamped_orchestrator_output():
    from dihiggs.app import adaptive_artifacts

    parsed = adaptive_artifacts.parse_run_dir_from_orchestrator_output(
        "[2026-02-19 22:38:11] [PATH] run_dir = /tmp/foo"
    )

    assert parsed == Path("/tmp/foo")
