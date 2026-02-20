def test_make_lam1_bins_deterministic_and_cover_range():
    from dihiggs.app.adaptive_policy import make_lam1_bins

    bins1 = make_lam1_bins(lam1_min=0.0, lam1_max=3.0, n_bins=3)
    bins2 = make_lam1_bins(lam1_min=0.0, lam1_max=3.0, n_bins=3)

    assert bins1 == bins2
    assert [b["id"] for b in bins1] == ["lam1_bin_000", "lam1_bin_001", "lam1_bin_002"]
    assert bins1[0]["lam1_lo"] == 0.0
    assert bins1[-1]["lam1_hi"] == 3.0
    assert bins1[0]["lam1_hi"] == bins1[1]["lam1_lo"]
    assert bins1[1]["lam1_hi"] == bins1[2]["lam1_lo"]


def test_beta_binomial_posterior_mean_matches_closed_form():
    from dihiggs.app.adaptive_policy import beta_binomial_posterior_mean

    assert beta_binomial_posterior_mean(successes=0, trials=0, alpha=1.0, beta=1.0) == 0.5
    assert beta_binomial_posterior_mean(successes=1, trials=1, alpha=1.0, beta=1.0) == 2.0 / 3.0
    assert beta_binomial_posterior_mean(successes=8, trials=10, alpha=1.0, beta=1.0) == 9.0 / 12.0


def test_budget_allocation_from_fake_task_summary():
    from dihiggs.app.adaptive_policy import allocate_budget, make_lam1_bins

    bins = make_lam1_bins(lam1_min=0.0, lam1_max=3.0, n_bins=3)
    bin_ids = [b["id"] for b in bins]
    hot_id, cold_id, missing_id = bin_ids

    successes_by_bin = {
        hot_id: 8,
        cold_id: 0,
    }
    trials_by_bin = {
        hot_id: 10,
        cold_id: 10,
    }

    total_budget = 17
    floor_points = 3

    budgets1 = allocate_budget(
        bins=bins,
        successes_by_bin=successes_by_bin,
        trials_by_bin=trials_by_bin,
        total_budget=total_budget,
        floor_points=floor_points,
        alpha=1.0,
        beta=1.0,
    )
    budgets2 = allocate_budget(
        bins=bins,
        successes_by_bin=successes_by_bin,
        trials_by_bin=trials_by_bin,
        total_budget=total_budget,
        floor_points=floor_points,
        alpha=1.0,
        beta=1.0,
    )

    assert budgets1 == budgets2
    assert set(budgets1.keys()) == set(bin_ids)
    assert sum(budgets1.values()) == total_budget
    assert all(budgets1[bid] >= floor_points for bid in bin_ids)
    assert budgets1[hot_id] >= budgets1[cold_id]
    assert budgets1[missing_id] >= floor_points
