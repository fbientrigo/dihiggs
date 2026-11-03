from dihiggs.app.orchestrator.aligned_constraints import (
    bfb_and_global_minimum, unitarity,
)


def test_general_scattering_blocks_keep_all_twelve_eigenvalues():
    lam = {f"lambda{i}": "0" for i in range(1,8)}
    lam["lambda1"] = "1"
    result = unitarity(lam, precision=50)
    assert result["eigenvalue_count"] == 12
    assert result["status"] == "PASS"
    assert result["tolerance"]["a0"] == "Lambda/(16*pi)"


def test_exact_charge_breaking_negative_witness():
    z = {f"Z{i}": "0" for i in range(1,8)}
    z.update(Z1="1", Z2="1", Z4="-2")
    bfb, global_min = bfb_and_global_minimum(z, {})
    assert bfb["status"] == "FAIL"
    assert global_min["status"] == "FAIL"
    assert bfb["evidence"]["x"] != "0" or bfb["evidence"]["y"] != "0"


def test_exact_charge_interior_negative_witness():
    z = {f"Z{i}": "0" for i in range(1,8)}
    z.update(Z1="1", Z2="1", Z3="-2")
    bfb, _ = bfb_and_global_minimum(z, {})
    assert bfb["status"] == "FAIL"
    assert bfb["evidence"]["u"] != "0" and bfb["evidence"]["b"] != "0"
    assert bfb["evidence"]["x"] == "0" and bfb["evidence"]["y"] == "0"


def test_exact_complex_negative_witness():
    z = {f"Z{i}": "0" for i in range(1,8)}
    z.update(Z1="1", Z2="1", Z5="2")
    bfb, _ = bfb_and_global_minimum(z, {})
    assert bfb["status"] == "FAIL"
    assert bfb["evidence"]["y"] != "0"


def test_marginal_degenerate_case_remains_inconclusive():
    z = {f"Z{i}": "0" for i in range(1,8)}
    bfb, global_min = bfb_and_global_minimum(z, {})
    assert bfb["status"] == "INCONCLUSIVE"
    assert global_min["status"] == "INCONCLUSIVE"


def test_unitarity_is_strict_and_near_boundary_is_not_forced():
    import mpmath as mp
    lam = {f"lambda{i}": "0" for i in range(1,8)}
    with mp.workdps(80):
        lam["lambda1"] = mp.nstr(8*mp.pi/3, 80)
    result = unitarity(lam, precision=80)
    assert result["status"] == "INCONCLUSIVE"
