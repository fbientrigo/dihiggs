from fractions import Fraction

import sympy as sp

from dihiggs.app.orchestrator.physical_certification import (
    RealPotential,
    certify_bfb,
    certify_global_minimum,
    certify_marginal_full_stability,
    certify_runaway,
    certify_tree_unitarity,
    enumerate_stationary_candidates,
    rounded_input_potential,
)


def point(*, Y=(0, 0, 0), Z=(0, 0, 0, 0, 0, 0, 0)):
    return RealPotential(tuple(map(sp.Rational, Y)), tuple(map(sp.Rational, Z)))


def test_exact_bfb_witnesses_cover_neutral_real_complex_and_charge_directions():
    real = certify_bfb(point(Z=(1, 1, -2, 0, 0, 0, 0)))
    complex_ = certify_bfb(point(Z=(1, 1, 0, -2, 0, 0, 0)))
    charged = certify_bfb(point(Z=(1, 1, 0, -1, -1, 0, 0)))
    assert [real["status"], complex_["status"], charged["status"]] == ["FAIL"] * 3
    assert all(sp.Rational(result["evidence"]["quartic_value"]) < 0
               for result in (real, complex_, charged))


def test_slemma_proves_strict_and_marginal_nonnegative_quartics():
    strict = certify_bfb(point(Z=(2, 2, 0, 0, 0, 0, 0)))
    marginal = certify_bfb(point(Z=(0, 2, 0, 0, 0, 0, 0)))
    assert strict["status"] == marginal["status"] == "PASS"
    assert strict["evidence"]["strict_certificate"] is True
    assert marginal["evidence"]["strict_certificate"] is False


def test_slemma_finds_narrow_strict_certificate_near_small_lambda1():
    result = certify_bfb(point(Z=("1e-28", 1, 0, 0, 0, 0, 0)))
    assert result["status"] == "PASS"
    assert result["evidence"]["strict_certificate"] is True


def test_corrected_uniform_C_gate_for_stable_marginal_case():
    p = point(Y=(0, -1, 0), Z=(0, 2, 0, 0, 0, 0, 0))
    result = certify_marginal_full_stability(p)
    assert result["status"] == "PASS"
    assert result["evidence"]["uniform_C"] == "1"


def test_flat_quartic_with_negative_mixed_scale_ray_is_not_called_bounded():
    p = point(Y=(0, 0, -1), Z=(0, 2, 0, 0, 0, 0, 0))
    assert certify_bfb(p)["status"] == "PASS"  # quartic alone is nonnegative
    result = certify_runaway(p)
    assert result["status"] == "FAIL"
    assert result["evidence"]["V_of_n"] == "-n**2"
    assert certify_global_minimum(p, (1, 0, 0, 0))["status"] == "FAIL"


def test_negative_quartic_ray_proves_full_unboundedness_without_stationary_scan():
    p=point(Z=(-1,1,0,0,0,0,0))
    result=certify_global_minimum(p,(0,0,0,0))
    assert result["status"]=="FAIL"
    assert result["evidence"]["V_of_n"]=="-n**2/32"
    assert result["evidence"]["stationary_enumeration"]["status"]=="NOT_TESTED"


def test_general_neutral_real_polynomial_is_exact_and_degree_at_most_four():
    p = point(Y=(-3, 2, 1), Z=(2, 4, 3, 5, 7, 11, 13))
    result = enumerate_stationary_candidates(p)
    r = sp.symbols("r", real=True)
    P = -3 + 2*r + 2*r**2
    Q = 1 + 22*r + 15*r**2 + 26*r**3 + 2*r**4
    expected = sp.Poly(sp.expand(2*Q*sp.diff(P,r)-P*sp.diff(Q,r)),r)
    actual_expr=sp.sympify(result["neutral_real_critical_polynomial"],locals={"r":r})
    assert sp.Poly(actual_expr,r) == expected
    assert result["neutral_real_critical_degree"] <= 4
    assert result["status"] in {"PASS", "INCONCLUSIVE"}


def test_aligned_pilot_has_complete_exact_strata_and_global_certificate():
    from dihiggs.app.orchestrator.aligned_lambda7 import build_aligned_lambda7_family
    family=build_aligned_lambda7_family(tan_beta="10",lambda1_ref="1",mh_GeV="125.13",
        mphi_GeV="200",mA_GeV="250",mHp_GeV="250",v_GeV="246",lambda6="1e-10")
    p=RealPotential.from_mapping(family.values)
    enum=enumerate_stationary_candidates(p)
    target=(-Fraction(p.Y[0])/Fraction(p.Z[0]),0,0,0)
    global_result=certify_global_minimum(p,target)
    assert enum["status"] == "PASS" and enum["complete"] is True
    assert global_result["status"] == "PASS"
    assert "decomposition" in global_result["evidence"]


def test_complex_fixed_multiplier_stationary_witness_is_retained():
    p = point(Y=(-3, -3, 0), Z=(2, 2, 0, 2, 1, 0, 0))
    result = enumerate_stationary_candidates(p)
    complex_points = [c for c in result["candidates"] if c["stratum"].startswith("neutral-complex")]
    assert complex_points
    assert sp.simplify(sp.sympify(complex_points[0]["u"])*sp.sympify(complex_points[0]["b"])
                       - sp.sympify(complex_points[0]["x"])**2
                       - sp.sympify(complex_points[0]["y"])**2) == 0


def test_degenerate_stationary_continuum_stays_inconclusive_and_finds_lower_point():
    # V=(u-1)^2+(b-1)^2 up to an irrelevant constant; the interior stationary
    # set has x^2+y^2<=u*b and energy -2, below the EW target energy -1.
    p = point(Y=(-2, -2, 0), Z=(2, 2, 0, 0, 0, 0, 0))
    enum = enumerate_stationary_candidates(p)
    assert enum["status"] == "INCONCLUSIVE"
    assert any("singular affine family" in reason for reason in enum["unresolved_reasons"])
    result = certify_global_minimum(p, (1, 0, 0, 0))
    assert result["status"] == "FAIL"
    assert result["evidence"]["lower_points"]
    assert "stationary_enumeration" in result["evidence"]


def test_exact_general_tree_unitarity_threshold_and_degenerate_boundary():
    below = {f"lambda{i}": "0" for i in range(1,8)}
    below["lambda1"] = "8.37"
    above = dict(below, lambda1="8.38")
    assert certify_tree_unitarity(below)["status"] == "PASS"
    assert certify_tree_unitarity(above)["status"] == "FAIL"
    degenerate = dict(below, lambda1="0")
    assert certify_tree_unitarity(degenerate)["status"] == "PASS"
    asymmetric={"lambda1":"1","lambda2":"2","lambda3":"0.3","lambda4":"0.4",
                "lambda5":"0.5","lambda6":"0.7","lambda7":"-0.8"}
    exact=certify_tree_unitarity(asymmetric)
    from dihiggs.app.orchestrator.aligned_constraints import unitarity
    assert exact["status"] == unitarity(asymmetric)["status"]


def test_rounded_input_model_uses_sm_getter_v2_and_exact_tadpole_identities():
    inputs = {f"lambda{i}": "0.5" for i in range(1,8)}
    inputs.update({"lambda4":"0", "lambda5":"0", "lambda6":"0", "lambda7":"0",
                   "m12_sq_GeV2": "1000.00000000001", "tan_beta": "3.2", "v_GeV": "246"})
    result = rounded_input_potential(inputs)
    assert result["inputs_binary64_rational"]["SM_GF_binary64"] is not None
    assert result["stationary"] is True
    assert result["status"] == "PASS"
    assert result["potential"]["v2_GeV2"] != str(sp.Rational(246)**2)


def test_nonaligned_empty_doublet_axis_is_enumerated_with_field_tadpole_only():
    p = point(Y=(-5, -1, 0), Z=(0, 2, 7, 0, 0, 0, 0))
    result = enumerate_stationary_candidates(p)
    assert any(c["stratum"] == "u=0 axis" and c["b"] == "1/2"
               and c["V_GeV4"] == "-1/4" for c in result["candidates"])


def test_nonstationary_global_target_fails_exactly():
    p = point(Y=(0, 0, 0), Z=(2, 2, 0, 0, 0, 0, 0))
    result = certify_global_minimum(p, (1, 0, 0, 0))
    assert result["status"] == "FAIL"
    assert "tadpoles" in result["evidence"]["reason"]
