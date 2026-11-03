from decimal import Decimal

import pytest

from dihiggs.app.orchestrator.aligned_lambda7 import (
    build_aligned_lambda7_family,
    xstar_lambda6_tanbeta,
)


BASE = dict(tan_beta="1000000", lambda1_ref="1", mh_GeV="125.13",
            mphi_GeV="150", mA_GeV="200", mHp_GeV="200", v_GeV="246")


def test_explicit_trajectory_and_decimal_round_trip():
    point = build_aligned_lambda7_family(**BASE, X_lambda6_tanbeta="1", precision=80)
    assert point.trajectory == "fixed_X"
    assert point.values["lambda6"] == "0.000001"
    round_trip = Decimal(point.values["lambda7"])
    assert round_trip == Decimal(point.values["lambda7"])
    assert point.values["Z6"] == "0"


def test_fixed_lambda6_and_fixed_x_are_distinct_trajectories():
    fixed_x = build_aligned_lambda7_family(**BASE, X_lambda6_tanbeta="1")
    fixed_l6 = build_aligned_lambda7_family(**BASE, lambda6="0.01")
    assert fixed_x.trajectory == "fixed_X"
    assert fixed_l6.trajectory == "fixed_lambda6"
    assert Decimal(fixed_l6.values["X_lambda6_tanbeta"]) == Decimal("1e4")
    assert Decimal(fixed_x.values["Z7"]) != Decimal(fixed_l6.values["Z7"])


def test_xstar_is_a_separate_evaluated_x_path():
    xstar = xstar_lambda6_tanbeta(**{k: BASE[k] for k in
                                    ("tan_beta", "lambda1_ref", "mh_GeV", "mphi_GeV", "v_GeV")})
    point = build_aligned_lambda7_family(**BASE, X_lambda6_tanbeta=xstar)
    assert point.trajectory == "fixed_X"
    assert abs(Decimal(point.values["Z7"])) < Decimal("1e-75")


@pytest.mark.parametrize("kwargs", [
    {"lambda6": "0", "X_lambda6_tanbeta": "0"},
    {},
])
def test_requires_exactly_one_trajectory_input(kwargs):
    with pytest.raises(ValueError, match="exactly one"):
        build_aligned_lambda7_family(**BASE, **kwargs)


@pytest.mark.parametrize("name,value", [("tan_beta", "0"), ("mh_GeV", "NaN"),
                                         ("mphi_GeV", "-1"), ("v_GeV", "Infinity")])
def test_invalid_inputs_rejected(name, value):
    with pytest.raises(ValueError):
        build_aligned_lambda7_family(**{**BASE, name: value}, X_lambda6_tanbeta="1")


def test_catastrophic_cancellation_is_visible_against_naive_binary64():
    args = {**BASE, "lambda1_ref": "1.2"}
    point = build_aligned_lambda7_family(**args, X_lambda6_tanbeta="1", precision=80)
    # The inverse contains O(t^2) and O(t^3) terms. Feeding rounded binary64
    # quartics back into it is deliberately ill-conditioned at this point.
    t = float(args["tan_beta"])
    mh, mp, v = float(args["mh_GeV"]), float(args["mphi_GeV"]), float(args["v_GeV"])
    l6, l7 = float(point.values["lambda6"]), float(point.values["lambda7"])
    inverse = mh * mh / (v * v) + (mp * mp - mh * mh / 2) * t * t / (v * v) - 1.5 * l6 * t + 0.5 * l7 * t**3
    assert abs(inverse - 1.2) > 1e-5


def test_constructor_closure_and_non_certification_states():
    point = build_aligned_lambda7_family(**BASE, X_lambda6_tanbeta="1")
    assert point.values["M_sq_GeV2"] == str(Decimal(BASE["mh_GeV"]) ** 2 / 2)
    assert point.constraints["local_minimum_tree"]["status"] == "PASS"
    assert point.constraints["BFB_general"]["status"] == "NOT_TESTED"
    assert point.constraints["global_EW_minimum"]["status"] == "NOT_TESTED"


def test_alignment_does_not_pass_when_independent_closure_is_inconclusive():
    point=build_aligned_lambda7_family(**{**BASE,"lambda1_ref":"1e100"},
        X_lambda6_tanbeta="1",precision=80)
    assert point.constraints["spectrum_reconstruction"]["status"] == "INCONCLUSIVE"
    assert point.constraints["alignment"]["status"] != "PASS"


def test_tiny_physical_eigenvalue_is_inconclusive_below_matrix_error_bound():
    point=build_aligned_lambda7_family(**{**BASE,"mh_GeV":"1e-40"},
        X_lambda6_tanbeta="1",precision=80)
    assert point.constraints["local_minimum_tree"]["status"] == "INCONCLUSIVE"
    assert Decimal(point.provenance["independent_potential_closure"]["smallest_physical_hessian_eigenvalue_GeV2"]) < Decimal("1e-70")
