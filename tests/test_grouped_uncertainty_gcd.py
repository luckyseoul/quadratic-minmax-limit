from __future__ import annotations

import pytest

from e1_gmin_m4_grouped_uncertainty_gcd import (
    even_support_counterexample_constraints,
    homogeneous_gcd_reduction,
    odd_support_grouped_uncertainty_theorem,
    projection_partition_factor_orders,
)


def test_local_odd_factors_detect_exactly_the_silent_partitions():
    silent = projection_partition_factor_orders(8, 0, (2, 2, 4))
    assert silent["factor_orders_by_class"] == (1, 1, 1, 1, 3, 3, 3, 3)
    assert silent["silent_direction"] is True
    assert silent["common_odd_factor"] is True

    radial = projection_partition_factor_orders(8, 2, (2, 4))
    assert radial["radial_factor_order"] == 3
    assert radial["silent_direction"] is True
    assert radial["common_odd_factor"] is True

    nonsilent = projection_partition_factor_orders(8, 2, (3, 3))
    assert nonsilent["factor_orders_by_class"] == (
        3,
        3,
        2,
        2,
        2,
        2,
        2,
        2,
    )
    assert nonsilent["silent_direction"] is False
    assert nonsilent["common_odd_factor"] is False


def test_homogeneous_gcd_formula_has_consistent_degrees():
    out = homogeneous_gcd_reduction(31, 8)
    assert out["proved"] is True
    assert out["detecting_form_degree"] == 15
    assert out["product_degree"] == 120
    assert out["exact_silent_gcd_formula"] == "z=deg(gcd_i odd(F_i))"
    assert out["odd_support_bound_proved"] is False
    assert out["even_support_bound_proved"] is False
    assert out["residual_ii_closed"] is False


def test_odd_support_branch_is_proved_but_even_branch_remains_open():
    for s in (1, 3, 5, 7):
        out = odd_support_grouped_uncertainty_theorem(31, s)
        assert out["proved"] is True
        assert out["conclusion"] == "z<=s"
    with pytest.raises(ValueError):
        odd_support_grouped_uncertainty_theorem(31, 8)


def test_s8_z9_counterexample_forces_eight_dependent_sextics():
    out = even_support_counterexample_constraints(31, 8, 9)
    assert out["proved"] is True
    assert out["necessary_common_squarefree_factor_degree"] == 9
    assert out["detecting_form_degree"] == 15
    assert out["quotient_form_degree_at_most"] == 6
    assert out["quotient_form_space_dimension_at_most"] == 7
    assert out["number_of_quotient_forms"] == 8
    assert out["quotient_forms_forced_linearly_dependent"] is True
    assert out["counterexample_excluded"] is False
    assert out["grouped_uncertainty_even_branch"] == "OPEN"
    assert out["residual_ii_closed"] is False


def test_invalid_partition_and_branch_inputs_fail_closed():
    with pytest.raises(ValueError):
        projection_partition_factor_orders(8, 2, (2, 2))
    with pytest.raises(ValueError):
        homogeneous_gcd_reduction(9, 8)
    with pytest.raises(ValueError):
        even_support_counterexample_constraints(31, 7, 9)
    with pytest.raises(ValueError):
        even_support_counterexample_constraints(31, 8, 8)
