from fractions import Fraction

import pytest

from e1_gmin_m4_mobius_global_incidence_barrier import (
    global_mobius_incidence_barrier,
    theorem_record,
)


@pytest.mark.parametrize("p", (31, 43, 47, 59, 71, 83, 103))
def test_all_prime_global_incidence_barrier(p):
    result = global_mobius_incidence_barrier(p)
    m = (p + 1) // 2
    assert result["proved"]
    assert result["physical_supports_pairwise_disjoint"]
    assert result["sum_of_halves_is_ternary"]
    assert result["sigma_equals_S_for_the_witness"]
    assert result["sigma_lower_bound_in_witness"] == 2 * m
    assert not result["standalone_global_incidence_bound_can_contradict_demand"]
    assert not result["endpoint_target_constructed"]
    assert not result["residual_ii_closed"]
    assert Fraction(result["expected_total_pair_overlap"]) < Fraction(3, 4)
    assert int(result["physical_pair_overlap_candidate_count"]) == (
        6 * p * p - 24 * p + 26
    )


def test_least_order_uses_exact_boundary_probability_gap():
    result = global_mobius_incidence_barrier(31)
    assert result["m"] == 16
    assert result["expected_total_pair_overlap"] == "10096/14415"
    assert Fraction(result["mcdiarmid_exponent"]) > Fraction(4, 3)
    assert "m=16" in result["probability_proof"]


def test_input_gate_and_record_scope():
    for bad in (29, 35, 41, True):
        with pytest.raises(ValueError):
            global_mobius_incidence_barrier(bad)
    record = theorem_record()
    assert record["proved_all_claimed_statements"]
    assert not record["residual_ii_closed"]
