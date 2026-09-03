from __future__ import annotations

import pytest

from e1_gmin_m4_mobius_half_intersections import (
    branch_c_pairwise_cancellation_bound,
    known_one_origin_witness,
    sharp_two_cancellation_witness,
    theorem_record,
    two_cancellation_locus_theorem,
    two_half_intersection_candidates,
)
from e1_gmin_m4_inversion_antisymmetric_radon import _negative_edge
from e1_gmin_m4_mobius_half_symmetric import (
    mobius_parameter_edges,
)


def test_known_common_origin_choice_has_exactly_one_cancellation():
    out = known_one_origin_witness(31)
    assert out["cancelled_inversion_orbits"] == 1
    assert out["parameters"] == [0, 0]
    assert out["two_trade_sum_is_ternary"] is True
    assert out["proved"] is True


def test_uniform_pair_bound_is_sharp_at_two_opposite_orbits():
    out = sharp_two_cancellation_witness(31)
    assert out["normalized_common_auxiliary"] == "(2/3)*(X+Y)"
    assert out["opposite_orientation_parameters"] == [[0, 0], [29, 29]]
    assert out["cancelled_inversion_orbits"] == 2
    assert out["same_orientation_shared_orbits"] == 0
    assert out["two_trade_sum_is_ternary"] is True
    assert out["pairwise_bound_is_sharp"] is True
    assert out["proved"] is True


def test_p31_sharp_witness_matches_the_canonical_mobius_halves():
    p = 31
    two_thirds = 2 * pow(3, -1, p) % p
    common_auxiliary = (two_thirds, two_thirds)
    first = set(
        mobius_parameter_edges(
            p, (1, 0), common_auxiliary, center=1
        ).values()
    )
    second = set(
        mobius_parameter_edges(
            p, (0, 1), common_auxiliary, center=1
        ).values()
    )
    negative_second = {_negative_edge(p, edge) for edge in second}
    assert len(first) == len(second) == p - 1
    assert len(first & second) == 0
    assert len(first & negative_second) == 2


def test_candidate_classifier_uses_only_four_forced_matchings():
    inverse_two = pow(2, -1, 31)
    out = two_half_intersection_candidates(
        31,
        q=inverse_two,
        r=inverse_two,
        A=3 * inverse_two,
        B=3 * inverse_two,
    )
    assert len(out["attempted_candidates"]) == 4
    assert out["shared_inversion_orbits"] == 2
    assert out["same_orientation_shared_orbits"] == 0
    assert out["opposite_orientation_shared_orbits"] == 2
    assert out["uniform_total_bound"] == 4
    assert out["uniform_opposite_orientation_bound"] == 2
    assert out["proved"] is True


def test_two_cancellation_locus_is_rigid_not_a_free_family():
    out = two_cancellation_locus_theorem(31)
    assert out["combined_factorization"] == (
        "(2*q-1)*(q+1)=0 with q=r"
    )
    assert out["unique_normalized_parameters"] == {
        "q": 16,
        "r": 16,
        "A": 17,
        "B": 17,
    }
    assert out["same_orientation_matches_at_unique_point"] == 0
    assert out["free_parameter_after_two_cancellations"] is False
    assert out["greedy_pairing_from_free_locus_available"] is False
    assert out["proved"] is True


def test_p31_branch_c_pairwise_bound_does_not_close_the_ray():
    lower = branch_c_pairwise_cancellation_bound(31, 68)
    assert lower["raw_selected_orbit_occurrences"] == 480
    assert lower["target_edge_count"] == 261
    assert lower["required_cancellation_units"] == 110
    assert lower["trade_pair_count"] == 120
    assert lower["pairwise_cancellation_upper_bound"] == 240
    assert lower["pairwise_bound_rules_out_required_cancellation"] is False
    assert lower["forced_fixed_weight_must_be_odd_if_feasible"] is True
    assert lower["fixed_edge_objective_evaluated_by_intersection_theorem"] is False

    upper = branch_c_pairwise_cancellation_bound(31, 177)
    assert upper["target_edge_count"] == 479
    assert upper["required_cancellation_units"] == 1
    assert upper["pairwise_cancellation_upper_bound"] == 240
    assert upper["pairwise_bound_rules_out_required_cancellation"] is False


def test_theorem_record_stops_before_fixed_edge_and_boolean_feasibility():
    out = theorem_record(31)
    assert out["proved"] is True
    assert out["fixed_edge_objective_closed"] is False
    assert out["residual_ii_closed"] is False
    assert out["status"] == (
        "PAIRWISE INTERSECTION THEOREM PROVED; "
        "FIXED-EDGE OBJECTIVE AND SYMMETRIC FIBRE OPEN"
    )


def test_parameter_guards():
    with pytest.raises(ValueError):
        two_half_intersection_candidates(9, 0, 0, 1, 1)
    with pytest.raises(ValueError):
        two_half_intersection_candidates(31, 0, 0, 0, 1)
    with pytest.raises(ValueError):
        branch_c_pairwise_cancellation_bound(31, 67)
