from __future__ import annotations

import pytest

from e1_gmin_m4_mobius_half_symmetric import mobius_parameter_edges
from e1_gmin_m4_prescribed_center_common_block import (
    alpha_equal_beta_anchor_hall_theorem,
    anchor_graph_pseudoforest_profile,
    dependent_anchor_line_cover_theorem,
    p31_prescribed_center_anchor_obstruction,
    theorem_record,
)


def test_equal_square_intrinsic_criterion_is_an_anchor_slot_problem():
    out = alpha_equal_beta_anchor_hall_theorem(31)
    assert out["proved"] is True
    assert out["half_count"] == 16
    assert out["eta_plus_reduction"] == "L(a_2)=-j/2"
    assert out["eta_minus_reduction"] == "L(a_1)=+j/2"
    assert out["criterion_is_necessary_and_sufficient"] is True
    assert out["anchor_necessary_condition"] == (
        "|union_(i in P) A_i| >= |P|-1 for every P subset I"
    )


def test_p31_prescribed_nonzero_hard_centers_have_deficiency_two():
    out = p31_prescribed_center_anchor_obstruction()
    assert out["proved"] is True
    assert out["common_direction_K"] == [1, 4]
    assert out["K_is_Paley_opposite"] is True
    assert out["hard_direction_count"] == 16
    assert out["all_prescribed_hard_centers_nonzero"] is True
    assert len(out["prescribed_centers"]) == 16
    assert out["four_anchor_q_values"] == [0, 2, 11, 14]
    assert out["six_distinct_zero_q_values"] == [1, 7, 8, 21, 22, 28]
    assert out["anchor_family_size"] == 6
    assert out["anchor_union_size"] == 4
    assert out["anchor_deficiency"] == 2
    assert out["fixed_K_equal_square_saturated_cover_impossible"] is True
    assert out["another_K_or_unequal_square_cover_excluded"] is False
    assert out["full_dual_support_excluded"] is False
    assert out["residual_ii_closed"] is False


def test_anchor_SDR_is_exactly_the_pseudoforest_condition():
    tree_and_cycle = anchor_graph_pseudoforest_profile(
        [(0, 1), (1, 2), (3, 4), (4, 5), (5, 3)]
    )
    assert tree_and_cycle["proved_graph_equivalence"] is True
    assert tree_and_cycle["pseudoforest"] is True
    assert tree_and_cycle["edge_to_incident_vertex_SDR_exists"] is True
    assert sorted(
        component["cycle_rank"]
        for component in tree_and_cycle["components"]
    ) == [0, 1]

    bicycle = anchor_graph_pseudoforest_profile(
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]
    )
    assert bicycle["pseudoforest"] is False
    assert bicycle["edge_to_incident_vertex_SDR_exists"] is False
    assert bicycle["components"] == [
        {
            "vertex_count": 4,
            "edge_count": 5,
            "cycle_rank": 2,
            "at_most_unicyclic": False,
        }
    ]


def test_good_dependent_anchor_line_gives_a_one_way_cover_theorem():
    out = dependent_anchor_line_cover_theorem(31)
    assert out["proved"] is True
    assert out["h"] == 15
    assert out["hard_half_count"] == 16
    assert out["dependent_half_is_single"] is True
    assert out["dependent_half_unique_common_parameter"] is True
    assert out["other_half_count"] == 15
    assert out["zero_points_distinct"] is True
    assert out["forbidden_triples_are_distinct_APs"] is True
    assert out["anchor_graph_simple"] is True
    assert out["anchor_SDR_iff_pseudoforest"] is True
    assert out["remaining_point_count_after_anchor_SDR"] == 16
    assert out["free_slot_matching_exists"] is True
    assert out["cover_profile"] == (
        "h doubled halves and one dependent singleton"
    )
    assert out["covered_midpoint_class_count"] == 31
    assert out["some_good_base_line_proved"] is False
    assert out["converse_claimed"] is False
    assert out["mutual_ternarity_proved"] is False
    assert out["full_dual_support_containment_proved"] is False
    assert out["residual_ii_closed"] is False


def test_dependent_half_supplies_an_arbitrary_unique_singleton():
    p = 31
    inverse_two = pow(2, -1, p)
    direction = (1, 0)
    center = 2
    desired_midpoint = (1, 3)
    auxiliary = (28, 1)  # M(1,3)=0
    edges = mobius_parameter_edges(p, direction, auxiliary, center)
    assert set(edges[0]) == {(0, 0), (2, 6)}

    common_parameters = []
    for parameter, edge in edges.items():
        midpoint = tuple(
            (edge[0][coordinate] + edge[1][coordinate]) * inverse_two % p
            for coordinate in (0, 1)
        )
        difference = tuple(
            (edge[1][coordinate] - edge[0][coordinate]) * inverse_two % p
            for coordinate in (0, 1)
        )
        if midpoint[0] ** 2 % p == 1 and difference[0] ** 2 % p == 1:
            common_parameters.append(parameter)
            assert midpoint in {desired_midpoint, (30, 28)}
            assert difference in {desired_midpoint, (30, 28)}
    assert common_parameters == [0]


def test_record_keeps_the_global_cover_question_open():
    out = theorem_record(31)
    assert out["proved_all_claimed_statements"] is True
    assert out["arbitrary_prescribed_centers_always_cover_fixed_K"] is False
    assert out["existence_after_varying_K_excluded"] is False
    assert out["residual_ii_closed"] is False


def test_guard():
    with pytest.raises(ValueError):
        alpha_equal_beta_anchor_hall_theorem(9)
    with pytest.raises(ValueError):
        dependent_anchor_line_cover_theorem(23)
    with pytest.raises(ValueError):
        anchor_graph_pseudoforest_profile([(0, 0)])
    with pytest.raises(ValueError):
        anchor_graph_pseudoforest_profile([(0, 1), (1, 0)])
