import pytest

from e1_gmin_m4_prop15732 import (
    block_triangle_syzygy_row,
    cleared_transition_exactness_row,
    endpoint_square_class_product_barrier_row,
    near_pairing_tangent_barrier_row,
    p31_hard_direction_component_barrier,
    potential_walk_row,
    proposition_15732,
    rich_direction_first_jet_row,
    weight_two_selection_span_row,
)


def test_cleared_transition_is_an_exact_potential_in_both_residues():
    c1 = cleared_transition_exactness_row(31)
    assert c1["c"] == 1
    assert c1["transition_quotient_degree"] == 3
    assert c1["repair_size_k"] == 22
    assert c1["envelope_degree_d"] == 22
    assert c1["vertex_potential_degree"] == 66
    assert c1["cleared_edge_degree"] == 66
    assert c1["c1_gauge"]["cycle_sum_still_a_coboundary"] is True
    assert c1["independent_linear_holonomy_obstruction"] is False
    assert all(c1["checks"].values())

    c2 = cleared_transition_exactness_row(41)
    assert c2["c"] == 2
    assert c2["transition_quotient_degree"] == 2
    assert c2["repair_size_k"] == 29
    assert c2["envelope_degree_d"] == 28
    assert c2["vertex_potential_degree"] == 86
    assert c2["cleared_edge_degree"] == 86
    assert c2["proved"] is True


def test_formal_closed_walk_telescopes_but_open_walk_does_not():
    triangle = potential_walk_row(("A", "B", "C", "A"))
    assert triangle["potential_coefficients_after_collection"] == {}
    assert triangle["telescopes_to_zero"] is True
    assert triangle["proved"] is True

    square = potential_walk_row(("A", "B", "C", "D", "A"))
    assert square["telescopes_to_zero"] is True

    open_walk = potential_walk_row(("A", "B", "C"))
    assert open_walk["potential_coefficients_after_collection"] == {"A": -1, "C": 1}
    assert open_walk["telescopes_to_zero"] is False
    assert open_walk["proved"] is True


def test_block_triangle_identity_has_the_residue_specific_degree():
    c1 = block_triangle_syzygy_row(31)
    assert c1["quotient_degree"] == 3
    assert c1["identity_degree"] == 6
    assert c1["identity_is_automatic_coboundary"] is True
    assert c1["endpoint_contradiction_from_identity_alone"] is False
    assert c1["proved"] is True

    c2 = block_triangle_syzygy_row(41)
    assert c2["quotient_degree"] == 2
    assert c2["identity_degree"] == 5
    assert "J(4,2)" in c2["support"]
    assert c2["proved"] is True


def test_rich_direction_quotient_has_a_nonzero_gauge_invariant_first_jet():
    for p in (31, 41, 43, 47):
        row = rich_direction_first_jet_row(p)
        assert row["edge_quotient_value"] == "Q_(a,z)(q)=0"
        assert row["first_jet_nonzero"] is True
        assert row["c1_first_jet_gauge_invariant"] is True
        assert row["phase_or_lift_formula_currently_proved"] is False
        assert all(row["checks"].values())
        assert row["proved"] is True


def test_near_pairing_directions_are_far_below_component_root_count():
    p31 = near_pairing_tangent_barrier_row(31)
    assert p31["near_pairing_profile"] == {
        "empty_fibres": 14,
        "singleton_fibres": 2,
        "double_fibres": 15,
    }
    assert p31["surviving_A_secant_floor"] == 5
    assert p31["A_tangent_ceiling"] == 12
    assert p31["envelope_degree"] == 22
    assert p31["roots_needed_to_force_direction_component"] == 23
    assert p31["root_deficit_at_least"] == 11
    assert p31["known_tangent_dual_points_are_distinct"] is True
    assert p31["intersection_multiplicity_used"] is False
    assert p31["direction_component_forced_by_known_distinct_roots"] is False

    for p in (17, 19, 23, 29, 31, 37, 41, 43, 47, 53):
        row = near_pairing_tangent_barrier_row(p)
        assert row["A_tangent_ceiling"] == row["R"] + 2
        assert row["root_deficit_at_least"] == row["R"] + 1
        assert all(row["checks"].values())


def test_every_p31_block_row_has_the_15728_component_barrier():
    for y in range(6):
        row = p31_hard_direction_component_barrier(y)
        assert row["15_728_nonrich_Paley_hard_direction_floor"] == 4 + y
        assert row["tangents_per_such_direction_at_most"] == 12
        assert row["envelope_degree"] == 22
        assert row["component_forced_in_any_such_direction"] is False
        assert row["proved"] is True


def test_pair_selection_masks_lose_triple_but_recover_quadruple_parity():
    triple = weight_two_selection_span_row(3)
    assert triple["generator_masks"] == ["110", "101", "011"]
    assert triple["span_rank"] == 2
    assert triple["span_size"] == 4
    assert triple["span_is_even_weight_subspace"] is True
    assert triple["full_block_mask"] == "111"
    assert triple["full_block_mask_in_span"] is False

    quadruple = weight_two_selection_span_row(4)
    assert quadruple["span_rank"] == 3
    assert quadruple["span_size"] == 8
    assert quadruple["span_is_even_weight_subspace"] is True
    assert quadruple["full_block_mask"] == "1111"
    assert quadruple["full_block_mask_in_span"] is True


def test_c2_trisecant_mask_barrier_has_the_exact_sole_target_exception():
    for y in range(7):
        row = endpoint_square_class_product_barrier_row(41, y)
        assert row["c"] == 2
        assert row["trisecants_x"] >= 1
        assert row["full_product_on_a_trisecant_recoverable_modulo_squares"] is False
        assert row["product_bridge_blocked_for_every_rich_target"] is (
            row["trisecants_x"] >= 2
        )
        assert row["product_bridge_blocked_for_some_rich_target"] is True
        assert row["sole_trisecant_target_exception"] is (
            row["trisecants_x"] == 1
        )
        assert row["proved"] is True

    all_quadruples = endpoint_square_class_product_barrier_row(31, 5)
    assert all_quadruples["c"] == 1
    assert all_quadruples["trisecants_x"] == 0
    assert all_quadruples["product_bridge_blocked_for_every_rich_target"] is False
    assert all_quadruples["product_bridge_blocked_for_some_rich_target"] is False


def test_proposition_package_is_a_proved_barrier_not_endpoint_progress():
    row = proposition_15732()
    assert row["prop"] == "15.732"
    assert row["result_status"] == "proved method barrier"
    assert row["finite_configuration_search_used"] is False
    assert row["nontrivial_linear_cycle_obstruction_available"] is False
    assert row["direction_component_from_15_728_profiles_available"] is False
    assert row["product_over_repairs_square_class_bridge_universally_available"] is False
    assert row["rich_direction_first_jet_proved"] is True
    assert row["phase_bridge_proved"] is False
    assert row["endpoint_excluded"] is False
    assert row["p_plus_one_shell_closed"] is False
    assert row["non_walsh_residual_ii_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True


def test_parameter_validation_is_strict():
    with pytest.raises(ValueError):
        potential_walk_row(("A",))
    with pytest.raises(ValueError):
        potential_walk_row(("A", ""))
    for invalid in (1, True, 2.0):
        with pytest.raises(ValueError):
            weight_two_selection_span_row(invalid)
