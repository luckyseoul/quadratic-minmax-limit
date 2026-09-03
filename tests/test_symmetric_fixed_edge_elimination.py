from __future__ import annotations

import pytest

from e1_gmin_m4_symmetric_fixed_edge_elimination import (
    directionwise_parallel_slices,
    exact_elimination_replay,
    exact_fixed_word_design_replay,
    exact_mobius_midpoint_replay,
    fixed_edge_elimination_theorem,
    fixed_word_block_basis_theorem,
    forced_fixed_word_parity_theorem,
    hamming_slice_identity,
    mobius_midpoint_direction_theorem,
    orbit_fixed_word,
    p31_mobius_cancellation_parity_ladder,
    theorem_record,
)


def test_symbolic_block_elimination_is_all_prime():
    for p in (3, 5, 7, 11, 31, 43):
        h = (p - 1) // 2
        d = p + 1
        out = fixed_edge_elimination_theorem(p)
        assert out["proved"] is True
        assert out["block_form"] == "R_plus=[[A,2B],[0,C]]"
        assert out["fixed_antipodal_variables"] == d * h
        assert out["compatible_fixed_residue_rank"] == d * h
        assert out["fixed_map_mod2_rank"] == d * h
        assert out["fixed_map_mod2_isomorphism"] is True
        assert out["remaining_boolean_variables"] == (
            "unused nonfixed inversion orbits only"
        )
        assert out["restricted_symmetric_fibre_closed"] is False
        assert out["residual_ii_closed"] is False


def test_small_exact_maps_have_the_claimed_blocks_and_unique_parity_solution():
    for p in (3, 7):
        h = (p - 1) // 2
        d = p + 1
        out = exact_elimination_replay(p)
        assert out["proved"] is True
        assert out["role"] == (
            "fail-when-wrong small-prime replay, not theorem evidence"
        )
        assert out["lower_left_block_zero"] is True
        assert out["upper_right_block_even"] is True
        assert out["fixed_map_mod2_rank"] == d * h
        assert out["fixed_binary_vector_recovered_uniquely"] is True
        assert out["explicit_inverse_formula_recovers_fixed_vector"] is True
        assert out["fixed_remainder_even"] is True
        assert out["divided_target_reconstructed_exactly"] is True
        assert out["hamming_slice"]["proved"] is True
        assert out["directionwise_parallel_slices"][
            "all_direction_slices_feasible"
        ] is True
        assert out["directionwise_parallel_slices"][
            "global_hamming_slice_recovered_by_summing_directions"
        ] is True


def test_hamming_slice_counts_single_fixed_and_double_orbit_edges_exactly():
    out = hamming_slice_identity(
        used_orbit_count=9,
        fixed_edge_bits=(1, 0, 1, 1),
        unused_double_orbit_bits=(0, 1, 1, 0, 1),
    )
    assert out["fixed_antipodal_edge_weight"] == 3
    assert out["unused_double_orbit_weight"] == 3
    assert out["physical_edge_count"] == 18
    assert out["slice_left"] == out["slice_right"] == 6
    assert out["proved"] is True


def test_fixed_word_parity_detects_odd_used_parallel_orbit_count():
    feasible = forced_fixed_word_parity_theorem(31, 481, 480, 16)
    assert feasible["graph_target_fixed_word_parity"] == 1
    assert feasible["used_nonparallel_orbit_count"] == 464
    assert feasible["central_remainder_fixed_word_parity"] == 1
    assert feasible["hamming_slice_numerator_even_automatically"] is False
    assert feasible["hamming_slice_numerator_parity"] == 0
    assert feasible["hamming_slice_parity_feasible"] is True
    assert feasible["compact_fixed_cell_word_may_be_omitted"] is False
    assert feasible["parity_excludes_symmetric_completion"] is False

    excluded = forced_fixed_word_parity_theorem(31, 481, 479, 15)
    assert excluded["used_nonparallel_orbit_count"] == 464
    assert excluded["central_remainder_fixed_word_parity"] == 1
    assert excluded["hamming_slice_numerator_parity"] == 1
    assert excluded["hamming_slice_parity_feasible"] is False
    assert excluded["parity_excludes_symmetric_completion"] is True

    feasible_even_remainder = forced_fixed_word_parity_theorem(
        31, 481, 479, 14
    )
    assert feasible_even_remainder["central_remainder_fixed_word_parity"] == 0
    assert feasible_even_remainder["hamming_slice_parity_feasible"] is True

    with pytest.raises(ValueError, match="cannot exceed"):
        forced_fixed_word_parity_theorem(31, 4, 5, 0)
    with pytest.raises(ValueError, match="parallel.*cannot exceed"):
        forced_fixed_word_parity_theorem(31, 4, 4, 5)


def test_p31_mobius_cancellation_parity_is_automatic_and_minimum_is_rigid():
    for t, zero_phi_cancellations in ((68, 8), (162, 4), (177, 1)):
        kappa = 178 - t
        out = p31_mobius_cancellation_parity_ladder(
            t, kappa, zero_phi_cancellations
        )
        assert out["graph_edge_count"] == 125 + 2 * t
        assert out["used_orbit_count"] == out["graph_edge_count"] - 1
        assert out["used_zero_phi_orbit_count"] == (
            16 - 2 * zero_phi_cancellations
        )
        assert out["used_nonzero_phi_orbit_count"] % 2 == 0
        assert out["used_zero_phi_orbit_count_parity"] == 0
        assert out["used_nonzero_phi_orbit_count_parity"] == 0
        assert out["fixed_word_weight_parity"] == 1
        assert out["hamming_slice_numerator_parity"] == 0
        assert out["hamming_slice_parity_is_automatic_for_sixteen_halves"]
        assert out["at_minimum_cancellation"] is True
        assert out["conditional_fixed_word_weight_if_completion"] == 1
        assert out["conditional_unused_double_orbit_count_if_completion"] == 0
        assert out["parity_excludes_mobius_completion"] is False
        assert out["residual_ii_closed"] is False

    above = p31_mobius_cancellation_parity_ladder(177, 2)
    assert above["remaining_edge_capacity"] == 3
    assert above["fixed_word_weight_parity"] == 1
    assert above["conditional_fixed_word_weight_if_completion"] is None
    assert above["conditional_unused_double_orbit_count_if_completion"] is None

    below = p31_mobius_cancellation_parity_ladder(68, 109)
    assert below["remaining_edge_capacity"] == -1
    assert below["support_size_feasible"] is False


def test_p31_mobius_cancellation_parity_validates_the_split_ledger():
    with pytest.raises(ValueError, match="68<=t<=177"):
        p31_mobius_cancellation_parity_ladder(67, 111)
    with pytest.raises(ValueError, match="between 0 and 240"):
        p31_mobius_cancellation_parity_ladder(177, 241)
    with pytest.raises(ValueError, match="between 0 and 8"):
        p31_mobius_cancellation_parity_ladder(177, 9, 9)
    with pytest.raises(ValueError, match="cannot exceed total"):
        p31_mobius_cancellation_parity_ladder(177, 1, 2)
    with pytest.raises(ValueError, match="cannot exceed 232"):
        p31_mobius_cancellation_parity_ladder(177, 233, 0)


def test_parallel_coordinates_pin_one_constant_weight_slice_per_direction():
    p = 7
    target = (5, 4, 3, 2, 7, 6, 1, 8)
    used = (1, 0, 1, 0, 1, 0, 1, 0)
    fixed = (0, 0, 0, 0, 0, 0, 0, 0)
    out = directionwise_parallel_slices(p, target, used, fixed)
    assert out["slice_numerators"] == [4, 4, 2, 2, 6, 6, 0, 8]
    assert out["selected_unused_double_orbits"] == [2, 2, 1, 1, 3, 3, 0, 4]
    assert out["all_direction_slices_feasible"] is True
    assert out["global_hamming_slice_recovered_by_summing_directions"] is True
    assert out["normalization_reason"] == (
        "for an edge parallel to L, epsilon_L*tau_e=1"
    )

    parity_failure = directionwise_parallel_slices(
        p,
        (4, 4, 3, 2, 7, 6, 1, 8),
        used,
        fixed,
    )
    assert parity_failure["all_direction_slices_feasible"] is False
    assert parity_failure["slice_integral"][0] is False


def test_one_orbit_toggles_zero_or_one_affine_p_block():
    transverse = orbit_fixed_word(7, ((1, 6), (1, 1)))
    assert transverse["midpoint"] == [1, 0]
    assert transverse["midpoint_parallel_to_difference"] is False
    assert transverse["fixed_word_weight"] == 7
    assert transverse["independent_of_selected_orbit_side"] is True
    assert transverse["proved"] is True

    negative = orbit_fixed_word(7, ((6, 1), (6, 6)))
    assert negative["fixed_word_support"] == transverse["fixed_word_support"]

    parallel = orbit_fixed_word(7, ((3, 0), (6, 0)))
    assert parallel["midpoint_parallel_to_difference"] is True
    assert parallel["fixed_word_support"] == []
    assert parallel["fixed_word_weight"] == 0
    assert parallel["proved"] is True


def test_paired_affine_line_words_form_an_orthogonal_binary_basis():
    for p in (3, 5, 7, 11, 31):
        h = (p - 1) // 2
        classes = (p + 1) * h
        out = fixed_word_block_basis_theorem(p)
        assert out["proved"] is True
        assert out["paired_affine_line_block_types"] == classes
        assert out["points_per_block"] == p
        assert out["blocks_per_point"] == p
        assert out["binary_gram_identity"] == "M*M^T=I, hence M^T*M=I"
        assert out["block_vectors_form_basis"] is True
        assert out["zero_word_orbits"] == classes * h
        assert out["multiplicity_per_nonzero_block_type"] == p * h
        assert out["disjoint_C_kernel_lifts_per_block_type"] == h
        assert out["columns_per_C_kernel_lift"] == p


def test_small_exact_fixed_words_match_the_formula_and_gram_counts():
    for p in (3, 7):
        h = (p - 1) // 2
        classes = (p + 1) * h
        out = exact_fixed_word_design_replay(p)
        assert out["proved"] is True
        assert out["role"] == (
            "fail-when-wrong small-prime replay, not theorem evidence"
        )
        assert out["direct_inverse_matches_affine_line_formula"] is True
        assert out["nonzero_block_types"] == classes
        assert out["exact_pair_intersection_numbers"] is True
        assert out["binary_gram_is_identity"] is True
        assert out["binary_block_rank"] == classes
        assert out["disjoint_C_kernel_lifts_per_block_type"] == h
        assert out["columns_per_C_kernel_lift"] == p
        assert out["exact_C_kernel_lifts"] is True


def test_one_mobius_half_has_two_midpoints_in_each_of_h_directions():
    for p in (3, 5, 7, 11, 31):
        out = mobius_midpoint_direction_theorem(p)
        assert out["proved"] is True
        assert out["midpoint_slope_M_over_L"] == "1-1/(t+1)^2"
        assert out["vertical_midpoint_direction_occurs"] is False
        assert out["distinct_midpoint_directions"] == (p - 1) // 2
        assert out["parameters_per_midpoint_direction"] == 2
        assert out["one_half_hits_any_midpoint_direction_at_most"] == 2

    for p in (3, 7):
        replay = exact_mobius_midpoint_replay(p)
        assert replay["proved"] is True
        assert replay["midpoint_slope_formula_holds"] is True
        assert replay["distinct_midpoint_directions"] == (p - 1) // 2
        assert replay["multiplicity_set"] == [2]
        assert replay["role"] == (
            "fail-when-wrong small-prime replay, not theorem evidence"
        )


def test_theorem_record_preserves_the_open_reduced_fibre():
    out = theorem_record(31)
    assert out["proved_all"] is True
    assert out["proved"]["block_triangular_form"] is True
    assert out["proved"]["unique_fixed_edge_parity_vector"] is True
    assert out["proved"]["explicit_fixed_edge_inverse"] is True
    assert out["proved"]["divided_target_equivalence"] is True
    assert out["proved"]["exact_hamming_slice"] is True
    assert out["proved"]["forced_fixed_word_parity_identity"] is True
    assert out["proved"]["p31_mobius_cancellation_parity_ladder"] is True
    assert out["p31_mobius_cancellation_parity"][
        "conditional_fixed_word_weight_if_completion"
    ] == 1
    assert out["p31_mobius_cancellation_parity"][
        "conditional_unused_double_orbit_count_if_completion"
    ] == 0
    assert out["proved"]["directionwise_parallel_slices"] is True
    assert out["proved"]["per_orbit_affine_fixed_word"] is True
    assert out["proved"]["fixed_word_blocks_form_binary_basis"] is True
    assert out["proved"]["disjoint_C_kernel_lifts"] is True
    assert out["proved"]["mobius_midpoint_direction_multiplicity_two"] is True
    assert out["proved"]["support_weight_coset_reduction"] is True
    assert out["proved"]["reduced_unused_orbit_fibre_nonempty"] is False
    assert out["proved"]["residual_ii_closed"] is False


def test_parameter_and_binary_guards():
    with pytest.raises(ValueError):
        fixed_edge_elimination_theorem(9)
    with pytest.raises(ValueError):
        exact_elimination_replay(5)
    with pytest.raises(ValueError):
        exact_elimination_replay(11)
    with pytest.raises(ValueError):
        exact_fixed_word_design_replay(5)
    with pytest.raises(ValueError):
        exact_mobius_midpoint_replay(5)
    with pytest.raises(ValueError):
        orbit_fixed_word(7, ((1, 0), (6, 0)))
    with pytest.raises(ValueError):
        hamming_slice_identity(-1, (), ())
    with pytest.raises(ValueError):
        hamming_slice_identity(0, (2,), ())
    with pytest.raises(ValueError):
        hamming_slice_identity(0, (), (True,))
    with pytest.raises(ValueError):
        directionwise_parallel_slices(7, (0,), (0,), (0,))
    with pytest.raises(ValueError):
        directionwise_parallel_slices(5, (0,) * 6, (0,) * 6, (0,) * 6)
