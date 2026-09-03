import pytest

from e1_gmin_m4_conic_odd_radon import (
    conic_reduction_constants,
    nonequianharmonic_constant_fiber_no_go,
    nonequianharmonic_score_three_compact_candidates,
    p31_equianharmonic_witness_certificate,
    scaled_family_exceptional_row_obstruction,
    tangent_conic_target,
    theorem_record,
)


def test_all_prime_conic_peeling_and_character_bounds_use_exact_constants():
    for p, b in ((31, 7), (43, 9), (47, 11), (59, 14)):
        row = conic_reduction_constants(p, b)
        assert row["proved"]
        assert row["support_at_most_3m"]
        assert row["peeling_contradiction"]
        assert row["no_constant_character_bound_strict"]
        assert row["one_constant_character_bound_strict"]
        assert row["beta_nonzero_excluded_by_l1"]
        assert row["normal_form_Omega_points"] == p - 2


def test_conic_reduction_rejects_parameters_outside_its_proved_scope():
    with pytest.raises(ValueError):
        conic_reduction_constants(29, 7)
    with pytest.raises(ValueError):
        conic_reduction_constants(31, 8)
    with pytest.raises(ValueError):
        conic_reduction_constants(41, 7)


def test_tangent_conic_parameterization_has_exact_orbit_support():
    for p, k in ((31, 11), (43, 13), (47, 2)):
        target = tangent_conic_target(p, k)
        assert len(target) == p - 2
        assert set(target.values()) <= {-1, 1}
        assert tangent_conic_target(p, -k) == target


def test_nonequianharmonic_score_three_candidates_and_global_counting_no_go():
    for p in (31, 43, 47, 59):
        r = (p - 3) // 4
        candidates = nonequianharmonic_score_three_compact_candidates(p, 2)
        assert len(candidates) <= 2
        for b in (0, r // 2, r):
            row = nonequianharmonic_constant_fiber_no_go(p, b, 2)
            assert row["proved"]
            assert row["q_cubed"] != 1
            assert row["score_three_compact_candidate_count"] <= 2
            assert row["total_atom_count"] <= row["total_atom_count_upper_bound"]
            assert row["score_upper_bound"] < row["target_score"]
            assert row["target_score"] == p - 2

    with pytest.raises(ValueError):
        nonequianharmonic_constant_fiber_no_go(31, 7, 11)


def test_p31_equianharmonic_witness_replays_edges_and_all_odd_channels():
    row = p31_equianharmonic_witness_certificate()
    assert row["proved"]
    assert row["k_squared_mod_p"] == row["minus_three_mod_p"]
    assert row["q_order"] == 3
    assert row["q_satisfies_equianharmonic_polynomial"]
    assert row["ae_atoms_are_q_cycle_triangles"]
    assert row["edge_orbit_replay_exact"]
    assert row["target_support"] == row["target_l1"] == 29
    assert row["odd_channel_count"] == 105
    assert row["all_odd_channels_zero"]
    assert row["degree_six"] == [11, 19, 10]
    assert row["degree_eight"] == [12, 11, 23, 6]
    assert not row["degree_six_and_eight_both_zero"]
    assert not row["central_signed_chain"]


def test_scaled_family_obstruction_and_paley_half_scope_are_separated():
    for p in (31, 43, 47):
        row = scaled_family_exceptional_row_obstruction(p)
        assert row["proved"]
        assert row["projective_direction_count"] > row["identity_degree"]
        assert not row["constant_nonzero_quadratic_character_on_P1_possible"]
        assert not row["actual_signed_Paley_exceptional_row_proved"]
        assert "Paley type" in row["paley_half_exception"]


def test_theorem_record_keeps_the_global_and_boolean_gates_open():
    record = theorem_record()
    proved = record["proved"]
    assert proved["conic_containing_low_weight_word_is_fully_conic_supported"]
    assert proved["high_intersection_irreducible_conic_is_triangle_tangent"]
    assert proved["p31_b7_equianharmonic_odd_zero_atom_witness_exists"]
    assert proved["nonequianharmonic_constant_branch_is_excluded"]
    assert proved["constant_branch_forces_q_cubed_equals_one"]
    assert proved["constant_branch_forces_p_congruent_7_mod_12"]
    assert not proved["p31_witness_degree_six_and_eight_both_zero"]
    assert not proved["nonstar_constant_branch_is_excluded"]
    assert not proved["actual_signed_Paley_exceptional_row_proved"]
    assert proved["paley_half_norm_form_coordination_is_algebraically_compatible"]
    assert not proved["common_Fp_atom_lift_constructed"]
    assert not proved["Boolean_lift_constructed"]
    assert not proved["residual_ii_closed"]
    assert record["L_status"] == "OPEN"
