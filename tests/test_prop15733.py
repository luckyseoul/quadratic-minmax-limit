import pytest

from e1_gmin_m4_prop15733 import (
    hard_b2_coefficient_row,
    hard_mean30_b2_upgrade,
    infinity_star_profile_options,
    odd_g_norm_obstruction,
    opposite_type_final_contradiction_row,
    p31_baseline_coefficient_rules,
    p31_block_direction_upgrade_row,
    phase_zero_filter_rows,
    pre_phase_zero_candidate_rows,
    proposition_15733,
    short_s8_parallel_rows,
    short_type_collapse_rows,
    surviving_global_rows,
)


def test_baseline_coefficient_offsets_have_the_required_one_unit_gaps():
    rules = p31_baseline_coefficient_rules()
    assert rules["hard_phase_one_b2"]["offset"] == 4
    assert rules["hard_phase_one_b30"]["offset"] == 3
    assert rules["opposite_phase_zero_b2"]["offset"] == 4
    assert rules["opposite_phase_zero_b30"]["offset"] == 5
    assert rules["proved"] is True


def test_all_fifteen_hard_mean30_directions_are_b2():
    row = hard_mean30_b2_upgrade()
    assert row["previous_b2_floor"] == 14
    assert row["mean_30_direction_count"] == 15
    assert row["equal_means_force_equal_parallel_counts"] is True
    assert row["coexistence_possible"] is False
    assert row["all_mean_30_directions_have_b2"] is True
    assert row["proved"] is True


def test_hard_b2_coefficient_row_has_exact_transverse_capacity():
    for I, parallel, g in ((28, 6, 2), (60, 4, 4), (92, 2, 6)):
        row = hard_b2_coefficient_row(I, parallel)
        assert row["g"] == g
        assert row["g_integral_nonnegative"] is True
        assert row["transverse_edge_capacity"] == 121 - 15 * g
        assert row["capacity_formula_121_minus_15g"] == 121 - 15 * g
        assert row["proved"] is True

    corrupt = hard_b2_coefficient_row(28, 5)
    assert corrupt["g_integral_nonnegative"] is False
    assert corrupt["proved"] is False


def test_odd_g_entrywise_norm_is_impossible():
    row = odd_g_norm_obstruction()
    assert row["pair_count"] == 465
    assert row["maximum_zero_pairs"] == 240
    assert row["base_entrywise_norm_floor"] == 225
    assert row["corrected_entrywise_norm_floor"] == 224
    assert row["transverse_capacity_ceiling"] == 106
    assert row["contradiction"] is True
    assert row["g_must_be_even"] is True
    assert row["proved"] is True


def test_short_type_argument_forces_s_equal_eight():
    rows = short_type_collapse_rows()
    assert [row["s"] for row in rows] == [0, 2, 4, 6, 8]
    assert [row["minimum_direction_mean"] for row in rows] == [16, 12, 8, 4, 32]
    assert [row["total_parallel_excess_above_minimum"] for row in rows] == [
        8,
        10,
        12,
        14,
        0,
    ]
    assert [row["s"] for row in rows if not row["excluded"]] == [8]


def test_short_s8_rows_all_close_by_boundary_or_phase_zero_congruence():
    rows = short_s8_parallel_rows()
    assert [(row["P"], row["I"]) for row in rows] == [
        (0, 124),
        (1, 108),
        (2, 92),
        (3, 76),
        (4, 60),
        (5, 44),
        (6, 28),
        (7, 12),
    ]
    assert rows[0]["boundary_support_excluded"] is True
    assert all(row["phase_zero_floor_compatible_b"] == [0, 2, 30] for row in rows)
    assert all(row["excluded"] for row in rows)
    for row in rows:
        if not row["boundary_support_excluded"]:
            assert row["b2_congruence_I_minus_4"] is False
            assert row["b30_congruence_I_minus_5"] is False
            assert row["all_opposite_directions_have_b0"] is True
            assert row["required_high_hard_direction_b"] == 42
            assert row["contradiction"] is True


def test_aggregate_and_boundary_bounds_leave_exactly_twelve_rows():
    rows = pre_phase_zero_candidate_rows()
    assert [(row["I"], row["P"], row["r"]) for row in rows] == [
        (0, 4, 0),
        (2, 2, 0),
        (4, 0, 0),
        (28, 6, 1),
        (30, 4, 1),
        (32, 2, 1),
        (34, 0, 1),
        (60, 4, 2),
        (62, 2, 2),
        (64, 0, 2),
        (92, 2, 3),
        (94, 0, 3),
    ]
    assert all(abs(row["signed_cell_aggregate"]) <= row["transverse_capacity"] for row in rows)


def test_phase_zero_residue_leaves_only_three_global_rows():
    rows = phase_zero_filter_rows()
    assert [(row["I"], row["P"], row["phase_zero_u0"]) for row in rows] == [
        (0, 4, 4),
        (2, 2, 6),
        (4, 0, 8),
        (28, 6, 0),
        (30, 4, 2),
        (32, 2, 4),
        (34, 0, 6),
        (60, 4, 0),
        (62, 2, 2),
        (64, 0, 4),
        (92, 2, 0),
        (94, 0, 2),
    ]
    survivors = [row for row in rows if not row["excluded"]]
    assert [(row["I"], row["P"], row["r"]) for row in survivors] == [
        (28, 6, 1),
        (60, 4, 2),
        (92, 2, 3),
    ]
    assert all(row["forced_zero_quotient_direction_mean"] <= 16 for row in rows if row["excluded"])
    assert all(row["nonzero_integral_lift_floor"] == 28 for row in rows)


def test_three_rows_exhaust_finite_edges_in_the_hard_type():
    rows = surviving_global_rows()
    assert [(row["I"], row["baseline_hard_parallel_count_P"], row["r"]) for row in rows] == [
        (28, 6, 1),
        (60, 4, 2),
        (92, 2, 3),
    ]
    for row in rows:
        assert row["hard_parallel_total"] == row["finite_edge_count"]
        assert row["every_finite_selected_edge_has_hard_sign"] is True
        assert row["every_opposite_direction_parallel_count"] == 0
        assert row["every_opposite_direction_mean"] == 32
        assert all(row["checks"].values())
        assert row["proved"] is True


def test_optional_infinity_star_profiles_are_exact():
    assert infinity_star_profile_options(1)["infinity_endpoint_count_histograms"] == [
        {0: 3, 1: 28}
    ]
    assert infinity_star_profile_options(2)["infinity_endpoint_count_histograms"] == [
        {0: 1, 2: 30},
        {1: 2, 2: 29},
    ]
    assert infinity_star_profile_options(3)["infinity_endpoint_count_histograms"] == [
        {2: 1, 3: 30}
    ]
    for r in (1, 2, 3):
        row = infinity_star_profile_options(r)
        assert row["deficit_sum"] == 4 - r
        assert row["negative_deficit_impossible"] is True
        assert row["proved"] is True


def test_each_surviving_row_forces_the_impossible_b42_direction():
    for I in (28, 60, 92):
        row = opposite_type_final_contradiction_row(I)
        assert row["phase_zero_floor_compatible_b"] == [0, 2, 30]
        assert row["phase_zero_b2_congruence_holds"] is False
        assert row["phase_zero_b30_congruence_holds"] is False
        assert row["all_opposite_directions_have_b0"] is True
        assert row["global_sum_b"] == 72
        assert row["hard_baseline_b_sum"] == 30
        assert row["required_high_hard_direction_b"] == 42
        assert row["maximum_possible_even_b"] == 30
        assert row["contradiction"] is True
        assert row["proved"] is True


def test_nonrich_hard_direction_floor_improves_to_five_plus_y():
    for y in range(6):
        row = p31_block_direction_upgrade_row(y)
        assert row["hard_b2_direction_count"] == 15
        assert row["nonrich_hard_b2_direction_floor"] == 5 + y
        assert row["proved"] is True


def test_proposition_package_closes_only_the_p31_endpoint():
    row = proposition_15733()
    assert row["prop"] == "15.733"
    assert row["result_status"] == "proved p=31 endpoint exclusion"
    assert row["finite_configuration_search_used"] is False
    assert row["p31_R10_endpoint_excluded"] is True
    assert row["p31_first_possible_positive_slack"] == 11
    assert row["first_unexcluded_endpoint_prime"] == 37
    assert row["endpoint_all_primes_closed"] is False
    assert row["p_plus_one_shell_closed"] is False
    assert row["non_walsh_residual_ii_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True


def test_parameter_validation_is_strict():
    for values in ((-1, 2), (28, -1), (True, 2), (28, False), (124, 2)):
        with pytest.raises(ValueError):
            hard_b2_coefficient_row(*values)
    for invalid in (0, 4, True):
        with pytest.raises(ValueError):
            infinity_star_profile_options(invalid)
    with pytest.raises(ValueError):
        opposite_type_final_contradiction_row(30)
