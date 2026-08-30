from e1_gmin_m4_prop15708 import (
    p17_slack_twenty_four_arithmetic_reduction,
    p17_slack_twenty_four_exclusion,
)


def test_slack_twenty_four_arithmetic_reduction_closes_142_profiles():
    row = p17_slack_twenty_four_arithmetic_reduction()
    assert row["proved_analytically"] is True
    assert row["slack_twenty_four_profiles_before"] == 151
    assert row["undetermined_direction_histogram_before"] == {
        0: 40,
        1: 47,
        2: 45,
        3: 19,
    }
    assert row["profiles_excluded_by_global_sign_identity"] == 142
    assert row["slack_twenty_four_profiles_after_arithmetic"] == 9
    assert row["profile_count_after_arithmetic"] == 878
    assert row["remaining_pair_slack_histogram_after_arithmetic"][24] == 9


def test_slack_twenty_four_rigid_direction_split_is_exact():
    row = p17_slack_twenty_four_arithmetic_reduction()
    assert row["residue_split"] == {"u0=0,u1=8": 142, "u0=8,u1=8": 9}
    assert row["phase_one_rigid_b2_lower_bound"] == 8
    assert row["easy_rigid_b0_lower_bound_histogram"] == {
        2: 3,
        3: 40,
        4: 59,
        5: 40,
    }
    assert row["easy_infinity_degree_candidate"] == 68
    assert row["hard_undetermined_direction_histogram"] == {2: 3, 3: 6}
    assert row["hard_infinity_degree_candidate"] == 4
    assert row["hard_forced_gauges_and_parallel_counts"] == {
        "phase_zero_b16": {"gauge": 1, "parallel_count": 7},
        "phase_one_b2": {"gauge": 0, "parallel_count": 0},
    }
    assert row["p17_second_all_finite_endpoint_closed"] is False
    assert row["top_level_gates_changed"] is False


def test_unique_even_fibre_identity_closes_the_nine_hard_profiles():
    row = p17_slack_twenty_four_exclusion()
    assert row["uses_solver"] is False
    assert row["uses_new_arc_classification"] is False
    assert row["hard_rigid_phase_zero_b16_lower_bounds"] == [
        3,
        3,
        3,
        3,
        2,
        2,
        2,
        2,
        2,
    ]
    assert row["profiles_excluded_by_unique_even_fibre_identity"] == 9
    assert row["profiles_excluded_here"] == 151
    assert row["slack_twenty_four_profiles_after"] == 0
    assert row["profile_count_after"] == 869
    assert 24 not in row["remaining_pair_slack_histogram"]
    certificate = row["unique_even_fibre_certificate"]
    assert certificate["phase_zero_finite_edges"] == 64
    assert certificate["phase_one_finite_edges"] == 1
    assert certificate["nonnegative_count_upper_bound"] == -3
    assert certificate["contradiction"] is True
    assert len(row["remaining_profile_indices"]) == 869
