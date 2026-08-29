from e1_gmin_m4_prop15689 import (
    p19_complete_arc_spectrum,
    p19_low_slack_geometric_exclusion,
    p19_low_slack_profile_ledger,
)


def test_p19_complete_arc_spectrum_has_the_needed_gap():
    row = p19_complete_arc_spectrum()
    assert row["complete_arc_sizes"] == [10, 11, 12, 13, 14, 20]
    assert row["no_complete_arc_sizes"] == [15, 16, 17, 18, 19]
    assert row["unique_size_twenty_arc"] == "nondegenerate conic"


def test_p19_low_slack_profile_block_is_exact():
    row = p19_low_slack_profile_ledger()
    assert row["profile_counts_by_slack"] == {0: 54, 4: 37, 8: 25, 12: 13}
    assert row["excluded_profile_count"] == 129
    assert row["slack_zero_profiles_with_at_least_three_undetermined"] == 25
    assert row["slack_zero_profiles_with_one_or_two_undetermined"] == 29
    assert row["small_t_arc_minimum_high_nonundetermined_b"] >= 6
    assert row["minimum_undetermined_directions"] == {8: 2, 12: 2}


def test_p19_low_slack_conic_reduction_leaves_fourteen_profiles():
    row = p19_low_slack_geometric_exclusion()
    assert row["profile_count_before"] == 143
    assert row["profile_count_excluded"] == 129
    assert row["profile_count_after"] == 14
    assert row["remaining_pair_slack_histogram"] == {
        16: 7,
        20: 4,
        24: 1,
        28: 1,
        32: 1,
    }
    assert row["p19_second_all_finite_endpoint_closed"] is False
    assert row["top_level_gates_changed"] is False
