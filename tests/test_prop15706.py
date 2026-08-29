from e1_gmin_m4_prop15706 import (
    p17_slack_zero_global_sign_certificate,
    p17_slack_zero_profile_exclusion,
)


def test_opposite_quadratic_types_force_impossible_infinity_degree():
    row = p17_slack_zero_global_sign_certificate()
    assert row["proved"] is True
    assert row["uses_solver"] is False
    assert row["opposite_type_comparison"] == "17*I=4+72*(g_+ + g_-)"
    assert row["infinity_degree_candidates_in_range"] == [68]
    assert row["possible_affine_odd_boundary_sizes"] == [66, 68, 70]
    assert row["required_affine_odd_boundary_size"] == 16


def test_both_slack_zero_profiles_are_excluded():
    row = p17_slack_zero_profile_exclusion()
    assert row["proved_analytically"] is True
    assert row["profile_count_before"] == 641
    assert row["profiles_excluded_here"] == 2
    assert row["profile_count_after"] == 639
    assert row["remaining_slack_zero_profiles"] == 0
    assert row["remaining_profiles_of_slack_at_least_twenty"] == 639
    assert min(row["remaining_pair_slack_histogram"]) == 20
