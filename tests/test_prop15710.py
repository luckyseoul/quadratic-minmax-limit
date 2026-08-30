from e1_gmin_m4_prop15710 import p17_phase_one_b16_global_sign_reduction


def test_complementary_global_sign_sweep_leaves_nineteen_profiles():
    row = p17_phase_one_b16_global_sign_reduction()
    assert row["proved_analytically"] is True
    assert row["uses_solver"] is False
    assert row["profile_count_before"] == 321
    assert row["profiles_excluded_by_b0_b16_comparison"] == 270
    assert row["profiles_excluded_by_b16_b16_comparison"] == 32
    assert row["profiles_excluded_here"] == 302
    assert row["profile_count_after"] == 19
    assert row["remaining_pair_slack_histogram"] == {
        96: 3,
        100: 4,
        104: 4,
        108: 3,
        112: 3,
        116: 1,
        128: 1,
    }
    assert row["remaining_residue_pair_histogram"] == {
        "u0=0,u1=0": 5,
        "u0=7,u1=0": 9,
        "u0=8,u1=0": 5,
    }


def test_complementary_global_sign_arithmetic_is_exact():
    row = p17_phase_one_b16_global_sign_reduction()
    assert row["rigid_phase_one_b16_lower_bound"] == 9
    assert row["rigid_phase_zero_b0_lower_bound_histogram"] == {
        1: 16,
        2: 40,
        3: 69,
        4: 90,
        5: 51,
        6: 4,
    }
    assert row["rigid_phase_zero_b16_lower_bound_histogram"] == {
        1: 4,
        2: 4,
        3: 4,
        4: 4,
        5: 4,
        6: 4,
        7: 4,
        8: 4,
    }
    assert row["b0_b16_contradiction"] == {
        "infinity_degree": 60,
        "forced_gauge_sum": 14,
        "minimum_gauges": [8, 7],
        "minimum_gauge_sum": 15,
    }
    assert row["b16_b16_contradiction"] == {
        "infinity_degree": 68,
        "forced_gauge_sum": 16,
        "minimum_gauges": [9, 8],
        "minimum_gauge_sum": 17,
    }
    assert row["p17_second_all_finite_endpoint_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert len(row["remaining_profile_indices"]) == 19
