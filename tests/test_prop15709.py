from e1_gmin_m4_prop15709 import p17_phase_one_residue_eight_exclusion


def test_all_phase_one_residue_eight_profiles_are_excluded():
    row = p17_phase_one_residue_eight_exclusion()
    assert row["proved_analytically"] is True
    assert row["uses_solver"] is False
    assert row["uses_new_arc_classification"] is False
    assert row["profile_count_before"] == 869
    assert row["profiles_excluded_by_global_gauge_identity"] == 334
    assert row["profiles_excluded_by_unique_even_fibre_identity"] == 214
    assert row["profiles_excluded_here"] == 548
    assert row["profile_count_after"] == 321
    assert row["minimum_remaining_pair_slack"] == 96
    assert row["all_phase_one_residue_eight_profiles_excluded"] is True
    assert row["all_survivors_have_phase_one_residue_zero"] is True
    assert row[
        "historical_orbiter_uncovered_slack_sixteen_profiles_received"
    ] == 74
    assert row[
        "historical_orbiter_uncovered_slack_sixteen_profiles_excluded_here"
    ] == 74
    assert row["slack_sixteen_block_closed_after_rigid_anchor_sweep"] is True


def test_full_rigid_anchor_sweep_has_exact_histograms():
    row = p17_phase_one_residue_eight_exclusion()
    assert row["rigid_phase_one_b2_lower_bound"] == 8
    assert row["rigid_phase_zero_b0_lower_bound_histogram"] == {
        2: 5,
        3: 60,
        4: 178,
        5: 91,
    }
    assert row["rigid_phase_zero_b16_lower_bound_histogram"] == {
        2: 4,
        3: 30,
        4: 36,
        5: 36,
        6: 36,
        7: 36,
        8: 36,
    }
    assert row["remaining_residue_pair_histogram"] == {
        "u0=0,u1=0": 275,
        "u0=7,u1=0": 9,
        "u0=8,u1=0": 37,
    }
    assert sum(row["remaining_pair_slack_histogram"].values()) == 321
    assert min(row["remaining_pair_slack_histogram"]) == 96
    assert row["p17_second_all_finite_endpoint_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert len(row["remaining_profile_indices"]) == 321
