from e1_gmin_m4_prop15709 import p17_phase_one_residue_eight_exclusion


def test_all_phase_one_residue_eight_profiles_are_excluded():
    row = p17_phase_one_residue_eight_exclusion()
    assert row["proved_analytically"] is True
    assert row["uses_solver"] is False
    assert row["uses_new_arc_classification"] is False
    assert row["profile_count_before"] == 507
    assert row["profiles_excluded_by_global_gauge_identity"] == 66
    assert row["profiles_excluded_by_unique_even_fibre_identity"] == 214
    assert row["profiles_excluded_here"] == 280
    assert row["profile_count_after"] == 227
    assert row["minimum_remaining_pair_slack"] == 96
    assert row["all_phase_one_residue_eight_profiles_excluded"] is True
    assert row["all_survivors_have_phase_one_residue_zero"] is True


def test_full_rigid_anchor_sweep_has_exact_histograms():
    row = p17_phase_one_residue_eight_exclusion()
    assert row["rigid_phase_one_b2_lower_bound"] == 8
    assert row["rigid_phase_zero_b0_lower_bound_histogram"] == {
        3: 10,
        4: 27,
        5: 29,
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
        "u0=0,u1=0": 181,
        "u0=7,u1=0": 9,
        "u0=8,u1=0": 37,
    }
    assert sum(row["remaining_pair_slack_histogram"].values()) == 227
    assert min(row["remaining_pair_slack_histogram"]) == 96
    assert row["p17_second_all_finite_endpoint_closed"] is False
    assert row["top_level_gates_changed"] is False
