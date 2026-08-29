from e1_gmin_m4_prop15713 import p7_positive_infinity_plus_seven_direction_reduction


def test_szonyi_reduces_positive_infinity_plus_seven_ledger():
    row = p7_positive_infinity_plus_seven_direction_reduction()
    assert row["proved_analytically"] is True
    assert row["uses_solver"] is False
    assert row["one_type_profile_count"] == 35
    assert row["projected_b_profile_count_before"] == 1217
    assert row["projected_b_profiles_excluded_here"] == 208
    assert row["projected_b_profile_count_after"] == 1009
    assert row["counts_residue_quotient_labelled_states"] is False
    assert row["positive_p7_infinity_plus_seven_closed"] is False


def test_exact_undetermined_direction_histograms():
    row = p7_positive_infinity_plus_seven_direction_reduction()
    assert row["undetermined_direction_histogram_before"] == {
        0: 217, 1: 300, 2: 280, 3: 210, 4: 126, 5: 56, 6: 21, 7: 6, 8: 1
    }
    assert row["excluded_undetermined_direction_histogram"] == {
        4: 126, 5: 56, 6: 21, 7: 4, 8: 1
    }
    assert row["remaining_undetermined_direction_histogram"] == {
        0: 217, 1: 300, 2: 280, 3: 210, 7: 2
    }
    assert row["negative_p7_infinity_plus_seven_changed"] is False
    assert row["top_level_gates_changed"] is False
